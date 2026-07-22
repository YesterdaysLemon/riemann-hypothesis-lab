"""Explicit prime-power transition searches over finite Weil matrices.

Version 2 is an orchestration layer over :mod:`riemann_lab.weil_search`.
It adds exact prime-power stencil contracts and per-cell policies without
changing the audited v1 matrix, attempt, or cell state machines.  Every
artifact remains exploratory and leaves the Riemann Hypothesis unresolved.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from .artifacts import content_sha256
from . import weil_search as v1


PLAN_SCHEMA = "rh-lab/weil-transition-plan/v2"
RUN_SCHEMA = "rh-lab/weil-transition-run/v2"
INDEX_SCHEMA = "rh-lab/weil-transition-index/v2"
CONTRACT_SCHEMA = "rh-lab/weil-transition-cell-contract/v2"
ENGINE_ID = "finite-weil-prime-transition-search-v2"

FROZEN_Q7_Q9_BATCH = {
    "batch_id": "weil-transition-q7-q9-v2",
    "cell_count": "81",
    "prime_powers": [
        {"q": "7", "prime": "7", "exponent": "1"},
        {"q": "8", "prime": "2", "exponent": "3"},
        {"q": "9", "prime": "3", "exponent": "2"},
    ],
    "stencil_exponents": ["4", "6", "8", "10"],
    "degree_ladders": {
        "7": ["12", "16", "20"],
        "8": ["16", "20", "24"],
        "9": ["16", "20", "24"],
    },
    "stencil_shape": (
        "four exact below/above pairs plus one unparameterized exact point "
        "for every q and degree"
    ),
}

EXPLORATORY_LIMITATION = (
    "This is a bounded finite-dimensional search at explicitly selected "
    "prime-power transition stencils. Positive cells do not prove RH, and "
    "negative candidates remain quarantined pending independent mathematical "
    "and backend audits. The global Riemann Hypothesis status is UNRESOLVED."
)


class WeilTransitionError(ValueError):
    """Base class for v2 transition-search failures."""


class WeilTransitionPlanError(WeilTransitionError):
    """Raised when a transition plan is malformed or noncanonical."""


class WeilTransitionCheckpointError(WeilTransitionError):
    """Raised when checkpoint state is corrupt or incompatible."""


class WeilTransitionVerificationError(WeilTransitionError):
    """Raised when replay of transition evidence fails."""


def _parse_integer(value: Any, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool):
        raise WeilTransitionPlanError(f"{name} must be an integer")
    if isinstance(value, int):
        parsed = value
    elif isinstance(value, str):
        try:
            parsed = int(value)
        except ValueError as exc:
            raise WeilTransitionPlanError(f"{name} must be an integer") from exc
        if value != str(parsed):
            raise WeilTransitionPlanError(
                f"{name} must use canonical decimal syntax"
            )
    else:
        raise WeilTransitionPlanError(f"{name} must be an integer")
    if parsed < minimum:
        raise WeilTransitionPlanError(f"{name} must be at least {minimum}")
    return parsed


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    for divisor in range(3, isqrt(value) + 1, 2):
        if value % divisor == 0:
            return False
    return True


def _is_prime_power(value: int) -> bool:
    """Return whether an integer is p^r for a prime p and r >= 1."""

    if value < 2:
        return False
    for prime in range(2, isqrt(value) + 1):
        if value % prime or not _is_prime(prime):
            continue
        remainder = value
        while remainder % prime == 0:
            remainder //= prime
        return remainder == 1
    return _is_prime(value)


def _expected_cutoff(q: int, position: str, stencil_exponent: int | None) -> Fraction:
    if position == "exact":
        if stencil_exponent is not None:
            raise WeilTransitionPlanError("exact transition cells cannot carry j")
        return Fraction(q, 1)
    if position not in {"below", "above"}:
        raise WeilTransitionPlanError("position must be below, exact, or above")
    if stencil_exponent is None:
        raise WeilTransitionPlanError("neighbor transition cells require j")
    power = 1 << stencil_exponent
    if position == "below":
        return Fraction(q * power, power + 1)
    return Fraction(q * (power + 1), power)


def _formula(position: str) -> str:
    return {
        "below": "q*2^j/(2^j+1)",
        "exact": "q",
        "above": "q*(2^j+1)/2^j",
    }[position]


def _validate_no_crossing(q: int, cutoff: Fraction) -> None:
    if cutoff < q:
        candidates = range((cutoff.numerator + cutoff.denominator - 1) // cutoff.denominator, q)
    elif cutoff > q:
        candidates = range(q + 1, cutoff.numerator // cutoff.denominator + 1)
    else:
        return
    if any(_is_prime_power(value) for value in candidates):
        raise WeilTransitionPlanError(
            "stencil neighbor reaches or crosses another prime-power transition"
        )


def _policy_from_record(raw: Mapping[str, Any]) -> v1.WeilSearchPolicy:
    if not isinstance(raw, Mapping):
        raise WeilTransitionPlanError("cell policy must be an object")
    witness = raw.get("witness_policy")
    if not isinstance(witness, Mapping):
        raise WeilTransitionPlanError("cell witness policy must be an object")
    try:
        policy = v1.WeilSearchPolicy(
            attempt_bits=tuple(raw.get("attempt_bits", ())),
            confirmation_bits=raw.get("confirmation_bits"),
            eigenpair_count=witness.get("eigenpair_count"),
            scale_bits=tuple(witness.get("scale_bits", ())),
        )
    except (TypeError, v1.WeilSearchPlanError) as exc:
        raise WeilTransitionPlanError("invalid per-cell v1 search policy") from exc
    if dict(raw) != policy.to_record():
        raise WeilTransitionPlanError("cell policy is not canonical v1 policy")
    return policy


@dataclass(frozen=True)
class WeilTransitionCellContract:
    """One exact stencil point, matrix input, and v1 numerical policy."""

    q: int
    prime: int
    exponent: int
    position: str
    stencil_exponent: int | None
    cutoff_numerator: int
    cutoff_denominator: int
    degree: int
    policy: v1.WeilSearchPolicy

    def __post_init__(self) -> None:
        if not isinstance(self.policy, v1.WeilSearchPolicy):
            raise WeilTransitionPlanError("policy must be a WeilSearchPolicy")
        q = _parse_integer(self.q, "q", 2)
        prime = _parse_integer(self.prime, "prime", 2)
        exponent = _parse_integer(self.exponent, "prime-power exponent", 1)
        degree = _parse_integer(self.degree, "degree", 0)
        if not _is_prime(prime) or prime**exponent != q:
            raise WeilTransitionPlanError("transition must satisfy q=p^r with p prime")
        j = self.stencil_exponent
        if j is not None:
            j = _parse_integer(j, "stencil exponent j", 1)
        expected = _expected_cutoff(q, self.position, j)
        supplied = Fraction(
            _parse_integer(self.cutoff_numerator, "cutoff numerator", 1),
            _parse_integer(self.cutoff_denominator, "cutoff denominator", 1),
        )
        if supplied != expected:
            raise WeilTransitionPlanError("cutoff does not match exact stencil formula")
        if supplied <= 1:
            raise WeilTransitionPlanError("cutoff must be strictly greater than one")
        _validate_no_crossing(q, supplied)
        object.__setattr__(self, "q", q)
        object.__setattr__(self, "prime", prime)
        object.__setattr__(self, "exponent", exponent)
        object.__setattr__(self, "stencil_exponent", j)
        object.__setattr__(self, "cutoff_numerator", supplied.numerator)
        object.__setattr__(self, "cutoff_denominator", supplied.denominator)
        object.__setattr__(self, "degree", degree)

    @property
    def cutoff(self) -> Fraction:
        return Fraction(self.cutoff_numerator, self.cutoff_denominator)

    @property
    def sort_key(self) -> tuple[Fraction, int, int, int, int, int, int]:
        return (
            self.cutoff,
            self.degree,
            self.cutoff_numerator,
            self.cutoff_denominator,
            self.q,
            -1 if self.stencil_exponent is None else self.stencil_exponent,
            {"below": 0, "exact": 1, "above": 2}[self.position],
        )

    @property
    def contract_id(self) -> str:
        location = (
            "exact"
            if self.position == "exact"
            else f"j-{self.stencil_exponent}-{self.position}"
        )
        return f"q-{self.q}-{location}-n-{self.degree}"

    @property
    def v1_cell(self) -> v1.WeilSearchCell:
        return v1.WeilSearchCell(
            self.cutoff_numerator, self.cutoff_denominator, self.degree
        )

    @property
    def matrix_input_sha256(self) -> str:
        return self.v1_cell.to_record()["input_sha256"]

    def to_record(self) -> dict[str, Any]:
        transition = {
            "q": str(self.q),
            "prime": str(self.prime),
            "exponent": str(self.exponent),
            "position": self.position,
            "stencil_exponent": (
                None if self.stencil_exponent is None else str(self.stencil_exponent)
            ),
            "formula": _formula(self.position),
        }
        matrix_input = self.v1_cell.to_record()
        body = {
            "schema": CONTRACT_SCHEMA,
            "contract_id": self.contract_id,
            "transition": transition,
            "matrix_input": matrix_input,
            "matrix_input_sha256": self.matrix_input_sha256,
            "policy": self.policy.to_record(),
            "numerical_engine_id": v1.ENGINE_ID,
            "v1_plan_binding": "cell_contract_sha256",
        }
        return {**body, "cell_contract_sha256": content_sha256(body)}


def _contract_from_record(raw: Mapping[str, Any]) -> WeilTransitionCellContract:
    if not isinstance(raw, Mapping):
        raise WeilTransitionPlanError("plan cell must be an object")
    if raw.get("schema") != CONTRACT_SCHEMA:
        raise WeilTransitionPlanError("unexpected cell contract schema")
    transition = raw.get("transition")
    matrix_input = raw.get("matrix_input")
    if not isinstance(transition, Mapping) or not isinstance(matrix_input, Mapping):
        raise WeilTransitionPlanError("cell transition and matrix input are required")
    cutoff = matrix_input.get("cutoff_c")
    if not isinstance(cutoff, Mapping):
        raise WeilTransitionPlanError("matrix cutoff must be an object")
    position = transition.get("position")
    j_raw = transition.get("stencil_exponent")
    j = None if j_raw is None else _parse_integer(j_raw, "stencil exponent j", 1)
    policy = _policy_from_record(raw.get("policy", {}))
    contract = WeilTransitionCellContract(
        q=_parse_integer(transition.get("q"), "q", 2),
        prime=_parse_integer(transition.get("prime"), "prime", 2),
        exponent=_parse_integer(transition.get("exponent"), "prime-power exponent", 1),
        position=str(position),
        stencil_exponent=j,
        cutoff_numerator=_parse_integer(cutoff.get("numerator"), "cutoff numerator", 1),
        cutoff_denominator=_parse_integer(cutoff.get("denominator"), "cutoff denominator", 1),
        degree=_parse_integer(matrix_input.get("degree"), "degree", 0),
        policy=policy,
    )
    canonical = contract.to_record()
    if dict(raw) != canonical:
        raise WeilTransitionPlanError(
            "cell contract contains noncanonical, missing, or unknown fields"
        )
    return contract


def _validate_stencils(contracts: Sequence[WeilTransitionCellContract]) -> None:
    groups: dict[tuple[int, int], list[WeilTransitionCellContract]] = {}
    matrix_hashes: set[str] = set()
    contract_ids: set[str] = set()
    for contract in contracts:
        if contract.matrix_input_sha256 in matrix_hashes:
            raise WeilTransitionPlanError("duplicate numerical matrix input")
        if contract.contract_id in contract_ids:
            raise WeilTransitionPlanError("duplicate cell contract id")
        matrix_hashes.add(contract.matrix_input_sha256)
        contract_ids.add(contract.contract_id)
        groups.setdefault((contract.q, contract.degree), []).append(contract)
    for key, group in groups.items():
        exact = [item for item in group if item.position == "exact"]
        below = {
            item.stencil_exponent for item in group if item.position == "below"
        }
        above = {
            item.stencil_exponent for item in group if item.position == "above"
        }
        if len(exact) != 1 or not below or below != above:
            raise WeilTransitionPlanError(
                f"incomplete exact below/exact/above stencil for q,N={key}"
            )
        if len(group) != 1 + 2 * len(below):
            raise WeilTransitionPlanError("duplicate stencil position")


def _validate_frozen_q7_q9_batch(
    contracts: Sequence[WeilTransitionCellContract],
) -> None:
    expected: set[tuple[int, int, int, int, str, int | None]] = set()
    for q, prime, exponent, degrees in (
        (7, 7, 1, (12, 16, 20)),
        (8, 2, 3, (16, 20, 24)),
        (9, 3, 2, (16, 20, 24)),
    ):
        for degree in degrees:
            expected.add((q, prime, exponent, degree, "exact", None))
            for j in (4, 6, 8, 10):
                expected.add((q, prime, exponent, degree, "below", j))
                expected.add((q, prime, exponent, degree, "above", j))
    actual = {
        (
            item.q,
            item.prime,
            item.exponent,
            item.degree,
            item.position,
            item.stencil_exponent,
        )
        for item in contracts
    }
    if len(contracts) != 81 or actual != expected:
        raise WeilTransitionPlanError(
            "frozen q7-q9 batch axes or 81-cell cardinality changed"
        )
    schedules = {
        (7, 12): ((192, 384), 768),
        (7, 16): ((192, 384, 768), 1536),
        (7, 20): ((384, 768), 1536),
        (8, 16): ((192, 384, 768), 1536),
        (8, 20): ((384, 768), 1536),
        (8, 24): ((384, 768), 1536),
        (9, 16): ((192, 384, 768), 1536),
        (9, 20): ((384, 768), 1536),
        (9, 24): ((384, 768), 1536),
    }
    for item in contracts:
        attempt_bits, confirmation_bits = schedules[(item.q, item.degree)]
        if (
            item.policy.attempt_bits != attempt_bits
            or item.policy.confirmation_bits != confirmation_bits
            or item.policy.eigenpair_count != 3
            or item.policy.scale_bits != (8, 12, 16, 24, 32, 48)
        ):
            raise WeilTransitionPlanError(
                "frozen q7-q9 per-cell precision or witness policy changed"
            )


def canonicalize_transition_plan(raw: Mapping[str, Any]) -> dict[str, Any]:
    """Validate and canonically order an explicit v2 transition plan."""

    if not isinstance(raw, Mapping) or raw.get("schema") != PLAN_SCHEMA:
        raise WeilTransitionPlanError("unexpected transition plan schema")
    if raw.get("classification", "EXPLORATORY") != "EXPLORATORY":
        raise WeilTransitionPlanError("transition plans must remain EXPLORATORY")
    if raw.get("hypothesis_status", "UNRESOLVED") != "UNRESOLVED":
        raise WeilTransitionPlanError("transition plans cannot resolve RH")
    allowed_fields = {
        "schema",
        "classification",
        "hypothesis_status",
        "cells",
        "engine",
        "backend_contract",
        "frozen_batch",
        "limitation",
        "plan_sha256",
    }
    if set(raw) - allowed_fields:
        raise WeilTransitionPlanError("transition plan has unknown fields")
    raw_cells = raw.get("cells")
    if not isinstance(raw_cells, list) or not raw_cells:
        raise WeilTransitionPlanError("cells must be a nonempty list")
    contracts = sorted(
        (_contract_from_record(cell) for cell in raw_cells),
        key=lambda contract: contract.sort_key,
    )
    _validate_stencils(contracts)
    frozen_batch = raw.get("frozen_batch")
    if frozen_batch is not None:
        if not isinstance(frozen_batch, Mapping) or dict(frozen_batch) != FROZEN_Q7_Q9_BATCH:
            raise WeilTransitionPlanError("frozen q7-q9 batch metadata changed")
        _validate_frozen_q7_q9_batch(contracts)
    backend = raw.get("backend_contract", v1._backend_record())
    if not isinstance(backend, Mapping):
        raise WeilTransitionPlanError("backend contract must be an object")
    backend_record = {
        "python_flint": str(backend.get("python_flint", "")),
        "flint": str(backend.get("flint", "")),
    }
    if set(backend) != {"python_flint", "flint"}:
        raise WeilTransitionPlanError("backend contract has unknown fields")
    if not all(backend_record.values()):
        raise WeilTransitionPlanError("backend versions are required")
    engine_record = {
        "algorithm_id": ENGINE_ID,
        "numerical_engine_id": v1.ENGINE_ID,
        "matrix": "A=P-R-S",
        "cell_order": (
            "exact rational cutoff ascending globally, then degree, "
            "numerator, denominator"
        ),
        "candidate_behavior": "finish bounded plan; never stop early",
        "parallelism": "serial; FLINT ctx.prec is shared mutable state",
    }
    supplied_engine = raw.get("engine")
    if supplied_engine is not None and (
        not isinstance(supplied_engine, Mapping)
        or dict(supplied_engine) != engine_record
    ):
        raise WeilTransitionPlanError("transition engine contract changed")
    if raw.get("limitation", EXPLORATORY_LIMITATION) != EXPLORATORY_LIMITATION:
        raise WeilTransitionPlanError("transition limitation changed")
    body = {
        "schema": PLAN_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "cells": [contract.to_record() for contract in contracts],
        "engine": engine_record,
        "backend_contract": backend_record,
        "limitation": EXPLORATORY_LIMITATION,
    }
    if frozen_batch is not None:
        body["frozen_batch"] = FROZEN_Q7_Q9_BATCH
    plan_sha256 = content_sha256(body)
    if raw.get("plan_sha256") not in {None, plan_sha256}:
        raise WeilTransitionPlanError("transition plan hash mismatch")
    return {**body, "plan_sha256": plan_sha256}


def load_transition_plan(path: Path) -> dict[str, Any]:
    """Load and canonicalize a v2 transition plan."""

    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise WeilTransitionPlanError(f"cannot read transition plan: {path}") from exc
    return canonicalize_transition_plan(raw)


def _contracts_from_plan(plan: Mapping[str, Any]) -> list[WeilTransitionCellContract]:
    return [_contract_from_record(record) for record in plan["cells"]]


def _with_payload_hash(payload: Mapping[str, Any]) -> dict[str, Any]:
    body = dict(payload)
    body.pop("payload_sha256", None)
    return {**body, "payload_sha256": content_sha256(body)}


def _require_payload_hash(
    artifact: Mapping[str, Any], error_type: type[WeilTransitionError]
) -> None:
    supplied = artifact.get("payload_sha256")
    body = {key: value for key, value in artifact.items() if key != "payload_sha256"}
    if supplied != content_sha256(body):
        raise error_type("payload hash mismatch")


def _runtime_backend(plan: Mapping[str, Any]) -> None:
    if plan["backend_contract"] != v1._backend_record():
        raise WeilTransitionPlanError(
            "installed python-flint/FLINT versions do not match the plan"
        )


def _run_manifest(plan: Mapping[str, Any]) -> dict[str, Any]:
    return _with_payload_hash(
        {
            "schema": RUN_SCHEMA,
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "plan": {"path": "plan.json", "plan_sha256": plan["plan_sha256"]},
            "backend": v1._backend_record(),
            "engine_id": ENGINE_ID,
            "numerical_engine_id": v1.ENGINE_ID,
            "candidate_behavior": "finish bounded plan; never stop early",
            "limitation": EXPLORATORY_LIMITATION,
        }
    )


def _load_json(path: Path, error_type: type[WeilTransitionError]) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise error_type(f"cannot read checkpoint: {path}") from exc
    if not isinstance(value, dict):
        raise error_type("checkpoint must contain an object")
    return value


def _safe_artifact_path(root: Path, relative: Any) -> Path:
    if not isinstance(relative, str) or not relative:
        raise WeilTransitionVerificationError("artifact path is invalid")
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise WeilTransitionVerificationError("artifact path escapes checkpoint")
    return path


def _discovery_attempt(artifact: Mapping[str, Any]) -> Mapping[str, Any]:
    discoveries = [
        attempt
        for attempt in artifact.get("attempts", [])
        if attempt.get("decision") == "NEGATIVE_WITNESS_DISCOVERED"
    ]
    if len(discoveries) != 1:
        raise WeilTransitionCheckpointError(
            "v1 candidate must have exactly one negative discovery"
        )
    return discoveries[0]


def _validate_transition_evidence(
    contract: WeilTransitionCellContract,
    evidence: Mapping[str, Any],
    precision_bits: int,
    error_type: type[WeilTransitionError],
) -> None:
    """Bind target p,r,q metadata to the regenerated prime transcript."""

    transcript = evidence.get("prime_power_transcript")
    if not isinstance(transcript, list):
        raise error_type("prime-power transcript is missing")
    target = [
        record
        for record in transcript
        if isinstance(record, Mapping) and record.get("power") == str(contract.q)
    ]
    expected_count = 0 if contract.position == "below" else 1
    if len(target) != expected_count:
        raise error_type("target prime-power endpoint presence changed")
    if not target:
        return
    record = target[0]
    if record.get("prime") != str(contract.prime) or record.get(
        "exponent"
    ) != str(contract.exponent):
        raise error_type("target q is not represented by its bound p^r tuple")
    try:
        with v1._working_precision(precision_bits):
            expected_weight = (
                v1.arb(contract.prime).log() / v1.arb(contract.q).sqrt()
            )
            expected_location = v1.arb(contract.q).log()
            recorded_weight = v1.arb_from_dyadic(record["weight"]["dyadic"])
            recorded_location = v1.arb_from_dyadic(
                record["log_location"]["dyadic"]
            )
            valid = recorded_weight.overlaps(
                expected_weight
            ) and recorded_location.overlaps(expected_location)
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise error_type("target prime-power transcript is malformed") from exc
    if not valid:
        raise error_type(
            "target transcript weight must be log(p)/sqrt(q), at log(q)"
        )


def _dedicated_confirmation(
    contract: WeilTransitionCellContract,
    v1_artifact: Mapping[str, Any],
) -> dict[str, Any]:
    """Confirm a v1 candidate at the policy's dedicated v2 precision.

    This artifact is deliberately outside the v1 attempt list, whose state
    machine is terminal after its own confirmation attempt.
    """

    if v1_artifact.get("status") != "NEGATIVE_CANDIDATE_QUARANTINED":
        raise WeilTransitionCheckpointError(
            "dedicated confirmation requires a v1 quarantined candidate"
        )
    discovery = _discovery_attempt(v1_artifact)
    candidate = v1_artifact.get("quarantined_candidate")
    if not isinstance(candidate, Mapping) or not isinstance(
        candidate.get("witness"), list
    ):
        raise WeilTransitionCheckpointError("v1 candidate witness is missing")
    witness = [int(value) for value in candidate["witness"]]
    precision_bits = contract.policy.confirmation_bits
    base = {
        "schema": "rh-lab/weil-transition-dedicated-confirmation/v2",
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "contract_id": contract.contract_id,
        "cell_contract_sha256": contract.to_record()["cell_contract_sha256"],
        "matrix_input_sha256": contract.matrix_input_sha256,
        "discovery_attempt_payload_sha256": discovery["payload_sha256"],
        "witness": [str(value) for value in witness],
        "precision_bits": str(precision_bits),
        "limitation": EXPLORATORY_LIMITATION,
    }
    terminal = v1_artifact["attempts"][-1]
    if (
        terminal.get("kind") == "CONFIRMATION"
        and terminal.get("precision_bits") == str(precision_bits)
        and terminal.get("decision") == "NEGATIVE_WITNESS_CONFIRMED"
    ):
        return _with_payload_hash(
            {
                **base,
                "mode": "REUSED_V1_TERMINAL_CONFIRMATION",
                "source_attempt_payload_sha256": terminal["payload_sha256"],
                "witness_evaluation": terminal["witness_evaluation"],
                "precision_containment": terminal["precision_containment"],
                "decision": "V2_NEGATIVE_CANDIDATE_CONFIRMED",
            }
        )

    with v1._working_precision(precision_bits):
        try:
            components = v1._build_components(contract.v1_cell)
        except ArithmeticError:
            return _with_payload_hash(
                {
                    **base,
                    "mode": "DEDICATED_REGENERATION",
                    "matrix_construction": {
                        "classification": "INCONCLUSIVE",
                        "failure_code": "SOURCE_ORACLE_OR_QUADRATURE_INCONCLUSIVE",
                    },
                    "decision": "V2_DEDICATED_CONFIRMATION_INCONCLUSIVE",
                }
            )
        evidence = v1._matrix_evidence(components)
        _validate_transition_evidence(
            contract,
            evidence,
            precision_bits,
            WeilTransitionCheckpointError,
        )
        evaluation = v1._stable_evaluation(
            v1.evaluate_integer_rayleigh(components["q"], witness), precision_bits
        )
        contained, failed_entry = v1._precision_containment(
            discovery["matrix_evidence"], evidence
        )
    decision = (
        "V2_NEGATIVE_CANDIDATE_CONFIRMED"
        if evaluation["classification"] == "NEGATIVE" and contained
        else "V2_DEDICATED_CONFIRMATION_INCONCLUSIVE"
    )
    return _with_payload_hash(
        {
            **base,
            "mode": "DEDICATED_REGENERATION",
            "matrix_construction": {"classification": "COMPLETE"},
            "matrix_evidence": evidence,
            "witness_evaluation": evaluation,
            "precision_containment": {
                "all_higher_precision_components_contained": contained,
                "first_failure": failed_entry,
            },
            "decision": decision,
        }
    )


def _derive_v1_artifact(
    contract: WeilTransitionCellContract,
    artifact: Mapping[str, Any],
    *,
    replay: bool,
) -> dict[str, Any]:
    try:
        v1._require_payload_hash(artifact, v1.WeilSearchCheckpointError)
        if artifact.get("schema") != v1.CELL_SCHEMA:
            raise v1.WeilSearchCheckpointError("unexpected v1 cell schema")
        if artifact.get("classification") != "EXPLORATORY":
            raise v1.WeilSearchCheckpointError("v1 cell was improperly promoted")
        if artifact.get("hypothesis_status") != "UNRESOLVED":
            raise v1.WeilSearchCheckpointError("v1 cell changed RH status")
        if artifact.get("cell") != contract.v1_cell.to_record():
            raise v1.WeilSearchCheckpointError("v1 matrix input changed")
        if artifact.get("policy") != contract.policy.to_record():
            raise v1.WeilSearchCheckpointError("v1 cell policy changed")
        contract_hash = contract.to_record()["cell_contract_sha256"]
        if artifact.get("plan_sha256") != contract_hash:
            raise v1.WeilSearchCheckpointError("v1 cell contract binding changed")
        attempts = artifact.get("attempts")
        if not isinstance(attempts, list) or not attempts:
            raise v1.WeilSearchCheckpointError("v1 cell attempts are missing")
        v1._validate_attempt_hashes(attempts)
        if replay:
            for attempt in attempts:
                v1._replay_attempt(contract.v1_cell, attempt)
                if attempt.get("matrix_construction", {}).get(
                    "classification"
                ) == "COMPLETE":
                    _validate_transition_evidence(
                        contract,
                        attempt["matrix_evidence"],
                        int(attempt["precision_bits"]),
                        WeilTransitionVerificationError,
                    )
        derived = v1._search_weil_cell(
            contract.v1_cell,
            contract.policy,
            existing_attempts=attempts,
            plan_sha256=contract_hash,
        )
    except (
        AttributeError,
        IndexError,
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        error_type = (
            WeilTransitionVerificationError
            if replay
            else WeilTransitionCheckpointError
        )
        raise error_type("v1 cell artifact failed derivation") from exc
    if dict(artifact) != derived:
        error_type = (
            WeilTransitionVerificationError
            if replay
            else WeilTransitionCheckpointError
        )
        raise error_type("v1 cell status does not follow from its attempts")
    return derived


def _v2_status(
    artifact: Mapping[str, Any], confirmation: Mapping[str, Any] | None
) -> str:
    status = artifact["status"]
    if status != "NEGATIVE_CANDIDATE_QUARANTINED":
        if confirmation is not None:
            raise WeilTransitionCheckpointError(
                "noncandidate v1 cell cannot have dedicated confirmation"
            )
        return str(status)
    if confirmation is None:
        raise WeilTransitionCheckpointError(
            "v1 negative candidate lacks dedicated v2 confirmation"
        )
    if confirmation.get("decision") == "V2_NEGATIVE_CANDIDATE_CONFIRMED":
        return "NEGATIVE_CANDIDATE_QUARANTINED"
    if confirmation.get("decision") == "V2_DEDICATED_CONFIRMATION_INCONCLUSIVE":
        return "NEGATIVE_OBSERVATION_QUARANTINED_INCONCLUSIVE"
    raise WeilTransitionCheckpointError("unknown dedicated confirmation decision")


def _confirmation_for_cell(
    contract: WeilTransitionCellContract,
    artifact: Mapping[str, Any],
    path: Path,
) -> dict[str, Any] | None:
    if artifact["status"] != "NEGATIVE_CANDIDATE_QUARANTINED":
        if path.exists():
            raise WeilTransitionCheckpointError(
                "unexpected dedicated confirmation for noncandidate cell"
            )
        return None
    expected = _dedicated_confirmation(contract, artifact)
    if path.exists():
        stored = _load_json(path, WeilTransitionCheckpointError)
        _require_payload_hash(stored, WeilTransitionCheckpointError)
        if stored != expected:
            raise WeilTransitionCheckpointError(
                "dedicated confirmation no longer regenerates"
            )
        return stored
    v1._atomic_write_json(path, expected)
    return expected


def _index_payload(
    plan: Mapping[str, Any],
    contracts: Sequence[WeilTransitionCellContract],
    completed: Mapping[
        str, tuple[Mapping[str, Any], Mapping[str, Any] | None]
    ],
) -> dict[str, Any]:
    references: list[dict[str, Any]] = []
    status_values: list[str] = []
    for contract in contracts:
        pair = completed.get(contract.contract_id)
        if pair is None:
            continue
        artifact, confirmation = pair
        status = _v2_status(artifact, confirmation)
        status_values.append(status)
        reference: dict[str, Any] = {
            "contract_id": contract.contract_id,
            "cell_contract_sha256": contract.to_record()[
                "cell_contract_sha256"
            ],
            "matrix_input_sha256": contract.matrix_input_sha256,
            "path": f"cells/{contract.contract_id}.json",
            "v1_cell_payload_sha256": artifact["payload_sha256"],
            "v1_status": artifact["status"],
            "status": status,
            "dedicated_confirmation": None,
        }
        if confirmation is not None:
            reference["dedicated_confirmation"] = {
                "path": f"confirmations/{contract.contract_id}.json",
                "payload_sha256": confirmation["payload_sha256"],
                "decision": confirmation["decision"],
            }
        references.append(reference)
    missing = [
        contract.contract_id
        for contract in contracts
        if contract.contract_id not in completed
    ]
    statuses = {
        status: str(status_values.count(status))
        for status in (
            "FINITE_POSITIVE_CERTIFIED",
            "NEGATIVE_CANDIDATE_QUARANTINED",
            "NEGATIVE_OBSERVATION_QUARANTINED_INCONCLUSIVE",
            "INCONCLUSIVE_MAX_PRECISION",
        )
    }
    if statuses["NEGATIVE_CANDIDATE_QUARANTINED"] != "0":
        conclusion = "QUARANTINED_NEGATIVE_CANDIDATE_IN_FINITE_BATCH"
    elif statuses["NEGATIVE_OBSERVATION_QUARANTINED_INCONCLUSIVE"] != "0":
        conclusion = (
            "QUARANTINED_CONFLICTING_NEGATIVE_OBSERVATION_IN_FINITE_BATCH"
        )
    elif missing:
        conclusion = "INCOMPLETE_BOUNDED_TRANSITION_SEARCH"
    else:
        conclusion = (
            "NO_CERTIFIED_NEGATIVE_WITNESS_FOUND_IN_FINITE_TRANSITION_BATCH"
        )
    return _with_payload_hash(
        {
            "schema": INDEX_SCHEMA,
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "conclusion": conclusion,
            "plan": {"path": "plan.json", "plan_sha256": plan["plan_sha256"]},
            "run": {"path": "run.json"},
            "backend": v1._backend_record(),
            "cells": references,
            "missing_contract_ids": missing,
            "progress": {
                "planned_cells": str(len(contracts)),
                "completed_cells": str(len(completed)),
                "status_counts": statuses,
            },
            "candidate_behavior": "finish bounded plan; never stop early",
            "limitation": EXPLORATORY_LIMITATION,
        }
    )


def run_weil_transition_search(
    plan: Mapping[str, Any],
    checkpoint_dir: Path,
    *,
    resume: bool = False,
    max_cells: int | None = None,
) -> dict[str, Any]:
    """Run or resume a deterministic explicit transition batch.

    A quarantined candidate does not terminate the batch. ``max_cells`` only
    bounds newly completed v1 numerical cells in this invocation.
    """

    canonical_plan = canonicalize_transition_plan(plan)
    _runtime_backend(canonical_plan)
    if max_cells is not None:
        if isinstance(max_cells, bool) or not isinstance(max_cells, int):
            raise TypeError("max_cells must be an integer or None")
        if max_cells < 0:
            raise ValueError("max_cells must be nonnegative")
    root = checkpoint_dir.resolve()
    plan_path = root / "plan.json"
    run_path = root / "run.json"
    expected_run = _run_manifest(canonical_plan)
    if root.exists() and any(root.iterdir()):
        if not resume:
            raise WeilTransitionCheckpointError(
                "checkpoint already exists; pass resume=True"
            )
        try:
            stored_plan_raw = _load_json(
                plan_path, WeilTransitionCheckpointError
            )
            stored_plan = canonicalize_transition_plan(stored_plan_raw)
        except WeilTransitionPlanError as exc:
            raise WeilTransitionCheckpointError(
                "stored transition plan is invalid"
            ) from exc
        if stored_plan_raw != stored_plan:
            raise WeilTransitionCheckpointError(
                "stored transition plan is not canonical"
            )
        if stored_plan != canonical_plan:
            raise WeilTransitionCheckpointError("checkpoint plan changed")
        stored_run = _load_json(run_path, WeilTransitionCheckpointError)
        _require_payload_hash(stored_run, WeilTransitionCheckpointError)
        if stored_run != expected_run:
            raise WeilTransitionCheckpointError(
                "checkpoint plan, backend, or engine does not match"
            )
    else:
        if resume:
            raise WeilTransitionCheckpointError(
                "resume requested without a checkpoint"
            )
        root.mkdir(parents=True, exist_ok=True)
        v1._atomic_write_json(plan_path, canonical_plan)
        v1._atomic_write_json(run_path, expected_run)

    contracts = _contracts_from_plan(canonical_plan)
    completed: dict[
        str, tuple[dict[str, Any], dict[str, Any] | None]
    ] = {}
    for contract in contracts:
        cell_path = root / "cells" / f"{contract.contract_id}.json"
        if not cell_path.exists():
            continue
        artifact = _load_json(cell_path, WeilTransitionCheckpointError)
        artifact = _derive_v1_artifact(contract, artifact, replay=False)
        confirmation = _confirmation_for_cell(
            contract,
            artifact,
            root / "confirmations" / f"{contract.contract_id}.json",
        )
        completed[contract.contract_id] = (artifact, confirmation)

    newly_completed = 0
    for contract in contracts:
        if contract.contract_id in completed:
            continue
        if max_cells is not None and newly_completed >= max_cells:
            break
        attempt_directory = root / "attempts" / contract.contract_id
        try:
            attempts = v1._load_attempts(attempt_directory)
        except v1.WeilSearchCheckpointError as exc:
            raise WeilTransitionCheckpointError(
                "v1 attempt checkpoint is invalid"
            ) from exc

        def save_attempt(attempt: dict[str, Any]) -> None:
            sequence = int(attempt["sequence"])
            v1._atomic_write_json(
                attempt_directory / f"{sequence:04d}.json", attempt
            )

        try:
            artifact = v1._search_weil_cell(
                contract.v1_cell,
                contract.policy,
                existing_attempts=attempts,
                attempt_sink=save_attempt,
                plan_sha256=contract.to_record()["cell_contract_sha256"],
            )
        except v1.WeilSearchError as exc:
            raise WeilTransitionCheckpointError(
                "v1 numerical cell search failed"
            ) from exc
        cell_path = root / "cells" / f"{contract.contract_id}.json"
        v1._atomic_write_json(cell_path, artifact)
        confirmation = _confirmation_for_cell(
            contract,
            artifact,
            root / "confirmations" / f"{contract.contract_id}.json",
        )
        completed[contract.contract_id] = (artifact, confirmation)
        newly_completed += 1

    index = _index_payload(canonical_plan, contracts, completed)
    v1._atomic_write_json(root / "index.json", index)
    return index


def verify_weil_transition_search(
    index: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
) -> dict[str, Any]:
    """Regenerate evidence and semantically derive every referenced v1 cell."""

    root = checkpoint_dir.resolve()
    if isinstance(index, Path):
        index_record = _load_json(index, WeilTransitionVerificationError)
    elif isinstance(index, Mapping):
        index_record = dict(index)
    else:
        raise TypeError("index must be a mapping or Path")
    _require_payload_hash(index_record, WeilTransitionVerificationError)
    if index_record.get("schema") != INDEX_SCHEMA:
        raise WeilTransitionVerificationError("unexpected transition index schema")
    if index_record.get("classification") != "EXPLORATORY":
        raise WeilTransitionVerificationError("transition index was promoted")
    if index_record.get("hypothesis_status") != "UNRESOLVED":
        raise WeilTransitionVerificationError("transition index changed RH status")
    plan_reference = index_record.get("plan")
    if not isinstance(plan_reference, Mapping):
        raise WeilTransitionVerificationError("index plan reference is missing")
    plan_path = _safe_artifact_path(root, plan_reference.get("path"))
    try:
        stored_plan_raw = _load_json(plan_path, WeilTransitionVerificationError)
        plan = canonicalize_transition_plan(stored_plan_raw)
    except WeilTransitionPlanError as exc:
        raise WeilTransitionVerificationError("stored plan is invalid") from exc
    if plan != stored_plan_raw:
        raise WeilTransitionVerificationError("stored plan is not canonical")
    if plan_reference.get("plan_sha256") != plan["plan_sha256"]:
        raise WeilTransitionVerificationError("index plan hash changed")
    try:
        _runtime_backend(plan)
    except WeilTransitionPlanError as exc:
        raise WeilTransitionVerificationError(str(exc)) from exc
    run_reference = index_record.get("run")
    if not isinstance(run_reference, Mapping):
        raise WeilTransitionVerificationError("index run reference is missing")
    run_path = _safe_artifact_path(root, run_reference.get("path"))
    run = _load_json(run_path, WeilTransitionVerificationError)
    _require_payload_hash(run, WeilTransitionVerificationError)
    if run != _run_manifest(plan):
        raise WeilTransitionVerificationError("run manifest changed")

    contracts = _contracts_from_plan(plan)
    by_id = {contract.contract_id: contract for contract in contracts}
    references = index_record.get("cells")
    if not isinstance(references, list):
        raise WeilTransitionVerificationError("index cells must be a list")
    seen: set[str] = set()
    completed: dict[
        str, tuple[dict[str, Any], dict[str, Any] | None]
    ] = {}
    replayed_attempts = 0
    replayed_integer_evaluations = 0
    for reference in references:
        if not isinstance(reference, Mapping):
            raise WeilTransitionVerificationError("cell reference must be an object")
        contract_id = reference.get("contract_id")
        if contract_id not in by_id or contract_id in seen:
            raise WeilTransitionVerificationError(
                "unknown or duplicate cell contract reference"
            )
        seen.add(str(contract_id))
        contract = by_id[str(contract_id)]
        canonical_contract = contract.to_record()
        if reference.get("cell_contract_sha256") != canonical_contract[
            "cell_contract_sha256"
        ]:
            raise WeilTransitionVerificationError("cell contract hash changed")
        if reference.get("matrix_input_sha256") != contract.matrix_input_sha256:
            raise WeilTransitionVerificationError("matrix input hash changed")
        cell_path = _safe_artifact_path(root, reference.get("path"))
        artifact = _load_json(cell_path, WeilTransitionVerificationError)
        artifact = _derive_v1_artifact(contract, artifact, replay=True)
        replayed_attempts += len(artifact["attempts"])
        for attempt in artifact["attempts"]:
            if attempt.get("kind") == "SEARCH":
                replayed_integer_evaluations += len(
                    attempt.get("candidate_generation", {}).get("tested", [])
                )
            elif attempt.get("witness_evaluation") is not None:
                replayed_integer_evaluations += 1
        if reference.get("v1_cell_payload_sha256") != artifact["payload_sha256"]:
            raise WeilTransitionVerificationError("v1 cell reference hash changed")
        if reference.get("v1_status") != artifact["status"]:
            raise WeilTransitionVerificationError("v1 cell status changed")

        confirmation_reference = reference.get("dedicated_confirmation")
        confirmation: dict[str, Any] | None = None
        if artifact["status"] == "NEGATIVE_CANDIDATE_QUARANTINED":
            if not isinstance(confirmation_reference, Mapping):
                raise WeilTransitionVerificationError(
                    "candidate lacks dedicated confirmation reference"
                )
            confirmation_path = _safe_artifact_path(
                root, confirmation_reference.get("path")
            )
            confirmation = _load_json(
                confirmation_path, WeilTransitionVerificationError
            )
            _require_payload_hash(confirmation, WeilTransitionVerificationError)
            expected_confirmation = _dedicated_confirmation(contract, artifact)
            if confirmation != expected_confirmation:
                raise WeilTransitionVerificationError(
                    "dedicated confirmation evidence changed on replay"
                )
            if confirmation_reference.get("payload_sha256") != confirmation[
                "payload_sha256"
            ]:
                raise WeilTransitionVerificationError(
                    "dedicated confirmation reference hash changed"
                )
            if confirmation_reference.get("decision") != confirmation["decision"]:
                raise WeilTransitionVerificationError(
                    "dedicated confirmation decision changed"
                )
            if confirmation.get("mode") == "DEDICATED_REGENERATION" and confirmation.get(
                "witness_evaluation"
            ) is not None:
                replayed_integer_evaluations += 1
        elif confirmation_reference is not None:
            raise WeilTransitionVerificationError(
                "noncandidate has dedicated confirmation reference"
            )
        try:
            status = _v2_status(artifact, confirmation)
        except WeilTransitionCheckpointError as exc:
            raise WeilTransitionVerificationError(str(exc)) from exc
        if reference.get("status") != status:
            raise WeilTransitionVerificationError("v2 cell status changed")
        completed[contract.contract_id] = (artifact, confirmation)

    expected_index = _index_payload(plan, contracts, completed)
    if index_record != expected_index:
        raise WeilTransitionVerificationError(
            "transition index does not match referenced artifacts"
        )
    return {
        "classification": "REPRODUCED_EXPLORATORY_TRANSITION_SEARCH",
        "hypothesis_status": "UNRESOLVED",
        "conclusion": index_record["conclusion"],
        "verified_cells": str(len(completed)),
        "replayed_v1_attempts": str(replayed_attempts),
        "replayed_integer_witness_evaluations": str(
            replayed_integer_evaluations
        ),
        "same_backend_replay_only": True,
    }


__all__ = [
    "CONTRACT_SCHEMA",
    "ENGINE_ID",
    "FROZEN_Q7_Q9_BATCH",
    "INDEX_SCHEMA",
    "PLAN_SCHEMA",
    "RUN_SCHEMA",
    "WeilTransitionCellContract",
    "WeilTransitionCheckpointError",
    "WeilTransitionError",
    "WeilTransitionPlanError",
    "WeilTransitionVerificationError",
    "canonicalize_transition_plan",
    "load_transition_plan",
    "run_weil_transition_search",
    "verify_weil_transition_search",
]
