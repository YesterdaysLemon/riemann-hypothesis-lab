"""Certificate-first finite natural-dilate Nyman distance searches.

This module is deliberately an orchestration layer over :mod:`riemann_lab.nyman`.
Approximate linear algebra may propose a coefficient vector, but every stored
coefficient and bound is an exact dyadic and the only terminal gate consists
of the two core interval certificates evaluated on those stored values.

The finite experiment is exploratory.  No finite list of natural-dilate
distances proves or disproves the Riemann Hypothesis.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from fractions import Fraction
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Mapping, Sequence

import flint
from flint import arb_mat, ctx

from .artifacts import content_sha256
from .balls import arb_to_dyadic
from . import nyman as core


PLAN_SCHEMA = "rh-lab/nyman-natural-plan/v1"
CELL_CONTRACT_SCHEMA = "rh-lab/nyman-natural-cell-contract/v1"
RUN_SCHEMA = "rh-lab/nyman-natural-run/v1"
CELL_SCHEMA = "rh-lab/nyman-natural-cell/v1"
INDEX_SCHEMA = "rh-lab/nyman-natural-index/v1"
KERNEL_SCHEMA = "rh-lab/nyman-natural-shared-kernel/v1"
ENGINE_ID = "natural-dilate-distance-search-v1"

FROZEN_N_VALUES = (8, 16, 32, 64, 128, 256)
FROZEN_GENERATION_BITS = 768
FROZEN_REPLAY_BITS = 1536
FROZEN_COEFFICIENT_GRID_BITS = 256
FROZEN_BOUND_GRID_BITS = 128
FROZEN_TARGET_LOWER_SLACK_BITS = 120

EXPLORATORY_LIMITATION = (
    "These are rigorous bounds for six finite-dimensional natural-dilate "
    "distances only. Finite decay does not prove that the limiting distance "
    "is zero, and a positive finite lower bound does not disprove RH. The "
    "global Riemann Hypothesis status is UNRESOLVED."
)

FROZEN_EXPERIMENT = {
    "experiment_id": "nyman-natural-v1",
    "n_values": [str(value) for value in FROZEN_N_VALUES],
    "cell_count": str(len(FROZEN_N_VALUES)),
    "generation_precision_bits": str(FROZEN_GENERATION_BITS),
    "replay_precision_bits": str(FROZEN_REPLAY_BITS),
    "coefficient_grid": f"2^-{FROZEN_COEFFICIENT_GRID_BITS}",
    "bound_grid": f"2^-{FROZEN_BOUND_GRID_BITS}",
    "target_lower_slack": f"2^-{FROZEN_TARGET_LOWER_SLACK_BITS}",
}


class NymanSearchError(ValueError):
    """Base class for natural-distance search failures."""


class NymanSearchPlanError(NymanSearchError):
    """Raised when a natural-distance search plan is malformed."""


class NymanSearchCheckpointError(NymanSearchError):
    """Raised when checkpoint state is missing, corrupt, or incompatible."""


class NymanSearchVerificationError(NymanSearchError):
    """Raised when stored finite-distance evidence does not reproduce."""


def _parse_integer(
    value: Any,
    name: str,
    minimum: int = 0,
    *,
    error_type: type[NymanSearchError] = NymanSearchPlanError,
) -> int:
    if isinstance(value, bool):
        raise error_type(f"{name} must be an integer")
    if isinstance(value, int):
        parsed = value
    elif isinstance(value, str):
        try:
            parsed = int(value)
        except ValueError as exc:
            raise error_type(f"{name} must be an integer") from exc
        if value != str(parsed):
            raise error_type(f"{name} must use canonical decimal syntax")
    else:
        raise error_type(f"{name} must be an integer")
    if parsed < minimum:
        raise error_type(f"{name} must be at least {minimum}")
    return parsed


def _parse_signed_integer(
    value: Any,
    name: str,
    *,
    error_type: type[NymanSearchError],
) -> int:
    if isinstance(value, bool):
        raise error_type(f"{name} must be an integer")
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        try:
            parsed = int(value)
        except ValueError as exc:
            raise error_type(f"{name} must be an integer") from exc
        if value == str(parsed):
            return parsed
    raise error_type(f"{name} must use canonical signed decimal syntax")


def _backend_record() -> dict[str, str]:
    return {
        "python_flint": str(flint.__version__),
        "flint": str(flint.__FLINT_VERSION__),
    }


def _with_payload_hash(payload: Mapping[str, Any]) -> dict[str, Any]:
    body = dict(payload)
    body.pop("payload_sha256", None)
    return {**body, "payload_sha256": content_sha256(body)}


def _require_payload_hash(
    artifact: Mapping[str, Any], error_type: type[NymanSearchError]
) -> None:
    supplied = artifact.get("payload_sha256")
    body = {key: value for key, value in artifact.items() if key != "payload_sha256"}
    if supplied != content_sha256(body):
        raise error_type("payload hash mismatch")


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _dyadic_record(numerator: int, denominator_exponent: int) -> dict[str, str]:
    return {
        "numerator": str(numerator),
        "denominator_exponent": str(denominator_exponent),
    }


def _parse_dyadic_record(
    raw: Any,
    name: str,
    *,
    expected_exponent: int,
    error_type: type[NymanSearchError],
) -> tuple[int, Fraction]:
    if not isinstance(raw, Mapping) or set(raw) != {
        "numerator",
        "denominator_exponent",
    }:
        raise error_type(f"{name} must be a canonical dyadic object")
    numerator = _parse_signed_integer(
        raw["numerator"], f"{name} numerator", error_type=error_type
    )
    exponent = _parse_integer(
        raw["denominator_exponent"],
        f"{name} denominator exponent",
        0,
        error_type=error_type,
    )
    if exponent != expected_exponent:
        raise error_type(f"{name} is not on the frozen dyadic grid")
    return numerator, Fraction(numerator, 1 << exponent)


def _shared_kernel_contract(max_n: int) -> dict[str, Any]:
    dilates = [str(value) for value in range(1, max_n + 1)]
    body = {
        "max_n": str(max_n),
        "dilates": dilates,
        "dilates_sha256": content_sha256(dilates),
        "normalization": "G_ab=A(a/b)/a; b_a=(log(a)+1-EulerGamma)/a",
        "prefix_rule": (
            "each cell uses the leading principal N-by-N Gram block and "
            "leading N target entries of the shared max-N kernel"
        ),
    }
    return {**body, "kernel_contract_sha256": content_sha256(body)}


@dataclass(frozen=True)
class NymanSearchCellContract:
    """One finite natural-dilate prefix and its shared-kernel binding."""

    n: int
    shared_kernel_n: int

    def __post_init__(self) -> None:
        n = _parse_integer(self.n, "N", 1)
        shared = _parse_integer(self.shared_kernel_n, "shared kernel N", n)
        object.__setattr__(self, "n", n)
        object.__setattr__(self, "shared_kernel_n", shared)

    @property
    def cell_id(self) -> str:
        return f"n-{self.n:04d}"

    def to_record(self) -> dict[str, Any]:
        kernel = _shared_kernel_contract(self.shared_kernel_n)
        body = {
            "schema": CELL_CONTRACT_SCHEMA,
            "cell_id": self.cell_id,
            "n": str(self.n),
            "dimension": str(self.n),
            "dilates": [str(value) for value in range(1, self.n + 1)],
            "shared_kernel_binding": {
                "kernel_contract_sha256": kernel["kernel_contract_sha256"],
                "shared_kernel_n": str(self.shared_kernel_n),
                "prefix_length": str(self.n),
            },
        }
        return {**body, "cell_contract_sha256": content_sha256(body)}


def _contract_from_record(raw: Any) -> NymanSearchCellContract:
    if not isinstance(raw, Mapping):
        raise NymanSearchPlanError("cell contract must be an object")
    if raw.get("schema") != CELL_CONTRACT_SCHEMA:
        raise NymanSearchPlanError("unexpected cell contract schema")
    allowed = {
        "schema",
        "cell_id",
        "n",
        "dimension",
        "dilates",
        "shared_kernel_binding",
        "cell_contract_sha256",
    }
    if set(raw) != allowed:
        raise NymanSearchPlanError("cell contract fields changed")
    n = _parse_integer(raw.get("n"), "N", 1)
    binding = raw.get("shared_kernel_binding")
    if not isinstance(binding, Mapping):
        raise NymanSearchPlanError("shared-kernel binding must be an object")
    shared = _parse_integer(binding.get("shared_kernel_n"), "shared kernel N", n)
    contract = NymanSearchCellContract(n=n, shared_kernel_n=shared)
    if dict(raw) != contract.to_record():
        raise NymanSearchPlanError("cell contract is not canonical")
    return contract


def _policy_record(
    generation_bits: int,
    replay_bits: int,
    coefficient_grid_bits: int,
    bound_grid_bits: int,
    slack_bits: int,
) -> dict[str, Any]:
    return {
        "generation_precision_bits": str(generation_bits),
        "replay_precision_bits": str(replay_bits),
        "coefficient_grid": {
            "denominator_exponent": str(coefficient_grid_bits),
            "quantum": f"2^-{coefficient_grid_bits}",
            "rounding": "exact-rational-round-ties-to-even",
        },
        "bound_grid": {
            "denominator_exponent": str(bound_grid_bits),
            "quantum": f"2^-{bound_grid_bits}",
            "upper_rounding": "least-grid-point-strictly-above-Arb-upper-endpoint",
        },
        "target_lower_slack": {
            "bits": str(slack_bits),
            "value": f"2^-{slack_bits}",
            "rule": "L=U-2^-target_lower_slack_bits",
        },
        "terminal_gate": (
            "direct exact-dyadic E(c)<=U certificate and positive-definite "
            "fixed-order interval LDL certificate for augmented K_N(L)"
        ),
        "failed_target_behavior": "INCONCLUSIVE; no adaptive weakening",
    }


def _parse_policy(raw: Any) -> tuple[int, int, int, int, int]:
    if not isinstance(raw, Mapping):
        raise NymanSearchPlanError("policy must be an object")
    generation = _parse_integer(raw.get("generation_precision_bits"), "generation bits", 64)
    replay = _parse_integer(raw.get("replay_precision_bits"), "replay bits", generation + 1)
    coefficient = raw.get("coefficient_grid")
    bound = raw.get("bound_grid")
    slack = raw.get("target_lower_slack")
    if not isinstance(coefficient, Mapping) or not isinstance(bound, Mapping) or not isinstance(slack, Mapping):
        raise NymanSearchPlanError("dyadic grid policy objects are required")
    coefficient_bits = _parse_integer(
        coefficient.get("denominator_exponent"), "coefficient grid bits", 1
    )
    bound_bits = _parse_integer(bound.get("denominator_exponent"), "bound grid bits", 1)
    slack_bits = _parse_integer(slack.get("bits"), "target lower slack bits", 1)
    expected = _policy_record(generation, replay, coefficient_bits, bound_bits, slack_bits)
    if dict(raw) != expected:
        raise NymanSearchPlanError("policy is not canonical")
    if slack_bits > bound_bits:
        raise NymanSearchPlanError("target lower slack must lie on the bound grid")
    return generation, replay, coefficient_bits, bound_bits, slack_bits


def _engine_record(max_n: int) -> dict[str, Any]:
    kernel = _shared_kernel_contract(max_n)
    return {
        "algorithm_id": ENGINE_ID,
        "basis": "natural dilates rho_a(x)={1/(a*x)}, a=1..N",
        "objective": "E(c)=1-2*b^T*c+c^T*G*c",
        "candidate_generation": (
            "approximate Arb solve followed by exact ties-to-even dyadic rounding; "
            "the solve is not sign-bearing"
        ),
        "cell_order": "N ascending",
        "candidate_behavior": "run every planned N; never stop early",
        "parallelism": "serial; FLINT ctx.prec is shared mutable state",
        "shared_kernel": kernel,
        "strict_improvement_role": (
            "U_2N<L_N is an index diagnostic only and never a cell-validity gate"
        ),
    }


def canonicalize_nyman_plan(raw: Mapping[str, Any]) -> dict[str, Any]:
    """Validate, order, and hash a natural-distance search plan."""

    if not isinstance(raw, Mapping) or raw.get("schema") != PLAN_SCHEMA:
        raise NymanSearchPlanError("unexpected Nyman plan schema")
    allowed = {
        "schema",
        "classification",
        "hypothesis_status",
        "cells",
        "policy",
        "engine",
        "backend_contract",
        "frozen_experiment",
        "limitation",
        "plan_sha256",
    }
    if set(raw) - allowed:
        raise NymanSearchPlanError("Nyman plan has unknown fields")
    if raw.get("classification", "EXPLORATORY") != "EXPLORATORY":
        raise NymanSearchPlanError("Nyman plans must remain EXPLORATORY")
    if raw.get("hypothesis_status", "UNRESOLVED") != "UNRESOLVED":
        raise NymanSearchPlanError("Nyman plans cannot resolve RH")
    raw_cells = raw.get("cells")
    if not isinstance(raw_cells, list) or not raw_cells:
        raise NymanSearchPlanError("cells must be a nonempty list")
    contracts = sorted((_contract_from_record(item) for item in raw_cells), key=lambda item: item.n)
    if len({item.n for item in contracts}) != len(contracts):
        raise NymanSearchPlanError("duplicate N cell")
    frozen = raw.get("frozen_experiment")
    if frozen is not None:
        if not isinstance(frozen, Mapping) or dict(frozen) != FROZEN_EXPERIMENT:
            raise NymanSearchPlanError("frozen Nyman experiment metadata changed")
        if tuple(item.n for item in contracts) != FROZEN_N_VALUES:
            raise NymanSearchPlanError("frozen Nyman N sequence changed")
    max_n = max(item.n for item in contracts)
    if any(item.shared_kernel_n != max_n for item in contracts):
        raise NymanSearchPlanError("all cells must bind to one max-N shared kernel")
    policy_values = _parse_policy(raw.get("policy"))
    engine = _engine_record(max_n)
    if raw.get("engine") is not None and (
        not isinstance(raw.get("engine"), Mapping) or dict(raw["engine"]) != engine
    ):
        raise NymanSearchPlanError("Nyman engine contract changed")
    backend = raw.get("backend_contract", _backend_record())
    if not isinstance(backend, Mapping) or set(backend) != {"python_flint", "flint"}:
        raise NymanSearchPlanError("backend contract must name python-flint and FLINT")
    backend_record = {"python_flint": str(backend["python_flint"]), "flint": str(backend["flint"])}
    if not all(backend_record.values()):
        raise NymanSearchPlanError("backend versions are required")
    if frozen is not None:
        if policy_values != (
            FROZEN_GENERATION_BITS,
            FROZEN_REPLAY_BITS,
            FROZEN_COEFFICIENT_GRID_BITS,
            FROZEN_BOUND_GRID_BITS,
            FROZEN_TARGET_LOWER_SLACK_BITS,
        ):
            raise NymanSearchPlanError("frozen Nyman arithmetic policy changed")
    if raw.get("limitation", EXPLORATORY_LIMITATION) != EXPLORATORY_LIMITATION:
        raise NymanSearchPlanError("Nyman limitation changed")
    body: dict[str, Any] = {
        "schema": PLAN_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "cells": [item.to_record() for item in contracts],
        "policy": _policy_record(*policy_values),
        "engine": engine,
        "backend_contract": backend_record,
        "limitation": EXPLORATORY_LIMITATION,
    }
    if frozen is not None:
        body["frozen_experiment"] = FROZEN_EXPERIMENT
    plan_sha256 = content_sha256(body)
    if raw.get("plan_sha256") not in {None, plan_sha256}:
        raise NymanSearchPlanError("Nyman plan hash mismatch")
    return {**body, "plan_sha256": plan_sha256}


def load_nyman_plan(path: Path) -> dict[str, Any]:
    raw = _load_json(path, NymanSearchPlanError, label="Nyman plan")
    return canonicalize_nyman_plan(raw)


def _policy_values(plan: Mapping[str, Any]) -> tuple[int, int, int, int, int]:
    return _parse_policy(plan["policy"])


def _contracts(plan: Mapping[str, Any]) -> list[NymanSearchCellContract]:
    return [_contract_from_record(item) for item in plan["cells"]]


@contextmanager
def _working_precision(precision_bits: int) -> Iterable[None]:
    previous = ctx.prec
    try:
        ctx.prec = precision_bits
        yield
    finally:
        ctx.prec = previous


def _dyadic_fraction(mantissa: Any, exponent: Any) -> Fraction:
    m = int(mantissa)
    e = int(exponent)
    return Fraction(m << e, 1) if e >= 0 else Fraction(m, 1 << -e)


def _arb_upper_fraction(value: Any) -> Fraction:
    midpoint = _dyadic_fraction(*value.mid().man_exp())
    radius = _dyadic_fraction(*value.rad().man_exp())
    return midpoint + radius


def _least_grid_point_strictly_above(value: Any, exponent: int) -> tuple[int, Fraction]:
    scaled = _arb_upper_fraction(value) * (1 << exponent)
    numerator = scaled.numerator // scaled.denominator + 1
    return numerator, Fraction(numerator, 1 << exponent)


def _validate_system(system: Any, expected_n: int) -> None:
    if not isinstance(system, core.NaturalSystem):
        raise NymanSearchCheckpointError("core returned an unexpected natural system")
    if tuple(system.dilates) != tuple(range(1, expected_n + 1)):
        raise NymanSearchCheckpointError("shared natural system dilates changed")
    if len(system.gram) != expected_n or any(len(row) != expected_n for row in system.gram):
        raise NymanSearchCheckpointError("shared natural Gram dimension changed")
    if len(system.target) != expected_n:
        raise NymanSearchCheckpointError("shared natural target dimension changed")


def _prefix_system(system: Any, n: int) -> Any:
    _validate_system(system, len(system.dilates))
    prefix = core.prefix_natural_system(system, n)
    _validate_system(prefix, n)
    if prefix.precision_bits != system.precision_bits:
        raise NymanSearchCheckpointError("shared-kernel prefix precision changed")
    return prefix


def _kernel_prefix_sha256(system: Any, n: int) -> str:
    """Commit to the exact Arb enclosures in one leading kernel prefix."""

    prefix = _prefix_system(system, n)
    body = {
        "precision_bits": str(prefix.precision_bits),
        "dilates": [str(value) for value in prefix.dilates],
        "target": [arb_to_dyadic(value) for value in prefix.target],
        "gram_upper_triangle": [
            {
                "row": str(row),
                "column": str(column),
                "value": arb_to_dyadic(prefix.gram[row][column]),
            }
            for row in range(n)
            for column in range(row, n)
        ],
    }
    return content_sha256(body)


def _kernel_manifest(
    system: Any,
    contracts: Sequence[NymanSearchCellContract],
    precision_bits: int,
) -> dict[str, Any]:
    max_n = max(item.shared_kernel_n for item in contracts)
    _validate_system(system, max_n)
    if system.precision_bits != precision_bits:
        raise NymanSearchCheckpointError("shared kernel precision changed")
    return _with_payload_hash(
        {
            "schema": KERNEL_SCHEMA,
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "precision_bits": str(precision_bits),
            "max_n": str(max_n),
            "kernel_contract_sha256": _shared_kernel_contract(max_n)[
                "kernel_contract_sha256"
            ],
            "core_system_content_sha256": core.natural_system_content_sha256(
                system
            ),
            "construction": (
                "one canonical full max-N kernel; each listed commitment is "
                "hashed from an exact leading in-memory prefix"
            ),
            "prefixes": [
                {
                    "n": str(contract.n),
                    "core_system_content_sha256": (
                        core.natural_system_content_sha256(
                            _prefix_system(system, contract.n)
                        )
                    ),
                    "prefix_kernel_sha256": _kernel_prefix_sha256(
                        system, contract.n
                    ),
                }
                for contract in contracts
            ],
        }
    )


def _kernel_cell_binding(
    kernel: Mapping[str, Any], contract: NymanSearchCellContract
) -> dict[str, Any]:
    matches = [
        item
        for item in kernel.get("prefixes", ())
        if isinstance(item, Mapping) and item.get("n") == str(contract.n)
    ]
    if len(matches) != 1:
        raise NymanSearchCheckpointError("shared kernel prefix commitment is missing")
    return {
        "path": "kernels/generation.json",
        "payload_sha256": kernel["payload_sha256"],
        "kernel_contract_sha256": contract.to_record()[
            "shared_kernel_binding"
        ]["kernel_contract_sha256"],
        "core_system_content_sha256": matches[0][
            "core_system_content_sha256"
        ],
        "prefix_kernel_sha256": matches[0]["prefix_kernel_sha256"],
        "prefix_length": str(contract.n),
    }


def _rounded_numerators(values: Sequence[Any], exponent: int) -> tuple[int, ...]:
    rounded = core.round_arb_vector_to_dyadic(values, exponent)
    if not isinstance(rounded, Sequence) or isinstance(rounded, (str, bytes)):
        raise NymanSearchCheckpointError("core dyadic rounding returned no vector")
    numerators: list[int] = []
    for value in rounded:
        if isinstance(value, bool):
            raise NymanSearchCheckpointError("core dyadic numerator is not an integer")
        if isinstance(value, int):
            numerators.append(value)
            continue
        if isinstance(value, Fraction):
            scaled = value * (1 << exponent)
            if scaled.denominator != 1:
                raise NymanSearchCheckpointError("core rounded value is off the coefficient grid")
            numerators.append(scaled.numerator)
            continue
        raise NymanSearchCheckpointError("core returned an unsupported dyadic value")
    return tuple(numerators)


def _candidate_record(
    numerators: Sequence[int],
    coefficient_exponent: int,
    upper_numerator: int,
    lower_numerator: int,
    bound_exponent: int,
) -> dict[str, Any]:
    body = {
        "coefficients": {
            "numerators": [str(value) for value in numerators],
            "denominator_exponent": str(coefficient_exponent),
        },
        "upper_bound": _dyadic_record(upper_numerator, bound_exponent),
        "lower_bound": _dyadic_record(lower_numerator, bound_exponent),
    }
    return {**body, "candidate_sha256": content_sha256(body)}


def _parse_candidate(
    raw: Any,
    contract: NymanSearchCellContract,
    policy: tuple[int, int, int, int, int],
    *,
    error_type: type[NymanSearchError],
) -> tuple[tuple[int, ...], Fraction, Fraction]:
    if not isinstance(raw, Mapping) or set(raw) != {
        "coefficients",
        "upper_bound",
        "lower_bound",
        "candidate_sha256",
    }:
        raise error_type("candidate must be a canonical exact-dyadic object")
    coefficient = raw.get("coefficients")
    if not isinstance(coefficient, Mapping) or set(coefficient) != {
        "numerators",
        "denominator_exponent",
    }:
        raise error_type("coefficient vector must use one dyadic denominator")
    _, _, coefficient_bits, bound_bits, slack_bits = policy
    exponent = _parse_integer(
        coefficient.get("denominator_exponent"),
        "coefficient denominator exponent",
        0,
        error_type=error_type,
    )
    if exponent != coefficient_bits:
        raise error_type("coefficient vector is off the frozen grid")
    raw_numerators = coefficient.get("numerators")
    if not isinstance(raw_numerators, list) or len(raw_numerators) != contract.n:
        raise error_type("coefficient vector dimension changed")
    numerators = tuple(
        _parse_signed_integer(
            value, "coefficient numerator", error_type=error_type
        )
        for value in raw_numerators
    )
    _, upper = _parse_dyadic_record(
        raw.get("upper_bound"),
        "upper bound",
        expected_exponent=bound_bits,
        error_type=error_type,
    )
    _, lower = _parse_dyadic_record(
        raw.get("lower_bound"),
        "lower bound",
        expected_exponent=bound_bits,
        error_type=error_type,
    )
    if upper - lower != Fraction(1, 1 << slack_bits):
        raise error_type("stored lower target does not equal U minus frozen slack")
    body = {key: value for key, value in raw.items() if key != "candidate_sha256"}
    if raw.get("candidate_sha256") != content_sha256(body):
        raise error_type("candidate hash mismatch")
    return numerators, upper, lower


def _validate_core_certificate(record: Any, side: str) -> str:
    if not isinstance(record, Mapping):
        raise NymanSearchCheckpointError(f"core {side} certificate must be an object")
    _require_payload_hash(record, NymanSearchCheckpointError)
    if record.get("classification") != "EXPLORATORY" or record.get("hypothesis_status") != "UNRESOLVED":
        raise NymanSearchCheckpointError(f"core {side} certificate classification changed")
    allowed = {
        "upper": {"UPPER_BOUND_CERTIFIED", "INCONCLUSIVE"},
        "lower": {"LOWER_BOUND_CERTIFIED", "INCONCLUSIVE"},
    }[side]
    decision = record.get("decision")
    if decision not in allowed:
        raise NymanSearchCheckpointError(f"core {side} certificate decision changed")
    checks = record.get("checks")
    if not isinstance(checks, Mapping):
        raise NymanSearchCheckpointError(f"core {side} certificate checks are missing")
    check_name = {
        "upper": "energy_at_most_claimed_upper_bound",
        "lower": "fixed_order_interval_ldlt_positive",
    }[side]
    expected = decision != "INCONCLUSIVE"
    if checks.get(check_name) is not expected:
        raise NymanSearchCheckpointError(f"core {side} decision disagrees with its check")
    return str(decision)


def _certify_candidate(
    contract: NymanSearchCellContract,
    candidate: Mapping[str, Any],
    policy: tuple[int, int, int, int, int],
    precision_bits: int,
    system: Any,
) -> tuple[dict[str, Any], dict[str, Any]]:
    numerators, upper, lower = _parse_candidate(
        candidate, contract, policy, error_type=NymanSearchCheckpointError
    )
    coefficient_bits = policy[2]
    dilates = tuple(range(1, contract.n + 1))
    _validate_system(system, contract.n)
    if system.precision_bits != precision_bits:
        raise NymanSearchCheckpointError("certificate prefix precision changed")
    upper_certificate = core.certify_dyadic_upper_bound(
        dilates,
        numerators,
        coefficient_bits,
        upper,
        precision_bits=precision_bits,
        system=system,
    )
    lower_certificate = core.certify_augmented_lower_bound(
        dilates,
        lower,
        precision_bits=precision_bits,
        system=system,
    )
    _validate_core_certificate(upper_certificate, "upper")
    _validate_core_certificate(lower_certificate, "lower")
    return dict(upper_certificate), dict(lower_certificate)


def _status_from_certificates(upper: Mapping[str, Any], lower: Mapping[str, Any]) -> str:
    upper_decision = _validate_core_certificate(upper, "upper")
    lower_decision = _validate_core_certificate(lower, "lower")
    if upper_decision == "UPPER_BOUND_CERTIFIED" and lower_decision == "LOWER_BOUND_CERTIFIED":
        return "FINITE_DISTANCE_BRACKET_CERTIFIED"
    return "INCONCLUSIVE"


def _generate_cell(
    contract: NymanSearchCellContract,
    plan: Mapping[str, Any],
    shared_system: Any,
    kernel: Mapping[str, Any],
) -> dict[str, Any]:
    policy = _policy_values(plan)
    generation_bits, replay_bits, coefficient_bits, bound_bits, slack_bits = policy
    prefix = _prefix_system(shared_system, contract.n)
    with _working_precision(generation_bits):
        matrix = arb_mat([list(row) for row in prefix.gram])
        right = arb_mat([[value] for value in prefix.target])
        try:
            solution = matrix.solve(right)
        except Exception as exc:  # pragma: no cover - backend-specific failure text
            raise NymanSearchCheckpointError("approximate natural-system solve failed") from exc
        numerators = _rounded_numerators(
            tuple(solution[index, 0] for index in range(contract.n)),
            coefficient_bits,
        )
        coefficients = core.dyadic_coefficients(numerators, coefficient_bits)
        energy = core.evaluate_natural_distance(prefix, coefficients)
        upper_numerator, _ = _least_grid_point_strictly_above(energy, bound_bits)
    slack_numerator = 1 << (bound_bits - slack_bits)
    candidate = _candidate_record(
        numerators,
        coefficient_bits,
        upper_numerator,
        upper_numerator - slack_numerator,
        bound_bits,
    )
    upper_certificate, lower_certificate = _certify_candidate(
        contract, candidate, policy, generation_bits, prefix
    )
    status = _status_from_certificates(upper_certificate, lower_certificate)
    body = {
        "schema": CELL_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "plan_sha256": plan["plan_sha256"],
        "cell_contract": contract.to_record(),
        "candidate": candidate,
        "generation": {
            "precision_bits": str(generation_bits),
            "solver_role": "approximate-untrusted-candidate-generator-only",
            "shared_kernel_binding": _kernel_cell_binding(kernel, contract),
            "upper_certificate": upper_certificate,
            "lower_certificate": lower_certificate,
        },
        "replay_precision_bits": str(replay_bits),
        "status": status,
        "terminal_statement": (
            "L_N < d_N^2 <= U_N" if status == "FINITE_DISTANCE_BRACKET_CERTIFIED" else None
        ),
        "limitation": EXPLORATORY_LIMITATION,
    }
    return _with_payload_hash(body)


def _derive_cell(
    contract: NymanSearchCellContract,
    plan: Mapping[str, Any],
    artifact: Mapping[str, Any],
    shared_system: Any,
    kernel: Mapping[str, Any],
    *,
    error_type: type[NymanSearchError],
) -> dict[str, Any]:
    _require_payload_hash(artifact, error_type)
    if artifact.get("schema") != CELL_SCHEMA:
        raise error_type("unexpected Nyman cell schema")
    if artifact.get("classification") != "EXPLORATORY" or artifact.get("hypothesis_status") != "UNRESOLVED":
        raise error_type("Nyman cell was improperly promoted")
    if artifact.get("plan_sha256") != plan["plan_sha256"]:
        raise error_type("Nyman cell plan binding changed")
    if artifact.get("cell_contract") != contract.to_record():
        raise error_type("Nyman cell contract changed")
    if artifact.get("limitation") != EXPLORATORY_LIMITATION:
        raise error_type("Nyman cell limitation changed")
    policy = _policy_values(plan)
    generation_bits, replay_bits, _, _, _ = policy
    if artifact.get("replay_precision_bits") != str(replay_bits):
        raise error_type("Nyman cell replay precision changed")
    generation = artifact.get("generation")
    if not isinstance(generation, Mapping) or set(generation) != {
        "precision_bits",
        "solver_role",
        "shared_kernel_binding",
        "upper_certificate",
        "lower_certificate",
    }:
        raise error_type("Nyman generation evidence changed")
    if generation.get("precision_bits") != str(generation_bits):
        raise error_type("Nyman generation precision changed")
    if generation.get("solver_role") != "approximate-untrusted-candidate-generator-only":
        raise error_type("approximate solve was given a certificate role")
    if generation.get("shared_kernel_binding") != _kernel_cell_binding(
        kernel, contract
    ):
        raise error_type("shared-kernel prefix binding changed")
    _parse_candidate(artifact.get("candidate"), contract, policy, error_type=error_type)
    prefix = _prefix_system(shared_system, contract.n)
    try:
        reproduced_upper, reproduced_lower = _certify_candidate(
            contract,
            artifact["candidate"],
            policy,
            generation_bits,
            prefix,
        )
    except NymanSearchCheckpointError as exc:
        raise error_type(str(exc)) from exc
    if generation.get("upper_certificate") != reproduced_upper:
        raise error_type("stored direct upper certificate does not reproduce")
    if generation.get("lower_certificate") != reproduced_lower:
        raise error_type("stored augmented lower certificate does not reproduce")
    expected_status = _status_from_certificates(reproduced_upper, reproduced_lower)
    if artifact.get("status") != expected_status:
        raise error_type("Nyman cell status does not follow from its certificates")
    expected_statement = "L_N < d_N^2 <= U_N" if expected_status == "FINITE_DISTANCE_BRACKET_CERTIFIED" else None
    if artifact.get("terminal_statement") != expected_statement:
        raise error_type("Nyman terminal statement changed")
    expected_body = {key: value for key, value in artifact.items() if key != "payload_sha256"}
    # Reject unknown fields by rebuilding the exact allowed shell around the
    # reproduced evidence.  Candidate validation above handles its inner shape.
    allowed = {
        "schema",
        "classification",
        "hypothesis_status",
        "plan_sha256",
        "cell_contract",
        "candidate",
        "generation",
        "replay_precision_bits",
        "status",
        "terminal_statement",
        "limitation",
    }
    if set(expected_body) != allowed:
        raise error_type("Nyman cell fields changed")
    return dict(artifact)


def _atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def _load_json(
    path: Path,
    error_type: type[NymanSearchError],
    *,
    label: str = "Nyman artifact",
) -> dict[str, Any]:
    def reject_duplicate_keys(
        pairs: list[tuple[str, Any]],
    ) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                raise error_type(f"duplicate JSON object key: {key}")
            value[key] = item
        return value

    def reject_nonstandard_constant(value: str) -> Any:
        raise error_type(f"nonstandard JSON constant: {value}")

    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=reject_duplicate_keys,
            parse_constant=reject_nonstandard_constant,
        )
    except error_type:
        raise
    except (OSError, json.JSONDecodeError) as exc:
        raise error_type(f"cannot read {label}: {path}") from exc
    if not isinstance(value, dict):
        raise error_type(f"{label} must contain an object")
    return value


def _runtime_backend(plan: Mapping[str, Any]) -> None:
    if plan["backend_contract"] != _backend_record():
        raise NymanSearchPlanError("installed python-flint/FLINT versions do not match the plan")


def _run_manifest(plan: Mapping[str, Any]) -> dict[str, Any]:
    return _with_payload_hash(
        {
            "schema": RUN_SCHEMA,
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "plan": {"path": "plan.json", "plan_sha256": plan["plan_sha256"]},
            "backend": _backend_record(),
            "engine_id": ENGINE_ID,
            "candidate_behavior": "run every planned N; never stop early",
            "limitation": EXPLORATORY_LIMITATION,
        }
    )


def _strict_improvement_diagnostics(
    contracts: Sequence[NymanSearchCellContract],
    completed: Mapping[str, Mapping[str, Any]],
    policy: tuple[int, int, int, int, int],
) -> list[dict[str, Any]]:
    by_n = {item.n: item for item in contracts}
    diagnostics: list[dict[str, Any]] = []
    for smaller_n in sorted(by_n):
        larger_n = 2 * smaller_n
        if larger_n not in by_n:
            continue
        smaller = completed.get(by_n[smaller_n].cell_id)
        larger = completed.get(by_n[larger_n].cell_id)
        decision = "UNAVAILABLE_IN_INCOMPLETE_BATCH"
        comparison: dict[str, Any] | None = None
        if smaller is not None and larger is not None:
            _, _, smaller_lower = _parse_candidate(
                smaller["candidate"], by_n[smaller_n], policy, error_type=NymanSearchCheckpointError
            )
            _, larger_upper, _ = _parse_candidate(
                larger["candidate"], by_n[larger_n], policy, error_type=NymanSearchCheckpointError
            )
            comparison = {
                "u_2n": _fraction_record(larger_upper),
                "l_n": _fraction_record(smaller_lower),
            }
            if (
                smaller.get("status")
                == "FINITE_DISTANCE_BRACKET_CERTIFIED"
                and larger.get("status")
                == "FINITE_DISTANCE_BRACKET_CERTIFIED"
            ):
                decision = (
                    "STRICT_IMPROVEMENT_CERTIFIED"
                    if larger_upper < smaller_lower
                    else "STRICT_IMPROVEMENT_NOT_CERTIFIED"
                )
            else:
                decision = "UNAVAILABLE_UNCERTIFIED_BRACKET"
        diagnostics.append(
            {
                "smaller_n": str(smaller_n),
                "larger_n": str(larger_n),
                "inequality": "U_2N < L_N",
                "decision": decision,
                "exact_comparison": comparison,
                "cell_validity_gate": False,
            }
        )
    return diagnostics


def _index_payload(
    plan: Mapping[str, Any],
    contracts: Sequence[NymanSearchCellContract],
    completed: Mapping[str, Mapping[str, Any]],
    kernel: Mapping[str, Any] | None,
) -> dict[str, Any]:
    references: list[dict[str, Any]] = []
    statuses: list[str] = []
    for contract in contracts:
        artifact = completed.get(contract.cell_id)
        if artifact is None:
            continue
        statuses.append(str(artifact["status"]))
        references.append(
            {
                "cell_id": contract.cell_id,
                "cell_contract_sha256": contract.to_record()["cell_contract_sha256"],
                "path": f"cells/{contract.cell_id}.json",
                "payload_sha256": artifact["payload_sha256"],
                "candidate_sha256": artifact["candidate"]["candidate_sha256"],
                "status": artifact["status"],
            }
        )
    missing = [item.cell_id for item in contracts if item.cell_id not in completed]
    counts = {
        "FINITE_DISTANCE_BRACKET_CERTIFIED": str(statuses.count("FINITE_DISTANCE_BRACKET_CERTIFIED")),
        "INCONCLUSIVE": str(statuses.count("INCONCLUSIVE")),
    }
    if missing:
        conclusion = "INCOMPLETE_BOUNDED_NYMAN_SEARCH"
    elif counts["INCONCLUSIVE"] != "0":
        conclusion = "FINITE_DISTANCE_BATCH_CONTAINS_INCONCLUSIVE_CELLS"
    else:
        conclusion = "ALL_FINITE_DISTANCE_BRACKETS_CERTIFIED"
    if completed and kernel is None:
        raise NymanSearchCheckpointError(
            "completed cells require a shared generation kernel"
        )
    kernel_reference = None
    if kernel is not None:
        kernel_reference = {
            "path": "kernels/generation.json",
            "payload_sha256": kernel["payload_sha256"],
            "precision_bits": kernel["precision_bits"],
            "kernel_contract_sha256": kernel["kernel_contract_sha256"],
            "core_system_content_sha256": kernel[
                "core_system_content_sha256"
            ],
        }
    return _with_payload_hash(
        {
            "schema": INDEX_SCHEMA,
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "conclusion": conclusion,
            "plan": {"path": "plan.json", "plan_sha256": plan["plan_sha256"]},
            "run": {"path": "run.json"},
            "backend": _backend_record(),
            "generation_kernel": kernel_reference,
            "cells": references,
            "missing_cell_ids": missing,
            "progress": {
                "planned_cells": str(len(contracts)),
                "completed_cells": str(len(completed)),
                "status_counts": counts,
            },
            "strict_improvement_diagnostics": _strict_improvement_diagnostics(
                contracts, completed, _policy_values(plan)
            ),
            "strict_improvement_role": "separate exact diagnostic; never a cell-validity gate",
            "limitation": EXPLORATORY_LIMITATION,
        }
    )


def run_nyman_search(
    plan: Mapping[str, Any],
    checkpoint_dir: Path,
    *,
    resume: bool = False,
    max_cells: int | None = None,
) -> dict[str, Any]:
    """Run or resume the complete bounded natural-distance experiment."""

    canonical_plan = canonicalize_nyman_plan(plan)
    _runtime_backend(canonical_plan)
    if max_cells is not None:
        if isinstance(max_cells, bool) or not isinstance(max_cells, int):
            raise TypeError("max_cells must be an integer or None")
        if max_cells < 0:
            raise ValueError("max_cells must be nonnegative")
    root = checkpoint_dir.resolve()
    expected_run = _run_manifest(canonical_plan)
    if root.exists() and any(root.iterdir()):
        if not resume:
            raise NymanSearchCheckpointError("checkpoint already exists; pass resume=True")
        stored_plan_raw = _load_json(root / "plan.json", NymanSearchCheckpointError)
        try:
            stored_plan = canonicalize_nyman_plan(stored_plan_raw)
        except NymanSearchPlanError as exc:
            raise NymanSearchCheckpointError("stored Nyman plan is invalid") from exc
        if stored_plan_raw != stored_plan or stored_plan != canonical_plan:
            raise NymanSearchCheckpointError("checkpoint Nyman plan changed")
        stored_run = _load_json(root / "run.json", NymanSearchCheckpointError)
        _require_payload_hash(stored_run, NymanSearchCheckpointError)
        if stored_run != expected_run:
            raise NymanSearchCheckpointError("checkpoint backend, plan, or engine changed")
    else:
        if resume:
            raise NymanSearchCheckpointError("resume requested without a checkpoint")
        root.mkdir(parents=True, exist_ok=True)
        _atomic_write_json(root / "plan.json", canonical_plan)
        _atomic_write_json(root / "run.json", expected_run)

    contracts = _contracts(canonical_plan)
    existing_contracts = [
        contract
        for contract in contracts
        if (root / "cells" / f"{contract.cell_id}.json").exists()
    ]
    pending = [
        contract for contract in contracts if contract not in existing_contracts
    ]
    if max_cells is not None:
        pending = pending[:max_cells]
    needs_kernel = bool(existing_contracts or pending)
    shared_system: Any | None = None
    kernel: dict[str, Any] | None = None
    if needs_kernel:
        generation_bits = _policy_values(canonical_plan)[0]
        max_n = max(item.shared_kernel_n for item in contracts)
        with _working_precision(generation_bits):
            shared_system = core.build_natural_system(
                tuple(range(1, max_n + 1))
            )
        _validate_system(shared_system, max_n)
        kernel = _kernel_manifest(
            shared_system, contracts, generation_bits
        )
        kernel_path = root / "kernels" / "generation.json"
        if kernel_path.exists():
            stored_kernel = _load_json(
                kernel_path, NymanSearchCheckpointError
            )
            _require_payload_hash(
                stored_kernel, NymanSearchCheckpointError
            )
            if stored_kernel != kernel:
                raise NymanSearchCheckpointError(
                    "stored shared generation kernel does not reproduce"
                )
        else:
            _atomic_write_json(kernel_path, kernel)

    completed: dict[str, dict[str, Any]] = {}
    for contract in existing_contracts:
        path = root / "cells" / f"{contract.cell_id}.json"
        assert shared_system is not None and kernel is not None
        artifact = _load_json(path, NymanSearchCheckpointError)
        completed[contract.cell_id] = _derive_cell(
            contract,
            canonical_plan,
            artifact,
            shared_system,
            kernel,
            error_type=NymanSearchCheckpointError,
        )

    if pending:
        assert shared_system is not None and kernel is not None
        for contract in pending:
            artifact = _generate_cell(
                contract, canonical_plan, shared_system, kernel
            )
            _atomic_write_json(root / "cells" / f"{contract.cell_id}.json", artifact)
            completed[contract.cell_id] = artifact

    index = _index_payload(canonical_plan, contracts, completed, kernel)
    _atomic_write_json(root / "index.json", index)
    return index


def _safe_artifact_path(root: Path, raw_path: Any) -> Path:
    if not isinstance(raw_path, str) or not raw_path or "\\" in raw_path:
        raise NymanSearchVerificationError("invalid Nyman artifact path")
    relative = Path(raw_path)
    if relative.is_absolute() or any(part in {"", ".", ".."} for part in relative.parts):
        raise NymanSearchVerificationError("invalid Nyman artifact path")
    resolved = (root / relative).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise NymanSearchVerificationError("Nyman artifact path escapes checkpoint") from exc
    return resolved


def verify_nyman_search(
    index: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
    *,
    replay_precision_bits: int = FROZEN_REPLAY_BITS,
) -> dict[str, Any]:
    """Replay stored ``c``, ``U``, and ``L`` without regenerating candidates."""

    root = checkpoint_dir.resolve()
    if isinstance(index, Path):
        index_record = _load_json(index, NymanSearchVerificationError)
    elif isinstance(index, Mapping):
        index_record = dict(index)
    else:
        raise TypeError("index must be a mapping or Path")
    _require_payload_hash(index_record, NymanSearchVerificationError)
    if index_record.get("schema") != INDEX_SCHEMA:
        raise NymanSearchVerificationError("unexpected Nyman index schema")
    if index_record.get("classification") != "EXPLORATORY" or index_record.get("hypothesis_status") != "UNRESOLVED":
        raise NymanSearchVerificationError("Nyman index was improperly promoted")
    plan_reference = index_record.get("plan")
    if not isinstance(plan_reference, Mapping):
        raise NymanSearchVerificationError("Nyman index plan reference is missing")
    plan_path = _safe_artifact_path(root, plan_reference.get("path"))
    stored_plan = _load_json(plan_path, NymanSearchVerificationError)
    try:
        plan = canonicalize_nyman_plan(stored_plan)
        _runtime_backend(plan)
    except NymanSearchPlanError as exc:
        raise NymanSearchVerificationError(str(exc)) from exc
    if stored_plan != plan:
        raise NymanSearchVerificationError("stored Nyman plan is not canonical")
    if plan_reference.get("plan_sha256") != plan["plan_sha256"]:
        raise NymanSearchVerificationError("Nyman index plan hash changed")
    requested = _parse_integer(
        replay_precision_bits,
        "replay precision bits",
        64,
        error_type=NymanSearchVerificationError,
    )
    if requested != _policy_values(plan)[1]:
        raise NymanSearchVerificationError("replay must use the frozen requested precision")
    run_reference = index_record.get("run")
    if not isinstance(run_reference, Mapping):
        raise NymanSearchVerificationError("Nyman index run reference is missing")
    run = _load_json(
        _safe_artifact_path(root, run_reference.get("path")),
        NymanSearchVerificationError,
    )
    _require_payload_hash(run, NymanSearchVerificationError)
    if run != _run_manifest(plan):
        raise NymanSearchVerificationError("Nyman run manifest changed")

    contracts = _contracts(plan)
    by_id = {item.cell_id: item for item in contracts}
    references = index_record.get("cells")
    if not isinstance(references, list):
        raise NymanSearchVerificationError("Nyman index cells must be a list")
    generation_bits = _policy_values(plan)[0]
    max_n = max(item.shared_kernel_n for item in contracts)
    kernel_reference = index_record.get("generation_kernel")
    generation_system: Any | None = None
    reproduced_kernel: dict[str, Any] | None = None
    if references:
        with _working_precision(generation_bits):
            generation_system = core.build_natural_system(
                tuple(range(1, max_n + 1))
            )
        _validate_system(generation_system, max_n)
        reproduced_kernel = _kernel_manifest(
            generation_system, contracts, generation_bits
        )
        expected_kernel_reference = {
            "path": "kernels/generation.json",
            "payload_sha256": reproduced_kernel["payload_sha256"],
            "precision_bits": str(generation_bits),
            "kernel_contract_sha256": reproduced_kernel[
                "kernel_contract_sha256"
            ],
            "core_system_content_sha256": reproduced_kernel[
                "core_system_content_sha256"
            ],
        }
        if kernel_reference != expected_kernel_reference:
            raise NymanSearchVerificationError(
                "Nyman generation-kernel reference changed"
            )
        stored_kernel = _load_json(
            _safe_artifact_path(root, kernel_reference["path"]),
            NymanSearchVerificationError,
        )
        _require_payload_hash(
            stored_kernel, NymanSearchVerificationError
        )
        if stored_kernel != reproduced_kernel:
            raise NymanSearchVerificationError(
                "stored shared generation kernel does not reproduce"
            )
    elif kernel_reference is not None:
        raise NymanSearchVerificationError(
            "empty Nyman index cannot reference a generation kernel"
        )

    completed: dict[str, dict[str, Any]] = {}
    seen: set[str] = set()
    for reference in references:
        if not isinstance(reference, Mapping):
            raise NymanSearchVerificationError("Nyman cell reference must be an object")
        cell_id = reference.get("cell_id")
        if cell_id not in by_id or cell_id in seen:
            raise NymanSearchVerificationError("unknown or duplicate Nyman cell reference")
        seen.add(str(cell_id))
        contract = by_id[str(cell_id)]
        if reference.get("cell_contract_sha256") != contract.to_record()["cell_contract_sha256"]:
            raise NymanSearchVerificationError("Nyman cell contract hash changed")
        expected_path = f"cells/{contract.cell_id}.json"
        if reference.get("path") != expected_path:
            raise NymanSearchVerificationError("Nyman cell path changed")
        artifact = _load_json(
            _safe_artifact_path(root, reference.get("path")),
            NymanSearchVerificationError,
        )
        assert generation_system is not None
        assert reproduced_kernel is not None
        artifact = _derive_cell(
            contract,
            plan,
            artifact,
            generation_system,
            reproduced_kernel,
            error_type=NymanSearchVerificationError,
        )
        if reference.get("payload_sha256") != artifact["payload_sha256"]:
            raise NymanSearchVerificationError("Nyman cell payload reference changed")
        if reference.get("candidate_sha256") != artifact["candidate"]["candidate_sha256"]:
            raise NymanSearchVerificationError("Nyman candidate reference changed")
        if reference.get("status") != artifact["status"]:
            raise NymanSearchVerificationError("Nyman cell status reference changed")
        completed[contract.cell_id] = artifact

    expected_index = _index_payload(plan, contracts, completed, reproduced_kernel)
    if index_record != expected_index:
        raise NymanSearchVerificationError("Nyman index does not match referenced artifacts")

    replay_counts = {
        "FINITE_DISTANCE_BRACKET_CERTIFIED": 0,
        "INCONCLUSIVE": 0,
    }
    replay_kernel_sha256: str | None = None
    if completed:
        with _working_precision(requested):
            replay_system = core.build_natural_system(
                tuple(range(1, max_n + 1))
            )
        _validate_system(replay_system, max_n)
        replay_kernel = _kernel_manifest(
            replay_system, contracts, requested
        )
        replay_kernel_sha256 = str(replay_kernel["payload_sha256"])
        for contract in contracts:
            artifact = completed.get(contract.cell_id)
            if artifact is None:
                continue
            replay_prefix = _prefix_system(replay_system, contract.n)
            try:
                replay_upper, replay_lower = _certify_candidate(
                    contract,
                    artifact["candidate"],
                    _policy_values(plan),
                    requested,
                    replay_prefix,
                )
                replay_status = _status_from_certificates(
                    replay_upper, replay_lower
                )
            except NymanSearchCheckpointError as exc:
                raise NymanSearchVerificationError(str(exc)) from exc
            replay_counts[replay_status] += 1
            if (
                artifact["status"]
                == "FINITE_DISTANCE_BRACKET_CERTIFIED"
                and replay_status != "FINITE_DISTANCE_BRACKET_CERTIFIED"
            ):
                raise NymanSearchVerificationError(
                    "terminal Nyman bracket failed requested-precision replay"
                )
    return {
        "classification": "REPRODUCED_EXPLORATORY_NYMAN_SEARCH",
        "hypothesis_status": "UNRESOLVED",
        "conclusion": index_record["conclusion"],
        "verified_cells": str(len(completed)),
        "replay_precision_bits": str(requested),
        "replay_status_counts": {key: str(value) for key, value in replay_counts.items()},
        "generation_kernel_payload_sha256": (
            reproduced_kernel["payload_sha256"]
            if reproduced_kernel is not None
            else None
        ),
        "replay_kernel_payload_sha256": replay_kernel_sha256,
        "stored_candidates_replayed_without_regeneration": True,
        "same_backend_replay_only": True,
    }


__all__ = [
    "CELL_CONTRACT_SCHEMA",
    "CELL_SCHEMA",
    "ENGINE_ID",
    "EXPLORATORY_LIMITATION",
    "FROZEN_BOUND_GRID_BITS",
    "FROZEN_COEFFICIENT_GRID_BITS",
    "FROZEN_EXPERIMENT",
    "FROZEN_GENERATION_BITS",
    "FROZEN_N_VALUES",
    "FROZEN_REPLAY_BITS",
    "FROZEN_TARGET_LOWER_SLACK_BITS",
    "INDEX_SCHEMA",
    "KERNEL_SCHEMA",
    "NymanSearchCellContract",
    "NymanSearchCheckpointError",
    "NymanSearchError",
    "NymanSearchPlanError",
    "NymanSearchVerificationError",
    "PLAN_SCHEMA",
    "RUN_SCHEMA",
    "canonicalize_nyman_plan",
    "load_nyman_plan",
    "run_nyman_search",
    "verify_nyman_search",
]
