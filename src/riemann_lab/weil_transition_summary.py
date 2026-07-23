"""Compact deterministic summaries for completed v2 transition searches.

The summary path checks the canonical checkpoint structure and stored
certificate state, but deliberately does not rebuild a Weil matrix or replay
any numerical calculation.  It is therefore an integrity and reporting layer,
not an independent reproduction of the finite-cell calculations.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping

from .artifacts import content_sha256
from . import weil_search as v1
from . import weil_transition as transition


SUMMARY_SCHEMA = "rh-lab/weil-transition-summary/v2"
FROZEN_PLAN_SHA256 = (
    "33185881cc619fda99b7835f5e71422b0c5cf745cb5411162638a8e08dd05336"
)
SUMMARY_LIMITATION = (
    f"{transition.EXPLORATORY_LIMITATION} This compact summary checks "
    "canonical checkpoint records and their content hashes without replaying "
    "the numerical calculations."
)


class WeilTransitionSummaryError(ValueError):
    """Raised when a completed checkpoint cannot be summarized safely."""


class WeilTransitionSummaryVerificationError(WeilTransitionSummaryError):
    """Raised when a supplied summary does not canonically regenerate."""


def _load_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise WeilTransitionSummaryError(f"cannot read checkpoint object: {path}") from exc
    if not isinstance(value, dict):
        raise WeilTransitionSummaryError("checkpoint object must be a JSON object")
    return value


def _safe_path(root: Path, relative: Any) -> Path:
    if not isinstance(relative, str) or not relative:
        raise WeilTransitionSummaryError("checkpoint reference path is invalid")
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise WeilTransitionSummaryError("checkpoint reference escapes its root")
    return path


def _canonical_integer(value: Any, name: str) -> int:
    if not isinstance(value, str):
        raise WeilTransitionSummaryError(f"{name} must be a canonical decimal string")
    try:
        parsed = int(value)
    except ValueError as exc:
        raise WeilTransitionSummaryError(f"{name} is not an integer") from exc
    if value != str(parsed):
        raise WeilTransitionSummaryError(f"{name} is not canonical")
    return parsed


def _dyadic_value(mantissa: Any, exponent: Any, name: str) -> Fraction:
    coefficient = _canonical_integer(mantissa, f"{name} mantissa")
    power = _canonical_integer(exponent, f"{name} exponent")
    if power >= 0:
        return Fraction(coefficient * (1 << power), 1)
    return Fraction(coefficient, 1 << -power)


def _dyadic_bounds(record: Any) -> tuple[Fraction, Fraction]:
    if not isinstance(record, Mapping):
        raise WeilTransitionSummaryError("Rump real enclosure must be an object")
    if set(record) != {"display", "dyadic", "is_exact"}:
        raise WeilTransitionSummaryError("Rump real enclosure has noncanonical fields")
    if not isinstance(record.get("display"), str) or not isinstance(
        record.get("is_exact"), bool
    ):
        raise WeilTransitionSummaryError("Rump real enclosure metadata is invalid")
    encoded = record.get("dyadic")
    if not isinstance(encoded, Mapping) or set(encoded) != {
        "mid_mantissa",
        "mid_exponent",
        "radius_mantissa",
        "radius_exponent",
    }:
        raise WeilTransitionSummaryError("Rump real enclosure has invalid dyadic data")
    midpoint = _dyadic_value(
        encoded["mid_mantissa"], encoded["mid_exponent"], "Rump midpoint"
    )
    radius = _dyadic_value(
        encoded["radius_mantissa"], encoded["radius_exponent"], "Rump radius"
    )
    if radius < 0:
        raise WeilTransitionSummaryError("Rump enclosure radius must be nonnegative")
    if record["is_exact"] is not (radius == 0):
        raise WeilTransitionSummaryError(
            "Rump enclosure exactness flag disagrees with its dyadic radius"
        )
    return midpoint - radius, midpoint + radius


def _normalized_enclosure(record: Mapping[str, Any]) -> dict[str, Any]:
    lower, upper = _dyadic_bounds(record)
    midpoint = (lower + upper) / 2
    radius = (upper - lower) / 2
    return {
        "midpoint": _fraction_record(midpoint),
        "radius": _fraction_record(radius),
        "lower_bound": _fraction_record(lower),
        "upper_bound": _fraction_record(upper),
        "is_exact": radius == 0,
    }


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _terminal_attempts(artifact: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    attempts = artifact.get("attempts")
    if not isinstance(attempts, list) or not attempts:
        raise WeilTransitionSummaryError("completed cell has no v1 attempts")
    if not all(isinstance(attempt, Mapping) for attempt in attempts):
        raise WeilTransitionSummaryError("v1 attempt must be an object")
    for attempt in attempts:
        _validate_attempt_fields(attempt)
    last = attempts[-1]
    terminal = (
        last.get("kind") == "SEARCH"
        and last.get("decision")
        in {"FINITE_POSITIVE_CERTIFIED", "INCONCLUSIVE_MAX_PRECISION"}
    ) or (
        last.get("kind") == "CONFIRMATION"
        and last.get("decision")
        in {"NEGATIVE_WITNESS_CONFIRMED", "NEGATIVE_CONFIRMATION_FAILED"}
    )
    if not terminal:
        raise WeilTransitionSummaryError(
            "cell attempts are not terminal; numerical continuation is forbidden"
        )
    return attempts


_ATTEMPT_BASE_FIELDS = {
    "schema",
    "classification",
    "kind",
    "cell",
    "precision_bits",
    "limitation",
    "sequence",
    "payload_sha256",
}
_SEARCH_COMPLETE_FIELDS = _ATTEMPT_BASE_FIELDS | {
    "matrix_construction",
    "matrix_evidence",
    "positive_definiteness",
    "secondary_rump_spectrum",
    "candidate_generation",
    "decision",
}
_SEARCH_CONSTRUCTION_FAILURE_FIELDS = _ATTEMPT_BASE_FIELDS | {
    "matrix_construction",
    "decision",
}
_CONFIRMATION_COMPLETE_FIELDS = _ATTEMPT_BASE_FIELDS | {
    "matrix_construction",
    "matrix_evidence",
    "witness_evaluation",
    "precision_containment",
    "decision",
}
_CONFIRMATION_CONSTRUCTION_FAILURE_FIELDS = _ATTEMPT_BASE_FIELDS | {
    "matrix_construction",
    "witness",
    "decision",
}


def _validate_attempt_fields(attempt: Mapping[str, Any]) -> None:
    """Reject fields that no canonical v1 attempt constructor can emit."""

    kind = attempt.get("kind")
    construction = attempt.get("matrix_construction")
    if construction == {"classification": "COMPLETE"}:
        expected = (
            _SEARCH_COMPLETE_FIELDS
            if kind == "SEARCH"
            else _CONFIRMATION_COMPLETE_FIELDS
            if kind == "CONFIRMATION"
            else None
        )
    elif construction == {
        "classification": "INCONCLUSIVE",
        "failure_code": "SOURCE_ORACLE_OR_QUADRATURE_INCONCLUSIVE",
    }:
        expected = (
            _SEARCH_CONSTRUCTION_FAILURE_FIELDS
            if kind == "SEARCH"
            else _CONFIRMATION_CONSTRUCTION_FAILURE_FIELDS
            if kind == "CONFIRMATION"
            else None
        )
    else:
        raise WeilTransitionSummaryError(
            "v1 attempt matrix-construction state is noncanonical"
        )
    if expected is None or set(attempt) != expected:
        raise WeilTransitionSummaryError("v1 attempt has noncanonical fields")
    if attempt.get("limitation") != v1.EXPLORATORY_LIMITATION:
        raise WeilTransitionSummaryError("v1 attempt limitation changed")


def _canonical_cell(
    contract: transition.WeilTransitionCellContract,
    artifact: dict[str, Any],
) -> dict[str, Any]:
    """Derive stored terminal state without allowing a numerical factory call."""

    _terminal_attempts(artifact)
    try:
        derived = transition._derive_v1_artifact(contract, artifact, replay=False)
    except (transition.WeilTransitionError, v1.WeilSearchError) as exc:
        raise WeilTransitionSummaryError(
            f"cell {contract.contract_id} is not canonical"
        ) from exc
    return derived


def _negative_evaluation(evaluation: Any) -> bool:
    return isinstance(evaluation, Mapping) and (
        evaluation.get("classification") == "NEGATIVE"
        and evaluation.get("sign") == "NEGATIVE"
        and evaluation.get("upper_bound_is_negative") is True
        and evaluation.get("lower_bound_is_positive") is False
    )


def _canonical_confirmation(
    contract: transition.WeilTransitionCellContract,
    artifact: Mapping[str, Any],
    confirmation: dict[str, Any],
) -> dict[str, Any]:
    try:
        transition._require_payload_hash(
            confirmation, transition.WeilTransitionCheckpointError
        )
    except transition.WeilTransitionError as exc:
        raise WeilTransitionSummaryError(
            f"confirmation for {contract.contract_id} has an invalid hash"
        ) from exc

    candidate = artifact.get("quarantined_candidate")
    discoveries = [
        attempt
        for attempt in artifact["attempts"]
        if attempt.get("decision") == "NEGATIVE_WITNESS_DISCOVERED"
    ]
    if (
        artifact.get("status") != "NEGATIVE_CANDIDATE_QUARANTINED"
        or not isinstance(candidate, Mapping)
        or not isinstance(candidate.get("witness"), list)
        or len(discoveries) != 1
    ):
        raise WeilTransitionSummaryError(
            "dedicated confirmation is not attached to one v1 candidate"
        )
    witness = [str(value) for value in candidate["witness"]]
    contract_record = contract.to_record()
    base = {
        "schema": "rh-lab/weil-transition-dedicated-confirmation/v2",
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "contract_id": contract.contract_id,
        "cell_contract_sha256": contract_record["cell_contract_sha256"],
        "matrix_input_sha256": contract.matrix_input_sha256,
        "discovery_attempt_payload_sha256": discoveries[0]["payload_sha256"],
        "witness": witness,
        "precision_bits": str(contract.policy.confirmation_bits),
        "limitation": transition.EXPLORATORY_LIMITATION,
    }
    for key, value in base.items():
        if confirmation.get(key) != value:
            raise WeilTransitionSummaryError(
                f"confirmation for {contract.contract_id} changed {key}"
            )

    mode = confirmation.get("mode")
    if mode == "REUSED_V1_TERMINAL_CONFIRMATION":
        terminal = artifact["attempts"][-1]
        expected = transition._with_payload_hash(
            {
                **base,
                "mode": mode,
                "source_attempt_payload_sha256": terminal["payload_sha256"],
                "witness_evaluation": terminal["witness_evaluation"],
                "precision_containment": terminal["precision_containment"],
                "decision": "V2_NEGATIVE_CANDIDATE_CONFIRMED",
            }
        )
        if confirmation != expected:
            raise WeilTransitionSummaryError(
                "reused dedicated confirmation is not canonical"
            )
        return confirmation

    if mode != "DEDICATED_REGENERATION":
        raise WeilTransitionSummaryError("unknown dedicated-confirmation mode")
    construction = confirmation.get("matrix_construction")
    if not isinstance(construction, Mapping):
        raise WeilTransitionSummaryError("dedicated matrix-construction state is missing")
    if construction.get("classification") != "COMPLETE":
        allowed = set(base) | {
            "mode",
            "matrix_construction",
            "decision",
            "payload_sha256",
        }
        if (
            set(confirmation) != allowed
            or dict(construction)
            != {
                "classification": "INCONCLUSIVE",
                "failure_code": "SOURCE_ORACLE_OR_QUADRATURE_INCONCLUSIVE",
            }
            or confirmation.get("decision")
            != "V2_DEDICATED_CONFIRMATION_INCONCLUSIVE"
        ):
            raise WeilTransitionSummaryError(
                "incomplete dedicated confirmation is not canonical"
            )
        return confirmation

    allowed = set(base) | {
        "mode",
        "matrix_construction",
        "matrix_evidence",
        "witness_evaluation",
        "precision_containment",
        "decision",
        "payload_sha256",
    }
    if set(confirmation) != allowed:
        raise WeilTransitionSummaryError(
            "complete dedicated confirmation has noncanonical fields"
        )
    if dict(construction) != {"classification": "COMPLETE"} or not isinstance(
        confirmation.get("matrix_evidence"), Mapping
    ):
        raise WeilTransitionSummaryError(
            "complete dedicated matrix evidence is malformed"
        )
    evaluation = confirmation.get("witness_evaluation")
    if not isinstance(evaluation, Mapping) or evaluation.get("witness") != witness:
        raise WeilTransitionSummaryError("dedicated confirmation changed its witness")
    classification = evaluation.get("classification")
    sign = evaluation.get("sign")
    lower_positive = evaluation.get("lower_bound_is_positive")
    upper_negative = evaluation.get("upper_bound_is_negative")
    valid_evaluation = (
        classification == "NEGATIVE"
        and sign == "NEGATIVE"
        and lower_positive is False
        and upper_negative is True
    ) or (
        classification == "POSITIVE_FOR_WITNESS"
        and sign == "POSITIVE"
        and lower_positive is True
        and upper_negative is False
    ) or (
        classification == "INCONCLUSIVE"
        and sign == "INCONCLUSIVE"
        and lower_positive is False
        and upper_negative is False
    )
    if not valid_evaluation:
        raise WeilTransitionSummaryError(
            "dedicated integer-witness sign fields are inconsistent"
        )
    containment = confirmation.get("precision_containment")
    if (
        not isinstance(containment, Mapping)
        or set(containment)
        != {"all_higher_precision_components_contained", "first_failure"}
        or not isinstance(
            containment.get("all_higher_precision_components_contained"), bool
        )
        or (
            containment.get("all_higher_precision_components_contained") is True
            and containment.get("first_failure") is not None
        )
    ):
        raise WeilTransitionSummaryError("dedicated precision containment is missing")
    expected_decision = (
        "V2_NEGATIVE_CANDIDATE_CONFIRMED"
        if _negative_evaluation(evaluation)
        and containment.get("all_higher_precision_components_contained") is True
        else "V2_DEDICATED_CONFIRMATION_INCONCLUSIVE"
    )
    if confirmation.get("decision") != expected_decision:
        raise WeilTransitionSummaryError(
            "dedicated confirmation decision does not follow from stored evidence"
        )
    return confirmation


def _witness_evaluation_counts(
    attempts: list[Mapping[str, Any]], confirmation: Mapping[str, Any] | None
) -> Counter[str]:
    counts: Counter[str] = Counter()
    for attempt in attempts:
        if attempt.get("kind") == "SEARCH":
            generation = attempt.get("candidate_generation")
            tested = generation.get("tested", []) if isinstance(generation, Mapping) else []
            if not isinstance(tested, list):
                raise WeilTransitionSummaryError("tested witnesses must be a list")
            counts["v1_search"] += len(tested)
        elif attempt.get("kind") == "CONFIRMATION" and attempt.get(
            "witness_evaluation"
        ) is not None:
            counts["v1_confirmation"] += 1
    if (
        confirmation is not None
        and confirmation.get("mode") == "DEDICATED_REGENERATION"
        and confirmation.get("witness_evaluation") is not None
    ):
        counts["v2_dedicated_confirmation"] += 1
    return counts


def _candidate_reference(
    reference: Mapping[str, Any], artifact: Mapping[str, Any]
) -> dict[str, Any]:
    return {
        "contract_id": reference["contract_id"],
        "index_status": reference["status"],
        "v1_status": reference["v1_status"],
        "cell": {
            "path": reference["path"],
            "payload_sha256": artifact["payload_sha256"],
        },
        "dedicated_confirmation": reference["dedicated_confirmation"],
    }


def _validated_rump_enclosures(
    spectrum: Mapping[str, Any],
    contract: transition.WeilTransitionCellContract,
) -> tuple[list[tuple[Mapping[str, Any], Fraction, Fraction]], bool] | None:
    rump_failure = {
        "classification": "INCONCLUSIVE",
        "failure_code": "RUMP_EIGENVALUE_ISOLATION_FAILED",
    }
    if dict(spectrum) == rump_failure:
        return None
    if "failure_code" in spectrum:
        raise WeilTransitionSummaryError("Rump diagnostic failure record is malformed")
    required_spectrum_fields = {
        "algorithm",
        "classification",
        "eigenvalues",
        "ordering",
        "smallest_eigenvalue",
    }
    if set(spectrum) != required_spectrum_fields:
        raise WeilTransitionSummaryError("Rump diagnostic record is noncanonical")
    if spectrum.get("algorithm") != "flint-arb_mat-eig-rump" or spectrum.get(
        "ordering"
    ) != "ascending by pairwise-separated real enclosures":
        raise WeilTransitionSummaryError("Rump diagnostic metadata is noncanonical")
    classification = spectrum.get("classification")
    if classification not in {"POSITIVE_SPECTRUM", "INCONCLUSIVE"}:
        raise WeilTransitionSummaryError("Rump diagnostic classification is invalid")

    eigenvalues = spectrum.get("eigenvalues")
    expected_dimension = 2 * contract.degree + 1
    if not isinstance(eigenvalues, list) or len(eigenvalues) != expected_dimension:
        raise WeilTransitionSummaryError(
            "Rump eigenvalue list does not match the matrix dimension"
        )
    eigenvalue_fields = {
        "index",
        "real",
        "imaginary",
        "imaginary_contains_zero",
        "real_part_is_positive",
    }
    enclosures: list[tuple[Mapping[str, Any], Fraction, Fraction]] = []
    real_bounds: list[tuple[Fraction, Fraction]] = []
    all_positive = True
    for index, eigenvalue in enumerate(eigenvalues):
        if not isinstance(eigenvalue, Mapping) or set(eigenvalue) != eigenvalue_fields:
            raise WeilTransitionSummaryError("Rump eigenvalue record is noncanonical")
        if eigenvalue.get("index") != str(index):
            raise WeilTransitionSummaryError("Rump eigenvalue index is noncanonical")
        real_lower, real_upper = _dyadic_bounds(eigenvalue["real"])
        imaginary_lower, imaginary_upper = _dyadic_bounds(eigenvalue["imaginary"])
        imaginary_contains_zero = imaginary_lower <= 0 <= imaginary_upper
        real_is_positive = imaginary_contains_zero and real_lower > 0
        if eigenvalue.get("imaginary_contains_zero") is not imaginary_contains_zero:
            raise WeilTransitionSummaryError(
                "Rump imaginary-zero flag disagrees with its exact bounds"
            )
        if eigenvalue.get("real_part_is_positive") is not real_is_positive:
            raise WeilTransitionSummaryError(
                "Rump positivity flag disagrees with its exact bounds"
            )
        real_bounds.append((real_lower, real_upper))
        enclosures.append((eigenvalue, real_lower, real_upper))
        all_positive = all_positive and real_is_positive

    separated_order = all(
        left[1] < right[0]
        for left, right in zip(real_bounds, real_bounds[1:])
    )
    expected_classification = (
        "POSITIVE_SPECTRUM"
        if all_positive and separated_order
        else "INCONCLUSIVE"
    )
    if classification != expected_classification:
        raise WeilTransitionSummaryError(
            "Rump classification does not follow from its exact enclosures"
        )

    smallest = spectrum.get("smallest_eigenvalue")
    smallest_fields = {
        "ordered_index",
        "real",
        "imaginary",
        "separated_from_all_others",
    }
    if not isinstance(smallest, Mapping) or set(smallest) != smallest_fields:
        raise WeilTransitionSummaryError(
            "Rump smallest-eigenvalue record is noncanonical"
        )
    if smallest.get("ordered_index") != "0":
        raise WeilTransitionSummaryError("Rump smallest eigenvalue index changed")
    if (
        smallest.get("real") != eigenvalues[0]["real"]
        or smallest.get("imaginary") != eigenvalues[0]["imaginary"]
    ):
        raise WeilTransitionSummaryError(
            "Rump smallest eigenvalue is not bound to eigenvalue zero"
        )
    smallest_is_separated = all(
        real_bounds[0][1] < other_lower
        for other_lower, _ in real_bounds[1:]
    )
    if smallest.get("separated_from_all_others") is not smallest_is_separated:
        raise WeilTransitionSummaryError(
            "Rump smallest-eigenvalue separation flag changed"
        )
    if classification == "POSITIVE_SPECTRUM" and not (real_bounds[0][0] > 0):
        raise WeilTransitionSummaryError(
            "POSITIVE_SPECTRUM has a nonpositive smallest lower bound"
        )
    return enclosures, separated_order


def _rump_attempt_reference(
    contract: transition.WeilTransitionCellContract,
    artifact: Mapping[str, Any],
    attempt: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "contract_id": contract.contract_id,
        "cell_payload_sha256": artifact["payload_sha256"],
        "attempt_payload_sha256": attempt["payload_sha256"],
        "attempt_sequence": attempt["sequence"],
        "precision_bits": attempt["precision_bits"],
    }


def _consider_rump(
    current: tuple[tuple[Any, ...], dict[str, Any]] | None,
    *,
    contract: transition.WeilTransitionCellContract,
    artifact: Mapping[str, Any],
    attempt: Mapping[str, Any],
) -> tuple[
    tuple[tuple[Any, ...], dict[str, Any]] | None,
    dict[str, Any],
]:
    if attempt.get("kind") != "SEARCH":
        raise WeilTransitionSummaryError("Rump diagnostics require a search attempt")
    reference = _rump_attempt_reference(contract, artifact, attempt)
    if "secondary_rump_spectrum" not in attempt:
        if attempt.get("matrix_construction", {}).get("classification") == "COMPLETE":
            raise WeilTransitionSummaryError(
                "complete search attempt lacks a Rump diagnostic record"
            )
        return current, {
            "kind": "MATRIX_CONSTRUCTION_INCONCLUSIVE_WITHOUT_SPECTRUM",
            "reference": reference,
        }
    spectrum = attempt["secondary_rump_spectrum"]
    if not isinstance(spectrum, Mapping):
        raise WeilTransitionSummaryError("Rump diagnostic record must be an object")
    validated = _validated_rump_enclosures(spectrum, contract)
    if validated is None:
        return current, {
            "kind": "FAILURE_SENTINEL",
            "reference": {
                **reference,
                "failure_code": "RUMP_EIGENVALUE_ISOLATION_FAILED",
            },
        }

    enclosures, pairwise_separated = validated
    sequence = _canonical_integer(attempt.get("sequence"), "attempt sequence")
    for eigenvalue, lower, upper in enclosures:
        eigenvalue_index = _canonical_integer(
            eigenvalue.get("index"), "Rump eigenvalue index"
        )
        key = (
            lower,
            upper,
            contract.sort_key,
            sequence,
            eigenvalue_index,
        )
        record = {
            **reference,
            "stored_eigenvalue_index": eigenvalue["index"],
            "spectrum_classification": spectrum.get("classification"),
            "spectrum_pairwise_separated_in_stored_order": pairwise_separated,
            "real_enclosure": _normalized_enclosure(eigenvalue["real"]),
        }
        if current is None or key < current[0]:
            current = key, record
    return current, {
        "kind": "USABLE_SPECTRUM",
        "eigenvalue_enclosures": len(enclosures),
        "pairwise_separated": pairwise_separated,
        "reference": {
            **reference,
            "spectrum_classification": spectrum.get("classification"),
        },
    }


def _count_record(counter: Counter[str]) -> dict[str, str]:
    def key(item: tuple[str, int]) -> tuple[int, int | str]:
        label = item[0]
        try:
            return 0, int(label)
        except ValueError:
            return 1, label

    return {label: str(count) for label, count in sorted(counter.items(), key=key)}


def generate_weil_transition_summary(checkpoint_dir: Path) -> dict[str, Any]:
    """Generate a compact summary of a complete frozen 81-cell v2 run.

    Stored terminal attempts are semantically derived, but no numerical
    evidence is rebuilt.  Use :func:`verify_weil_transition_search` when a
    numerical replay is required.
    """

    root = checkpoint_dir.resolve()
    plan_path = root / "plan.json"
    run_path = root / "run.json"
    index_path = root / "index.json"

    raw_plan = _load_object(plan_path)
    try:
        plan = transition.canonicalize_transition_plan(raw_plan)
    except transition.WeilTransitionPlanError as exc:
        raise WeilTransitionSummaryError("stored transition plan is invalid") from exc
    if raw_plan != plan:
        raise WeilTransitionSummaryError("stored transition plan is not canonical")
    if (
        plan.get("frozen_batch") != transition.FROZEN_Q7_Q9_BATCH
        or len(plan["cells"]) != 81
    ):
        raise WeilTransitionSummaryError("summary requires the frozen 81-cell batch")
    if plan.get("plan_sha256") != FROZEN_PLAN_SHA256:
        raise WeilTransitionSummaryError("frozen transition plan hash changed")
    runtime_backend = v1._backend_record()
    if plan.get("backend_contract") != runtime_backend:
        raise WeilTransitionSummaryError(
            "frozen plan backend does not match the current backend"
        )

    run = _load_object(run_path)
    try:
        transition._require_payload_hash(run, transition.WeilTransitionCheckpointError)
    except transition.WeilTransitionError as exc:
        raise WeilTransitionSummaryError("run manifest hash is invalid") from exc
    if run != transition._run_manifest(plan):
        raise WeilTransitionSummaryError("run manifest is not canonical for the plan")
    if run.get("backend") != runtime_backend:
        raise WeilTransitionSummaryError("run backend binding changed")

    index = _load_object(index_path)
    try:
        transition._require_payload_hash(index, transition.WeilTransitionCheckpointError)
    except transition.WeilTransitionError as exc:
        raise WeilTransitionSummaryError("transition index hash is invalid") from exc
    if index.get("schema") != transition.INDEX_SCHEMA:
        raise WeilTransitionSummaryError("unexpected transition index schema")
    if index.get("classification") != "EXPLORATORY":
        raise WeilTransitionSummaryError("transition index classification changed")
    if index.get("hypothesis_status") != "UNRESOLVED":
        raise WeilTransitionSummaryError("transition index changed hypothesis status")
    if index.get("backend") != runtime_backend:
        raise WeilTransitionSummaryError("index backend binding changed")

    contracts = transition._contracts_from_plan(plan)
    by_id = {contract.contract_id: contract for contract in contracts}
    references = index.get("cells")
    if not isinstance(references, list) or len(references) != 81:
        raise WeilTransitionSummaryError("summary requires exactly 81 index references")
    if index.get("missing_contract_ids") != []:
        raise WeilTransitionSummaryError("summary requires a complete transition index")

    status_counts: Counter[str] = Counter()
    q_counts: Counter[str] = Counter()
    position_counts: Counter[str] = Counter()
    degree_counts: Counter[str] = Counter()
    precision_counts: Counter[str] = Counter()
    evaluation_counts: Counter[str] = Counter()
    total_attempts = 0
    candidate_references: list[dict[str, Any]] = []
    conflict_references: list[dict[str, Any]] = []
    rump_minimum: tuple[tuple[Any, ...], dict[str, Any]] | None = None
    rump_coverage: Counter[str] = Counter()
    rump_failure_references: list[dict[str, Any]] = []
    rump_unseparated_references: list[dict[str, Any]] = []
    rump_construction_failure_references: list[dict[str, Any]] = []
    completed: dict[
        str, tuple[dict[str, Any], dict[str, Any] | None]
    ] = {}
    seen: set[str] = set()

    for reference in references:
        if not isinstance(reference, Mapping):
            raise WeilTransitionSummaryError("index cell reference must be an object")
        contract_id = reference.get("contract_id")
        if not isinstance(contract_id, str) or contract_id not in by_id or contract_id in seen:
            raise WeilTransitionSummaryError("unknown or duplicate cell reference")
        seen.add(contract_id)
        contract = by_id[contract_id]
        canonical_contract = contract.to_record()
        if reference.get("cell_contract_sha256") != canonical_contract[
            "cell_contract_sha256"
        ]:
            raise WeilTransitionSummaryError("index cell-contract hash changed")
        if reference.get("matrix_input_sha256") != contract.matrix_input_sha256:
            raise WeilTransitionSummaryError("index matrix-input hash changed")

        artifact = _canonical_cell(contract, _load_object(_safe_path(root, reference.get("path"))))
        if reference.get("v1_cell_payload_sha256") != artifact["payload_sha256"]:
            raise WeilTransitionSummaryError("index v1 cell hash changed")
        if reference.get("v1_status") != artifact["status"]:
            raise WeilTransitionSummaryError("index v1 status changed")

        confirmation_reference = reference.get("dedicated_confirmation")
        if "dedicated_confirmation" not in reference:
            raise WeilTransitionSummaryError(
                "index cell reference is missing dedicated_confirmation"
            )
        confirmation: dict[str, Any] | None = None
        if confirmation_reference is not None:
            if not isinstance(confirmation_reference, Mapping):
                raise WeilTransitionSummaryError("confirmation reference must be an object")
            confirmation = _canonical_confirmation(
                contract,
                artifact,
                _load_object(_safe_path(root, confirmation_reference.get("path"))),
            )
            if confirmation_reference.get("payload_sha256") != confirmation[
                "payload_sha256"
            ] or confirmation_reference.get("decision") != confirmation["decision"]:
                raise WeilTransitionSummaryError("confirmation reference changed")
        elif artifact["status"] == "NEGATIVE_CANDIDATE_QUARANTINED":
            raise WeilTransitionSummaryError("v1 candidate lacks dedicated confirmation")

        try:
            status = transition._v2_status(artifact, confirmation)
        except transition.WeilTransitionError as exc:
            raise WeilTransitionSummaryError("v2 status cannot be derived") from exc
        if reference.get("status") != status:
            raise WeilTransitionSummaryError("index v2 status changed")

        attempts = _terminal_attempts(artifact)
        terminal_precision = (
            confirmation["precision_bits"] if confirmation is not None else attempts[-1]["precision_bits"]
        )
        _canonical_integer(terminal_precision, "terminal precision")
        total_attempts += len(attempts)
        evaluation_counts.update(_witness_evaluation_counts(attempts, confirmation))
        status_counts[status] += 1
        q_counts[str(contract.q)] += 1
        position_counts[contract.position] += 1
        degree_counts[str(contract.degree)] += 1
        precision_counts[terminal_precision] += 1

        compact_reference = _candidate_reference(reference, artifact)
        if status == "NEGATIVE_CANDIDATE_QUARANTINED":
            candidate_references.append(compact_reference)
        if (
            status == "NEGATIVE_OBSERVATION_QUARANTINED_INCONCLUSIVE"
            or (
                isinstance(artifact.get("quarantined_candidate"), Mapping)
                and artifact["quarantined_candidate"].get("status")
                == "CONFLICTING_NEGATIVE_REPLAY_QUARANTINED"
            )
        ):
            conflict_references.append(compact_reference)

        for attempt in attempts:
            if attempt.get("kind") != "SEARCH":
                continue
            rump_coverage["search_attempts"] += 1
            rump_minimum, rump_observation = _consider_rump(
                rump_minimum,
                contract=contract,
                artifact=artifact,
                attempt=attempt,
            )
            if rump_observation["kind"] == "USABLE_SPECTRUM":
                rump_coverage["usable_spectra"] += 1
                rump_coverage["usable_eigenvalue_enclosures"] += rump_observation[
                    "eigenvalue_enclosures"
                ]
                if rump_observation["pairwise_separated"]:
                    rump_coverage["pairwise_separated_usable_spectra"] += 1
                else:
                    rump_coverage["unseparated_usable_spectra"] += 1
                    rump_unseparated_references.append(rump_observation["reference"])
            elif rump_observation["kind"] == "FAILURE_SENTINEL":
                rump_coverage["failure_sentinels"] += 1
                rump_failure_references.append(rump_observation["reference"])
            else:
                rump_coverage[
                    "matrix_construction_inconclusive_without_spectrum"
                ] += 1
                rump_construction_failure_references.append(
                    rump_observation["reference"]
                )
        completed[contract_id] = (
            {
                "status": artifact["status"],
                "payload_sha256": artifact["payload_sha256"],
            },
            None
            if confirmation is None
            else {
                "decision": confirmation["decision"],
                "payload_sha256": confirmation["payload_sha256"],
            },
        )

    if seen != set(by_id):
        raise WeilTransitionSummaryError("complete index does not cover the frozen plan")
    expected_index = transition._index_payload(plan, contracts, completed)
    if index != expected_index:
        raise WeilTransitionSummaryError("transition index is not canonical")
    if rump_minimum is None:
        raise WeilTransitionSummaryError("completed batch has no Rump diagnostic interval")

    for status in (
        "FINITE_POSITIVE_CERTIFIED",
        "NEGATIVE_CANDIDATE_QUARANTINED",
        "NEGATIVE_OBSERVATION_QUARANTINED_INCONCLUSIVE",
        "INCONCLUSIVE_MAX_PRECISION",
    ):
        status_counts.setdefault(status, 0)

    total_evaluations = sum(evaluation_counts.values())
    body = {
        "schema": SUMMARY_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "conclusion": index["conclusion"],
        "source": {
            "plan": {"path": "plan.json", "plan_sha256": plan["plan_sha256"]},
            "run": {"path": "run.json", "payload_sha256": run["payload_sha256"]},
            "index": {
                "path": "index.json",
                "payload_sha256": index["payload_sha256"],
            },
        },
        "frozen_batch": transition.FROZEN_Q7_Q9_BATCH,
        "counts": {
            "cells": "81",
            "by_status": _count_record(status_counts),
            "by_q": _count_record(q_counts),
            "by_position": _count_record(position_counts),
            "by_degree": _count_record(degree_counts),
            "by_terminal_precision_bits": _count_record(precision_counts),
            "terminal_precision_definition": (
                "dedicated v2 confirmation precision when present; otherwise "
                "the final v1 attempt precision"
            ),
            "total_v1_attempts": str(total_attempts),
            "stored_candidate_evaluation_records": {
                "definition": (
                    "hash-bound stored records counted from candidate-generation "
                    "tested arrays or confirmation fields; this summary neither "
                    "replays their numerical values nor validates that recorded "
                    "witness entries are canonical integers"
                ),
                "v1_search": str(evaluation_counts["v1_search"]),
                "v1_confirmation": str(evaluation_counts["v1_confirmation"]),
                "v2_dedicated_confirmation": str(
                    evaluation_counts["v2_dedicated_confirmation"]
                ),
                "total": str(total_evaluations),
            },
        },
        "quarantine_references": {
            "negative_candidates": candidate_references,
            "conflicting_negative_observations": conflict_references,
        },
        "rump_diagnostic": {
            "role": "DIAGNOSTIC_ONLY",
            "scope": (
                "all stored Rump eigenvalue enclosures in usable spectra from "
                "SEARCH attempts; this is not an ordered-eigenvalue or RH claim"
            ),
            "coverage": {
                "search_attempts": str(rump_coverage["search_attempts"]),
                "usable_spectra": str(rump_coverage["usable_spectra"]),
                "usable_eigenvalue_enclosures": str(
                    rump_coverage["usable_eigenvalue_enclosures"]
                ),
                "pairwise_separated_usable_spectra": str(
                    rump_coverage["pairwise_separated_usable_spectra"]
                ),
                "unseparated_usable_spectra": str(
                    rump_coverage["unseparated_usable_spectra"]
                ),
                "failure_sentinels": str(rump_coverage["failure_sentinels"]),
                "matrix_construction_inconclusive_without_spectrum": str(
                    rump_coverage[
                        "matrix_construction_inconclusive_without_spectrum"
                    ]
                ),
                "failure_references": rump_failure_references,
                "unseparated_references": rump_unseparated_references,
                "matrix_construction_inconclusive_references": (
                    rump_construction_failure_references
                ),
            },
            "least_lower_endpoint_among_usable_stored_enclosures": {
                "selection_rule": (
                    "least exact lower endpoint, then least exact upper endpoint, "
                    "then canonical cell, attempt, and stored eigenvalue order"
                ),
                **rump_minimum[1],
            },
        },
        "summary_method": {
            "numerical_replay_performed": False,
            "stored_terminal_state_semantically_derived": True,
            "content_hashes_checked": [
                "plan_sha256",
                "run.payload_sha256",
                "index.payload_sha256",
                "cell.payload_sha256",
                "attempt.payload_sha256",
                "dedicated_confirmation.payload_sha256_when_present",
            ],
        },
        "limitation": SUMMARY_LIMITATION,
    }
    return {**body, "payload_sha256": content_sha256(body)}


def verify_weil_transition_summary(
    summary: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
) -> dict[str, Any]:
    """Regenerate a compact summary and reject any supplied difference."""

    if isinstance(summary, Path):
        supplied = _load_object(summary)
    elif isinstance(summary, Mapping):
        supplied = dict(summary)
    else:
        raise TypeError("summary must be a mapping or Path")
    body = {key: value for key, value in supplied.items() if key != "payload_sha256"}
    if supplied.get("payload_sha256") != content_sha256(body):
        raise WeilTransitionSummaryVerificationError("summary payload hash mismatch")
    expected = generate_weil_transition_summary(checkpoint_dir)
    if supplied != expected:
        raise WeilTransitionSummaryVerificationError(
            "summary does not canonically regenerate from the checkpoint"
        )
    return {
        "classification": "REPRODUCED_EXPLORATORY_SUMMARY",
        "hypothesis_status": "UNRESOLVED",
        "summary_payload_sha256": supplied["payload_sha256"],
        "plan_sha256": supplied["source"]["plan"]["plan_sha256"],
        "index_payload_sha256": supplied["source"]["index"]["payload_sha256"],
        "verified_cells": supplied["counts"]["cells"],
        "numerical_replay_performed": False,
    }


__all__ = [
    "FROZEN_PLAN_SHA256",
    "SUMMARY_LIMITATION",
    "SUMMARY_SCHEMA",
    "WeilTransitionSummaryError",
    "WeilTransitionSummaryVerificationError",
    "generate_weil_transition_summary",
    "verify_weil_transition_summary",
]
