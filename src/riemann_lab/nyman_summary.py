"""Compact deterministic summaries for frozen natural-distance evidence.

This module validates the canonical checkpoint structure and the exact hashes
of the stored finite-dimensional Nyman evidence.  It deliberately does not
rebuild an Arb kernel, solve a Gram system, or replay either certificate.  The
result is therefore an integrity and reporting layer, not an independent
numerical reproduction.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from .artifacts import content_sha256
from . import nyman as core
from . import nyman_search as search


SUMMARY_SCHEMA = "rh-lab/nyman-natural-summary/v1"
FROZEN_PLAN_SHA256 = (
    "27ebd77ca7bfae0ff898bd33fcd739c1c069a79b23b9bff4f03f84ac91f1540a"
)
SUMMARY_LIMITATION = (
    f"{search.EXPLORATORY_LIMITATION} This compact summary checks canonical "
    "checkpoint records and their content hashes without rebuilding the "
    "768-bit generation kernel or rerunning the 1536-bit numerical verifier. "
    "Stored energy displays are reporting metadata, not a new certificate."
)

_HASH_PATTERN = re.compile(r"[0-9a-f]{64}")
_KERNEL_CONSTRUCTION = (
    "one canonical full max-N kernel; each listed commitment is hashed from "
    "an exact leading in-memory prefix"
)
_PROVENANCE_CONSTRUCTION = (
    "canonical-rational-vasyunin-with-reciprocity-intersection"
)


class NymanSummaryError(ValueError):
    """Raised when a Nyman checkpoint cannot be summarized safely."""


class NymanSummaryVerificationError(NymanSummaryError):
    """Raised when a supplied summary does not canonically regenerate."""


def _load_object(
    path: Path,
    error_type: type[NymanSummaryError] = NymanSummaryError,
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
        raise error_type(f"cannot read Nyman checkpoint object: {path}") from exc
    if not isinstance(value, dict):
        raise error_type("Nyman checkpoint object must be a JSON object")
    return value


def _safe_path(root: Path, relative: Any) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise NymanSummaryError("invalid Nyman checkpoint reference path")
    candidate = Path(relative)
    if candidate.is_absolute() or any(
        part in {"", ".", ".."} for part in candidate.parts
    ):
        raise NymanSummaryError("invalid Nyman checkpoint reference path")
    resolved = (root / candidate).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise NymanSummaryError("Nyman checkpoint reference escapes its root") from exc
    return resolved


def _require_fields(record: Any, expected: set[str], name: str) -> Mapping[str, Any]:
    if not isinstance(record, Mapping) or set(record) != expected:
        raise NymanSummaryError(f"{name} has noncanonical fields")
    return record


def _require_hash(value: Any, name: str) -> str:
    if not isinstance(value, str) or _HASH_PATTERN.fullmatch(value) is None:
        raise NymanSummaryError(f"{name} must be a lowercase SHA-256 digest")
    return value


def _require_payload_hash(record: Mapping[str, Any], name: str) -> None:
    supplied = _require_hash(record.get("payload_sha256"), f"{name} payload hash")
    body = {key: value for key, value in record.items() if key != "payload_sha256"}
    if supplied != content_sha256(body):
        raise NymanSummaryError(f"{name} payload hash mismatch")


def _canonical_integer(value: Any, name: str) -> int:
    if not isinstance(value, str):
        raise NymanSummaryError(f"{name} must be a canonical decimal string")
    try:
        parsed = int(value)
    except ValueError as exc:
        raise NymanSummaryError(f"{name} is not an integer") from exc
    if value != str(parsed):
        raise NymanSummaryError(f"{name} is not canonical")
    return parsed


def _dyadic_value(mantissa: Any, exponent: Any, name: str) -> Fraction:
    coefficient = _canonical_integer(mantissa, f"{name} mantissa")
    power = _canonical_integer(exponent, f"{name} exponent")
    if power >= 0:
        return Fraction(coefficient * (1 << power), 1)
    return Fraction(coefficient, 1 << -power)


def _arb_bounds(record: Any, name: str) -> tuple[Fraction, Fraction]:
    value = _require_fields(record, {"display", "dyadic", "is_exact"}, name)
    if not isinstance(value["display"], str) or not value["display"]:
        raise NymanSummaryError(f"{name} display must be a nonempty string")
    if not isinstance(value["is_exact"], bool):
        raise NymanSummaryError(f"{name} exactness flag must be Boolean")
    encoded = _require_fields(
        value["dyadic"],
        {
            "mid_mantissa",
            "mid_exponent",
            "radius_mantissa",
            "radius_exponent",
        },
        f"{name} dyadic enclosure",
    )
    midpoint = _dyadic_value(
        encoded["mid_mantissa"], encoded["mid_exponent"], f"{name} midpoint"
    )
    radius = _dyadic_value(
        encoded["radius_mantissa"],
        encoded["radius_exponent"],
        f"{name} radius",
    )
    if radius < 0:
        raise NymanSummaryError(f"{name} radius must be nonnegative")
    if value["is_exact"] is not (radius == 0):
        raise NymanSummaryError(f"{name} exactness flag disagrees with its radius")
    return midpoint - radius, midpoint + radius


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _bound_record(raw: Mapping[str, Any], value: Fraction) -> dict[str, Any]:
    return {
        "dyadic": dict(raw),
        "exact_fraction": _fraction_record(value),
    }


def _canonical_dyadic_record(value: Fraction) -> dict[str, str]:
    denominator = value.denominator
    if denominator <= 0 or denominator & (denominator - 1):
        raise NymanSummaryError("derived bracket width is not dyadic")
    return {
        "numerator": str(value.numerator),
        "denominator_exponent": str(denominator.bit_length() - 1),
    }


def _validate_kernel(
    kernel: dict[str, Any],
    plan: Mapping[str, Any],
    contracts: Sequence[search.NymanSearchCellContract],
) -> dict[str, Mapping[str, Any]]:
    _require_fields(
        kernel,
        {
            "schema",
            "classification",
            "hypothesis_status",
            "precision_bits",
            "max_n",
            "kernel_contract_sha256",
            "core_system_content_sha256",
            "construction",
            "prefixes",
            "payload_sha256",
        },
        "generation kernel",
    )
    _require_payload_hash(kernel, "generation kernel")
    generation_bits = search._policy_values(plan)[0]
    max_n = max(contract.shared_kernel_n for contract in contracts)
    expected_contract_hash = plan["engine"]["shared_kernel"][
        "kernel_contract_sha256"
    ]
    expected_metadata = {
        "schema": search.KERNEL_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "precision_bits": str(generation_bits),
        "max_n": str(max_n),
        "kernel_contract_sha256": expected_contract_hash,
        "construction": _KERNEL_CONSTRUCTION,
    }
    for key, expected in expected_metadata.items():
        if kernel.get(key) != expected:
            raise NymanSummaryError(f"generation kernel changed {key}")
    _require_hash(kernel.get("core_system_content_sha256"), "kernel core-system hash")

    prefixes = kernel.get("prefixes")
    if not isinstance(prefixes, list) or len(prefixes) != len(contracts):
        raise NymanSummaryError("generation kernel prefix list changed")
    by_n: dict[str, Mapping[str, Any]] = {}
    for prefix, contract in zip(prefixes, contracts, strict=True):
        item = _require_fields(
            prefix,
            {"n", "core_system_content_sha256", "prefix_kernel_sha256"},
            "generation kernel prefix",
        )
        if item.get("n") != str(contract.n):
            raise NymanSummaryError("generation kernel prefix order changed")
        _require_hash(
            item.get("core_system_content_sha256"),
            "prefix core-system hash",
        )
        _require_hash(item.get("prefix_kernel_sha256"), "prefix kernel hash")
        by_n[str(contract.n)] = item
    if by_n[str(max_n)]["core_system_content_sha256"] != kernel[
        "core_system_content_sha256"
    ]:
        raise NymanSummaryError("max-N prefix and full kernel hashes disagree")
    return by_n


def _scope(contract: search.NymanSearchCellContract) -> dict[str, Any]:
    return {
        "dilates": [str(value) for value in range(1, contract.n + 1)],
        "dimension": str(contract.n),
        "finite_only": True,
    }


def _kernel_provenance(
    contract: search.NymanSearchCellContract,
    generation_bits: int,
    prefix: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "construction": _PROVENANCE_CONSTRUCTION,
        "normalization": core.NYMAN_NORMALIZATION,
        "system_content_sha256": prefix["core_system_content_sha256"],
        "precision_bits": str(generation_bits),
        "source_dilates": [
            str(value) for value in range(1, contract.shared_kernel_n + 1)
        ],
        "source_dimension": str(contract.shared_kernel_n),
        "scope_is_leading_principal_prefix": True,
        "supplied_system_reused": True,
    }


def _validate_upper_certificate(
    record: Any,
    contract: search.NymanSearchCellContract,
    candidate: Mapping[str, Any],
    upper: Fraction,
    generation_bits: int,
    prefix: Mapping[str, Any],
) -> tuple[str, str]:
    certificate = _require_fields(
        record,
        {
            "schema",
            "classification",
            "hypothesis_status",
            "precision_bits",
            "normalization",
            "limitation",
            "decision",
            "scope",
            "kernel_provenance",
            "coefficients",
            "claimed_upper_bound",
            "energy",
            "margin",
            "checks",
            "payload_sha256",
        },
        "upper certificate",
    )
    _require_payload_hash(certificate, "upper certificate")
    expected_preamble = {
        "schema": core.NYMAN_UPPER_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "precision_bits": str(generation_bits),
        "normalization": core.NYMAN_NORMALIZATION,
        "limitation": core.NYMAN_LIMITATION,
        "scope": _scope(contract),
        "kernel_provenance": _kernel_provenance(
            contract, generation_bits, prefix
        ),
        "claimed_upper_bound": _fraction_record(upper),
    }
    for key, expected in expected_preamble.items():
        if certificate.get(key) != expected:
            raise NymanSummaryError(f"upper certificate changed {key}")
    expected_coefficients = {
        "encoding": "signed-integers-over-common-power-of-two",
        **candidate["coefficients"],
    }
    if certificate.get("coefficients") != expected_coefficients:
        raise NymanSummaryError("upper certificate changed its exact coefficients")

    _arb_bounds(certificate["energy"], "stored energy")
    margin_lower, _ = _arb_bounds(certificate["margin"], "stored upper margin")
    certified = margin_lower >= 0
    strict = margin_lower > 0
    checks = _require_fields(
        certificate["checks"],
        {
            "coefficients_are_exact_dyadics",
            "energy_at_most_claimed_upper_bound",
            "energy_strictly_below_claimed_upper_bound",
            "direct_energy_evaluation_only",
            "approximate_solve_used_as_evidence",
        },
        "upper certificate checks",
    )
    expected_checks = {
        "coefficients_are_exact_dyadics": True,
        "energy_at_most_claimed_upper_bound": certified,
        "energy_strictly_below_claimed_upper_bound": strict,
        "direct_energy_evaluation_only": True,
        "approximate_solve_used_as_evidence": False,
    }
    if dict(checks) != expected_checks:
        raise NymanSummaryError("upper certificate decision checks are inconsistent")
    expected_decision = "UPPER_BOUND_CERTIFIED" if certified else "INCONCLUSIVE"
    if certificate.get("decision") != expected_decision:
        raise NymanSummaryError("upper certificate decision is inconsistent")
    return expected_decision, str(certificate["energy"]["display"])


def _validate_ldlt(record: Any, dimension: int) -> tuple[str, int]:
    ldlt = _require_fields(
        record,
        {"classification", "fixed_order", "failed_pivot_index", "pivots"},
        "stored augmented LDLT",
    )
    if ldlt.get("fixed_order") != [str(index) for index in range(dimension)]:
        raise NymanSummaryError("stored augmented LDLT order changed")
    classification = ldlt.get("classification")
    if classification not in {"POSITIVE_DEFINITE", "INCONCLUSIVE"}:
        raise NymanSummaryError("stored augmented LDLT classification changed")
    pivots = ldlt.get("pivots")
    if not isinstance(pivots, list) or not pivots:
        raise NymanSummaryError("stored augmented LDLT pivots are missing")
    lower_endpoints: list[Fraction] = []
    for index, pivot in enumerate(pivots):
        item = _require_fields(pivot, {"index", "value"}, "stored LDLT pivot")
        if item.get("index") != str(index):
            raise NymanSummaryError("stored LDLT pivot index changed")
        lower, _ = _arb_bounds(item["value"], "stored LDLT pivot")
        lower_endpoints.append(lower)
    failed = ldlt.get("failed_pivot_index")
    if classification == "POSITIVE_DEFINITE":
        if failed is not None or len(pivots) != dimension or not all(
            value > 0 for value in lower_endpoints
        ):
            raise NymanSummaryError("positive-definite LDLT record is inconsistent")
    else:
        failed_index = _canonical_integer(failed, "failed LDLT pivot index")
        if (
            failed_index < 0
            or failed_index >= dimension
            or len(pivots) != failed_index + 1
            or not all(value > 0 for value in lower_endpoints[:failed_index])
            or lower_endpoints[failed_index] > 0
        ):
            raise NymanSummaryError("inconclusive LDLT record is inconsistent")
    return str(classification), len(pivots)


def _validate_lower_certificate(
    record: Any,
    contract: search.NymanSearchCellContract,
    lower: Fraction,
    generation_bits: int,
    prefix: Mapping[str, Any],
) -> tuple[str, int]:
    certificate = _require_fields(
        record,
        {
            "schema",
            "classification",
            "hypothesis_status",
            "precision_bits",
            "normalization",
            "limitation",
            "decision",
            "scope",
            "kernel_provenance",
            "claimed_lower_bound",
            "augmented_matrix",
            "augmented_ldlt",
            "checks",
            "payload_sha256",
        },
        "lower certificate",
    )
    _require_payload_hash(certificate, "lower certificate")
    expected_preamble = {
        "schema": core.NYMAN_LOWER_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "precision_bits": str(generation_bits),
        "normalization": core.NYMAN_NORMALIZATION,
        "limitation": core.NYMAN_LIMITATION,
        "scope": _scope(contract),
        "kernel_provenance": _kernel_provenance(
            contract, generation_bits, prefix
        ),
        "claimed_lower_bound": _fraction_record(lower),
        "augmented_matrix": {
            "block_form": "[[G,-b],[-b^T,1-L]]",
            "quadratic_identity": "[c,1]^T M [c,1]=E(c)-L",
            "dimension": str(contract.n + 1),
        },
    }
    for key, expected in expected_preamble.items():
        if certificate.get(key) != expected:
            raise NymanSummaryError(f"lower certificate changed {key}")

    classification, pivot_count = _validate_ldlt(
        certificate["augmented_ldlt"], contract.n + 1
    )
    certified = classification == "POSITIVE_DEFINITE"
    checks = _require_fields(
        certificate["checks"],
        {
            "fixed_order_interval_ldlt_positive",
            "approximate_solve_used_as_evidence",
        },
        "lower certificate checks",
    )
    if dict(checks) != {
        "fixed_order_interval_ldlt_positive": certified,
        "approximate_solve_used_as_evidence": False,
    }:
        raise NymanSummaryError("lower certificate decision checks are inconsistent")
    expected_decision = "LOWER_BOUND_CERTIFIED" if certified else "INCONCLUSIVE"
    if certificate.get("decision") != expected_decision:
        raise NymanSummaryError("lower certificate decision is inconsistent")
    return expected_decision, pivot_count


def _validate_cell(
    artifact: dict[str, Any],
    contract: search.NymanSearchCellContract,
    plan: Mapping[str, Any],
    kernel: Mapping[str, Any],
    prefix: Mapping[str, Any],
) -> tuple[dict[str, Any], Fraction, Fraction, str, int]:
    _require_fields(
        artifact,
        {
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
            "payload_sha256",
        },
        "Nyman cell",
    )
    _require_payload_hash(artifact, "Nyman cell")
    generation_bits, replay_bits, _, _, _ = search._policy_values(plan)
    expected_shell = {
        "schema": search.CELL_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "plan_sha256": plan["plan_sha256"],
        "cell_contract": contract.to_record(),
        "replay_precision_bits": str(replay_bits),
        "limitation": search.EXPLORATORY_LIMITATION,
    }
    for key, expected in expected_shell.items():
        if artifact.get(key) != expected:
            raise NymanSummaryError(f"Nyman cell changed {key}")

    try:
        _, upper, lower = search._parse_candidate(
            artifact["candidate"],
            contract,
            search._policy_values(plan),
            error_type=search.NymanSearchCheckpointError,
        )
    except search.NymanSearchError as exc:
        raise NymanSummaryError(f"cell {contract.cell_id} candidate is invalid") from exc
    if not (0 <= lower < upper <= 1):
        raise NymanSummaryError("stored finite-distance bracket is out of range")

    generation = _require_fields(
        artifact["generation"],
        {
            "precision_bits",
            "solver_role",
            "shared_kernel_binding",
            "upper_certificate",
            "lower_certificate",
        },
        "Nyman generation evidence",
    )
    if generation.get("precision_bits") != str(generation_bits):
        raise NymanSummaryError("Nyman generation precision changed")
    if generation.get("solver_role") != "approximate-untrusted-candidate-generator-only":
        raise NymanSummaryError("approximate solve was given a certificate role")
    try:
        expected_binding = search._kernel_cell_binding(kernel, contract)
    except search.NymanSearchError as exc:
        raise NymanSummaryError("shared kernel binding cannot be derived") from exc
    if generation.get("shared_kernel_binding") != expected_binding:
        raise NymanSummaryError("Nyman shared-kernel prefix binding changed")

    upper_decision, energy_display = _validate_upper_certificate(
        generation["upper_certificate"],
        contract,
        artifact["candidate"],
        upper,
        generation_bits,
        prefix,
    )
    lower_decision, pivot_count = _validate_lower_certificate(
        generation["lower_certificate"],
        contract,
        lower,
        generation_bits,
        prefix,
    )
    expected_status = (
        "FINITE_DISTANCE_BRACKET_CERTIFIED"
        if upper_decision == "UPPER_BOUND_CERTIFIED"
        and lower_decision == "LOWER_BOUND_CERTIFIED"
        else "INCONCLUSIVE"
    )
    expected_statement = (
        "L_N < d_N^2 <= U_N"
        if expected_status == "FINITE_DISTANCE_BRACKET_CERTIFIED"
        else None
    )
    if artifact.get("status") != expected_status:
        raise NymanSummaryError("Nyman cell status does not follow from stored evidence")
    if artifact.get("terminal_statement") != expected_statement:
        raise NymanSummaryError("Nyman cell terminal statement changed")
    return artifact, lower, upper, energy_display, pivot_count


def _count_record(counter: Counter[str]) -> dict[str, str]:
    return {key: str(counter[key]) for key in sorted(counter)}


def generate_nyman_summary(checkpoint_dir: Path) -> dict[str, Any]:
    """Generate a compact structural summary of the complete frozen Nyman run.

    No natural system is rebuilt and no numerical certificate is replayed.
    Use :func:`riemann_lab.nyman_search.verify_nyman_search` for that stronger,
    same-backend numerical reproduction.
    """

    root = checkpoint_dir.resolve()
    raw_plan = _load_object(root / "plan.json")
    try:
        plan = search.canonicalize_nyman_plan(raw_plan)
    except search.NymanSearchPlanError as exc:
        raise NymanSummaryError("stored Nyman plan is invalid") from exc
    if raw_plan != plan:
        raise NymanSummaryError("stored Nyman plan is not canonical")
    if (
        plan.get("frozen_experiment") != search.FROZEN_EXPERIMENT
        or tuple(int(cell["n"]) for cell in plan["cells"])
        != search.FROZEN_N_VALUES
    ):
        raise NymanSummaryError("summary requires the frozen six-cell experiment")
    if plan.get("plan_sha256") != FROZEN_PLAN_SHA256:
        raise NymanSummaryError("frozen Nyman plan hash changed")
    runtime_backend = search._backend_record()
    if plan.get("backend_contract") != runtime_backend:
        raise NymanSummaryError("frozen plan backend does not match the current backend")
    contracts = search._contracts(plan)

    run = _load_object(root / "run.json")
    _require_fields(
        run,
        {
            "schema",
            "classification",
            "hypothesis_status",
            "plan",
            "backend",
            "engine_id",
            "candidate_behavior",
            "limitation",
            "payload_sha256",
        },
        "Nyman run manifest",
    )
    _require_payload_hash(run, "Nyman run manifest")
    if run != search._run_manifest(plan):
        raise NymanSummaryError("Nyman run manifest is not canonical for the plan")

    index = _load_object(root / "index.json")
    _require_fields(
        index,
        {
            "schema",
            "classification",
            "hypothesis_status",
            "conclusion",
            "plan",
            "run",
            "backend",
            "generation_kernel",
            "cells",
            "missing_cell_ids",
            "progress",
            "strict_improvement_diagnostics",
            "strict_improvement_role",
            "limitation",
            "payload_sha256",
        },
        "Nyman index",
    )
    _require_payload_hash(index, "Nyman index")
    if index.get("schema") != search.INDEX_SCHEMA:
        raise NymanSummaryError("unexpected Nyman index schema")
    if index.get("classification") != "EXPLORATORY":
        raise NymanSummaryError("Nyman index classification changed")
    if index.get("hypothesis_status") != "UNRESOLVED":
        raise NymanSummaryError("Nyman index changed hypothesis status")
    if index.get("limitation") != search.EXPLORATORY_LIMITATION:
        raise NymanSummaryError("Nyman index limitation changed")
    if index.get("backend") != runtime_backend:
        raise NymanSummaryError("Nyman index backend binding changed")
    if index.get("plan") != {
        "path": "plan.json",
        "plan_sha256": plan["plan_sha256"],
    } or index.get("run") != {"path": "run.json"}:
        raise NymanSummaryError("Nyman index source references changed")

    references = index.get("cells")
    if not isinstance(references, list) or len(references) != len(contracts):
        raise NymanSummaryError("summary requires exactly six Nyman cell references")
    if index.get("missing_cell_ids") != []:
        raise NymanSummaryError("summary requires a complete Nyman index")

    kernel_reference = _require_fields(
        index.get("generation_kernel"),
        {
            "path",
            "payload_sha256",
            "precision_bits",
            "kernel_contract_sha256",
            "core_system_content_sha256",
        },
        "generation kernel reference",
    )
    if kernel_reference.get("path") != "kernels/generation.json":
        raise NymanSummaryError("generation kernel path changed")
    kernel = _load_object(_safe_path(root, kernel_reference["path"]))
    prefixes = _validate_kernel(kernel, plan, contracts)
    expected_kernel_reference = {
        "path": "kernels/generation.json",
        "payload_sha256": kernel["payload_sha256"],
        "precision_bits": kernel["precision_bits"],
        "kernel_contract_sha256": kernel["kernel_contract_sha256"],
        "core_system_content_sha256": kernel["core_system_content_sha256"],
    }
    if dict(kernel_reference) != expected_kernel_reference:
        raise NymanSummaryError("generation kernel reference changed")

    completed: dict[str, dict[str, Any]] = {}
    cell_summaries: list[dict[str, Any]] = []
    status_counts: Counter[str] = Counter()
    total_coefficients = 0
    total_pivots = 0
    seen: set[str] = set()
    by_id = {contract.cell_id: contract for contract in contracts}
    reference_fields = {
        "cell_id",
        "cell_contract_sha256",
        "path",
        "payload_sha256",
        "candidate_sha256",
        "status",
    }
    for reference in references:
        item = _require_fields(reference, reference_fields, "Nyman cell reference")
        cell_id = item.get("cell_id")
        if not isinstance(cell_id, str) or cell_id not in by_id or cell_id in seen:
            raise NymanSummaryError("unknown or duplicate Nyman cell reference")
        seen.add(cell_id)
        contract = by_id[cell_id]
        expected_path = f"cells/{contract.cell_id}.json"
        if item.get("path") != expected_path:
            raise NymanSummaryError("Nyman cell path changed")
        if item.get("cell_contract_sha256") != contract.to_record()[
            "cell_contract_sha256"
        ]:
            raise NymanSummaryError("Nyman cell-contract hash changed")
        artifact, lower, upper, energy_display, pivot_count = _validate_cell(
            _load_object(_safe_path(root, item["path"])),
            contract,
            plan,
            kernel,
            prefixes[str(contract.n)],
        )
        if item.get("payload_sha256") != artifact["payload_sha256"]:
            raise NymanSummaryError("Nyman cell payload reference changed")
        if item.get("candidate_sha256") != artifact["candidate"][
            "candidate_sha256"
        ]:
            raise NymanSummaryError("Nyman candidate reference changed")
        if item.get("status") != artifact["status"]:
            raise NymanSummaryError("Nyman cell status reference changed")

        width = upper - lower
        status_counts[str(artifact["status"])] += 1
        total_coefficients += contract.n
        total_pivots += pivot_count
        completed[cell_id] = artifact
        generation = artifact["generation"]
        cell_summaries.append(
            {
                "cell_id": cell_id,
                "n": str(contract.n),
                "status": artifact["status"],
                "terminal_statement": artifact["terminal_statement"],
                "bounds": {
                    "lower": _bound_record(
                        artifact["candidate"]["lower_bound"], lower
                    ),
                    "upper": _bound_record(
                        artifact["candidate"]["upper_bound"], upper
                    ),
                    "width": {
                        "dyadic": _canonical_dyadic_record(width),
                        "exact_fraction": _fraction_record(width),
                    },
                },
                "stored_energy_display": energy_display,
                "decisions": {
                    "upper": generation["upper_certificate"]["decision"],
                    "lower": generation["lower_certificate"]["decision"],
                },
                "hashes": {
                    "cell_contract_sha256": item["cell_contract_sha256"],
                    "cell_payload_sha256": item["payload_sha256"],
                    "candidate_sha256": item["candidate_sha256"],
                    "upper_certificate_payload_sha256": generation[
                        "upper_certificate"
                    ]["payload_sha256"],
                    "lower_certificate_payload_sha256": generation[
                        "lower_certificate"
                    ]["payload_sha256"],
                    "prefix_core_system_content_sha256": prefixes[str(contract.n)][
                        "core_system_content_sha256"
                    ],
                    "prefix_kernel_sha256": prefixes[str(contract.n)][
                        "prefix_kernel_sha256"
                    ],
                },
            }
        )

    if seen != set(by_id):
        raise NymanSummaryError("complete index does not cover the frozen plan")
    try:
        expected_index = search._index_payload(plan, contracts, completed, kernel)
    except search.NymanSearchError as exc:
        raise NymanSummaryError("Nyman index cannot be structurally derived") from exc
    if index != expected_index:
        raise NymanSummaryError("Nyman index is not canonical for its artifacts")

    diagnostics = index["strict_improvement_diagnostics"]
    decision_counts = Counter(str(item["decision"]) for item in diagnostics)
    for status in ("FINITE_DISTANCE_BRACKET_CERTIFIED", "INCONCLUSIVE"):
        status_counts.setdefault(status, 0)

    body = {
        "schema": SUMMARY_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "conclusion": index["conclusion"],
        "source": {
            "plan": {
                "path": "plan.json",
                "plan_sha256": plan["plan_sha256"],
            },
            "run": {
                "path": "run.json",
                "payload_sha256": run["payload_sha256"],
            },
            "index": {
                "path": "index.json",
                "payload_sha256": index["payload_sha256"],
            },
        },
        "frozen_experiment": search.FROZEN_EXPERIMENT,
        "generation_kernel": {
            "path": "kernels/generation.json",
            "precision_bits": kernel["precision_bits"],
            "max_n": kernel["max_n"],
            "payload_sha256": kernel["payload_sha256"],
            "kernel_contract_sha256": kernel["kernel_contract_sha256"],
            "core_system_content_sha256": kernel[
                "core_system_content_sha256"
            ],
            "prefixes": [
                {
                    "n": prefix["n"],
                    "core_system_content_sha256": prefix[
                        "core_system_content_sha256"
                    ],
                    "prefix_kernel_sha256": prefix["prefix_kernel_sha256"],
                }
                for prefix in kernel["prefixes"]
            ],
        },
        "counts": {
            "cells": str(len(cell_summaries)),
            "by_status": _count_record(status_counts),
            "stored_exact_dyadic_coefficients": str(total_coefficients),
            "stored_ldlt_pivots": str(total_pivots),
            "strict_improvement_diagnostics": str(len(diagnostics)),
            "strict_improvement_by_decision": _count_record(decision_counts),
        },
        "cells": cell_summaries,
        "strict_improvement": {
            "role": index["strict_improvement_role"],
            "diagnostics": diagnostics,
        },
        "summary_method": {
            "numerical_replay_performed": False,
            "generation_kernel_rebuilt": False,
            "stored_numerical_values_treated_as_hash_bound_records": True,
            "content_hashes_checked": [
                "plan.plan_sha256",
                "run.payload_sha256",
                "index.payload_sha256",
                "generation_kernel.payload_sha256",
                "cell.payload_sha256",
                "candidate.candidate_sha256",
                "upper_certificate.payload_sha256",
                "lower_certificate.payload_sha256",
            ],
        },
        "limitation": SUMMARY_LIMITATION,
    }
    return {**body, "payload_sha256": content_sha256(body)}


def verify_nyman_summary(
    summary: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
) -> dict[str, Any]:
    """Regenerate a compact summary and reject every supplied difference."""

    if isinstance(summary, Path):
        supplied = _load_object(summary, NymanSummaryVerificationError)
    elif isinstance(summary, Mapping):
        supplied = dict(summary)
    else:
        raise TypeError("summary must be a mapping or Path")
    body = {key: value for key, value in supplied.items() if key != "payload_sha256"}
    if supplied.get("payload_sha256") != content_sha256(body):
        raise NymanSummaryVerificationError("summary payload hash mismatch")
    expected = generate_nyman_summary(checkpoint_dir)
    if supplied != expected:
        raise NymanSummaryVerificationError(
            "summary does not canonically regenerate from the checkpoint"
        )
    return {
        "classification": "REPRODUCED_EXPLORATORY_NYMAN_SUMMARY",
        "hypothesis_status": "UNRESOLVED",
        "summary_payload_sha256": supplied["payload_sha256"],
        "plan_sha256": supplied["source"]["plan"]["plan_sha256"],
        "index_payload_sha256": supplied["source"]["index"]["payload_sha256"],
        "verified_cells": supplied["counts"]["cells"],
        "numerical_replay_performed": False,
        "generation_kernel_rebuilt": False,
    }


__all__ = [
    "FROZEN_PLAN_SHA256",
    "NymanSummaryError",
    "NymanSummaryVerificationError",
    "SUMMARY_LIMITATION",
    "SUMMARY_SCHEMA",
    "generate_nyman_summary",
    "verify_nyman_summary",
]
