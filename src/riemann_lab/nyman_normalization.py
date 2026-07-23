"""Tracked normalization evidence for the finite Nyman experiment.

The bundle in this module compares the closed rational autocorrelation formula
used by :mod:`riemann_lab.nyman` with the independently implemented truncated
integral in :mod:`riemann_lab.nyman_oracle`.  It also compares Arb's
``1-EulerGamma`` value with the oracle's exact-harmonic-sum enclosure.

This is a normalization audit, not an independent high-precision backend and
not evidence that resolves the Riemann Hypothesis.
"""

from __future__ import annotations

from contextlib import contextmanager
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import flint
from flint import arb, ctx

from .artifacts import content_sha256
from .balls import arb_record
from . import nyman
from . import nyman_oracle as oracle


BUNDLE_SCHEMA = "rh-lab/nyman-normalization-audit-bundle/v1"
HARMONIC_AUDIT_SCHEMA = "rh-lab/nyman-harmonic-normalization-audit/v1"
FROZEN_BUNDLE_ID = "nyman-natural-normalization-v1"
FORMULA_EVALUATOR_LABEL = (
    "riemann_lab.nyman.autocorrelation_a rational Vasyunin formula"
)

NORMALIZATION_LIMITATION = (
    "The truncated-integral oracle has only about 12 bits of resolving power "
    "because its generic tail is 1/4096. The formula and oracle use the same "
    "Arb/FLINT backend, so this is not a clean-room or alternate-backend "
    "verification. A passing finite normalization audit has no implication "
    "for RH; the global hypothesis status remains UNRESOLVED."
)


class NymanNormalizationError(ValueError):
    """Base class for malformed normalization evidence."""


class NymanNormalizationVerificationError(NymanNormalizationError):
    """Raised when a stored bundle does not exactly regenerate."""


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
    artifact: Mapping[str, Any],
    error_type: type[NymanNormalizationError],
    label: str,
) -> None:
    supplied = artifact.get("payload_sha256")
    body = {key: value for key, value in artifact.items() if key != "payload_sha256"}
    if supplied != content_sha256(body):
        raise error_type(f"{label} payload hash mismatch")


@contextmanager
def _clean_precision(precision_bits: int) -> Iterable[None]:
    """Run one bundle from a clean FLINT cache and restore global precision."""

    previous_precision = ctx.prec
    ctx.cleanup()
    ctx.prec = precision_bits
    try:
        yield
    finally:
        ctx.prec = previous_precision
        ctx.cleanup()


def _frozen_scope() -> dict[str, Any]:
    return {
        "ordered_supplied_ratio_pairs": [
            {"numerator": str(numerator), "denominator": str(denominator)}
            for numerator, denominator in oracle.FROZEN_RATIO_PAIRS
        ],
        "ratio_count": str(len(oracle.FROZEN_RATIO_PAIRS)),
        "truncation": str(oracle.FROZEN_TRUNCATION),
        "harmonic_last_index": str(oracle.FROZEN_HARMONIC_M),
        "precision_bits": str(oracle.FROZEN_ORACLE_PRECISION_BITS),
        "generic_tail_upper_bound": {
            "numerator": "1",
            "denominator": str(oracle.FROZEN_TRUNCATION),
        },
    }


def _harmonic_audit() -> dict[str, Any]:
    precision_bits = oracle.FROZEN_ORACLE_PRECISION_BITS
    harmonic = oracle.one_minus_euler_gamma_oracle(
        oracle.FROZEN_HARMONIC_M,
        precision_bits,
    )
    core_value = arb(1) - arb.const_euler()
    overlaps = bool(core_value.overlaps(harmonic.enclosure))
    return _with_payload_hash(
        {
            "schema": HARMONIC_AUDIT_SCHEMA,
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "audit_outcome": (
                "CORE_ONE_MINUS_EULER_GAMMA_OVERLAPS_HARMONIC_ORACLE"
                if overlaps
                else "CORE_ONE_MINUS_EULER_GAMMA_DISJOINT_FROM_HARMONIC_ORACLE"
            ),
            "core_expression": "1-arb.const_euler()",
            "core_enclosure": arb_record(core_value),
            "independent_harmonic_oracle": harmonic.to_record(),
            "comparison": "rigorous Arb enclosure overlap",
            "precision_bits": str(precision_bits),
            "limitation": (
                "This checks one normalization constant using a harmonic-sum "
                "bound. Both sides are enclosed with the same Arb backend, "
                "and the result has no implication for RH."
            ),
        }
    )


def generate_nyman_normalization_bundle() -> dict[str, Any]:
    """Generate the canonical frozen fourteen-ratio normalization bundle."""

    precision_bits = oracle.FROZEN_ORACLE_PRECISION_BITS
    with _clean_precision(precision_bits):
        formula_audit = oracle.audit_frozen_ratio_normalization(
            nyman.autocorrelation_a,
            evaluator_label=FORMULA_EVALUATOR_LABEL,
            precision_bits=precision_bits,
        )
        harmonic_audit = _harmonic_audit()

        formula_passed = formula_audit.get("audit_outcome") == (
            "ALL_14_CONSISTENT_WITH_TRUNCATION_ORACLE"
        )
        harmonic_passed = harmonic_audit.get("audit_outcome") == (
            "CORE_ONE_MINUS_EULER_GAMMA_OVERLAPS_HARMONIC_ORACLE"
        )
        passed = formula_passed and harmonic_passed
        payload = {
            "schema": BUNDLE_SCHEMA,
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "bundle_id": FROZEN_BUNDLE_ID,
            "audit_outcome": (
                "NORMALIZATION_AUDIT_PASSED_EXPLORATORY"
                if passed
                else "NORMALIZATION_AUDIT_FAILED_EXPLORATORY"
            ),
            "frozen_scope": _frozen_scope(),
            "formula_vs_truncated_integral": formula_audit,
            "one_minus_euler_gamma": harmonic_audit,
            "backend": _backend_record(),
            "cache_policy": (
                "clean FLINT cache before and after the complete deterministic "
                "128-bit bundle generation"
            ),
            "independence": {
                "truncated_integral_implementation": (
                    "riemann_lab.nyman_oracle; exact rational breakpoint "
                    "partition and piece antiderivative; no nyman formula import"
                ),
                "formula_implementation": (
                    "riemann_lab.nyman.autocorrelation_a; rational Vasyunin formula"
                ),
                "harmonic_implementation": (
                    "exact Fraction harmonic sum plus an Arb-enclosed logarithm"
                ),
                "shared_arithmetic_backend": "python-flint/Arb",
                "clean_room_or_alternate_backend": False,
            },
            "limitation": NORMALIZATION_LIMITATION,
        }
        return _with_payload_hash(payload)


def _load_bundle(
    artifact: Mapping[str, Any] | Path,
) -> dict[str, Any]:
    if isinstance(artifact, Path):
        def reject_duplicate_keys(
            pairs: list[tuple[str, Any]],
        ) -> dict[str, Any]:
            value: dict[str, Any] = {}
            for key, item in pairs:
                if key in value:
                    raise NymanNormalizationVerificationError(
                        f"duplicate JSON object key: {key}"
                    )
                value[key] = item
            return value

        def reject_nonstandard_constant(value: str) -> Any:
            raise NymanNormalizationVerificationError(
                f"nonstandard JSON constant: {value}"
            )

        try:
            loaded = json.loads(
                artifact.read_text(encoding="utf-8"),
                object_pairs_hook=reject_duplicate_keys,
                parse_constant=reject_nonstandard_constant,
            )
        except NymanNormalizationVerificationError:
            raise
        except (OSError, json.JSONDecodeError) as exc:
            raise NymanNormalizationVerificationError(
                f"cannot read normalization bundle: {artifact}"
            ) from exc
        if not isinstance(loaded, dict):
            raise NymanNormalizationVerificationError(
                "normalization bundle must contain an object"
            )
        return loaded
    if isinstance(artifact, Mapping):
        return dict(artifact)
    raise TypeError("artifact must be a mapping or Path")


def _validate_bundle_structure(bundle: Mapping[str, Any]) -> None:
    _require_payload_hash(
        bundle,
        NymanNormalizationVerificationError,
        "normalization bundle",
    )
    allowed = {
        "schema",
        "classification",
        "hypothesis_status",
        "bundle_id",
        "audit_outcome",
        "frozen_scope",
        "formula_vs_truncated_integral",
        "one_minus_euler_gamma",
        "backend",
        "cache_policy",
        "independence",
        "limitation",
        "payload_sha256",
    }
    if set(bundle) != allowed:
        raise NymanNormalizationVerificationError(
            "normalization bundle fields changed"
        )
    if bundle.get("schema") != BUNDLE_SCHEMA:
        raise NymanNormalizationVerificationError(
            "unexpected normalization bundle schema"
        )
    if bundle.get("classification") != "EXPLORATORY":
        raise NymanNormalizationVerificationError(
            "normalization bundle was improperly promoted"
        )
    if bundle.get("hypothesis_status") != "UNRESOLVED":
        raise NymanNormalizationVerificationError(
            "normalization bundle changed RH status"
        )
    if bundle.get("bundle_id") != FROZEN_BUNDLE_ID:
        raise NymanNormalizationVerificationError(
            "normalization bundle id changed"
        )
    if bundle.get("frozen_scope") != _frozen_scope():
        raise NymanNormalizationVerificationError(
            "normalization frozen scope changed"
        )
    if bundle.get("backend") != _backend_record():
        raise NymanNormalizationVerificationError(
            "normalization backend changed"
        )
    if bundle.get("limitation") != NORMALIZATION_LIMITATION:
        raise NymanNormalizationVerificationError(
            "normalization limitation changed"
        )

    formula = bundle.get("formula_vs_truncated_integral")
    harmonic = bundle.get("one_minus_euler_gamma")
    if not isinstance(formula, Mapping) or not isinstance(harmonic, Mapping):
        raise NymanNormalizationVerificationError(
            "normalization sub-audits must be objects"
        )
    _require_payload_hash(
        formula,
        NymanNormalizationVerificationError,
        "formula normalization audit",
    )
    _require_payload_hash(
        harmonic,
        NymanNormalizationVerificationError,
        "harmonic normalization audit",
    )
    if formula.get("classification") != "EXPLORATORY" or formula.get(
        "hypothesis_status"
    ) != "UNRESOLVED":
        raise NymanNormalizationVerificationError(
            "formula normalization audit was promoted"
        )
    if formula.get("evaluator_label") != FORMULA_EVALUATOR_LABEL:
        raise NymanNormalizationVerificationError(
            "formula evaluator binding changed"
        )
    comparisons = formula.get("comparisons")
    if not isinstance(comparisons, list) or len(comparisons) != len(
        oracle.FROZEN_RATIO_PAIRS
    ):
        raise NymanNormalizationVerificationError(
            "formula normalization comparison count changed"
        )
    if harmonic.get("schema") != HARMONIC_AUDIT_SCHEMA:
        raise NymanNormalizationVerificationError(
            "unexpected harmonic normalization schema"
        )
    if harmonic.get("classification") != "EXPLORATORY" or harmonic.get(
        "hypothesis_status"
    ) != "UNRESOLVED":
        raise NymanNormalizationVerificationError(
            "harmonic normalization audit was promoted"
        )
    harmonic_oracle = harmonic.get("independent_harmonic_oracle")
    if not isinstance(harmonic_oracle, Mapping):
        raise NymanNormalizationVerificationError(
            "harmonic oracle evidence is missing"
        )
    _require_payload_hash(
        harmonic_oracle,
        NymanNormalizationVerificationError,
        "harmonic oracle",
    )


def verify_nyman_normalization_bundle(
    artifact: Mapping[str, Any] | Path,
) -> dict[str, Any]:
    """Exactly regenerate and compare every field of a stored audit bundle."""

    stored = _load_bundle(artifact)
    _validate_bundle_structure(stored)
    regenerated = generate_nyman_normalization_bundle()
    if stored != regenerated:
        raise NymanNormalizationVerificationError(
            "normalization bundle does not exactly regenerate"
        )
    return {
        "classification": (
            "REPRODUCED_EXPLORATORY_NYMAN_NORMALIZATION_AUDIT"
        ),
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": stored["audit_outcome"],
        "verified_ratio_comparisons": str(len(oracle.FROZEN_RATIO_PAIRS)),
        "harmonic_comparison_reproduced": True,
        "harmonic_comparison_passed": stored["one_minus_euler_gamma"].get(
            "audit_outcome"
        )
        == "CORE_ONE_MINUS_EULER_GAMMA_OVERLAPS_HARMONIC_ORACLE",
        "payload_sha256": stored["payload_sha256"],
        "exact_regeneration": True,
        "same_arb_backend_not_clean_room": True,
    }


__all__ = [
    "BUNDLE_SCHEMA",
    "FORMULA_EVALUATOR_LABEL",
    "FROZEN_BUNDLE_ID",
    "HARMONIC_AUDIT_SCHEMA",
    "NORMALIZATION_LIMITATION",
    "NymanNormalizationError",
    "NymanNormalizationVerificationError",
    "generate_nyman_normalization_bundle",
    "verify_nyman_normalization_bundle",
]
