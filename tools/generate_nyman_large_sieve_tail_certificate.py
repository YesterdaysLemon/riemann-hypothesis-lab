"""Generate or independently replay the Fourier/Farey tail certificate."""

from __future__ import annotations

import argparse
from dataclasses import fields
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
from typing import Any, Mapping

import flint
import numpy as np

from riemann_lab.artifacts import content_sha256

if __package__:
    from tools.certify_nyman_balanced_prefix import (
        PrefixBound,
        _arb_from_fraction,
        certified_prefix_gain,
    )
    from tools.certify_nyman_balanced_tail import (
        ExactTailInputs,
        LargeSieveTailBound,
        build_exact_vectors,
        coefficient_sum,
        convolve_coefficients,
        harmonic_sum,
        large_sieve_tail_center_radius,
        require_balanced,
        strict_json_loads,
        tail_center_radius,
    )
    from tools.generate_nyman_balanced_tail_certificate import (
        TRUSTED_COMPUTING_BASE as V1_TRUSTED_COMPUTING_BASE,
        _ball_from_record,
        _dyadic_vector_record,
        _fraction_record,
        _full_gain,
        _parse_dyadic_vector,
        _parse_fraction,
        _prefix_record,
        _require_environment,
        _source_record,
        _strict_lower_proved,
        _wide_ball_record,
    )
else:
    from certify_nyman_balanced_prefix import (
        PrefixBound,
        _arb_from_fraction,
        certified_prefix_gain,
    )
    from certify_nyman_balanced_tail import (
        ExactTailInputs,
        LargeSieveTailBound,
        build_exact_vectors,
        coefficient_sum,
        convolve_coefficients,
        harmonic_sum,
        large_sieve_tail_center_radius,
        require_balanced,
        strict_json_loads,
        tail_center_radius,
    )
    from generate_nyman_balanced_tail_certificate import (
        TRUSTED_COMPUTING_BASE as V1_TRUSTED_COMPUTING_BASE,
        _ball_from_record,
        _dyadic_vector_record,
        _fraction_record,
        _full_gain,
        _parse_dyadic_vector,
        _parse_fraction,
        _prefix_record,
        _require_environment,
        _source_record,
        _strict_lower_proved,
        _wide_ball_record,
    )


SCHEMA = "rh-lab/nyman-large-sieve-tail-certificate/v1"
DEFAULT_OUTPUT = Path("results/nyman-large-sieve-tail-v1.json")
FROZEN_ARTIFACT_PAYLOAD_SHA256 = (
    "2ee75c9af3625872c8e7b9a59ea76e353859180d636ec7446cac09db99e25918"
)
CLAIMED_GAIN_LOWER_BOUND = Fraction(1, 5_000)
PRIOR_CLAIMED_GAIN_LOWER_BOUND = Fraction(7, 50_000)
CONFIG = {
    "n": 256,
    "multiplier_limit": 16,
    "shell_bits": 9,
    "z_bits": 9,
    "old_bits": 16,
    "cutoff": 1 << 22,
    "generation_block_size": 1 << 18,
    "generation_precision_bits": 256,
    "replay_block_size": 1 << 17,
    "replay_precision_bits": 384,
}
PRIOR_CUTOFF = 1 << 26
AS_OF = "2026-07-23"
STATEMENT = (
    "For one explicit exact-dyadic N=256 old Nyman vector and one explicit "
    "balanced K=16 correction supported through index 8192, the complete "
    "infinite weighted-L2 direct gain is greater than 1/5000. A "
    "Fourier/Farey large-sieve bound certifies the tail after 2^22 intervals."
)
LIMITATION = (
    "This proves one complete finite-dimensional contraction and one general "
    "periodic-tail lemma. It supplies no uniform all-scale contraction theorem "
    "and does not prove or disprove the Riemann Hypothesis."
)
REFERENCE = {
    "authors": "H. L. Montgomery and R. C. Vaughan",
    "title": "The large sieve",
    "journal": "Mathematika 20 (1973), 119-134",
    "doi": "10.1112/S0025579300004708",
    "url": "https://doi.org/10.1112/S0025579300004708",
    "used_result": (
        "Theorem 1 in dual form: a delta-separated frequency set has "
        "consecutive-interval square norm at most (H+delta^-1) times its "
        "Fourier coefficient energy."
    ),
}
METHOD = {
    "interval_identity": (
        "gain=sum_M[(2*g_M*q_M-g_M^2)/(M(M+1))"
        "-2*P1*g_M*log((M+1)/M)]"
    ),
    "fourier_coefficients": (
        "at reduced r/d, g_hat(r/d)=-S_d/(1-exp(-2*pi*i*r/d)), "
        "where S_d=sum_{d|n}a_n/n"
    ),
    "farey_spacing": (
        "all active reduced denominators are at most Q, so distinct "
        "frequencies, including zero, are separated by at least 1/(Q(Q-1))"
    ),
    "two_sided_interval_bound": (
        "apply the Montgomery-Vaughan upper large sieve to an interval and "
        "its cyclic complement in a common rational period"
    ),
    "difference_of_squares": (
        "2*g*v-g^2=v^2-(v-g)^2, giving C_LS=Q(Q-1)*(rho+tau), "
        "tau=rho+mu-2*nu"
    ),
    "tail_center": "(2*nu-mu)/(T+1)",
    "tail_radius": "C_LS/((T+1)(T+2))+abs(P1)*A0/(12*T^2)",
    "exact_means": "Jordan-J2 divisor aggregation and exact rational arithmetic",
    "prefix": (
        "exact int64 divisor recurrences, pinned binary64 reductions with an "
        "exact rational forward-error radius, and Arb log/log-gamma Abel "
        "compression"
    ),
    "selection_boundary": (
        "the previously frozen exploratory scout selected the direction; this "
        "certificate recomputes every vector and proof quantity from exact "
        "dyadic source data"
    ),
}
TRUSTED_COMPUTING_BASE = dict(V1_TRUSTED_COMPUTING_BASE)
EXPECTED_CHECKS = {
    "old_support_at_most_256": True,
    "added_support_at_most_8192": True,
    "added_coefficient_sum_zero": True,
    "added_harmonic_sum_zero": True,
    "new_periodic_mean_square_nonnegative": True,
    "large_sieve_bound_strictly_below_lcm_bound": True,
    "cutoff_is_sixteen_times_shorter_than_prior_certificate": True,
    "full_gain_strictly_above_claim": True,
    "resolves_rh": False,
}
TOP_LEVEL_FIELDS = {
    "schema",
    "classification",
    "hypothesis_status",
    "as_of",
    "statement",
    "configuration",
    "claimed_gain_lower_bound",
    "source",
    "vectors",
    "prefix",
    "tail",
    "comparison",
    "full_gain",
    "method",
    "reference",
    "trusted_computing_base",
    "checks",
    "limitation",
    "payload_sha256",
}


def _construct_vectors(root: Path) -> ExactTailInputs:
    return build_exact_vectors(
        root,
        n=CONFIG["n"],
        multiplier_limit=CONFIG["multiplier_limit"],
        shell_bits=CONFIG["shell_bits"],
        z_bits=CONFIG["z_bits"],
        old_bits=CONFIG["old_bits"],
    )


def _tail_record(bound: LargeSieveTailBound) -> dict[str, object]:
    record: dict[str, object] = {}
    for field in fields(LargeSieveTailBound):
        value = getattr(bound, field.name)
        record[field.name] = str(value) if isinstance(value, int) else _fraction_record(value)
    return record


def _verify_tail_record(record: object, replay: LargeSieveTailBound) -> None:
    expected_fields = {field.name for field in fields(LargeSieveTailBound)}
    if not isinstance(record, Mapping) or set(record) != expected_fields:
        raise ValueError("large-sieve tail record fields changed")
    for field in fields(LargeSieveTailBound):
        expected = getattr(replay, field.name)
        supplied = record[field.name]
        if isinstance(expected, int):
            if not isinstance(supplied, str) or str(int(supplied)) != supplied:
                raise ValueError(f"tail {field.name} is not a canonical integer")
            if int(supplied) != expected:
                raise ValueError(f"tail {field.name} changed")
        elif _parse_fraction(supplied, f"tail {field.name}") != expected:
            raise ValueError(f"tail {field.name} changed")


def _comparison_record(
    large_sieve: LargeSieveTailBound,
    lcm_bound: Any,
) -> dict[str, object]:
    return {
        "prior_certificate_cutoff": str(PRIOR_CUTOFF),
        "cutoff_reduction_factor": str(PRIOR_CUTOFF // CONFIG["cutoff"]),
        "prior_claimed_gain_lower_bound": _fraction_record(
            PRIOR_CLAIMED_GAIN_LOWER_BOUND
        ),
        "prior_lcm_periodic_partial_sum_bound": _fraction_record(
            lcm_bound.periodic_partial_sum_bound
        ),
        "large_sieve_periodic_partial_sum_bound": _fraction_record(
            large_sieve.periodic_partial_sum_bound
        ),
        "large_sieve_to_lcm_bound_ratio": _fraction_record(
            large_sieve.periodic_partial_sum_bound
            / lcm_bound.periodic_partial_sum_bound
        ),
    }


def _verify_prefix_record(record: object, replay: PrefixBound) -> None:
    prefix_fields = {
        "cutoff",
        "block_size",
        "block_count",
        "added_denominator",
        "old_denominator",
        "rational_midpoint",
        "rational_error",
        "rational_lower",
        "rational_upper",
        "computed_absolute_term_sum_upper",
        "logarithmic_gain",
        "prefix_gain",
        "last_added_step_numerator",
        "last_old_intercept_numerator",
        "maximum_added_step_numerator",
        "maximum_old_intercept_numerator",
        "maximum_interval_numerator",
    }
    if not isinstance(record, Mapping) or set(record) != prefix_fields:
        raise ValueError("certificate prefix record is missing")
    expected_blocks = (
        CONFIG["cutoff"] + CONFIG["generation_block_size"] - 1
    ) // CONFIG["generation_block_size"]
    integer_expectations = {
        "cutoff": CONFIG["cutoff"],
        "block_size": CONFIG["generation_block_size"],
        "block_count": expected_blocks,
        "added_denominator": replay.added_denominator,
        "old_denominator": replay.old_denominator,
        "last_added_step_numerator": replay.last_added_step_numerator,
        "last_old_intercept_numerator": replay.last_old_intercept_numerator,
        "maximum_added_step_numerator": replay.maximum_added_step_numerator,
        "maximum_old_intercept_numerator": replay.maximum_old_intercept_numerator,
        "maximum_interval_numerator": replay.maximum_interval_numerator,
    }
    for name, expected in integer_expectations.items():
        supplied = record[name]
        if not isinstance(supplied, str) or str(int(supplied)) != supplied:
            raise ValueError(f"stored prefix {name} is not canonical")
        if int(supplied) != expected:
            raise ValueError(f"stored prefix {name} changed")
    midpoint = _parse_fraction(record["rational_midpoint"], "stored midpoint")
    error = _parse_fraction(record["rational_error"], "stored error")
    lower = _parse_fraction(record["rational_lower"], "stored lower")
    upper = _parse_fraction(record["rational_upper"], "stored upper")
    if error < 0 or lower != midpoint - error or upper != midpoint + error:
        raise ValueError("stored rational prefix identities failed")
    if lower > replay.rational_lower or upper < replay.rational_upper:
        raise ValueError("stored rational prefix does not contain replay")
    stored_absolute_term_sum = _parse_fraction(
        record["computed_absolute_term_sum_upper"],
        "stored absolute term sum",
    )
    if stored_absolute_term_sum < 0:
        raise ValueError("stored absolute term sum is negative")
    if (
        stored_absolute_term_sum
        < replay.computed_absolute_term_sum_upper
    ):
        raise ValueError("stored absolute term sum does not contain replay")
    stored_prefix = _ball_from_record(record["prefix_gain"], "stored prefix gain")
    if not stored_prefix.contains(replay.prefix_gain):
        raise ValueError("stored prefix ball does not contain replay")
    stored_log = _ball_from_record(record["logarithmic_gain"], "stored log gain")
    if not stored_log.contains(replay.logarithmic_gain):
        raise ValueError("stored logarithmic ball does not contain replay")


def _verify_full_gain_record(
    record: object,
    replay_full: Any,
) -> Any:
    expected_fields = {
        "enclosure",
        "lower_endpoint",
        "upper_endpoint",
        "interval_display",
    }
    if not isinstance(record, Mapping) or set(record) != expected_fields:
        raise ValueError("stored full gain record changed")
    stored = _ball_from_record(record["enclosure"], "stored full gain")
    if record != _wide_ball_record(stored):
        raise ValueError("stored full gain is not canonical")
    stored_lower = _ball_from_record(record["lower_endpoint"], "stored lower")
    stored_upper = _ball_from_record(record["upper_endpoint"], "stored upper")
    if not stored_lower.overlaps(stored.lower()):
        raise ValueError("stored full-gain lower endpoint changed")
    if not stored_upper.overlaps(stored.upper()):
        raise ValueError("stored full-gain upper endpoint changed")
    threshold = _arb_from_fraction(CLAIMED_GAIN_LOWER_BOUND)
    if not _strict_lower_proved(stored, threshold):
        raise ValueError("stored full gain does not prove the claim")
    if not _strict_lower_proved(replay_full, threshold):
        raise ValueError("replayed full gain does not prove the claim")
    if not stored.contains(replay_full):
        raise ValueError("stored full-gain ball does not contain replay")
    return stored


def generate(root: Path) -> dict[str, object]:
    _require_environment()
    vectors = _construct_vectors(root)
    prefix = certified_prefix_gain(
        vectors.old_coefficients,
        vectors.added_coefficients,
        CONFIG["cutoff"],
        block_size=CONFIG["generation_block_size"],
        precision_bits=CONFIG["generation_precision_bits"],
    )
    tail = large_sieve_tail_center_radius(
        vectors.old_coefficients,
        vectors.added_coefficients,
        CONFIG["cutoff"],
    )
    old_tail = tail_center_radius(
        vectors.old_coefficients,
        vectors.added_coefficients,
        CONFIG["cutoff"],
    )
    full_gain = _full_gain(prefix, tail, CONFIG["generation_precision_bits"])
    if not _strict_lower_proved(
        full_gain,
        _arb_from_fraction(CLAIMED_GAIN_LOWER_BOUND),
    ):
        raise ArithmeticError("large-sieve full gain does not clear the claim")
    checks = {
        "old_support_at_most_256": max(vectors.old_coefficients) <= 256,
        "added_support_at_most_8192": max(vectors.added_coefficients) <= 8192,
        "added_coefficient_sum_zero": coefficient_sum(vectors.added_coefficients)
        == 0,
        "added_harmonic_sum_zero": harmonic_sum(vectors.added_coefficients) == 0,
        "new_periodic_mean_square_nonnegative": (
            tail.new_periodic_mean_square >= 0
        ),
        "large_sieve_bound_strictly_below_lcm_bound": (
            tail.periodic_partial_sum_bound
            < old_tail.periodic_partial_sum_bound
        ),
        "cutoff_is_sixteen_times_shorter_than_prior_certificate": (
            CONFIG["cutoff"] * 16 == PRIOR_CUTOFF
        ),
        "full_gain_strictly_above_claim": True,
        "resolves_rh": False,
    }
    if checks != EXPECTED_CHECKS:
        raise ArithmeticError("generated large-sieve checks did not close")
    body: dict[str, object] = {
        "schema": SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "as_of": AS_OF,
        "statement": STATEMENT,
        "configuration": {key: str(value) for key, value in CONFIG.items()},
        "claimed_gain_lower_bound": _fraction_record(CLAIMED_GAIN_LOWER_BOUND),
        "source": _source_record(vectors),
        "vectors": {
            "old": _dyadic_vector_record(vectors.old_coefficients),
            "shell": _dyadic_vector_record(vectors.shell),
            "multiplier": _dyadic_vector_record(vectors.multiplier),
            "added": _dyadic_vector_record(vectors.added_coefficients),
        },
        "prefix": _prefix_record(prefix),
        "tail": _tail_record(tail),
        "comparison": _comparison_record(tail, old_tail),
        "full_gain": _wide_ball_record(full_gain),
        "method": METHOD,
        "reference": REFERENCE,
        "trusted_computing_base": TRUSTED_COMPUTING_BASE,
        "checks": checks,
        "limitation": LIMITATION,
    }
    return {**body, "payload_sha256": content_sha256(body)}


def verify(root: Path, artifact_path: Path) -> dict[str, object]:
    _require_environment()
    supplied = strict_json_loads(artifact_path.read_text(encoding="utf-8"))
    if not isinstance(supplied, Mapping):
        raise ValueError("certificate must be a JSON object")
    if set(supplied) != TOP_LEVEL_FIELDS:
        raise ValueError("certificate top-level fields changed")
    body = {key: value for key, value in supplied.items() if key != "payload_sha256"}
    if supplied.get("payload_sha256") != content_sha256(body):
        raise ValueError("certificate canonical payload SHA-256 mismatch")
    if supplied.get("payload_sha256") != FROZEN_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("certificate is not the frozen large-sieve v1 payload")
    expected_scalars = {
        "schema": SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "as_of": AS_OF,
        "statement": STATEMENT,
        "limitation": LIMITATION,
    }
    for name, expected in expected_scalars.items():
        if supplied.get(name) != expected:
            raise ValueError(f"certificate {name} changed")
    if supplied.get("method") != METHOD:
        raise ValueError("certificate method changed")
    if supplied.get("reference") != REFERENCE:
        raise ValueError("certificate reference changed")
    if supplied.get("trusted_computing_base") != TRUSTED_COMPUTING_BASE:
        raise ValueError("certificate trusted-computing-base record changed")
    if supplied.get("configuration") != {
        key: str(value) for key, value in CONFIG.items()
    }:
        raise ValueError("certificate configuration changed")
    if _parse_fraction(
        supplied.get("claimed_gain_lower_bound"),
        "claimed gain lower bound",
    ) != CLAIMED_GAIN_LOWER_BOUND:
        raise ValueError("claimed gain lower bound changed")

    vectors = _construct_vectors(root)
    if supplied.get("source") != _source_record(vectors):
        raise ValueError("certificate source binding changed")
    stored_vectors = supplied.get("vectors")
    if not isinstance(stored_vectors, Mapping) or set(stored_vectors) != {
        "old",
        "shell",
        "multiplier",
        "added",
    }:
        raise ValueError("certificate vector bundle changed")
    parsed_old = _parse_dyadic_vector(stored_vectors["old"], "old vector")
    parsed_shell = _parse_dyadic_vector(stored_vectors["shell"], "shell")
    parsed_multiplier = _parse_dyadic_vector(
        stored_vectors["multiplier"],
        "multiplier",
    )
    parsed_added = _parse_dyadic_vector(stored_vectors["added"], "added vector")
    if parsed_old != vectors.old_coefficients:
        raise ValueError("stored old vector changed")
    if parsed_shell != vectors.shell:
        raise ValueError("stored shell changed")
    if parsed_multiplier != vectors.multiplier:
        raise ValueError("stored multiplier changed")
    if parsed_added != vectors.added_coefficients:
        raise ValueError("stored added vector changed")
    require_balanced(parsed_added, "stored added vector")
    if convolve_coefficients(parsed_shell, parsed_multiplier) != parsed_added:
        raise ValueError("stored convolution identity failed")

    replay_prefix = certified_prefix_gain(
        parsed_old,
        parsed_added,
        CONFIG["cutoff"],
        block_size=CONFIG["replay_block_size"],
        precision_bits=CONFIG["replay_precision_bits"],
    )
    replay_tail = large_sieve_tail_center_radius(
        parsed_old,
        parsed_added,
        CONFIG["cutoff"],
    )
    replay_old_tail = tail_center_radius(parsed_old, parsed_added, CONFIG["cutoff"])
    _verify_tail_record(supplied.get("tail"), replay_tail)
    if supplied.get("comparison") != _comparison_record(
        replay_tail,
        replay_old_tail,
    ):
        raise ValueError("certificate comparison record changed")
    _verify_prefix_record(supplied.get("prefix"), replay_prefix)
    replay_full = _full_gain(
        replay_prefix,
        replay_tail,
        CONFIG["replay_precision_bits"],
    )
    _verify_full_gain_record(supplied.get("full_gain"), replay_full)
    if supplied.get("checks") != EXPECTED_CHECKS:
        raise ValueError("certificate check ledger changed")
    return {
        "artifact_payload_sha256": supplied["payload_sha256"],
        "classification": supplied["classification"],
        "hypothesis_status": supplied["hypothesis_status"],
        "claimed_gain_lower_bound": "1/5000",
        "replay_full_gain_lower": replay_full.lower().str(20),
        "replay_full_gain_upper": replay_full.upper().str(20),
        "periodic_bound_reduction_factor_lower": str(
            replay_old_tail.periodic_partial_sum_bound
            // replay_tail.periodic_partial_sum_bound
        ),
        "replay_environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "python_flint": flint.__version__,
            "numpy": np.__version__,
        },
        "verified": True,
    }


def _raw_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    artifact_path = args.artifact
    if not artifact_path.is_absolute():
        artifact_path = root / artifact_path
    if args.verify:
        print(json.dumps(verify(root, artifact_path), sort_keys=True))
        return
    artifact = generate(root)
    artifact_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = artifact_path.with_name(f".{artifact_path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    temporary.replace(artifact_path)
    print(
        json.dumps(
            {
                "artifact": str(artifact_path),
                "payload_sha256": artifact["payload_sha256"],
                "raw_sha256": _raw_sha256(artifact_path),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
