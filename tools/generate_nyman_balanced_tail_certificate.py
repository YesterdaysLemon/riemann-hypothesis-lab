"""Generate or independently replay the balanced-multiplier tail certificate."""

from __future__ import annotations

import argparse
from dataclasses import fields
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
from typing import Any, Mapping

import flint
from flint import ctx
import numpy as np

from riemann_lab.artifacts import content_sha256
from riemann_lab.balls import arb_from_dyadic, arb_record
if __package__:
    from tools.certify_nyman_balanced_prefix import (
        CERTIFICATE_NUMPY_VERSION,
        PrefixBound,
        _arb_from_fraction,
        _radius_ball,
        certified_prefix_gain,
        common_dyadic_numerators,
        require_binary64_tcb,
    )
    from tools.certify_nyman_balanced_tail import (
        ExactTailInputs,
        TailBound,
        build_exact_vectors,
        coefficient_sum,
        convolve_coefficients,
        harmonic_sum,
        require_balanced,
        strict_json_loads,
        tail_center_radius,
    )
else:
    from certify_nyman_balanced_prefix import (
        CERTIFICATE_NUMPY_VERSION,
        PrefixBound,
        _arb_from_fraction,
        _radius_ball,
        certified_prefix_gain,
        common_dyadic_numerators,
        require_binary64_tcb,
    )
    from certify_nyman_balanced_tail import (
        ExactTailInputs,
        TailBound,
        build_exact_vectors,
        coefficient_sum,
        convolve_coefficients,
        harmonic_sum,
        require_balanced,
        strict_json_loads,
        tail_center_radius,
    )


SCHEMA = "rh-lab/nyman-balanced-full-tail-certificate/v1"
DEFAULT_OUTPUT = Path("results/nyman-balanced-full-tail-v1.json")
FROZEN_ARTIFACT_PAYLOAD_SHA256 = (
    "47adcf95031f0e51203c56332beafcfd049374b3446fbe212d4ae7a474958377"
)
CLAIMED_GAIN_LOWER_BOUND = Fraction(7, 50_000)
CONFIG = {
    "n": 256,
    "multiplier_limit": 16,
    "shell_bits": 9,
    "z_bits": 9,
    "old_bits": 16,
    "cutoff": 1 << 26,
    "generation_block_size": 1 << 20,
    "generation_precision_bits": 256,
    "replay_block_size": 1 << 19,
    "replay_precision_bits": 384,
}
AS_OF = "2026-07-23"
STATEMENT = (
    "For one explicit exact-dyadic N=256 old Nyman vector and one "
    "explicit balanced K=16 correction supported through index 8192, "
    "the complete infinite weighted-L2 direct gain is greater than 7/50000."
)
LIMITATION = (
    "This proves one complete finite-dimensional contraction for one explicit "
    "vector. It supplies no uniform all-scale estimate and does not prove or "
    "disprove the Riemann Hypothesis."
)
METHOD = {
    "interval_identity": (
        "gain=sum_M[(2*g_M*q_M-g_M^2)/(M(M+1))"
        "-2*P1*g_M*log((M+1)/M)]"
    ),
    "prefix": (
        "exact int64 divisor recurrences, pinned binary64 reductions with an "
        "exact rational forward-error radius, and Arb log/log-gamma Abel "
        "compression"
    ),
    "tail_center": "(2*nu-mu)/(T+1)",
    "tail_radius": "C_per/((T+1)(T+2))+abs(P1)*A0/(12*T^2)",
    "exact_means": "Jordan-J2 divisor aggregation of periodic covariance",
    "periodic_bound": "ordered-pair LCM sums and Abel summation",
    "selection_boundary": (
        "the exploratory scout selected the direction; Arb proves the first "
        "255 shell rounding bins, the final shell coefficient is imposed by "
        "exact zero-sum balance, and all certificate arithmetic is recomputed "
        "from exact dyadic coefficients"
    ),
}
TRUSTED_COMPUTING_BASE = {
    "python": "CPython >=3.11",
    "python_flint": "0.9.0",
    "numpy": CERTIFICATE_NUMPY_VERSION,
    "binary64": "IEEE-754 round-to-nearest-ties-to-even",
    "reduction": "NumPy add.reduce with an order-independent gamma bound",
}
EXPECTED_CHECKS = {
    "old_support_at_most_256": True,
    "added_support_at_most_8192": True,
    "added_coefficient_sum_zero": True,
    "added_harmonic_sum_zero": True,
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
    "full_gain",
    "method",
    "trusted_computing_base",
    "checks",
    "limitation",
    "payload_sha256",
}


def _require_environment() -> None:
    require_binary64_tcb()
    if platform.python_implementation() != "CPython" or sys.version_info < (3, 11):
        raise RuntimeError("the certificate requires CPython 3.11 or newer")
    if flint.__version__ != TRUSTED_COMPUTING_BASE["python_flint"]:
        raise RuntimeError("the certificate requires python-flint 0.9.0")


def _fraction_record(value: Fraction) -> dict[str, str]:
    exact = Fraction(value)
    return {
        "numerator": str(exact.numerator),
        "denominator": str(exact.denominator),
    }


def _parse_fraction(record: object, name: str) -> Fraction:
    if not isinstance(record, Mapping) or set(record) != {
        "numerator",
        "denominator",
    }:
        raise ValueError(f"{name} is not a canonical rational record")
    value = Fraction(int(record["numerator"]), int(record["denominator"]))
    if str(value.numerator) != record["numerator"] or str(value.denominator) != record[
        "denominator"
    ]:
        raise ValueError(f"{name} is not reduced")
    return value


def _dyadic_vector_record(values: Mapping[int, Fraction]) -> dict[str, Any]:
    denominator, numerators = common_dyadic_numerators(values)
    exponent = denominator.bit_length() - 1
    return {
        "denominator_exponent": str(exponent),
        "entries": [
            {"index": str(index), "numerator": str(numerator)}
            for index, numerator in sorted(numerators.items())
        ],
        "support_count": str(len(numerators)),
        "maximum_index": str(max(numerators, default=0)),
    }


def _parse_dyadic_vector(record: object, name: str) -> dict[int, Fraction]:
    if not isinstance(record, Mapping) or set(record) != {
        "denominator_exponent",
        "entries",
        "support_count",
        "maximum_index",
    }:
        raise ValueError(f"{name} has invalid vector fields")
    exponent = int(record["denominator_exponent"])
    if exponent < 0 or str(exponent) != record["denominator_exponent"]:
        raise ValueError(f"{name} has a negative denominator exponent")
    entries = record["entries"]
    if not isinstance(entries, list):
        raise ValueError(f"{name} entries must be a list")
    result: dict[int, Fraction] = {}
    previous = 0
    for entry in entries:
        if not isinstance(entry, Mapping) or set(entry) != {"index", "numerator"}:
            raise ValueError(f"{name} has an invalid sparse entry")
        index = int(entry["index"])
        numerator = int(entry["numerator"])
        if index <= previous or index < 1 or numerator == 0:
            raise ValueError(f"{name} entries are not canonical")
        if str(index) != entry["index"] or str(numerator) != entry["numerator"]:
            raise ValueError(f"{name} entry is not canonical")
        result[index] = Fraction(numerator, 1 << exponent)
        previous = index
    if int(record["support_count"]) != len(result):
        raise ValueError(f"{name} support count mismatch")
    if int(record["maximum_index"]) != max(result, default=0):
        raise ValueError(f"{name} maximum index mismatch")
    if record != _dyadic_vector_record(result):
        raise ValueError(f"{name} is not a canonical dyadic vector")
    return result


def _tail_record(bound: TailBound) -> dict[str, object]:
    return {
        field.name: (
            str(getattr(bound, field.name))
            if field.name == "cutoff"
            else _fraction_record(getattr(bound, field.name))
        )
        for field in fields(TailBound)
    }


def _verify_tail_record(record: object, replay: TailBound) -> None:
    if not isinstance(record, Mapping) or set(record) != {
        field.name for field in fields(TailBound)
    }:
        raise ValueError("tail record fields changed")
    for field in fields(TailBound):
        if field.name == "cutoff":
            if int(record[field.name]) != replay.cutoff:
                raise ValueError("tail cutoff changed")
        elif _parse_fraction(record[field.name], f"tail {field.name}") != getattr(
            replay,
            field.name,
        ):
            raise ValueError(f"tail {field.name} changed")


def _prefix_record(bound: PrefixBound) -> dict[str, object]:
    return {
        "cutoff": str(bound.cutoff),
        "block_size": str(bound.block_size),
        "block_count": str(bound.block_count),
        "added_denominator": str(bound.added_denominator),
        "old_denominator": str(bound.old_denominator),
        "rational_midpoint": _fraction_record(bound.rational_midpoint),
        "rational_error": _fraction_record(bound.rational_error),
        "rational_lower": _fraction_record(bound.rational_lower),
        "rational_upper": _fraction_record(bound.rational_upper),
        "computed_absolute_term_sum_upper": _fraction_record(
            bound.computed_absolute_term_sum_upper
        ),
        "logarithmic_gain": arb_record(bound.logarithmic_gain),
        "prefix_gain": arb_record(bound.prefix_gain),
        "last_added_step_numerator": str(bound.last_added_step_numerator),
        "last_old_intercept_numerator": str(
            bound.last_old_intercept_numerator
        ),
        "maximum_added_step_numerator": str(
            bound.maximum_added_step_numerator
        ),
        "maximum_old_intercept_numerator": str(
            bound.maximum_old_intercept_numerator
        ),
        "maximum_interval_numerator": str(bound.maximum_interval_numerator),
    }


def _ball_from_record(record: object, name: str) -> Any:
    if not isinstance(record, Mapping) or set(record) != {
        "dyadic",
        "display",
        "is_exact",
    }:
        raise ValueError(f"{name} has invalid Arb fields")
    if not isinstance(record["dyadic"], dict):
        raise ValueError(f"{name} has no dyadic Arb enclosure")
    return arb_from_dyadic(record["dyadic"])


def _wide_ball_record(ball: Any) -> dict[str, object]:
    lower = ball.lower()
    upper = ball.upper()
    return {
        "enclosure": arb_record(ball),
        "lower_endpoint": arb_record(lower),
        "upper_endpoint": arb_record(upper),
        "interval_display": (
            lower.str(30) + " <= gain <= " + upper.str(30)
        ),
    }


def _full_gain(prefix: PrefixBound, tail: TailBound, bits: int) -> Any:
    previous_precision = ctx.prec
    ctx.prec = bits
    try:
        return (
            prefix.prefix_gain
            + _arb_from_fraction(tail.center)
            + _radius_ball(tail.radius, bits)
        )
    finally:
        ctx.prec = previous_precision


def _strict_lower_proved(ball: Any, threshold: Any) -> bool:
    """Return true only when Arb proves the ball's lower endpoint is above."""

    return bool(ball.lower() > threshold)


def _source_record(vectors: ExactTailInputs) -> dict[str, object]:
    return {
        "scout": {
            "path": "results/nyman-balanced-multiplier-scout-v1.json",
            "raw_sha256": vectors.scout_raw_sha256,
            "payload_sha256": vectors.scout_payload_sha256,
            "cell_payload_sha256": vectors.scout_cell_payload_sha256,
        },
        "candidate": {
            "path": vectors.source_candidate_path,
            "raw_sha256": vectors.source_candidate_raw_sha256,
            "payload_sha256": vectors.source_candidate_payload_sha256,
            "candidate_sha256": vectors.source_candidate_candidate_sha256,
        },
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
    tail = tail_center_radius(
        vectors.old_coefficients,
        vectors.added_coefficients,
        CONFIG["cutoff"],
    )
    full_gain = _full_gain(prefix, tail, CONFIG["generation_precision_bits"])
    if not _strict_lower_proved(
        full_gain,
        _arb_from_fraction(CLAIMED_GAIN_LOWER_BOUND),
    ):
        raise ArithmeticError("full-tail gain does not clear the frozen claim")

    checks = {
        "old_support_at_most_256": max(vectors.old_coefficients) <= 256,
        "added_support_at_most_8192": max(vectors.added_coefficients) <= 8192,
        "added_coefficient_sum_zero": coefficient_sum(
            vectors.added_coefficients
        )
        == 0,
        "added_harmonic_sum_zero": harmonic_sum(
            vectors.added_coefficients
        )
        == 0,
        "full_gain_strictly_above_claim": True,
        "resolves_rh": False,
    }
    if checks != EXPECTED_CHECKS:
        raise ArithmeticError("generated certificate checks did not close")

    body: dict[str, object] = {
        "schema": SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "as_of": AS_OF,
        "statement": STATEMENT,
        "configuration": {key: str(value) for key, value in CONFIG.items()},
        "claimed_gain_lower_bound": _fraction_record(
            CLAIMED_GAIN_LOWER_BOUND
        ),
        "source": _source_record(vectors),
        "vectors": {
            "old": _dyadic_vector_record(vectors.old_coefficients),
            "shell": _dyadic_vector_record(vectors.shell),
            "multiplier": _dyadic_vector_record(vectors.multiplier),
            "added": _dyadic_vector_record(vectors.added_coefficients),
        },
        "prefix": _prefix_record(prefix),
        "tail": _tail_record(tail),
        "full_gain": _wide_ball_record(full_gain),
        "method": METHOD,
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
        raise ValueError("certificate is not the frozen v1 payload")
    if supplied.get("schema") != SCHEMA:
        raise ValueError("certificate schema changed")
    if supplied.get("classification") != "CERTIFIED_FINITE":
        raise ValueError("certificate classification changed")
    if supplied.get("hypothesis_status") != "UNRESOLVED":
        raise ValueError("certificate RH status changed")
    if supplied.get("as_of") != AS_OF:
        raise ValueError("certificate date changed")
    if supplied.get("statement") != STATEMENT:
        raise ValueError("certificate statement changed")
    if supplied.get("limitation") != LIMITATION:
        raise ValueError("certificate limitation changed")
    if supplied.get("method") != METHOD:
        raise ValueError("certificate method changed")
    if supplied.get("trusted_computing_base") != TRUSTED_COMPUTING_BASE:
        raise ValueError("certificate trusted-computing-base record changed")
    expected_configuration = {key: str(value) for key, value in CONFIG.items()}
    if supplied.get("configuration") != expected_configuration:
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
    replay_tail = tail_center_radius(parsed_old, parsed_added, CONFIG["cutoff"])
    _verify_tail_record(supplied.get("tail"), replay_tail)

    stored_prefix = supplied.get("prefix")
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
    if not isinstance(stored_prefix, Mapping) or set(stored_prefix) != prefix_fields:
        raise ValueError("certificate prefix record is missing")
    if int(stored_prefix.get("cutoff", -1)) != CONFIG["cutoff"]:
        raise ValueError("stored prefix cutoff changed")
    if int(stored_prefix["block_size"]) != CONFIG["generation_block_size"]:
        raise ValueError("stored generation block size changed")
    expected_block_count = (
        CONFIG["cutoff"] + CONFIG["generation_block_size"] - 1
    ) // CONFIG["generation_block_size"]
    if int(stored_prefix["block_count"]) != expected_block_count:
        raise ValueError("stored generation block count changed")
    if int(stored_prefix["added_denominator"]) != replay_prefix.added_denominator:
        raise ValueError("stored added denominator changed")
    if int(stored_prefix["old_denominator"]) != replay_prefix.old_denominator:
        raise ValueError("stored old denominator changed")
    for name in (
        "last_added_step_numerator",
        "last_old_intercept_numerator",
        "maximum_added_step_numerator",
        "maximum_old_intercept_numerator",
        "maximum_interval_numerator",
    ):
        if int(stored_prefix[name]) != getattr(replay_prefix, name):
            raise ValueError(f"stored prefix {name} changed")
    stored_midpoint = _parse_fraction(
        stored_prefix["rational_midpoint"],
        "stored rational midpoint",
    )
    stored_error = _parse_fraction(
        stored_prefix["rational_error"],
        "stored rational error",
    )
    stored_lower = _parse_fraction(
        stored_prefix.get("rational_lower"),
        "stored rational lower",
    )
    stored_upper = _parse_fraction(
        stored_prefix.get("rational_upper"),
        "stored rational upper",
    )
    if stored_error < 0 or stored_lower != stored_midpoint - stored_error:
        raise ValueError("stored rational lower identity failed")
    if stored_upper != stored_midpoint + stored_error:
        raise ValueError("stored rational upper identity failed")
    if _parse_fraction(
        stored_prefix["computed_absolute_term_sum_upper"],
        "stored absolute term sum",
    ) < 0:
        raise ValueError("stored absolute term sum is negative")
    if stored_lower > replay_prefix.rational_lower:
        raise ValueError("stored rational lower does not contain replay")
    if stored_upper < replay_prefix.rational_upper:
        raise ValueError("stored rational upper does not contain replay")
    stored_prefix_ball = _ball_from_record(
        stored_prefix.get("prefix_gain"),
        "stored prefix gain",
    )
    if not stored_prefix_ball.contains(replay_prefix.prefix_gain):
        raise ValueError("stored prefix ball does not contain replay")
    stored_logarithmic_ball = _ball_from_record(
        stored_prefix["logarithmic_gain"],
        "stored logarithmic gain",
    )
    if not stored_logarithmic_ball.contains(replay_prefix.logarithmic_gain):
        raise ValueError("stored logarithmic ball does not contain replay")

    stored_full_record = supplied.get("full_gain")
    if not isinstance(stored_full_record, Mapping) or set(stored_full_record) != {
        "enclosure",
        "lower_endpoint",
        "upper_endpoint",
        "interval_display",
    }:
        raise ValueError("stored full gain record changed")
    stored_full_ball = _ball_from_record(
        stored_full_record["enclosure"],
        "stored full gain",
    )
    stored_lower_endpoint = _ball_from_record(
        stored_full_record["lower_endpoint"],
        "stored full gain lower endpoint",
    )
    stored_upper_endpoint = _ball_from_record(
        stored_full_record["upper_endpoint"],
        "stored full gain upper endpoint",
    )
    if not stored_lower_endpoint.overlaps(stored_full_ball.lower()):
        raise ValueError("stored full-gain lower endpoint changed")
    if not stored_upper_endpoint.overlaps(stored_full_ball.upper()):
        raise ValueError("stored full-gain upper endpoint changed")
    replay_full_ball = _full_gain(
        replay_prefix,
        replay_tail,
        CONFIG["replay_precision_bits"],
    )
    threshold = _arb_from_fraction(CLAIMED_GAIN_LOWER_BOUND)
    if not _strict_lower_proved(stored_full_ball, threshold):
        raise ValueError("stored full gain does not prove the claim")
    if not _strict_lower_proved(replay_full_ball, threshold):
        raise ValueError("replayed full gain does not prove the claim")
    if not stored_full_ball.contains(replay_full_ball):
        raise ValueError("stored full-gain ball does not contain replay")

    checks = supplied.get("checks")
    if checks != EXPECTED_CHECKS:
        raise ValueError("certificate check ledger changed")
    return {
        "artifact_payload_sha256": supplied["payload_sha256"],
        "classification": supplied["classification"],
        "hypothesis_status": supplied["hypothesis_status"],
        "claimed_gain_lower_bound": "7/50000",
        "replay_full_gain_lower": replay_full_ball.lower().str(20),
        "replay_full_gain_upper": replay_full_ball.upper().str(20),
        "replay_environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "python_flint": flint.__version__,
            "numpy": np.__version__,
            "byteorder": sys.byteorder,
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
