"""Rigorous finite verification of zeta zeros using FLINT/Arb."""

from __future__ import annotations

from typing import Any

import flint
from flint import acb, arb, ctx

from .artifacts import content_sha256
from .balls import arb_from_dyadic, arb_record


FLINT_ZERO_DOCS = "https://flintlib.org/doc/acb_dirichlet.html"
ZERO_CLAIM = (
    "The isolated critical-line Hardy-Z roots in the stated prefix exhaust the "
    "nontrivial zeta zeros through the final separator, as certified by Turing counts."
)
ZERO_LIMITATION = (
    "This is a finite, same-backend calibration far below the published record. "
    "The FLINT critical-line routine is not an off-line-zero search and this result "
    "cannot prove the Riemann Hypothesis."
)
ZERO_ROUTINE = "flint.acb.zeta_zeros and flint.arb.zeta_nzeros"


class ZeroCertificateError(ValueError):
    """Raised when a serialized zero certificate fails replay."""


def certify_critical_line_zeros(
    first_index: int = 1,
    count: int = 100,
    precision_bits: int = 192,
    count_block_size: int = 1000,
) -> dict[str, Any]:
    """Certify consecutive positive-ordinate zeta zeros on the critical line.

    FLINT's numbered zero isolation uses Hardy-Z isolation together with
    Turing's method.  We additionally place a separator enclosure after each
    requested block and ask ``zeta_nzeros`` for the total number of zeros below
    it.  This makes completeness checks explicit in the artifact.
    """

    for name, value in (("first_index", first_index), ("count", count)):
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"{name} must be an integer")
        if value < 1:
            raise ValueError(f"{name} must be positive")
    if precision_bits < 64:
        raise ValueError("precision_bits must be at least 64")
    if count_block_size < 1:
        raise ValueError("count_block_size must be positive")

    previous_precision = ctx.prec
    try:
        ctx.prec = precision_bits
        # One look-ahead zero supplies a separator beyond the requested range.
        isolated = acb.zeta_zeros(first_index, count + 1)
        selected = isolated[:count]
        half = arb(1) / 2

        embedded_on_critical_line = all(
            zero.real.is_exact() and zero.real == half for zero in selected
        )
        positive_imaginary = all(zero.imag > 0 for zero in selected)
        strictly_ordered = all(
            left.imag < right.imag for left, right in zip(selected, selected[1:])
        )

        records = [
            {
                "index": str(first_index + offset),
                "real": {
                    "numerator": "1",
                    "denominator": "2",
                    "exact": zero.real.is_exact() and zero.real == half,
                },
                "imaginary": arb_record(zero.imag),
            }
            for offset, zero in enumerate(selected)
        ]

        block_checks: list[dict[str, Any]] = []
        last_index = first_index + count - 1
        if first_index == 1:
            block_end = min(count_block_size, count)
            while block_end <= count:
                left = isolated[block_end - 1].imag.upper()
                right = isolated[block_end].imag.lower()
                separator = (left + right) / 2
                count_ball = separator.zeta_nzeros()
                unique_count = count_ball.unique_fmpz()
                between_adjacent_roots = (
                    selected[block_end - 1].imag < separator
                    and separator < isolated[block_end].imag
                )
                block_checks.append(
                    {
                        "last_index": str(block_end),
                        "separator": arb_record(separator),
                        "total_zero_count": (
                            str(unique_count) if unique_count is not None else None
                        ),
                        "count_certified": (
                            unique_count is not None
                            and int(unique_count) == block_end
                            and between_adjacent_roots
                        ),
                        "between_adjacent_root_balls": between_adjacent_roots,
                    }
                )
                if block_end == count:
                    break
                block_end = min(block_end + count_block_size, count)

        complete_prefix = (
            first_index == 1
            and bool(block_checks)
            and int(block_checks[-1]["last_index"]) == count
            and all(check["count_certified"] for check in block_checks)
        )
        finite_range_certified = (
            embedded_on_critical_line
            and positive_imaginary
            and strictly_ordered
            and (complete_prefix or first_index > 1)
        )
        transcript_sha256 = content_sha256(records)
        payload = {
            "schema": "rh-lab/zero-certificate/v1",
            "classification": (
                "CERTIFIED_FINITE" if finite_range_certified else "INCONCLUSIVE"
            ),
            "claim": ZERO_CLAIM,
            "scope": {
                "first_index": str(first_index),
                "last_index": str(last_index),
                "positive_imaginary_ordinate": True,
                "multiplicity_counted": True,
            },
            "precision_bits": str(precision_bits),
            "backend": {
                "python_flint": flint.__version__,
                "flint": flint.__FLINT_VERSION__,
                "routine": ZERO_ROUTINE,
            },
            "checks": {
                "hardy_z_roots_embedded_on_critical_line": embedded_on_critical_line,
                "imaginary_parts_positive": positive_imaginary,
                "imaginary_balls_strictly_ordered": strictly_ordered,
                "complete_prefix_by_turing_counts": complete_prefix,
                "blocks": block_checks,
            },
            "zeros": records,
            "zeros_sha256": transcript_sha256,
            "source": FLINT_ZERO_DOCS,
            "limitation": ZERO_LIMITATION,
        }
        return {**payload, "payload_sha256": content_sha256(payload)}
    finally:
        ctx.prec = previous_precision


def verify_zero_certificate(
    artifact: dict[str, Any], replay_precision_bits: int = 384
) -> dict[str, Any]:
    """Reject corrupted records and replay the complete transcript.

    This is a separate invocation of the same FLINT backend, not an independent
    mathematical implementation. It is therefore a consistency-reproduction
    gate, not the clean-room audit required for a global resolution claim.
    """

    supplied_hash = artifact.get("payload_sha256")
    payload = {key: value for key, value in artifact.items() if key != "payload_sha256"}
    if supplied_hash != content_sha256(payload):
        raise ZeroCertificateError("payload hash mismatch")
    if artifact.get("schema") != "rh-lab/zero-certificate/v1":
        raise ZeroCertificateError("unexpected schema")
    if artifact.get("classification") != "CERTIFIED_FINITE":
        raise ZeroCertificateError("artifact is not a finite certificate")
    if artifact.get("claim") != ZERO_CLAIM:
        raise ZeroCertificateError("unexpected claim text")
    if artifact.get("limitation") != ZERO_LIMITATION:
        raise ZeroCertificateError("unexpected limitation text")
    if artifact.get("source") != FLINT_ZERO_DOCS:
        raise ZeroCertificateError("unexpected source")
    expected_backend = {
        "python_flint": flint.__version__,
        "flint": flint.__FLINT_VERSION__,
        "routine": ZERO_ROUTINE,
    }
    if artifact.get("backend") != expected_backend:
        raise ZeroCertificateError("backend manifest does not match this verifier")

    records = artifact.get("zeros")
    if not isinstance(records, list) or not records:
        raise ZeroCertificateError("zero transcript is empty")
    if artifact.get("zeros_sha256") != content_sha256(records):
        raise ZeroCertificateError("zero transcript hash mismatch")

    try:
        first_index = int(artifact["scope"]["first_index"])
        last_index = int(artifact["scope"]["last_index"])
        original_bits = int(artifact["precision_bits"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ZeroCertificateError("invalid scope or precision") from exc
    if len(records) != last_index - first_index + 1:
        raise ZeroCertificateError("transcript length does not match scope")
    if artifact.get("scope", {}).get("positive_imaginary_ordinate") is not True:
        raise ZeroCertificateError("positive-ordinate scope flag is missing")
    if artifact.get("scope", {}).get("multiplicity_counted") is not True:
        raise ZeroCertificateError("multiplicity scope flag is missing")
    if replay_precision_bits <= original_bits:
        raise ZeroCertificateError("replay precision must exceed original precision")

    checks = artifact.get("checks", {})
    required_true = (
        "hardy_z_roots_embedded_on_critical_line",
        "imaginary_parts_positive",
        "imaginary_balls_strictly_ordered",
    )
    if any(checks.get(key) is not True for key in required_true):
        raise ZeroCertificateError("generation checks were not all conclusive")

    previous_precision = ctx.prec
    try:
        ctx.prec = replay_precision_bits
        decoded: list[Any] = []
        for offset, record in enumerate(records):
            expected_index = first_index + offset
            if record.get("index") != str(expected_index):
                raise ZeroCertificateError("zero indices are not consecutive")
            real = record.get("real", {})
            if real != {"numerator": "1", "denominator": "2", "exact": True}:
                raise ZeroCertificateError(
                    f"zero {expected_index} has invalid embedded real part"
                )
            try:
                imaginary = record["imaginary"]
                ball = arb_from_dyadic(imaginary["dyadic"])
                displayed = arb(imaginary["display"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ZeroCertificateError(
                    f"zero {expected_index} has invalid imaginary enclosure"
                ) from exc
            if imaginary.get("is_exact") is not ball.is_exact():
                raise ZeroCertificateError(
                    f"zero {expected_index} has inconsistent exactness metadata"
                )
            if not displayed.contains(ball):
                raise ZeroCertificateError(
                    f"zero {expected_index} display does not contain its dyadic ball"
                )
            if not ball.is_finite() or not ball > 0:
                raise ZeroCertificateError(
                    f"zero {expected_index} is not positive finite"
                )
            if decoded and not decoded[-1] < ball:
                raise ZeroCertificateError("zero balls overlap or are out of order")
            decoded.append(ball)

        replayed_zeros = acb.zeta_zeros(first_index, len(records) + 1)
        if len(replayed_zeros) != len(records) + 1:
            raise ZeroCertificateError("higher-precision batch replay is incomplete")
        half = arb(1) / 2
        for offset, (stored, replayed) in enumerate(zip(decoded, replayed_zeros)):
            if offset == len(records):
                break
            index = first_index + offset
            if not (replayed.real.is_exact() and replayed.real == half):
                raise ZeroCertificateError(
                    f"Hardy-Z replay was not embedded at 1/2 for index {index}"
                )
            if not stored.contains(replayed.imag):
                raise ZeroCertificateError(
                    f"stored ball does not contain replayed root at index {index}"
                )

        reproduced_blocks: list[str] = []
        if first_index == 1:
            if checks.get("complete_prefix_by_turing_counts") is not True:
                raise ZeroCertificateError("prefix completeness was not certified")
            blocks = checks.get("blocks")
            if not isinstance(blocks, list) or not blocks:
                raise ZeroCertificateError("prefix certificate has no count blocks")
            previous_end = 0
            for block in blocks:
                expected_count = int(block["last_index"])
                if not previous_end < expected_count <= len(records):
                    raise ZeroCertificateError("count blocks are not strictly increasing")
                if block.get("count_certified") is not True:
                    raise ZeroCertificateError("stored count block is not certified")
                if block.get("between_adjacent_root_balls") is not True:
                    raise ZeroCertificateError("stored separator placement is not certified")
                if block.get("total_zero_count") != str(expected_count):
                    raise ZeroCertificateError("stored total count contradicts block end")
                separator_record = block["separator"]
                separator = arb_from_dyadic(separator_record["dyadic"])
                displayed_separator = arb(separator_record["display"])
                if separator_record.get("is_exact") is not separator.is_exact():
                    raise ZeroCertificateError("separator exactness metadata is inconsistent")
                if not displayed_separator.contains(separator):
                    raise ZeroCertificateError("separator display misses its dyadic ball")
                lower_root = decoded[expected_count - 1]
                upper_root = (
                    decoded[expected_count]
                    if expected_count < len(records)
                    else replayed_zeros[len(records)].imag
                )
                if not (lower_root < separator and separator < upper_root):
                    raise ZeroCertificateError(
                        f"separator is not between roots at index {expected_count}"
                    )
                unique_count = separator.zeta_nzeros().unique_fmpz()
                if unique_count is None or int(unique_count) != expected_count:
                    raise ZeroCertificateError(
                        f"Turing count replay failed at index {expected_count}"
                    )
                reproduced_blocks.append(str(expected_count))
                previous_end = expected_count
            if previous_end != len(records):
                raise ZeroCertificateError("final count block does not close the prefix")
    except (KeyError, TypeError, ValueError) as exc:
        if isinstance(exc, ZeroCertificateError):
            raise
        raise ZeroCertificateError("invalid block certificate") from exc
    finally:
        ctx.prec = previous_precision

    return {
        "classification": "REPRODUCED",
        "original_precision_bits": str(original_bits),
        "replay_precision_bits": str(replay_precision_bits),
        "reproduced_block_counts": reproduced_blocks,
        "reproduced_zero_count": str(len(records)),
        "zeros_sha256": artifact["zeros_sha256"],
    }
