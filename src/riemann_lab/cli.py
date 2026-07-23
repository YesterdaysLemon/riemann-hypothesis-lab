"""Command-line entrypoint for certificate generation and verification."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from .artifacts import write_json
from .claims import ClaimValidationError, load_and_validate_registry
from .lagarias import (
    LagariasCertificateError,
    certify_lagarias_range,
    verify_lagarias_certificate,
)
from .nyman_search import (
    NymanSearchError,
    load_nyman_plan,
    run_nyman_search,
    verify_nyman_search,
)
from .nyman_normalization import (
    NymanNormalizationError,
    generate_nyman_normalization_bundle,
    verify_nyman_normalization_bundle,
)
from .nyman_summary import (
    NymanSummaryError,
    generate_nyman_summary,
    verify_nyman_summary,
)
from .nyman_rebound import (
    NymanReboundError,
    generate_nyman_rebound_audit,
    verify_nyman_rebound_audit,
)
from .nyman_mobius import (
    NymanMobiusError,
    generate_nyman_mobius_audit,
    verify_nyman_mobius_audit,
)
from .nyman_beta2 import (
    FROZEN_REPLAY_BITS as NYMAN_BETA2_REPLAY_BITS,
    NymanBeta2Error,
    generate_nyman_beta2_audit,
    propose_nyman_beta2_candidate,
    verify_nyman_beta2_audit,
)
from .zeros import (
    ZeroCertificateError,
    certify_critical_line_zeros,
    verify_zero_certificate,
)
from .weil import (
    WeilCertificateError,
    certify_weil_matrix,
    verify_weil_certificate,
)
from .weil_audits import (
    WeilAuditVerificationError,
    certify_degree_nesting_audit,
    certify_parity_audit,
    verify_degree_nesting_audit,
    verify_parity_audit,
)
from .weil_search import (
    WeilSearchCell,
    WeilSearchError,
    load_search_plan,
    run_weil_search,
    verify_weil_search,
)
from .weil_transition import (
    WeilTransitionError,
    load_transition_plan,
    run_weil_transition_search,
    verify_weil_transition_search,
)
from .weil_transition_summary import (
    WeilTransitionSummaryError,
    generate_weil_transition_summary,
    verify_weil_transition_summary,
)


def _integer_witness(value: str) -> tuple[int, ...]:
    """Parse a comma-separated exact integer vector for an audit command."""

    try:
        witness = tuple(int(item.strip()) for item in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "witness entries must be comma-separated integers"
        ) from exc
    if not witness or any(not item.strip() for item in value.split(",")):
        raise argparse.ArgumentTypeError(
            "witness entries must be comma-separated integers"
        )
    if not any(witness):
        raise argparse.ArgumentTypeError("witness must be nonzero")
    return witness


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rh-lab",
        description="Generate and validate certificate-first RH research artifacts.",
    )
    commands = parser.add_subparsers(dest="command", required=True)

    zeros = commands.add_parser("zeros", help="certify consecutive zeta zeros")
    zeros.add_argument("--first", type=int, default=1)
    zeros.add_argument("--count", type=int, default=100)
    zeros.add_argument("--bits", type=int, default=192)
    zeros.add_argument("--block-size", type=int, default=1000)
    zeros.add_argument("--output", type=Path, required=True)

    lagarias = commands.add_parser(
        "lagarias", help="certify a finite range of the Lagarias inequality"
    )
    lagarias.add_argument("--limit", type=int, default=10_000)
    lagarias.add_argument("--bits", type=int, default=192)
    lagarias.add_argument("--output", type=Path, required=True)

    weil = commands.add_parser(
        "weil", help="certify the frozen c=5/2, degree-4 finite Weil matrix"
    )
    weil.add_argument("--bits", type=int, default=192)
    weil.add_argument("--output", type=Path, required=True)

    weil_search = commands.add_parser(
        "weil-search",
        help="run or resume a checkpointed exploratory finite Weil grid",
    )
    weil_search.add_argument(
        "--plan", type=Path, default=Path("plans/weil-grid-v1.json")
    )
    weil_search.add_argument("--checkpoint-dir", type=Path, required=True)
    weil_search.add_argument("--resume", action="store_true")
    weil_search.add_argument(
        "--max-cells",
        type=int,
        help="complete at most this many new cells before checkpointing",
    )

    transition_search = commands.add_parser(
        "weil-transition-search",
        help="run or resume the exploratory prime-power transition batch",
    )
    transition_search.add_argument(
        "--plan",
        type=Path,
        default=Path("plans/weil-transition-q7-q9-v2.json"),
    )
    transition_search.add_argument(
        "--checkpoint-dir", type=Path, required=True
    )
    transition_search.add_argument("--resume", action="store_true")
    transition_search.add_argument(
        "--max-cells",
        type=int,
        help="complete at most this many new cells before checkpointing",
    )

    nyman_search = commands.add_parser(
        "nyman-search",
        help="run or resume the exploratory natural-dilate distance batch",
    )
    nyman_search.add_argument(
        "--plan",
        type=Path,
        default=Path("plans/nyman-natural-v1.json"),
    )
    nyman_search.add_argument("--checkpoint-dir", type=Path, required=True)
    nyman_search.add_argument("--resume", action="store_true")
    nyman_search.add_argument(
        "--max-cells",
        type=int,
        help="complete at most this many new cells before checkpointing",
    )

    nyman_normalization = commands.add_parser(
        "nyman-normalization-audit",
        help="generate the frozen exploratory Nyman normalization audit",
    )
    nyman_normalization.add_argument("--output", type=Path, required=True)

    nyman_summary = commands.add_parser(
        "summarize-nyman",
        help="derive a compact structural summary of the frozen Nyman batch",
    )
    nyman_summary.add_argument("--checkpoint-dir", type=Path, required=True)
    nyman_summary.add_argument("--output", type=Path, required=True)

    nyman_rebound = commands.add_parser(
        "nyman-rebound-audit",
        help="generate the frozen exploratory Nyman forced-rebound audit",
    )
    nyman_rebound.add_argument("--summary", type=Path, required=True)
    nyman_rebound.add_argument("--checkpoint-dir", type=Path, required=True)
    nyman_rebound.add_argument("--output", type=Path, required=True)

    nyman_mobius = commands.add_parser(
        "nyman-mobius-core-tail-audit",
        help="generate the frozen log-tapered Mobius core/tail audit",
    )
    nyman_mobius.add_argument("--summary", type=Path, required=True)
    nyman_mobius.add_argument("--checkpoint-dir", type=Path, required=True)
    nyman_mobius.add_argument("--output", type=Path, required=True)

    nyman_beta2_candidate = commands.add_parser(
        "propose-nyman-beta2-n512",
        help="propose the frozen exact N=512 dyadic witness (not a certificate)",
    )
    nyman_beta2_candidate.add_argument("--output", type=Path, required=True)

    nyman_beta2 = commands.add_parser(
        "nyman-beta2-n512-audit",
        help="certify the finite N=256 to N=512 contraction",
    )
    nyman_beta2.add_argument("--candidate", type=Path, required=True)
    nyman_beta2.add_argument("--summary", type=Path, required=True)
    nyman_beta2.add_argument("--checkpoint-dir", type=Path, required=True)
    nyman_beta2.add_argument("--output", type=Path, required=True)

    parity = commands.add_parser(
        "weil-parity-audit",
        help="generate an exploratory reversal-parity audit for one Weil cell",
    )
    parity.add_argument("--cutoff-numerator", type=int, required=True)
    parity.add_argument("--cutoff-denominator", type=int, default=1)
    parity.add_argument("--degree", type=int, required=True)
    parity.add_argument("--bits", type=int, default=192)
    parity.add_argument(
        "--witness",
        type=_integer_witness,
        action="append",
        default=[],
        help="optional comma-separated full-basis integer vector; repeatable",
    )
    parity.add_argument("--output", type=Path, required=True)

    nesting = commands.add_parser(
        "weil-nesting-audit",
        help="generate an exploratory fixed-cutoff degree-nesting audit",
    )
    nesting.add_argument("--cutoff-numerator", type=int, required=True)
    nesting.add_argument("--cutoff-denominator", type=int, default=1)
    nesting.add_argument("--lower-degree", type=int, required=True)
    nesting.add_argument("--higher-degree", type=int, required=True)
    nesting.add_argument("--bits", type=int, default=192)
    nesting.add_argument("--replay-bits", type=int, default=384)
    nesting.add_argument(
        "--witness",
        type=_integer_witness,
        action="append",
        default=[],
        help="optional comma-separated lower-basis integer vector; repeatable",
    )
    nesting.add_argument("--output", type=Path, required=True)

    transition_summary = commands.add_parser(
        "summarize-weil-transition",
        help="derive a compact non-numerical summary of a complete v2 batch",
    )
    transition_summary.add_argument(
        "--checkpoint-dir", type=Path, required=True
    )
    transition_summary.add_argument("--output", type=Path, required=True)

    verify_zeros = commands.add_parser(
        "verify-zeros", help="validate and replay a zero certificate"
    )
    verify_zeros.add_argument("--artifact", type=Path, required=True)
    verify_zeros.add_argument("--bits", type=int, default=384)

    verify_lagarias = commands.add_parser(
        "verify-lagarias", help="validate and replay a Lagarias certificate"
    )
    verify_lagarias.add_argument("--artifact", type=Path, required=True)

    verify_weil = commands.add_parser(
        "verify-weil", help="validate and replay a finite Weil certificate"
    )
    verify_weil.add_argument("--artifact", type=Path, required=True)
    verify_weil.add_argument("--bits", type=int, default=384)

    verify_weil_search_parser = commands.add_parser(
        "verify-weil-search",
        help="validate and replay a checkpointed exploratory Weil search",
    )
    verify_weil_search_parser.add_argument("--index", type=Path, required=True)
    verify_weil_search_parser.add_argument(
        "--checkpoint-dir", type=Path, required=True
    )

    verify_transition = commands.add_parser(
        "verify-weil-transition-search",
        help="validate and replay an exploratory prime-power transition batch",
    )
    verify_transition.add_argument("--index", type=Path, required=True)
    verify_transition.add_argument(
        "--checkpoint-dir", type=Path, required=True
    )

    verify_nyman = commands.add_parser(
        "verify-nyman-search",
        help="validate and replay an exploratory natural-dilate distance batch",
    )
    verify_nyman.add_argument("--index", type=Path, required=True)
    verify_nyman.add_argument("--checkpoint-dir", type=Path, required=True)
    verify_nyman.add_argument("--bits", type=int, default=1536)

    verify_nyman_normalization = commands.add_parser(
        "verify-nyman-normalization-audit",
        help="exactly regenerate a frozen Nyman normalization audit",
    )
    verify_nyman_normalization.add_argument(
        "--artifact", type=Path, required=True
    )

    verify_nyman_summary_parser = commands.add_parser(
        "verify-nyman-summary",
        help="regenerate a compact Nyman summary from its evidence tree",
    )
    verify_nyman_summary_parser.add_argument(
        "--summary", type=Path, required=True
    )
    verify_nyman_summary_parser.add_argument(
        "--checkpoint-dir", type=Path, required=True
    )

    verify_nyman_rebound = commands.add_parser(
        "verify-nyman-rebound-audit",
        help="exactly regenerate the frozen Nyman forced-rebound audit",
    )
    verify_nyman_rebound.add_argument("--artifact", type=Path, required=True)
    verify_nyman_rebound.add_argument("--summary", type=Path, required=True)
    verify_nyman_rebound.add_argument(
        "--checkpoint-dir", type=Path, required=True
    )

    verify_nyman_mobius = commands.add_parser(
        "verify-nyman-mobius-core-tail-audit",
        help="regenerate at 256 bits and replay the Mobius audit at 512 bits",
    )
    verify_nyman_mobius.add_argument("--artifact", type=Path, required=True)
    verify_nyman_mobius.add_argument("--summary", type=Path, required=True)
    verify_nyman_mobius.add_argument(
        "--checkpoint-dir", type=Path, required=True
    )

    verify_nyman_beta2 = commands.add_parser(
        "verify-nyman-beta2-n512-audit",
        help="regenerate and replay the finite N=512 contraction audit",
    )
    verify_nyman_beta2.add_argument("--artifact", type=Path, required=True)
    verify_nyman_beta2.add_argument("--candidate", type=Path, required=True)
    verify_nyman_beta2.add_argument("--summary", type=Path, required=True)
    verify_nyman_beta2.add_argument(
        "--checkpoint-dir", type=Path, required=True
    )
    verify_nyman_beta2.add_argument(
        "--bits", type=int, default=NYMAN_BETA2_REPLAY_BITS
    )

    verify_parity = commands.add_parser(
        "verify-weil-parity-audit",
        help="validate and replay an exploratory Weil parity audit",
    )
    verify_parity.add_argument("--artifact", type=Path, required=True)

    verify_nesting = commands.add_parser(
        "verify-weil-nesting-audit",
        help="validate and replay an exploratory Weil degree-nesting audit",
    )
    verify_nesting.add_argument("--artifact", type=Path, required=True)

    verify_summary = commands.add_parser(
        "verify-weil-transition-summary",
        help="regenerate a compact v2 summary from its evidence tree",
    )
    verify_summary.add_argument("--summary", type=Path, required=True)
    verify_summary.add_argument(
        "--checkpoint-dir", type=Path, required=True
    )

    claims = commands.add_parser("verify-claims", help="validate the claim ledger")
    claims.add_argument(
        "--registry", type=Path, default=Path("claims/registry.json")
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)

    if args.command == "zeros":
        artifact = certify_critical_line_zeros(
            first_index=args.first,
            count=args.count,
            precision_bits=args.bits,
            count_block_size=args.block_size,
        )
        write_json(args.output, artifact)
        print(json.dumps({
            "classification": artifact["classification"],
            "output": str(args.output),
            "payload_sha256": artifact["payload_sha256"],
        }, sort_keys=True))
        return 0 if artifact["classification"] == "CERTIFIED_FINITE" else 2

    if args.command == "lagarias":
        artifact = certify_lagarias_range(args.limit, args.bits)
        write_json(args.output, artifact)
        print(json.dumps({
            "classification": artifact["classification"],
            "output": str(args.output),
            "payload_sha256": artifact["payload_sha256"],
        }, sort_keys=True))
        return 0 if artifact["classification"] == "CERTIFIED_FINITE" else 2

    if args.command == "weil":
        artifact = certify_weil_matrix(precision_bits=args.bits)
        write_json(args.output, artifact)
        print(json.dumps({
            "classification": artifact["classification"],
            "output": str(args.output),
            "payload_sha256": artifact["payload_sha256"],
        }, sort_keys=True))
        return 0 if artifact["classification"] == "CERTIFIED_FINITE" else 2

    if args.command == "weil-search":
        try:
            plan = load_search_plan(args.plan)
            index = run_weil_search(
                plan,
                args.checkpoint_dir,
                resume=args.resume,
                max_cells=args.max_cells,
            )
        except (WeilSearchError, OSError, ValueError, TypeError) as exc:
            print(f"Weil search rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": index["classification"],
                    "conclusion": index["conclusion"],
                    "index": str(args.checkpoint_dir / "index.json"),
                    "payload_sha256": index["payload_sha256"],
                    "progress": index["progress"],
                },
                sort_keys=True,
            )
        )
        return 0

    if args.command == "weil-transition-search":
        try:
            plan = load_transition_plan(args.plan)
            index = run_weil_transition_search(
                plan,
                args.checkpoint_dir,
                resume=args.resume,
                max_cells=args.max_cells,
            )
        except (WeilTransitionError, OSError, TypeError, ValueError) as exc:
            print(f"Weil transition search rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": index["classification"],
                    "conclusion": index["conclusion"],
                    "hypothesis_status": index["hypothesis_status"],
                    "index": str(args.checkpoint_dir / "index.json"),
                    "payload_sha256": index["payload_sha256"],
                    "progress": index["progress"],
                },
                sort_keys=True,
            )
        )
        return 0

    if args.command == "nyman-search":
        try:
            plan = load_nyman_plan(args.plan)
            index = run_nyman_search(
                plan,
                args.checkpoint_dir,
                resume=args.resume,
                max_cells=args.max_cells,
            )
        except (NymanSearchError, OSError, TypeError, ValueError) as exc:
            print(f"Nyman search rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": index["classification"],
                    "conclusion": index["conclusion"],
                    "hypothesis_status": index["hypothesis_status"],
                    "index": str(args.checkpoint_dir / "index.json"),
                    "payload_sha256": index["payload_sha256"],
                    "progress": index["progress"],
                },
                sort_keys=True,
            )
        )
        return 0

    if args.command == "nyman-normalization-audit":
        try:
            artifact = generate_nyman_normalization_bundle()
            write_json(args.output, artifact)
        except (NymanNormalizationError, OSError, TypeError, ValueError) as exc:
            print(f"Nyman normalization audit rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": artifact["classification"],
                    "audit_outcome": artifact["audit_outcome"],
                    "hypothesis_status": artifact["hypothesis_status"],
                    "output": str(args.output),
                    "payload_sha256": artifact["payload_sha256"],
                },
                sort_keys=True,
            )
        )
        return 0 if artifact["audit_outcome"] == (
            "NORMALIZATION_AUDIT_PASSED_EXPLORATORY"
        ) else 2

    if args.command == "summarize-nyman":
        try:
            artifact = generate_nyman_summary(args.checkpoint_dir)
            write_json(args.output, artifact)
        except (NymanSummaryError, OSError, TypeError, ValueError) as exc:
            print(f"Nyman summary rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": artifact["classification"],
                    "conclusion": artifact["conclusion"],
                    "hypothesis_status": artifact["hypothesis_status"],
                    "output": str(args.output),
                    "payload_sha256": artifact["payload_sha256"],
                    "counts": artifact["counts"],
                },
                sort_keys=True,
            )
        )
        return 0

    if args.command == "nyman-rebound-audit":
        try:
            artifact = generate_nyman_rebound_audit(
                args.summary,
                args.checkpoint_dir,
            )
            write_json(args.output, artifact)
        except (NymanReboundError, OSError, TypeError, ValueError) as exc:
            print(f"Nyman rebound audit rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": artifact["classification"],
                    "audit_outcome": artifact["audit_outcome"],
                    "hypothesis_status": artifact["hypothesis_status"],
                    "output": str(args.output),
                    "payload_sha256": artifact["payload_sha256"],
                    "scaled_declines": artifact[
                        "five_exact_scaled_declines"
                    ]["count"],
                },
                sort_keys=True,
            )
        )
        return 0

    if args.command == "nyman-mobius-core-tail-audit":
        try:
            artifact = generate_nyman_mobius_audit(
                args.summary,
                args.checkpoint_dir,
            )
            write_json(args.output, artifact)
        except (NymanMobiusError, OSError, TypeError, ValueError) as exc:
            print(f"Nyman Mobius audit rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": artifact["classification"],
                    "audit_outcome": artifact["audit_outcome"],
                    "hypothesis_status": artifact["hypothesis_status"],
                    "output": str(args.output),
                    "payload_sha256": artifact["payload_sha256"],
                    "certified_cells": str(
                        len(artifact["generation"]["records"])
                    ),
                },
                sort_keys=True,
            )
        )
        return 0

    if args.command == "propose-nyman-beta2-n512":
        try:
            artifact = propose_nyman_beta2_candidate()
            write_json(args.output, artifact)
        except (NymanBeta2Error, OSError, TypeError, ValueError) as exc:
            print(f"Nyman beta=2 candidate rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": artifact["classification"],
                    "hypothesis_status": artifact["hypothesis_status"],
                    "output": str(args.output),
                    "payload_sha256": artifact["payload_sha256"],
                    "solver_role": artifact["solver_role"],
                },
                sort_keys=True,
            )
        )
        return 0

    if args.command == "nyman-beta2-n512-audit":
        try:
            artifact = generate_nyman_beta2_audit(
                args.candidate,
                args.summary,
                args.checkpoint_dir,
            )
            write_json(args.output, artifact)
        except (NymanBeta2Error, OSError, TypeError, ValueError) as exc:
            print(f"Nyman beta=2 audit rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": artifact["classification"],
                    "audit_outcome": artifact["audit_outcome"],
                    "hypothesis_status": artifact["hypothesis_status"],
                    "output": str(args.output),
                    "payload_sha256": artifact["payload_sha256"],
                },
                sort_keys=True,
            )
        )
        return 0

    if args.command == "weil-parity-audit":
        try:
            artifact = certify_parity_audit(
                WeilSearchCell(
                    args.cutoff_numerator,
                    args.cutoff_denominator,
                    args.degree,
                ),
                args.bits,
                witnesses=args.witness,
            )
            write_json(args.output, artifact)
        except (OSError, TypeError, ValueError) as exc:
            print(f"Weil parity audit rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": artifact["classification"],
                    "decision": artifact["decision"],
                    "hypothesis_status": artifact["hypothesis_status"],
                    "output": str(args.output),
                    "payload_sha256": artifact["payload_sha256"],
                },
                sort_keys=True,
            )
        )
        return 0 if artifact["decision"] != "AUDIT_FAILED" else 2

    if args.command == "weil-nesting-audit":
        try:
            lower = WeilSearchCell(
                args.cutoff_numerator,
                args.cutoff_denominator,
                args.lower_degree,
            )
            higher = WeilSearchCell(
                args.cutoff_numerator,
                args.cutoff_denominator,
                args.higher_degree,
            )
            artifact = certify_degree_nesting_audit(
                lower,
                higher,
                args.bits,
                args.replay_bits,
                negative_witnesses=args.witness,
            )
            write_json(args.output, artifact)
        except (OSError, TypeError, ValueError) as exc:
            print(f"Weil nesting audit rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": artifact["classification"],
                    "decision": artifact["decision"],
                    "hypothesis_status": artifact["hypothesis_status"],
                    "output": str(args.output),
                    "payload_sha256": artifact["payload_sha256"],
                },
                sort_keys=True,
            )
        )
        return 0 if artifact["decision"] != "AUDIT_FAILED" else 2

    if args.command == "summarize-weil-transition":
        try:
            artifact = generate_weil_transition_summary(args.checkpoint_dir)
            write_json(args.output, artifact)
        except (WeilTransitionSummaryError, OSError) as exc:
            print(f"Weil transition summary rejected: {exc}")
            return 2
        print(
            json.dumps(
                {
                    "classification": artifact["classification"],
                    "conclusion": artifact["conclusion"],
                    "hypothesis_status": artifact["hypothesis_status"],
                    "output": str(args.output),
                    "payload_sha256": artifact["payload_sha256"],
                },
                sort_keys=True,
            )
        )
        return 0

    if args.command == "verify-zeros":
        try:
            artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
            result = verify_zero_certificate(artifact, args.bits)
        except (ZeroCertificateError, OSError, json.JSONDecodeError) as exc:
            print(f"zero certificate rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-lagarias":
        try:
            artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
            result = verify_lagarias_certificate(artifact)
        except (LagariasCertificateError, OSError, json.JSONDecodeError) as exc:
            print(f"Lagarias certificate rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-weil":
        try:
            artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
            result = verify_weil_certificate(artifact, args.bits)
        except (WeilCertificateError, OSError, json.JSONDecodeError) as exc:
            print(f"Weil certificate rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-weil-search":
        try:
            result = verify_weil_search(args.index, args.checkpoint_dir)
        except (WeilSearchError, OSError, json.JSONDecodeError) as exc:
            print(f"Weil search artifact rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-weil-transition-search":
        try:
            result = verify_weil_transition_search(
                args.index,
                args.checkpoint_dir,
            )
        except (WeilTransitionError, OSError, json.JSONDecodeError) as exc:
            print(f"Weil transition search rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-nyman-search":
        try:
            result = verify_nyman_search(
                args.index,
                args.checkpoint_dir,
                replay_precision_bits=args.bits,
            )
        except (NymanSearchError, OSError, json.JSONDecodeError) as exc:
            print(f"Nyman search artifact rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-nyman-normalization-audit":
        try:
            result = verify_nyman_normalization_bundle(args.artifact)
        except (
            NymanNormalizationError,
            OSError,
            json.JSONDecodeError,
        ) as exc:
            print(f"Nyman normalization audit rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-nyman-summary":
        try:
            result = verify_nyman_summary(
                args.summary,
                args.checkpoint_dir,
            )
        except (
            NymanSummaryError,
            OSError,
            json.JSONDecodeError,
        ) as exc:
            print(f"Nyman summary rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-nyman-rebound-audit":
        try:
            result = verify_nyman_rebound_audit(
                args.artifact,
                args.summary,
                args.checkpoint_dir,
            )
        except (
            NymanReboundError,
            OSError,
            json.JSONDecodeError,
        ) as exc:
            print(f"Nyman rebound audit rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-nyman-mobius-core-tail-audit":
        try:
            result = verify_nyman_mobius_audit(
                args.artifact,
                args.summary,
                args.checkpoint_dir,
            )
        except (
            NymanMobiusError,
            OSError,
            json.JSONDecodeError,
        ) as exc:
            print(f"Nyman Mobius audit rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-nyman-beta2-n512-audit":
        try:
            result = verify_nyman_beta2_audit(
                args.artifact,
                args.candidate,
                args.summary,
                args.checkpoint_dir,
                replay_precision_bits=args.bits,
            )
        except (
            NymanBeta2Error,
            OSError,
            json.JSONDecodeError,
        ) as exc:
            print(f"Nyman beta=2 audit rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-weil-parity-audit":
        try:
            artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
            result = verify_parity_audit(artifact)
        except (
            WeilAuditVerificationError,
            OSError,
            json.JSONDecodeError,
        ) as exc:
            print(f"Weil parity audit rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-weil-nesting-audit":
        try:
            artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
            result = verify_degree_nesting_audit(artifact)
        except (
            WeilAuditVerificationError,
            OSError,
            json.JSONDecodeError,
        ) as exc:
            print(f"Weil nesting audit rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "verify-weil-transition-summary":
        try:
            result = verify_weil_transition_summary(
                args.summary,
                args.checkpoint_dir,
            )
        except (WeilTransitionSummaryError, OSError, json.JSONDecodeError) as exc:
            print(f"Weil transition summary rejected: {exc}")
            return 2
        print(json.dumps(result, sort_keys=True))
        return 0

    try:
        registry = load_and_validate_registry(args.registry)
    except (ClaimValidationError, OSError, json.JSONDecodeError) as exc:
        print(f"claim ledger rejected: {exc}")
        return 2
    print(
        json.dumps(
            {
                "claims": len(registry["claims"]),
                "hypothesis_status": registry["hypothesis"]["status"],
                "registry": str(args.registry),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
