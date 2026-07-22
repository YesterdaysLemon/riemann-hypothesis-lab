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
from .weil_search import (
    WeilSearchError,
    load_search_plan,
    run_weil_search,
    verify_weil_search,
)


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
