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

    verify_zeros = commands.add_parser(
        "verify-zeros", help="validate and replay a zero certificate"
    )
    verify_zeros.add_argument("--artifact", type=Path, required=True)
    verify_zeros.add_argument("--bits", type=int, default=384)

    verify_lagarias = commands.add_parser(
        "verify-lagarias", help="validate and replay a Lagarias certificate"
    )
    verify_lagarias.add_argument("--artifact", type=Path, required=True)

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
