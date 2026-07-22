# Contributing

This project welcomes mathematics, implementations, failed approaches, test
cases, and adversarial audits. It does not welcome unsupported solution claims.

Every contribution must assign each result exactly one classification:

- `PROVED`: a global theorem with all assumptions and infinite limiting steps
  discharged;
- `DISPROVED`: a rigorous counterexample to RH or to a proved equivalent;
- `CERTIFIED_FINITE`: a finite statement checked with exact or enclosing
  arithmetic;
- `CONDITIONAL`: a theorem whose unproved assumptions are listed;
- `EXPLORATORY`: numerical or symbolic evidence without a certificate; or
- `FAILED`: a retained null result, bug, or invalidated idea.

Pull requests that touch a headline claim must include the artifact, its
reproduction command, an assumption manifest, and an independent audit. A
finite calculation must state its finite scope and explicitly say that it does
not prove RH.

Before opening a pull request, run:

```console
python -m pip install -e ".[dev]"
pytest
rh-lab verify-claims
```

For a proposed proof, include an axiom audit (`#print axioms` in Lean or its
equivalent) and remove every placeholder such as `sorry`, `admit`, or an
unproved custom axiom. For a proposed off-line zero, include an interval box
disjoint from `Re(s)=1/2`, a rigorous argument-principle count, and an
independent reproduction using a separately implemented verifier.
