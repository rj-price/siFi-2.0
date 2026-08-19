# End-to-end baseline

`design_TUB3_fragment.json` and `offtarget_TUB3_fragment.json` are the port's own whole-pipeline output for the
committed query, pinned by `tests/test_end_to_end.py`. They are *not* golden fixtures in the sense of
`tests/golden/`: those were captured from the running Python 2 original and carry the fidelity guarantee, while
these were produced by the port and only pin it against itself. A full original-vs-port comparison was never
available — the original's design mode blocks on a Qt dialog and it shipped Windows `.exe` binaries.

They exist so a whole-pipeline change is visible as a diff. `PLAN.md` Phase 6 moves these numbers
deliberately, one attributable commit at a time; each such commit regenerates both files with the commands
below.

Regenerate (from the repo root, with a scratch directory `$TMP`):

```sh
conda run -n sifi2 sifi db build --fasta tests/data/reference.fasta --name testdb --db-location "$TMP/dbs"
conda run -n sifi2 sifi offtarget --query tests/data/query.fasta --db testdb \
    --db-location "$TMP/dbs" --outdir "$TMP/offtarget" -q
conda run -n sifi2 sifi design --query tests/data/query.fasta --db testdb --all-targets-main \
    --db-location "$TMP/dbs" --outdir "$TMP/design" -q
cp "$TMP/offtarget/TUB3_fragment.json" tests/baseline/offtarget_TUB3_fragment.json
cp "$TMP/design/TUB3_fragment.json"    tests/baseline/design_TUB3_fragment.json
```

Everything here runs at the CLI defaults (`--mismatches 0`, `--sirna-size 21`), which is why the record counts
differ from `tests/data/pipeline_design.json` — that one is real Python 2 output captured at `--mismatches 2`.
