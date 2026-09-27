# lava plugin

Migrated from the legacy `lava_container` (single-locus univ/bivar/pcor/multireg)
and `lava_scan_container` (multiple-locus scan) wrappers in nodes-io. One
directory = one plugin family = one git-able unit. Both nodes share one image
and one panel binding, so they are one plugin with two `[[nodes]]` entries
(the migration doc's family rule).

## Layout

- `manifest.toml` — node kinds `lava` and `lava_scan`: params, ports, panel,
  image provenance
- `scripts/lava.sh` — single-locus Rscript program (`scripts/*.sh` is the
  family convention; the content is R, executed as
  `Rscript /work/.autonomics/script`)
- `scripts/lava_scan.sh` — multiple-locus scan Rscript program
- `Dockerfile` — image build recipe (moved verbatim from
  `containers/lava/`; the image installs the official LAVA 0.1.5 package
  from the upstream GitHub tarball at the pinned commit)
- `LAVA/` — vendored upstream checkout (reference material and the
  `vignettes/data` fixtures used by live integration tests). The embedded
  `.git` was stripped on migration so the plugin commits it as plain files
  (migration doc pitfall 7); upstream history lives at
  github.com/josefin-werme/LAVA.

## Install

```sh
export AUTONOMICS_PLUGIN_ROOT=/mnt/projects/node-plugins
cargo test -p container-plugin --test lava_migration   # golden parity
```

After publishing, declare the source in `~/.autonomics/plugins.toml` with a
pinned `rev`, exactly like the other wave families.

## Panel decision: one static panel, tutorial variant is a DSL gap

The legacy wrappers carried a `panel_id` param selecting between two catalog
packages, both mounted at `/panels/lava_ref` but with different filename
prefixes:

| panel | ref.prefix the legacy script rendered |
|---|---|
| `wjixiang/catalog-lava-ref-ukb-eur` (default, production) | `/panels/lava_ref/lava-ukb-v1.1` |
| `wjixiang/catalog-lava-ref-1000g-test` (tutorial, migration validation only) | `/panels/lava_ref/g1000_test` |

The v0 plugin DSL has **no parameter-driven panel selection** (family panels
are static; `from_param` panels are a recorded future extension), and a node
cannot express "either panel" any other way. Mounting both panels to preserve
run-time choice is impossible for this family: the loader rejects two panels
on one mount path, and repathing one panel would change the `ref.prefix`
string the script passes to `process.input` anyway. Splitting the family into
two plugins per panel would contradict the share-panels family rule without
helping, because the choice is per *spec*, not per node.

Chosen design: bind the **default panel statically**
(`wjixiang/catalog-lava-ref-ukb-eur` at `/panels/lava_ref`), which is what
both legacy factories advertised through `data_bundles()` and what the
backlog records as the production binding. The `panel_id`,
`mount_path_override`, and `ref_prefix_override` params are dropped; the
script hard-codes `ref.prefix = "/panels/lava_ref/lava-ukb-v1.1"`, the exact
string the legacy wrapper rendered for this panel. The tutorial panel variant
(`catalog-lava-ref-1000g-test` / `g1000_test`), which existed for migration
regression against the official vignette, returns when `from_param` panel
selection lands; until then a regression run can point
`AUTONOMICS_LAVA_IT_SOURCE` at the vendored vignettes and use the legacy-era
image contract, or temporarily edit the single `[[panels]]` entry.

## Port decision: the dominant legacy shape, statically

The legacy factories computed ports per spec:
`port_layout_for(sample_overlap, phenotypes.len())` — input.info, loci, an
optional sample-overlap port, then one sumstats port per phenotype. The v0
plugin DSL has **no dynamic ports**: every declared input port is required,
and unconnected declared ports fail DAG validation. A superset layout is
therefore not "compatible with more specs" — it would force dummy files onto
ports a given run does not use, because the manifest DSL cannot mark input
ports optional.

Chosen design: bless the **dominant shape** both wrappers' own tests, the
live integration test, and the official vignette all use —
`sample_overlap = true` with **two phenotypes**:

- ports 0-1: input.info and loci tables (fixed in every legacy shape)
- port 2: sample.overlap table (always required now)
- ports 3-4: sumstats for `phenotypes[1]` and `phenotypes[2]`
- `phenotypes` is pinned with `min_len = max_len = 2` so the param list and
  the port count cannot drift apart

Outputs keep the legacy count for each kind: three for `lava`
(`lava.tsv`, `lava.RDS`, `lava.log`), four for `lava_scan` (`univ`/`bivar`
TSVs, RDS, log).

Consequences (recorded as the dynamic-ports gap, to revisit when the DSL
grows optional/variadic input ports):

- single-phenotype `univ` runs are not expressible (wire a second phenotype —
  `run.univ` then reports both — or wait for optional ports);
- overlap-free runs (`sample_overlap = false` in the legacy spec) are not
  expressible;
- scans over 3+ phenotypes are not expressible (`run.univ.bivar` needs two;
  larger scans wait for variadic ports).

## Script parity

Both scripts are the legacy wrapper's R programs with every parameter moved
from Rust string-building into environment variables (the mvmr/ldsc plugin
conventions): booleans render `"true"`/`"false"` and convert with
`as.logical`, string arrays render space-joined and rebuild with
`strsplit(..., fixed = TRUE)`, optional params render `""` and test with
`nzchar`. The official LAVA API tokens are identical to the legacy scripts:
`process.input -> read.loci -> process.locus -> run.*` on the single-locus
node, and one `process.input` plus a `seq_len(nrow(loci))` loop with
per-locus `tryCatch` on the scan node. The legacy `validate()` rules the v0
param DSL cannot express (enum membership for `analysis`, target-shape per
analysis, duplicate/empty ID rejection) run inside the container with the
legacy error messages; the numeric and array-length bounds ARE manifest
bounds and fail at compile time. `param.lim` stays 1.25, or 1.5 for
`run.multireg`, selected in the script exactly like the legacy
`container_spec()`.

## Golden test

`crates/container-plugin/tests/lava_migration.rs` compiles both nodes and
asserts field-for-field parity with the legacy `container_spec` output:
image, outputs, formats, panel bundle, network, read-only rootfs, pull
policy, timeouts, env defaults, and semantic script markers. Known deltas,
each deliberate:

- artifact prefixes follow the kind rename
  (`/artifacts/lava_container` -> `/artifacts/lava`,
  `/artifacts/lava_scan_container` -> `/artifacts/lava_scan`), matching the
  ldsc wave precedent;
- `panel_bundles` is the static default panel instead of the spec-selected
  one (panel gap above);
- analysis dispatch lives in one script switch instead of four Rust-built
  scripts (semantic, not byte, parity);
- the input-count guard expects 5 (the static contract) instead of
  `2 + sample_overlap + phenotypes.len()`.
