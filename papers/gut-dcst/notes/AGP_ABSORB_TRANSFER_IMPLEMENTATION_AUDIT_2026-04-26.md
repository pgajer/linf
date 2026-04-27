# AGP absorb transfer implementation audit

Date: 2026-04-26

## Scope

This audit checks the **logic of the corrected frozen-transfer implementation**
in:

`/Users/pgajer/current_projects/gut_microbiome/scripts/run_external_case_control_validation_frozen_agp.R`

It does **not** reinterpret the previously generated transfer outputs, because
those outputs were produced under the older copy-forward rule.

## Correction made

The transfer reference is now built from the AGP assignment table using the
full realized child labels at each nominal depth, rather than using only the
suffix components and then copying the parent label across all remaining
depths.

As a result:

- repeated labels are allowed only when the AGP assignment table itself
  realizes that repeated label at the deeper nominal depth;
- a lineage can stay flat at depth `d` and then split again at depth `d + 1`;
- the script no longer freezes all later depths once it encounters a stable
  lineage state.

## Logic checks

### Check 1: stable depth-2 lineage can split again at depth 3

Toy AGP hierarchy:

- depth 1: `A`
- depth 2: `A`
- depth 3: `A > X` or `A > Y`

Toy sample support:

- `A = 10`
- `X = 3`
- `Y = 2`

Expected behavior:

- depth 1: `A`
- depth 2: `A`
- depth 3: `A > X`
- depth 4: `A > X`

Observed from the corrected implementation:

`A | A | A__X | A__X`

Interpretation:

- the stable depth-2 lineage was preserved;
- the later split at depth 3 remained reachable;
- this is the behavior the previous copy-forward implementation could not
  represent correctly.

### Check 2: parent cannot repeat if AGP never realizes that repetition

Toy AGP hierarchy:

- depth 1: `P`
- depth 2: `P > Q1` or `P > Q2`
- no realized `depth2 = P`

Toy sample support:

- `P = 10`
- `Q1 = 0`
- `Q2 = 0`

Expected behavior:

- depth 1: `P`
- depth 2: `NA`

Observed from the corrected implementation:

`P | NA | NA | NA`

Interpretation:

- the script does **not** fabricate a repeated parent label at depth 2 when AGP
  never realizes that repeated parent there.

### Check 3: real Halfvarson sample now follows the intended stable-then-split rule

Validation sample:

- `ERR1746351`

Under the previous procedural copy-forward rule, this sample was effectively
frozen as:

- `Gemmiger sp.`
- `Gemmiger sp.`
- `Gemmiger sp.`
- `Gemmiger sp.`

The corrected implementation assigns:

- `Gemmiger sp.`
- `Gemmiger sp.`
- `Gemmiger sp.__Faecalibacterium prausnitzii`
- `Gemmiger sp.__Faecalibacterium prausnitzii__Agathobacter faecis`

Interpretation:

- `Gemmiger sp.` is a legitimate stable lineage at depth 2 in the AGP absorb
  reference;
- the lineage then splits again at depth 3;
- the corrected transfer logic now reproduces that intended hierarchy behavior.

## Conclusion

The corrected implementation now matches the intended AGP absorb transfer rule:

- validation samples are compared against realized AGP lineage sets at each
  nominal depth;
- repeated labels are valid only when those repeated labels are realized by AGP
  at that depth;
- stable lineage states can persist for one or more depths and still split
  later if the AGP hierarchy realizes deeper descendants.
