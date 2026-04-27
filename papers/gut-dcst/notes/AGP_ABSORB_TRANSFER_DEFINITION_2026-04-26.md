# AGP absorb transfer definition

Date: 2026-04-26

## Purpose

This note defines the intended semantics of AGP absorb lineage-set transfer to
external validation samples. It is meant to replace the looser "carry-forward"
language used in earlier notes.

## Objects

Let:

- `L_d` be the set of retained AGP absorb lineage sets at nominal depth `d`;
- `a_d(x)` be the AGP absorb assignment of AGP sample `x` at depth `d`;
- `s` be an external validation sample after taxon-label harmonization into the
  AGP label space.

The authoritative retained lineage sets are the distinct values already present
in the AGP assignment table columns:

- `depth1_absorb`
- `depth2_absorb`
- `depth3_absorb`
- `depth4_absorb`

for the AGP run used as the frozen reference.

## Valid lineage sets at each nominal depth

A lineage set is valid at depth `d` if and only if it is realized by at least
one AGP sample in `depthd_absorb`.

This definition is crucial because under the absorb policy a lineage can remain
stable across adjacent depths. For example, if some AGP samples satisfy:

- `depth1_absorb = s1`
- `depth2_absorb = s1`
- `depth3_absorb = s1 > s2`

then `s1` is a legitimate retained depth-2 lineage set, even though it is only
length 1 as a taxonomic chain.

More generally:

- a lineage may remain unchanged across multiple nominal depths;
- if the AGP assignment table repeats the same label at depth `d`, that label
  is a valid depth-`d` lineage set;
- such repetition is part of the absorb hierarchy itself, not a transfer
  artifact.

## Parent-child consistency

For `d > 1`, a candidate depth-`d` lineage set for an external sample is valid
only if it is a realized AGP child of the sample's assigned depth-`d-1`
lineage set.

Operationally, the frozen transfer reference is a depth-indexed parent-child
graph built directly from adjacent AGP assignment columns:

- parent column: `depth(d-1)_absorb`
- child column: `depthd_absorb`

The child may be:

- a more specific lineage set, such as `s1 > s2`, or
- the same lineage repeated, such as `s1`.

Both are legitimate if they occur in the AGP assignment table.

## Transfer rule

For each validation sample `s` and each nominal depth `d`, assign `s` to the
closest valid AGP absorb lineage set reachable at that depth.

In the current frozen-transfer implementation, "closest" is operationalized
greedily within the AGP parent-child graph:

1. At depth 1, compare all realized AGP depth-1 lineage sets.
2. At depth `d > 1`, compare only the realized AGP children of the previously
   assigned depth-`d-1` lineage set.
3. Score each candidate child by the abundance of its last lineage component in
   the validation sample after AGP taxon harmonization.
4. Break ties by the child support observed in the AGP reference.

This means that a repeated label across depths is valid only when that repeated
label is itself a realized AGP lineage set at the deeper nominal depth.

## What is not allowed

The transfer should **not** fabricate a deeper nominal assignment simply by
copying a shallower parent label when that repeated label is not realized by
AGP at the deeper depth.

So, if AGP contains:

- depth 1: `p`
- depth 2 children: `p > q1`, `p > q2`

but never `depth2_absorb = p`, then a validation sample assigned to `p` at
depth 1 must **not** be reported as `p` again at depth 2.

In that situation the sample either:

- matches one of the valid depth-2 children, or
- remains unassigned at that nominal depth under the current greedy rule.

## Correct interpretation of repeated labels

A repeated transferred label such as:

- depth 1: `Blautia sp.`
- depth 2: `Blautia sp.`
- depth 3: `Blautia sp.`

is legitimate if and only if:

- `Blautia sp.` is realized by AGP in `depth2_absorb` and `depth3_absorb`.

This is an absorb-consistent stable lineage, not a sample-loss artifact.

## Suggested manuscript wording

Suggested concise wording for paper-facing use:

> Frozen AGP transfer assigned each validation sample to the closest retained
> AGP absorb lineage set at each nominal depth. Because absorb lineages may
> remain stable across adjacent depths, repeated labels at deeper levels were
> treated as valid only when those same repeated lineage sets were realized in
> the AGP reference hierarchy.
