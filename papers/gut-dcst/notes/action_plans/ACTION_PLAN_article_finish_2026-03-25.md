# Action Plan: Gut Application Paper Finishing Pass

Date: 2026-03-25
Status: Completed

## Objective

Finish the current preprint draft by converting the remaining planned items into
reproducible manuscript assets and integrating them into the paper.

## Workstreams

1. Build the missing main-text tables.
   - Table 1: top depth-1 DCSTs in the full AGP cohort.
   - Table 2: strongest adjusted depth-1 association for each phenotype with at
     least one adjusted-significant signal.
   - Table 3: disease-wise comparison of univariate, adjusted, and
     contamination-aware counts, together with representative retained and
     attenuated examples.
   - Table 4: SILVA versus GG2 agreement summary, including both overall
     agreement metrics and representative remappings.

2. Integrate the tables into the manuscript.
   - Place the shared-overview tables in Section 3.
   - Add the external-validation table to Section 5.
   - Update surrounding prose so the tables are explained and not merely
     inserted.

3. Calibrate claims across disease sections.
   - Keep IBD, IBS, seasonal allergy, CDI, and kidney disease as the strongest
     disease-wise sections.
   - Soften umbrella or weaker sections, especially autoimmune disease, lung
     disease, migraine, acid reflux, diabetes, and CVD.
   - Make the distinction between discovery-screen phenotypes and
     disease-specific mechanistic literature more explicit.

4. Add standard preprint back matter.
   - Data and code availability.
   - Funding.
   - Competing interests.
   - Acknowledgments.

5. Rebuild and visually verify the PDF draft.
   - Render the updated manuscript.
   - Check the title page, the new tables, and the back-matter pages for layout
     issues.

## Deliverables

- Updated manuscript source:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`
- Updated manuscript PDF:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf`
- Reproducible table builder:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/scripts/build_manuscript_tables.py`
- Table assets:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/TABLE_1_top_depth1_dcsts.*`
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/TABLE_2_headline_adjusted_associations.*`
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/TABLE_3_univariate_vs_adjusted_summary.*`
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/TABLE_4_taxonomy_concordance_summary.*`

## Execution Notes

- Built the missing manuscript tables with
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/scripts/build_manuscript_tables.py`.
- Inserted Tables 1-5 into the draft manuscript and updated the disease-wise
  prose to soften umbrella or weaker claims.
- Added Data and Code Availability, Funding, Competing Interests, and
  Acknowledgments.
- Rebuilt and visually checked the PDF, including the new overview tables, the
  compact validation table, and the back-matter pages.
