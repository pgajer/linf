# AI mSystems Review Prompt System for the dCST Methods Paper

This prompt system is for reviewing the combined methods-focused dCST paper:

- `/Users/pgajer/current_projects/linf/papers/dcst-methods/build/dcst_methods_paper.pdf`

against the main source papers it draws from:

- `/Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf`
- `/Users/pgajer/current_projects/CT_clearance/docs/ct_dcst_manuscript.pdf`
- `/Users/pgajer/current_projects/Linf_paper_arXiv_submission/Linf_paper.pdf`

The main submission target is `mSystems`, so the review lens should be:

- microbial systems biology
- method contribution
- cross-application usefulness
- interpretability
- transportability / portability
- validation and reproducibility

rather than primarily:

- disease-specialist novelty
- clinical biomarker readiness
- mechanistic causality claims

## Core Goal

The paper should be judged mainly as:

- a dCST framework paper
- a systems-microbiology methods paper
- a cross-application representation paper

with gut and CT acting as supporting applications rather than the manuscript
trying to stand or fall as a purely IBD or CT disease paper.

## Best Use

Use this in stages:

1. overall `mSystems` fit and reviewer read
2. method-contribution audit
3. source-selection and cross-paper synthesis audit
4. AI-note / prose-naturalness audit
5. flow and readability audit
6. overclaim and framing audit
7. final submission-readiness priorities

## Core Reviewer Frame

Use this frame before any of the prompts below if you want tighter output:

```text
Act as an expert microbial systems biology reviewer for mSystems.

Assume you are experienced in:
- microbiome methods
- microbial community state representation
- computational microbiology
- cross-cohort and cross-application validation
- scientific writing and reviewer standards for systems-biology journals

Your job is not to be agreeable. Your job is to determine:
1. whether this manuscript is a good fit for mSystems
2. whether the dCST framework contribution is clear and substantial
3. whether the chosen applications support the framework well
4. whether the manuscript reads like a polished human-written research paper
5. what major and minor issues need to be fixed before submission

Important constraints:
- Distinguish clearly between direct observations from the manuscript and your own inference.
- Do not invent missing analyses, results, or literature.
- Cite the relevant page, section, figure, table, or supplement item whenever possible.
- When comparing the combined methods paper to the source papers, point to the relevant section or figure in the source paper as well.
- Separate major concerns from minor comments.
- If something is fixable, explain how.
- If something is strong, explain why.
```

## Prompt 1: Overall mSystems Reviewer Read

Use this first.

```text
Please read these documents carefully:

- /Users/pgajer/current_projects/linf/papers/dcst-methods/build/dcst_methods_paper.pdf
- /Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf
- /Users/pgajer/current_projects/CT_clearance/docs/ct_dcst_manuscript.pdf
- /Users/pgajer/current_projects/Linf_paper_arXiv_submission/Linf_paper.pdf

Act as an expert microbial systems biology reviewer for mSystems.

I want your overall impression of the combined methods paper, with the source papers available for comparison.

Please tell me:
1. whether the combined paper feels suitable for mSystems
2. what the manuscript is trying to contribute
3. whether the central message is clear and convincing
4. what the strongest aspects are
5. what the main weaknesses are
6. what mSystems reviewers would most likely criticize first

Focus especially on:
- whether the paper now works better as a methods/framework paper than the gut-only paper
- whether the CT and gut applications complement each other well
- whether the manuscript is framed at the right level for mSystems

Return your answer in this structure:

1. Overall impression
2. Journal-fit assessment
   - strong fit / plausible fit / weak fit
   - why
3. Strongest aspects
4. Major concerns
   For each:
   - issue
   - where it appears
   - why it matters for mSystems
   - what would improve it
   - confidence: high / medium / low
5. Minor comments
6. Highest-value revisions
```

## Prompt 2: Method Contribution Audit

```text
Read the combined methods paper as an mSystems reviewer focused on the method contribution.

Please assess whether the manuscript makes a clear and compelling case for dCSTs as a methodological contribution.

Evaluate:
1. what is genuinely new here
2. whether the framework is explained clearly enough
3. whether the manuscript explains why dCSTs are useful beyond one application domain
4. whether the paper demonstrates enough methodological generality
5. whether the contribution feels substantial enough for mSystems

Focus especially on:
- determinism
- sample-level representation
- hierarchy and interpretability
- portability across cohorts or applications
- what dCSTs add beyond clustering-defined CST/enterotype-style summaries

Return:
1. What feels genuinely novel
2. What feels incremental or underdeveloped
3. Where the methods explanation is unclear
4. What would make the methods contribution more compelling
```

## Prompt 3: Source-Selection And Synthesis Audit

This prompt is specifically for your question about whether the selected results
from the CT paper, gut paper, and a small amount of the Linf paper are
appropriate.

```text
Compare the combined methods paper against these source documents:

- /Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf
- /Users/pgajer/current_projects/CT_clearance/docs/ct_dcst_manuscript.pdf
- /Users/pgajer/current_projects/Linf_paper_arXiv_submission/Linf_paper.pdf

Your task is to evaluate whether the selection, extraction, and combination of material is appropriate for the new combined methods manuscript.

Please assess:
1. whether the chosen gut results are the right ones for a methods-focused mSystems paper
2. whether the chosen CT results are the right ones for a methods-focused mSystems paper
3. whether the small amount of Linf-paper content is appropriate and well integrated
4. whether any imported results feel unnecessary, too detailed, or distracting
5. whether any important supporting results from the source papers are missing
6. whether the balance between gut and CT is right
7. whether the combined paper preserves the strongest logic of the originals without importing too much application-specific baggage

For each concern, return:
- issue
- where it appears in the combined paper
- where the relevant source material appears
- why the current selection is strong / weak / incomplete / excessive
- what I should add, remove, compress, or reframe

Then give:
1. What was selected well
2. What was selected poorly
3. What important material is missing
4. The 5 most important source-selection improvements
```

## Prompt 4: AI-Generated Notes / Prose Naturalness Audit

This one directly addresses your concern about text that may still read like AI
notes rather than a polished paper.

```text
Read the combined methods paper as a careful scientific editor and mSystems reviewer.

I want you to identify any places that read more like:
- AI-generated notes
- outline residue
- scaffold prose
- internal planning language
- reviewer-response language
- stitched-together transitional text

rather than like a polished human-written research manuscript.

Please look for signs such as:
- repetitive or overly symmetrical sentence structure
- meta-writing about what the manuscript is doing
- language that sounds like planning notes rather than final prose
- transitions that feel mechanical or over-signposted
- vague but polished-sounding generalities
- sections that sound patched together from different source papers

For each issue, return:
- the problematic wording or passage type
- where it appears
- why it sounds artificial, note-like, or non-native to a polished paper
- whether the issue is minor, moderate, or serious
- how I should revise it

Then summarize:
1. whether the manuscript overall reads like a human-written paper
2. which sections are strongest stylistically
3. which sections most need prose cleanup
```

## Prompt 5: Flow And Readability Audit

```text
Read the combined methods paper as an mSystems reviewer and evaluate the flow of ideas.

Please assess:
1. whether the progression of ideas is smooth
2. whether the move from framework to gut application to CT application feels natural
3. whether the manuscript repeats itself unnecessarily
4. whether the discussion circles back cleanly to the main methods claim
5. whether the paper is easy to read for an mSystems audience
6. where the writing could be made smoother, leaner, or better paced

For each flow problem, return:
- issue
- where it appears
- why it disrupts the reading experience
- how I should revise it

Then give:
1. Flow strengths
2. Flow weaknesses
3. The 5 highest-value edits to make the narrative smoother
```

## Prompt 6: Overclaim And Framing Audit

```text
Identify places where the combined methods paper may be framed too strongly or in the wrong way for mSystems.

Focus on:
- method claims
- transportability / portability claims
- reproducibility claims
- disease or clinical interpretation
- generality across domains
- biological significance

Be especially alert to places where:
- the paper sounds too much like a clinical disease paper
- the applications overshadow the dCST framework contribution
- portability is stated more strongly than the evidence warrants
- the manuscript implies mechanistic or causal biology that is not really shown
- the combined paper inherits overclaim from the source papers without reframing it

For each issue, return:
- the claim or wording
- where it appears
- why it is vulnerable from an mSystems perspective
- whether it is a serious evidentiary problem or mainly a framing problem
- a more defensible alternative phrasing
```

## Prompt 7: Final Submission-Readiness Review

Use this after one or more earlier prompts.

```text
Please give a final mSystems-oriented submission-readiness review of the combined dCST methods paper.

Answer these questions directly:
1. Is the paper suitable for submission to mSystems?
2. If yes, what are the main major issues that must still be addressed?
3. What are the most important minor issues?
4. Is the selected material from the gut paper, CT paper, and Linf paper appropriate?
5. Does any part of the paper still read like AI-generated notes rather than a polished manuscript?
6. Is the flow of ideas as smooth as it should be?
7. Is the paper easy to read for the mSystems audience?
8. What else most needs to be addressed before submission?

Return:
1. Go / plausible with revision / not yet suitable
2. Major issues
3. Minor issues
4. Source-selection judgment
5. Prose-naturalness judgment
6. Flow/readability judgment
7. Top 5 next actions
```

## Fast One-Prompt Version

If you want one prompt only, use this:

```text
Please read these documents carefully:

- /Users/pgajer/current_projects/linf/papers/dcst-methods/build/dcst_methods_paper.pdf
- /Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf
- /Users/pgajer/current_projects/CT_clearance/docs/ct_dcst_manuscript.pdf
- /Users/pgajer/current_projects/Linf_paper_arXiv_submission/Linf_paper.pdf

Act as an expert microbial systems biology reviewer for mSystems.

I want you to review the combined methods paper from the standpoint of mSystems and answer these questions:

1. Is the paper suitable for submission to mSystems?
2. If yes, what are the major and minor issues that need to be addressed before submission?
3. Is the selection of results from the CT paper, gut paper, and small amount of Linf-paper content appropriate?
4. Does any part of the paper read like AI-generated notes rather than a polished human-written research paper?
5. Is the flow of ideas smooth, or can it be made smoother?
6. Is the paper easy to read for the mSystems audience?
7. What else should be addressed?

Focus especially on:
- whether the paper now works as a stronger methods/framework paper than the gut-only manuscript
- whether the dCST methodological contribution is clear and substantial
- whether the chosen applications support the framework well
- whether the framing is right for mSystems
- whether the prose and structure feel submission-ready

Important constraints:
- Distinguish clearly between direct observations and inference.
- Do not invent missing analyses or results.
- Cite page, section, figure, table, or supplement item whenever possible.
- When comparing to the source papers, cite the relevant source location too.
- Separate major concerns from minor comments.
- Be specific and constructive.

Return:
1. Overall impression
2. Journal-fit assessment
3. Major concerns
4. Minor concerns
5. Source-selection audit
6. AI-note / prose-naturalness audit
7. Flow/readability audit
8. Final submission-readiness judgment
9. Top 5 revisions
```

## Recommendation

Yes, I think reviewing this paper from the `mSystems` perspective is the right
move.

The combined manuscript should be judged mainly on whether it convincingly
establishes dCSTs as a deterministic, interpretable, reusable representation
framework across applications, not on whether either the gut or CT application
alone would carry the paper as a disease-specific manuscript.
