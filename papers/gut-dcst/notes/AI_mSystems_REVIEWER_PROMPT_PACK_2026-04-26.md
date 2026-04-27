# AI mSystems Reviewer Prompt Pack

This prompt pack is for AI review of the gut dCST paper from the standpoint of
an `mSystems` reviewer rather than a disease-specialist IBD reviewer.

Primary files:

- `/Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper_supplement.pdf`

## Why This Perspective

This paper currently looks strongest as:

- a microbial systems / microbiome methods paper;
- a deterministic representation and transportability paper;
- a cross-cohort validation and interpretability paper;
- a gut-microbiome application of a more general framework.

That makes `mSystems` a more natural review lens than a purely clinical IBD
lens. The core question is not "is this a definitive IBD biology paper?" but
rather:

- does the dCST framework produce meaningful systems-level microbiome insight?
- is the method contribution clear and nontrivial?
- is the validation architecture convincing?
- are the claims about transportability, interpretability, and reproducibility
  appropriate?

## Best Use

Use a staged workflow:

1. mSystems fit and overall impression
2. method and contribution audit
3. validation and reproducibility audit
4. overclaim and framing audit
5. figure/supplement clarity audit
6. highest-value revision priorities

## Core Reviewer Frame

Use this before any of the prompts below if you want tighter output:

```text
Act as an expert microbial systems biology reviewer for mSystems.

Assume you are experienced in:
- microbiome data analysis
- microbial community state representation
- computational and statistical microbiology
- cross-cohort validation
- high-dimensional biological data interpretation

Your job is not to be agreeable. Your job is to determine:
1. whether this manuscript is genuinely interesting for mSystems
2. whether the method contribution is clear and substantial
3. whether the biological application supports the systems-level claims
4. where the paper is strong
5. where it is vulnerable to reviewer criticism
6. what revisions would most improve its chances at mSystems

Important constraints:
- Distinguish clearly between direct observations from the manuscript and your own inference.
- Do not invent missing analyses, results, or literature.
- Cite the relevant page, section, figure, table, or supplement item whenever possible.
- Separate major concerns from minor comments.
- If something is fixable, explain how.
- If something is strong, explain why.
```

## Prompt 1: mSystems Fit And Overall Read

Use this first.

```text
Please read these two documents carefully:

- /Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf
- /Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper_supplement.pdf

Act as an expert microbial systems biology reviewer for mSystems.

Read this as an mSystems reviewer and tell me:
1. what the manuscript is trying to contribute
2. whether the paper feels like a good fit for mSystems
3. what the strongest parts are
4. what the main weaknesses are
5. what an mSystems reviewer would most likely criticize

Focus especially on:
- whether the dCST framework itself feels novel and useful
- whether the gut application demonstrates systems-level value
- whether the manuscript feels more like a framework paper than a disease-association paper
- whether the cross-cohort validation is convincing enough for this journal

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

This is one of the most important prompts for this paper.

```text
Read the manuscript and supplement as an mSystems reviewer focused on the method contribution.

I want you to assess whether the manuscript makes a clear and compelling case for dCSTs as a methodological contribution.

Please evaluate:
1. what is actually new here methodologically
2. whether the manuscript explains the method clearly enough for a systems-microbiology audience
3. whether the contribution is merely descriptive or whether it offers a reusable analytical framework
4. whether the manuscript clearly explains why dCSTs are useful beyond this one IBD application
5. whether the method's strengths and limits are stated honestly

Focus especially on:
- determinism
- representation and hierarchy
- portability across cohorts
- interpretability
- what the framework adds beyond more standard clustering, ordination, or association summaries

Return:
1. What feels genuinely novel
2. What feels incremental or underdeveloped
3. Where the method explanation is unclear
4. What would make the method contribution much more compelling to mSystems reviewers
```

## Prompt 3: Validation And Reproducibility Audit

```text
Evaluate the manuscript's validation and reproducibility logic as an mSystems reviewer.

Please assess whether the manuscript convincingly supports claims about:
- transferability
- external validation
- reproducibility across cohorts
- where portability holds and where it breaks

Focus especially on whether the paper clearly separates:
- discovery-stage screening
- external validation
- rebuilt-cohort rediscovery
- direct transfer of AGP-derived labels
- within-cohort follow-up analyses

Please answer:
1. Is the validation architecture conceptually strong?
2. Is it communicated clearly enough?
3. Are the reproducibility claims appropriately calibrated?
4. Where does the manuscript overstate transfer, portability, or generality?
5. What additional framing or analyses would make the validation story more convincing?

Return:
1. Validation strengths
2. Validation weaknesses
3. Overstated reproducibility/transfer claims
4. The 3 to 5 most valuable validation-related revisions
```

## Prompt 4: Systems-Level Biological Insight Audit

```text
Read the manuscript as an mSystems reviewer who wants both computational novelty and biologically meaningful systems-level insight.

Evaluate:
- whether the dCST framework yields real microbiome insight rather than only a new labeling scheme
- whether the biological interpretations are appropriately linked to the method output
- whether the paper shows enough about ecological structure, hierarchy, and cross-cohort behavior to justify publication in mSystems
- whether the gut application demonstrates general value or feels too narrowly disease-framed

Please identify:
1. where the systems-level insight feels strongest
2. where the paper feels too descriptive
3. where the biology is underinterpreted
4. where the biology is overinterpreted
5. what would improve the balance between computational and biological contribution
```

## Prompt 5: Overclaim And Framing Audit

```text
Identify places where the manuscript may be framed too much like a clinical disease paper or where the claims go beyond what mSystems reviewers would consider well supported.

Focus on:
- disease interpretation
- clinical significance
- transportability
- generalizability
- method claims
- reproducibility claims

Be especially alert to places where:
- the manuscript sounds like a clinical biomarker paper
- the disease application overshadows the framework contribution
- portability is stated more strongly than the evidence warrants
- the paper implies causal or mechanistic disease insight that is not really shown

For each issue, return:
- the claim or wording
- where it appears
- why it is vulnerable from an mSystems reviewer perspective
- whether it is a serious evidentiary problem or mainly a framing problem
- a more defensible alternative phrasing
```

## Prompt 6: Title, Abstract, And Narrative Positioning Audit

```text
Read the title, abstract, introduction, and discussion as an mSystems reviewer.

I want to know whether the paper is positioned correctly for that journal.

Please evaluate:
- whether the title foregrounds the right contribution
- whether the abstract sounds too clinical, too descriptive, or appropriately systems-focused
- whether the introduction clearly motivates the method problem
- whether the discussion lands at the right level of generality
- whether the manuscript clearly says what the paper is and is not claiming

Return:
1. What currently works in the positioning
2. What feels mispositioned for mSystems
3. Specific revisions to make the paper read more like a strong mSystems paper
4. A short proposed rewrite of the central claim in more defensible mSystems-oriented language
```

## Prompt 7: Figure And Supplement Burden Audit

```text
Review the figures, tables, and supplement as an mSystems reviewer.

Please identify:
1. which figures best support the central framework contribution
2. which figures mainly support the disease application
3. whether the balance between method, validation, and biological application is right
4. whether too much essential argument is currently pushed into the supplement
5. whether the reader can understand the main message without doing too much interpretive work

Return:
1. Figure/table strengths
2. Figure/table weaknesses
3. Supplement burden problems
4. Specific presentation changes that would improve clarity and journal fit
```

## Prompt 8: Most Likely mSystems Reviewer Objections

```text
If this manuscript were reviewed at mSystems, what are the 5 most likely reviewer objections?

Please make these realistic reviewer objections, not generic ones.

For each:
- write the objection as a reviewer might actually phrase it
- explain whether it is fair
- say whether it threatens journal fit, methodological credibility, validation credibility, or only framing
- explain how I should revise the manuscript to blunt or resolve it
```

## Prompt 9: Highest-Value Revision Priorities

```text
Assume this manuscript is fundamentally worth improving for mSystems, but the authors can only make a limited number of changes before submission.

What are the 5 revisions that would most improve:
- journal fit
- methodological clarity
- reviewer confidence
- systems-level biological value
- overall publishability at mSystems

For each revision:
- say exactly what should change
- explain why it matters
- classify it as mainly:
  - framing
  - methods explanation
  - validation logic
  - figure/supplement presentation
  - biological interpretation

Then rank them from most important to least important.
```

## Prompt 10: Fast One-Prompt Version

If you only want one prompt, use this:

```text
Please read these two documents carefully:

- /Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf
- /Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper_supplement.pdf

Act as an expert microbial systems biology reviewer for mSystems.

Your job is not to be agreeable. Read this as an mSystems reviewer and tell me:
- whether this feels like a good fit for mSystems
- what is genuinely strong
- what is weak or vulnerable
- whether the dCST methodology is presented clearly and compellingly
- whether the validation and portability claims are appropriately calibrated
- whether the gut disease application supports the framework without overclaiming clinical significance
- what 5 changes would most improve the paper before submission

Important constraints:
- Distinguish clearly between direct observations and inference.
- Do not invent missing analyses or results.
- Cite page, section, figure, table, or supplement item whenever possible.
- Separate major concerns from minor comments.
- Be specific and constructive.

Return:
1. Overall impression
2. Journal-fit assessment
3. Strongest aspects
4. Major concerns
5. Minor comments
6. Overclaim/framing audit
7. Highest-value revisions
```

## Recommendation

Yes, I think getting review from the `mSystems` perspective is the right move.

The strongest version of this manuscript is not "we have a definitive IBD
paper." It is closer to:

- "we have a deterministic community-state framework,"
- "we show how it behaves across gut cohorts,"
- "we test where transport works and where it breaks,"
- "and we use IBD/IBS and related outcomes as a strong application space."

That is exactly the kind of framing where an `mSystems` review lens is more
useful than a purely disease-specialist lens.
