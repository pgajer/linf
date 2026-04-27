# AI IBD Expert Prompt Pack

This prompt pack is for expert-style AI reading of the gut DCST IBD paper and
its supplement:

- `/Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper_supplement.pdf`

The goal is not generic praise. The goal is to get domain-serious feedback from
an IBD-clinician-scientist perspective.

## Best Use

Do not ask for only "impressions." That usually produces polite but shallow
feedback. Use a staged workflow:

1. overall IBD expert read
2. skeptical clinical-audience audit
3. validation and overclaim audit
4. figure/supplement clarity audit
5. highest-value revision priorities

## Core Reviewer Frame

Use this frame before any of the prompts below if you want tighter output:

```text
Act as a PhD medical doctor, gastroenterologist, and IBD expert who is also an experienced scientific reviewer.

Your job is not to be agreeable. Your job is to identify:
1. what is genuinely strong
2. what is weak, unclear, or vulnerable to criticism
3. where the clinical or biological interpretation is too strong
4. what would most improve the manuscript for an expert IBD audience

Important constraints:
- Distinguish clearly between direct observations from the manuscript and your own inference.
- Do not invent missing analyses, results, or literature.
- Cite the relevant page, section, figure, table, or supplement item whenever possible.
- Separate major concerns from minor comments.
- If something is fixable, explain how.
- If something is persuasive, explain why.
```

## Prompt 1: Overall IBD Expert Read

Use this first.

```text
Please read these two documents carefully:

- /Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf
- /Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper_supplement.pdf

Act as a PhD medical doctor, gastroenterologist, and IBD expert who is also an experienced scientific reviewer.

Read this as an IBD specialist and tell me:
1. what the manuscript is trying to show
2. whether the central message would feel interesting and credible to an IBD audience
3. what the strongest parts are
4. what the main weaknesses are
5. what IBD reviewers would be most likely to criticize

Focus especially on:
- clinical relevance to IBD
- whether the framing matches what clinicians and IBD researchers would find convincing
- whether the manuscript explains clearly why these dominance states matter
- whether the external validation strategy feels persuasive

Return your answer in this structure:

1. Overall impression
2. Strongest aspects
3. Major concerns
   For each:
   - issue
   - where it appears
   - why it matters for an IBD audience
   - what would improve it
   - confidence: high / medium / low
4. Minor comments
5. Highest-value revisions
```

## Prompt 2: Skeptical IBD Reviewer Audit

Use this after Prompt 1.

```text
Read the manuscript and supplement again as a skeptical IBD reviewer.

Assume you are knowledgeable about:
- Crohn disease and ulcerative colitis heterogeneity
- stool vs mucosal microbiome interpretation
- treatment effects and disease activity confounding
- cohort heterogeneity across adult, pediatric, and clinically annotated IBD studies

Your task is to identify the places where a strong IBD reviewer would say:
- "this is interesting, but not yet clinically convincing"
- "this interpretation is too broad for the evidence shown"
- "this needs better caveating or clearer framing"

Focus especially on:
- the use of AGP self-reported IBD as the discovery-stage screen
- whether the jump from AGP screening to clinically annotated external IBD cohorts feels conceptually solid
- whether the manuscript sufficiently distinguishes pooled IBD, Crohn disease, and ulcerative colitis
- whether cohort differences could explain some of the transfer pattern
- whether the claims are stronger than the evidence on portability, reproducibility, or biological meaning

Return:
1. The 5 most likely reviewer objections
2. For each objection:
   - where it appears
   - whether the issue is major or mainly framing
   - how I would revise the manuscript to reduce the objection
3. Which objections are probably overstated and should not distract me
```

## Prompt 3: Validation Architecture Audit

This one is especially important for this paper.

```text
Evaluate the manuscript's validation logic as an IBD expert.

I want you to assess whether the manuscript clearly and convincingly separates:
- exploratory AGP screening
- external validation in clinician-annotated IBD cohorts
- rebuilt-cohort rediscovery of disease-linked structure
- direct AGP-derived dCST label portability
- within-disease follow-up analyses such as Crohn location and UC extent

Please answer:
1. Is the validation architecture conceptually sound?
2. Is it explained clearly enough for an expert but non-methods-focused IBD reader?
3. Where does the logic feel strongest?
4. Where does it feel confusing, under-motivated, or vulnerable to criticism?
5. Does the manuscript overstate exact label transfer when the stronger result may actually be rediscovery of related structure?

Return:
1. Validation strengths
2. Validation weaknesses
3. Places where the logic should be clarified or reframed
4. The 3 revisions that would most improve credibility of the validation story
```

## Prompt 4: Clinical Relevance And Biological Plausibility Audit

```text
Read the manuscript as an IBD clinician-scientist and evaluate the biological and clinical plausibility of the findings.

Focus especially on:
- whether the manuscript connects the reported dCST patterns to known IBD biology in a convincing way
- whether taxa-level interpretation is careful enough
- whether the manuscript distinguishes hypothesis generation from clinically actionable insight
- whether the Crohn disease vs ulcerative colitis interpretation is appropriately cautious
- whether the discussion of stool-based states, mucosal validation, and subtype structure is balanced

Please identify:
1. places where the biological interpretation feels strong and well grounded
2. places where it feels speculative or overstated
3. missing caveats or missing IBD context
4. missing literature or conceptual framing that an IBD reviewer might expect

Return:
1. Well-grounded biological/clinical interpretations
2. Weak or overstated interpretations
3. Missing context or caveats
4. Specific wording changes that would make the discussion more credible
```

## Prompt 5: Overclaim And Wording Audit

```text
Identify every place where the manuscript's claims may go beyond the evidence shown.

Focus on:
- clinical significance
- biological interpretation
- portability / generalizability
- reproducibility
- disease-subtype claims
- translational relevance

For each issue, return:
- the claim or wording
- where it appears
- why it may be too strong
- whether this is a serious evidentiary problem or mainly a framing problem
- a more defensible alternative phrasing

Be especially alert to places where:
- self-reported AGP findings are implicitly treated as stronger than they are
- transfer and rediscovery are blurred together
- pooled IBD results are allowed to stand in for subtype-specific insight
- within-disease follow-up findings are presented more confidently than warranted
```

## Prompt 6: Figure And Supplement Clarity Audit

```text
Review the manuscript figures, tables, and supplement as an expert IBD reader.

Please identify:
1. which figures best support the paper's central message
2. which figures are hard to interpret without too much effort
3. whether the supplement is carrying too much essential argument
4. places where figure legends, table framing, or supplement organization could be improved
5. whether an IBD audience would clearly understand the main takeaway from Figures 2 and 3 and the corresponding supplement material

Return:
1. Figure/table strengths
2. Figure/table clarity problems
3. Supplement problems
4. Specific presentation changes that would improve expert readability
```

## Prompt 7: Title, Abstract, And Framing Audit

```text
Read the title, abstract, introduction, and discussion as an IBD expert.

Evaluate:
- whether the title is appropriately ambitious
- whether the abstract accurately reflects what the manuscript truly shows
- whether the introduction motivates the problem in a way that feels important to IBD readers
- whether the discussion lands at the right level of biological and clinical interpretation

Return:
1. What currently works in the framing
2. What feels too broad, too technical, or not persuasive enough
3. Specific edits that would make the framing stronger for an IBD audience
4. A short proposed rewrite of the central claim in more defensible language
```

## Prompt 8: Highest-Value Revision Priorities

Use this near the end after one or more earlier prompts.

```text
Assume the manuscript is fundamentally worth improving, but the authors can only make a limited number of revisions before submission.

Based on the paper and supplement, what are the 5 changes that would most improve:
- credibility
- clarity
- persuasiveness to an IBD audience
- reviewer reception

For each revision:
- say exactly what should change
- explain why it matters
- say whether it is mainly:
  - framing
  - analysis presentation
  - figure/supplement presentation
  - clinical context
  - discussion/interpretation

Then rank them from most important to least important.
```

## Prompt 9: What Will Reviewers Attack First?

```text
If this paper were sent to a good but skeptical microbiome/IBD journal reviewer, what are the first 5 criticisms they would be most likely to raise?

Please make these realistic reviewer criticisms, not generic ones.

For each:
- write the criticism as a reviewer might actually phrase it
- explain whether the criticism is fair
- say how I should revise the paper to blunt or resolve it
```

## Prompt 10: Fast Version

If you only want one prompt, use this:

```text
Please read these two documents carefully:

- /Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf
- /Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper_supplement.pdf

Act as a PhD medical doctor, gastroenterologist, and IBD expert who is also an experienced scientific reviewer.

Your job is not to be agreeable. Read this as an expert IBD reviewer and tell me:
- what is genuinely strong
- what is weak or vulnerable
- where the interpretation goes beyond the evidence
- what an IBD audience may find unconvincing or underexplained
- what 5 changes would most improve the manuscript before submission

Focus especially on:
- the logic of AGP discovery versus external clinician-annotated IBD validation
- whether transfer, rediscovery, and portability are clearly distinguished
- whether pooled IBD, Crohn disease, and ulcerative colitis are interpreted appropriately
- whether the biological and clinical framing feels credible

Important constraints:
- Distinguish clearly between direct observations and inference.
- Do not invent missing analyses or results.
- Cite section, figure, table, supplement item, or page whenever possible.
- Separate major concerns from minor comments.
- Be specific and constructive.

Return:
1. Overall impression
2. Strongest aspects
3. Major concerns
4. Minor comments
5. Overclaim audit
6. Highest-value revisions
```

## Should This Become A Skill?

Probably not yet.

For this paper, a note-local prompt pack is the right size. Turn this into a
Codex skill only if you find yourself repeatedly doing the same workflow across
multiple biomedical papers, for example:

- paper review as a clinician expert
- overclaim audit
- clinical relevance audit
- figure/supplement clarity audit
- revision-priority ranking

If that becomes a recurring pattern across projects, then a reusable skill such
as `biomedical-manuscript-review` or `clinical-paper-audit` would make sense.
For one paper or one manuscript family, a prompt pack is cleaner and easier to
evolve.
