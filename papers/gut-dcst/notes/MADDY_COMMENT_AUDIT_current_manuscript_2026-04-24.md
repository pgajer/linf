# Audit of Maddy's Word Comments Against the Current Manuscript

Date: 2026-04-24

Reviewed files:

- commented draft: `/Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper_Built_2026-04-20_081555_EDT_maa.docx`
- current manuscript: `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`

Scope:

- This audit covers the explicit Word comment balloons authored by Madeline Alizadeh in the commented DOCX.
- I did not treat every tracked insertion/deletion as a separate comment item unless it also had a comment balloon.

High-level summary:

- `14` explicit comments were audited.
- `14` are addressed in the current manuscript after the final wording pass completed on 2026-04-24.

Bottom line:

- The major structural critiques were addressed.
  - The manuscript no longer relies on AGP as the main disease-facing IBD evidence.
  - The main Results are now centered on external clinically annotated IBD cohorts.
  - IBS and autoimmune analyses were removed from the main Results narrative and pushed to the supplement.
- The remaining wording issues identified in the first audit pass were also addressed.
  - generic `phenotype` wording was replaced with `condition`, `outcome`, `case definition`, or literal reference to the AGP `Phenotype` field only where needed;
  - the unclear sentence about low-resolution retained-feature labels was rewritten;
  - residual shorthand that could sound more clinical than the data support was tightened.

## Comment-by-comment audit

### Comment 1

- Anchor in commented DOCX: `Self-reported`
- Comment text:
  `how is this defined? were they confirmed diagnosis previously and then self-reported? this might be a huge point of contention if you get an IBD expert reviewing as diagnostic methodology is an important component of methods in IBD based studies`
- Status: `Addressed`
- How it was addressed:
  - The manuscript now explicitly downgrades AGP from disease-defining evidence to structure-learning context.
  - The Abstract now states that AGP is used `primarily as a large ecological structure-learning cohort rather than as a source of diagnosis-grade IBD case definitions`.
  - The Methods now define the AGP outcome fields as fixed substring matches from the AGP `Phenotype` metadata field and explicitly state that they are not diagnosis-grade clinical case definitions.
- Current manuscript evidence:
  - Abstract lines around 180-191.
  - Introduction lines around 287-295.
  - Methods lines around 500-519.

### Comment 2

- Anchor in commented DOCX:
  `IBS and autoimmune disease showed weaker but directionally related dominance patterns.`
- Comment text:
  `from a clinical sense, it doesn’t make sense to group these together. i think the focus should either be a 1:1 or IBD:IBS, or an analysis that incorporates the other autoimmune data (and this can be much more flexible to the degree to which you include it or focus on it).`
- Status: `Addressed`
- How it was addressed:
  - The current Abstract no longer groups IBS and autoimmune disease as a shared secondary biological result.
  - Instead, both are moved out of the main narrative and explicitly described as exploratory supplement material.
  - The main text is now IBD-centered, with external CD/UC datasets carrying the disease-facing story.
- Current manuscript evidence:
  - Abstract lines around 200-204.
  - Introduction lines around 298-305.
  - Results lines around 437-440.

### Comment 4

- Anchor in commented DOCX: `Introduction`
- Comment text:
  `we have to be careful with this. many people are in a very deep remission, but inflammatory biomarkers don’t always equate to clinical observance of inflammation, so it’s hard to make any assumptions if we aren’t directly provided this information by those running a cohort.`
- Status: `Addressed`
- How it was addressed:
  - The manuscript now explicitly states that AGP does not provide diagnosis-grade disease activity or treatment information.
  - The disease-facing claims are now anchored in external clinically annotated cohorts rather than in AGP alone.
  - The Discussion also explicitly notes the lack of disease-activity information and warns against overinterpretation.
  - Residual shorthand that could imply measured inflammation was tightened in the final wording pass.
- Current manuscript evidence:
  - Introduction lines around 289-295.
  - Methods lines around 515-519.
  - Discussion lines around 649-653.

### Comment 40

- Anchor in commented DOCX:
  `but they are not used as primary biological evidence for the IBD interpretation.`
- Comment text:
  `im not quite sure what this sentence means. is this meant in a clinical diagnostic context?`
- Status: `Addressed`
- How it was addressed:
  - The original vague sentence was rewritten in plain language.
  - The manuscript now states directly that these low-resolution or potentially non-gut labels are kept for accounting but are not interpreted as disease-relevant gut states in the IBD results.
- Current manuscript evidence:
  - Results lines around 347-354.

### Comment 41

- Anchor in commented DOCX: `phenotype`
- Comment text:
  `phenotype has a very specific meaning in IBD, usually meaning you’re comparing inflammatory, stricturing, penetrating, and stricturing+penetrating disease. I don’t think that’s what is meant here`
- Status: `Addressed`
- How it was addressed:
  - The generic uses of `phenotype` were replaced with `condition`, `outcome`, `case definition`, or more literal phrasing.
  - The only remaining `Phenotype` instance is the literal name of the AGP metadata field.
- Current manuscript evidence:
  - Abstract lines around 180-191.
  - Figure 2 caption lines around 317-320.
  - Results lines around 356-357.
  - Methods lines around 500-514 and 549-553.
  - Discussion lines around 646-649.

### Comment 42

- Anchor in commented DOCX: `phenotype`
- Comment text:
  `do you mean depthwise condition screen?`
- Status: `Addressed`
- How it was addressed:
  - The Results text now uses `exploratory depthwise AGP condition screen`.
- Current manuscript evidence:
  - Results lines around 356-357.

### Comment 45

- Anchor in commented DOCX:
  `non-IBD examples are summarized in the supplement`
- Comment text:
  `IBS could be included here!`
- Status: `Addressed`
- How it was addressed:
  - IBS is now explicitly named as exploratory supplement material rather than left ambiguous.
- Current manuscript evidence:
  - Results lines around 370-372.
  - Results lines around 437-440.

### Comment 47

- Anchor in commented DOCX: `Self-reported IBD`
- Comment text:
  `again, need to clarify.`
- Status: `Addressed`
- How it was addressed:
  - The original AGP-self-reported-IBD subsection was removed.
  - The manuscript now repeatedly clarifies that AGP IBD is self-reported, not subtype-resolved, and not the main source of admissible disease inference.
- Current manuscript evidence:
  - Introduction lines around 289-295.
  - Methods lines around 512-519.
  - Discussion lines around 611-619 and 642-645.

### Comment 48

- Anchor in commented DOCX:
  `complete posterior support for enrichment in IBD`
- Comment text:
  `this should probably be given a bit more explanation as to what it means`
- Status: `Addressed`
- How it was addressed:
  - That exact language no longer appears in the main manuscript.
  - The Methods now define the Bayesian summaries in more explicit statistical terms:
    posterior median odds ratio, credible interval, and posterior probability that the odds ratio exceeded 1.
- Current manuscript evidence:
  - Methods lines around 560-569.

### Comment 49

- Anchor in commented DOCX:
  `Morganella contained 22 IBD samples`
- Comment text:
  `shouldnt it be the other day around, that 22 IBD samples contained Morganella?`
- Status: `Addressed`
- How it was addressed:
  - The count sentence is no longer in the main manuscript.
  - The rare-state example was simplified and de-emphasized.
- Current manuscript evidence:
  - The old count sentence has been removed from the current main text.
  - Rare enterobacterial states are now mentioned only in a more compressed exploratory way in the Abstract and Discussion.

### Comment 50

- Anchor in commented DOCX:
  `or as evidence of disease severity`
- Comment text:
  `you don’t have data to even make this connection.you would need endoscopic and disease activity indices which are standard (SESCD/Rutgeerts, or Mayo are the endoscopic scores to look for, and CDAI, full Mayo are the total DAI)`
- Status: `Addressed`
- How it was addressed:
  - The disease-severity suggestion is no longer present.
  - The current manuscript explicitly says AGP does not support disease-activity or detailed treatment inference.
- Current manuscript evidence:
  - Introduction lines around 289-291.
  - Methods lines around 515-519.
  - Discussion lines around 649-653.

### Comment 51

- Anchor in commented DOCX:
  `a Bacteroides-first retained-feature context`
- Comment text:
  `this should be explained more clearer for readers`
- Status: `Addressed`
- How it was addressed:
  - The specific dense explanatory paragraph was removed from the main Results.
  - The concept is now introduced more cleanly upfront through the early dCST definition and the dominance-lineage explanation in the Introduction.
- Current manuscript evidence:
  - Introduction lines around 255-271.
  - The old `Bacteroides-first retained-feature context` sentence no longer appears in the main Results.

### Comment 52

- Anchor in commented DOCX: `AGP phenotype`
- Comment text:
  `again, would call these conditions, as phenotype has a specific meaning in crohns`
- Status: `Addressed`
- How it was addressed:
  - The AGP methods language now consistently uses `condition`, `outcome`, `condition-by-state`, and `condition-specific`.
  - The only remaining `Phenotype` token is the literal AGP field name.
- Current manuscript evidence:
  - Methods lines around 500-514.
  - Methods lines around 551-553 and 579-581.
  - Discussion lines around 646-649.

### Comment 53

- Anchor in commented DOCX:
  `, acid reflux, or seasonal allergies`
- Comment text:
  `there isn’t a biological relevance to compare these conditions based on what’s being looked at, and i would remove from the manuscript based on simplicity sake, and make IBS a supplement`
- Status: `Addressed`
- How it was addressed:
  - IBS, autoimmune disease, acid reflux, and seasonal allergies were removed from the main Results narrative.
  - They remain only as exploratory AGP screens summarized in the supplement and briefly acknowledged in Methods or framing sentences.
- Current manuscript evidence:
  - Introduction lines around 298-305.
  - Results lines around 370-372 and 437-440.

## Overall conclusion

Maddy's major critiques were successfully incorporated into the new manuscript architecture.

The current version is much closer to what her review was pushing toward:

- AGP for hierarchy learning rather than diagnosis-grade inference;
- external clinically annotated IBD cohorts as the main disease-facing evidence;
- CD and UC separated where metadata allow it;
- IBS and other non-IBD analyses moved out of the main story.

After the final wording pass, I consider the explicit Word comments closed in the current manuscript. Any remaining future changes would be optional manuscript-development refinements rather than unresolved comment debt from this review round.
