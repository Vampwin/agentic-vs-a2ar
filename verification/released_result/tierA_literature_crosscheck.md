# Tier-A literature cross-check

**Searched 2026-08-21 via PubMed.** Both hits and non-hits are recorded. This
search was performed *after* the candidates were identified computationally and
is therefore contextual, not confirmatory. It does not validate any compound.

An earlier version of this project asserted a flavone/A2AR precedent from
memory, without a search, and that claim was formally withdrawn (ADR-AIIM-073)
because its A2A data turned out to be rat-receptor. This document exists so
that the same mistake cannot recur: every statement below is tied to a
retrieved record, and the rat/human distinction is stated wherever it applies.

## The three candidates

| ID | Formula | SMILES | Identification |
|---|---|---|---|
| NPLIB01185 | C19H12O6 | `COc1c2ccoc2cc2oc(-c3ccc4c(c3)OCO4)cc(=O)c12` | methoxy-furanoflavone with a methylenedioxyphenyl B-ring |
| NPLIB01643 | C18H12O4 | `COc1ccc(-c2cc(=O)c3cc4ccoc4cc3o2)cc1` | methoxyphenyl furanoflavone (linear furan fusion) |
| NPLIB03056 | C18H12O4 | `COc1c2ccoc2cc2oc(-c3ccccc3)cc(=O)c12` | 3-methoxy-2-phenyl furanoflavone — **consistent with karanjin** |

All three are furanoflavones recorded in LOTUS as *Pongamia pinnata*
constituents.

## Hit 1 — flavones do bind adenosine receptors, in the micromolar range

**Ji XD, Melman N, Jacobson KA. *Interactions of flavonoids and other
phytochemicals with adenosine receptors.* J Med Chem 1996;39(3):781-8.**
[10.1021/jm950661k](https://doi.org/10.1021/jm950661k)

Retrieved facts:

- Flavone derivatives bind three adenosine receptor subtypes in the
  **micromolar** range.
- **Galangin: Ki = 1 µM at rat A1 and rat A2A**; 3 µM at human A3.
- A2A affinities in that study were measured at **rat brain A2A receptors**
  using [3H]CGS21680 — **not human A2A**.
- Structure–activity finding, directly relevant here: *"the hydroxyl groups of
  naturally occurring flavones are not essential for affinity at adenosine
  receptors."* The three candidates are **methoxylated, not hydroxylated**, so
  this observation is the reason their scaffold class is not excluded a priori.
- The authors' own conclusion: flavones "may provide leads for the development
  of novel adenosine antagonists."

**How this may and may not be used.** It establishes that the flavone scaffold
class has measurable, micromolar adenosine-receptor affinity, and that
methylation of the hydroxyls does not abolish it. It does **not** establish
anything about these three compounds, and its A2A numbers are rat-receptor
while this study's ligand-based model is trained on human-receptor ChEMBL
data. The species mismatch is precisely what invalidated the earlier withdrawn
claim and must be stated in the manuscript wherever this reference is cited.

## Hit 2 — one candidate matches a well-characterised natural product

**Singh A, Bhatt G, Gujre N, *et al.* *Karanjin.* Phytochemistry
2021;183:112641.**
[10.1016/j.phytochem.2020.112641](https://doi.org/10.1016/j.phytochem.2020.112641)

Retrieved facts:

- Karanjin is a **furanoflavonoid first isolated from *Pongamia pinnata***,
  IUPAC 3-methoxy-2-phenylfuro-(2,3-h-chrome-4-ol), molecular formula C18H12O4.
- Reported activities: anti-diabetic, anti-cancer, anti-inflammatory,
  anti-hyperglycaemic, antioxidant, anti-colitis, anti-ulcer and
  **anti-Alzheimer**.
- Toxicological studies cited as showing karanjin **non-toxic at physiological
  conditions** in vitro and in vivo.

NPLIB03056's formula (C18H12O4) and connectivity are consistent with karanjin.
**This is a structural correspondence noted from formula and SMILES, not an
authenticated identification**; confirming it would require comparison against
an authentic standard.

Relevance: if the correspondence holds, one of the three candidates is a
known, isolable, non-toxic natural product with reported CNS-relevant activity
— which matters for feasibility of any follow-up, and not at all for whether
it binds A2AR.

## Non-hits — recorded because absence of evidence is the finding

| Query | Result |
|---|---|
| `Pongamia pinnata furanoflavone karanjin adenosine receptor` | **0 records** |
| `karanjin pharmacology neuroprotective anti-inflammatory Pongamia` | 0 records |
| `furanoflavone adenosine receptor antagonist` | **0 records** |
| `pongapin lanceolatin Pongamia pinnata furanoflavonoid biological activity` | 0 records |

**No published record links any *Pongamia pinnata* furanoflavone, or any
furanoflavone at all, to an adenosine receptor.** The three candidates are
therefore not rediscoveries of a known adenosine ligand. That cuts both ways:
there is no prior support for them, and equally no prior report excluding them.

## What the manuscript may state

Permissible:

- The flavone scaffold class has documented micromolar affinity at adenosine
  receptors, with A2A values measured in rat, and hydroxyl groups are not
  required for that affinity (Ji et al. 1996).
- One candidate is structurally consistent with karanjin, a characterised
  *Pongamia pinnata* furanoflavonoid reported as non-toxic (Singh et al. 2021).
- No furanoflavone has been reported at any adenosine receptor.

Not permissible:

- Any statement that these compounds are known or expected A2AR ligands.
- Any use of the galangin rat-A2A Ki as support for a human-A2A prediction.
- Any claim that the karanjin correspondence is a confirmed identification.
