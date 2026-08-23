"""
build_references_v6.py -- attach a numbered reference list to the v6 manuscript.

Two jobs, kept separate so each can be audited:

  1. REFS   the bibliography. Every entry records where its metadata came from.
            `verified: pubmed` means the record was retrieved from PubMed in
            this project and the volume/pages/DOI are the retrieved ones.
            `verified: v4` means it was carried over from the v4 submission's
            own reference list unchanged. `verified: author_action` means the
            claim it supports is grounded in the project record but the
            bibliographic record was never captured -- those are emitted with
            an explicit [AUTHOR ACTION] note rather than a guessed citation.

  2. CITES  the insertions. Each is an exact (anchor, replacement) pair applied
            to the manuscript text. A missing anchor is a hard error: silently
            skipping one would leave a claim uncited while the script reported
            success.

Numbering is by order of first appearance, computed after insertion, which is
what Elsevier's numbered style requires.

This is a ONE-SHOT transformation, not an idempotent one: once it has run, the
anchors it looks for have been replaced by the citations it inserted, and a
second run reports every anchor as missing. That is the intended behaviour --
the alternative, silently matching already-cited text, would double-cite. The
pre-citation manuscript is kept at
`30_manuscript/CURRENT_v5/.manuscript_AIIM_v6.md.pre-refs.bak`; to change the
reference list or a citation placement, edit this script, restore that file,
and run again. Note that the backup predates any manuscript edit made after
the citations were inserted, so check `git diff`-style before restoring.

References deliberately NOT carried over from v4, with reasons:
  DiffDock, ESM         no generative pose prediction is reported in v6
  OpenMM, ff14SB, TIP3P,
  MMPBSA.py, interaction
  entropy               superseded or unused: v6's membrane protocol is
                        ff19SB/Lipid21/OPC under GROMACS, and reports no
                        end-state free energies. GAFF and AM1-BCC ARE carried
                        over -- they parameterise the ligand in that protocol.
  PLIP                  checked: the interaction fingerprint in
                        `ifp_rescore_benchmark.py` is implemented in-house and
                        is purely geometric. PLIP was never called, so citing
                        it would misattribute the method.
  Antonioli 2013        adenosine in immunity/cancer; v6's framing is
                        Parkinson's disease only.

Usage:
    python build_references_v6.py [--check]      # --check = report, do not write
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
MS = ROOT / "30_manuscript/CURRENT_v5/manuscript_AIIM_v6.md"

# --------------------------------------------------------------------------
# The bibliography. Keys are the markers used in CITES below.
# --------------------------------------------------------------------------
REFS: dict[str, dict[str, str]] = {
    "schwarzschild2006": {"verified": "v4", "text":
        "Schwarzschild MA, Agnati L, Fuxe K, Chen J-F, Morelli M. Targeting adenosine "
        "A2A receptors in Parkinson's disease. Trends Neurosci. 2006;29(11):647–654. "
        "doi:10.1016/j.tins.2006.09.004"},
    "ohno2023": {"verified": "pubmed", "text":
        "Ohno Y, Suzuki M, Asada H, et al. In vitro pharmacological profile of KW-6356, "
        "a novel adenosine A2A receptor antagonist/inverse agonist. Mol Pharmacol. "
        "2023;103(6):311–324. doi:10.1124/molpharm.122.000633"},
    "kondo2015": {"verified": "v4", "text":
        "Kondo T, Mizuno Y; Japanese Istradefylline Study Group. A long-term study of "
        "istradefylline safety and efficacy in patients with Parkinson disease. Clin "
        "Neuropharmacol. 2015;38(2):41–46. doi:10.1097/WNF.0000000000000073"},
    "liu2012": {"verified": "v4", "text":
        "Liu W, Chun E, Thompson AA, et al. Structural basis for allosteric regulation "
        "of GPCRs by sodium ions. Science. 2012;337(6091):232–236. "
        "doi:10.1126/science.1219218"},
    "kapoor2023": {"verified": "v4", "text":
        "Kapoor S, Narayanan A. Leakage and the reproducibility crisis in "
        "machine-learning-based science. Patterns. 2023;4(9):100804. "
        "doi:10.1016/j.patter.2023.100804"},
    "roberts2021": {"verified": "v4", "text":
        "Roberts M, Driggs D, Thorpe M, et al. Common pitfalls and recommendations for "
        "using machine learning to detect and prognosticate for COVID-19 using chest "
        "radiographs and CT scans. Nat Mach Intell. 2021;3(3):199–217. "
        "doi:10.1038/s42256-021-00307-0"},
    "degrave2021": {"verified": "v4", "text":
        "DeGrave AJ, Janizek JD, Lee SI. AI for radiographic COVID-19 detection selects "
        "shortcuts over signal. Nat Mach Intell. 2021;3(7):610–619. "
        "doi:10.1038/s42256-021-00338-7"},
    "wynants2020": {"verified": "v4", "text":
        "Wynants L, Van Calster B, Collins GS, et al. Prediction models for diagnosis "
        "and prognosis of covid-19: systematic review and critical appraisal. BMJ. "
        "2020;369:m1328. doi:10.1136/bmj.m1328"},
    "claudecode": {"verified": "v4", "text":
        "Anthropic. Claude Code: agentic coding in the terminal. 2025. "
        "https://claude.ai/code [accessed 30 June 2026]."},
    "react2023": {"verified": "v4", "text":
        "Yao S, Zhao J, Yu D, Du N, Shafran I, Narasimhan K, Cao Y. ReAct: synergizing "
        "reasoning and acting in language models. In: Proceedings of the 11th "
        "International Conference on Learning Representations (ICLR); 2023."},
    "jia2011": {"verified": "web", "text":
        "Jia Y, Harman M. An analysis and survey of the development of mutation "
        "testing. IEEE Trans Softw Eng. 2011;37(5):649–678. doi:10.1109/TSE.2010.62"},
    "chen2018": {"verified": "web", "text":
        "Chen TY, Kuo F-C, Liu H, Poon P-L, Towey D, Tse TH, Zhou ZQ. Metamorphic "
        "testing: a review of challenges and opportunities. ACM Comput Surv. "
        "2018;51(1):Article 4. doi:10.1145/3143561"},
    "snplib": {"verified": "v4", "text":
        "Isolation source papers for the SNP-Lib compound collection, compiled from the "
        "Sritularak group's publication record (Faculty of Pharmaceutical Sciences, "
        "Chulalongkorn University), 2009–2025. The 38 primary reports are listed in "
        "full in Supplementary Reference S1."},
    "mysinger2012": {"verified": "v4", "text":
        "Mysinger MM, Carchia M, Irwin JJ, Shoichet BK. Directory of Useful Decoys, "
        "Enhanced (DUD-E): better ligands and decoys for better benchmarking. J Med "
        "Chem. 2012;55(14):6582–6594. doi:10.1021/jm300687e"},
    "fourches2010": {"verified": "pubmed", "text":
        "Fourches D, Muratov E, Tropsha A. Trust, but verify: on the importance of "
        "chemical structure curation in cheminformatics and QSAR modeling research. "
        "J Chem Inf Model. 2010;50(7):1189–1204. doi:10.1021/ci100176x"},
    "fourches2016": {"verified": "pubmed", "text":
        "Fourches D, Muratov E, Tropsha A. Trust, but verify II: a practical guide to "
        "chemogenomics data curation. J Chem Inf Model. 2016;56(7):1243–1252. "
        "doi:10.1021/acs.jcim.6b00129"},
    "rdkit": {"verified": "v4", "text":
        "Landrum G, et al. RDKit: open-source cheminformatics. Version 2025.09.5. "
        "Zenodo. doi:10.5281/zenodo.5548890"},
    "baell2010": {"verified": "v4", "text":
        "Baell JB, Holloway GA. New substructure filters for removal of pan assay "
        "interference compounds (PAINS) from screening libraries and for their "
        "exclusion in bioassays. J Med Chem. 2010;53(7):2719–2740. "
        "doi:10.1021/jm901137j"},
    "lipinski1997": {"verified": "v4", "text":
        "Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and computational "
        "approaches to estimate solubility and permeability in drug discovery and "
        "development settings. Adv Drug Deliv Rev. 1997;23(1–3):3–25. "
        "doi:10.1016/S0169-409X(96)00423-1"},
    "veber2002": {"verified": "v4", "text":
        "Veber DF, Johnson SR, Cheng H-Y, Smith BR, Ward KW, Kopple KD. Molecular "
        "properties that influence the oral bioavailability of drug candidates. J Med "
        "Chem. 2002;45(12):2615–2623. doi:10.1021/jm020017n"},
    "brenk2008": {"verified": "v4", "text":
        "Brenk R, Schipani A, James D, Krasowski A, Gilbert IH, Frearson J, Wyatt PG. "
        "Lessons learnt from assembling screening libraries for drug discovery for "
        "neglected diseases. ChemMedChem. 2008;3(3):435–444. "
        "doi:10.1002/cmdc.200700139"},
    "bickerton2012": {"verified": "v4", "text":
        "Bickerton GR, Paolini GV, Besnard J, Muresan S, Hopkins AL. Quantifying the "
        "chemical beauty of drugs. Nat Chem. 2012;4(2):90–98. doi:10.1038/nchem.1243"},
    "rutz2022": {"verified": "pubmed", "text":
        "Rutz A, Sorokina M, Galgonek J, et al. The LOTUS initiative for open knowledge "
        "management in natural products research. eLife. 2022;11:e70780. "
        "doi:10.7554/eLife.70780"},
    "lotus_v11": {"verified": "artifact", "text":
        "LOTUS initiative. LOTUS frozen snapshot, version 11, 13 April 2026. Zenodo. "
        "doi:10.5281/zenodo.19360665 (CC BY 4.0)."},
    "sorokina2021": {"verified": "pubmed", "text":
        "Sorokina M, Merseburger P, Rajan K, Yirik MA, Steinbeck C. COCONUT online: "
        "Collection of Open Natural Products database. J Cheminform. 2021;13(1):2. "
        "doi:10.1186/s13321-020-00478-9"},
    "coconut2": {"verified": "author_action", "text":
        "**[AUTHOR ACTION]** COCONUT 2.0 database publication. The statement it "
        "supports — that the COCONUT NP-likeness distribution is multimodal, "
        "\"suggesting the coexistence of both natural products and non-natural "
        "products, or rather synthetic compounds, within the dataset\" — is a verbatim "
        "quotation recorded in this project's decision log (ADR-AIIM-079) at the time "
        "the library source was chosen, but the bibliographic record of the COCONUT 2.0 "
        "paper was not captured then and could not be recovered from PubMed. Supply the "
        "full citation, or delete the sentence in §2.3 and cite only [sorokina2021] for "
        "the resource itself."},
    "mendez2019": {"verified": "v4", "text":
        "Mendez D, Gaulton A, Bento AP, et al. ChEMBL: towards direct deposition of "
        "bioassay data. Nucleic Acids Res. 2019;47(D1):D930–D940. "
        "doi:10.1093/nar/gky1075"},
    "tropsha2010": {"verified": "pubmed", "text":
        "Tropsha A. Best practices for QSAR model development, validation, and "
        "exploitation. Mol Inform. 2010;29(6–7):476–488. doi:10.1002/minf.201000061"},
    "rogers2010": {"verified": "v4", "text":
        "Rogers D, Hahn M. Extended-connectivity fingerprints. J Chem Inf Model. "
        "2010;50(5):742–754. doi:10.1021/ci100050t"},
    "pedregosa2011": {"verified": "v4", "text":
        "Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in "
        "Python. J Mach Learn Res. 2011;12:2825–2830."},
    "bemis1996": {"verified": "v4", "text":
        "Bemis GW, Murcko MA. The properties of known drugs. 1. Molecular frameworks. "
        "J Med Chem. 1996;39(15):2887–2893. doi:10.1021/jm9602928"},
    "rucker2007": {"verified": "pubmed", "text":
        "Rücker C, Rücker G, Meringer M. y-Randomization and its variants in QSPR/QSAR. "
        "J Chem Inf Model. 2007;47(6):2345–2357. doi:10.1021/ci700157b"},
    "norinder2014": {"verified": "pubmed", "text":
        "Norinder U, Carlsson L, Boyer S, Eklund M. Introducing conformal prediction in "
        "predictive modeling. A transparent and flexible alternative to applicability "
        "domain determination. J Chem Inf Model. 2014;54(6):1596–1603. "
        "doi:10.1021/ci5001168"},
    "norinder2015": {"verified": "pubmed", "text":
        "Norinder U, Carlsson L, Boyer S, Eklund M. Introducing conformal prediction in "
        "predictive modeling for regulatory purposes. A transparent and flexible "
        "alternative to applicability domain determination. Regul Toxicol Pharmacol. "
        "2015;71(2):279–284. doi:10.1016/j.yrtph.2014.12.021"},
    "sun2017": {"verified": "pubmed", "text":
        "Sun J, Carlsson L, Ahlberg E, Norinder U, Engkvist O, Chen H. Applying "
        "Mondrian cross-conformal prediction to estimate prediction confidence on large "
        "imbalanced bioactivity data sets. J Chem Inf Model. 2017;57(7):1591–1598. "
        "doi:10.1021/acs.jcim.7b00159"},
    "trott2010": {"verified": "v4", "text":
        "Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking "
        "with a new scoring function, efficient optimization, and multithreading. "
        "J Comput Chem. 2010;31(2):455–461. doi:10.1002/jcc.21334"},
    "eberhardt2021": {"verified": "v4", "text":
        "Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: new "
        "docking methods, expanded force field, and Python bindings. J Chem Inf Model. "
        "2021;61(8):3891–3898. doi:10.1021/acs.jcim.1c00203"},
    "warren2006": {"verified": "pubmed", "text":
        "Warren GL, Andrews CW, Capelli A-M, et al. A critical assessment of docking "
        "programs and scoring functions. J Med Chem. 2006;49(20):5912–5931. "
        "doi:10.1021/jm050362n"},
    "riniker2015": {"verified": "v4", "text":
        "Riniker S, Landrum GA. Better informed distance geometry: using what we know "
        "to improve conformation generation. J Chem Inf Model. 2015;55(12):2562–2574. "
        "doi:10.1021/acs.jcim.5b00654"},
    "halgren1996": {"verified": "v4", "text":
        "Halgren TA. Merck molecular force field. I. Basis, form, scope, "
        "parameterization, and performance of MMFF94. J Comput Chem. "
        "1996;17(5–6):490–519. "
        "doi:10.1002/(SICI)1096-987X(199604)17:5/6<490::AID-JCC1>3.0.CO;2-P"},
    "truchon2007": {"verified": "v4", "text":
        "Truchon J-F, Bayly CI. Evaluating virtual screening methods: good and bad "
        "metrics for the \"early recognition\" problem. J Chem Inf Model. "
        "2007;47(2):488–508. doi:10.1021/ci600426e"},
    "durant2002": {"verified": "pubmed", "text":
        "Durant JL, Leland BA, Henry DR, Nourse JG. Reoptimization of MDL keys for use "
        "in drug discovery. J Chem Inf Comput Sci. 2002;42(6):1273–1280. "
        "doi:10.1021/ci010132r"},
    "gobbi1998": {"verified": "web", "text":
        "Gobbi A, Poppinger D. Genetic optimization of combinatorial libraries. "
        "Biotechnol Bioeng. 1998;61(1):47–54. "
        "doi:10.1002/(SICI)1097-0290(199824)61:1<47::AID-BIT9>3.0.CO;2-Z"},
    "stiefl2006": {"verified": "pubmed", "text":
        "Stiefl N, Watson IA, Baumann K, Zaliani A. ErG: 2D pharmacophore descriptions "
        "for scaffold hopping. J Chem Inf Model. 2006;46(1):208–220. "
        "doi:10.1021/ci050457y"},
    "ji1996": {"verified": "pubmed", "text":
        "Ji XD, Melman N, Jacobson KA. Interactions of flavonoids and other "
        "phytochemicals with adenosine receptors. J Med Chem. 1996;39(3):781–788. "
        "doi:10.1021/jm950661k"},
    "singh2021": {"verified": "pubmed", "text":
        "Singh A, Bhatt G, Gujre N, Mitra S, Swaminathan R, Limaye AM, Rangan L. "
        "Karanjin. Phytochemistry. 2021;183:112641. "
        "doi:10.1016/j.phytochem.2020.112641"},
    "tian2020": {"verified": "pubmed", "text":
        "Tian C, Kasavajhala K, Belfon KAA, et al. ff19SB: amino-acid-specific protein "
        "backbone parameters trained against quantum mechanics energy surfaces in "
        "solution. J Chem Theory Comput. 2020;16(1):528–552. "
        "doi:10.1021/acs.jctc.9b00591"},
    "dickson2022": {"verified": "pubmed", "text":
        "Dickson CJ, Walker RC, Gould IR. Lipid21: complex lipid membrane simulations "
        "with AMBER. J Chem Theory Comput. 2022;18(3):1726–1736. "
        "doi:10.1021/acs.jctc.1c01217"},
    "izadi2014": {"verified": "web", "text":
        "Izadi S, Anandakrishnan R, Onufriev AV. Building water models: a different "
        "approach. J Phys Chem Lett. 2014;5(21):3863–3871. doi:10.1021/jz501780a"},
    "wang2004": {"verified": "v4", "text":
        "Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. Development and testing of "
        "a general amber force field. J Comput Chem. 2004;25(9):1157–1174. "
        "doi:10.1002/jcc.20035"},
    "jakalian2002": {"verified": "v4", "text":
        "Jakalian A, Jack DB, Bayly CI. Fast, efficient generation of high-quality "
        "atomic charges. AM1-BCC model: II. Parameterization and validation. J Comput "
        "Chem. 2002;23(16):1623–1641. doi:10.1002/jcc.10128"},
    "abraham2015": {"verified": "web", "text":
        "Abraham MJ, Murtola T, Schulz R, Páll S, Smith JC, Hess B, Lindahl E. "
        "GROMACS: high performance molecular simulations through multi-level "
        "parallelism from laptops to supercomputers. SoftwareX. 2015;1–2:19–25. "
        "doi:10.1016/j.softx.2015.06.001"},
    "schottverdugo2019": {"verified": "web", "text":
        "Schott-Verdugo S, Gohlke H. PACKMOL-Memgen: a simple-to-use, generalized "
        "workflow for membrane-protein–lipid-bilayer system building. J Chem Inf "
        "Model. 2019;59(6):2522–2528. doi:10.1021/acs.jcim.9b00269"},
    "ng2014": {"verified": "pubmed", "text":
        "Ng HW, Laughton CA, Doughty SW. Molecular dynamics simulations of the "
        "adenosine A2a receptor in POPC and POPE lipid bilayers: effects of membrane "
        "on protein behavior. J Chem Inf Model. 2014;54(2):573–581. "
        "doi:10.1021/ci400463z"},
    "mcnutt2025": {"verified": "pubmed", "text":
        "McNutt AT, Li Y, Meli R, Aggarwal R, Koes DR. GNINA 1.3: the next increment in "
        "molecular docking with deep learning. J Cheminform. 2025;17(1):28. "
        "doi:10.1186/s13321-025-00973-x"},
    "corso2023": {"verified": "arxiv+iclr+dblp", "text":
        "Corso G, Stärk H, Jing B, Barzilay R, Jaakkola TS. DiffDock: diffusion steps, "
        "twists, and turns for molecular docking. In: International Conference on "
        "Learning Representations (ICLR 2023). arXiv:2210.01776"},
    "corso2024": {"verified": "arxiv+iclr+openreview", "text":
        "Corso G, Deng A, Fry B, Polizzi N, Barzilay R, Jaakkola T. Deep confident steps "
        "to new pockets: strategies for docking generalization. In: International "
        "Conference on Learning Representations (ICLR 2024). arXiv:2402.18396"},
}

# --------------------------------------------------------------------------
# Citation insertions: (anchor, replacement). Anchors are exact substrings of
# the manuscript with newlines normalised to single spaces at match time.
# --------------------------------------------------------------------------
CITES: list[tuple[str, str]] = [
    # ---- Introduction ----
    ("interest, and the adenosine A2A receptor (A2AR) is the most clinically advanced among them:",
     "interest, and the adenosine A2A receptor (A2AR) is the most clinically advanced "
     "among them {schwarzschild2006}:"),
    ('experiencing "OFF" episodes [REF: Ohno 2023 Mol Pharmacol, 10.1124/molpharm.122.000633].',
     'experiencing "OFF" episodes {ohno2023,kondo2015}.'),
    ("receptor with multiple co-crystal structures, which makes it unusually well suited to "
     "structure-based computational screening.",
     "receptor with multiple co-crystal structures {liu2012}, which makes it unusually "
     "well suited to structure-based computational screening."),
    ("can misdirect that allocation while producing no visible symptom.",
     "can misdirect that allocation while producing no visible symptom. The parallel "
     "with machine learning in clinical prediction is exact, and instructive: there, "
     "systematic review found that the overwhelming majority of published models "
     "carried methodological faults invisible in the reported result "
     "{wynants2020,roberts2021}, that models learned confounders rather than signal "
     "{degrave2021}, and that data leakage propagated across whole literatures "
     "{kapoor2023}."),
    # ---- 2.1 ----
    ("operating under nine sequential gates, G0 (input provenance)",
     "{claudecode,react2023} operating under nine sequential gates, G0 (input provenance)"),
    # ---- 2.2 ----
    ("giving 133 in total, and detection is measured as the complement of the silent "
     "error escape rate (SEER)",
     "giving 133 in total. The design follows established mutation-testing practice "
     "{jia2011}, and detection is measured as the complement of the silent error escape "
     "rate (SEER)"),
    ("It needs no oracle, which is what makes it usable on agent-written code:",
     "It needs no oracle {chen2018}, which is what makes it usable on agent-written code:"),
    # ---- 2.3 ----
    ("The screening library comprises 236 natural products isolated from Thai orchid "
     "species (SNP-Lib).",
     "The screening library comprises 236 natural products isolated from Thai orchid "
     "species (SNP-Lib) {snplib}."),
    ("Benchmark sets are drawn from ChEMBL A2AR bioactivity records and from DUD-E "
     "[REF: Mysinger 2012, 10.1021/jm300687e].",
     "Benchmark sets are drawn from ChEMBL A2AR bioactivity records {mendez2019} and "
     "from DUD-E {mysinger2012}."),
    ("Curation follows established cheminformatics practice [REF: Fourches 2010, "
     "10.1021/ci100176x; Fourches 2016, 10.1021/acs.jcim.6b00129]:",
     "Curation follows established cheminformatics practice {fourches2010,fourches2016}, "
     "implemented in RDKit {rdkit}:"),
    ("The expansion library derives from LOTUS (version 11, 13 April 2026; Zenodo DOI "
     "10.5281/zenodo.19360665; CC BY 4.0), selected over COCONUT after the latter's own "
     "publication reported synthetic-compound contamination in its aggregated records.",
     "The expansion library derives from LOTUS {rutz2022,lotus_v11} (version 11, "
     "13 April 2026; CC BY 4.0), selected over COCONUT {sorokina2021} after the latter's "
     "own publication reported synthetic-compound contamination in its aggregated "
     "records {coconut2}. Property filtering applied conventional drug-likeness and "
     "screening-artefact criteria {lipinski1997,veber2002,brenk2008,baell2010,bickerton2012}."),
    # ---- 2.4 ----
    ("[REF: Tropsha 2010, 10.1002/minf.201000061], model selection used nested",
     "{tropsha2010}, model selection used nested"),
    ("Following Rücker et al. (10.1021/ci700157b), the activity labels were permuted",
     "Following Rücker et al. {rucker2007}, the activity labels were permuted"),
    ("alternative to conventional applicability-domain determination [REF: "
     "Norinder 2015, 10.1016/j.yrtph.2014.12.021].",
     "alternative to conventional applicability-domain determination "
     "{norinder2014,norinder2015}."),
    # ---- 2.5 ----
    ("Docking used AutoDock Vina [REF: Eberhardt 2021, PMID 34278794].",
     "Docking used AutoDock Vina {trott2010,eberhardt2021}."),
    ("the threshold conventionally used in large comparative assessments of "
     "binding-mode prediction [REF: Warren 2006, 10.1021/jm050362n].",
     "the threshold conventionally used in large comparative assessments of "
     "binding-mode prediction {warren2006}."),
    ("Ligand conformers were generated with ETKDGv3.",
     "Ligand conformers were generated with ETKDGv3 {riniker2015} and minimised with "
     "MMFF94 {halgren1996}."),
    ("as ROC-AUC with bootstrap confidence intervals, enrichment factor and BEDROC20.",
     "as ROC-AUC with bootstrap confidence intervals, enrichment factor and BEDROC20 "
     "{truchon2007}."),
    # ---- 4. Discussion ----
    ("alternative to applicability-domain determination [REF: Norinder 2015].",
     "alternative to applicability-domain determination {norinder2015}."),
    ("whether the findings replicate under an independent engine such as GNINA [REF: "
     "McNutt 2025, 10.1186/s13321-025-00973-x] or a generative pose predictor was not "
     "tested.",
     "whether the findings replicate under an independent engine such as GNINA "
     "{mcnutt2025} or a generative pose predictor was not tested."),
]

# Citations that must be attached to text this script also has to locate by a
# looser rule, because the surrounding prose was written in several passes.
LOOSE_CITES: list[tuple[str, str]] = [
    ("one fixed rule was applied identically in eight",
     "one fixed rule was applied identically in eight",
     "{rogers2010,durant2002,gobbi1998,stiefl2006}"),
    ("scaffold-grouped cross-validation", "scaffold-grouped cross-validation",
     "{bemis1996,pedregosa2011}"),
    ("Ji et al. (J Med Chem 1996;39:781-8)", "Ji et al. {ji1996}", None),
    ("(Singh et al., Phytochemistry 2021;183:112641)", "{singh2021}", None),
    ("class-conditional (Mondrian) inductive conformal prediction",
     "class-conditional (Mondrian) inductive conformal prediction {sun2017}", None),
]


def _replace_flexible(text: str, anchor: str, repl: str) -> tuple[str, bool]:
    """Replace `anchor`, allowing any run of whitespace where it has whitespace.

    The manuscript is hard-wrapped, so an anchor written on one line here has
    to match text broken across two there. Matching on a whitespace-flexible
    regex is the only way to keep the anchors in this file readable.
    """
    pat = re.compile(r"\s+".join(re.escape(w) for w in anchor.split()))
    m = pat.search(text)
    if not m:
        return text, False
    return text[:m.start()] + repl + text[m.end():], True


def apply_citations(text: str) -> tuple[str, list[str]]:
    problems: list[str] = []
    for anchor, repl in CITES:
        text, ok = _replace_flexible(text, anchor, repl)
        if not ok:
            problems.append(f"ANCHOR NOT FOUND: {anchor[:70]}...")
    for anchor, repl, suffix in LOOSE_CITES:
        full = repl if suffix is None else f"{repl} {suffix}"
        text, ok = _replace_flexible(text, anchor, full)
        if not ok:
            problems.append(f"LOOSE ANCHOR NOT FOUND: {anchor[:70]}...")
    return text, problems


def number_and_render(text: str) -> tuple[str, list[tuple[int, str, str]]]:
    """Replace {key,key} markers with [n], numbered by order of first appearance."""
    order: list[str] = []
    marker = re.compile(r"\{([a-z0-9_,]+)\}")

    def keys_of(m: re.Match) -> list[str] | None:
        """The manuscript legitimately contains brace expressions that are not
        citations -- the conformal prediction sets are written {inactive} and
        {inactive, active}. A match is a citation only if every comma-separated
        token is a known reference key; anything else is left alone."""
        ks = [k.strip() for k in m.group(1).split(",")]
        return ks if ks and all(k in REFS for k in ks) else None

    def first_pass(m: re.Match) -> str:
        ks = keys_of(m)
        if ks:
            for k in ks:
                if k not in order:
                    order.append(k)
        return m.group(0)

    marker.sub(first_pass, text)
    num = {k: i + 1 for i, k in enumerate(order)}

    def second_pass(m: re.Match) -> str:
        ks = keys_of(m)
        if not ks:
            return m.group(0)
        ns = sorted(num[k] for k in ks)
        # Collapse runs of three or more into a range, per Elsevier style.
        out, i = [], 0
        while i < len(ns):
            j = i
            while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
                j += 1
            out.append(f"{ns[i]}-{ns[j]}" if j - i >= 2
                       else ",".join(str(x) for x in ns[i:j + 1]))
            i = j + 1
        return "[" + ",".join(out) + "]"

    text = marker.sub(second_pass, text)
    listing = [(num[k], k, REFS[k]["text"]) for k in order]
    unused = [k for k in REFS if k not in num]
    if unused:
        print(f"NOTE: {len(unused)} defined references are never cited and are "
              f"omitted from the list: {', '.join(sorted(unused))}")
    return text, listing


def build_section(listing: list[tuple[int, str, str]]) -> str:
    lines = [
        "## References",
        "",
        "*Provenance of every entry is recorded in "
        "`10_new_system/scripts/build_references_v6.py`: `pubmed` = record retrieved "
        "from PubMed for this manuscript; `v4` = carried unchanged from the previous "
        "submission's reference list; `web` = verified against the publisher's DOI "
        "record; `author_action` = the supported claim is in the project record but "
        "the bibliographic entry could not be recovered and must be supplied.*",
        "",
    ]
    for n, _key, txt in listing:
        lines.append(f"[{n}] {txt}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    text = MS.read_text(encoding="utf-8")
    text, problems = apply_citations(text)
    if problems:
        print("PROBLEMS:")
        for p in problems:
            print("  ", p)
        return 1

    text, listing = number_and_render(text)
    section = build_section(listing)

    # A numbered-style Elsevier article carries References after the
    # Declarations block and before the tables and figures.
    anchor = "\n## Tables\n"
    if anchor not in text:
        raise SystemExit("could not find the Tables heading to insert before")
    text = text.replace(anchor, "\n---\n\n" + section + "\n---\n" + anchor, 1)

    if args.check:
        print("check only; nothing written")
    else:
        MS.write_text(text, encoding="utf-8")
        print(f"written -> {MS}")

    print(f"\ncited references : {len(listing)}")
    by_src: dict[str, int] = {}
    for _n, k, _t in listing:
        v = REFS[k]["verified"]
        by_src[v] = by_src.get(v, 0) + 1
    for k, v in sorted(by_src.items()):
        print(f"  {k:15s} {v}")
    n_marks = len(re.findall(r"\[\d+(?:[,-]\d+)*\]", text))
    print(f"in-text citation markers : {n_marks}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
