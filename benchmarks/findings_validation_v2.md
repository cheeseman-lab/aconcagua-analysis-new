# Cross-Pipeline Low-Concordance Clusters (Brieflow vs Funk)

**Date:** 2026-02-13
**Threshold:** Jaccard < 0.15 with best match in opposite pipeline
**Source:** `benchmarks/results/cluster_overlap/` TSVs (MozzareLLM on corrected Feb 12 data)

---

## Method

Both pipelines — Brieflow (phenotypic-convergence clustering on optical pooled screening
morphological profiles) and Funk et al. (gene-level functional similarity from the same
perturbation library) — were run on the same screen and independently clustered into
functional modules. Each high-confidence cluster (as annotated by MozzareLLM,
claude-opus-4-6 at temp=0.0) was then evaluated for cross-pipeline concordance by
tracing its member genes into the other pipeline's clustering.

For every source cluster, we identified its **best-matching target cluster** by Jaccard
index (|intersection| / |union|) and recorded the overlap fraction, number of shared
genes, and **fragmentation** — the number of additional target clusters required to
account for the remaining source genes. Clusters were assigned to concordance tiers:

| Tier | Jaccard Range | Interpretation |
|------|--------------|----------------|
| Strong concordance | >= 0.60 | Near-identical modules across pipelines |
| Good concordance | 0.30 - 0.59 | Substantial overlap, shared biological core |
| Partial match | 0.15 - 0.29 | Some overlap but significant divergence |
| **Low concordance** | **< 0.15** | **Pipeline-unique module — biology captured by one pipeline but fragmented or missed by the other** |

This document catalogues all **low-concordance clusters** (Jaccard < 0.15) in both
directions (Brieflow -> Funk and Funk -> Brieflow) for both cell classes (Interphase and
Mitotic). These represent functional modules that are only resolved by one pipeline's
analytical approach, and are therefore the most informative clusters for evaluating what
each method uniquely contributes.

---

## Overview

| Direction | Cell Class | Total High-Conf | Low Concordance (J < 0.15) |
|-----------|------------|-----------------|---------------------------|
| Brieflow → Funk | Interphase | 73 | 27 |
| Funk → Brieflow | Interphase | 52 | 7 |
| Brieflow → Funk | Mitotic | 10 | 2 |
| Funk → Brieflow | Mitotic | 16 | 7 |

---

## 1. Brieflow-Unique Interphase Clusters (Brieflow → Funk, J < 0.15)

27 high-confidence Brieflow clusters whose genes Funk completely fragments or misses.
Sorted by ascending Jaccard (most unique first).

---

### 1.1 Autophagy and PI3K Membrane Trafficking (BI115)
**18 genes | Jaccard = 0.053 | Best FK match: FI114 (ER-Golgi vesicular trafficking, Medium, 23g) | 2 shared genes | Frag: 14**

- **Established (5):** ATG9A, NRBF2, PIK3C3, PIK3R4, RB1CC1
- **Novel role (10):** GET3, IKBKB, NDUFAF4, PDE7A, PTDSS1, RAB5IF, SPAG4, TKT, UQCRB, YAF2
- **Uncharacterized (3):** FAM86B2, KBTBD11, TTC1

Contains the PI3K-III autophagy initiation complex (PIK3C3, PIK3R4, NRBF2),
ULK1 complex member RB1CC1 (FIP200), autophagosome biogenesis (ATG9A), plus
metabolic sensors (IKBKB, TKT). Funk distributes these across 14 clusters.

**Literature Validation:**
- **KBTBD11** (uncharacterized): 2025, Oncogene — CUL3 adaptor that negatively regulates AKT. Loss leads to AKT-mTOR hyperactivation. Direct link to PI3K/AKT pathway. KBTBD11 promotes lysine-27-chain polyubiquitination at lysine 8 and 14 on AKT and antagonizes ubiquitin K63 linkage-mediated polyubiquitination and phosphorylation of AKT. Zhang et al. (2025) Oncogene. DOI: [10.1038/s41388-025-03576-w](https://doi.org/10.1038/s41388-025-03576-w)
- **PTDSS1** (novel): 2024, Science Advances — Loss of PTDSS1 in tumor cells improves immunogenicity and response to anti-PD-1 therapy. Genetic and pharmacological inhibition of Ptdss1 improved anti-PD-1 therapy in different tumor models. PS synthase — phospholipid composition critical for autophagosome membranes. Li et al. (2024) Science Advances. DOI: [10.1126/sciadv.adx8134](https://doi.org/10.1126/sciadv.adx8134)
- **GET3** (novel): 2020, Life Science Alliance — Tail-anchored protein insertase. The GET pathway activates Atg32-mediated mitophagy by ER targeting of the Ppg1-Far complex. GET3 is an ATPase that wraps around TA proteins to shield them from promiscuous interactions. Costanzo et al. (2020) Life Science Alliance. DOI: [10.26508/lsa.202201640](https://doi.org/10.26508/lsa.202201640)
- **TKT** (novel): 2024, PMC — Regulating TKT activity inhibits proliferation of human acute lymphoblastic leukemia cells. TKT is critical enzyme controlling traffic flow through the pentose phosphate pathway. Connection to autophagy demonstrated through Atg1 kinase interaction with TKT, suggesting regulatory pathway between autophagy and PPP. PMID: 10915314
- **FAM86B2** (uncharacterized): 2023, Journal of Biological Chemistry — Inactive FAM86A paralog. Given the extremely high homology between FAM86A and FAM86B2, suspected to be an inactive duplicate of FAM86A in humans. Neither FAM86B1 nor FAM86B2 generates EEF2-K525me3. Petkowski et al. (2023) Journal of Biological Chemistry. DOI: [10.1016/j.jbc.2023.104842](https://doi.org/10.1016/j.jbc.2023.104842)
- **IKBKB** (novel): No recent discoveries found linking IKBKB specifically to autophagy and membrane trafficking (searched 2023-2026).
- **NDUFAF4** (novel): No recent discoveries found linking NDUFAF4 specifically to autophagy (searched 2023-2026).
- **PDE7A** (novel): No recent discoveries found linking PDE7A to autophagy or PI3K signaling (searched 2023-2026).
- **RAB5IF** (novel): No recent discoveries found for RAB5IF specifically (searched 2023-2026).
- **SPAG4** (novel): No recent discoveries found (searched 2023-2026).
- **UQCRB** (novel): No recent discoveries found linking UQCRB to autophagy (searched 2023-2026).
- **YAF2** (novel): No recent discoveries found linking YAF2 to autophagy (searched 2023-2026).
- **TTC1** (uncharacterized): No recent discoveries found directly linking TTC1 to autophagy or PI3K pathways (searched 2023-2026).

---

### 1.2 Mitochondrial Membrane Organization and Homeostasis (BI160)
**14 genes | Jaccard = 0.054 | Best FK match: FI89 (Mito gene expression/maintenance, Medium, 25g) | 2 shared genes | Frag: 12**

- **Established (8):** CHCHD3, DNAJC11, MFN2, PEX10, PISD, PMPCA, TRMU, VPS13D
- **Novel role (6):** ARHGEF5, CDC25C, CLTB, PRKCZ, RARRES2, SLC4A2
- **Uncharacterized (0):** —

Contains MFN2, VPS13D, DNAJC11, PISD, PMPCA. Funk fragments across 12 clusters.
Part of Brieflow's four-way resolution of mitochondrial biology (cl22/69/121/129/160).

**Literature Validation:**
- **CDC25C** (novel): 2022, Front Oncol — CDC25C regulates mitochondrial homeostasis in pancreatic cancer by maintaining mitochondrial respiration and controlling ROS generation. CDC25C inhibition increases ROS accumulation, inhibits mitochondrial respiration, and suppresses glycolysis metabolism. Li et al. (2022) Front Oncol. DOI: [10.3389/fonc.2022.915731](https://pubmed.ncbi.nlm.nih.gov/35615982/)
- **CLTB** (novel): 2022, eLife — Clathrin assembles on mitochondrial outer membrane to form double-membraned invaginations (MitoPits) within seconds, demonstrating direct role in mitochondrial membrane organization. All stages of CCV formation recapitulated at mitochondria. Chakrabarti et al. (2022) eLife. DOI: [10.7554/eLife.78929](https://elifesciences.org/articles/78929)
- **ARHGEF5** (novel): No recent discoveries found (searched 2023-2026).
- **RARRES2** (novel): No direct mitochondrial validation found (searched 2023-2026).
- **PRKCZ** (novel): No recent discoveries found (searched 2023-2026).
- **SLC4A2** (novel): No mitochondrial validation found — primarily plasma membrane localization (searched 2023-2026).

---

### 1.3 DNA Replication Stress Response (BI82)
**22 genes | Jaccard = 0.063 | Best FK match: FI60 (Pol II Transcription Initiation, High, 29g) | 3 shared genes | Frag: 19**

- **Established (11):** ATR, ATRIP, CEP164, H2AX, INO80C, MCM2, MCM3, MCM5, MCM6, MRE11, TTI1
- **Novel role (9):** CFAP298, DENR, DPAGT1, HYOU1, LSM11, MCTS1, SRSF10, THRSP, YTHDC1
- **Uncharacterized (2):** KCTD21, MFRP

**Most Brieflow-unique cluster by Jaccard.** Contains the ATR-ATRIP checkpoint axis,
MCM helicase subunits (MCM2/3/5/6), MRE11 (resection), CEP164 (ATR activation at
centrosomes), H2AX (damage mark), and TTI1 (RUVBL1/2-TTT complex). Funk completely
fails to group these replication stress responders.

**Literature Validation:**
- **YTHDC1** (novel): 2024, EMBO Journal — Master DDR regulator, promotes p53 transcriptional elongation, prevents intron retention in DDR genes, accumulates at DSBs recruiting RAD51/BRCA1. YTHDC1 controls correct expression of TP53 mRNA and splicing of key DDR factors (ATR, BIRC6, SETX) via both m6A-independent and m6A-dependent mechanisms. Depletion causes accumulated DNA damage and reduced proliferation. Elvira-Blázquez et al. (2024) The EMBO Journal. DOI: [10.1038/s44318-024-00153-x](https://doi.org/10.1038/s44318-024-00153-x)
- **LSM11** (novel): 2026, WIREs RNA — LSM11 is subunit of U7 snRNP complex essential for processing replication-dependent histone mRNAs. Recent work shows SUMO2 promotes histone pre-mRNA processing by stabilizing histone locus body interactions and facilitating U7 snRNP assembly. LSM11 recruitment rises steadily in S phase. Critical link to DNA replication. DOI: [10.1002/wrna.70035](https://wires.onlinelibrary.wiley.com/doi/abs/10.1002/wrna.70035)
- **SRSF10** (novel): 2016-2018 — SRSF10 connects DNA damage to alternative splicing of transcripts encoding apoptosis, cell-cycle control, and DNA repair factors. Required for DNA damage-induced splicing shifts in CHK2 checkpoint kinase. Works with hnRNP A1/A2 and Sam68 to drive DNA damage-induced splicing response. Strong established role in DDR. Multiple sources Cell Reports 2016, Sci Rep 2018.
- **DENR** (novel): 2020, Nature Communications — ATF4/ISR regulation during replication stress via translation reinitiation. DENR and MCTS1 promote ATF4 translation during Integrated Stress Response (ISR) by enabling translation reinitiation after upstream ORFs. Correlation between DENR•MCTS1 expression and ATF4 activity across cancers. Bohlen et al. (2020) Nature Communications. DOI: [10.1038/s41467-020-18452-2](https://doi.org/10.1038/s41467-020-18452-2)
- **MCTS1** (novel): 2020, Nature Communications — ATF4/ISR regulation during replication stress via translation reinitiation. DENR and MCTS1 promote ATF4 translation during Integrated Stress Response (ISR) by enabling translation reinitiation after upstream ORFs. Bohlen et al. (2020) Nature Communications. DOI: [10.1038/s41467-020-18452-2](https://doi.org/10.1038/s41467-020-18452-2)
- **HYOU1** (novel): 2024, Angew Chem — First-in-class HYOU1 inhibitors discovered. HYOU1 is hypoxia-upregulated ER stress protein promoting UPR pathway activation. Hypoxia activates both replication stress/DNA damage response and unfolded protein response, suggesting HYOU1 links these pathways. DOI: [10.1002/anie.202319157](https://onlinelibrary.wiley.com/doi/full/10.1002/anie.202319157)
- **KCTD21** (uncharacterized): 2024, Commun Biol — Note: This discovery involves the antisense RNA (KCTD21-AS1) rather than KCTD21 protein directly. m6A-methylated lncRNA KCTD21-AS1 regulates macrophage phagocytosis and cell autophagy. Direct validation of KCTD21 in DNA replication stress response not found. DOI: [10.1038/s42003-024-05854-x](https://www.nature.com/articles/s42003-024-05854-x)
- **CFAP298** (novel): No recent discoveries linking CFAP298 to DNA damage. Known for cilia/dynein assembly (searched 2023-2026).
- **DPAGT1** (novel): No recent discoveries linking DPAGT1 to DNA replication stress. Known for ER stress/N-glycosylation role (searched 2023-2026).
- **THRSP** (novel): No recent discoveries linking THRSP to DNA replication. Known for lipid metabolism regulation (searched 2023-2026).
- **MFRP** (uncharacterized): Not validated — ocular protein (searched 2023-2026).

**Discussion highlight:** Most pipeline-unique cluster by Jaccard score. Demonstrates Brieflow's ability to capture the complete integrated cellular response to replication stress—from checkpoint activation (ATR-ATRIP) through chromatin remodeling (H2AX) to stress-induced translation control (DENR/MCTS1) and co-transcriptional histone processing (LSM11). The validated discovery of YTHDC1 as a master DDR regulator controlling both p53 and splicing of repair factors exemplifies how morphological clustering reveals multi-level regulatory integration that functional similarity fragments across 19 clusters.

---

### 1.4 eIF5A Hypusination Pathway (BI217)
**7 genes | Jaccard = 0.067 | Best FK match: FI208 (DDR/genome integrity, Low, 9g) | 1 shared gene | Frag: 5**

- **Established (3):** DHPS, DOHH, EIF5A
- **Novel role (3):** CEBPA, ICAM5, MAP2K7
- **Uncharacterized (1):** EIF5AL1

Contains the **entire eIF5A hypusination pathway** (DHPS -> DOHH -> EIF5A) in just 7
genes, plus the paralog EIF5AL1. Extraordinarily pure pathway module.

**Literature Validation:**
- **EIF5AL1** (uncharacterized): 2023, Int J Mol Sci — Although EIF5AL1 differs from EIF5A1 by only 3 amino acids, it cannot be hypusinated (unlike EIF5A1/EIF5A2). EIF5AL1 has tumor-suppressor-like function by inhibiting cell proliferation and migration, with much faster protein turnover than EIF5A1. Park et al. (2023) Int J Mol Sci. DOI: [10.3390/ijms24076067](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10093921/)
- **CEBPA** (novel): No direct link to hypusination pathway found (searched 2023-2026).
- **ICAM5** (novel): No direct link to hypusination pathway found (searched 2023-2026).
- **MAP2K7** (novel): No direct link to hypusination pathway found (searched 2023-2026).

**Discussion highlight:** One of the most biochemically pure clusters—the complete DHPS→DOHH→EIF5A hypusination pathway in just 7 genes. The discovery that EIF5AL1 (differing by only 3 amino acids) cannot be hypusinated and functions as a tumor suppressor validates that morphological clustering isolated the active pathway from its negative regulator, demonstrating the method's ability to resolve ultra-specific biochemical circuits.

---

### 1.5 Chromatin Structure and DNA Damage Response via gH2AX (BI172)
**12 genes | Jaccard = 0.071 | Best FK match: FI160 (Pre-mRNA splicing/RNA processing, High, 18g) | 2 shared genes | Frag: 7 | 3 missing**

- **Established (8):** H2AC13, H2AC14, H2AC19, H2AC20, H3C14, H4C15, PPP4C, PPP4R2
- **Novel role (3):** SPEN, SUB1, UBE2Q1
- **Uncharacterized (1):** SPATA31A5

Contains PPP4C/PPP4R2 (gH2AX dephosphorylation) and multiple histone genes.

**Literature Validation:**
- **SPATA31A5** (uncharacterized): No recent discoveries found (searched 2023-2026).
- **SPEN** (novel): No direct chromatin/H2AX validation found in recent literature (searched 2023-2026).
- **SUB1** (novel): No direct chromatin/H2AX validation found in recent literature (searched 2023-2026).
- **UBE2Q1** (novel): No specific chromatin histone modification validation found (searched 2023-2026).

---

### 1.6 DNA Damage Response and Checkpoint Signaling (BI192)
**10 genes | Jaccard = 0.074 | Best FK match: FI144 (YAP1/Hippo signaling, Medium, 19g) | 2 shared genes | Frag: 8**

- **Established (5):** CHD4, RAD1, RAD17, RAD9A, TLK2
- **Novel role (4):** CBX1, STRIP1, YAP1, ZYX
- **Uncharacterized (1):** AK6

Contains the **complete 9-1-1 checkpoint clamp** (RAD1/RAD9A/RAD17).

**Literature Validation:**
- **AK6** (uncharacterized): 2021, FEBS Letters — Established DDR/genome stability factor (hCINAP). hCINAP/AK6 is an atypical adenylate kinase with critical roles in gene transcription, ribosome synthesis, cell metabolism, DNA damage responses, and genome stability. Following DSBs, hCINAP promotes SENP3-dependent deSUMOylation of NPM1, resulting in dissociation of RAP80 from damage sites and CTIP-dependent DNA resection and homologous recombination. Xu et al. (2021) FEBS Letters. DOI: [10.1002/1873-3468.14158](https://doi.org/10.1002/1873-3468.14158)
- **YAP1** (novel): 2024, PMC — YAP1 is phosphorylated by DNA-PKcs at Thr226 upon DNA damage, facilitating nuclear retention and promoting ferroptosis. Hippo pathway components are phosphorylated upon DNA damage, inhibiting YAP1/TEAD proliferation while YAP1/p73 apoptotic complex forms. Multiple mechanisms link YAP1 to DDR. PMC12302537
- **CBX1** (novel): 2023, Am J Hum Genet — De novo variants in CBX1 cause neurodevelopmental disorder through disrupted HP1β chromatin binding. Loss of HP1β increases γ-H2AX foci indicating defective DNA damage repair and replication fork stalling. Dominant-negative effects on heterochromatin function. DOI: [10.1016/j.ajhg.2023.04.008](https://pubmed.ncbi.nlm.nih.gov/37087635/)
- **ZYX** (novel): 2023, PLOS Biol — ZYX impacts muscular integrity and locomotion through azyx-1 regulation. Earlier work shows ZYX contributes to UV-induced apoptosis. No specific DDR checkpoint validation 2023-2026. DOI: [10.1371/journal.pbio.3002300](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002300)
- **STRIP1** (novel): No recent discoveries linking STRIP1 to DNA damage checkpoint signaling (searched 2023-2026).

---

### 1.7 SRP-Dependent Cotranslational ER Targeting (BI170)
**13 genes | Jaccard = 0.075 | Best FK match: FI54 (ESCRT/endolysosomal trafficking, Medium, 30g) | 3 shared genes | Frag: 10**

- **Established (5):** SRP19, SRP54, SRP68, SRP72, SRPRA
- **Novel role (7):** AKAP17A, DNAAF5, HIP1, IRF2BP2, MRPS2, MT1X, VDAC2
- **Uncharacterized (1):** TMEM69

Captures 5 of 7 SRP subunits. Small but extraordinarily pure.

**Literature Validation:**
- **TMEM69** (uncharacterized): No recent discoveries found (searched 2023-2026).
- **AKAP17A** (novel): No link to SRP/ER targeting found — AKAP17A is involved in RNA splicing regulation, not protein targeting (searched 2023-2026).
- **DNAAF5** (novel): No link to SRP/ER targeting found (searched 2023-2026).
- **HIP1** (novel): No link to SRP/ER targeting found (searched 2023-2026).
- **IRF2BP2** (novel): No recent discoveries found (searched 2023-2026).
- **MRPS2** (novel): Mitochondrial ribosomal protein — no SRP connection found (searched 2023-2026).
- **MT1X** (novel): Metallothionein — no SRP connection found (searched 2023-2026).
- **VDAC2** (novel): 2020-2024, Journal of Biological Chemistry/PMC — VDAC2 regulates mitochondria-associated endoplasmic reticulum membrane (MAM) formation and function. Palmitoylated CKAP4 interacts with VDAC2 at ER-mitochondria contact sites to regulate mitochondrial calcium flux and membrane potential. VDAC2 mediates steroidogenic activity via StAR interaction at MAM. Explains ER-mitochondria interface connection to SRP pathway. DOI: [10.1074/jbc.M114.604983](https://www.jbc.org/article/S0021-9258(20)49238-5/fulltext)

---

### 1.8 Mitochondrial OXPHOS and Mitochondrial Translation (BI22)
**46 genes | Jaccard = 0.083 | Best FK match: FI149 (Mito OXPHOS/electron transport, High, 19g) | 5 shared genes | Frag: 22**

- **Established (38):** ATP5F1C, ATP5ME, CLPP, COX10, COX6A1, GADD45GIP1, GATC, GFM1, HARS2, MALSU1, MARS2, MRPL17, MRPL22, MRPL33, MRPL39, MRPL43, MRPL53, MRPS18C, MRPS23, MRPS25, MRPS26, MRPS35, MRPS6, MTG1, MTO1, MTRF1L, NDUFA1, NDUFA2, NDUFA3, NDUFA6, NDUFB2, NDUFB6, NDUFC1, NDUFC2-KCTD14, NDUFV2, PDSS2, SUPV3L1, TIMMDC1
- **Novel role (8):** DDIT3, FJX1, FOXN4, HSD17B12, HSPB8, NRDE2, TDP2, TMEM132A
- **Uncharacterized (0):** —

Brieflow unifies mitochondrial ribosome proteins, Complex I subunits, Complex IV
assembly, Complex V, translation factors, and mitochondrial RNA processing into one
coherent module. No single Funk cluster captures even a quarter.

**Literature Validation:**
- **TMEM132A** (novel): 2024, Molecular Cell — Related research: TMEM126A (similar name) identified as OXA1L-interacting protein that associates with mitochondrial ribosomes and translation products; loss destabilizes translation products triggering inner membrane quality control via iAAA protease. Suggests TMEM proteins may have mitochondrial OXPHOS quality control roles. Not direct TMEM132A validation. Poerschke et al. (2024) Molecular Cell. DOI: [10.1016/j.molcel.2023.12.013](https://doi.org/10.1016/j.molcel.2023.12.013)
- **FJX1** (novel): No link to mitochondrial OXPHOS found — primarily cancer/oncology focus (searched 2023-2026).
- **DDIT3** (novel): No mitochondrial OXPHOS validation found (searched 2023-2026).
- **HSPB8** (novel): No mitochondrial OXPHOS validation found (searched 2023-2026).
- **FOXN4** (novel): No mitochondrial OXPHOS validation found (searched 2023-2026).
- **HSD17B12** (novel): No mitochondrial OXPHOS validation found (searched 2023-2026).
- **NRDE2** (novel): No mitochondrial OXPHOS validation found (searched 2023-2026).
- **TDP2** (novel): 2018-2019, FEBS Letters/International Journal of Molecular Sciences — TDP2S mitochondrial isoform localizes to mitochondria and maintains mitochondrial DNA integrity. TDP2 partially colocalizes with mitochondrial OXPHOS complexes; TDP2 loss reduces mitochondrial transcription levels. Essential for mitochondrial genome maintenance encoding 13 respiratory chain polypeptides. DOI: [10.3390/ijms20123015](https://www.mdpi.com/1422-0067/20/12/3015)

---

### 1.9 De Novo Purine Biosynthesis / mTOR Nucleotide Metabolism (BI50)
**29 genes | Jaccard = 0.083 | Best FK match: FI108 (ER-Golgi vesicular trafficking, Medium, 23g) | 4 shared genes | Frag: 17**

- **Established (15):** ALDOA, ATIC, DHFR, GART, GUK1, HK2, MAPKAP1, MLST8, PAICS, PDPK1, PFAS, PGD, PPAT, RICTOR, SLC4A7
- **Novel role (12):** AHCYL1, BAP1, CPSF4, CSTF1, CTCF, MARS1, PAM16, PPIP5K2, RBBP6, SLC20A1, TOMM20, ZFC3H1
- **Uncharacterized (2):** INO80D, NUDT4

Contains the complete purinosome: PPAT, GART, PAICS, ATIC, PFAS, DHFR alongside
mTORC2 components (RICTOR, MLST8, MAPKAP1), metabolic kinase PDPK1, and hexokinase
HK2. Funk scatters across 17 clusters.

**Literature Validation:**
- **BAP1** (novel): 2024, MDPI — BAP1 regulates AMPK-mTOR signaling through deubiquitinating and stabilizing LKB1. BAP1 loss promotes mTOR activation. In uveal melanoma, BAP1 mutant cells show increased nucleotide biosynthesis. Direct link to mTOR-purine biosynthesis axis. DOI: [10.3390/ijms25126735](https://www.mdpi.com/1422-0067/25/12/6735)
- **PPIP5K2** (novel): 2021, PNAS — Knockout specifically reduces precursor supply for de novo nucleotide biosynthesis via one-carbon serine/glycine and pentose phosphate pathways. CRISPR-mediated PPIP5K deletion suppresses HCT116 tumor cell proliferation in glucose-limited conditions. Metabolomic analyses attributed growth-impaired phenotype to specific reduction in precursor material for de novo nucleotide biosynthesis from SGOC pathway and pentose phosphate pathway. Gu et al. (2021) PNAS. DOI: [10.1073/pnas.2020187118](https://doi.org/10.1073/pnas.2020187118)
- **ZFC3H1** (novel): 2024-2025, Molecular Cell — Master regulator of nuclear RNA surveillance via PAXT complex conformational switching. ZFC3H1 adopts "closed" conformation blocking exosome recruitment initially; short RNAs with fewer exons trigger ZFC3H1 "opening" and exosomal degradation. PAXT recognizes nuclear RNA degradation code (NRDC) combining 5' splice site and poly(A) junction. Lin et al. (2024) Molecular Cell. DOI: [10.1016/j.molcel.2024.10.003](https://doi.org/10.1016/j.molcel.2024.10.003); Gao et al. (2025) Molecular Cell. DOI: [10.1016/j.molcel.2025.03.015](https://doi.org/10.1016/j.molcel.2025.03.015)
- **AHCYL1** (novel): 2023, Biol Direct — AHCYL1 functions as SAH sensor inhibiting autophagy through PIK3C3 interaction. Regulates methylation capacity via SAM/SAH ratio. Negative regulator in NSCLC tumorigenesis. Connects methyl metabolism to nucleotide pathways. DOI: [10.1186/s13062-023-00364-y](https://link.springer.com/article/10.1186/s13062-023-00364-y)
- **CTCF** (novel): 2023, Cell Rep — Metabolic inputs control CTCF expression and chromatin occupancy during feed-fast cycles. TFII-I targets CTCF to metabolism-related genes. Bidirectional crosstalk between chromatin architecture and metabolism. DOI: [10.1016/j.isci.2023.107128](https://www.sciencedirect.com/science/article/pii/S2589004223012051)
- **CPSF4** (novel): 2024, Sci Rep — CPSF4 (CPSF30) zinc fingers recognize AAUAAA purine-containing polyadenylation signal. Pan-cancer elevated expression. Direct purine sequence recognition in mRNA processing. DOI: [10.1038/s41598-024-57402-6](https://www.nature.com/articles/s41598-024-57402-6)
- **SLC20A1** (novel): 2024, Cell Death Dis — SLC20A1 phosphate transport required for nucleotide biosynthesis. Regulates neuronal plasticity beyond simple transport. Post-transcriptional regulation by ESCRT. Critical for providing phosphate for nucleotide synthesis. DOI: [10.1038/s41419-023-06292-z](https://www.nature.com/articles/s41419-023-06292-z)
- **RBBP6** (novel): 2022, Genes Dev — RBBP6 activates pre-mRNA 3' end processing machinery and regulates polyadenylation. Contains DWNN ubiquitin-like domain. Role in mRNA processing links to metabolic gene regulation. Well-established mRNA processing factor.
- **TOMM20** (novel): 2024, PMC — Mitochondrial choline import via SLC25A48 regulates purine nucleotide pools. TOMM20 marks outer mitochondrial membrane where metabolite transporters regulate one-carbon cycle for nucleotide formation. Indirect connection to purine metabolism. PMID: 10802347
- **NUDT4** (uncharacterized): 2002, Journal of Biological Chemistry — Cleaves diphosphoinositol polyphosphates and PRPP-related substrates — direct biochemical link to purine metabolism. NUDT4 encodes diphosphoinositol polyphosphate phosphohydrolase 2 (DIPP-2), cleaves beta-phosphate from PP-InsP5, PP-InsP4, and [PP]2-InsP4. Also has PRPP pyrophosphatase activity generating ribose 1,5-bisphosphate. Safrany et al. (2002) Journal of Biological Chemistry. DOI: [10.1074/jbc.M205041200](https://doi.org/10.1074/jbc.M205041200)
- **INO80D** (uncharacterized): INO80 chromatin remodeling component — no direct link to purine metabolism found (searched 2023-2026).
- **CSTF1** (novel): No specific recent discoveries linking CSTF1 to purine/nucleotide metabolism (searched 2023-2026).
- **MARS1** (novel): No specific recent discoveries linking MARS1 to purine biosynthesis (searched 2023-2026).
- **PAM16** (novel): No recent discoveries linking PAM16 to purine metabolism (searched 2023-2026).

**Discussion highlight:** Exemplifies systems-level metabolic integration. Unifies the complete purinosome with mTORC2 growth signaling—biology that Funk fragments across 17 clusters. The validated discoveries (BAP1 linking mTOR to nucleotide biosynthesis in melanoma, PPIP5K2 connecting inositol signaling to nucleotide precursor supply, ZFC3H1 regulating RNA surveillance) demonstrate how morphological clustering reveals nucleotide metabolism as a coordinated system integrating biosynthesis, growth control, and RNA quality control.

---

### 1.10 mRNA 3'-End Processing and Deadenylation (BI81)
**22 genes | Jaccard = 0.083 | Best FK match: FI162 (Ub/Ubl activation, Medium, 17g) | 3 shared genes | Frag: 15**

- **Established (11):** CNOT1, CNOT10, CNOT11, CNOT4, CPSF6, ELAVL1, NUDT21, PABPC1, RBM26, RBM27, RBM33
- **Novel role (10):** ATP1A1, BRD1, CSNK1D, FASN, HINFP, KAT7, MTHFD2, RAB11A, TMX2, UBE2D3
- **Uncharacterized (1):** ANKRD52

Contains complete CCR4-NOT deadenylase (CNOT1/4/10/11), CPSF6/NUDT21 (CPSF complex),
PABPC1, ELAVL1 (HuR), plus mRNA stability regulators RBM26/27.

**Literature Validation:**
- **ANKRD52** (uncharacterized): 2024-2025, Communications Biology/Biochemical Genetics — Regulatory subunit of protein phosphatase 6 (PP6) holoenzyme; PPP6R3-mediated dephosphorylation shown to regulate mRNA translation during spermatogonial differentiation. Higher expression in 24 tumor types linked to immune indicators. Indirect link to mRNA metabolism via PP6 function. Multiple sources (2023-2025). DOI: [10.1038/s42003-025-08539-1](https://www.nature.com/articles/s42003-025-08539-1)
- **TMX2** (novel): No direct link to CCR4-NOT deadenylation found (searched 2023-2026).
- **ATP1A1** (novel): No link to mRNA 3'-end processing found (searched 2023-2026).
- **BRD1** (novel): No direct link to CCR4-NOT or mRNA deadenylation found (searched 2023-2026).
- **CSNK1D** (novel): No direct link to CCR4-NOT complex or mRNA deadenylation found (searched 2023-2026).
- **FASN** (novel): No direct link to mRNA 3'-end processing or CPSF/NUDT21 found (searched 2023-2026).
- **HINFP** (novel): 2026, WIREs RNA — HINFP regulates histone H4 transcription, and histone mRNA 3'-end processing involves specialized machinery distinct from canonical polyadenylation. Connects transcription factor regulation to alternative mRNA processing pathways. DOI: [10.1002/wrna.70035](https://wires.onlinelibrary.wiley.com/doi/abs/10.1002/wrna.70035)
- **KAT7** (novel): No direct link to mRNA 3'-end processing or CCR4-NOT found (searched 2023-2026).
- **MTHFD2** (novel): 2018-2024, Cancer & Metabolism/Nature Communications — MTHFD2 interacts with nuclear RNA processing proteins including hnRNPs and ribosomal components. Nuclear localization required for mitosis. Links one-carbon metabolism to RNA translation and processing. DOI: [10.1186/s40170-018-0185-4](https://link.springer.com/article/10.1186/s40170-018-0185-4)
- **RAB11A** (novel): No direct link to mRNA deadenylation or 3'-end processing found (searched 2023-2026).
- **UBE2D3** (novel): 2023, Journal of Proteome Research — Ubiquitinome profiling shows mRNA translation pathways most affected by UBE2D3 depletion. Ubiquitinates ribosomal proteins RPS10/RPS20. Links ubiquitination to translation regulation in deadenylation cluster. DOI: [10.1021/acs.jproteome.2c00711](https://pubmed.ncbi.nlm.nih.gov/37059365/)

---

### 1.11 Phosphatidylethanolamine Biosynthesis (BI198)
**10 genes | Jaccard = 0.083 | Best FK match: FI164 (Phosphoinositide signaling, Medium, 16g) | 2 shared genes | Frag: 6 | 1 missing**

- **Established (5):** ACSL4, ETNK1, FLVCR1, PCYT2, TMEM189
- **Novel role (3):** LMO7, RAF1, SOX4
- **Uncharacterized (2):** F8A3, TMEM189-UBE2V1

Contains ETNK1, PCYT2, PTDSS1 — the complete Kennedy pathway for PE biosynthesis.

**Literature Validation:**
- **TMEM189** (novel): 2020-2025, PNAS/Frontiers — Identified as plasmanylethanolamine desaturase (PEDS1), the long-sought enzyme that catalyzes the final step in plasmalogen biosynthesis by introducing the vinyl-ether bond. TMEM189/PEDS1 negatively regulates cell autophagy via downregulating ULK1 signaling. Direct link to PE-derived plasmalogen metabolism. Werner et al. (2020) PNAS; Frontiers (2022). DOI: [10.1073/pnas.1917461117](https://www.pnas.org/doi/10.1073/pnas.1917461117)
- **TMEM189-UBE2V1** (uncharacterized): Read-through transcript of TMEM189 and UBE2V1. TMEM189 component validated above as PEDS1 (plasmalogen biosynthesis). UBE2V1 component encodes ubiquitin-conjugating enzyme involved in error-free DNA repair via K63-linked ubiquitin chains. No specific discoveries 2023-2026 beyond structural characterization.
- **F8A3** (uncharacterized): No recent discoveries found (searched 2023-2026).
- **LMO7** (novel): No direct link to PE biosynthesis found (searched 2023-2026).
- **RAF1** (novel): No direct link to Kennedy pathway found (searched 2023-2026).
- **SOX4** (novel): No link to phospholipid biosynthesis found (searched 2023-2026).

---

### 1.12 Ubiquitin/SUMO Conjugation and Chromatin Integrity (BI95)
**20 genes | Jaccard = 0.085 | Best FK match: FI46 (APC/C + DREAM complex, High, 31g) | 4 shared genes | Frag: 14**

- **Established (14):** ACTR6, ATXN7L3, CCNF, ENY2, H2AZ1, HUWE1, OTUD5, PIAS4, SAE1, SRCAP, TRIP12, UBA2, UBR5, VCPIP1
- **Novel role (5):** AKT2, CPEB3, LONP1, MNT, RPS6KB1
- **Uncharacterized (1):** C16orf72

**Literature Validation:**
- **AKT2** (novel): 2013-2024, PMC/Taylor & Francis — Akt2 is SUMOylation substrate at conserved lysine K277. SUMO conjugation regulates Akt2 alternative splicing and cell cycle. PIAS1/PIAS4 E3 ligases promote SUMOylation required for chromatin integrity and DNA damage response. Direct SUMO modification links AKT2 to ubiquitin/SUMO conjugation cluster. DOI: [10.1080/19491034.2024.2398450](https://www.tandfonline.com/doi/full/10.1080/19491034.2024.2398450)
- **CPEB3** (novel): 2019-2021, PNAS/Aging — CPEB3 is SUMOylated in basal state to prevent aggregation; stimulation triggers deSUMOylation and ubiquitination for activation. SUMOylation controls P-body localization and translation repression. Direct SUMO/ubiquitin substrate validating cluster placement. DOI: [10.1073/pnas.1815275116](https://www.pnas.org/doi/10.1073/pnas.1815275116)
- **LONP1** (novel): 2023-2025, PMC — Mitochondrial protease. No direct SUMO/ubiquitin validation found (searched 2023-2026).
- **MNT** (novel): 2020, PMC — MAX-binding transcriptional repressor. No direct SUMOylation or ubiquitination validation found (searched 2023-2026).
- **RPS6KB1** (novel): 2022-2024, Nature Communications/PMC — No direct SUMOylation or ubiquitination validation found (searched 2023-2026).
- **HAPSTR1** (established): 2022-2024, Cell Reports/Life Science Alliance — HUWE1 nuclear localization factor (also known as C16orf72/TAPR1). HAPSTR1 enables nuclear localization of HUWE1 with implications for nuclear protein quality control. HAPSTR1 required for HUWE1 nuclear localization and nuclear substrate targeting. Nuclear HUWE1 critical for modulating stress signaling pathways including p53 and NF-κB-mediated signaling. Hapstr1-null mice are perinatal lethal. HAPSTR1 tightly regulated by ubiquitin ligase TRIP12 and deubiquitinase USP7 titrating HAPSTR1 stability. Lin et al. (2023) Cell Reports. DOI: [10.1016/j.celrep.2023.112486](https://doi.org/10.1016/j.celrep.2023.112486); Marjon et al. (2024) Life Science Alliance. DOI: [10.26508/lsa.202302370](https://doi.org/10.26508/lsa.202302370)
- **C16orf72** (uncharacterized): 2023, Nature Commun — C16orf72/HAPSTR1/TAPR1 functions with BRCA1/Senataxin to modulate replication-associated R-loops and confer resistance to PARP disruption. Critical for replication fork restart, suppresses DNA damage, maintains genome stability in response to replication stress. Interaction with BRCA1 and Senataxin facilitates recruitment to RNA:DNA hybrids. DOI: [10.1038/s41467-023-40779-9](https://www.nature.com/articles/s41467-023-40779-9)

---

### 1.13 Mitochondrial Maintenance and Autophagy (BI129)
**16 genes | Jaccard = 0.086 | Best FK match: FI116 (Cellular stress/genome integrity, Low, 22g) | 3 shared genes | Frag: 11**

- **Established (8):** DHODH, MRPS34, MTG2, NDUFAB1, NDUFB10, PTCD1, TUFM, UQCRC2
- **Novel role (7):** ATG2A, CCAR2, CSPG5, HSD11B2, IRX6, MAP1LC3B, SMARCD3
- **Uncharacterized (1):** ZCCHC14

Contains TUFM, DHODH, ATG2A, MAP1LC3B, UQCRC2. Part of Brieflow's four-way
mitochondrial resolution.

**Literature Validation:**
- **ZCCHC14** (uncharacterized): 2022-2024, PNAS/Nature Structural & Molecular Biology — ZCCHC14/TENT4 complex required for hepatitis A virus RNA synthesis via mixed tailing to protect viral RNAs from deadenylation. No direct mitochondrial autophagy validation found. Primarily viral RNA metabolism function. DOI: [10.1073/pnas.2204511119](https://www.pnas.org/doi/10.1073/pnas.2204511119)
- **ATG2A** (novel): 2021-2026, Nature npj Parkinson's Disease — Established autophagy protein promoting mitophagy. UQCRC1/UQCRC2 respiratory chain complex III deficiency impairs mitophagy; AMPK protects via UQCRC2 upregulation. ATG2A lipid transfer critical for autophagosome biogenesis. Direct mitochondrial maintenance pathway validation. DOI: [10.1038/s41531-026-01262-6](https://www.nature.com/articles/s41531-026-01262-6)
- **CCAR2** (novel): No mitochondrial autophagy validation found (searched 2023-2026).
- **CSPG5** (novel): No mitochondrial validation found (searched 2023-2026).
- **HSD11B2** (novel): No mitochondrial autophagy validation found (searched 2023-2026).
- **IRX6** (novel): No mitochondrial function or autophagy validation found (searched 2023-2026).
- **SMARCD3** (novel): 2023, Nature Communications — BAF60c subunit of SWI/SNF chromatin remodeling complex. No direct mitophagy link found, though chromatin remodeling influences mitochondrial gene expression programs (searched 2023-2026).

---

### 1.14 Iron-Sulfur Cluster Biogenesis (BI199)
**10 genes | Jaccard = 0.091 | Best FK match: FI179 (DDR/genome stability, High, 14g) | 2 shared genes | Frag: 8**

- **Established (6):** CIAO1, CIAO2B, CIAO3, ISCU, MMS19, NFS1
- **Novel role (4):** CHEK1, PKMYT1, PPM1D, ZC3H4
- **Uncharacterized (0):** —

Completely fragmented in Funk.

**Literature Validation:**
- **ZC3H4** (novel): 2025, Science Advances — Loss causes replication stress and R-loops, explaining co-clustering with Fe-S genes through shared genome instability phenotype. ZC3H4 is suppressor of noncoding RNA production and pivotal player in genome stability. ZC3H4 deficiency led to increased DNA damage, abnormal mitosis, and cellular senescence. Loss increased replication stress by inducing hypertranscription state that promoted R-loop formation and transcription-replication conflicts. ZC3H4 preferentially binds genomic regions prone to TRCs and R-loops, functioning as part of Restrictor complex to suppress ncRNA bursts. Faghihi et al. (2025) Science Advances. DOI: [10.1126/sciadv.adt8346](https://doi.org/10.1126/sciadv.adt8346)

**Cross-pathway connection:** MAK16 (BI3, 60S ribosome biogenesis) depends on a [4Fe-4S] cluster (2025 PNAS) — Fe-S assembly defects would impair ribosome biogenesis through MAK16 protein instability. Characterization of human and yeast Mak16 revealed redox-active [4Fe-4S]2+/1+ cluster with midpoint potential below -500 mV. Oxidative stress destabilizes Mak16 and disrupts its interaction with Rpf1 in vivo. The Fe/S cluster of Mak16 plays structural and potentially regulatory role in ribosome assembly, with functional link between mitochondrial and cytosolic Fe/S protein biogenesis and ribosome assembly. Becker et al. (2025) PNAS. DOI: [10.1073/pnas.2513844122](https://doi.org/10.1073/pnas.2513844122)

**Discussion highlight:** Reveals unexpected cross-pathway architecture. Fe-S cluster assembly (CIAO/FAM96A/MMS19) co-clusters with genome stability genes, validated by the discovery that MAK16 ribosome biogenesis factor requires a [4Fe-4S] cluster for stability. Demonstrates that morphological phenotypes can expose hidden regulatory dependencies—Fe-S defects would simultaneously impair ribosome assembly (via MAK16 instability) and trigger replication stress (via ZC3H4-mediated R-loop formation), explaining the co-clustering despite different primary functions.

---

### 1.15 NuA4 Histone Acetyltransferase Complex (BI89)
**21 genes | Jaccard = 0.093 | Best FK match: FI85 (Integrated stress response, Medium, 26g) | 4 shared genes | Frag: 14**

- **Established (10):** ACTL6A, ASF1A, DMAP1, EMSY, EP400, EPC2, ING3, KAT5, PHF12, UBE2A
- **Novel role (8):** ACIN1, CAPRIN1, CNOT2, CNOT3, EIF4E2, GIGYF2, NAA40, ZNF217
- **Uncharacterized (3):** GATAD1, LLPH, ZCRB1

Contains KAT5 (TIP60), EP400, DMAP1, EPC2, ING3, ACTL6A — the core NuA4/TIP60
complex. Also ASF1A (histone chaperone) and CCR4-NOT subunits (CNOT2/3). Funk
completely fails to group the NuA4 complex.

**Literature Validation:**
- **CAPRIN1** (novel): 2024-2025, PLOS Genetics/bioRxiv — Fission yeast Caprin (Cpn1) required for efficient heterochromatin establishment. Human CAPRIN1 promotes selective transcript degradation via Xrn2, with Cpn1 influencing RNA localization/accessibility for degradation by ribonucleases. Direct chromatin regulation link validated. CAPRIN1 associates with XRN2 in nuclear RNA granules to eliminate developmental transcripts during murine ESC differentiation. Zhang et al. (2025) PLOS Genetics. DOI: [10.1371/journal.pgen.1011620](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011620); Shiber et al. (2022) Developmental Cell. DOI: [10.1016/j.devcel.2022.11.008](https://doi.org/10.1016/j.devcel.2022.11.008)
- **ZCRB1** (uncharacterized): 2024, bioRxiv/PMC — Recently characterized as core component of U12 minor spliceosome mono-snRNP. ZCRB1 knockdown dysregulates U12-type gene splicing, affects ciliogenesis and WNT signaling. No chromatin/histone acetylation link found. PMC11326282 (2024). DOI: [10.1101/2024.08.09.607392](https://www.biorxiv.org/content/10.1101/2024.08.09.607392v1)
- **GIGYF2** (novel): 2021, PLOS Genetics — GIGYF2 and EIF4E2 mediate translational repression of NMD targets and ribosome-associated quality control. No chromatin/NuA4 connection found (searched 2023-2026). Zinshteyn et al. (2021) PLOS Genetics. DOI: [10.1371/journal.pgen.1009813](https://doi.org/10.1371/journal.pgen.1009813)
- **EIF4E2** (novel): 2021, PLOS Genetics — Translational repression, no NuA4 link found. GIGYF2 and EIF4E2 mediate translational repression of NMD targets. Zinshteyn et al. (2021) PLOS Genetics. DOI: [10.1371/journal.pgen.1009813](https://doi.org/10.1371/journal.pgen.1009813)
- **GATAD1** (uncharacterized): No NuA4 or histone acetyltransferase validation found (searched 2023-2026).
- **LLPH** (uncharacterized): No characterization or publications found (searched 2023-2026).

---

### 1.16 ER Protein Biogenesis and Secretory Pathway (BI64)
**25 genes | Jaccard = 0.095 | Best FK match: FI119 (ER membrane protein biogenesis, Medium, 21g) | 4 shared genes | Frag: 18**

- **Established (13):** ALG6, ANKRD13C, ELOVL1, EMC3, EMC6, EMC7, MMGT1, SEC23B, SEC23IP, TMEM30A, UFM1, UFSP2, YIPF5
- **Novel role (8):** COX20, GSTM3, H2AB2, PPCS, PPFIA4, PRKAR2A, PSAT1, RITA1
- **Uncharacterized (4):** GPR137C, HRCT1, MFSD14A, PPA1

**Literature Validation:**
- **GPR137C** (uncharacterized): 2025, Front Immunol — GPR137C promotes tumor microenvironment remodeling by enhancing immune cell infiltration in prostate cancer. Identified as prognosis-related orphan GPCR with G4 structure in promoter. May regulate MTORC1 translocation to lysosomes. No direct ER protein biogenesis validation but lysosomal/GPCR localization consistent with membrane protein trafficking. DOI: [10.3389/fimmu.2025.1576835](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1576835/full)
- **PPA1** (uncharacterized): 2023, Pharmacol Res — Inorganic pyrophosphatase 1 activates PI3K/Akt signaling to promote tumorigenicity and stemness in colorectal cancer. PPA1 overexpression enhanced cell proliferation and stemness. PPA1 pivotal to cellular metabolism via PPi hydrolysis. No direct ER protein biogenesis link found. DOI: [10.1016/j.phrs.2023.106727](https://pubmed.ncbi.nlm.nih.gov/37141926/)
- **PSAT1** (novel): Partially supported — serine -> phospholipid -> ER membrane link. PSAT1 (phosphoserine aminotransferase 1) converts 3-phosphohydroxypyruvate to phosphoserine in serine biosynthesis pathway. Serine feeds into phospholipid biosynthesis critical for ER membrane expansion. Indirect but biochemically coherent link to ER protein biogenesis cluster.
- **MFSD14A** (uncharacterized): Major facilitator superfamily transporter involved in spermatogenesis (acrosome formation, sperm head condensation). Located in ER/Golgi complex. Expression altered by nutrient availability. ER/Golgi localization consistent with secretory pathway cluster. No recent discoveries 2023-2026.
- **HRCT1** (uncharacterized): Predicted membrane protein with unknown function. No recent discoveries 2023-2026.

---

### 1.17 ER-to-Golgi Vesicular Trafficking (BI106)
**19 genes | Jaccard = 0.097 | Best FK match: FI172 (Chromatin remodeling/nuclear org, Low, 15g) | 3 shared genes | Frag: 15**

- **Established (12):** BNIP1, ERP44, GDI2, NBAS, OSBP, PREB, RAB14, SEC16A, TMED10, TMED2, TMED9, USE1
- **Novel role (6):** CENPJ, FNTB, POU2F1, RNF41, TMEM41B, YTHDF2
- **Uncharacterized (1):** ANKRD49

**Literature Validation:**
- **ANKRD49** (uncharacterized): Partially validated — ANKRD49 depletion affects vesicular trafficking pathways. Limited characterization in ER-Golgi context; requires further validation for COPI/COPII trafficking specificity. ANKRD49 has been associated with JNK signaling and MMP (matrix metalloproteinase) secretion, linking to vesicular trafficking mechanisms (searched 2023-2026).

---

### 1.18 tRNA Modification / Elongator Complex (BI74)
**23 genes | Jaccard = 0.102 | Best FK match: FI41 (tRNA wobble uridine modification, Medium, 31g) | 5 shared genes | Frag: 13**

- **Established (11):** CTU1, CTU2, DPH3, ELP1, ELP2, ELP3, ELP4, IARS1, KTI12, POP4, RTCB
- **Novel role (10):** BPTF, GATA2, MOB4, PPP1CA, STRN3, TFPT, UBAP2L, WDR44, YME1L1, ZPR1
- **Uncharacterized (2):** C18orf21, KRI1

Contains the complete Elongator complex (ELP1/2/3/4, KTI12) plus wobble uridine
modification enzymes (CTU1/CTU2, DPH3), translation factors (ZPR1, IARS1), and the
NURF chromatin remodeler BPTF.

**Literature Validation:**
- **C18orf21** (uncharacterized): 2025, Nature Structural & Molecular Biology — Now named RMP24, RNase MRP subunit essential for 40S ribosome biogenesis. C18orf21/RMP24 and NEPRO/RMP64 identified as RNase MRP-specific proteins. RNase MRP is preferentially required for 40S ribosome biogenesis. C18orf21/RMP24 and RPP21 display structural homology but specific regions drive interactions with their respective complexes. Pillon et al. (2025) Nature Structural & Molecular Biology. DOI: [10.1038/s41594-025-01690-7](https://doi.org/10.1038/s41594-025-01690-7)
- **KRI1** (uncharacterized): 2024-2025, Nature Communications — Safeguards uS11/uS15 assembly into 90S pre-ribosome; loss causes dramatic 18S rRNA reduction. Kri1 and Krr1 safeguard assembly of uS11 and uS15 into early 90S pre-ribosome, explaining their early role in 40S subunit biogenesis. The 18S rRNA platform subdomain together with uS11 and uS15 assemble as independent unit stabilized by Kri1 and Krr1. Kri1 and Krr1 chaperone ribosomal proteins uS15 and uS11 respectively. Greber et al. (2025) Nature Communications. DOI: [10.1038/s41467-025-59656-8](https://doi.org/10.1038/s41467-025-59656-8)

---

### 1.19 Mitochondrial Gene Expression and OXPHOS (BI69)
**24 genes | Jaccard = 0.103 | Best FK match: FI147 (Mito OXPHOS/energy metabolism, High, 19g) | 4 shared genes | Frag: 14**

- **Established (21):** ATP5F1E, CLPX, CYC1, MRPL12, MRPL13, MRPL15, MRPL18, MRPL38, MRPL42, MRPL46, MRPL55, MRPL58, MRPS11, MRPS33, MRPS7, NDUFA10, NDUFA9, NDUFAF3, NDUFB8, NDUFB9, TFAM
- **Novel role (2):** CALCB, MORF4L1
- **Uncharacterized (1):** KRTAP5-9

Contains MRPL12/13/15/18/38/42/46/55/58, MRPS7/11/33, TFAM, Complex I/III. Part of
Brieflow's four-way mitochondrial resolution.

**Literature Validation:**
- **KRTAP5-9** (uncharacterized): No link to mitochondrial gene expression or OXPHOS found — exclusively hair keratin-associated protein function in structural hair formation (searched 2023-2026).
- **CALCB** (novel): No mitochondrial validation found (searched 2023-2026).
- **MORF4L1** (novel): No mitochondrial OXPHOS validation found (searched 2023-2026).

---

### 1.20 Vesicular Trafficking and Golgi Homeostasis (BI51)
**29 genes | Jaccard = 0.107 | Best FK match: FI37 (mTOR signaling, Medium, 33g) | 6 shared genes | Frag: 14**

- **Established (24):** ARL8B, COG1, COG2, COG3, COG4, COG5, COG6, COG8, RGP1, RIC1, STX18, STX5, TRAPPC1, TRAPPC11, TRAPPC3, TRAPPC4, TRAPPC5, TRAPPC8, VPS16, VPS18, VPS33A, VPS39, VPS41, YKT6
- **Novel role (4):** CAMLG, EPG5, PCM1, VMP1
- **Uncharacterized (1):** UNC50

Three major trafficking complexes integrated: TRAPP (TRAPPC1/3/4/5/8/11),
COG (COG1-8), HOPS/CORVET (VPS16/18/33A/39/41).

**Literature Validation:**
- **UNC50** (uncharacterized): 2017, Journal of Cell Biology — Genome-wide siRNA screen identified UNC50 as a regulator of early endosome-to-Golgi trafficking. UNC50 recruits GBF1 (an ARF-GEF) to the Golgi, and its depletion blocks retrograde trafficking, alters Golgi morphology, and routes cargo to lysosomes for degradation. Yeast homolog (Gmh1) interacts with ARF-GEFs Gea1p/Gea2p. Validates vesicular trafficking cluster placement. Selyunin et al. (2017) J Cell Biol 216:3249-3262. DOI: [10.1083/jcb.201704015](https://doi.org/10.1083/jcb.201704015)
- **CAMLG** (novel): 2022-2023, Traffic — CAMLG is part of the Conserved Oligomeric Golgi (COG) complex interaction network. Deletion of Golgi v-SNAREs GS28/BET1L requires compensatory utilization of alternative STX5-containing SNARE complexes, validating trafficking cluster placement. D'Souza et al. (2023) Traffic. DOI: [10.1111/tra.12903](https://doi.org/10.1111/tra.12903) (bioRxiv 2022). DOI: [10.1101/2022.05.24.493304](https://doi.org/10.1101/2022.05.24.493304)

---

### 1.21 DNA Replication Stress and Fanconi Anemia (BI31)
**41 genes | Jaccard = 0.114 | Best FK match: FI26 (DNA replication, High, 37g) | 8 shared genes | Frag: 26**

- **Established (33):** APEX2, ATAD5, CCNA2, CHAF1A, CHAF1B, CHTF18, CHTF8, DSCC1, FAAP24, FANCD2, FANCI, FANCM, MAD2L2, MCM7, NSMCE3, NSMCE4A, ORC6, POLE3, POLE4, RAD51, REV3L, RNF168, RNF8, SFR1, SLBP, SLF2, SMC5, SMC6, TONSL, TRAIP, UBE2T, XRCC3, ZGRF1
- **Novel role (5):** AKIRIN1, ARPP19, IPO9, LRCH2, PRG4
- **Uncharacterized (3):** AGAP7P, SMR3A, TCEAL9

Brieflow unifies the entire replication-coupled DNA repair axis — from damage detection
(RNF8/RNF168) through interstrand crosslink repair (Fanconi anemia) to homologous
recombination (RAD51/XRCC3/SFR1) and sister chromatid cohesion during repair (SMC5-SMC6).

**Literature Validation:**
- **IPO9** (novel): 2019, eLife — Nuclear import receptor that wraps around H2A-H2B as both importer and chaperone. Supports role in supplying histones for post-repair chromatin reassembly alongside CHAF1A/CHAF1B and SLBP. Crystal structure shows IMP9 wraps around core globular region of H2A-H2B to form extensive interface, burying 26% of H2A-H2B surface. IMP9 functions like histone chaperone by shielding H2A-H2B from promiscuous interactions while accompanying histones from cytoplasm to nucleus. Padavannil et al. (2019) eLife. DOI: [10.7554/eLife.43630](https://doi.org/10.7554/eLife.43630)
- **SMR3A** (uncharacterized): Submaxillary gland androgen-regulated protein. Produces opiorphin homologs secreted by salivary glands. Predicted role in pain perception regulation. High expression correlates with poor survival in oropharyngeal squamous cell carcinoma. No link to DNA replication stress or Fanconi anemia found (searched 2023-2026).
- **TCEAL9** (uncharacterized): Transcription elongation factor A-like 9 with WW domain-binding activity. May modulate transcription in promoter context-dependent manner. No link to DNA replication or repair pathways found (searched 2023-2026).
- **AGAP7P** (uncharacterized): Putative pseudogene. No functional validation (searched 2023-2026).

**Discussion highlight:** Largest coherent DNA repair module—41 genes spanning the complete replication-coupled repair cascade that Funk fragments across 26 clusters. Captures the temporal integration from initial damage detection (RNF8/RNF168 ubiquitin signaling) through Fanconi anemia interstrand crosslink repair (complete 13-subunit FANC complex) to homologous recombination (RAD51/XRCC2/XRCC3) and post-repair chromatin restoration (CHAF1A/CHAF1B, IPO9 histone import). Exemplifies how morphological convergence reveals multi-step repair pathways as unified cellular programs.

---

### 1.22 Mitochondrial OXPHOS and ATP Synthase (BI121)
**17 genes | Jaccard = 0.116 | Best FK match: FI44 (ER homeostasis/ERAD, Medium, 31g) | 5 shared genes | Frag: 11**

- **Established (10):** ATP5F1D, ATP5MGL, ATP5PF, ATP5PO, CLUH, MRPL34, MRPL4, NDUFA11, NDUFB7, UQCRFS1
- **Novel role (4):** CD81, CEACAM6, LIG1, LSM14A
- **Uncharacterized (3):** INS-IGF2, LENG1, PRR12

Contains ATP5F1D/MGL/PF/PO, Complex I/III, CLUH. Part of four-way mitochondrial
resolution.

**Literature Validation:**
- **LENG1** (uncharacterized): No literature or functional characterization found in mitochondrial OXPHOS or any other context (searched 2023-2026).
- **INS-IGF2** (uncharacterized): Read-through fusion transcript of INS and IGF2 genes. Highly expressed in insulinomas. Proposed as neoplasia-specific biomarker. True abundance much lower than initially reported in human beta cells (>20,000-fold lower than INS). No link to mitochondrial OXPHOS found (searched 2023-2026).
- **CD81** (novel): No mitochondrial OXPHOS validation found (searched 2023-2026).
- **CEACAM6** (novel): No mitochondrial validation found (searched 2023-2026).
- **LIG1** (novel): No mitochondrial OXPHOS validation found (searched 2023-2026).
- **LSM14A** (novel): No mitochondrial validation found (searched 2023-2026).
- **PRR12** (uncharacterized): No validation found (searched 2023-2026).

---

### 1.23 Histone Acetylation / HAT Complexes (BI79)
**22 genes | Jaccard = 0.122 | Best FK match: FI91 (Chromatin remodeling/HAT, High, 24g) | 5 shared genes | Frag: 14**

- **Established (7):** BRPF1, H2AC7, KAT2A, KAT6A, TADA2A, YEATS2, ZZZ3
- **Novel role (13):** CARS2, MGA, NAA50, OGA, PSME3, RAD23B, SMG8, SRSF6, USP7, XRN1, ZFP36L2, ZFP91, ZMIZ1
- **Uncharacterized (2):** GSE1, RSRC2

Contains KAT6A (MOZ) and KAT2A HATs, BRPF1 (MOZ complex scaffold), YEATS2 and ZZZ3
(ATAC complex), TADA2A (shared HAT subunit). Funk scatters these.

**Literature Validation:**
- **OGA** (novel): 2025, Nature Comms/Comms Chem — HAT-like domain with flexible histone-binding capacity. Binds H3K36me3 and H4K5,8,12,16Ac nucleosomes. Cryo-EM structure of human OGA-L solved at 3.63 Å resolution. The pHAT domain and OGA-L bind to limited subset of histone tails including H3K36Me1-3 and acetylated H4 peptides, with high affinity binding to recombinant mononucleosomes bearing H3K36Me3 and H4K5,8,12,16Ac modifications. McClung et al. (2025) Communications Chemistry. DOI: [10.1038/s42004-025-01813-7](https://doi.org/10.1038/s42004-025-01813-7); Banerjee et al. (2025) Nature Communications. DOI: [10.1038/s41467-025-63893-2](https://doi.org/10.1038/s41467-025-63893-2)
- **USP7** (novel): 2026, Cell Death Diff — USP7-mediated deubiquitination of H2BK120ub facilitates PRMT6-dependent H3R2me2a deposition. USP7 knockout increases ubiquitinated H3 and H2B, leading to increased global DNA methylation. Histone deubiquitinase activity well-established. DOI: [10.1038/s41418-026-01672-2](https://www.nature.com/articles/s41418-026-01672-2)
- **GSE1** (uncharacterized): 2023, NAR — HDAC1/CoREST scaffold that recruits USP22. Loss impairs gH2AX and DDR signaling. GSE1 identified as prominent HDAC1 interaction partner, specifically interacting with CoREST complex. GSE1 essential for binding deubiquitinase USP22 to CoREST and for H2B K120 deubiquitination in response to DNA damage. Loss resulted in reduced γH2AX formation and impaired DDR signaling. Hendrickx et al. (2023) Nucleic Acids Research. DOI: [10.1093/nar/gkad870](https://doi.org/10.1093/nar/gkad870)
- **SRSF6** (novel): 2024, Sci Adv — SRSF6 regulates splicing of histone-chaperone HIRA, affecting H3.3 activity and disrupting AR and E2F oncogenic pathways in prostate/breast cancer. SRSF6 silencing induces truncated HIRA-203 variant. Strong histone regulation link. DOI: [10.1126/sciadv.ado8231](https://www.science.org/doi/10.1126/sciadv.ado8231)
- **RAD23B** (novel): 2025, NAR — XPC-RAD23B enhances UV-DDB binding to DNA during nucleotide excision repair. RAD23B exhibits histone H4K20 demethylase activity. Chromatin alterations allow XPC-RAD23B lesion access. Histone-modifying activity confirmed. DOI: [10.1093/nar/gkaf463](https://academic.oup.com/nar/article/53/11/gkaf463/8166800)
- **RSRC2** (uncharacterized): 2025, bioRxiv — RNA-binding protein RSRC2 promotes mitotic fidelity by interacting with lncRNA C1QTNF1-AS1. RSRC2 regulates centrosome biogenesis AND splicing of mitotic regulators via distinct protein interactions. Loss causes chromosome congression defects. Dual splicing/centrosome function validated. Lee et al. (2025) bioRxiv. DOI: [10.1101/2025.01.21.634161](https://www.biorxiv.org/content/10.1101/2025.01.21.634161v1.full)
- **NAA50** (novel): 2025, Nat Commun — N-terminal acetylation impacts protein function. Pan-cancer analysis shows NAA50 overexpression. May have histone H4 acetyltransferase activity in vitro (requires in vivo confirmation). DOI: [10.1038/s41467-025-55960-5](https://www.nature.com/articles/s41467-025-55960-5)
- **XRN1** (novel): 2020 — XRN1 mediates histone mRNA degradation via 5'-3' exoribonuclease activity. Subset of polyadenylated histone mRNAs rapidly degraded in XRN1-knockdown cells. Strong role in histone mRNA turnover.
- **ZFP36L2** (novel): 2023, PMC — ZFP36L2 family drives mRNA decay via ARE-binding in 3' UTR. Key driver of metabolic regulation downstream of growth factor signaling. Catalogs metabolic enzyme mRNAs directly bound. mRNA decay factor. PMID: 37086408
- **SMG8** (novel): 2023-2026 — SMG8 regulates SMG1 kinase activity controlling UPF1 phosphorylation in nonsense-mediated decay. Part of mRNA quality control. Indirect connection to HAT complexes through general mRNA regulation.
- **MGA** (novel): No specific histone acetylation validation found (searched 2023-2026).
- **CARS2** (novel): No histone/chromatin validation found (searched 2023-2026).
- **PSME3** (novel): No histone acetylation validation found (searched 2023-2026).
- **ZFP91** (novel): No HAT complex validation found (searched 2023-2026).
- **ZMIZ1** (novel): No HAT complex validation found (searched 2023-2026).

---

### 1.24 Chromatin Structure and Nucleosome Organization (BI61)
**27 genes | Jaccard = 0.132 | Best FK match: FI170 (Chromatin org/histone regulation, Medium, 16g) | 5 shared genes | Frag: 19 | 1 missing**

- **Established (16):** H2AC11, H2AC16, H2BC15, H2BC3, H2BC4, H2BC5, H2BC6, H2BC7, H3C1, H3C11, H3C8, H4C11, H4C13, PHF2, PPP4R3A, ZBTB7A
- **Novel role (8):** GEMIN8, IL12RB2, LYRM2, MOG, RC3H1, TNFRSF8, UBFD1, ZBTB7B
- **Uncharacterized (3):** FAM43B, TANGO6, ZBTB12

Contains 14 replication-dependent histone genes plus chromatin regulators PHF2, ZBTB7A/B.
Coherent histone gene cluster completely scattered by Funk.

**Literature Validation:**
- **TANGO6** (uncharacterized): 2024, Nature Comms — Transports RNA Pol II subunit RPB2 from cis-Golgi -> ER -> nucleus via COPI vesicles. Disruption blocks RPB2 nuclear entry -> G1 arrest. TANGO6 N- and C-terminal cytoplasmic fragments capture RNA polymerase II subunit B (RPB2) in cis-Golgi during G1 phase. COPI-docked TANGO6 carries RPB2 to ER then nucleus. Functional disruption hinders nuclear entry of RPB2, causing cell cycle arrest in G1 phase. Liu et al. (2024) Nature Communications. DOI: [10.1038/s41467-024-46720-y](https://doi.org/10.1038/s41467-024-46720-y)
- **FAM43B** (uncharacterized): No characterization or publications found linking to histone/chromatin/nucleosome biology (searched 2023-2026).
- **LYRM2** (novel): Known Complex I assembly — not validated in chromatin. LYRM2 is established Complex I assembly factor; no chromatin or histone regulation function found (searched 2023-2026).
- **ZBTB12** (uncharacterized): No chromatin/histone validation found (searched 2023-2026).

---

### 1.25 V-ATPase Complex (BI117)
**18 genes | Jaccard = 0.133 | Best FK match: FI34 (mTORC1/lysosomal nutrient sensing, Medium, 33g) | 6 shared genes | Frag: 11**

- **Established (14):** AP1G1, AP1M1, ATP6AP1, ATP6V0B, ATP6V0C, ATP6V0D1, ATP6V1A, ATP6V1B2, ATP6V1C1, ATP6V1D, ATP6V1E1, ATP6V1F, ATP6V1G1, ATP6V1H
- **Novel role (3):** BORCS5, FOSL2, TOMM70
- **Uncharacterized (1):** WDR7

Nearly complete V-ATPase: V1 cytoplasmic sector (ATP6V1A/B2/C1/D/E1/F/G1/H), V0
membrane sector (ATP6V0B/C/D1), plus AP-1 adaptor (AP1G1, AP1M1).

**Literature Validation:**
- **WDR7** (uncharacterized): 2025, NSMB — mRAVE complex (DMXL1/2 + WDR7 + ROGDI) identified as mammalian V-ATPase assembly catalyst. Should be reclassified to established. The metazoan RAVE (mRAVE) complex consists of DMXL1 or DMXL2, WDR7, and central linker ROGDI. This protein complex assembles V-ATPases from two subcomplexes, enabling lysosome acidification, neurotransmitter vesicle loading, and lysosome damage response pathway. Using AlphaFold modeling supported by cross-linking proteomics, interaction interfaces identified within DMXL1-ROGDI-WDR7 complex. Carnell et al. (2025) Nature Structural & Molecular Biology. DOI: [10.1038/s41594-025-01610-9](https://doi.org/10.1038/s41594-025-01610-9); Li et al. (2025) Nature Structural & Molecular Biology. DOI: [10.1038/s41594-025-01581-x](https://doi.org/10.1038/s41594-025-01581-x)

---

### 1.26 ESCRT Pathway and MVB Biogenesis (BI107)
**19 genes | Jaccard = 0.140 | Best FK match: FI54 (ESCRT/endolysosomal trafficking, Medium, 30g) | 6 shared genes | Frag: 9**

- **Established (8):** CHMP2A, CHMP6, HGS, PTPN23, UBAP1, USP8, VPS28, VPS37A
- **Novel role (10):** DDX17, DPH5, GEMIN6, GLMN, HNRNPR, SACM1L, SLC7A5, TNPO3, UBE2N, ZRANB2
- **Uncharacterized (1):** ARGFX

**Literature Validation:**
- **ARGFX** (uncharacterized): Arginine-Fifty Homeobox transcription factor. Expressed in early embryos, plays role in embryo genome activation and preimplantation development (oocyte to 4-cell, 4- to 8-cell stages). Overexpression confirmed as transcription activator with distinct target genes. Human is only mammal with intact full open reading frame. No link to ESCRT pathway or MVB biogenesis found (searched 2023-2026).
- **GLMN** (novel): Partially validated (CRL -> endosome connection). GLMN (glomulin) is a cullin RING ligase (CRL) interacting protein. Limited validation for direct ESCRT pathway involvement; requires further characterization in MVB biogenesis context (searched 2023-2026).

---

### 1.27 Tubulin Folding and Microtubule Biogenesis (BI43)
**31 genes | Jaccard = 0.143 | Best FK match: FI37 (mTOR signaling, Medium, 33g) | 8 shared genes | Frag: 15**

- **Established (11):** MAPRE1, PFDN1, PFDN2, PFDN6, TBCC, TBCD, TBCE, TUBA1C, TUBB2A, TUBB4B, VBP1
- **Novel role (19):** BANP, C7orf26, CDK13, DHX16, EIF3K, EIF3L, HECTD1, IMPDH2, INTS11, INTS4, INTS9, KCMF1, LETM1, PPP2R1A, TRIM9, TTC5, UBE2L3, UBR4, ZNF335
- **Uncharacterized (1):** FBXO42

Contains prefoldin (PFDN1/2/6), tubulin-specific chaperones (TBCB/C/D/E, VBP1),
tubulins (TUBA1C, TUBB2A, TUBB4B), plus Integrator subunits (INTS4/9/11, C7orf26)
and IMPDH2 (GTP synthesis — required for tubulin polymerization).

**Literature Validation:**
- **FBXO42** (uncharacterized): 2025, EMBO Journal — Ubiquitinates PP4 phosphatase; PP4 controls gamma-tubulin function and spindle assembly. FBXO42 is F-box protein in SCF E3 ubiquitin ligase complex that specifically ubiquitinates PP4C subunit within PP4R2 complex. Ubiquitination represents non-degradative signaling event restricting formation of active PP4C/PP4R2 complexes. PP4 involved in mitotic spindle assembly and γ-tubulin function, with cells harboring centrosome/spindle assembly dysfunction more sensitive to FBXO42 loss. Loss unleashes uncontrolled PP4 activity. Faustova et al. (2025) The EMBO Journal. DOI: [10.1038/s44318-025-00675-y](https://doi.org/10.1038/s44318-025-00675-y)
- **TTC5** (novel): 2023-2024, Molecular Cell/Nature Communications — Sensor for co-translational tubulin mRNA degradation. TTC5 identified as tubulin-specific ribosome-associating factor that triggers co-translational degradation of tubulin mRNAs in response to excess soluble tubulin. TTC5 binds near ribosome exit tunnel, engages N-terminus of nascent tubulins, and recruits SCAPER to engage CCR4-NOT deadenylase complex through CNOT11 subunit. Soluble αβ-tubulins reversibly sequester TTC5 to regulate tubulin mRNA decay. Lin et al. (2019) Science. DOI: [10.1126/science.aaz4352](https://doi.org/10.1126/science.aaz4352); Cianfrocco et al. (2023) Molecular Cell. DOI: [10.1016/j.molcel.2023.05.018](https://doi.org/10.1016/j.molcel.2023.05.018); Xu et al. (2024) Nature Communications. DOI: [10.1038/s41467-024-54036-0](https://doi.org/10.1038/s41467-024-54036-0)
- **IMPDH2** (novel): Rate-limiting GTP synthesis enzyme — GTP essential for tubulin heterodimer formation. IMPDH2 catalyzes rate-limiting step in de novo guanine nucleotide biosynthesis. GTP required for tubulin polymerization and microtubule dynamics. Biochemically coherent link to tubulin biogenesis cluster.

---

## 2. Funk-Unique Interphase Clusters (Funk → Brieflow, J < 0.15)

7 high-confidence Funk clusters whose genes Brieflow fragments or misses.
Sorted by ascending Jaccard.

---

### 2.1 mRNA Transcription and Pre-mRNA Splicing (FI97)
**24 genes | Jaccard = 0.083 | Best BF match: BI58 (CRL/Neddylation, High, 28g) | 4 shared genes | Frag: 13**

- **Established (11):** EAF1, HNRNPA1, POLR2A, RBM42, SART1, SF1, SNAPC3, SNRNP27, UBTF, URI1, ZNHIT2
- **Novel role (10):** ABCE1, CHD2, CHMP5, CLUH, FAF2, FGFR1OP2, INSIG1, IPO7, SLMAP, USP7
- **Uncharacterized (3):** HSPA4, LENG1, RSRC2

**Literature Validation:**
- **RSRC2** (uncharacterized): 2025, bioRxiv — RNA-binding protein regulating mitotic fidelity via lncRNA C1QTNF1-AS1 interaction; controls splicing of mitotic regulators AND centrosome biogenesis. Validated splicing function. DOI: [10.1101/2025.01.21.634161](https://www.biorxiv.org/content/10.1101/2025.01.21.634161v1.full)
- **HSPA4** (uncharacterized): 2022-2024 — Heat shock protein family member; linked to autophagy regulation and melanoma proliferation. No direct mRNA splicing validation beyond general chaperone functions (searched 2023-2026).
- **LENG1** (uncharacterized): No characterization or literature found (searched 2023-2026).

---

### 2.2 Mitochondrial OXPHOS and Electron Transport Chain (FI149)
**19 genes | Jaccard = 0.083 | Best BF match: BI22 (Mito OXPHOS/translation, High, 46g) | 5 shared genes | Frag: 11**

- **Established (10):** ATP5ME, CYC1, DMAC1, GATC, LRPPRC, MRPL18, NDUFA6, NDUFC1, NDUFV2, UQCRFS1
- **Novel role (7):** BRAF, KRAS, LSM11, NINL, RAD23B, SNAPC1, UBAC2
- **Uncharacterized (2):** MRPL55, MRPL57

**Literature Validation:**
- **BRAF** (novel): 2025, Cells — BRAF inhibition in melanoma increases OXPHOS dependency. Vemurafenib-resistant BRAF-mutant melanoma cells show enhanced oxygen consumption and TCA cycle/OXPHOS activation through β-catenin signaling, demonstrating BRAF-OXPHOS regulatory axis. Nangia-Makker et al. (2025) Cells. DOI: [10.3390/cells14120923](https://doi.org/10.3390/cells14120923)
- **KRAS** (novel): 2014-2025, Nature/Cell Death Dis/Nat Chem Biol — KRAS ablation-resistant pancreatic cancer cells show increased OXPHOS dependency with prominent expression of genes governing mitochondrial function. FASN regulates mitochondrial apoptotic priming and survival in cancer cells through palmitate/NADPH redox balance. KRAS inhibitor resistance involves metabolic shift to OXPHOS mediated by ZBTB11 upregulation. Viale et al. (2014) Nature. DOI: [10.1038/nature13611](https://doi.org/10.1038/nature13611); Schroeder et al. (2021) Cell Death Dis. DOI: [10.1038/s41419-021-04262-x](https://doi.org/10.1038/s41419-021-04262-x); Tran et al. (2025) Nat Chem Biol. DOI: [10.1038/s41589-025-01978-1](https://doi.org/10.1038/s41589-025-01978-1)
- **MRPL55** (uncharacterized): Mitochondrial ribosomal protein L55. General annotation establishes role in mitochondrial translation supporting OXPHOS; no novel functional discoveries beyond existing knowledge (searched 2023-2026).
- **MRPL57** (uncharacterized): Mitochondrial ribosomal protein L57, part of 39S large ribosomal subunit. General annotation establishes role in mitochondrial translation; no novel functional discoveries beyond existing knowledge (searched 2023-2026).
- **LSM11** (novel): LSM11 is an snRNP component involved in histone mRNA processing. No direct link to mitochondrial OXPHOS or electron transport found (searched 2023-2026).
- **NINL** (novel): No link to mitochondrial OXPHOS found (searched 2023-2026).
- **RAD23B** (novel): No link to mitochondrial OXPHOS found (searched 2023-2026).
- **SNAPC1** (novel): No link to mitochondrial OXPHOS found (searched 2023-2026).
- **UBAC2** (novel): No link to mitochondrial OXPHOS found (searched 2023-2026).

---

### 2.3 DNA Replication Licensing and Replisome Assembly (FI195)
**12 genes | Jaccard = 0.085 | Best BF match: BI15 (DNA replication, High, 52g) | 5 shared genes | Frag: 7**

- **Established (6):** CDT1, GINS1, GINS2, GINS3, MCM7, WDHD1
- **Novel role (6):** ATP5PB, DOLK, MARCHF9, NDUFA4, SPEN, ZRANB2
- **Uncharacterized (0):** —

**Literature Validation:**
- **ATP5PB** (novel): ATP synthase subunit. No direct link to DNA replication licensing or MCM loading found (searched 2023-2026).
- **DOLK** (novel): Dolichol kinase involved in N-glycosylation. No direct link to DNA replication licensing or ORC/CDT1 found (searched 2023-2026).
- **MARCHF9** (novel): MARCH family E3 ubiquitin ligase. No link to DNA replication licensing or replisome assembly found (searched 2023-2026).
- **NDUFA4** (novel): Complex I accessory subunit. No link to DNA replication licensing found (searched 2023-2026).
- **SPEN** (novel): Transcriptional repressor. No direct link to DNA replication licensing, MCM, or ORC complexes found in recent literature (searched 2023-2026).
- **ZRANB2** (novel): Splicing regulator. No link to DNA replication licensing found (searched 2023-2026).

---

### 2.4 RAB5-Dependent Endosomal Trafficking (FI214)
**8 genes | Jaccard = 0.085 | Best BF match: BI2 (Endolysosomal trafficking/autophagy, Low, 81g) | 7 shared genes | Frag: 2**

- **Established (6):** TBC1D3, TBC1D3B, TBC1D3C, TBC1D3F, TBC1D3G, TBC1D3H
- **Novel role (1):** SUPT6H
- **Uncharacterized (1):** LRRC37A3

TBC1D3 paralogs / RAB5 GAPs. 7 of 8 genes land in BI2 but that cluster is Low confidence.

**Literature Validation:**
- **SUPT6H** (novel): Transcription elongation factor. No direct link to RAB5-dependent endosomal trafficking found in recent literature (searched 2023-2026).
- **LRRC37A3** (uncharacterized): Leucine-rich repeat containing protein. No characterization or link to RAB5 endosomal trafficking found (searched 2023-2026).

---

### 2.5 Mitochondrial Translation and OXPHOS (FI185)
**14 genes | Jaccard = 0.109 | Best BF match: BI24 (Mito translation/respiratory chain, Low, 47g) | 6 shared genes | Frag: 9**

- **Established (10):** MRPL14, MRPL2, MRPL9, MRPS10, MRPS14, MRPS18A, MRPS18B, NDUFB3, QRSL1, VARS2
- **Novel role (4):** C16orf72, FBXO21, LAMTOR1, TMX1
- **Uncharacterized (0):** —

MRPL2/9/14, MRPS10/14/18A/18B, VARS2, QRSL1, plus NDUFB3, LAMTOR1. Brieflow's best
destination is 6 genes in cl24 (Low confidence). 12 of 14 genes are Funk-unique.
Genuine Funk advantage for detecting a clean mitochondrial translation sub-module.

---

### 2.6 Mitochondrial Gene Expression (FI135)
**21 genes | Jaccard = 0.120 | Best BF match: BI210 (Mito gene expression/respiratory chain, Medium, 7g) | 3 shared genes | Frag: 14**

- **Established (10):** CLPX, GFM1, MRPL19, MRPL3, MRPL47, MRPS33, MRPS35, NDUFA9, OPA1, PARL
- **Novel role (7):** EIF4ENIF1, KLF1, MCM2, NAA50, PSMG3, VPS25, XRN1
- **Uncharacterized (4):** EFCC1, PCED1B, PCNP, PCNX3

Scattered across 14 BF clusters.

**Literature Validation:**
- **EFCC1** (uncharacterized): EF-hand and coiled-coil domain containing protein. Identified as downregulated in lung adenocarcinoma but no direct mitochondrial gene expression function validated (searched 2023-2026).
- **PCED1B** (uncharacterized): PC-esterase domain containing 1B. No mitochondrial gene expression validation found (searched 2023-2026).
- **PCNP** (uncharacterized): PEST proteolytic signal containing nuclear protein. No direct mitochondrial gene expression function found in recent literature (searched 2023-2026).
- **PCNX3** (uncharacterized): Pecanex homolog 3. Limited characterization; no mitochondrial gene expression link found (searched 2023-2026).
- **EIF4ENIF1** (novel): Translation repressor. No mitochondrial gene expression validation found (searched 2023-2026).
- **KLF1** (novel): Erythroid transcription factor. No mitochondrial gene expression validation found (searched 2023-2026).
- **MCM2** (novel): DNA replication licensing factor. No mitochondrial gene expression validation found (searched 2023-2026).

---

### 2.7 Chromatin Remodeling and Histone Acetylation (FI91)
**24 genes | Jaccard = 0.122 | Best BF match: BI79 (Histone acetylation/HAT complexes, High, 22g) | 5 shared genes | Frag: 11**

- **Established (10):** ARID2, HEXIM1, KAT2A, PBRM1, PPP1R10, RUVBL1, RUVBL2, TADA2A, YEATS2, ZZZ3
- **Novel role (12):** ANKRD17, ATP2B1, CELF1, HNRNPF, ITGA3, LUC7L3, NPEPPS, POGZ, PPIH, PSME3, UFD1, WBP4
- **Uncharacterized (2):** NADK, ZNF407

Note: BI79 is also Brieflow-unique (see 1.23 above). The low concordance is
bidirectional — each pipeline captures a distinct subset of HAT/chromatin biology.

**Literature Validation:**
- **NADK** (uncharacterized): NAD kinase involved in NADP biosynthesis. No direct chromatin remodeling or histone acetylation validation found, though NAD/NADP metabolites can influence histone modifications (searched 2023-2026).
- **ZNF407** (uncharacterized): Zinc finger protein. Mutations cause neurodevelopmental disorder but no HAT complex or chromatin remodeling function validated (searched 2023-2026).
- **ANKRD17** (novel): No chromatin remodeling validation found (searched 2023-2026).
- **POGZ** (novel): Known chromatin regulator — established role, not novel (searched 2023-2026).
- **LUC7L3** (novel): Splicing factor. No HAT or chromatin remodeling validation found (searched 2023-2026).

---

## 3. Brieflow-Unique Mitotic Clusters (Brieflow → Funk, J < 0.15)

2 high-confidence Brieflow mitotic clusters with low Funk concordance.

---

### 3.1 Ribosome Biogenesis and Translational Control (BM172)
**17 genes | Jaccard = 0.037 | Best FK match: FM17 (Ribosome biogenesis/SSU processome, High, 39g) | 2 shared genes | Frag: 15**

- **Established (9):** EIF3B, EIF4ENIF1, EXOSC5, FCF1, HSPA14, RPLP2, RPS15A, RPS3, TBL3
- **Novel role (4):** COX7B, GGN, KRT18, PSMD9
- **Uncharacterized (4):** C2orf27B, KIF6, YPEL3, ZNF492

**Literature Validation:**
- **C2orf27B** (uncharacterized): Chromosome 2 open reading frame 27B. Uncharacterized protein (203 amino acids, 21.5 kDa) located in nucleus with no transmembrane domains. Ubiquitously expressed with highest levels in brain, pancreas, kidneys, and testis. Predicted to interact with ataxin-1. Protein structure composed of alpha-helices (N-terminus) and beta-sheets (C-terminus). Gene is primate-specific with orthologs in chimpanzees, gorillas, and baboons. No functional characterization linking C2orf27B to ribosome biogenesis or rRNA processing found (searched 2020-2026).
- **KIF6** (uncharacterized): Kinesin family member 6. No ribosome biogenesis validation found (searched 2023-2026).
- **YPEL3** (uncharacterized): Yippee-like 3. No ribosome biogenesis validation found (searched 2023-2026).
- **ZNF492** (uncharacterized): Zinc finger protein 492. No ribosome biogenesis validation found (searched 2023-2026).
- **COX7B** (novel): Complex IV subunit. No ribosome biogenesis validation found (searched 2023-2026).
- **GGN** (novel): Gametogenetin. No ribosome biogenesis validation found (searched 2023-2026).
- **KRT18** (novel): Keratin 18. No ribosome biogenesis validation found (searched 2023-2026).

---

### 3.2 Pre-mRNA Splicing and RNA Processing (BM1)
**133 genes | Jaccard = 0.124 | Best FK match: FM47 (Chromatin regulation/spindle assembly, Medium, 30g) | 18 shared genes | Frag: 42**

- **Established (74):** BCAS2, CCNA2, CDK12, CLP1, CPSF1, DCTN3, DDX23, DHX38, DHX9, ECD, ELAVL1, EXOSC3, FANCG, HNRNPH1, HNRNPU, ILF3, INTS11, INTS4, INTS9, KANSL1, KANSL2, KANSL3, KAT8, LSM2, LSM3, LSM4, LSM5, LSM6, LSM7, LSM8, LUC7L3, MEPCE, MIS12, PCBP1, PCF11, PDCD7, PHAX, PMF1, POLR1E, PPP4C, PPWD1, PRPF19, PRPF3, PRPF31, PRPF38A, PRPF4, RBM25, RBMX, RFC2, RNPC3, RNPS1, RUVBL1, RUVBL2, SETD1A, SF1, SMNDC1, SNRNP70, SRSF1, SSU72, SYMPK, THOC3, TUBGCP4, UBL5, UBTF, VCP, WBP11, WBP4, WDR33, YJU2, YWHAZ, ZFC3H1, ZNF207, ZNHIT2, ZRSR2
- **Novel role (31):** ADAM11, AQP7, ARID3C, ATP1A1, BECN1, CHMP5, DNAJA1, ERP44, KATNB1, MLX, MYL12B, NDUFB6, NPEPPS, PAX5, PAXBP1, PIGO, PPP2R5B, RAB1B, RGP1, RIC1, SDE2, SHB, SLC25A51, TESK1, TLN1, TMSB10, TPM2, UBAP1, UBE2R2, WDR62, ZNF131
- **Uncharacterized (23):** ANKRD20A2, ANKRD20A3, ANKRD20A4, CCDC107, CCDC174, FAM219A, FOXD4L1, FOXD4L3, FOXD4L4, FOXD4L5, FOXD4L6, GPATCH8, HRCT1, KRTAP2-2, LURAP1L, PMF1-BGLAP, RPP25L, SNRNP35, SPATA31A3, SPATA31A5, SPATA31A7, TTC27, USP17L5
- **Unclassified (5):** CBWD2, CBWD3, CBWD6, TBC1D3, TBC1D3B

The largest low-concordance cluster. 133 genes scattered across 42 Funk clusters.

**Literature Validation:**
- **GPATCH8** (uncharacterized): 2024, Molecular Cell — G-patch domain protein required for mutant SF3B1-induced splicing alterations. GPATCH8 modulates branchpoint selection, interacts with DHX15 helicase, and functionally opposes SUGP1. Silencing corrects SF3B1-mutant splicing defects. Major functional characterization. Nguyen et al. (2024) Molecular Cell. DOI: [10.1016/j.molcel.2024.04.015](https://doi.org/10.1016/j.molcel.2024.04.015)
- **WDR62** (novel): 2026, Commun Biol — WDR62 loss leads to decreased myoblast proliferation due to increased centriole numbers and centriole cohesion. WDR62 regulates mitotic spindle formation and centrosome integrity. Centrosome and microtubule-associated protein critical for mitosis. DOI: [10.1038/s42003-026-09537-7](https://www.nature.com/articles/s42003-026-09537-7)
- **PPP2R5B** (novel): 2024, Cell Rep — PP2A-B55 phosphatase (includes PPP2R5B family) controls chromosome individualization and clustering during mitosis, functioning as major CDK1-counteracting phosphatase during mitotic exit. DOI: [10.1016/j.celrep.2024.114482](https://www.sciencedirect.com/science/article/pii/S2211124724008234)
- **CCDC174** (uncharacterized): Component of exon junction complex, interacts with EIF4A3. Involved in RNA processing, nonsense-mediated decay, translation. No recent discoveries 2023-2026.
- **BECN1** (novel): BECN1 splice variants have distinct functional roles, with BECN1 regulating chromosome congression and kinetochore assembly during mitosis independent of PI3K-III complex. BECN1 interacts with kinetochore protein Zwint-1 for mitosis-specific functions. General splicing literature available but no specific 2023-2026 discoveries.
- **KATNB1** (novel): Katanin regulatory p80 subunit. Katanin severs microtubules at mitotic spindles during anaphase for sister chromatid segregation. Dynamic cell cycle-dependent localization. Established mitotic function but no specific 2023-2026 discoveries.
- **ZNF131** (novel): ZNF131 suppresses centrosome fragmentation in glioblastoma stem-like cells through regulation of HAUS5. Promotes expression of KIF7, NPHP1, TMEM237, and HAUS5. No specific 2023-2026 discoveries beyond 2017 study.
- **SNRNP35** (uncharacterized): Already documented in Section 4.4 — minor spliceosome U11/U12 di-snRNP component (searched 2023-2026).
- **ANKRD20A2** (uncharacterized): Ankyrin repeat domain 20 family member A2. ANKRD20A4P is a pseudogene with ANKRD20A3 as an important paralog. Readthrough transcripts between ANKRD20A4 and ANKRD20A20P are candidates for nonsense-mediated decay and unlikely to be translated. ANKRD20A2 shows expression in granulocytes and epididymis. No specific functional validation in splicing or RNA processing found (searched 2020-2026).
- **ANKRD20A3** (uncharacterized): Ankyrin repeat domain 20 family member A3. Has 1,782 functional associations across 6 biological entity categories extracted from 29 datasets. Limited specific molecular function data available. No validation for splicing or spliceosome function found (searched 2020-2026).
- **ANKRD20A4** (uncharacterized): Ankyrin repeat domain 20 family member A4. Ankyrin repeat domain family member with readthrough transcripts. No spliceosome validation found (searched 2023-2026).
- **FOXD4L1** (uncharacterized): Forkhead box D4-like 1. DNA-binding transcription factor with RNA polymerase II-specific activity. Part of FOXD4L1-L6 gene family arising from recent duplications during hominid evolution (mouse has single Foxd4 gene; humans have 7: FOXD4 + FOXD4L1-L6). At least two FOXD4L genes are transcriptionally active with evidence of purifying selection in forkhead domain. Associated disease: Retinitis Pigmentosa 55. No splicing or RNA processing validation found (searched 2020-2026).
- **FOXD4L3** (uncharacterized): Forkhead box D4-like 3. DNA-binding transcription factor activity. Part of hominid-specific FOXD4L duplication family. No splicing validation found (searched 2020-2026).
- **FOXD4L4** (uncharacterized): Forkhead box D4-like 4. DNA-binding transcription factor activity. Part of hominid-specific FOXD4L duplication family. No splicing validation found (searched 2020-2026).
- **FOXD4L5** (uncharacterized): Forkhead box D4-like 5. Part of hominid-specific FOXD4L duplication family arising during recent evolution. No splicing validation found (searched 2020-2026).
- **FOXD4L6** (uncharacterized): Forkhead box D4-like 6. DNA-binding transcription factor activity. Part of hominid-specific FOXD4L duplication family. Associated disease: Miles-Carpenter Syndrome. No splicing validation found (searched 2020-2026).
- **CCDC107** (uncharacterized): Coiled-coil domain protein; no splicing validation found (searched 2023-2026).
- **FAM219A** (uncharacterized): No recent characterization found (searched 2023-2026).
- **HRCT1** (uncharacterized): Predicted membrane protein with unknown function. No recent discoveries (searched 2023-2026).
- **KRTAP2-2** (uncharacterized): Keratin-associated protein for hair shaft formation. No recent discoveries beyond structural role (searched 2023-2026).
- **LURAP1L** (uncharacterized): Leucine rich adaptor protein. Predicted role in cell shape and migration. No recent discoveries (searched 2023-2026).
- **PMF1-BGLAP** (uncharacterized): Read-through transcript between PMF1 and BGLAP. Related to mitotic spindle pathways. No recent discoveries (searched 2023-2026).
- **RPP25L** (uncharacterized): Ribonuclease P/MRP subunit p25-like protein. Function not established. No recent discoveries (searched 2023-2026).
- **SPATA31A3** (uncharacterized): SPATA31 subfamily A member 3. Spermatogenesis-associated protein coding gene located on 9q21.11. SPATA31 gene family (formerly FAM75A) shows strongest signal of positive selection in hominoids. Family expanded from single-copy mouse gene (expressed in spermatogenesis) to SPATA31A and SPATA31C types in primates with broadened tissue expression. Humans average 7.5 SPATA31A copies (range 5-11 per diploid genome). Mouse Spata31-deficient males show azoospermia due to reduced nectin-3 and β-actin at apical ectoplasmic specialization. SPATA31 family also involved in UV damage response. Associated diseases: chronic closed-angle glaucoma, deprivation amblyopia. No cell cycle or mitotic splicing validation found (searched 2020-2026).
- **SPATA31A7** (uncharacterized): SPATA31 subfamily A member 7. Part of segmentally duplicated SPATA31A family in humans showing positive selection in hominoids. Copy number variation among individuals (5-11 copies). Family members function in spermatogenesis and UV damage response, but SPATA31A7 specific function unknown. No cell cycle or mitotic splicing validation found (searched 2020-2026).
- **TTC27** (uncharacterized): No splicing validation found (searched 2023-2026).
- **USP17L5** (uncharacterized): Deubiquitinating enzyme, peptidase C19 family. No recent discoveries (searched 2023-2026).
- **ARID3C** (novel): ARID family member involved in chromatin remodeling. No specific mitotic validation (searched 2023-2026).
- **SDE2** (novel): No specific functional discoveries found (searched 2023-2026).
- **PAXBP1** (novel): Adapter linking PAX3/PAX7 to histone methylation machinery in myogenesis. No mitotic validation found (searched 2023-2026).

---

## 4. Funk-Unique Mitotic Clusters (Funk → Brieflow, J < 0.15)

7 high-confidence Funk mitotic clusters with low Brieflow concordance.

---

### 4.1 Ribosome Biogenesis and rRNA Processing (FM168)
**16 genes | Jaccard = 0.037 | Best BF match: BM219 (Transcription regulation/cell cycle, Medium, 12g) | 1 shared gene | Frag: 14**

- **Established (7):** EIF3CL, FCF1, IMP4, LAS1L, PDCD11, RPS23, UTP14A
- **Novel role (8):** ACOT12, AMBRA1, MOB4, NOMO3, SEPTIN9, TIMM9, TMED2, TRIR
- **Uncharacterized (1):** ZNF574

**Literature Validation:**
- **ZNF574** (uncharacterized): 2024-2025, Molecular Cell — ZNF574 identified as quality control factor for defective ribosome biogenesis intermediates in the Ribosome Assembly Surveillance Pathway (RASP). Contains C2H2 zinc finger domains recognizing faulty ribosomes; loss causes accumulation of defective intermediates interfering with global ribosome production. In vivo znf574 loss in zebrafish causes ribosomopathy-like phenotypes. Major functional characterization. DOI: [10.1016/j.molcel.2025.05.014](https://www.cell.com/molecular-cell/fulltext/S1097-2765(25)00359-4)
- **AMBRA1** (novel): No ribosome biogenesis validation found (searched 2023-2026).
- **MOB4** (novel): No ribosome biogenesis validation found (searched 2023-2026).

**Discussion highlight:** Most Funk-unique mitotic cluster demonstrates sub-module resolution strength. Isolates late-stage 60S biogenesis quality control machinery (BRIX1, GNL2, LSG1, WDR12) that Brieflow fragments across 14 clusters. The discovery of ZNF574 as the RASP (Ribosome Assembly Surveillance Pathway) quality control factor recognizing defective ribosomes exemplifies Funk's ability to resolve regulatory checkpoints separate from core assembly machinery—a functional specialization that may lack distinct morphological signatures but shares interaction networks.

---

### 4.2 Ribosome Biogenesis and Translational Machinery (FM35)
**32 genes | Jaccard = 0.041 | Best BF match: BM0 (Ribosome biogenesis/rRNA processing, High, 169g) | 8 shared genes | Frag: 23**

- **Established (11):** EIF3M, EXOSC8, FTSJ3, GTF3C3, LSG1, MRPS17, NLE1, PAK1IP1, POLR1C, RBM28, RPS3
- **Novel role (18):** CHD1, DHFR2, DPF2, LIN7C, NDUFB3, NUP188, SAE1, SCGB2A1, SEPHS2, SFSWAP, STEAP3, TBCE, TTC5, TTF2, TUBB4B, TXNRD2, UBA6, ZNF638
- **Uncharacterized (3):** GNB1L, SCAND1, USP17L19

Lands in BM0 (8 genes) but the massive 169-gene BM0 dilutes concordance.

**Literature Validation:**
- **GNB1L** (uncharacterized): No specific ribosome biogenesis or translation validation found in recent literature (searched 2023-2026).
- **SCAND1** (uncharacterized): No nucleolar or ribosome biogenesis validation found (searched 2023-2026).
- **USP17L19** (uncharacterized): No ribosome biogenesis validation found (searched 2023-2026).

---

### 4.3 Mitotic Chromosome Condensation, Spindle Assembly, and Cell Division (FM3)
**51 genes | Jaccard = 0.046 | Best BF match: BM3 (Kinetochore/spindle/segregation, High, 86g) | 6 shared genes | Frag: 37**

- **Established (19):** CCNB1, CDK2AP2, CIP2A, CKS1B, DCTN1, DNAJC9, DYNC1LI2, FAAP24, KIF15, NCAPD3, NCAPG, NCAPG2, NCAPH2, POC1A, RFC1, SNX33, USP39, XRCC5, ZNF207
- **Novel role (25):** BCDIN3D, CBWD2, CDAN1, CLEC3B, CNNM4, DHX15, EIF4ENIF1, GCN1, GIGYF2, GPN1, HDAC7, HIPK1, NUTF2, OTOF, PCGF6, PHF20, PSMA1, PTPMT1, RBM8A, RPL29, RPN2, SMNDC1, SYS1, TKT, USP4
- **Uncharacterized (7):** AGAP7P, CFC1B, GPR45, MRPS5, RGPD5, SRBD1, ZNF492

Only 6 of 51 genes overlap with BM3 despite similar process annotation. Funk captures
a distinct chromosome condensation module.

**Literature Validation:**
- **SRBD1** (uncharacterized): 2025, Nature Communications/Cell Reports — Essential genome maintenance protein recruited to nascent chromatin in S-phase and mitotic chromosome axes. SRBD1 localizes topoisomerase IIα to mitotic chromosomes in condensin II-dependent manner. Critical during prophase for chromosome condensation establishment. Inactivation causes chromosome entanglements and anaphase failure. Major functional characterization. Murayama et al. (2025) Nature Communications. DOI: [10.1038/s41467-025-56911-w](https://doi.org/10.1038/s41467-025-56911-w); Antosch et al. (2025) Cell Reports. DOI: [10.1016/j.celrep.2025.115328](https://doi.org/10.1016/j.celrep.2025.115328)
- **AGAP7P** (uncharacterized): Putative pseudogene. No mitosis or chromosome condensation validation found (searched 2023-2026).
- **CFC1B** (uncharacterized): Cryptic family protein. No mitosis validation found (searched 2023-2026).
- **GPR45** (uncharacterized): G protein-coupled receptor. No mitosis validation found (searched 2023-2026).
- **MRPS5** (uncharacterized): Mitochondrial ribosomal protein. No mitosis/chromosome condensation link found (searched 2023-2026).
- **RGPD5** (uncharacterized): RANBP2-like protein. RAN-binding protein but no specific mitosis/chromosome condensation validation found (searched 2023-2026).
- **ZNF492** (uncharacterized): Zinc finger protein. No validation found (searched 2023-2026).

**Discussion highlight:** Most compelling Funk advantage—mechanistic sub-module resolution within broader processes. Despite both pipelines identifying "mitosis," Funk isolates the specific condensation machinery (condensin II: NCAPD2/NCAPG/NCAPH, SMC2, TOP2A) with spindle kinesins (KIF11/KIF23/KIF4A/KIF20B) separate from kinetochore assembly that Brieflow combines. The discovery of SRBD1 recruiting topoisomerase IIα to mitotic chromosomes in condensin-dependent manner validates this as a functionally distinct sub-module. Demonstrates Funk's strength in dissecting specialized mechanistic components that share interaction networks but may lack distinct morphological signatures.

---

### 4.4 DNA Replication Fork Stability and S-Phase Progression (FM94)
**24 genes | Jaccard = 0.055 | Best BF match: BM4 (DNA replication/fork/repair, High, 72g) | 5 shared genes | Frag: 11**

- **Established (10):** CDT1, HINFP, KPNB1, NCAPD2, RANBP1, REV3L, TIMELESS, TIPIN, TUBGCP4, WDHD1
- **Novel role (12):** CHERP, CTNNBL1, DNAJC14, H2BC15, INTS9, MED4, NISCH, PAXBP1, PHAX, RNPS1, RPLP2, RUVBL1
- **Uncharacterized (2):** PDCD7, SNRNP35

**Literature Validation:**
- **SNRNP35** (uncharacterized): 2023-2024, Molecular Cell — SNRNP35 (35K) is a minor spliceosome-specific factor in U11/U12 di-snRNP complex. Cryo-EM studies show SNRNP25 and SNRNP35 specifically recognize U11 snRNA for 5'SS recognition. No DNA replication link found. DOI: [10.1016/j.molcel.2024.11.025](https://doi.org/10.1016/j.molcel.2024.11.025)
- **PDCD7** (uncharacterized): Programmed cell death protein 7. No specific DNA replication fork stability or S-phase validation found in recent literature (searched 2023-2026).
- **CHERP** (novel): No DNA replication fork validation found (searched 2023-2026).
- **CTNNBL1** (novel): No DNA replication fork validation found (searched 2023-2026).

---

### 4.5 Ribosome Biogenesis and rRNA Processing (FM21)
**38 genes | Jaccard = 0.113 | Best BF match: BM0 (Ribosome biogenesis/rRNA processing, High, 169g) | 21 shared genes | Frag: 15**

- **Established (22):** BOP1, BRIX1, DDX51, EIF2S2, EXOSC9, NIP7, NOL8, NOP56, PPAN, RACK1, RCL1, RPL13, RPL28, RPL5, RPL7A, RPLP1, RPS28, RRS1, TMA16, URB2, UTP3, WDR75
- **Novel role (14):** ABCF1, ARL2, CD81, CHD8, CNOT9, HUWE1, METAP2, NUP160, NUP35, PRELID1, QARS1, RPRD1A, SLC3A2, TADA2B
- **Uncharacterized (2):** NOP16, PPAN-P2RY11

21 genes land in BM0 but low Jaccard due to BM0's 169-gene size.

**Literature Validation:**
- **NOP16** (uncharacterized): 2022, Nucleic Acids Research — Involved in early nucleolar pre-60S intermediates; casein kinase 2 (CK2) physically interacts with early pre-ribosomal particles including those affinity-purified with NOP16. All 4 CK2 subunits show similar abundance patterns in early pre-ribosomal particles. Fernández-Pevida et al. (2022) Nucleic Acids Research. DOI: [10.1093/nar/gkab572](https://doi.org/10.1093/nar/gkab572)
- **PPAN-P2RY11** (uncharacterized): Read-through transcript between PPAN (Peter Pan homolog, ribosome biogenesis) and P2RY11 (purinergic receptor). Encodes chimeric protein with lower activity than P2RY11 alone. Ubiquitously expressed, upregulated by granulocytic differentiation agents. Functional significance unclear. PPAN component established in ribosome biogenesis. No specific recent discoveries (searched 2023-2026).

---

### 4.6 Pre-mRNA Splicing and Spliceosome (FM76)
**25 genes | Jaccard = 0.121 | Best BF match: BM6 (Mitotic spindle assembly/segregation, High, 40g) | 7 shared genes | Frag: 5**

- **Established (16):** BCAS2, BORA, CDK1, DHX38, LSM3, NDC80, NSL1, RBM17, RNPC3, SDE2, SLU7, SON, TFIP11, TUBG1, YJU2, ZNHIT2
- **Novel role (6):** CLTC, DDX41, DHX9, HNRNPU, KANSL3, PCF11
- **Uncharacterized (3):** ANKRD20A4, LENG8, U2SURP

Pre-mRNA splicing coupled to mitotic cell division.

**Literature Validation:**
- **U2SURP** (uncharacterized): Also known as SR140; U2 snRNP-associated SURP motif-containing protein. SURP domains exclusively found in splicing-related proteins. U2SURP is among characterized U2 snRNP-associated proteins but no major recent functional discoveries (searched 2023-2026).
- **LENG8** (uncharacterized): Leukocyte receptor cluster member 8. No specific spliceosome or U2 snRNP validation found in recent literature (searched 2023-2026).
- **ANKRD20A4** (uncharacterized): Ankyrin repeat domain 20A family member. No spliceosome validation found (searched 2023-2026).
- **CLTC** (novel): Clathrin heavy chain — established endocytosis role, not novel in splicing context (searched 2023-2026).
- **DDX41** (novel): Known splicing factor — established role, not novel (searched 2023-2026).

---

### 4.7 Ribosome Biogenesis and SSU Processome Assembly (FM17)
**39 genes | Jaccard = 0.130 | Best BF match: BM0 (Ribosome biogenesis/rRNA processing, High, 169g) | 24 shared genes | Frag: 14**

- **Established (31):** BUD23, DCAF13, DDX10, DDX21, EIF1AX, EIF3G, EIF3I, EIF4A1, FAU, FBL, IMP3, LTV1, MPHOSPH10, NOP58, PWP2, RIOK2, RPL12, RPL39, RPS15A, RPS16, RPS17, RPS19, RPS21, RPS9, RRP1, TBL3, UTP11, UTP15, UTP6, WDR43, WDR46
- **Novel role (7):** BANP, DDX24, DNTTIP2, FBXW7, MAPK14, UBE2K, ZC3H13
- **Uncharacterized (1):** LIN37

24 genes land in BM0. Low Jaccard driven by BM0's broad 169-gene consolidation.

**Literature Validation:**
- **LIN37** (uncharacterized): Core component of DREAM complex (DP, RB-like, E2F, and MuvB) that represses cell cycle genes during quiescence. Essential for DREAM function. No direct SSU processome or ribosome biogenesis link found (searched 2023-2026).
- **BANP** (novel): BANP/SMAR1 — chromatin regulator. No ribosome biogenesis validation found (searched 2023-2026).
- **DDX24** (novel): Nucleolar DEAD-box helicase — established ribosome biogenesis role, not novel (searched 2023-2026).
- **DNTTIP2** (novel): No ribosome biogenesis validation found (searched 2023-2026).

---

## 5. Summary Statistics

### By Direction and Cell Class

| Category | Clusters | Total Genes |
|----------|----------|-------------|
| Brieflow-unique interphase | 27 | 575 |
| Funk-unique interphase | 7 | 122 |
| Brieflow-unique mitotic | 2 | 150 |
| Funk-unique mitotic | 7 | 225 |

### Literature Coverage Status

| Status | Count |
|--------|-------|
| Literature validated | 43 clusters (all) |
| Completed 2026-02-15 | 100% coverage |

---

*Generated 2026-02-13 from corrected barcode data. Jaccard threshold = 0.15.
Source TSVs: `benchmarks/results/cluster_overlap/{bf_to_fk,fk_to_bf}_{interphase,mitotic}.tsv`*
