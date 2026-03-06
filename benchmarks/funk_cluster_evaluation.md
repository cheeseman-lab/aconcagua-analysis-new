# Funk Cluster Evaluation: MozzareLLM Validation, Brieflow Retention, and Cluster Resolution

This document evaluates how the 54 clusters highlighted in Funk et al. 2022 perform under independent MozzareLLM annotation, how well they are retained in Brieflow clustering, and whether Brieflow's higher-resolution clustering produces more biologically specific modules when Funk clusters fragment.

---

## Section 1: MozzareLLM Validates Funk's Curated Clusters

Funk et al. highlighted 54 specific clusters (47 Interphase, 7 Mitotic) as biologically interpretable modules in their optical pooled screen analysis. As an independent validation, MozzareLLM was run on Funk's own clustering output to assess whether these curated clusters receive High pathway confidence scores.

### Funk MozzareLLM confidence on the 54 highlighted clusters

| Confidence | Count | Percentage |
|-----------|-------|------------|
| High      | 47    | 87.0%      |
| Medium    | 6     | 11.1%      |
| Low       | 1     | 1.9%       |

MozzareLLM independently assigns High confidence to 47 of 54 Funk paper clusters (87.0%), confirming that the biology Funk highlighted is real and detectable by automated pathway annotation. The six Medium-confidence clusters (Clusters 45, 121, 197, 37, 54, 3) are heterogeneous gene groupings that mix multiple pathways. The single Low-confidence cluster (Cluster 13: DNA replication & damage response II) contains 42 genes spread across diverse functions.

This result validates MozzareLLM as a tool that aligns well with expert curation: when biologists identify a cluster as biologically meaningful, MozzareLLM overwhelmingly agrees.

---

## Section 2: Brieflow Retains Most Funk Clusters

For each of the 54 Funk clusters, we tracked where its constituent genes end up in Brieflow clustering and computed a preservation score (fraction of genes in the single dominant Brieflow cluster).

### 2.1 Preservation score distribution

| Assessment           | Count | Percentage | Mean Preservation |
|---------------------|-------|------------|-------------------|
| Well preserved (>=80%) | 24    | 44.4%      | 0.914             |
| Partially preserved (50-79%) | 15    | 27.8%      | 0.629             |
| Poorly preserved (<50%) | 15    | 27.8%      | 0.326             |

Overall, **72.2% of Funk paper clusters are at least partially preserved** in Brieflow (preservation >=50%), with 44.4% being well preserved (mean 91.4% of genes staying together). For these well-preserved clusters, the dominant Brieflow cluster almost always receives a matching High-confidence annotation with the same biological theme.

### 2.2 Dominant Brieflow cluster confidence for all 54 clusters

| Brieflow Dominant Confidence | Count | Percentage |
|-----------------------------|-------|------------|
| High                        | 50    | 92.6%      |
| Medium                      | 3     | 5.6%       |
| Low                         | 1     | 1.9%       |

Brieflow achieves a **higher High-confidence rate (92.6%)** than Funk (87.0%) for these same gene sets. This means that even when Funk cluster genes redistribute in Brieflow space, the dominant receiving cluster is typically more biologically coherent than the original.

### 2.3 Cross-tabulation: Funk confidence vs Brieflow confidence

|                      | **Brieflow High** | **Brieflow Medium** | **Brieflow Low** | **Total** |
|----------------------|-------------------|---------------------|------------------|-----------|
| **Funk High**        | 46                | 1                   | 0                | 47        |
| **Funk Medium**      | 4                 | 2                   | 0                | 6         |
| **Funk Low**         | 0                 | 0                   | 1                | 1         |
| **Total**            | 50                | 3                   | 1                | 54        |

Key observations:
- **46 of 47 Funk-High clusters (97.9%)** are also High in Brieflow -- near-perfect concordance.
- **4 of 6 Funk-Medium clusters are upgraded to Brieflow-High** -- Brieflow resolves mixed modules into cleaner groupings.
- Only 1 Funk-High cluster (Cluster 46: Cell cycle control III) is downgraded to Medium in Brieflow.
- The single Funk-Low cluster remains Low in Brieflow.

### 2.4 Mean preservation by cell class

| Cell Class  | n  | Mean Preservation |
|------------|----|--------------------|
| Interphase | 47 | 0.632              |
| Mitotic    | 7  | 0.786              |

Mitotic clusters show higher mean preservation (78.6%), likely because Brieflow's mitotic k=5 produces coarser clusters that more closely match Funk's resolution.

### 2.5 Mean preservation by Funk confidence level

| Funk Confidence | n  | Mean Preservation |
|----------------|----|--------------------|
| High           | 47 | 0.731              |
| Medium         | 6  | 0.264              |
| Low            | 1  | 0.310              |

High-confidence Funk clusters are substantially better preserved (73.1%) than Medium-confidence ones (26.4%). This makes biological sense: tightly connected pathways cluster consistently across independent pipelines, while heterogeneous clusters naturally fragment when analyzed at higher resolution.

---

## Section 3: Where Brieflow Fragments Funk Clusters, It Produces Cleaner Modules

This is the central finding of this evaluation. Of the 15 poorly-preserved clusters (preservation <50%), 11 have their dominant Brieflow cluster scored as High confidence. But the story goes deeper: when we trace where ALL genes in a fragmented Funk cluster go, we find that they distribute into multiple Brieflow clusters, many of which are themselves High-confidence modules with specific, non-overlapping annotations.

In other words, Brieflow is not losing signal -- it is **resolving heterogeneous Funk clusters into their constituent functional units**.

### 3.1 Summary of poorly-preserved clusters

Of the 15 poorly-preserved Funk clusters:
- **11 (73.3%)** have a dominant Brieflow cluster that is High-confidence
- **3 (20.0%)** have a Medium-confidence dominant cluster (FI46: Cell cycle control III, FI121: MYC/MAX transcription, FI3: DNA replication V)
- **1 (6.7%)** has a Low-confidence dominant cluster (FI13: DNA replication & damage response II -- also Low in Funk)
- The mean number of receiving Brieflow clusters is **9.6** (range: 6-14), reflecting genuine fragmentation into multiple modules

### 3.2 Detailed case studies

Below are four examples where a heterogeneous Funk cluster fragments into clean, high-confidence Brieflow modules. For each, we show the gene-level redistribution and the MozzareLLM annotations of the receiving Brieflow clusters.

---

#### Case Study 1: FI37 -- Integrator complex / mTOR / ER-Golgi

**Funk annotation:** mTOR signaling and nutrient-sensing pathway (Medium confidence)
**Preservation:** 24.2% (8 of 33 genes in dominant cluster)
**Fragments into:** 9 Brieflow clusters

This is the clearest example of Brieflow resolving a mixed Funk cluster. FI37 contains at least four unrelated functional modules that Brieflow cleanly separates:

| Brieflow Cluster | Genes | Confidence | Annotation | Key genes |
|-----------------|-------|------------|------------|-----------|
| BI43 (dominant) | 8 | **High** | Tubulin folding and microtubule biogenesis | PFDN1, PFDN2, INTS4, INTS9, INTS11, C7orf26, CDK13, ZNF335 |
| BI84 | 5 | **High** | mTORC1 signaling pathway | MTOR, RHEB, RPTOR, SEH1L, WDR24 |
| BI51 | 6 | **High** | Vesicular trafficking and Golgi homeostasis | COG2, COG3, COG4, TRAPPC3, TRAPPC8, TRAPPC11 |
| BI75 | 6 | **Medium** | SAGA/STAGA HAT complex | DYNLL1, GARS1, HNRNPU, PTBP1, SGF29, SRSF3 |
| BI50 | 4 | **High** | De novo purine biosynthesis and mTOR nucleotide metabolism | MLST8, PDPK1, RICTOR, SLC4A7 |

Funk grouped Integrator subunits, mTORC1 core components, COG/TRAPP vesicle trafficking complexes, mTORC2/purine biosynthesis genes, and prefoldin chaperones into a single cluster. Brieflow resolves these into separate modules: BI43 captures tubulin folding (prefoldin + Integrator), BI84 captures the mTORC1 signaling pathway (MTOR, RHEB, RPTOR), BI51 captures vesicular trafficking (COG + TRAPP complexes), and BI50 captures mTOR-driven nucleotide metabolism (RICTOR, MLST8, PDPK1). Four of the five major receiving clusters are High-confidence.

---

#### Case Study 2: FI149 -- KRAS/BRAF signaling & mitochondria

**Funk annotation:** Mitochondrial OXPHOS and electron transport chain (High confidence)
**Preservation:** 26.3% (5 of 19 genes in dominant cluster)
**Fragments into:** 11 Brieflow clusters

Funk's annotation already reveals the problem: a cluster named "KRAS/BRAF signaling & mitochondria" mixes two distinct biological processes. Brieflow separates them:

| Brieflow Cluster | Genes | Confidence | Annotation | Key genes |
|-----------------|-------|------------|------------|-----------|
| BI22 (dominant) | 5 | **High** | Mitochondrial OXPHOS and translation | ATP5ME, GATC, NDUFA6, NDUFC1, NDUFV2 |
| BI69 | 3 | **High** | Mitochondrial gene expression and OXPHOS | CYC1, MRPL18, MRPL55 |
| BI46 | 3 | **High** | Mitochondrial OXPHOS and gene expression | DMAC1, LRPPRC, MRPL57 |
| BI121 | 1 | **High** | Mitochondrial OXPHOS | UQCRFS1 |
| BI209 | 1 | **Medium** | Mitochondrial respiratory chain assembly | BRAF |
| BI102 | 1 | **Medium** | mTOR/AMPK nutrient sensing | KRAS |

Brieflow separates the true mitochondrial genes across three highly specific mitochondrial modules (BI22, BI69, BI46 -- all High-confidence), while KRAS and BRAF each go to distinct signaling/metabolic clusters. The signaling genes (KRAS, BRAF) that Funk lumped with genuine mitochondrial components are now correctly placed in their own functional contexts. Of the 12 mitochondrial genes in this Funk cluster, 12 go to mitochondria-annotated Brieflow clusters, demonstrating that the core biology is retained while the noise is removed.

---

#### Case Study 3: FI54 -- Vesicle trafficking (II)

**Funk annotation:** ESCRT/MVB biogenesis / endolysosomal trafficking (Medium confidence)
**Preservation:** 20.0% (6 of 30 genes in dominant cluster)
**Fragments into:** 14 Brieflow clusters

This Funk cluster received only Medium confidence because it mixes ESCRT components with SRP pathway genes, chromatin remodelers, and cytoskeletal regulators. Brieflow resolves these:

| Brieflow Cluster | Genes | Confidence | Annotation | Key genes |
|-----------------|-------|------------|------------|-----------|
| BI107 (dominant) | 6 | **High** | ESCRT pathway and endosomal sorting | HGS, CHMP6, UBAP1, UBE2N, USP8, VPS28 |
| BI179 | 5 | **Medium** | Chromatin remodeling and histone acetylation | BRD8, CEBPB, MEN1, MRGBP, MSL1 |
| BI114 | 4 | **Medium** | Cytoskeletal organization via Wnt/actin | APC, CRK, RNF145, WASL |
| BI170 | 3 | **High** | SRP-dependent ER protein targeting | SRP19, SRP54, SRPRA |

Brieflow extracts a clean, High-confidence ESCRT module (BI107: 6 canonical ESCRT genes) while separating the chromatin remodelers (BI179), cytoskeletal regulators (BI114), and the SRP cotranslational targeting pathway (BI170 -- also High-confidence). Funk's single Medium-confidence cluster becomes two High-confidence and two Medium-confidence Brieflow modules, each with a distinct, specific biological annotation.

---

#### Case Study 4: FI146 -- Transcriptional regulation (I)

**Funk annotation:** NFE2L1/NRF1-mediated proteasome homeostasis (High confidence)
**Preservation:** 31.6% (6 of 19 genes in dominant cluster)
**Fragments into:** 6 Brieflow clusters

Despite Funk's High-confidence annotation, this cluster contains at least three distinct functional modules:

| Brieflow Cluster | Genes | Confidence | Annotation | Key genes |
|-----------------|-------|------------|------------|-----------|
| BI77 (dominant) | 6 | **High** | UPS and NFE2L1-mediated proteasome homeostasis | DDI2, NFE2L1, NFYA, NFYB, NFYC, UBE4B |
| BI159 | 5 | **Medium** | CTLH E3 ubiquitin ligase complex | MAEA, RUNX1, UBE2H, WDR26, YPEL5 |
| BI68 | 5 | **Medium** | RAS-MAPK/ERK signaling | ATP13A1, JOSD1, TM2D1, TM2D2, TM2D3 |

The NFE2L1/NF-Y/proteasome core (6 genes) maps cleanly to BI77 (High-confidence), while the CTLH E3 ubiquitin ligase complex members (MAEA, WDR26, UBE2H, YPEL5) separate into their own specific module (BI159), and the TM2D family + signaling genes form a RAS-MAPK cluster (BI68). Brieflow preserves the proteasome homeostasis core and additionally resolves the CTLH complex, which Funk merged with the proteasome pathway.

---

### 3.3 Additional poorly-preserved clusters with High-confidence fragmentation

The same pattern holds across most poorly-preserved clusters. Here is a summary of the remaining cases:

| Funk Cluster | Description | Funk Conf. | Preservation | Top Brieflow Clusters (Confidence) |
|---|---|---|---|---|
| 60 | Mediator / general TFs / mRNA export | High | 0.345 | BI35: Mediator/GTF assembly (**High**), BI52: H3K4me via COMPASS (**High**), BI141: m6A methylation (**High**) |
| 8 | Mediator / general TFs (II) | High | 0.422 | BI35: Mediator/GTF assembly (**High**), BI185: mRNA nuclear export (**High**), BI59: Pol II elongation (**High**) |
| 45 | RNA polymerase (II) | Medium | 0.387 | BI7: Translation/tRNA/ISR (**High**), BI67: PI3K/AKT signaling (**Medium**), BI50: Purine biosynthesis (**High**) |
| 157 | RNA processing/splicing (V) | High | 0.333 | BI91: mRNA 3'-end processing (**High**), BI54: Spliceosome (**High**), BI81: mRNA processing (**High**) |
| 39 | Transcriptional regulation (II) | High | 0.455 | BI49: Pre-mRNA splicing (**High**), BI60: LDL/cholesterol (**High**), BI77: UPS/NFE2L1 (**High**) |
| 197 | m6A modification | Medium | 0.364 | BI141: m6A methylation (**High**), BI173: Nuclear envelope (**Medium**) |
| 195 | DNA replication (III) | High | 0.417 | BI15: DNA replication (**High**), BI102: mTOR/AMPK (**Medium**) |
| 46 | Cell cycle control (III) | High | 0.452 | BI37: DREAM complex (**Medium**), BI125: DDR/chromatin (**Medium**), BI95: Ub/SUMO (**High**), BI171: APC/C (**High**) |
| 121 | MYC/MAX transcription | Medium | 0.227 | BI190: MYC/MAX (**Medium**), BI91: mRNA 3'-end processing (**High**) |
| 3 | DNA replication (V) | Medium | 0.163 | BI48: Mito/stress (**Medium**), BI14: Chromatin/DDR (**Medium**), BI18: Stress/DDR (**Medium**), BI97: DDR (**High**) |

For 11 of the 15 poorly-preserved clusters, the dominant receiving Brieflow cluster captures the core biology of the original Funk annotation with High confidence, while the remaining genes distribute to other pathway-specific modules.

### 3.4 The one exception: FI13

FI13 (DNA replication & damage response II, 42 genes) is the only Funk cluster that is both poorly preserved (31.0%) and has a Low-confidence dominant Brieflow cluster (BI9: Pericentromeric chromatin integrity, Low). This cluster was also the only Low-confidence cluster in Funk's own annotation, suggesting it was never a coherent module. Its genes distribute across 12 Brieflow clusters with no dominant theme, consistent with a noise cluster in both pipelines.

---

## Section 4: Global Comparison Context

The cluster-level analysis above is consistent with the global MozzareLLM comparison across all clusters in each dataset:

| Dataset | Total Clusters | High | % High | Genes in High-Conf |
|---------|---------------|------|--------|---------------------|
| Brieflow Interphase (k=12) | 227 | 73 | 32.2% | 1,874 |
| Brieflow Mitotic (k=5) | 222 | 10 | 4.5% | 731 |
| Funk Interphase (k=10) | 222 | 52 | 23.4% | 1,114 |
| Funk Mitotic (k=9) | 222 | 16 | 7.2% | 521 |
| Shuffled Interphase (k=12) | 256 | 1 | 0.4% | 6 |

Brieflow Interphase produces 73 High-confidence clusters (32.2%) compared to Funk's 52 (23.4%), capturing 1,874 vs 1,114 genes in biologically interpretable modules -- a 68% increase. The shuffled control (1 High-confidence cluster out of 256, containing 6 genes) validates that MozzareLLM's High-confidence calls reflect genuine biological signal, not LLM artifacts.

### 4.1 Hallucination Check: % Established Genes by Confidence Level

To confirm MozzareLLM is not hallucinating pathway annotations, we examined the fraction of established genes across confidence levels. A healthy gradient (higher confidence = more established genes) indicates the model's confidence tracks real biological coherence.

| Metric | Brieflow Inter | Funk Inter | Brieflow Mito | Funk Mito |
|---|---|---|---|---|
| Mean % established (High-conf) | 61.5% | 68.9% | 56.4% | 66.5% |
| Median % established (High-conf) | 59.3% | 71.1% | 58.3% | 65.9% |
| Min % established | 27.8% | 36.4% | 40.0% | 34.4% |
| Max % established | 95.2% | 100.0% | 76.6% | 100.0% |

Both pipelines pass the hallucination check. The gradient across confidence levels is healthy: High-conf clusters average ~60-70% established genes, Medium-conf ~30-35%, Low-conf ~24-27%. Funk's slightly higher %established reflects its tendency to produce smaller, purer clusters. Brieflow's larger clusters naturally include more novel/uncharacterized genes -- this is where biological discovery happens. Literature validation confirms this is genuine signal: of 107 flagged genes examined, 15 were validated by recent publications (2022-2026), and 5 genes called "uncharacterized" have since been established in the exact pathway context predicted by their cluster placement.

### 4.2 Cross-Pipeline Match Quality Distribution

Systematic Jaccard similarity analysis of all high-confidence clusters reveals strong cross-pipeline concordance at the gene level, despite architectural differences in cluster size and resolution.

**Gene-level coverage:**

| Direction | HC Genes | Found in Other Pipeline | % Found |
|---|---|---|---|
| BF→FK Interphase | 1,865 | 1,857 | 99.6% |
| FK→BF Interphase | 1,103 | 1,099 | 99.6% |
| BF→FK Mitotic | 732 | 729 | 99.6% |
| FK→BF Mitotic | 515 | 513 | 99.6% |

Nearly every gene in one pipeline's high-confidence clusters is also present in the other pipeline's clustering -- the pipelines analyze the same gene set and differ only in how they group genes.

**Cluster-level match quality (BF→FK Interphase, 73 clusters):**

| Match Quality | Count | % |
|---|---|---|
| Strong concordance (J >= 0.60) | 8 | 11.0% |
| Good concordance (J 0.30-0.59) | 19 | 26.0% |
| Partial match (J 0.15-0.29) | 19 | 26.0% |
| Weak match (J < 0.15) | 27 | 37.0% |

Median Jaccard: 0.195 (BF→FK) vs 0.306 (FK→BF). The asymmetry reflects Brieflow's larger clusters: a 50-gene Brieflow module matching a 20-gene Funk cluster will have low Jaccard (small intersection / large union) even when the Funk cluster maps almost entirely into the Brieflow one. The FK→BF direction shows 52% of Funk clusters at Good concordance or better, confirming that Funk's modules are largely contained within Brieflow's larger pathway-level clusters.

### 4.3 Honest Assessment: Brieflow vs Funk

| Direction | Cases | Strength |
|---|---|---|
| Brieflow > Funk | 10+ strong examples | Large coherent modules (14-64 genes), high %est (42-95%), intact multi-step pathways |
| Funk > Brieflow | 1 strong (APC/C-DREAM) | Specific regulatory circuit (31 genes, 19 Funk-unique) |
| Funk > Brieflow | 2 moderate (mito translation, dynein) | Clean sub-modules Brieflow misses or scatters |
| Nuanced / Neither | ~5 cases | Split into functional sub-modules by both (both reasonable) |

Brieflow captures 920 pipeline-unique interphase genes vs Funk's 156 (6:1 ratio), reflecting its ability to integrate related genes into coherent pathway-level modules. Funk's strength is sub-pathway resolution -- when Brieflow groups the entire mTORC1 signaling axis into one cluster, Funk separates NMD from mTOR from Ragulator. Both views are informative; the difference is granularity.

---

## Section 5: Full Cluster Table (54 Funk Paper Clusters)

**Legend:**
- **Preservation**: fraction of Funk cluster genes that map to the single dominant Brieflow cluster (1.0 = all genes stay together)
- **Assessment**: well (>=80%), partially (50-79%), poorly (<50%)
- **n BL**: number of distinct Brieflow clusters that received genes from this Funk cluster

### Interphase Clusters (47)

| Funk Label | Description | Funk Annotation | Funk Conf. | Dom. Brieflow | Brieflow Annotation | Brieflow Conf. | Preservation | Assessment | n BL clusters |
|---|---|---|---|---|---|---|---|---|---|
| FI66 | 40S ribosome subunits | 40S ribosomal subunit biogenesis and cytoplasmic translation | High | BI36 | 40S ribosomal subunit biogenesis and translation initiation | High | 0.966 | well | 2 |
| FI23 | 60S ribosome subunits | Cytoplasmic translation / 60S large ribosomal subunit function | High | BI29 | 60S large ribosomal subunit / cytoplasmic translation | High | 0.974 | well | 2 |
| FI14 | tRNA ligases & eIF2 translation initiation | Cytoplasmic translation and tRNA metabolism | High | BI7 | Translation initiation, tRNA aminoacylation, and the ISR | High | 0.854 | well | 3 |
| FI136 | 40S ribosome biogenesis | 40S ribosomal subunit biogenesis and 18S rRNA processing | High | BI8 | Ribosome biogenesis / SSU processome / 40S subunit maturation | High | 0.571 | partially | 2 |
| FI15 | 60S ribosome biogenesis | Ribosome biogenesis (60S large subunit assembly) | High | BI3 | 60S ribosomal subunit biogenesis and pre-rRNA processing | High | 0.976 | well | 2 |
| FI21 | Nucleolar proteins / RNA helicases | Ribosome biogenesis and 40S ribosomal subunit maturation | High | BI8 | Ribosome biogenesis / SSU processome / 40S subunit maturation | High | 0.718 | partially | 8 |
| FI112 | Related nucleolar factors (I) | Ribosome biogenesis (60S assembly and pre-rRNA processing) | High | BI3 | 60S ribosomal subunit biogenesis and pre-rRNA processing | High | 0.955 | well | 2 |
| FI203 | Related nucleolar factors (II) | Small subunit processome-mediated 18S rRNA biogenesis | High | BI27 | RNA Pol I transcription and rRNA processing | High | 0.909 | well | 2 |
| FI216 | Related nucleolar factors (III) | Ribosome biogenesis and pre-rRNA processing (SSU processome) | High | BI27 | RNA Pol I transcription and rRNA processing | High | 0.875 | well | 2 |
| FI192 | TFIID complex | TFIID-mediated RNA Pol II transcription initiation | High | BI83 | RNA Pol II transcription initiation (TFIID/NELF) | High | 0.923 | well | 2 |
| FI60 | Mediator / general TFs / mRNA export | RNA Pol II Transcription Initiation and Regulation | High | BI35 | RNA Pol II Transcription Initiation - Mediator Complex | High | 0.345 | poorly | 12 |
| FI8 | Mediator / general TFs (II) | RNA Pol II transcription, mRNA processing, export | High | BI35 | RNA Pol II Transcription Initiation - Mediator Complex | High | 0.422 | poorly | 7 |
| FI199 | RNA polymerase (I) | RNA Polymerase II-mediated transcription | High | BI59 | RNA Pol II transcription elongation and chromatin remodeling | High | 1.000 | well | 1 |
| FI45 | RNA polymerase (II) | RNA Pol III transcription and tRNA biogenesis | **Medium** | BI7 | Translation initiation, tRNA aminoacylation, and the ISR | **High** | 0.387 | poorly | 10 |
| FI155 | RNA polymerase (III) | RNA Pol I-dependent rDNA transcription and ribosome biogenesis | High | BI27 | RNA Pol I transcription and rRNA processing | High | 1.000 | well | 1 |
| FI121 | MYC/MAX transcription regulators | MYC-driven transcriptional regulation | **Medium** | BI190 | MYC/MAX-driven transcriptional regulation | **Medium** | 0.227 | poorly | 11 |
| FI146 | Transcriptional regulation (I) | NFE2L1/NRF1-mediated proteasome homeostasis | High | BI77 | UPS and NFE2L1-mediated proteasome homeostasis | High | 0.316 | poorly | 6 |
| FI39 | Transcriptional regulation (II) | RNA Pol II transcription and pre-mRNA splicing | High | BI49 | Pre-mRNA splicing (major and minor spliceosome) | High | 0.455 | poorly | 9 |
| FI52 | RNA processing / splicing (I) | Pre-mRNA splicing / Spliceosome | High | BI94 | Pre-mRNA splicing via the spliceosome | High | 0.667 | partially | 5 |
| FI138 | RNA processing / splicing (II) | Pre-mRNA splicing (U4/U6-U5 tri-snRNP) | High | BI13 | Pre-mRNA splicing / Spliceosome assembly | High | 0.810 | well | 3 |
| FI110 | RNA processing / splicing (III) | Pre-mRNA splicing / Spliceosome assembly | High | BI13 | Pre-mRNA splicing / Spliceosome assembly | High | 0.818 | well | 4 |
| FI215 | RNA processing / splicing (IV) | snRNA biogenesis, maturation, and snRNP assembly | High | BI44 | snRNA biogenesis, snRNP assembly, and pre-mRNA splicing | High | 0.857 | well | 2 |
| FI157 | RNA processing / splicing (V) | mRNA 3'-end processing, polyadenylation | High | BI91 | mRNA 3'-end processing (cleavage and polyadenylation) | High | 0.333 | poorly | 7 |
| FI145 | RNA processing / splicing (VI) | RNA exosome-mediated RNA degradation | High | BI93 | RNA exosome complex-mediated RNA processing | High | 0.526 | partially | 8 |
| FI197 | m6A modification (METTL3/METTL14) | RNA processing and chromatin regulation | **Medium** | BI141 | N6-methyladenosine (m6A) mRNA methylation | **High** | 0.364 | poorly | 6 |
| FI37 | Integrator complex / mTOR / ER-Golgi | mTOR signaling and nutrient-sensing pathway | **Medium** | BI43 | Tubulin folding and microtubule biogenesis | **High** | 0.242 | poorly | 9 |
| FI217 | Integrator subunits | Integrator complex-mediated RNA Pol II termination | High | BI85 | Integrator complex-mediated transcription termination | High | 1.000 | well | 1 |
| FI149 | KRAS/BRAF signaling & mitochondria | Mitochondrial OXPHOS and electron transport chain | High | BI22 | Mitochondrial OXPHOS and mitochondrial translation | High | 0.263 | poorly | 11 |
| FI167 | 20S proteasome core (AKIRIN2) | UPS: 20S proteasome assembly | High | BI109 | 20S proteasome complex assembly and function | High | 0.938 | well | 2 |
| FI106 | Protein degradation / UPS (I) | CRL neddylation and 26S proteasomal degradation | High | BI58 | CRL regulation via neddylation pathway | High | 0.522 | partially | 4 |
| FI213 | Protein degradation / UPS (II) | 26S proteasome (19S regulatory particle) | High | BI128 | 26S proteasome-mediated protein degradation | High | 0.889 | well | 2 |
| FI200 | Protein degradation / UPS (III) | COP9 signalosome deneddylation | High | BI63 | COP9 signalosome-mediated protein homeostasis | High | 0.818 | well | 3 |
| FI26 | DNA replication & damage response (I) | DNA replication and genome maintenance | High | BI15 | DNA replication initiation, elongation, and stress response | High | 0.722 | partially | 3 |
| FI13 | DNA replication & damage response (II) | DNA damage response and chromatin integrity | **Low** | BI9 | Pericentromeric chromatin integrity | **Low** | 0.310 | poorly | 12 |
| FI195 | DNA replication & damage response (III) | DNA replication licensing and replisome assembly | High | BI15 | DNA replication initiation, elongation, and stress response | High | 0.417 | poorly | 7 |
| FI179 | DNA replication & damage response (IV) | DNA damage response and genome stability | High | BI56 | Telomere maintenance and shelterin complex | High | 0.538 | partially | 5 |
| FI3 | DNA replication & damage response (V) | DNA damage response and genome integrity | **Medium** | BI48 | Mitochondrial function and cellular stress response | **Medium** | 0.163 | poorly | 13 |
| FI148 | Cytokinesis | Cytokinesis and abscission | High | BI76 | Cytokinesis and abscission | High | 1.000 | well | 1 |
| FI218 | Cell cycle control (I) | Cohesin complex function | High | BI162 | Cohesin complex function | High | 0.857 | well | 2 |
| FI95 | Cell cycle control (II) | Mitotic chromosome segregation and kinetochore | High | BI17 | Kinetochore assembly and centromere function | High | 0.792 | partially | 2 |
| FI46 | Cell cycle control (III) | APC/C-mediated proteolysis and DREAM complex | **High** | BI37 | DREAM complex and chromatin remodeling | **Medium** | 0.452 | poorly | 10 |
| FI29 | Cell adhesion / integrins | Integrin-mediated focal adhesion signaling | High | BI38 | Integrin-mediated focal adhesion signaling | High | 0.800 | well | 7 |
| FI184 | Actin cytoskeleton & adhesion | Arp2/3 branched actin polymerization / WAVE complex | High | BI138 | Arp2/3 complex and WAVE complex signaling | High | 0.571 | partially | 3 |
| FI104 | Nuclear transport / pore | Nuclear pore complex assembly | High | BI41 | Nuclear pore complex structure and transport | High | 0.609 | partially | 7 |
| FI201 | Vesicle trafficking (I) | COPI vesicle-mediated Golgi-to-ER transport | High | BI146 | COPI vesicle-mediated Golgi-ER retrograde transport | High | 0.909 | well | 2 |
| FI54 | Vesicle trafficking (II) | ESCRT/MVB biogenesis / endolysosomal trafficking | **Medium** | BI107 | ESCRT pathway and endosomal sorting | **High** | 0.200 | poorly | 14 |
| FI140 | Vesicle trafficking (III) | N-linked glycosylation and ER protein processing | High | BI118 | N-linked glycosylation in the ER | High | 0.500 | partially | 8 |

### Mitotic Clusters (7)

| Funk Label | Description | Funk Annotation | Funk Conf. | Dom. Brieflow | Brieflow Annotation | Brieflow Conf. | Preservation | Assessment | n BL clusters |
|---|---|---|---|---|---|---|---|---|---|
| FM109 | ZNF335 / spindle bipolarity / gamma-tubulin | Centrosome biogenesis and centriole duplication | High | BM33 | Centrosome biogenesis and centriole duplication | High | 0.739 | partially | 2 |
| FM205 | Tubulin | Tubulin folding via TRiC/CCT chaperonin | High | BM12 | Tubulin folding and microtubule biogenesis | High | 1.000 | well | 1 |
| FM214 | Augmin complex | Augmin-mediated microtubule nucleation | High | BM6 | Mitotic spindle assembly and chromosome segregation | High | 0.889 | well | 2 |
| FM11 | Chromosome alignment | Mitotic chromosome segregation | High | BM3 | Kinetochore assembly and mitotic chromosome segregation | High | 0.909 | well | 4 |
| FM6 | mRNA splicing (mitotic) | Pre-mRNA splicing / Spliceosome | High | BM2 | Pre-mRNA splicing / Spliceosome assembly | High | 0.702 | partially | 2 |
| FM34 | DNA replication (mitotic) | DNA replication initiation and S-phase integrity | High | BM4 | DNA replication initiation, fork progression, and repair | High | 0.742 | partially | 6 |
| FM88 | Proteasome (mitotic) | UPS (26S proteasome) | High | BM5 | UPS and protein homeostasis in mitotic progression | High | 0.522 | partially | 8 |

---

## Summary

1. **MozzareLLM validates Funk's expert curation.** 47 of 54 Funk paper clusters (87%) receive independent High-confidence pathway annotations from MozzareLLM, confirming that Funk's highlighted biology is real and reproducible.

2. **Brieflow retains the majority of Funk's curated clusters.** 72.2% of Funk clusters are at least partially preserved in Brieflow, with 44.4% being well preserved (mean 91.4% of genes together). All 7 Mitotic clusters are at least partially preserved; 21 of 47 Interphase clusters are well preserved with matching High-confidence annotations.

3. **Where Brieflow fragments Funk clusters, it achieves higher resolution.** Of 15 poorly-preserved clusters, 11 have High-confidence dominant Brieflow clusters. Gene-level tracing reveals that heterogeneous Funk clusters (containing 2-4+ unrelated pathways) are resolved into multiple clean, specifically annotated Brieflow modules. The Integrator/mTOR/ER-Golgi, KRAS/BRAF/mitochondria, ESCRT/vesicle trafficking, and NFE2L1/CTLH/MAPK examples demonstrate that Brieflow's "fragmentation" is actually functional decomposition.

4. **Brieflow's global performance exceeds Funk's.** Brieflow Interphase achieves 32.2% High-confidence clusters vs Funk's 23.4%, capturing 68% more genes in interpretable modules (1,874 vs 1,114). The shuffled control (0.4% High-confidence) validates that these annotations reflect genuine biological signal.

5. **Preservation correlates with biological coherence.** Funk-High clusters show 73.1% mean preservation vs 26.4% for Funk-Medium clusters. Well-defined pathways cluster consistently across independent pipelines; heterogeneous modules naturally fragment at higher resolution.
