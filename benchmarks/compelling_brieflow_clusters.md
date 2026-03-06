# Compelling High-Confidence Clusters from Brieflow MozzareLLM Analysis

**Screen**: Aconcagua OPS, DAPI/TUBULIN/GH2AX/PHALLOIDIN
**Conditions**: Interphase (k=12) and Mitotic (k=5)
**MozzareLLM model**: claude-opus-4-6
**Date**: 2026-02-18

---

## Overall Statistics

### Interphase (k=12)
| Metric | Value |
|--------|-------|
| Total clusters | 227 |
| High-confidence clusters | 73 (32.2%) |
| Medium-confidence clusters | 60 (26.4%) |
| Low-confidence clusters | 94 (41.4%) |

### Mitotic (k=5)
| Metric | Value |
|--------|-------|
| Total clusters | 222 |
| High-confidence clusters | 10 (4.5%) |
| Medium-confidence clusters | 62 (27.9%) |
| Low-confidence clusters | 150 (67.6%) |

### Combined
| Metric | Interphase | Mitotic | Combined |
|--------|-----------|---------|----------|
| High-confidence clusters | 73 | 10 | 83 |
| Total genes in high-conf | ~1,874 | ~731 | ~2,605 |

The dramatic difference in high-confidence rates (32.2% vs. 4.5%) reflects the biology: interphase is the default cellular state where most pathways are actively operating and thus detectable by morphological profiling, while mitotic cells represent a specialized transient state where fewer pathways produce distinctive phenotypic signatures.

---

## Summary Table: All High-Confidence Clusters

### Interphase High-Confidence Clusters (73 clusters)

| Cluster | Process | Genes | Estab. | Novel | Unchar. | %Estab. |
|---------|---------|-------|--------|-------|---------|---------|
| BI3 | 60S ribosomal subunit biogenesis / pre-rRNA processing | 76 | 55 | 12 | 9 | 72.4% |
| BI7 | Translation initiation, tRNA aminoacylation, and ISR | 65 | 41 | 22 | 2 | 63.1% |
| BI8 | SSU processome / 40S subunit maturation | 66 | 49 | 11 | 6 | 74.2% |
| BI13 | Pre-mRNA splicing / Spliceosome assembly | 54 | 40 | 10 | 4 | 74.1% |
| BI15 | DNA replication initiation, elongation, stress response | 52 | 46 | 5 | 1 | 88.5% |
| BI17 | Kinetochore assembly, centromere function, spindle org. | 52 | 43 | 8 | 1 | 82.7% |
| BI22 | Mitochondrial OXPHOS and mitochondrial translation | 48 | 38 | 8 | 2 | 79.2% |
| BI27 | RNA Pol I transcription and rRNA processing | 44 | 36 | 7 | 1 | 81.8% |
| BI29 | 60S large ribosomal subunit / cytoplasmic translation | 42 | 39 | 2 | 1 | 92.9% |
| BI31 | Replication stress / FA/HR DNA repair | 41 | 33 | 5 | 3 | 80.5% |
| BI32 | CDC42/Rho-dependent actin organization | 41 | 17 | 22 | 2 | 41.5% |
| BI35 | RNA Pol II transcription initiation (Mediator/GTFs) | 41 | 33 | 6 | 2 | 80.5% |
| BI36 | 40S ribosomal subunit / translation initiation | 40 | 38 | 1 | 1 | 95.0% |
| BI38 | Integrin-mediated focal adhesion signaling | 39 | 16 | 21 | 2 | 41.0% |
| BI41 | Nuclear pore complex structure and transport | 36 | 15 | 19 | 2 | 41.7% |
| BI43 | Tubulin folding / microtubule biogenesis | 31 | 11 | 19 | 1 | 35.5% |
| BI44 | snRNA biogenesis, snRNP assembly, splicing | 31 | 18 | 10 | 3 | 58.1% |
| BI46 | Mitochondrial OXPHOS and gene expression | 30 | 22 | 7 | 1 | 73.3% |
| BI49 | Pre-mRNA splicing (major and minor spliceosome) | 32 | 17 | 12 | 3 | 53.1% |
| BI50 | De novo purine biosynthesis / mTOR nucleotide metab. | 29 | 15 | 12 | 2 | 51.7% |
| BI51 | Vesicular trafficking / Golgi homeostasis (TRAPP/COG/HOPS) | 29 | 24 | 4 | 1 | 82.8% |
| BI52 | H3K4 methylation / COMPASS-MLL complexes | 29 | 20 | 8 | 1 | 69.0% |
| BI54 | Pre-mRNA splicing / spliceosome function | 29 | 16 | 9 | 4 | 55.2% |
| BI56 | Telomere maintenance / shelterin + DDR | 28 | 12 | 8 | 8 | 42.9% |
| BI58 | CRL regulation via neddylation pathway | 28 | 9 | 17 | 2 | 32.1% |
| BI59 | RNA Pol II transcription elongation | 29 | 24 | 4 | 1 | 82.8% |
| BI61 | Chromatin structure / nucleosome organization | 27 | 16 | 8 | 3 | 59.3% |
| BI63 | Splicing + COP9 signalosome protein homeostasis | 27 | 18 | 6 | 3 | 66.7% |
| BI64 | ER protein biogenesis / membrane insertion / secretory | 25 | 13 | 8 | 4 | 52.0% |
| BI67 | PI3K/AKT + INO80 chromatin remodeling | 25 | 15 | 8 | 2 | 60.0% |
| BI69 | Mitochondrial gene expression and OXPHOS | 24 | 21 | 2 | 1 | 87.5% |
| BI74 | tRNA modification / Elongator wobble uridine pathway | 23 | 11 | 10 | 2 | 47.8% |
| BI76 | Cytokinesis and abscission | 23 | 15 | 8 | 0 | 65.2% |
| BI77 | UPS and NFE2L1-mediated proteasome homeostasis | 23 | 10 | 12 | 1 | 43.5% |
| BI79 | HAT complexes (ATAC, MOZ/MORF) | 22 | 7 | 13 | 2 | 31.8% |
| BI81 | mRNA 3'-end processing / deadenylation (CCR4-NOT) | 22 | 11 | 10 | 1 | 50.0% |
| BI82 | DNA replication stress / MCM/ATR checkpoint | 22 | 11 | 9 | 2 | 50.0% |
| BI83 | RNA Pol II transcription initiation (TFIID/NELF) | 22 | 15 | 6 | 1 | 68.2% |
| BI84 | mTORC1 signaling pathway | 22 | 13 | 8 | 1 | 59.1% |
| BI85 | Integrator complex / Pol II pause-release | 22 | 9 | 12 | 1 | 40.9% |
| BI89 | NuA4 HAT complex / chromatin remodeling | 21 | 10 | 8 | 3 | 47.6% |
| BI91 | mRNA 3'-end cleavage and polyadenylation | 21 | 10 | 10 | 1 | 47.6% |
| BI93 | RNA exosome complex | 21 | 11 | 8 | 2 | 52.4% |
| BI94 | Pre-mRNA splicing (spliceosome) - purest cluster | 21 | 20 | 0 | 1 | 95.2% |
| BI95 | Ub/SUMO conjugation / chromatin DDR | 20 | 14 | 5 | 1 | 70.0% |
| BI101 | Dynein-dynactin microtubule transport | 20 | 15 | 3 | 2 | 75.0% |
| BI106 | ER-to-Golgi vesicular trafficking | 19 | 12 | 6 | 1 | 63.2% |
| BI107 | ESCRT pathway / MVB sorting | 19 | 8 | 10 | 1 | 42.1% |
| BI109 | 20S proteasome complex assembly | 19 | 17 | 2 | 0 | 89.5% |
| BI115 | Autophagy / PI3K-mediated membrane trafficking | 18 | 5 | 10 | 3 | 27.8% |
| BI117 | V-ATPase endolysosomal acidification | 18 | 14 | 3 | 1 | 77.8% |
| BI118 | N-linked glycosylation in ER | 17 | 9 | 7 | 1 | 52.9% |
| BI121 | Mitochondrial OXPHOS / gene expression | 17 | 10 | 4 | 3 | 58.8% |
| BI128 | 26S proteasome (19S regulatory particle) | 17 | 15 | 2 | 0 | 88.2% |
| BI129 | Mitochondrial function / OXPHOS | 16 | 8 | 7 | 1 | 50.0% |
| BI138 | Arp2/3 + WAVE branched actin polymerization | 16 | 8 | 7 | 1 | 50.0% |
| BI139 | Nuclear pore complex / nucleocytoplasmic transport | 16 | 12 | 3 | 1 | 75.0% |
| BI141 | m6A mRNA methylation | 16 | 6 | 9 | 1 | 37.5% |
| BI145 | Centrosome biogenesis / centriole duplication | 15 | 6 | 8 | 1 | 40.0% |
| BI146 | COPI Golgi-ER retrograde transport | 15 | 11 | 4 | 0 | 73.3% |
| BI160 | Mitochondrial membrane organization | 14 | 8 | 6 | 0 | 57.1% |
| BI162 | Cohesin complex / chromatin organization | 14 | 6 | 7 | 1 | 42.9% |
| BI163 | TRiC/CCT chaperonin / tubulin-actin folding | 14 | 10 | 3 | 1 | 71.4% |
| BI165 | Nuclear pore complex / transport | 13 | 7 | 6 | 0 | 53.8% |
| BI170 | SRP cotranslational protein targeting to ER | 13 | 5 | 7 | 1 | 38.5% |
| BI171 | APC/C ubiquitin-dependent degradation | 13 | 9 | 2 | 2 | 69.2% |
| BI172 | Histone/PP4 chromatin and DDR | 12 | 8 | 3 | 1 | 66.7% |
| BI185 | mRNA nuclear export (TREX/NXF1) | 11 | 10 | 1 | 0 | 90.9% |
| BI186 | Condensin I complex / chromatin architecture | 11 | 5 | 5 | 1 | 45.5% |
| BI192 | DDR / checkpoint signaling (9-1-1 clamp) | 10 | 5 | 4 | 1 | 50.0% |
| BI197 | N-glycosylation / ERAD | 10 | 7 | 3 | 0 | 70.0% |
| BI198 | PE biosynthesis / Kennedy pathway | 10 | 5 | 3 | 2 | 50.0% |
| BI199 | Iron-sulfur cluster biogenesis / CIA pathway | 10 | 6 | 4 | 0 | 60.0% |
| BI217 | eIF5A hypusination pathway | 7 | 3 | 3 | 1 | 42.9% |

### Mitotic High-Confidence Clusters (10 clusters)

| Cluster | Process | Genes | Estab. | Novel | Unchar. | %Estab. |
|---------|---------|-------|--------|-------|---------|---------|
| BM0 | Ribosome biogenesis and rRNA processing | 171 | 131 | 30 | 10 | 76.6% |
| BM1 | Pre-mRNA splicing and RNA processing | 130 | 76 | 31 | 23 | 58.5% |
| BM2 | Pre-mRNA splicing / Spliceosome (mitotic) | 99 | 58 | 27 | 14 | 58.6% |
| BM3 | Kinetochore / spindle / chromosome segregation | 86 | 50 | 28 | 8 | 58.1% |
| BM4 | DNA replication / fork progression / repair | 72 | 43 | 18 | 11 | 59.7% |
| BM5 | UPS / protein homeostasis in mitotic progression | 55 | 22 | 31 | 2 | 40.0% |
| BM6 | Mitotic spindle assembly / chromosome segregation | 40 | 21 | 15 | 4 | 52.5% |
| BM12 | Tubulin folding / microtubule biogenesis (mitotic) | 33 | 20 | 9 | 4 | 60.6% |
| BM7 | DNA replication stress response | -- | -- | -- | -- | -- |
| BM21 | Chromatin integrity | -- | -- | -- | -- | -- |

*Note: Mitotic clusters 7 and 21 are additional high-confidence clusters included in the count of 10. Detailed profiles focus on the 8 clusters with the strongest biological signals.*

---

## Interphase Highlights

### Part I: Clean Biological Machines -- Validated Clusters with Follow-Up Targets

These clusters capture nearly complete representations of known molecular machines, with 70-100% established gene membership. The key insight: when an unbiased morphological screen reconstructs a known complex almost perfectly, any mystery gene sitting among the established members becomes a very high-confidence functional prediction.

---

#### 1. BI94 -- The Purest Spliceosome (Interphase)

**Process**: Pre-mRNA splicing via the spliceosome
**Genes**: 21 | **Established**: 20 | **Novel**: 0 | **Uncharacterized**: 1 | **%Established**: 95.2%

| Gene | Role |
|------|------|
| SF3A1, SF3A3 | established (SF3A subcomplex) |
| SF3B1, SF3B2, SF3B3, SF3B5 | established (SF3B subcomplex) |
| SNRPA1, SNU13 | established (U2 snRNP) |
| PLRG1 | established (PRP19/CDC5L complex) |
| SLU7, DHX8, RBM22, CDC40 | established (catalytic step factors) |
| AQR, CWC22, CRNKL1, XAB2 | established (spliceosomal components) |
| BUD31, SMU1, ISY1 | established (additional splicing factors) |
| ISY1-RAB43 | uncharacterized (read-through fusion transcript) |

**Why this is compelling**: This is the single purest cluster in the entire dataset -- 20 of 21 genes are canonical spliceosome components spanning the SF3A, SF3B, U2 snRNP, PRP19, and catalytic step machinery. The screen has essentially reconstructed the spliceosome from morphological phenotypes alone. The sole uncharacterized gene, ISY1-RAB43, is a read-through fusion transcript whose functional characterization as a distinct entity is lacking. Its phenotypic co-clustering with core spliceosome factors provides strong evidence it functions within the splicing pathway.

---

#### 2. BI29 -- The Complete 60S Ribosome (Interphase)

**Process**: 60S large ribosomal subunit / cytoplasmic translation
**Genes**: 42 | **Established**: 39 | **Novel**: 2 | **Uncharacterized**: 1 | **%Established**: 92.9%

| Gene | Role |
|------|------|
| RPL3, RPL4, RPL5, RPL6, RPL7, RPL7A, RPL8, RPL9, RPL10, RPL10A, RPL11, RPL12, RPL13, RPL13A, RPL14, RPL15, RPL17, RPL18, RPL18A, RPL19, RPL21, RPL22, RPL23, RPL23A, RPL24, RPL26, RPL27, RPL27A, RPL28, RPL29, RPL30, RPL31, RPL32, RPL34, RPL35, RPL36, RPL37A, RPL38, RPLP0 | established (60S ribosomal proteins) |
| RPL17-C18orf32 | uncharacterized (read-through fusion) |
| AAMP | novel_role (angiogenesis/migration protein) |

**Why this is compelling**: This cluster contains 39 of the ~47 human 60S ribosomal proteins -- essentially a complete large ribosomal subunit assembled purely from phenotypic similarity in a CRISPR screen. The single outlier gene, **AAMP**, is an angiogenesis and migration protein with no previously known ribosomal function. Its co-clustering with a virtually complete ribosome is a strong prediction that AAMP has an unrecognized role in ribosome biogenesis, translational regulation, or a shared nucleolar stress response pathway. This is a top-priority follow-up target.

---

#### 3. BI109 -- 20S Proteasome Complex (Interphase)

**Process**: 20S proteasome complex assembly and function
**Genes**: 19 | **Established**: 17 | **Novel**: 2 | **Uncharacterized**: 0 | **%Established**: 89.5%

| Gene | Role |
|------|------|
| PSMA1-7 | established (all 7 alpha subunits) |
| PSMB1-7 | established (all 7 beta subunits) |
| POMP | established (proteasome assembly chaperone) |
| PSMG3 | established (proteasome assembly chaperone) |
| AKIRIN2 | established (nuclear import of 20S) |
| ORC3 | novel_role (origin recognition complex) |
| PPIH | novel_role (peptidyl-prolyl isomerase/splicing) |

**Why this is compelling**: All 14 core 20S proteasome subunits (7 alpha + 7 beta) cluster together along with two assembly chaperones and the recently characterized nuclear import factor AKIRIN2. The two outliers, **ORC3** (DNA replication origin licensing) and **PPIH** (cyclophilin-type isomerase involved in splicing), have no established proteasome roles. ORC3's presence is particularly surprising and suggests either a novel connection to proteasome biology or shared downstream phenotypic consequences. PPIH's chaperone activity provides a plausible mechanistic hypothesis for involvement in proteasome folding/assembly.

---

#### 4. BI128 -- 26S Proteasome 19S Regulatory Particle (Interphase)

**Process**: 26S proteasome regulatory particle
**Genes**: 17 | **Established**: 15 | **Novel**: 2 | **Uncharacterized**: 0 | **%Established**: 88.2%

| Gene | Role |
|------|------|
| PSMC1-6 | established (all 6 AAA-ATPase subunits) |
| PSMD1, PSMD2, PSMD3, PSMD6, PSMD7, PSMD8, PSMD11, PSMD12, PSMD14 | established (lid/base non-ATPase subunits) |
| SEM1 | established (lid subunit, dual-function) |
| VCP | established (functional partner) |

**Why this is compelling**: The screen independently separated the 20S core particle (cluster 109) from the 19S regulatory particle (cluster 128), each with near-complete subunit representation. The 19S contains all six AAA-ATPase base subunits and nine non-ATPase subunits. This is a remarkable demonstration of the screen's resolution power -- it distinguishes subcomplexes of the same molecular machine based on distinct morphological signatures.

---

#### 5. BI185 -- mRNA Nuclear Export (Interphase)

**Process**: mRNA nuclear export via TREX complex and TAP/NXF1 pathway
**Genes**: 11 | **Established**: 10 | **Novel**: 1 | **Uncharacterized**: 0 | **%Established**: 90.9%

| Gene | Role |
|------|------|
| THOC1, THOC2, THOC5, THOC7 | established (THO subcomplex) |
| ALYREF (THOC4) | established (mRNA export adaptor) |
| DDX39B | established (TREX helicase) |
| NXF1, NXT1 | established (TAP-p15 heterodimer) |
| GLE1 | established (mRNA export factor) |
| MCM3AP | established (TREX-2 component) |
| ACTR5 | novel_role (INO80 chromatin remodeling) |

**Why this is compelling**: Ten of eleven genes are canonical mRNA nuclear export factors, spanning the TREX complex (THO subcomplex + ALYREF + DDX39B), the export receptor (NXF1/NXT1), and the mRNA release factor GLE1. The sole outlier, **ACTR5**, is a core subunit of the INO80 chromatin remodeling complex. Its phenotypic co-clustering with the complete mRNA export machinery makes it a compelling candidate for a novel role linking chromatin remodeling to mRNP biogenesis or nuclear export -- a connection with mechanistic plausibility given INO80's role in transcription but no prior experimental support.

---

#### 6. BI15 -- DNA Replication Machinery (Interphase)

**Process**: DNA replication initiation, elongation, and replication stress response
**Genes**: 52 | **Established**: 46 | **Novel**: 5 | **Uncharacterized**: 1 | **%Established**: 88.5%

| Gene | Role |
|------|------|
| ORC1, CDC6, CDT1, GMNN | established (origin licensing) |
| CDC7, TICRR, CDC45, GINS1-4, MCM10 | established (origin firing) |
| PCNA, RFC1-5, RPA1-3 | established (replication fork) |
| POL alpha/delta/epsilon subunits | established (polymerases) |
| TIMELESS, TIPIN, DONSON | established (fork protection) |
| DTL, DDB1, SKP2, CKS1B, FBXO5, LRR1 | established (ubiquitin regulation) |
| RRM1, RRM2, DTYMK, DUT | established (nucleotide biosynthesis) |
| RTEL1, MMS22L | established (repair) |
| TIMM29 | novel_role (mitochondrial import factor) |
| HNRNPH3 | novel_role (RNA processing) |

**Why this is compelling**: This cluster captures virtually every stage of eukaryotic DNA replication -- origin licensing through firing, elongation, nucleotide supply, fork protection, and replication-coupled ubiquitin regulation. The single most compelling follow-up target is **TIMM29**, a mitochondrial protein import factor with absolutely no known nuclear or replication function. Its presence among 46 established replication factors is unexplained by current literature and represents a potentially important biological connection worth experimental investigation.

---

#### 7. BI146 -- COPI Retrograde Transport (Interphase)

**Process**: COPI vesicle-mediated Golgi-ER retrograde transport
**Genes**: 15 | **Established**: 11 | **Novel**: 4 | **Uncharacterized**: 0 | **%Established**: 73.3%

| Gene | Role |
|------|------|
| COPA, COPB1, COPB2, COPG1, COPZ1, ARCN1 | established (COPI coatomer subunits) |
| NSF, NAPA (alpha-SNAP) | established (vesicle fusion machinery) |
| GOSR2, SCFD1, SYS1 | established (SNARE/Golgi trafficking) |
| HMGCR, HMGCS1 | established (mevalonate pathway; COPI-dependent SREBP) |
| TNPO1 | novel_role (nuclear import receptor) |
| PABPN1 | novel_role (poly(A) binding protein) |

**Why this is compelling**: Six of seven COPI coatomer subunits cluster together with SNARE fusion machinery and Golgi trafficking factors. The mevalonate pathway enzymes HMGCR/HMGCS1 make biological sense since COPI trafficking is essential for SREBP-mediated cholesterol homeostasis. The follow-up targets **TNPO1** (nuclear import receptor) and **PABPN1** (poly(A) binding nuclear protein) have no established vesicular trafficking roles, yet they phenocopy COPI loss. These could reveal novel connections between nuclear processes and secretory pathway morphology.

---

#### 8. BI199 -- Iron-Sulfur Cluster Biogenesis (Interphase)

**Process**: Iron-sulfur cluster biogenesis and cytosolic Fe-S protein assembly
**Genes**: 10 | **Established**: 6 | **Novel**: 4 | **Uncharacterized**: 0 | **%Established**: 60.0%

| Gene | Role |
|------|------|
| NFS1, ISCU | established (mitochondrial ISC assembly) |
| CIAO1, CIAO2B, CIAO3, MMS19 | established (CIA complex) |
| CHEK1 | novel_role (DDR kinase) |
| PPM1D | novel_role (DDR phosphatase) |
| PKMYT1 | novel_role (CDK1 inhibitory kinase) |
| ZC3H4 | novel_role (lncRNA termination factor) |

**Why this is compelling**: This cluster captures both arms of Fe-S cluster biogenesis -- mitochondrial ISC assembly (NFS1, ISCU) and the cytosolic CIA complex (CIAO1, CIAO2B, CIAO3, MMS19). The co-clustering of CHEK1 and PPM1D (DDR signaling) makes mechanistic sense since Fe-S clusters are essential cofactors for DNA repair helicases like XPD and FANCJ. The most intriguing follow-up target is **ZC3H4**, an lncRNA transcription termination factor with no known Fe-S connection, whose presence here could reveal a novel regulatory axis.

---

### Part II: Novel Connections -- Discovery Clusters

These clusters reveal unexpected biological connections, higher fractions of uncharacterized genes, or cross-pathway phenotypic convergence that suggests new biology.

---

#### 9. BI32 -- Actin Cytoskeleton Meets TCA Cycle (Interphase)

**Process**: CDC42/Rho GTPase-dependent actin organization and polarized membrane trafficking
**Genes**: 41 | **Established**: 17 | **Novel**: 22 | **Uncharacterized**: 2 | **%Established**: 41.5%

| Gene | Role |
|------|------|
| CDC42, PFN1, MRTFA, SRF, CALD1, DBNL, CDC42EP4, FLII | established (Rho GTPase/actin core) |
| EXOC1, EXOC4, EXOC5, EXOC7, EXOC8 | established (exocyst complex) |
| GGPS1, PGGT1B | established (prenylation) |
| OGDH, DLST, DLD, SUCLG1, SDHB, SDHD, MRPS36 | novel_role (TCA cycle/OGDHC) |
| DET1 | novel_role (JUN degradation) |
| PON1, BASP1 | novel_role (unexpected) |
| TMEM209 | novel_role (nuclear/membrane) |

**Why this is compelling**: This cluster reveals a striking co-clustering of the complete CDC42/Rho-actin signaling axis with the exocyst complex AND multiple TCA cycle components (OGDH/DLST/DLD = the oxoglutarate dehydrogenase complex; SDHB/SDHD = succinate dehydrogenase). The mechanistic hypothesis is that OGDHC-mediated histone succinylation in the nucleus could regulate SRF/MRTFA target gene expression, connecting mitochondrial metabolic flux to actin cytoskeletal gene programs. This represents a potential metabolic-epigenetic-cytoskeletal regulatory axis that has not been previously described.

---

#### 10. BI38 -- Focal Adhesion Meets NSL Complex (Interphase)

**Process**: Integrin-mediated focal adhesion signaling and cell-ECM adhesion
**Genes**: 39 | **Established**: 16 | **Novel**: 21 | **Uncharacterized**: 2 | **%Established**: 41.0%

| Gene | Role |
|------|------|
| ITGB1, ITGAV, ITGB5 | established (integrins) |
| ILK, LIMS1 | established (IPP complex) |
| PTK2/FAK, PXN, TLN1, FERMT2, BCAR1 | established (focal adhesion) |
| CRKL, RAPGEF1, RAC1, ELMO2, TNS3 | established (Rho/adaptor signaling) |
| KANSL1, KANSL2, KANSL3, KAT8, MCRS1 | novel_role (NSL HAT complex) |
| SEC61A1, SEC61G | novel_role (ER translocon) |
| DDX3X, ATF4, COP1 | novel_role (stress/ubiquitin) |

**Why this is compelling**: The co-clustering of four NSL histone acetyltransferase complex members (KANSL1/2/3, KAT8) with canonical focal adhesion components (integrins, FAK, paxillin, talin) is a novel finding. The NSL complex acetylates histone H4K16, and its known role in microtubule cytoskeleton organization in non-ciliated cells could intersect with focal adhesion dynamics at cortical microtubule capture sites. This suggests the NSL complex may transcriptionally regulate integrin/focal adhesion gene expression programs -- a connection not previously established.

---

## Mitotic Highlights

### Part I: Clean Biological Machines -- Validated Mitotic Clusters with Follow-Up Targets

Mitotic clusters are fewer in number but uniquely valuable -- they capture biology specific to cell division (spindle assembly, kinetochore function, cytokinesis) that would not appear in interphase screens.

---

#### 1. BM3 -- Complete Kinetochore and Spindle Machinery (Mitotic)

**Process**: Kinetochore assembly, spindle-microtubule attachment, chromosome segregation
**Genes**: 86 | **Established**: 50 | **Novel**: 28 | **Uncharacterized**: 8 | **%Established**: 58.1%

| Gene | Role |
|------|------|
| CENPW, CENPC, CENPH, CENPM, CENPN, CENPT, CENPI, CENPK, CENPL | established (inner kinetochore) |
| NDC80, NUF2, SPC24, SPC25 | established (NDC80 complex) |
| SKA1, SKA2, SKA3 | established (SKA complex) |
| DSN1, NSL1 | established (MIS12 complex) |
| ANAPC1/2/4/5/10/11, CDC16, CDC20, CDC23, CDC26 | established (APC/C) |
| DYNC1H1, DYNC1I2, DYNC1LI1, DYNLRB1 | established (dynein motor) |
| DCTN1/2/4/5, ACTR1A, ACTR10 | established (dynactin) |
| RAD21, CDCA5, SGO1, NAA50 | established (cohesin regulators) |
| CDK1, CCNB1, MAD2L1BP, PAFAH1B1, SPDL1 | established (mitotic kinases/checkpoint) |
| KIAA1211L | uncharacterized (no known function) |
| RGPD2, RGPD5, RGPD8 | uncharacterized (RANBP2-like) |
| SSBP4 | uncharacterized |
| CD8B, CLEC3B, PLGLB2, ARID5A | novel_role (unexpected) |
| COPS4, TRIM43 | novel_role (plausible connections) |

**Why this is compelling**: This is the most complete mitotic kinetochore/spindle cluster in the dataset. It contains nearly the full inner kinetochore (CENP-C/H/I/K/L/M/N/T/W), the complete NDC80 and SKA complexes, the APC/C ubiquitin ligase, and the dynein-dynactin motor -- the core machinery for chromosome segregation. Among the uncharacterized genes, three RGPD family members (RGPD2/5/8) are RANBP2-like proteins that co-cluster, suggesting a previously unrecognized family-wide mitotic role. **KIAA1211L** is completely uncharacterized and becomes a high-confidence mitotic gene candidate by virtue of sitting among this machinery.

---

#### 2. BM12 -- Tubulin Folding Pipeline (Mitotic)

**Process**: Tubulin folding, microtubule biogenesis, and mitotic spindle assembly
**Genes**: 33 | **Established**: 20 | **Novel**: 9 | **Uncharacterized**: 4 | **%Established**: 60.6%

| Gene | Role |
|------|------|
| TUBA1A, TUBA1B, TUBA1C | established (alpha-tubulins) |
| TUBB, TUBB2A, TUBB4B | established (beta-tubulins) |
| TCP1, CCT2-8 | established (all 8 TRiC/CCT chaperonin subunits) |
| PFDN1 | established (prefoldin co-chaperone) |
| TBCC, TBCD, TBCE | established (tubulin cofactors) |
| ARL2, DLGAP5 | established (tubulin cofactor regulator, spindle protein) |
| CHCHD5 | uncharacterized |
| ELOA3D, ELOA3B | uncharacterized |
| C19orf53 | uncharacterized |
| CHKA | novel_role (phospholipid biosynthesis) |
| ZNF24 | novel_role (transcription factor) |

**Why this is compelling**: This cluster captures the entire tubulin production pipeline: raw polypeptides are delivered by prefoldin (PFDN1) to the TRiC/CCT chaperonin (all 8 subunits), then passed to tubulin-specific folding cofactors (TBCC/D/E) and the ARL2 GTPase for final alpha-beta heterodimer assembly. The mitotic context makes this cluster especially meaningful, as the massive demand for tubulin during spindle assembly amplifies phenotypic effects. **CHCHD5** is completely uncharacterized and becomes a compelling candidate for involvement in the tubulin folding pathway.

---

#### 3. BM6 -- Mitotic Spindle Assembly (Mitotic)

**Process**: Mitotic spindle assembly and chromosome segregation
**Genes**: 40 | **Established**: 21 | **Novel**: 15 | **Uncharacterized**: 4 | **%Established**: 52.5%

| Gene | Role |
|------|------|
| HAUS1-8 | established (complete augmin/HAUS complex) |
| TUBG1, TUBGCP3, TUBGCP6, MZT1 | established (gamma-TuRC) |
| TACC3, CLTC | established (inter-MT bridge) |
| NCAPD2, NCAPG, NCAPH, SMC2, SMC4 | established (condensin I) |
| LENG8 | uncharacterized |
| C19orf84, POTEH | uncharacterized |
| COMMD8 | novel_role (endosomal recycling) |
| ATP5F1A | novel_role (mitochondrial ATP synthase) |
| FAM32A, TBC1D3G, ERH | novel_role |

**Why this is compelling**: This cluster contains the complete augmin/HAUS complex (all 8 subunits), which is the primary machinery for branched microtubule nucleation during spindle assembly. It also captures gamma-TuRC components, the TACC3-clathrin inter-microtubule bridge, and all five condensin I subunits. The most surprising member is **ATP5F1A** (mitochondrial ATP synthase alpha subunit), whose co-clustering with spindle assembly genes in mitotic cells suggests a connection between mitochondrial ATP production and spindle function that warrants investigation. **LENG8** is completely uncharacterized and represents a strong candidate for a novel mitotic spindle role.

---

### Part II: Novel Connections -- Discovery Mitotic Clusters

---

#### 4. BM5 -- Ubiquitin-Proteasome System in Mitosis (Mitotic)

**Process**: Ubiquitin-proteasome system and protein homeostasis in mitotic progression
**Genes**: 55 | **Established**: 22 | **Novel**: 31 | **Uncharacterized**: 2 | **%Established**: 40.0%

| Gene | Role |
|------|------|
| PSMC1-6, PSMD2/4/13 | established (proteasome subunits) |
| PSMA2/5, PSMB2/7, POMP | established (20S core) |
| RBX1, CUL3, ARIH1, NEDD8, UBA3, KCTD10 | established (CRL machinery) |
| USP9X, NPLOC4 | established (deubiquitinase, VCP cofactor) |
| BTBD9 | uncharacterized (BTB-domain protein) |
| DDN | uncharacterized (dendrin) |
| CS | novel_role (citrate synthase -- metabolic-proteasomal crosstalk) |
| CNPY2 | novel_role (UPS/myosin regulation) |
| ARHGAP11A | novel_role (RhoGAP with emerging mitotic function) |
| CFL1, CAPZB, TUBGCP5 | novel_role (cytoskeletal regulators) |
| RAD51D, GINS1 | novel_role (DDR/replication) |

**Why this is compelling**: This mitotic cluster reveals how the ubiquitin-proteasome system operates during cell division. It contains proteasome subunits alongside cullin-RING E3 ligase machinery (RBX1, CUL3, NEDD8), showing that both degradation machinery and its substrate-recognition apparatus must be intact for mitotic progression. The most intriguing discovery is **BTBD9**, a BTB-domain protein with virtually no functional annotation. BTB domains predict CUL3 substrate adaptor function, and its co-clustering with CUL3 strongly suggests BTBD9 is an uncharacterized CUL3 adaptor with mitotic substrates. **CS** (citrate synthase) is also unexpected, potentially revealing metabolic-proteasomal crosstalk during cell division.

---

#### 5. BM4 -- DNA Replication with Novel FAM25 Family (Mitotic)

**Process**: DNA replication initiation, fork progression, and replication-coupled repair
**Genes**: 72 | **Established**: 43 | **Novel**: 18 | **Uncharacterized**: 11 | **%Established**: 59.7%

| Gene | Role |
|------|------|
| ORC1, ORC6, CDC6, CDT1 | established (origin licensing) |
| CDC7, CDC45, TICRR, MCM7, MCM10, GINS3 | established (firing) |
| POLA1, POLA2, PRIM1, POLE, POLE2, POLD3 | established (polymerases) |
| PCNA, RFC1, CHTF18 | established (clamp/loading) |
| RPA1, RPA2, RPA3 | established (ssDNA binding) |
| TIMELESS, TIPIN, CLSPN | established (fork protection) |
| CHAF1A, CHAF1B | established (chromatin assembly) |
| MMS22L, TONSL, TRAIP, BARD1, RTEL1, REV3L | established (repair) |
| FAM25A, FAM25C, FAM25G | uncharacterized (FAM25 family) |
| ANTXRL | uncharacterized |
| HSF1 | novel_role (heat shock master regulator) |
| CNOT1, CNOT11 | novel_role (CCR4-NOT deadenylase) |

**Why this is compelling**: Beyond the comprehensive replication machinery, the most striking finding is the co-clustering of three **FAM25 family members** (FAM25A/C/G), which are completely uncharacterized proteins with no known function. Their shared presence in a replication-focused mitotic cluster suggests this entire gene family may have an unrecognized role in DNA replication or replication-coupled processes. **HSF1** (the master heat shock transcription factor) appearing among replication factors is also surprising and suggests a direct replication role beyond its canonical stress response function.

---

## Key Takeaways

### Screen Validation
The Brieflow Aconcagua OPS screen demonstrates extraordinary resolving power. It independently reconstructs:
- Complete 60S ribosome (39/47 proteins in one cluster)
- Separate 20S and 19S proteasome subcomplexes
- The spliceosome with >95% purity
- The full DNA replication fork machinery
- Complete augmin/HAUS complex (all 8 subunits)
- The entire tubulin folding pipeline

### Top Follow-Up Targets from Clean Validated Clusters

These genes are embedded among near-complete molecular machines, making their functional predictions exceptionally high-confidence:

| Gene | Predicted Function | Cluster Context | Condition |
|------|-------------------|-----------------|-----------|
| AAMP | Ribosome biogenesis or translational regulation | 39 of 42 genes are 60S RPL proteins | Interphase |
| TIMM29 | Novel role in DNA replication | 46 established replication factors | Interphase |
| ACTR5 | mRNP biogenesis/nuclear export | 10 of 11 genes are TREX/export factors | Interphase |
| ISY1-RAB43 | Spliceosome function | 20 of 21 genes are core spliceosome | Interphase |
| BTBD9 | CUL3 substrate adaptor (mitotic) | BTB domain + CUL3 co-clustering | Mitotic |
| KIAA1211L | Kinetochore/spindle function | Among complete NDC80/SKA/APC-C | Mitotic |
| RGPD2/5/8 | Family-wide mitotic role | Three RANBP2-like proteins co-cluster | Mitotic |
| FAM25A/C/G | DNA replication-associated | Three family members + replication fork | Mitotic |
| CHCHD5 | Tubulin folding | Among complete TRiC + tubulin cofactors | Mitotic |
| ZC3H4 | Fe-S cluster biogenesis link | Among ISC + CIA pathway | Interphase |
| LENG8 | Mitotic spindle assembly | Among HAUS/gamma-TuRC/condensin | Mitotic |

### Novel Biological Connections

| Connection | Evidence | Clusters |
|-----------|---------|----------|
| TCA cycle --> Actin cytoskeleton | OGDHC enzymes cluster with CDC42/SRF/exocyst | BI32 |
| NSL HAT complex --> Focal adhesion | KANSL1/2/3 + KAT8 cluster with integrins/FAK | BI38 |
| Metabolic-proteasomal crosstalk in mitosis | Citrate synthase clusters with UPS machinery | BM5 |
| Heat shock --> DNA replication | HSF1 among replication factors | BM4 |
| Nuclear import --> Golgi transport | TNPO1 clusters with COPI coatomer | BI146 |

---

## Genes Now Established: Screen Predictions Validated by Recent Literature

MozzareLLM's gene classifications reflect the LLM's training data cutoff. Several genes labeled "uncharacterized" or "novel_role" have since been independently validated by publications from 2020-2025, confirming the screen's predictive power.

| Gene | Cluster | MozzareLLM Label | Recent Evidence |
|------|---------|-----------------|-----------------|
| NEPRO | BI8 (SSU processome) | uncharacterized | 2025 NSMB: RNase MRP subunit essential for 40S ribosome biogenesis |
| MAK16 | BI3 (60S biogenesis) | uncharacterized | 2025 PNAS: [4Fe-4S]-dependent 60S assembly factor |
| NOC4L | BI8 (SSU processome) | uncharacterized | 2021-2023 cryo-EM confirms in human SSU processome |
| NOL10 | BI8 (SSU processome) | uncharacterized | 2022 Nat Commun: 40S biogenesis + nucleolar scaffold |
| ZC3H18 | BI44 (snRNA biogenesis) | uncharacterized | 2023-2025: CBC-NEXT adaptor for snRNA processing |
| BTBD9 | BI58 (CRL/neddylation) | uncharacterized | 2020: Confirmed CUL3 substrate adaptor |
| MTBP | BI15 (DNA replication) | novel_role | 2022-2025: Established replication origin firing factor (Treslin-MTBP) |
| RBM42 | BI13 (spliceosome) | uncharacterized | 2023 Nat Commun: Confirmed splicing regulator |
| ECD | BI13 (spliceosome) | novel_role | 2025 bioRxiv: Direct U5 snRNA binding confirmed |
| ADNP | BI162 (cohesin) | novel_role | 2023: ADNP deficiency alters cohesin-insulated regions |
| ORC3 | BI109 (20S proteasome) | novel_role | 2024: Proteasome regulates ORC protein chromatin dynamics |

These validations demonstrate that when the screen places a gene in a specific pathway context, subsequent independent experiments frequently confirm the prediction. The true fraction of established genes in high-confidence clusters is higher than MozzareLLM reports.

---

*Report generated from MozzareLLM analysis of Aconcagua OPS screen (claude-opus-4-6). All clusters shown are High confidence unless otherwise noted. Gene role classifications (established/novel_role/uncharacterized) are assigned by MozzareLLM based on literature evidence relative to the cluster's dominant biological process.*
