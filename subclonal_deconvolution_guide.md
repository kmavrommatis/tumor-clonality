# Subclonal deconvolution in bulk NGS data: tools, statistics, and feasibility

This guide covers tools for subclonal deconvolution of tumors from bulk whole exome sequencing (WES) and whole genome sequencing (WGS), with attention to how each tool handles copy-number alterations (CNAs). It estimates how sequencing depth and tumor purity limit subclone detection, and describes four problems that affect real analyses: copy-number effects, errors in upstream purity and copy-number calls, the neutral tail of low-frequency mutations, and the trade-off between deeper sequencing and more samples.

## 1. Subclonal deconvolution tools and how they handle CNAs

Restricting an analysis to diploid regions discards much of the genome, so the table records how each tool treats copy number. Tools are listed by year of publication, followed by CliPP (a widely used preprint) and two benchmark studies. In the CNA handling column, *Yes* means the tool models copy number when it estimates cancer cell fractions (CCFs) or builds trees; *CN-centric* means it infers subclones from copy-number profiles rather than from SNVs; *No* means it should be run only on copy-neutral regions; and *Not stated* means the original abstract does not say. Every DOI was checked against its PubMed record, except Cloe (*Ann Appl Stat*, not indexed in PubMed), whose DOI comes from the publisher. BayClone appeared in conference proceedings without a DOI, so it links to PubMed instead.

| # | Tool | Publication | CNA handling | Algorithm & Requirements Summary |
| :--- | :--- | :--- | :--- | :--- |
| 1 | **ABSOLUTE** | Carter et al., *Nat Biotechnol* 2012 [DOI](https://doi.org/10.1038/nbt.2203) | Yes | Estimates purity and ploidy from segmented CN and SNV allele fractions, then converts each mutation to a CCF and calls it clonal or subclonal; needs CN segments and SNV counts. |
| 2 | **DPClust** (+ Battenberg) | Nik-Zainal et al., *Cell* 2012 [DOI](https://doi.org/10.1016/j.cell.2012.04.023); Bolli et al., *Nat Commun* 2014 [DOI](https://doi.org/10.1038/ncomms3997) | Yes | Dirichlet-process clustering of CN- and purity-corrected CCFs, single or multi-sample; needs SNV counts plus subclonal CN and purity, typically from WGS. |
| 3 | **THetA / THetA2** | Oesper et al., *Genome Biol* 2013 [DOI](https://doi.org/10.1186/gb-2013-14-7-r80); Oesper et al., *Bioinformatics* 2014 [DOI](https://doi.org/10.1093/bioinformatics/btu651) | CN-centric | Infers the mixture of tumor populations and their CN profiles from read depth, with THetA2 adding B-allele frequency (BAF) and WES support; needs tumor/normal depth and BAF. |
| 4 | **TrAp** | Strino et al., *Nucleic Acids Res* 2013 [DOI](https://doi.org/10.1093/nar/gkt641) | No | Solves for all subclonal trees consistent with aggregate mutation frequencies under sum constraints; needs cellular frequencies from one sample. |
| 5 | **PyClone** | Roth et al., *Nat Methods* 2014 [DOI](https://doi.org/10.1038/nmeth.2883) | Yes | Dirichlet-process beta-binomial MCMC that integrates over possible mutation genotypes given major/minor CN; needs ref/alt counts, allele-specific CN and purity, and supports multiple samples. |
| 6 | **SciClone** | Miller et al., *PLoS Comput Biol* 2014 [DOI](https://doi.org/10.1371/journal.pcbi.1003665) | No | Variational Bayesian beta/binomial mixture on VAFs from one or more samples; uses CN segments only to restrict analysis to copy-neutral, LOH-free regions. |
| 7 | **PhyloSub** | Jiao et al., *BMC Bioinformatics* 2014 [DOI](https://doi.org/10.1186/1471-2105-15-35) | No | Tree-structured stick-breaking prior with MCMC on VAF-derived frequencies (PhyloWGS's predecessor); assumes diploid loci. |
| 8 | **Clomial** | Zare et al., *PLoS Comput Biol* 2014 [DOI](https://doi.org/10.1371/journal.pcbi.1003703) | No | Binomial latent-clone model fit by EM to multi-sample allele counts; assumes CN-neutral loci. |
| 9 | **EXPANDS** | Andor et al., *Bioinformatics* 2014 [DOI](https://doi.org/10.1093/bioinformatics/btt622) | Yes | Builds a cell-frequency probability distribution for each mutation from allele frequency and local CN, then clusters these into subpopulations; needs SNV frequencies plus CN segments. |
| 10 | **CHAT** | Li & Li, *Genome Biol* 2014 [DOI](https://doi.org/10.1186/s13059-014-0473-4) | Yes | Estimates segment-level aneuploid fractions from depth and BAF, then adjusts SNV allele frequencies to CCFs; needs tumor/normal depth and BAF. |
| 11 | **cloneHD** | Fischer et al., *Cell Rep* 2014 [DOI](https://doi.org/10.1016/j.celrep.2014.04.055) | Yes | Hidden Markov model jointly inferring subclonal CN, BAF and SNV genotypes across samples; needs read-depth, BAF and SNV tracks, best suited to WGS. |
| 12 | **SubcloneSeeker** | Qiao et al., *Genome Biol* 2014 [DOI](https://doi.org/10.1186/s13059-014-0443-x) | Yes | Exhaustively enumerates subclone structures consistent with the prevalences of SNV and CNV clusters; needs pre-clustered prevalence estimates. |
| 13 | **TITAN** | Ha et al., *Genome Res* 2014 [DOI](https://doi.org/10.1101/gr.180281.114) | CN-centric | HMM segmenting subclonal CN and LOH and estimating the cellular prevalence of each event; needs tumor/normal WGS or WES depth and germline allelic ratios. |
| 14 | **PhyloWGS** | Deshwar et al., *Genome Biol* 2015 [DOI](https://doi.org/10.1186/s13059-015-0602-8) | Yes | Tree-structured stick-breaking MCMC that places SNVs and CNA segments jointly on a phylogeny with a CN correction to VAFs; needs SNV counts plus subclonal CN calls. |
| 15 | **LICHeE** | Popic et al., *Genome Biol* 2015 [DOI](https://doi.org/10.1186/s13059-015-0647-8) | No | Builds lineage trees from multi-sample VAF clusters under ancestry constraints; requires CN-neutral SNVs. |
| 16 | **CITUP** | Malikic et al., *Bioinformatics* 2015 [DOI](https://doi.org/10.1093/bioinformatics/btv003) | No | Combinatorial or quadratic-integer-programming search over tree topologies fitting multi-sample frequencies; diploid assumption. |
| 17 | **AncesTree** | El-Kebir et al., *Bioinformatics* 2015 [DOI](https://doi.org/10.1093/bioinformatics/btv261) | No | Integer linear program (ILP) for the VAF factorization problem under a perfect phylogeny; needs multi-sample VAFs at diploid loci. |
| 18 | **SCHISM** | Niknafs et al., *PLoS Comput Biol* 2015 [DOI](https://doi.org/10.1371/journal.pcbi.1004416) | Yes (via CCF input) | Hypothesis tests of lineage precedence between clusters, followed by a genetic-algorithm tree search; needs CN- and purity-adjusted CCFs. |
| 19 | **BayClone** | Sengupta et al., *Pac Symp Biocomput* 2015 (no DOI; [PubMed](https://pubmed.ncbi.nlm.nih.gov/25592605/)) | No | Bayesian feature-allocation model assigning SNVs to latent subclones; uses read counts at loci assumed diploid. |
| 20 | **Canopy** | Jiang et al., *PNAS* 2016 [DOI](https://doi.org/10.1073/pnas.1522203113) | Yes | Bayesian MCMC phylogeny jointly using SNA allele counts and allele-specific CNA segment ratios; needs single- or multi-sample SNV counts and allele-specific CN. |
| 21 | **SPRUCE** | El-Kebir et al., *Cell Syst* 2016 [DOI](https://doi.org/10.1016/j.cels.2016.07.004) | Yes | Enumerates multi-state perfect phylogenies in which SNVs and CNAs co-evolve; needs multi-sample VAFs and allele-specific CN. |
| 22 | **BubbleTree** | Zhu et al., *Nucleic Acids Res* 2016 [DOI](https://doi.org/10.1093/nar/gkv1102) | CN-centric | Jointly estimates purity, ploidy, allele-specific CN and clonality from a graph of depth ratio against BAF; needs tumor/normal depth, BAF and SNV calls. |
| 23 | **Cloe** | Marass et al., *Ann Appl Stat* 2016 [DOI](https://doi.org/10.1214/16-AOAS986) | No | Phylogenetic latent-feature model fit by MCMC to multi-sample VAFs; assumes CN-neutral loci. |
| 24 | **Treeomics** | Reiter et al., *Nat Commun* 2017 [DOI](https://doi.org/10.1038/ncomms14114) | No | Bayesian inference of whether each mutation is present per sample, followed by an ILP phylogeny for multi-region data; tolerant of CN effects but does not model them. |
| 25 | **ClonEvol** | Dang et al., *Ann Oncol* 2017 [DOI](https://doi.org/10.1093/annonc/mdx517) | Inherited | Bootstrap-based clonal ordering and visualization of pre-computed clusters; CN handling depends on the upstream clustering tool. |
| 26 | **CTPsingle** | Donmez et al., *J Comput Biol* 2017 [DOI](https://doi.org/10.1089/cmb.2016.0148) | Not stated | Infers purity and clonality from a single low-coverage sample (around 30x); needs SNV counts. |
| 27 | **QuantumClone** | Deveau et al., *Bioinformatics* 2018 [DOI](https://doi.org/10.1093/bioinformatics/bty016) | Yes | EM clustering correcting for genotype and CN, with BIC model selection, for joint multi-sample analysis; needs SNV counts, CN genotypes and purity. |
| 28 | **Sclust** | Cun et al., *Nat Protoc* 2018 [DOI](https://doi.org/10.1038/nprot.2018.033) | Yes | Infers purity, ploidy and subclonal CN, then clusters CN-corrected CCFs; needs tumor/normal BAMs (WGS preferred). |
| 29 | **CloneFinder** | Miura et al., *Bioinformatics* 2018 [DOI](https://doi.org/10.1093/bioinformatics/bty469) | No | Phylogeny-aware iterative decomposition of multi-region SNV frequencies into clone genotypes; best used on CN-neutral SNVs. |
| 30 | **TUSV** | Eaton et al., *Bioinformatics* 2018 [DOI](https://doi.org/10.1093/bioinformatics/bty270) | Yes | ILP deconvolution into a phylogeny that includes SNVs, CNAs and structural variants (SVs); needs multi-sample WGS SNV, CN and SV calls. |
| 31 | **CloneDeMix** | Tai et al., *PLoS One* 2018 [DOI](https://doi.org/10.1371/journal.pone.0206579) | CN-centric | Two-way mixture Poisson model on read depth inferring subclonal CN, mutational cellular prevalence and mutation order; needs tumor/normal depth. |
| 32 | **BAMSE** | Toosi et al., *BMC Bioinformatics* 2019 [DOI](https://doi.org/10.1186/s12859-019-2824-3) | No | Bayesian model selection over the number of subclones and tree topology from read counts; mainly assumes diploid loci. |
| 33 | **CALDER** | Myers et al., *Cell Syst* 2019 [DOI](https://doi.org/10.1016/j.cels.2019.05.010) | Yes (via input) | Mixed-ILP phylogeny that respects the time order of longitudinal samples; needs CN-adjusted counts from serial samples. |
| 34 | **GenoClone** | Zou, Jin & Au, *Brief Bioinform* 2019 [DOI](https://doi.org/10.1093/bib/bby084) | Not stated | Combines VAF with linkage between somatic and germline variants from paired-end or long reads to infer subclones; needs reads spanning both. |
| 35 | **B-SCITE** | Malikic et al., *Nat Commun* 2019 [DOI](https://doi.org/10.1038/s41467-019-10737-5) | Partial | Infers tumor phylogenies from combined single-cell and bulk sequencing data and is reported to stay robust when bulk VAFs are distorted by CN; needs both bulk and single-cell data. |
| 36 | **BaSiC** | Sun et al., *Biometrics* 2020 [DOI](https://doi.org/10.1111/biom.13198) | Not stated | Jointly models bulk and single-cell DNA-seq to infer intratumor heterogeneity; needs both data types. |
| 37 | **PyClone-VI** | Gillis & Roth, *BMC Bioinformatics* 2020 [DOI](https://doi.org/10.1186/s12859-020-03919-2) | Yes | Variational-inference version of PyClone's CN-aware mixture that scales to WGS-sized mutation sets; needs ref/alt counts, major/minor CN and purity. |
| 38 | **MOBSTER** | Caravagna et al., *Nat Genet* 2020 [DOI](https://doi.org/10.1038/s41588-020-0675-5); R package: *BMC Bioinformatics* 2020 [DOI](https://doi.org/10.1186/s12859-020-03863-1) | Yes (per karyotype) | Beta-mixture plus Pareto-tail model that accounts for neutral within-clone mutations, fitted separately for each karyotype; needs high-depth WGS counts plus allele-specific CN and purity. |
| 39 | **FastClone** | Xiao et al., *Nat Commun* 2020 [DOI](https://doi.org/10.1038/s41467-020-18169-2) | Yes | Fast probabilistic deconvolution assigning SNVs to subclones while modeling which CN allele they lie on; single-sample only, and needs SNV counts plus CN and purity. |
| 40 | **HATCHet** | Zaccaria & Raphael, *Nat Commun* 2020 [DOI](https://doi.org/10.1038/s41467-020-17967-y) | CN-centric | Jointly infers allele- and clone-specific CN across samples by clustering depth ratio and BAF; needs multi-sample WGS (or WES) BAMs. |
| 41 | **SVclone** | Cmero et al., *Nat Commun* 2020 [DOI](https://doi.org/10.1038/s41467-020-14351-8) | Yes | Estimates the CCF and copy number of SV breakpoints from both breakends; needs WGS BAMs, SV calls and CN. |
| 42 | **TrackSig** | Rubanova et al., *Nat Commun* 2020 [DOI](https://doi.org/10.1038/s41467-020-14352-7) | Yes | Orders mutations by CN-corrected CCF and detects changes in mutational-signature activity along that order; needs WGS SNVs with CN and purity. |
| 43 | **msphy-SCNAClonal** | Chu, Nie & Wang, *Front Genet* 2020 [DOI](https://doi.org/10.3389/fgene.2019.01374) | CN-centric | Subclonal reconstruction of somatic CN alterations with read-count-ratio bias correction and false-breakpoint filtering; needs tumor/normal NGS data. |
| 44 | **DeCiFer** | Satas et al., *Cell Syst* 2021 [DOI](https://doi.org/10.1016/j.cels.2021.07.006) | Yes | Clusters SNVs by descendant cell fraction, explicitly modeling SNVs gained or lost through CNAs, jointly across samples; needs multi-sample counts and allele-specific CN. |
| 45 | **CloneSig** | Abécassis et al., *Nat Commun* 2021 [DOI](https://doi.org/10.1038/s41467-021-24992-y) | Yes | Joint model of subclonal clustering and mutational signatures, with signature exposures specific to each clone; needs SNVs with sequence context plus CN and purity. |
| 46 | **SubMARine** | Sundermann et al., *PLoS Comput Biol* 2021 [DOI](https://doi.org/10.1371/journal.pcbi.1008400) | Yes (extension) | Polynomial-time partial clone tree summarizing all valid trees, with an extension for subclonal CNAs; needs pre-clustered subclone frequencies. |
| 47 | **CRUST** | Chattopadhyay et al., *Brief Bioinform* 2021 [DOI](https://doi.org/10.1093/bib/bbab292) | Yes | Semi-supervised multi-sample deconvolution with purity normalization and CN adjustment, able to estimate CN itself if none is supplied; works at about 50x. |
| 48 | **DEVOLUTION** | Andersson et al., *Commun Biol* 2021 [DOI](https://doi.org/10.1038/s42003-021-02637-6) | Yes | Subclonal deconvolution followed by phylogeny reconstruction that integrates copy number and sequencing information across tumor regions, provided each mutation's mutated clone fraction is known; accepts SNP-array, WES, WGS or targeted data. |
| 49 | **Unnamed AML method** | Stiehl & Marciniak-Czochra, *Front Physiol* 2021 [DOI](https://doi.org/10.3389/fphys.2021.596194) | No | Enumerates all clonal hierarchies compatible with the bulk VAFs of one sample; needs VAFs only. |
| 50 | **HDP bulk+single-cell model** | He et al., *Ann Appl Stat* 2021 [DOI](https://doi.org/10.1214/20-aoas1434) | Not stated | Hierarchical Dirichlet process mixture with Gibbs sampling for structured bulk and single-cell sampling designs; needs both data types. |
| 51 | **Pairtree** | Wintersinger et al., *Blood Cancer Discov* 2022 [DOI](https://doi.org/10.1158/2643-3230.BCD-21-0092) | Yes (via ω) | Computes pairwise relations between mutation clusters, then samples trees by MCMC, for many samples; needs counts plus a CN-derived variant read probability (ω) per mutation. |
| 52 | **PICTograph** | Zheng et al., *Bioinformatics* 2022 [DOI](https://doi.org/10.1093/bioinformatics/btac367) | Yes | Bayesian (JAGS) model integrating over tree space using CN-adjusted CCFs; needs SNV counts, CN and purity. |
| 53 | **MEDICC2** | Kaufmann et al., *Genome Biol* 2022 [DOI](https://doi.org/10.1186/s13059-022-02794-9) | CN-centric | Minimum-event distance phylogenies from allele-specific CN profiles, including whole-genome doubling; needs allele-specific CN segments. |
| 54 | **MAGOS** | Ahmadinejad et al., *Mol Biol Evol* 2022 [DOI](https://doi.org/10.1093/molbev/msac136) | Yes | Clustering error model accounting for how VAF variance depends on depth and mean VAF, working down to about 30x; accepts SNVs and CNVs from one or more samples. |
| 55 | **cfdna-wgs** | Herberts et al., *Nature* 2022 [DOI](https://doi.org/10.1038/s41586-022-04975-9) | Yes | K-means-based clustering of CN-corrected mutation fractions from deep ctDNA WGS; needs SNV counts, CN and ctDNA fraction. |
| 56 | **CONIPHER** | Grigoriadis et al., *Nat Protoc* 2024 [DOI](https://doi.org/10.1038/s41596-023-00913-9) | Yes | TRACERx framework combining CN-aware clustering, error correction (including mutation loss) and multi-sample tree building; needs multi-region counts, allele-specific CN and purity. |
| 57 | **Orchard** | Kulman et al., *PLoS Comput Biol* 2024 [DOI](https://doi.org/10.1371/journal.pcbi.1012653) | Yes (via ω) | Stochastic combinatorial search for large phylogenies, the successor to Pairtree; same inputs as Pairtree. |
| 58 | **RETCHER** | Wang et al., *Brief Bioinform* 2024 [DOI](https://doi.org/10.1093/bib/bbae516) | Yes | Identifies the number of mutated copies of each SNV under CN changes to correct CCFs, assesses how credible each cluster is, and infers phylogenies; needs multi-sample SNV and CN data. |
| 59 | **PhyClone** | Hurtado, Bouchard-Côté & Roth, *Bioinformatics* 2025 [DOI](https://doi.org/10.1093/bioinformatics/btaf344) | Yes | Bayesian sampler over phylogenies (forest-structured Chinese restaurant process with particle Gibbs) that models mutation loss and outliers; needs allele-specific CN and counts from one or more samples. |
| 60 | **MyClone** | Guo et al., *BMC Bioinformatics* 2025 [DOI](https://doi.org/10.1186/s12859-025-06252-8) | Yes (heuristic) | Variational-EM clustering of copy-neutral mutations, then reassignment of CNA-affected mutations to those clones by matching possible genotypes; designed for deep (around 2,000x) targeted or ctDNA data and weaker at WES/WGS depth. |
| 61 | **Sapling** | Qi, Schmidt & El-Kebir, *Algorithms Mol Biol* 2026 [DOI](https://doi.org/10.1186/s13015-025-00295-5) | Not stated | Infers a small set of backbone trees that summarize the space of plausible phylogenies from bulk read counts, and can also infer full trees. |
| 62 | **CliPP** *(preprint, not peer-reviewed)* | Jiang et al., *bioRxiv* 2021, updated 2026 [DOI](https://doi.org/10.1101/2021.03.31.437383) | Yes | Penalized-likelihood (pairwise fusion) clustering of cellular prevalences that needs no extensive priors; takes WGS/WES variant calls with read depth, allele counts and copy-number segments. |
| 63 | *Benchmark:* **Tanner et al.** | *Nat Commun* 2021 [DOI](https://doi.org/10.1038/s41467-021-26698-7) | n/a | Benchmarked complete pipelines on 80 simulated bulk WES datasets; Mutect2, then FACETS, then PyClone-VI was the best combination, and purity and purity-corrected depth drove accuracy. |
| 64 | *Benchmark:* **ICGC-TCGA DREAM Challenge (Salcedo et al.)** | *Nat Biotechnol* 2024 [DOI](https://doi.org/10.1038/s41587-024-02250-y) | n/a | Benchmarked 31 single-sample subclonal reconstruction algorithms on 51 simulated tumors across seven tasks; no single method won all tasks, and CN state, purity-adjusted depth and mappability affected performance. |

## 2. Feasibility of subclone detection

Whether a subclone can be detected depends on tumor purity, the fraction of tumor cells in the subclone, sequencing depth, and how many mutations the subclone carries. For a single sample, the tables estimate whether a subclone separates from the clonal (trunk) mutations as its own cluster, and roughly how many mutations it needs.

All values come from the model in Section 3, with these definitions:

* Expected VAF, for a heterozygous mutation in a diploid region, is $VAF = Purity \times CCF / 2$.
* Separation score ($S$) is the distance between the trunk and subclone VAFs divided by their combined binomial standard deviation at the stated depth (Section 3A).
* Retained is the fraction of the subclone's true mutations expected to pass a caller filter of $\ge 4$ alternate reads, calculated from the binomial distribution.
* Feasible is *No* if $S < 1.5$ or fewer than 20% of the subclone's mutations are retained, *Marginal* if $1.5 \le S < 2$, and *Yes* if $S \ge 2$.
* Biological mutations needed is a rough count of the detected mutations a mixture model needs before it opens a separate cluster (10-15 for $S \ge 4$, 20-30 for $3 \le S < 4$, 40-60 for $2 \le S < 3$, and 100 or more for marginal cases), divided by the retained fraction and rounded to the nearest 5.

These are best-case estimates. They assume heterozygous mutations in diploid regions, the same depth at every site, known purity, pure binomial noise, and a single subclone alongside the trunk. Real data are harder. Exome coverage varies widely between sites, overdispersion widens VAF distributions (Section 3B), retained mutations are biased toward higher VAF (Section 3C), and copy-number changes, purity errors and the neutral tail can create or hide clusters (Section 4). The mutation counts are rules of thumb, not the result of a formal power calculation. To estimate power for a specific study design, simulate data and run it through the pipeline you plan to use.

### 100x (typical WES)

| Tumor Purity | Subclone Fraction | Expected VAFs (Trunk / Subclone) | Separation Score | Subclonal Mutations Retained (≥4 alt reads) | Feasible? | Biological Mutations Needed |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **10%** | **20%** | 5.0% / 1.0% | 1.67 | 2% | No | N/A |
| **10%** | **40%** | 5.0% / 2.0% | 1.16 | 14% | No | N/A |
| **10%** | **60%** | 5.0% / 3.0% | 0.72 | 35% | No | N/A |
| **20%** | **20%** | 10.0% / 2.0% | 2.42 | 14% | No | N/A |
| **20%** | **40%** | 10.0% / 4.0% | 1.67 | 57% | Marginal | >175 |
| **20%** | **60%** | 10.0% / 6.0% | 1.05 | 86% | No | N/A |
| **30%** | **20%** | 15.0% / 3.0% | 3.03 | 35% | Yes | ≈55-85 |
| **30%** | **40%** | 15.0% / 6.0% | 2.10 | 86% | Yes | ≈45-70 |
| **30%** | **60%** | 15.0% / 9.0% | 1.31 | 98% | No | N/A |
| **40%** | **20%** | 20.0% / 4.0% | 3.59 | 57% | Yes | ≈35-55 |
| **40%** | **40%** | 20.0% / 8.0% | 2.48 | 96% | Yes | ≈40-60 |
| **40%** | **60%** | 20.0% / 12.0% | 1.55 | 100% | Marginal | >100 |
| **50%** | **20%** | 25.0% / 5.0% | 4.13 | 74% | Yes | ≈15-20 |
| **50%** | **40%** | 25.0% / 10.0% | 2.85 | 99% | Yes | ≈40-60 |
| **50%** | **60%** | 25.0% / 15.0% | 1.78 | 100% | Marginal | >100 |

### 200x (deep WES)

| Tumor Purity | Subclone Fraction | Expected VAFs (Trunk / Subclone) | Separation Score | Subclonal Mutations Retained (≥4 alt reads) | Feasible? | Biological Mutations Needed |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **10%** | **20%** | 5.0% / 1.0% | 2.36 | 14% | No | N/A |
| **10%** | **40%** | 5.0% / 2.0% | 1.64 | 57% | Marginal | >175 |
| **10%** | **60%** | 5.0% / 3.0% | 1.02 | 85% | No | N/A |
| **20%** | **20%** | 10.0% / 2.0% | 3.42 | 57% | Yes | ≈35-55 |
| **20%** | **40%** | 10.0% / 4.0% | 2.37 | 96% | Yes | ≈40-60 |
| **20%** | **60%** | 10.0% / 6.0% | 1.48 | 100% | No | N/A |
| **30%** | **20%** | 15.0% / 3.0% | 4.29 | 85% | Yes | ≈10-20 |
| **30%** | **40%** | 15.0% / 6.0% | 2.97 | 100% | Yes | ≈40-60 |
| **30%** | **60%** | 15.0% / 9.0% | 1.85 | 100% | Marginal | >100 |
| **40%** | **20%** | 20.0% / 4.0% | 5.08 | 96% | Yes | ≈10-15 |
| **40%** | **40%** | 20.0% / 8.0% | 3.51 | 100% | Yes | ≈20-30 |
| **40%** | **60%** | 20.0% / 12.0% | 2.20 | 100% | Yes | ≈40-60 |
| **50%** | **20%** | 25.0% / 5.0% | 5.83 | 99% | Yes | ≈10-15 |
| **50%** | **40%** | 25.0% / 10.0% | 4.03 | 100% | Yes | ≈10-15 |
| **50%** | **60%** | 25.0% / 15.0% | 2.52 | 100% | Yes | ≈40-60 |

A subclone present in 60% of tumor cells is often harder to detect than one present in 40%, because its VAF sits close to the trunk peak. At low purity the limit shifts from peak overlap to variant-caller attrition: at 10% purity, most subclonal mutations never reach the clustering step. For WGS at 30-60x, separation scores are roughly $\sqrt{100/N}$ times those in the 100x table (about 0.7 times at 50x), so a single WGS sample can resolve only well-separated subclones in high-purity tumors. WGS partly makes up for this by detecting far more mutations.

## 3. Statistical model of subclone detection

The numbers in Section 2 depend on the binomial spread of each mutation's VAF, on sequencing noise beyond the binomial, and on the loss of low-VAF mutations before clustering. They also depend on how clustering models decide how many clusters to report.

### A. Binomial variance and peak separation ($\Delta$VAF)

Tools model the variant read count at each site ($v_i$) against its total depth ($N_i$). For a heterozygous mutation in a diploid region, the expected VAF ($p_i$) is:
$$p_i = \frac{Purity \times CCF_i}{2}$$
where $CCF_i$ is the cancer cell fraction, the fraction of tumor cells carrying the mutation. Section 4A gives the formula for other copy-number states.

Reads are a random sample, so the variant read count follows a binomial distribution, $v_i \sim \text{Binomial}(N_i, p_i)$, and the standard deviation ($\sigma$) of a single mutation's observed VAF is:
$$\sigma = \sqrt{\frac{p_i(1-p_i)}{N_i}}$$

The separation score divides the distance between the expected VAFs of two clusters by their combined spread:
$$S = \frac{|VAF_{trunk} - VAF_{subclone}|}{\sqrt{\sigma_{trunk}^2 + \sigma_{subclone}^2}}$$
Above about 4 the clusters separate cleanly, and 10-15 detected mutations are usually enough. Between 2 and 4 the distributions partly overlap and more mutations are needed. Below about 2 the combined histogram usually shows a single peak.

$N_i$ is the depth at each site, and exome coverage is uneven, so many sites are sequenced well below the mean and the real separation is lower than the score at mean depth suggests. In the benchmark by Tanner et al. (*Nat Commun* 2021, [DOI](https://doi.org/10.1038/s41467-021-26698-7)), purity-corrected depth was one of the strongest predictors of accuracy. A low score does not make separation impossible, though. With $n$ mutations in a cluster, the standard error of the cluster's mean VAF is roughly $\sigma/\sqrt{n}$, so a mixture model can still resolve overlapping clusters if there are enough mutations. It does mean that many individual mutations will be assigned to the wrong cluster, and that the number of mutations needed rises steeply.

### B. Beta-binomial sequencing noise

Real sequencing data vary more than the binomial model predicts, which is called overdispersion. Causes include PCR duplicates and amplification bias, differences in capture efficiency between baits, reference bias against the alternate allele during mapping, and local copy-number heterogeneity that the segmentation misses. Several tools, including PyClone and PyClone-VI, offer a beta-binomial likelihood with a precision (overdispersion) parameter to absorb it:
$$v_i \sim \text{Beta-Binomial}(N_i, p_i, \rho)$$
Overdispersion widens each cluster, lowers the effective separation score and raises the number of mutations needed. Fitting a pure binomial model to overdispersed data carries the opposite risk: it tends to split one real cluster into several spurious ones.

### C. Variant caller attrition and detection bias

The somatic caller decides which mutations reach the clustering step. Mutect2 has no hard minimum number of alternate reads by default (it uses a likelihood-based score), but many pipelines add a minimum of 3-5 alternate reads, and every caller has low sensitivity when only a few reads support a variant. With a cutoff of $\ge 4$ alternate reads at 100x depth, the binomial model predicts the loss of about 14% of true mutations at 6% VAF (a mean of 6 alternate reads), 26% at 5% VAF, 43% at 4% VAF and 86% at 2% VAF.

The surviving mutations are also biased. They are the ones that happened to receive extra alternate reads, so the observed VAFs of a subclone near the detection limit are too high. The subclone moves toward the trunk, the gap between them narrows, and the subclone's CCF is overestimated. In multi-sample studies, force-call (genotype) the union of all mutations in every sample and record read counts even where a site was not called, because a mutation that was not called in a sample may still be present in it.

### D. Model complexity and cluster number

Every clustering method penalizes complexity in some way when it decides how many clusters the data support. Dirichlet-process models (PyClone, DPClust, PhyloWGS) use a concentration parameter that controls how readily a new cluster opens. Variational methods such as PyClone-VI fit a mixture with a sparsity-inducing prior and drop empty components. Others use BIC (QuantumClone) or a penalized likelihood (CliPP). Under all of these schemes, a few mutations at an intermediate VAF are cheaper to explain as noise around an existing cluster than as a new population, so a subclone is reported only once it has enough detected mutations. Because these settings can be tuned, results can change with the priors, the concentration parameter and even the random seed. Report the settings, and check that the inferred structure is stable across runs.

## 4. Complications in real data

The feasibility analysis assumes diploid, heterozygous mutations and correct purity. Four further issues decide whether an inferred subclonal structure is real.

### A. Copy-number effects

Most solid tumors are aneuploid, and copy number changes how VAF relates to CCF. For a mutation present on $m$ copies (its multiplicity) at a locus with total tumor copy number $CN_t$, and a normal copy number of 2:
$$VAF = \frac{Purity \times m \times CCF}{Purity \times CN_t + 2(1 - Purity)}$$

The table gives the expected VAF of a clonal mutation (CCF = 1) at 50% purity in several copy-number states, and the CCF that a naive diploid conversion ($CCF = 2 \times VAF / Purity$) would report:

| Copy-number state at locus | $CN_t$ | Multiplicity $m$ | Expected VAF | Naive diploid CCF |
| :--- | :--- | :--- | :--- | :--- |
| Diploid, heterozygous | 2 | 1 | 25.0% | 1.00 |
| Copy-neutral LOH (mutation duplicated) | 2 | 2 | 50.0% | 2.00 |
| Single-copy deletion of the other allele | 1 | 1 | 33.3% | 1.33 |
| Single-copy gain, mutation not on gained allele | 3 | 1 | 20.0% | 0.80 |
| Single-copy gain, mutation on gained allele | 3 | 2 | 40.0% | 1.60 |
| After whole-genome doubling, mutation post-doubling | 4 | 1 | 16.7% | 0.67 |
| After whole-genome doubling, mutation pre-doubling | 4 | 2 | 33.3% | 1.33 |

Every row is a clonal mutation, but a diploid model would put the last six rows at five different wrong CCFs. In a genome-doubled tumor, mutations acquired after the doubling form a large, tight cluster at half the VAF of those acquired before it, and a diploid model reports them as a subclone at CCF 0.67. Genome doubling is common, so this artifact is too. CNA-aware tools are needed whenever a substantial fraction of mutations lies outside diploid regions.

Multiplicity is not observed, so CNA-aware tools have to infer $m$ for each mutation. PyClone integrates over the possible genotypes, while other tools pick the most likely multiplicity. When $CN_t > 2$, mutations with different multiplicities and CCFs can produce the same VAF, so some ambiguity remains even when the copy number is correct.

A CNA present in only some tumor cells makes the local copy number a mixture. A mutation's CCF then depends on whether it arose before or after the CNA, and in which lineage. Only some tools, such as PhyloWGS, Canopy, SPRUCE and DeCiFer, model this; others require subclonal CNA regions to be excluded.

A deletion can also remove the allele carrying a mutation in a descendant lineage, so a clonal mutation looks subclonal or disappears from one region of a multi-region study. This breaks the infinite-sites assumption behind most tree-building methods. CONIPHER, PhyClone, DeCiFer and SPRUCE allow for mutation loss. With other tools, treat mutations in regions with subclonal losses with caution.

MOBSTER analyzes mutations separately by karyotype (for example 1:1, 2:1 and 2:2), and the same split makes a good quality check: in each karyotype, the clonal cluster should sit where theory predicts.

### B. Upstream errors in purity, ploidy and copy-number calls

Every CCF depends on the purity and copy-number estimates supplied as input, so an upstream error affects every mutation.

Purity and ploidy are often ambiguous. A copy-number caller can fit the same data almost equally well with different combinations, such as a near-diploid tumor at one purity and a near-tetraploid tumor at another. Picking the wrong solution rescales every CCF and can create or remove whole clusters. Smaller purity errors move the clonal peak. If purity is underestimated, clonal mutations appear to have CCFs above 1; if it is overestimated, they appear subclonal and may be reported as a large subclone.

Copy-number calls from WES are noisier than calls from WGS. Off-target reads, uneven capture and fewer heterozygous SNPs make segmentation and allele-specific calls less reliable from exomes. In the Tanner et al. WES benchmark, the best pipeline combined Mutect2, FACETS and PyClone-VI, so the choice of copy-number caller mattered. Other sources of error include tumor cells contaminating the matched normal (which removes true mutations or lowers their VAFs), poorly mappable or repetitive regions, sex chromosomes, and segmentation errors near breakpoints.

In the ICGC-TCGA DREAM benchmark (Salcedo et al., *Nat Biotechnol* 2024, [DOI](https://doi.org/10.1038/s41587-024-02250-y)), copy-number state, purity-adjusted depth and read mappability affected the performance of most methods. Before clustering, confirm that the clonal mutation peak lies at CCF ≈ 1 in each major karyotype. Compare purity estimates from at least two independent methods (such as ABSOLUTE, ASCAT, FACETS or Sequenza) with the purity implied by the mutation VAFs, and look at VAF histograms for each copy-number state. A quality-control tool such as CNAqc (Antonello et al., *Genome Biol* 2024) tests whether the copy-number calls and mutation VAFs agree.

### C. The neutral tail

Cells in every expanding clone keep acquiring new mutations, each carried by a small fraction of cells. Under neutral growth these form a long tail of low-frequency mutations whose number keeps rising toward lower VAF. The tail follows a power law, with the cumulative number of mutations above frequency $f$ roughly proportional to $1/f$. It is real biology, but it does not correspond to a distinct subpopulation.

Standard mixture models assume every mutation belongs to a discrete cluster, so they tend to fit the tail with one or more spurious low-CCF subclones. Deep WGS makes this worse because it detects many tail mutations, and the caller's detection limit cuts the tail off so that it looks like a hump just above the threshold. MOBSTER (Caravagna et al., *Nat Genet* 2020, [DOI](https://doi.org/10.1038/s41588-020-0675-5)) models the tail as a Pareto component alongside beta-distributed clusters. Its authors showed that ignoring the tail causes systematic errors, especially across multiple samples. Treat the lowest-CCF cluster in any analysis as suspect until other evidence supports it, particularly if it sits just above the detection limit.

### D. Multi-region and longitudinal sampling

A single bulk sample can only separate subclones whose CCFs differ in that sample. Two distinct subclones with similar prevalence in one region merge into one cluster, and one sample usually cannot tell a linear tree from a branching one. With several regions or time points, each mutation has a vector of CCFs, and clusters that overlap in one sample are often well separated in another. Multi-sample data also reveal region-specific subclones and constrain the tree topology much more tightly.

Sequencing several samples at standard depth can therefore be more informative than sequencing one sample very deeply. The MAGOS authors (Ahmadinejad et al., *Mol Biol Evol* 2022, [DOI](https://doi.org/10.1093/molbev/msac136)) made this case on cost grounds.

Multi-region designs have their own requirements. Force-call every mutation in every sample (Section 3C), so that a missing mutation reflects biology rather than detection limits. Estimate purity and copy number separately for each sample, since both often differ between regions. Use methods built for multiple samples, such as PyClone-VI, DPClust, Pairtree, CONIPHER, PhyClone or DeCiFer. Watch for sampling bias too: a mutation that is clonal in every sampled region may still be subclonal in the tumor as a whole (the "illusion of clonality"), and low-purity regions add little information.

## 5. Intuitive guide for bioinformaticians

*This section explains the ideas above without the statistics, and what they mean for your pipeline.*

### 1. Coverage is focus: the blurry eye chart

Think of a histogram of your VAFs as an eye chart. The trunk clone is a large letter at the top, the subclone a smaller letter further down, and sequencing depth is your glasses. At 40x WGS the glasses are poor. Every mutation's VAF is blurred across the chart, and if the trunk and subclone are close together their blur merges into one blob. The algorithm cannot sharpen the image: if you cannot see two humps in the histogram, the clustering tool will usually report one clone. At 200x or 300x the glasses are strong, mutations sit close to their true VAFs, the overlap shrinks, and the algorithm can separate the two groups. The depth that counts is the depth at each site after allowing for purity. A tumor at 30% purity sequenced at 200x carries as much tumor information as a pure tumor at 60x.

### 2. Why a subclone needs enough mutations

Clustering tools follow a form of Occam's razor and prefer the simplest explanation. Given a large peak of 100 mutations at 25% VAF and three stray mutations at 10%, the algorithm effectively asks: *"Is it more likely that a separate population of cancer cells exists just for these three mutations, or that these are noisy reads of the main tumor?"* It almost always picks the second answer. Depending on how far the subclone sits from the trunk, it typically needs from about 10 to more than 50 detected mutations before the algorithm accepts it. These settings can be tuned, so check that your result does not depend on a particular prior or random seed.

### 3. The "Death Valley" of variant callers

Variant callers are built to avoid false positives, and most pipelines discard mutations with only a few supporting reads, so the caller can remove much of a subclone before clustering starts. If your subclone sits at 6% VAF and average depth is 100x, a typical subclonal mutation has 6 variant reads, and about 1 in 7 will fall below a 4-read cutoff by chance. At 4% VAF almost half are lost, and at 2% nearly all. The survivors are the ones that happened to get extra reads, so the subclone also looks closer to the trunk than it really is. If the clustering tool needs 20 detected mutations to accept a subclone at 4% VAF, the subclone must carry about 35-40 real mutations for enough of them to survive.

### 4. Copy number moves the humps

A mutation's VAF depends on how many cells carry it, how many copies of the chromosome those cells have, and how many of those copies carry the mutation. In a genome-doubled tumor, clonal mutations form two humps, one at half the VAF of the other, and a diploid model calls the lower hump a subclone. Before trusting a subclone, check which copy-number states its mutations come from. A "subclone" made up mostly of mutations from gained, lost or doubled regions is very likely a copy-number artifact.

### 5. Purity and copy-number errors carry through

Every CCF is calculated from the purity and copy-number estimates you supply. If purity is overestimated, the whole distribution shifts left and the clonal peak is reported as a large subclone. Check that the clonal peak lands at CCF ≈ 1 before you interpret anything below it.

### 6. A low bump is not always a subclone

Growing tumors constantly pick up private mutations carried by few cells. Deep sequencing detects them as a long low-VAF tail, which clustering tools tend to split off as a fake minor subclone, especially just above the caller's detection limit. Treat the lowest cluster with suspicion, or use a tool such as MOBSTER that models the tail.

### 7. More samples can beat more depth

Two subclones that look identical in one biopsy often look different in another. If the budget allows, sequencing two or three regions at standard depth often reveals more structure than sequencing one region twice as deep, and it is the only reliable way to resolve branching evolution. Force-call all mutations in all samples so that a mutation that went undetected is not counted as absent.

### Summary

Subclonal deconvolution can only find what the data support. With tumor purity below about 30% and standard 100x WES, only well-separated subclones that carry many mutations are detectable. Purer samples, deeper sequencing, more regions, and checked purity and copy-number calls improve results far more than a more sophisticated algorithm. Whatever tool you use, interpret its clusters in light of copy number, purity and the neutral tail before treating them as cell populations.
