# EpiPsychosis_IGF2
This repo contains code used to generate results for the manuscript:

**Differential DNA methylation of an enhancer at the IGF2 locus affects dopamine synthesis in patients with major psychosis***.
Pai S., Li P., Killinger B., Marshall L., Jia P., Liao J., Petronis A., Szabó P.E., and Labrie V.

bioRXiv preprint: https://doi.org/10.1101/296756  

# Prerequisites
Most of the analysis uses R and BioConductor packages. In addition, the standalone
packages dataExplore (http://github.com/shraddhapai/dataExplore) and IlluminaEPICtools
(contained in `IlluminaEPICTools` in this repo) are required to be installed.

# Case-Control EPIC microarrays
Location: `CaseControlEpic/` folder

## General
* Array preprocessing QC, normalization, MDS: `processInput.R`
* Post-processing QC and data exploration, data cleaning, merging with processed genetic data: `RunMe.R` (Supplementary Figure 3,4)
* Call differentially-methylated probes and plotting results: 
  * Main result reported in paper (all samples): `callDMPs.R` (Supplementary Figure 6a)
  * CEU samples: `callDMPs_CEU.R` 
  * IGF2 confirmation with CEU males: `dmp_getMeanBetaDiff_CEUmale.R`
* Volcano plot in Figure 1A: `volcanoPlot.R`
* Pathway analysis and Figure 1B: Call to `dmp_pathwayORA` within `callDMPs.R`; Enrichment Map created in [Cytoscape](https://cytoscape.org/).
* Compare genotypes with genetic data: Supplementary Fiure 5a: `getEPICsnpgeno.R` and `match_genoEPIC.R`

## Methylation at IGF2 locus
* Figure 2a: `targetView.R`
* Figure 2b: `dmp_plotGroupMeans.R`
* Supplementary Figure 8: Case/control IGF2 methylation plotted with lifestyle variables: `relateIGF2_confounds.R`

## meQTL analysis
* Pre-imputation processing:  `preimpute_check.sh`, `vcf4impute.sh`
* Extracts genotypes in DMP window: `CaseControlSNP/prepareDMPwins.R`, `CaseControlSNP/getSNPsInRanges.sh`
* meQTL analysis: Supplementary Table 8 `meQTL_useLMER.R`
* Figure 1c: `meQTL_plotCis.R`

# Case-Control genotyping
Code in `CaseControlSNP/`
* Quality control: `fixPlink_coords.R`, `plinkQC.sh`
* Compute lambda: `computeLambda.R` 

# Case-Control targeted bisulfite sequencing
Code in `SeqCap1/` and `SeqCap2/`
Processing pipeline followed manufacturer's instructions from this online document: http://netdocs.roche.com/DDM/Effective/07187009001_RNG_SeqCap-EZ_TchNote_Eval-data_v2.1.pdf
* Processing: 
  * QC of fastq files: `fastQC_batch.R`
  * Trim reads: `trim.sh`
  * Count reads: `countreads.sh`
  * Align reads: `align.sh`
  * Post-alignment processing: `postalign_runner.sh` > `postalign.sh`
* Validation of IGF2 locus: Figure 2c,d: `IGF2_view.R`
* Correlate IGF2 methylation with EPIC microarrays: Supplementary Figure 9: `corrEPIC_IGF2.R`

# Case-Control transcriptomes
Code in `CaseControlRNA/`
* Processing: `SCZ_RNAseq.bash`
* Differential expression and pathway analysis: `diffEx_simple.R`. Enrichment Map created in Cytoscape.
* Volcano plot:Supplementary Figure 7a: `volcanoPlot.R`
* Correlation of expression with BrainSpan gene sets: Supplementary Figure 7b: `diffEx_annotate.R`
* Comparing genotype with SNP arrays: `RNAseq_callVar.sh`, `makePlink.sh`, `matchGeno.R`

# Mouse transcriptomes
Code in `MouseRNA/`
* Processing: `../CaseControlRNA/SCZ_RNAseq.bash`
* Differential expression: Figure 4a, Supplementary Figure 14. `diffEx.R`
* Pathway analysis: Figure 4b, Supplementary Figure 15. `runPathway_FC.R` and `runPathway_STR.R`. Enrichment map created in Cytoscape. 

# Mouse Synaptosomes
* Pathway analysis: Supplementary Figure 17. `diffEx_testAllProt.R`
  
