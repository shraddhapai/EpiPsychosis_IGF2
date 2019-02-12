# EpiPsychosis_IGF2
This repo contains code used to generate results for the manuscript:

**Differential DNA methylation of an enhancer at the IGF2 locus affects dopamine synthesis in patients with major psychosis***.
Pai S., Li P., Killinger B., Marshall L., Jia P., Liao J., Petronis A., Szabó P.E., and Labrie V.
[ Publication info ] 

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

## Methylation at IGF2 locus
* Figure 2a: `targetView.R`
* Figure 2b: `dmp_plotGroupMeans.R`
* Supplementary Figure 8: Case/control IGF2 methylation plotted with lifestyle variables: `relateIGF2_confounds.R`

# Case-Control genotyping

# Case-Control targeted bisulfite sequencing

# Case-Control transcriptomes





  
