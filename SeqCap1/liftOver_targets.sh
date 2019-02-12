#!/bin/bash

# liftOver Roche targets from hg19 to hg38
LIFTOVER=/home/shraddhapai/software/kent_utilities/liftOver
HG19toHG38=/home/shraddhapai/software/kent_utilities/hg19ToHg38.over.chain.gz

primaryFile=/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/input/OID43737_hg19_161128_design_deliverables/OID43737_hg19_161128_primary_targets.bed
outDir=`dirname $primaryFile`
baseF=`basename $primaryFile .bed`
outF=${outDir}/${baseF}
$LIFTOVER $primaryFile $HG19toHG38 ${outF}.hg38 ${outF}.unmapped

captureFile=/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/input/OID43737_hg19_161128_design_deliverables/OID43737_hg19_161128_capture_targets.bed
outDir=`dirname $captureFile`
baseF=`basename $captureFile .bed`
outF=${outDir}/${baseF}
$LIFTOVER $captureFile $HG19toHG38 ${outF}.hg38 ${outF}.unmapped
