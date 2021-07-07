options(stringsAsFactors = F)
options(scipen = 20)
options(bitmapType='cairo')
suppressMessages(library(dplyr))
suppressMessages(library(parallel))
suppressMessages(library(ggplot2))
suppressMessages(library(GenomicRanges))
source("functions.R")
source("helpers.R")

BOWTIE2.DEFAULTS <- " --xeq --dovetail"
BOWTIE2 <- "/opt/apps/bowtie2/2.3.4.1/bowtie2"
SAMTOOLS <- "/bar/cfan/anaconda2/envs/jupyter/bin/samtools"
BWA.TARGET <- "~/software/bwa/bwa_0.7.16a-r1181/bwa-0.7.16/bwa"
