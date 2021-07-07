# simulate genome:
mini1.blocks.df <- data.frame(name = c("A", "B", "S", "A", "C"),
                              length = c(500, 500, 1000000, 500, 500))
mini1 <- simulate.genome(genome.name = "mini1", out.dir = "genomes/mini1/", blocks.df = mini1.blocks.df, 
                         base.seed = 42)

mini1.reads.se <- simulate.reads(fa = mini1$seq.file, start = 1001001, read.length = 75, bPairEnd = F, 
                                # insert.size = 500, 
                                n.reads = 100, interval = 1, 
                                out.dir = "mini1_se/fastq/", root.name = "mini1_se")

# mostly default parameters of bowtie2
bowtie2.wrapper(fastq.or.pair = mini1.reads.se %>% unlist(), genome.index = mini1$index, 
                out.bam = "mini1_se/bam/mini1_se_aln_1.bam", threads = 1, a = F, k = NULL,
                X = 500, report.unaligned = T, preset = "", defaults = "")

# default parameters of bwa
bwa.mem(fastq.vec = mini1.reads.se %>% unlist(), ref.fa = mini1$seq.file, 
        out.bam = "mini1_se/bam/mini1_se_aln_2.bam", threads = 1)

# simulate pair-end reads
mini1.reads.pe <- simulate.reads(fa = mini1$seq.file, start = 1, read.length = 75, bPairEnd = T, 
                                 insert.size = 500, 
                                 n.reads = 100, interval = 1, 
                                 out.dir = "mini1_pe/fastq/", root.name = "mini1_pe")
# bowtie2 alignment
bowtie2.wrapper(fastq.or.pair = mini1.reads.pe %>% unlist(), genome.index = mini1$index, 
                out.bam = "mini1_pe/bam/mini1_pe_aln_1.bam", threads = 1, a = T, k = NULL,
                X = 500, report.unaligned = T)

Rsamtools::bamFlagAsBitMatrix(as.integer(c(353, 97, 145))) %>% 
.[, c("isFirstMateRead", "isSecondaryAlignment")] %>% 
  `rownames<-`(as.character(c(353, 97, 145)))

# bowtie2 alignment, relax -X
bowtie2.wrapper(fastq.or.pair = mini1.reads.pe %>% unlist(), genome.index = mini1$index, 
                out.bam = "mini1_pe/bam/mini1_pe_aln_2.bam", threads = 1, a = T, k = NULL,
                X = 700, report.unaligned = T)

Rsamtools::bamFlagAsBitMatrix(as.integer(c(99, 147))) %>% 
  .[, c("isFirstMateRead", "isSecondaryAlignment")] %>% 
  `rownames<-`(as.character(c(99, 147)))

