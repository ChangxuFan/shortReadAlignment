mini1.reads.pe <- simulate.reads(fa = mini1$seq.file, start = 1, read.length = 75, bPairEnd = F, 
                                 insert.size = 500, 
                                 n.reads = 100, interval = 1, 
                                 out.dir = "mini1_pe/fastq/", root.name = "mini1_se")