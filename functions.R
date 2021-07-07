simulate.genome <- function(genome.name, out.dir, blocks.df, 
                            samtools = SAMTOOLS,
                            base.seed = 42, 
                            bwa = "/opt/apps/bwa/0.7.15/bwa") {
  # the format of blocks.df: name and length. name will also be used to generate seed
  blocks.df$seed <- base.seed + factor(blocks.df$name) %>% as.numeric()
  if (is.null(blocks.df$mutation.rate)) {
    blocks.df$mutation.rate <- 0
  }
  # generate sequence:
  seq <- lapply(seq_along(blocks.df$seed), function(i) {
    set.seed(seed = blocks.df$seed[i])
    seq.sub <- sample(c("A", "T", "C", "G"), size = blocks.df$length[i], replace = T)
    if (blocks.df$mutation.rate[i] > 0) {
      pos.mutate <- sample(1:length(seq.sub), size = blocks.df$mutation.rate[i] * length(seq.sub), 
                           replace = F)
      seq.sub[pos.mutate] <- "G"
    }
    return(seq.sub)
  }) %>% unlist()
  dir.create(out.dir, showWarnings = F, recursive = T)
  out.fa <- paste0(out.dir, "/", genome.name, ".fa")
  seqinr::write.fasta(sequences = seq, names = genome.name, file.out = out.fa)
  cmd <- paste0(samtools, " faidx ", out.fa)
  print(cmd); system(cmd)
  
  cmd <- paste0(bwa, " index ", out.fa)
  print(cmd); system(cmd)
  index <- bowtie2.index(fa = out.fa, threads = 4, 
                                         log.file = paste0(out.dir, "/bowtie2-build.log"))
  
  blocks.df$end <- cumsum(blocks.df$length)
  blocks.df$start <- c(1, 1+blocks.df$end[-length(blocks.df$end)])
  blocks.df$seqnames <- genome.name
  gr <- blocks.df %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  bed.file <-  paste0(out.dir, "/", genome.name, "_construct.bed")
  trash <- utilsFanc::write.zip.fanc(df = blocks.df[, c("seqnames", "start", "end", "name")], out.file = bed.file,
                                     bed.shift = T)
  res <- list(genome.name = genome.name, seq = seq, seq.file = out.fa, index = index,
              gr = gr, bed.file = bed.file, out.dir = out.dir, blocks.df = blocks.df)
  return(res)
}

simulate.reads <- function(fa, start, read.length, bPairEnd = T,
                           insert.size, n.reads, interval = 1,
                           out.dir, root.name, zip = T) {
  if (file.exists(fa))
    fa <- seqinr::read.fasta(fa, forceDNAtolower = F)[[1]] %>% as.character()
  
  ref.length <- length(fa)
  
  starts <- seq(from = start, by = interval, length.out = n.reads)
  df <- data.frame(qname = 1:length(starts) %>% paste0("seq", .), 
                   start.1 = starts, end.1 = starts + read.length - 1)
  df$read.1 <- df %>% split(f = 1:nrow(df)) %>% 
    sapply(function(x) return(fa[x$start.1:x$end.1] %>% paste0(collapse = ""))) %>% 
    `names<-`(NULL)
  df <- df %>% filter(start.1 <= ref.length, end.1 <= ref.length)
  if (bPairEnd == T) {
    df$start.2 <- df$start.1 + read.length + insert.size - 1
    df$end.2 <- df$start.2 + read.length -1 
    df$read.2 <- df %>% split(f = 1:nrow(df)) %>% 
      sapply(function(x) return(fa[x$start.2:x$end.2] %>% rev() %>% 
                                  paste0(collapse = "") %>% chartr("ATGC","TACG",.))) %>% 
      `names<-`(NULL) 
    df <- df %>% filter(start.2 <= ref.length, end.2 <= ref.length)
  }
  
  ### now write things out
  res <- list()
  
  if (bPairEnd == T) {
    res$read.1 <- write.fastq.fanc(df = df, seq.col = "read.1", name.col = "qname", zip = zip, 
                                   file.out = paste0(out.dir, "/", root.name, "_1.fastq"))
    res$read.2 <- write.fastq.fanc(df = df, seq.col = "read.2", name.col = "qname", zip = zip, 
                                   file.out = paste0(out.dir, "/", root.name, "_2.fastq"))
  } else {
    res$read.1 <- write.fastq.fanc(df = df, seq.col = "read.1", name.col = "qname", zip = zip, 
                                   file.out = paste0(out.dir, "/", root.name, ".fastq"))
  }
  return(res)
}


write.fastq.fanc <- function(df, seq.col, name.col, quality.col = NULL, zip, file.out) {
  df[, name.col] <- df[, name.col] %>% paste0("@", .)
  df$third <- "+"
  # browser()
  if (is.null(quality.col)) {
    df$quality <- sapply(df[, seq.col] %>% nchar(), 
                         function(x) return(paste0(rep("C", x), collapse = "")))
    quality.col <- "quality"
  }
  
  df <- df[, c(name.col, seq.col, "third", quality.col)]
  dir.create(dirname(file.out), showWarnings = F, recursive = T)
  if (grepl("gz$", file.out))
    file.out <- file.out %>% sub(".gz$", "", .)
  write.table(x = df, file = file.out, quote = F, sep = "\n" ,row.names = F, col.names = F)
  if (zip == T) {
    cmd <- paste0("gzip -f ", file.out)
    print(cmd); system(cmd)
    file.out <- paste0(file.out, ".gz")
  }
  
  return(file.out)
}

bowtie2.wrapper <- function (fastq.or.pair, genome.index, align.out.dir, out.bam = NULL, 
          threads = 6, a = F, k = 10, X = 10000, mm = F, report.unaligned = F, 
          preset = " --very-sensitive ", defaults = BOWTIE2.DEFAULTS, 
          add.params = "", run = T, log.file = NULL, bowtie2 = BOWTIE2, 
          samtools = SAMTOOLS) 
{
  if (is.null(out.bam)) {
    root.name <- trim.fastq(fastq.or.pair[1], T)
    align.out.dir <- normalizePath(align.out.dir, mustWork = F)
    out.bam <- paste0(align.out.dir, "/", root.name, ".bam")
  }
  if (length(fastq.or.pair) == 1) 
    cmd.fastq <- paste0(" -U ", fastq.or.pair)
  else if (length(fastq.or.pair) == 2) 
    cmd.fastq <- paste0(" -1 ", fastq.or.pair[1], " -2 ", 
                        fastq.or.pair[2])
  else stop("must provide one fastq or a pair of fastqs")
  if (a == T) {
    report <- " -a "
  }
  else if (!is.null(k)) {
    report <- paste0(" -k ", k)
  }
  else {
    report <- ""
  }
  if (mm == T) {
    mm <- " --mm "
  }
  else {
    mm <- ""
  }
  if (report.unaligned == F) {
    unal <- " --no-unal "
  }
  else {
    unal <- ""
  }
  dir.create(dirname(out.bam), showWarnings = F, recursive = T)
  cmd <- paste0(bowtie2, " -x ", genome.index, " ", cmd.fastq, 
                " ", preset, " ", defaults, " ", report, " -p ", threads, 
                " ", mm, " ", unal, " -X ", X, " ", add.params)
  if (!is.null(log.file)) {
    dir.create(dirname(log.file), showWarnings = F, recursive = T)
    cmd <- paste0(cmd, " 2>", log.file, " ")
  }
  cmd <- paste0(cmd, " | ", samtools, " sort - -O bam -m 2G -o ", 
                out.bam, " -@ ", threads)
  cmd.exec.fanc(cmd = cmd, stdout.file = NULL, intern = F, 
                           run = run)
  system(paste0(samtools, " index ", out.bam))
  if (!file.exists(out.bam)) 
    stop(paste0(out.bam, " failed to generate"))
  else return(out.bam)
}

bwa.mem <- function (fastq.vec, ref.fa, out.bam, threads, samtools = SAMTOOLS, 
          bwa = BWA.TARGET, mapped.only = F, add.params = "", 
          log.file = NULL) {
  dir.create(dirname(out.bam), showWarnings = F, recursive = T)
  cmd <- paste0(bwa, " mem -t ", threads, " ", ref.fa, " ", 
                paste0(fastq.vec, collapse = " "), " ", add.params)
  if (!is.null(log.file)) 
    cmd <- paste0(cmd, " 2>", log.file)
  if (mapped.only == T) 
    cmd <- paste0(cmd, " | ", samtools, " view -h -@ ", threads, 
                  " -F 0x4 - ")
  cmd <- paste0(cmd, " | ", samtools, " sort - -O bam -o ", 
                out.bam, " -@ ", threads)
  utilsFanc::cmd.exec.fanc(cmd, stdout.file = NULL, run = T, 
                           intern = F)
  if (!file.exists(out.bam)) 
    stop(paste0(out.bam, " failed to generate"))
  cmd <- paste0(samtools, " index ", out.bam)
  cmd.exec.fanc(cmd, run = T, intern = F)
  return(out.bam)
}