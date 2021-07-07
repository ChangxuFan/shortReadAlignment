cmd.exec.fanc <- function (cmd, stdout.file = NULL, intern, run = T) {
  if (!is.null(stdout.file)) {
    system("mkdir -p " %>% paste0(dirname(stdout.file)))
    cmd <- paste0(cmd, " 1>", stdout.file, " 2>&1")
  }
  cat(cmd)
  cat("\n")
  if (run == T) {
    res <- system(cmd, intern = intern)
    if (intern == T) 
      return(res)
    else {
      if (res != 0) 
        stop("the following command returned non-zero exit code: " %>% 
               paste0(cmd))
      return(res)
    }
  }
  return(cmd)
}

bowtie2.index <- function (fa, root.name = NULL, out.dir = NULL, threads, bowtie2.build = "/opt/apps/bowtie2/2.3.4.1/bowtie2-build", 
                           log.file = NULL) {
  if (is.null(root.name)) 
    root.name <- basename(fa) %>% sub(".fn*a.*$", "", .)
  if (is.null(out.dir)) 
    out.dir <- paste0(dirname(fa), "/bowtie2/")
  dir.create(out.dir, showWarnings = F, recursive = T)
  cmd <- paste0(bowtie2.build, " --threads ", threads, " ", 
                fa, " ", out.dir, "/", root.name)
  cmd.exec.fanc(cmd = cmd, stdout.file = log.file, 
                intern = F, run = T)
  return(paste0(out.dir, "/", root.name))
}


insert.name <- function (name, insert, ext, trim.dir = F) {
  if (trim.dir == T) 
    name <- basename(name)
  inserted <- sub(paste0("(.*)(", ext, ")"), paste0("\\1_", 
                                                    insert, "\\2"), name)
  return(inserted)
}

bam.filter <- function (in.bam, out.bam = NULL, creat.index = F,  
          filter.expression, remove.mate = F, debug = F, thread, 
          samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools", 
          sambamba = "~/software/sambamba/sambamba") {
  if (is.null(out.bam)) 
    out.bam <- insert.name(in.bam, insert = "filtered", ext = ".bam", 
                           trim.dir = F)
  out.bam.tmpt <- insert.name(out.bam, insert = "tempt", ext = ".bam", 
                              trim.dir = F)
  filter.cmd <- paste0(sambamba, " view -f bam -F ", "\"", 
                       filter.expression, "\"", " -t ", thread, " ", in.bam, 
                       " > ", out.bam.tmpt)
  print(filter.cmd)
  system(filter.cmd)
  if (remove.mate == T) {
    bam.nsort <- insert.name(name = out.bam.tmpt, insert = "nsort", 
                             ext = ".bam", trim.dir = F)
    cmd <- paste0(samtools, " sort -n ", " -@ ", thread, 
                  " -o ", bam.nsort, " ", out.bam.tmpt)
    print(cmd)
    system(cmd)
    bam.nsort.fixmate <- insert.name(name = bam.nsort, insert = "fixmate", 
                                     ext = ".bam", trim.dir = F)
    cmd <- paste0(samtools, " fixmate ", bam.nsort, " ", 
                  bam.nsort.fixmate)
    print(cmd)
    system(cmd)
    cmd <- paste0(samtools, " sort ", " -@ ", thread, " ", 
                  bam.nsort.fixmate, " | ", sambamba, " view -f bam -F ", 
                  "\"", "paired", "\"", " -t ", thread, " ", "/dev/stdin", 
                  " > ", out.bam.tmpt)
    print(cmd)
    system(cmd)
    if (debug == F) {
      system(paste0("rm ", bam.nsort, " ", bam.nsort.fixmate))
    }
  }
  cmd <- paste0("mv ", out.bam.tmpt, " ", out.bam)
  print(cmd)
  system(cmd)
  if (creat.index == T) {
    cmd <- paste0(samtools, " index ", out.bam)
    print(cmd)
    system(cmd)
  }

  if (file.exists(out.bam)) 
    return(out.bam)
  else (stop(paste0("bam file was not successfully filtered for: ", 
                    in.bam)))
  return(out.bam)
}