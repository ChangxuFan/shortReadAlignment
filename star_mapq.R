# first, filter alignments for primary alignments only:

bam.filter(in.bam = "star_mapq/Lin_neg_star_nkc.bam",
           filter.expression = "not secondary_alignment", remove.mate = F, 
           creat.index = T, thread = 1)

star.reads <- Rsamtools::scanBam(file = "star_mapq/Lin_neg_star_nkc_filtered.bam")[[1]] %>%
  as.data.frame()

star.reads %>% ggplot(aes(x = mapq)) +
  geom_density(adjust = 0.5) +
  xlim(c(0, 270)) +
  ggsave("star_mapq/distro.png")

star.reads$mapq %>% table()
# 0     1     3   255 
# 89   238  1990 45625 