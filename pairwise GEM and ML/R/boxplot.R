library(magrittr) 
library(dplyr)    
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggprism)

result_path = '..\\model_performance_v2.csv'
score <- read.table(result_path,sep=',', header=TRUE, fileEncoding="UTF-8-BOM")


summary <- score %>% group_by(classifier) %>% summarise_each(funs(mean, sd))
summary <- data.frame(summary)

write.csv(summary,"..\\model_performance_summary.csv")

pwc <- score %>%
  pairwise_t_test(
    F1 ~ classifier, paired = FALSE,
    p.adjust.method = "BH"
  )
pwc

pwc <- pwc %>% add_xy_position(x = "classifier")
ggboxplot(score, x = "classifier", y = "F1", color = "classifier", 
          add = "jitter", legend = "none",  short.panel.labs = FALSE, width = 0.4) +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(score$F1), linetype = 2)+
  stat_pvalue_manual(pwc, hide.ns = TRUE)

ggsave("images\\ML_performance_f1score_v2.pdf")


result_path = '..\\count.csv'
score <- read.csv(result_path,sep=',', header=TRUE, fileEncoding="UTF-8-BOM")
rownames(score) <- score$patient.id
score <- score %>% select(-(patient.id))

singleGEM <- ggboxplot(score, x = "label", y = "Single.GEM", color = "label", 
          add = "jitter", legend = "none",  short.panel.labs = FALSE,
          palette = c("#F79B9B", "#9BB9F4"), width = 0.2) +
          rotate_x_text(angle = 0)

ggsave("images/singleGEM.pdf", plot = singleGEM)


pairGEM <- ggboxplot(score, x = "label", y = "Pairwise.GEM", color = "label", 
          add = "jitter", legend = "none",  short.panel.labs = FALSE,
          palette = c("#F79B9B", "#9BB9F4"), width = 0.2) +
          rotate_x_text(angle = 0)
ggsave("images/PairGEM.pdf", plot = pairGEM)

sp <- ggscatter(score, x = "Single.GEM", y = "Pairwise.GEM",
                color = "label",
                palette = c("#F79B9B", "#9BB9F4"),
                size = 3)
sp
ggsave("images/scatter.pdf", plot = sp)

