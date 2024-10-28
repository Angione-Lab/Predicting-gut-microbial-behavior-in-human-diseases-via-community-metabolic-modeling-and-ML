library(tidyverse)
library(umap)
library(rstatix)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw(18))

#UMAP
growth_rate = read.csv('../community_extracted_exRxns/community_growth_rates.csv')
growth_rate <- growth_rate[!duplicated(growth_rate$Model), ]

rownames(growth_rate) <- growth_rate$Model
growth_rate <- growth_rate %>% select(-c("Model"))

growth_rate_meta <- growth_rate %>%
  select(c(Patient, TypeOfInteraction))


growth_rate_meta$ID <- rownames(growth_rate)


growth_rate <- growth_rate %>% 
  select(-c(Patient, TypeOfInteraction))%>%
  mutate(ID=row_number()) 

custom.config <- umap.defaults
custom.config$n_neighbors <- 10
custom.config$min_dist <- 0.5

set.seed(142)
umap_fit <- growth_rate %>%
  select(where(is.numeric)) %>%
  scale() %>% 
  umap(config=custom.config)



umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  # rename(UMAP1="V1",
  #        UMAP2="V2") %>%
  mutate(ID=rownames(umap_fit$layout)) %>%
  inner_join(growth_rate_meta, by="ID")


umap_df %>%
  ggplot(aes(x = V1, 
             y = V2, 
             color = Patient,
             # shape = TypeOfInteraction
             ))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.text.y = element_text(size = 18), 
        plot.title = element_text(size=20,  hjust = 0.5))+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "Growth rate UMAP") +
  # scale_color_manual(values = c(" Neutralism" = "red",
  #                               " Amensalism" = "blue",
  #                               " Competition" = "green"),
  #                    )

  scale_color_manual(values = c("COVID" = "#F79B9B",
                                "Healthy" = "#9BB9F4"))
  
     

ggsave("images/Community_growth_rate_UMAP.pdf")


#%%
#boxplot
comm_gr_rate = read.csv('../community_extracted_exRxns/community_growth_rates.csv')
comm_gr_rate <- comm_gr_rate[!duplicated(comm_gr_rate$Model), ]

rownames(comm_gr_rate) <- comm_gr_rate$Model
comm_gr_rate <- comm_gr_rate %>% select(-c("Model"))

interactions <- comm_gr_rate %>% select("GRSpeciesAFull", "GRSpeciesBFull", "GRASolo" , "GRBSolo", "TypeOfInteraction")
sep_into_cols <- colnames(interactions %>% select(-c('TypeOfInteraction')))


pf <- interactions %>%
  pivot_longer(
    sep_into_cols,
    names_to = "model_type",
    values_to = "Growth_rate",
    values_drop_na = TRUE
  )


pf$Growth_rate <- as.numeric(pf$Growth_rate)

pwc <- pf %>%
  group_by(`model_type`) %>%
  pairwise_t_test(
    Growth_rate ~ TypeOfInteraction, paired = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc

pwc <- pwc %>% add_xy_position(x = "model_type")

p <- ggboxplot(pf, x = "model_type", y = "Growth_rate", color = "TypeOfInteraction", alpha= 1,
               #palette = c("#F79B9B", "#9BB9F4"),
               add = "mean_se",   short.panel.labs = TRUE,
               ggtheme = theme_bw(base_size = 20),
               title = "Community interaction type and growth rate" ) +
  geom_hline(yintercept = mean(pf$Growth_rate), linetype = 2)+
  stat_pvalue_manual(pwc, hide.ns = TRUE)

p <- p + theme(axis.text.x = element_text(angle = 45, hjust=1), 
               plot.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "top",
               axis.text.y = element_text(size = 12), 
               plot.title = element_text(size=18,  hjust = 0.5)) #+ stat_pvalue_manual(pwc, hide.ns = TRUE)

#p <-p + coord_flip()
print(p)


ggsave(
  "images/Community_interaction_type_box_plot.pdf",
  plot = p,
  device = NULL,
  #path = 'Results/metabolic_analysis/',
  scale = 1,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  bg = NULL
)


