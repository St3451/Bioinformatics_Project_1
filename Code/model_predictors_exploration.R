#### Model building predictor variables exploration and plotting model coefficients 

library(tidyverse)
library(corrplot)
library(RColorBrewer)
library(gridExtra)

## Data exploration
sim_matrix <- read_delim("output_cos_sim_samples_signatures_exomes.txt", ",")
metadata <- read_delim("output_metadata_case_samples.txt", ",")

metadata$Response[is.na(metadata$Response)] = "1"

save_plot <- function(filename, plot, width = 20, height = 12, dpi = 300){
  ggsave(filename,
         plot,
         width = width,
         height = height,
         dpi = dpi) 
}

metadata %>% 
  gather(key = model_features, value = value,
         Age, TMB) %>%
  ggplot(aes(x=Response, y = value, fill = Response))+
  geom_boxplot()+
  facet_grid(~model_features, scales="free") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  theme_bw(base_size = 15) +
  ylab("Value Age TMB") + 
  theme(legend.position='none') -> plot1

  
metadata %>%
  ggplot(aes(x = Gender, fill = Response)) +
  geom_bar(color = "black", position = "dodge") +
  labs(fill = "Response",
       title = "asd") +
  xlab("Gender") + 
  ylab("Frequency") + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  theme_bw(base_size = 15) -> plot2


metadata %>% 
  mutate(Gender = TMB) %>%
  gather(key = model_features, value = value,
         Age, Gender) %>%
  ggplot(aes(x=Response, y = value, fill = Response))+
  geom_boxplot()+
  facet_grid(~model_features, scales="free") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  theme_bw(base_size = 15) +
  ylab("Value Age TMB") + 
  theme(legend.position='none') -> plot1



plot12 <- grid.arrange(plot1, plot2, ncol = 2, widths = 5:2)


save_plot("plot_predictor_metadata3.png", plot12,
          width = 8, height = 4)


## Correlogram case samples
df <- bind_cols(sim_matrix, metadata[,c("Age", "TMB")])
colnames(df) <- gsub("S_", "s. ", c(colnames(sim_matrix), "Age", "TMB"))

corr <- cor(df[,2:ncol(df)])
corrplot(corr,
         order = "hclust",
         type = "upper",
         col=brewer.pal(n=8, name="RdYlBu"))


## Correlogram complete cohort
sim_matrix <- read_delim("output_cos_sim_samples_signatures_exomes_all_samples.txt", ",")
metadata <- read_delim("output_metadata_all_samples.txt", ",")
TMB <- read_delim("TMB_all_samples.txt", ",")

df <- bind_cols(sim_matrix, metadata$Age, TMB$TMB)
colnames(df) <- gsub("S_", "s. ", c(colnames(sim_matrix), "Age", "TMB"))

df %>% filter(!is.na(Age)) -> df

corr <- cor(df[,2:ncol(df)])
corrplot(corr,
         order = "hclust",
         type = "upper",
         col=brewer.pal(n=8, name="RdYlBu"))

## Model coefficients
filename = "rforest_features_importance_100"
logistic_weights <- read.table(paste(filename, ".txt", sep = ""), sep = "\t", header = TRUE)

logistic_weights %>% 
  as_tibble() %>% 
  rename(Predictors = X) %>%
  arrange(desc(abs(Log_weights))) %>%
  mutate(Predictors = factor(Predictors,
                             levels = Predictors)) %>%
  mutate(Type = ifelse(substr(Predictors, 1, 2) %in% "S_", "Signature", 
                       ifelse(Predictors == "Intercept", "Intercept",
                              "Metadata"))) -> logistic_weights
logistic_weights

save_plot <- function(filename, plot, width = 20, height = 12, dpi = 300){
  ggsave(paste("results_modeling/", 
               filename, 
               sep = ""), 
         plot,
         width = width,
         height = height,
         dpi = dpi) 
}

logistic_weights %>%
  ggplot(aes(x = Predictors, y = Log_weights, fill = Type)) +
  geom_bar(stat = "identity", col = "black") + theme_bw() +
  labs(title = "Random forest features important",
       fill = "Variable type") +
  xlab("Predictor variables") + 
  ylab("Coefficient") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust=1,
                                   size = 7),
        plot.title = element_text(size = 25,
                                  hjust = 0.5)) -> logistic_weights_barplot

save_plot(paste(filename, "_barplot_new01.png", sep = ""),
          logistic_weights_barplot,
          width = 10, height = 6)
logistic_weights_barplot
