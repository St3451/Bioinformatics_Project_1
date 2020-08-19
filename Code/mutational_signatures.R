#!/usr/bin/env Rscript

#### 0 Load libraries and set flags

# Use the command line flag to specify if case samples only should be analysed
library(argparse)
parser <- ArgumentParser()
parser$add_argument("-c", "--case_only", action = "store_true",
                    help = "Specify this flag to use case samples only")
args <- parser$parse_args()
use_case_samples_only <- args$case_only

if (use_case_samples_only == FALSE){
  writeLines("Mutational signatures profile pipeline using all samples")
} else{
  writeLines("Mutational signatures profile pipeline using case samples only")
}

# Load other libraries
writeLines("\nLoad libraries")
suppressPackageStartupMessages({
  library(MutationalPatterns)
  library(BSgenome)
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  library(ref_genome, character.only = TRUE)
  library(gridExtra)
  library(NMF)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(tidyverse)
  library(ggrepel)
  library(viridis)
  library(RColorBrewer)
  library(pheatmap)
  library(ggbiplot)
  library(ggforce)
})

#### 1. Set paths
# Set working directory to the folder where the vcf files are located
setwd("/ngc/projects/gm_ext/stepel/data/samples_vcf_files_benchmark/")
# Set path where the signatures sets are saved
signatures_path <- "/ngc/projects/gm_ext/stepel/data/signatures_sets/"
# Set path for the outputs of the scripts
if (use_case_samples_only == TRUE){
  results_path <- "/ngc/projects/gm_ext/stepel/results/results_case_samples/"
} else{
  results_path <- "/ngc/projects/gm_ext/stepel/results/results_all_samples/"
}
# Set path for the file containing sample metadata
metadata_path <- "/ngc/projects/gm_ext/stepel/data/metadata_updated.txt"
# Set path for TBM.txt file
TBM_path <- "/ngc/projects/gm_ext/stepel/data/TMB/TMB.csv"
# Set path for sequencing date
dates_path <- "/ngc/projects/gm_ext/stepel/data/sequencing_dates.csv"


##### 2. Load VCF data and metadata
writeLines("\nLoad VCF files and metadata")

#### 2.1 VCF data
print("Load VCF data")
# Load vcf data and create GRange object
vcf_files <- list.files(path = ".", pattern = ".vcf")
if (use_case_samples_only == TRUE){
  vcf_files <- vcf_files[grep("case", vcf_files)]
}
sample_names <- str_extract(vcf_files, ".+?(?=_)")
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
summary(vcfs)


#### 2.2 Metadata
### Load metadata
print("Load metadata")
metadata <- read_tsv(metadata_path)
metadata$Subtype[metadata$Subtype == "Endometrial_+_adeno_fra_tuba"] <- "Endometrial_cancer"
metadata$Subtype[metadata$Subtype == "Endometrial + adeno fra tuba"] <- "Endometrial_cancer"
#metadata$Subtype[metadata$Subtype == "Ovarian cancer"] <- "Ovarian_cancer"
#metadata$Subtype[metadata$Subtype == "Bile duct_cancer"] <- "Bile_duct_cancer"
metadata$Subtype[metadata$Subtype == "Ovarian"] <- "Ovarian cancer"
metadata$Subtype[metadata$Subtype == "Cervical"] <- "Cervical cancer"
metadata$Subtype[metadata$Subtype == "Vulvovaginal"] <- "Vulvovaginal cancer"
metadata$Subtype[metadata$Subtype == "Pancreas/breast/ovarian"] <- "Multiple cancer"
metadata$Subtype <- gsub("_", " ", metadata$Subtype)
metadata$Subtype <- gsub("cancer", "c.", metadata$Subtype)

metadata_sorted_by_samplesID <- metadata %>% arrange(SampleID)
metadata_sorted_by_samplesID$Gender[is.na(metadata_sorted_by_samplesID$Gender)] = "NA"
metadata_sorted_by_samplesID$Response[is.na(metadata_sorted_by_samplesID$Response)] = "NA"
metadata_sorted_by_samplesID$Response[metadata_sorted_by_samplesID$Response %in% c("SD", "PR")] = "1"
metadata_sorted_by_samplesID$Response[metadata_sorted_by_samplesID$Response == "PD"] = "0"

if (use_case_samples_only == TRUE){
  metadata_sorted_by_samplesID %>% filter(Group == "Case") -> metadata_sorted_by_samplesID  
  print(metadata_sorted_by_samplesID[,1:9])
  print(metadata_sorted_by_samplesID[,9:12])
}

tumortype_count <- count(metadata_sorted_by_samplesID$Subtype) %>% 
  add_row(x = "Total", freq = sum(count(metadata_sorted_by_samplesID$Subtype)$freq))
colnames(tumortype_count) <- c("Tumor", "Count")
print("Tumor types:")
tumortype_count
tumortype_count %>% write_tsv(paste(results_path, "tumortype_count.txt", sep = ""))


#### 2.3 Find differences between metadata and vcf files
## Count no of samples in metadata and vcf files
paste("Case samples in metadata:", nrow(metadata %>% filter(Group == "Case")))
case_vcf_files <- vcf_files[grep("case", vcf_files)]
paste("Case samples in vcf files:", length(case_vcf_files))

if (use_case_samples_only == FALSE){
  paste("Control samples in metadata:", nrow(metadata %>% filter(Group == "Control" | Group == "Benchmark")))
  control_vcf_files <- vcf_files[grep("control", vcf_files)]
  paste("Control samples in vcf files:", length(control_vcf_files))
  paste("Total samples:", length(case_vcf_files) + length(control_vcf_files))
}

## Find missing files that are present as vcf but not in metadata
missing_samples <- setdiff(sample_names, metadata_sorted_by_samplesID$SampleID)
no_missing_samples <- length(missing_samples)
# Generate a regex pattern object with the samples ID that differ between the two files
c <- 0
for (sample in missing_samples){
  if (c == 0){
    regex_pattern <- sample
  } else{
    regex_pattern <- paste(regex_pattern, sample, sep = "|")
  }
  c <- c + 1
}
# If there are difference between the samples ID, some samples are missing
if (no_missing_samples > 0){
  paste("There are", no_missing_samples, "missing samples")
  # Check if the missing samples are present in the vcf files
  if (length(grep(regex_pattern, vcf_files)) == 0){
    print("Vcf files missing")  
    print(missing_samples)
  # If not, the missing samples are in the metadata file  
  } else{
    paste("Metadata missing samples:", vcf_files[grep(regex_pattern, vcf_files)])
  } 
# If there are no difference, there are no missing file
} else{
    paste("There are no missing files in metadata")
}

#### 2.4 Load sequencing dates and TMB
print("Sequencing dates and TMB")

TMB <- read.csv(TBM_path)
sequencing_dates <- read_csv(dates_path, col_types = "cc")
colnames(sequencing_dates) <- c("Date_yymmdd", "Sample")
sequencing_dates %>%
  distinct() %>% 
  mutate(Sample_number = gsub("\\D", "", Sample)) %>%
  arrange(Sample_number) %>%
  select(Sample, Date_yymmdd) %>%
  mutate(Date_yymmdd = as.numeric(Date_yymmdd)) -> sequencing_dates

metadata_dates_TMB <- merge(sequencing_dates, TMB, by = "Sample", all = TRUE) %>%
  filter(Sample %in% metadata_sorted_by_samplesID$SampleID)

# Check if the samples_names of the dates and TBM match with the original samples names
extra_missing_samples <- setdiff(sample_names, metadata_dates_TMB$Sample)
no_extra_missing_samples <- length(extra_missing_samples)
if (no_extra_missing_samples > 0){
  paste("There are", no_extra_missing_samples, "different samples in the sequencing dates and TBM")
} else{
  paste("The sample names in sequencing dates and TBM files match the vcf files and the original metadata")
}

## Data exploratory analysis
discrete_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                      "#66A61E", "504", "#A6761D", "#E6AB02",
                      "lightskyblue", "#666666", "#0072B2", "#CC79A7",
                      "darkseagreen2", "chartreuse")
save_plot <- function(path_in_results, plot, width = 20, height = 12, dpi = 300){
  ggsave(paste(results_path, 
               path_in_results, 
               sep = ""), 
         plot,
         width = width,
         height = height,
         dpi = dpi) 
}

if (use_case_samples_only == FALSE){
  # No samples
  temp_df <- rbind(metadata_sorted_by_samplesID, 
                   transform(metadata_sorted_by_samplesID, Group = "Total"))
  temp_df$Group[temp_df$Group == "Case"] <- "Treated"
  temp_df$Group[temp_df$Group == "Benchmark"] <- "Validation"
  
  count(temp_df$Group[temp_df$Group != "Control"]) %>%
    arrange(freq) %>%
    mutate(x = factor(x,
                      levels = x)) %>%
    ggplot(aes(x = x, y = freq, fill = x)) +
    geom_bar(color = "black", stat = "identity") +
    labs(title = "Number of samples in the three datasets",
         fill = "Cancer type") +
    xlab("Cancer type") + 
    ylab("Frequency") +
    scale_fill_manual(values = c("YellowGreen", "coral", "deepskyblue2")) +
    scale_y_continuous(limits = c(0, 113), breaks = scales::pretty_breaks(n = 15)) +
    theme_bw() +
    theme(legend.position='none') +
    coord_flip() -> no_samples1
  
  count(metadata_sorted_by_samplesID$Subtype) %>%
    arrange(freq) %>%
    mutate(x = factor(x,
                      levels = x)) %>%
    ggplot(aes(x = x, y = freq, fill = x)) +
    geom_bar(color = "black", stat = "identity") +
    labs(title = "Complete cohort",
         fill = "Cancer type") +
    xlab("Cancer type") + 
    ylab("Frequency") +
    scale_fill_manual(values = discrete_palette[c(10, 4, 13, 1,
                                                  7, 8, 3, 6,
                                                  2, 5, 9, 12, 
                                                  11)]) +
    scale_y_continuous(limits = c(0, 34),
                       breaks = scales::pretty_breaks(n = 20)) +
    theme_bw() +
    theme(legend.position='bottom') +
    coord_flip() -> no_samples2
  
  ## RELEVANT
  count(metadata_sorted_by_samplesID$Subtype) %>%
    arrange(freq) %>%
    mutate(x = factor(x,
                      levels = x)) %>%
    ggplot(aes(x = x, y = freq, fill = x)) +
    geom_bar(color = "black", stat = "identity") +
    labs(title = "Complete cohort cancer types",
         fill = "Cancer type") +
    xlab("Cancer type") + 
    ylab("Frequency") +
    scale_fill_manual(values = discrete_palette[c(10, 4, 13, 1,
                                                  7, 8, 3, 6,
                                                  2, 5, 9, 12, 
                                                  11)]) +
    scale_y_continuous(limits = c(0, 34),
                       breaks = scales::pretty_breaks(n = 20)) +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5, size = 25)) +
    coord_flip() -> relevant_data1
  
  count(metadata_sorted_by_samplesID$Subtype[metadata_sorted_by_samplesID$Group == "Case"]) %>%
    arrange(freq) %>%
    mutate(x = factor(x,
                      levels = x)) %>%
    ggplot(aes(x = x, y = freq, fill = x)) +
    geom_bar(color = "black", stat = "identity") +
    labs(title = "Case subgroup cancer types",
         fill = "Cancer type") +
    xlab("Cancer type") + 
    ylab("Frequency") +
    scale_fill_manual(values = discrete_palette[c(3, 8, 10, 6, 
                                                  12, 11)]) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    theme_bw() +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5, size = 25)) +
    coord_flip() -> relevant_data2
  save_plot("relevant_results/data1.png", relevant_data1,
            width = 10, height = 5)
  save_plot("relevant_results/data2.png", relevant_data2,
            width = 7, height = 3)
  ###
  
  count(metadata_sorted_by_samplesID$Subtype[metadata_sorted_by_samplesID$Group == "Case"]) %>%
    arrange(freq) %>%
    mutate(x = factor(x,
                      levels = x)) %>%
    ggplot(aes(x = x, y = freq, fill = x)) +
    geom_bar(color = "black", stat = "identity") +
    labs(title = "Treatment group",
         fill = "Cancer type") +
    xlab("Cancer type") + 
    ylab("Frequency") +
    scale_fill_manual(values = discrete_palette[c(3, 8, 10, 6, 
                                                  12, 11)]) +
    scale_y_continuous(limits = c(0, 34),
                       breaks = scales::pretty_breaks(n = 20)) +
    theme_bw() +
    theme(legend.position='none') +
    coord_flip() -> no_samples3
  
  count(metadata_sorted_by_samplesID$Subtype[metadata_sorted_by_samplesID$Group == "Benchmark"]) %>%
    arrange(freq) %>%
    mutate(x = factor(x,
                      levels = x)) %>%
    ggplot(aes(x = x, y = freq, fill = x)) +
    geom_bar(color = "black", stat = "identity") +
    labs(title = "Validation group",
         fill = "Cancer type") +
    xlab("Cancer type") + 
    ylab("Frequency") +
    scale_y_continuous(limits = c(0, 34),
                       breaks = scales::pretty_breaks(n = 20)) +
    scale_fill_manual(values = discrete_palette[c(3, 12, 4, 7)])  +
    theme_bw() +
    theme(legend.position='none') +
    coord_flip() -> no_samples4
  
  metadata_sorted_by_samplesID %>%
    filter(!is.na(Age) & !is.na(Gender) & Gender != "NA") %>%
    ggplot(aes(x = Gender, fill = Gender)) +
    geom_bar(color = "black") +
    labs(title = "Number of samples by gender") +
    ylab("Frequency") +
    theme_bw() +
    theme(legend.position='none') +
    scale_y_continuous(limits = c(0, 113), breaks = scales::pretty_breaks(n = 15)) +
    coord_flip() -> no_samples5
  
  no_samples3_4 <- grid.arrange(no_samples3, no_samples4, ncol = 1)
  no_samples2_3_4 <- grid.arrange(no_samples2, no_samples3_4, ncol = 1)
  no_samples1_5 <- grid.arrange(no_samples1, no_samples5, ncol = 1, heights = 3:2)
  no_samples2_3_4_1 <- grid.arrange(no_samples1_5, no_samples2_3_4, ncol = 2,
                                    widths = 2:3,
                                    top = textGrob("Number of samples", 
                                                   gp = gpar(fontsize = 20)))
  
  save_plot("no_samples.png", no_samples2_3_4_1)
  
  # Age
  metadata_sorted_by_samplesID %>% 
    filter(!is.na(Age)) %>% 
    arrange(Age) %>%
    mutate(SampleID = factor(SampleID,
                             levels = SampleID)) %>%
    ggplot(aes(x = SampleID, y = Age, fill = Subtype)) +
    geom_bar(colour="black", stat = "identity") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust=1)) + 
    scale_fill_manual(values = discrete_palette[c(1:3, 5:6,
                                                  8:13)]) +
    labs(title = "Age by sample ID",
         fill = "Cancer type") +
    xlab("Sample ID") + 
    coord_cartesian(ylim=c(30, 80)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) -> age1
  
  metadata_sorted_by_samplesID %>% 
    filter(!is.na(Age)) %>%
    ggplot(aes(x = Subtype, y = Age, fill = Subtype)) +
    geom_boxplot(colour="black") +
    theme_bw() + 
    xlab("Cancer type") +
    scale_fill_manual(values = discrete_palette[c(1:3, 5:6,
                                                  8:13)]) +
    labs(title = "Age by cancer type") +
    theme(legend.position='none') +
    coord_cartesian(ylim=c(30, 80)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) -> age2
  
  metadata_sorted_by_samplesID %>% 
    filter(!is.na(Age) & !is.na(Gender) & Gender != "NA") %>%
    ggplot(aes(x = Gender, y = Age, fill = Gender)) +
    geom_boxplot(colour="black") +
    theme_bw() +
    labs(title = "Age by gender",
         fill = "Gender") +
    theme(legend.position='none') +
    ylab("") +
    coord_cartesian(ylim=c(30, 80)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) -> age3
  
  temp_df %>% 
    filter(!is.na(Age) & Group != "Control" & Group != "Validation") %>%
    ggplot(aes(x = Group, y = Age, fill = Group)) +
    geom_boxplot(colour = "black")  +
    labs(title = "Age total and treated patients",
         fill = "Cancer type") +
    theme_bw() + 
    theme(legend.position='none') +
    ylab("") +
    scale_fill_manual(values = c("deepskyblue2", "coral", "YellowGreen")) +
    coord_cartesian(ylim=c(30, 80)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) -> age4
  
  age3_4 <- grid.arrange(age3, age4, ncol = 2)
  age2_3_4 <- grid.arrange(age2, age3_4, ncol = 2, widths = 3:2)
  age_plots <- grid.arrange(age1, age2_3_4, ncol = 1,
                            top = textGrob("Patients age", 
                                           gp = gpar(fontsize = 20)))
  save_plot("age_samples.png", age_plots)
  
  # TMB
  metadata_dates_TMB %>%
    mutate(Subtype = metadata_sorted_by_samplesID$Subtype,
           Gender = metadata_sorted_by_samplesID$Gender,
           Group = metadata_sorted_by_samplesID$Group,
           Age = metadata_sorted_by_samplesID$Age) -> metadata_dates_TMB 
  
  metadata_dates_TMB %>% 
    arrange(TMB) %>%
    mutate(Sample = factor(Sample,
                        levels = Sample)) %>%
    ggplot(aes(x = Sample, y = TMB, fill = Subtype)) +
    geom_bar(colour="black", stat = "identity") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust=1)) + 
    scale_fill_manual(values = discrete_palette) +
    labs(title = "TMB by sample ID",
         fill = "Cancer type") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 20)) -> tmb1
  
  metadata_dates_TMB %>%
    ggplot(aes(x = Subtype, y = TMB, fill = Subtype)) +
    geom_boxplot(colour="black") +
    theme_bw() + 
    xlab("Cancer type") +
    scale_fill_manual(values = discrete_palette) +
    labs(title = "TMB by cancer type") +
    theme(legend.position='none') +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) -> tmb2
  
  metadata_dates_TMB %>% 
    filter(!is.na(Age) & !is.na(Gender) & Gender != "NA") %>%
    ggplot(aes(x = Gender, y = TMB, fill = Gender)) +
    geom_boxplot(colour="black") +
    theme_bw() +
    labs(title = "TMB by gender",
         fill = "Gender") +
    ylab("") +
    theme(legend.position='none') +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) -> tmb3
  
  temp_df2 <- rbind(metadata_dates_TMB, 
                   transform(metadata_dates_TMB, Group = "Total"))
  temp_df2$Group[temp_df2$Group == "Case"] <- "Treated"
  temp_df2$Group[temp_df2$Group == "Benchmark"] <- "Validation"
  
  temp_df2 %>% 
    filter(Group != "Control") %>%
    ggplot(aes(x = Group, y = TMB, fill = Group)) +
    geom_boxplot(colour = "black")  +
    labs(title = "TMB by groups",
         fill = "Cancer type") +
    theme_bw() + 
    ylab("") +
    theme(legend.position='none') +
    scale_fill_manual(values = c("deepskyblue2", "coral", "YellowGreen")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) -> tmb4
  
  tmb3_4 <- grid.arrange(tmb3, tmb4, ncol = 2, widths = 2:3)
  tmb2_3_4 <- grid.arrange(tmb2, tmb3_4, ncol = 2, widths = 3:2)
  tmb_plots <- grid.arrange(tmb1, tmb2_3_4, ncol = 1,
                             top = textGrob("Patients TMB", 
                                            gp = gpar(fontsize = 20)))
  save_plot("tmb_samples.png", tmb_plots)
  
  # Sequencing dates
  metadata_dates_TMB$Date_yymmdd <- 
    paste("20", metadata_dates_TMB$Date_yymmdd, sep = "")
  metadata_dates_TMB$Date_yymmdd <-
    as.Date(metadata_dates_TMB$Date_yymmdd, 
            format("%Y%m%d"))
  lab_dates <- pretty(metadata_dates_TMB$Date_yymmdd)
  
  metadata_dates_TMB %>% 
    mutate(Subtype = metadata_sorted_by_samplesID$Subtype) %>%
    filter(!is.na(Date_yymmdd)) %>%
    arrange(Date_yymmdd) %>%
    mutate(Sample = factor(Sample,
                           levels = Sample)) %>%
    ggplot(aes(x = Sample, y = Date_yymmdd, fill = Subtype)) +
    geom_bar(colour="black", stat = "identity") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust=1)) + 
    scale_fill_manual(values = discrete_palette[c(1:3, 5:6,
                                                  8:13)]) +
    labs(title = "Sequencing date by sample ID",
         fill = "Cancer type") +
    ylab("Sequencing date") +
    coord_cartesian(ylim = as.Date(c("2014-01-01", "2020-01-01"))) +
    scale_y_date(breaks = pretty_breaks(10)) -> date1
  
  metadata_dates_TMB %>% 
    mutate(Subtype = metadata_sorted_by_samplesID$Subtype) %>%
    filter(!is.na(Date_yymmdd)) %>%
    ggplot(aes(x = Subtype, y = Date_yymmdd, fill = Subtype)) +
    geom_boxplot(colour="black") +
    theme_bw() + 
    xlab("Cancer type") +
    ylab("Sequencing date") +
    scale_fill_manual(values = discrete_palette[c(1:3, 5:6,
                                                  8:13)]) +
    labs(title = "Sequencing date by cancer type") +
    theme(legend.position='none') +
    coord_cartesian(ylim = as.Date(c("2014-01-01", "2020-01-01"))) +
    scale_y_date(breaks = scales::pretty_breaks(n = 10)) -> date2
  
  metadata_dates_TMB %>% 
    filter(!is.na(Date_yymmdd) & !is.na(Gender) & Gender != "NA") %>%
    ggplot(aes(x = Gender, y = Date_yymmdd, fill = Gender)) +
    geom_boxplot(colour="black") +
    theme_bw() +
    labs(title = "Sequencing date by gender",
         fill = "Gender") +
    ylab("") +
    theme(legend.position='none') +
    coord_cartesian(ylim = as.Date(c("2014-01-01", "2020-01-01"))) +
    scale_y_date(breaks = scales::pretty_breaks(n = 10)) -> date3
  
  temp_df2 <- rbind(metadata_dates_TMB, 
                    transform(metadata_dates_TMB, Group = "Total"))
  temp_df2$Group[temp_df2$Group == "Case"] <- "Treated"
  temp_df2$Group[temp_df2$Group == "Benchmark"] <- "Validation"
  
  temp_df2 %>% 
    filter(Group != "Control") %>%
    ggplot(aes(x = Group, y = Date_yymmdd, fill = Group)) +
    geom_boxplot(colour = "black") +
    labs(title = "Sequencing date by groups",
         fill = "Cancer type") +
    theme_bw() + 
    ylab("") +
    theme(legend.position='none') +
    scale_fill_manual(values = c("deepskyblue2", "coral", "YellowGreen")) +
    coord_cartesian(ylim = as.Date(c("2014-01-01", "2020-01-01"))) +
    scale_y_date(breaks = scales::pretty_breaks(n = 10)) -> date4
  
  date3_4 <- grid.arrange(date3, date4, ncol = 2, widths = 2:3)
  date2_3_4 <- grid.arrange(date2, date3_4, ncol = 2, widths = 3:2)
  date_plots <- grid.arrange(date1, date2_3_4, ncol = 1,
                             top = textGrob("Patients sequencing date", 
                                            gp = gpar(fontsize = 20)))
  save_plot("date_samples.png", date_plots)
}

#metadata_sorted_by_samplesID$Subtype <- gsub(" ", "_", metadata_sorted_by_samplesID$Subtype)
#metadata_sorted_by_samplesID$Subtype <- gsub("cancer", "c" ,metadata_sorted_by_samplesID$Subtype)

#### 3. Mutation chracteristics
## Make a 96 trinucleodide mutation count matrix (96 mutational profile)
writeLines("\nMutation chracteristics")
print("Generating 96 mutational profile count matrix")
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

## Mutation spectrum
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)

# General spectrum
print("Generating spectrum")
general_spectrum <- plot_spectrum(type_occurrences) + 
  ggtitle("General mutation spectrum") + 
  theme(plot.title = element_text(hjust = 0.5))
general_spectrum_ct <- plot_spectrum(type_occurrences, CT = TRUE) +
  ggtitle("General mutation spectrum with C to T CpG distinction") + 
  theme(plot.title = element_text(hjust = 0.5))
general_spectrum_arranged <- grid.arrange(general_spectrum, general_spectrum_ct)
ggsave(paste(results_path, "COSMIC/spectrum/general_spectrum.png", sep = ""), general_spectrum_arranged)

# Spectrum divided by tumor subtypes
spectrum_by_types <- plot_spectrum(type_occurrences, by = metadata_sorted_by_samplesID$Subtype, 
                                   CT = TRUE, legend = TRUE) + 
  ggtitle("Mutation spectrum divided by tumor types") +
  theme_bw(base_size = 15) + 
  theme(plot.title = element_text(hjust = 0.5, size = 25))
save_plot("COSMIC/spectrum/spectrum_by_types.png", spectrum_by_types, width = 12, height = 6)
save_plot("relevant_results/spectrum_by_types.png", spectrum_by_types, width = 12, height = 8)

# Spectrum divided by sample (with added information of tumor tipes from metadata)
concatenate_two_vectors <- function(vector1, vector2){
  if (length(vector1) == length(vector2)){
    concatenated_vector <- c()
    for (n in seq(length(vector1))){
      concatenated_element <- (paste(vector1[n], vector2[n], sep = " - "))
      concatenated_vector <- c(concatenated_vector, concatenated_element)
    }
  } else{
    print("Error: the two vectors must be of the same length")
  }
  return(concatenated_vector)
}
sample_names_and_subtypes <- concatenate_two_vectors(sample_names, metadata_sorted_by_samplesID$Subtype)

spectrum_by_samples <- plot_spectrum(type_occurrences, by = sample_names_and_subtypes, 
                                     CT = TRUE, legend = TRUE) + 
  ggtitle("Mutation spectrum for each sample") + 
  theme(plot.title = element_text(hjust = 0.5))
save_plot("COSMIC/spectrum/spectrum_by_samples.png", spectrum_by_samples)


#### 4. Mutational signatures
writeLines("\nSignaures characteristics")

### 4.1 Load signatures sets
print("Load signatures sets")

## Exomes cancer signatures
# Load data
cancer_signatures_exomes <- read_csv(paste(signatures_path, 
                                           "sigProfiler_exome_SBS_signatures.csv", sep = ""))
# Create SomaticMutationType column to sort the mutation types
cancer_signatures_exomes %>% 
  mutate(SomaticMutationType = paste(substr(SubType, 1, 1), "[", Type, "]", substr(SubType, 3, 3), sep=""),
         .after = SubType) -> cancer_signatures_exomes
print(paste("Dimension exomes set", dim(cancer_signatures_exomes)))
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures_exomes$SomaticMutationType)
# Reorder cancer signatures dataframe
cancer_signatures_exomes = as.data.frame(cancer_signatures_exomes[as.vector(new_order),])
# Add trinucletiode changes names as row.names
row.names(cancer_signatures_exomes) = cancer_signatures_exomes$SomaticMutationType
# Keep only 96 contributions of the signatures in matrix
cancer_signatures_exomes = as.matrix(select(cancer_signatures_exomes, starts_with("SBS")))
dim(cancer_signatures_exomes)
colnames(cancer_signatures_exomes) <- gsub("\\D{3}", "S.", colnames(cancer_signatures_exomes))


### 4.2 Signatures characteristics
writeLines("\nSignatures characteristics")

## Original mutation profile
# Functions to save 96 profile plot divided in batches 
first_capital_letter <- function(word){
  word <- tolower(word)
  substr(word, 1, 1) <- toupper(substr(word, 1, 1))
  return(word)
}
save_96_profile_plots <- function(cancer_signatures, plot_name, old_signatures=FALSE){
  # Define the indexes intervals for the different plots
  if (old_signatures == FALSE){
    signatures_batch_intervals <- list(1:20, 21:40, 41:dim(cancer_signatures)[2])
  }
  else{
    signatures_batch_intervals <- list(1:dim(cancer_signatures)[2])
  }
  # Save a plot for each batch
  print(paste("Saving", plot_name))
  for (batch_interval in signatures_batch_intervals){
    # Create plot object
    name <- paste(plot_name, "_", batch_interval[1], "_to_", batch_interval[length(batch_interval)], sep = "")
    plot_96_profile <- plot_96_profile(cancer_signatures[, batch_interval], condensed = TRUE, ymax = 0.3) + 
      ggtitle(gsub("_", " ", first_capital_letter(name))) + 
      theme(plot.title = element_text(hjust = 0.5))
    # Save plot
    save_plot(path_in_results = paste("COSMIC/96_profile_plots/", name, ".png", sep=""), plot = plot_96_profile)
  }
}
# Save to server
save_96_profile_plots(cancer_signatures_exomes, "exomes_signatures_96_profile")
# Plot most interesting signatures
top_selected_sig <- c("S.1", "S.3", "S.4", "S.5", "S.6", "S.7a", "S.7b", 
                      "S.11", "S.13", "S.15", "S.18", "S.25", "S.29",
                      "S.30", "S.35", "S.36", "S.38", "S.40", "S.43", 
                      "S.44", "S.45","S.55", "S.56", "S.60") 

plot_96_profile(cancer_signatures_exomes[,top_selected_sig], condensed = TRUE, ymax = 0.3) +
  ggtitle("Most relevant exomes COSMIC signatures 96 profiles") +
  theme(plot.title = element_text(hjust = 0.5, size = 25)) -> plot_96_profile_selected
save_plot(path_in_results = paste("relevant_results/96_profile_selected.png", sep=""), 
          plot = plot_96_profile_selected, height = 14, width= 9)


## Hierarchically clustering
print("Running signatures hierarchical clustering")
# Exomes
hclust_SBS_exomes <- cluster_signatures(cancer_signatures_exomes, method = "average")
cosmic_order_exomes = colnames(cancer_signatures_exomes)[hclust_SBS_exomes$order]
save_plot("COSMIC/hierarchical_clustering_signatures/hierarchical_clustering_exomes_signatures.png",
          plot(hclust_SBS_exomes))


#### 5. Pairwise similarity between samples original mutational profiles and COSMIC signatures 
writeLines("\nSimilarity between original mutational profiles and COSMIC sets (exomes, whole genomes, old set)")

## PCA 
print("Generating PCA plots")
# Signatures as input, color as signature average cosine similarity across samples
pca_with_metadata_signatures <- function(matrix,
                                         title = "",
                                         added_col = "average"){
  labels = colnames(matrix)
  # Transpose and extract average cosine similarity for each signature
  if (added_col == "average"){
    matrix %>% t -> matrix_t
    matrix_t %>% as_tibble() %>% 
      mutate(added_col = apply(., 1, mean)) %>%
      select(added_col) -> added_column
    fill_legend = "Average"
  } 
  if (added_col == "variance"){
    matrix %>% t -> matrix_t
    matrix_t %>% as_tibble() %>% 
      mutate(added_col = apply(., 1, var)) %>%
      select(added_col) -> added_column
    fill_legend = "Variance"
  }
  # PCA
  pca <- prcomp(matrix_t,
                center = TRUE)
  percent_variance <- summary(pca)$importance["Proportion of Variance",] * 100
  # Plot
  as_tibble(pca$x) %>% 
    ggplot(aes(x=PC1, y=PC2, 
               label = colnames(matrix),
               fill = added_column$added_col)) + 
    geom_point() + 
    geom_label_repel(segment.alpha = 0.5) +
    labs(title = title,
         fill = fill_legend) +
    xlab(label = paste("PC1", percent_variance[1], "%")) +
    ylab(label = paste("PC2", percent_variance[2], "%")) + 
    scale_fill_distiller(palette = "Spectral") +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5, size = 25))
}
# Samples as input, addition of metadata information
pca_with_metadata_samples <- function(matrix, 
                                      title = "", 
                                      colored_by = "tumor", 
                                      labels = sample_names){
  # Assign additional metadata information
  if (colored_by == "tumor"){
    metadata_for_color = metadata_sorted_by_samplesID$Subtype 
    label_color = "Tumor type"
    scale_fill = scale_fill_manual(values = discrete_palette)
  } 
  if (colored_by == "age"){
    metadata_for_color = metadata_sorted_by_samplesID$Age
    label_color = "Age"
    scale_fill = scale_fill_viridis() 
  }
  if (colored_by == "TMB"){
    metadata_for_color = metadata_dates_TMB$TMB
    label_color = "TMB"
    scale_fill = scale_fill_distiller(palette = "Spectral") 
  }
  if (colored_by == "seq_date"){
    metadata_for_color = as.numeric(metadata_dates_TMB$Date_yymmdd)
    label_color = "Sequencing date"
    scale_fill =  scale_fill_viridis(breaks = as.numeric(lab_dates), 
                                     labels = lab_dates)
  }
  if (colored_by == "response"){
    metadata_for_color = metadata_sorted_by_samplesID$Response
    label_color = "Response"
    scale_fill = scale_fill_manual(values = c("dodgerblue", "coral", "darkgoldenrod1")) 
  }
  if (colored_by == "group"){
    metadata_for_color = metadata_sorted_by_samplesID$Group
    label_color = "Group"
    scale_fill = scale_fill_manual(values = discrete_palette)
  }
  
  # PCA
  pca <- prcomp(matrix,
                center = TRUE)
  percent_variance <- summary(pca)$importance["Proportion of Variance",] * 100
  # Plot
  as_tibble(pca$x) %>% 
    ggplot(aes(x=PC1, y=PC2, 
               label = labels,
               fill = metadata_for_color,
               shape = metadata_sorted_by_samplesID$Gender)) + 
    geom_point() + 
    scale_shape(solid = TRUE) +
    geom_label_repel(segment.alpha = 0.5) +
    labs(title = title,
         fill = label_color,
         shape = "Gender") +
    xlab(label = paste("PC1", percent_variance[1], "%")) +
    ylab(label = paste("PC2", percent_variance[2], "%")) + 
    scale_fill +
    scale_shape_manual(values = c(21,24,23)) +
    scale_size(guide="none") +
    theme_bw(base_size = 15) + 
    theme(plot.title = element_text(hjust = 0.5, size = 25))
}
# Extract pairwise cosine similarity
cos_sim_samples_signatures_exomes = cos_sim_matrix(mut_mat, cancer_signatures_exomes)

pca_plot_cos_sim_exomes_sig_input <- 
  pca_with_metadata_signatures(cos_sim_samples_signatures_exomes,
                               title = "PCA plot of pairwise cosine similarity of original profile VS exomes COSMIC signatures (signatures as input, colored by average cosine similarity)")

pca_plot_cos_sim_exomes_samples_input_col_subtype <- 
  pca_with_metadata_samples(cos_sim_samples_signatures_exomes,
                            "PCA plot of pairwise cosine similarity of sample original profile VS exomes COSMIC signatures (samples as input, colored by tumor type)")

pca_plot_cos_sim_exomes_samples_input_col_age_labeled_subtypes <- 
  pca_with_metadata_samples(cos_sim_samples_signatures_exomes,
                            "PCA plot of pairwise cosine similarity of sample original profile VS exomes COSMIC signatures (samples as input, colored by patient age, labeled by tumor type)",
                            colored_by = "age",
                            labels = metadata_sorted_by_samplesID$Subtype)

if (use_case_samples_only == FALSE){
  pca_plot_cos_sim_exomes_samples_input_col_TMB_labeled_subtypes <- 
    pca_with_metadata_samples(cos_sim_samples_signatures_exomes,
                              "PCA plot of pairwise cosine similarity of sample original profile VS exomes COSMIC signatures (samples as input, colored by TMB, labeled by tumor type)",
                              colored_by = "TMB",
                              labels = metadata_sorted_by_samplesID$Subtype)
  
  pca_plot_cos_sim_exomes_samples_input_col_date_labeled_subtypes <- 
    pca_with_metadata_samples(cos_sim_samples_signatures_exomes,
                              "PCA plot of pairwise cosine similarity of sample original profile VS exomes COSMIC signatures (samples as input, colored by sequencing date, labeled by tumor type)",
                              colored_by = "seq_date",
                              labels = metadata_sorted_by_samplesID$Subtype)
  pca_plot_cos_sim_exomes_samples_input_col_group_labeled_subtype <- 
    pca_with_metadata_samples(cos_sim_samples_signatures_exomes,
                              "PCA plot of pairwise cosine similarity of sample original profile VS exomes COSMIC signatures (samples as input, colored by group, labeled by tumor type)",
                              colored_by = "group",
                              labels = metadata_sorted_by_samplesID$Subtype)
} else{
  pca_plot_cos_sim_exomes_samples_input_col_subtype_labeled_cancergroup <- 
    pca_with_metadata_samples(cos_sim_samples_signatures_exomes,
                              "PCA plot of pairwise cosine similarity of sample original profile VS exomes COSMIC signatures (samples as input, colored by tumor type, labeled by cancer group)",
                              labels = metadata_sorted_by_samplesID$Cancer_group)
  pca_plot_cos_sim_exomes_samples_input_col_subtype_labeled_HRD <- 
    pca_with_metadata_samples(cos_sim_samples_signatures_exomes,
                              "PCA plot of pairwise cosine similarity of sample original profile VS exomes COSMIC signatures (samples as input, colored by tumor type, labeled by somatic HRD mutation)",
                              labels = metadata_sorted_by_samplesID$Somatic_HRD_mt)
  pca_plot_cos_sim_exomes_samples_input_col_response_labeled_cancergroup <- 
    pca_with_metadata_samples(cos_sim_samples_signatures_exomes,
                              "PCA plot of pairwise cosine similarity of sample original profile VS exomes COSMIC signatures (samples as input, colored by response, labeled by tumor type)",
                              labels = metadata_sorted_by_samplesID$Subtype,
                              colored_by = "response")
}
# Save
print("Saving exomes set PCA plots")
save_plot("COSMIC/cosine_similarity_original_vs_cosmic/pca/pca_plot_cos_sim_exomes_sig_input.png", 
          pca_plot_cos_sim_exomes_sig_input)
save_plot("COSMIC/cosine_similarity_original_vs_cosmic/pca/pca_plot_cos_sim_exomes_samples_input_col_subtype.png", 
          pca_plot_cos_sim_exomes_samples_input_col_subtype)
save_plot("COSMIC/cosine_similarity_original_vs_cosmic/pca/pca_plot_cos_sim_exomes_samples_input_col_age_labeled_subtypes.png", 
          pca_plot_cos_sim_exomes_samples_input_col_age_labeled_subtypes)
if (use_case_samples_only == FALSE){
  save_plot("COSMIC/cosine_similarity_original_vs_cosmic/pca/pca_plot_cos_sim_exomes_samples_input_col_TMB_labeled_subtypes.png", 
            pca_plot_cos_sim_exomes_samples_input_col_TMB_labeled_subtypes)
  save_plot("COSMIC/cosine_similarity_original_vs_cosmic/pca/pca_plot_cos_sim_exomes_samples_input_col_date_labeled_subtypes.png", 
            pca_plot_cos_sim_exomes_samples_input_col_date_labeled_subtypes)
  save_plot("COSMIC/cosine_similarity_original_vs_cosmic/pca/pca_plot_cos_sim_exomes_samples_input_col_group_labeled_subtype.png", 
            pca_plot_cos_sim_exomes_samples_input_col_group_labeled_subtype)
} else {
  save_plot("COSMIC/cosine_similarity_original_vs_cosmic/pca/pca_plot_cos_sim_exomes_samples_input_col_subtype_labeled_cancergroup.png", 
            pca_plot_cos_sim_exomes_samples_input_col_subtype_labeled_cancergroup)
  save_plot("COSMIC/cosine_similarity_original_vs_cosmic/pca/pca_plot_cos_sim_exomes_samples_input_col_subtype_labeled_HRD.png", 
            pca_plot_cos_sim_exomes_samples_input_col_subtype_labeled_HRD)
  save_plot("COSMIC/cosine_similarity_original_vs_cosmic/pca/pca_plot_cos_sim_exomes_samples_input_col_response_labeled_cancergroup.png", 
            pca_plot_cos_sim_exomes_samples_input_col_response_labeled_cancergroup)
}

## PCA benchmark data

if (use_case_samples_only == FALSE){
  benchmark_samples <- metadata_sorted_by_samplesID %>% 
    filter(Group == "Benchmark") 
  benchmark_samples$Annotation[benchmark_samples$Annotation == "HRD/BRCAness profile"] <- "HRD/BRCA p."
  
  benchmark_sample_index <- rownames(cos_sim_samples_signatures_exomes) %in% benchmark_samples$SampleID
  cos_sim_benchmark_signatures_exome <- cos_sim_samples_signatures_exomes[benchmark_sample_index,]
  
  pca_benchmark <- function(matrix,
                            label = "annotation",
                            title = ""){
    pca <- prcomp(matrix,
                  center = TRUE)
    percent_variance <- summary(pca)$importance["Proportion of Variance",] * 100
    
    if (label == "annotation"){
      labels = benchmark_samples$Annotation
    }
    if (label == "sample"){
      labels = benchmark_samples$SampleID
    }
    # Plot
    as_tibble(pca$x) %>% 
      ggplot(aes(x=PC1, y=PC2, 
                 label = labels,
                 fill = benchmark_samples$Subtype)) + 
      geom_point() +
      geom_label_repel(segment.alpha = 0.5) +
      xlab(label = paste("PC1", percent_variance[1], "%")) +
      ylab(label = paste("PC2", percent_variance[2], "%")) + 
      scale_fill_manual(values = c(discrete_palette[3], discrete_palette[4], 
                                   discrete_palette[7], discrete_palette[12])) +
      labs(title = title,
           fill = "Tumor type") +
      theme_bw(base_size = 15) + 
      theme(plot.title = element_text(hjust = 0.34, size = 25))
  }
  pca_plot_cos_sim_exomes_samples_input_benchmark <- 
    pca_benchmark(cos_sim_benchmark_signatures_exome, 
                  title = "PCA plot of pairwise cosine similarity of original profile VS exomes COSMIC signatures (signatures as input, benchmark data, labeled by sample ID)",
                  label = "sample")
  save_plot("benchmark/COSMIC/cosine_similarity_original_vs_cosmic/pca_cos_sim_exomes_samples_input_benchmark.png", 
            pca_plot_cos_sim_exomes_samples_input_benchmark)
  
  pca_plot_cos_sim_exomes_samples_input_benchmark_annotation <- 
    pca_benchmark(cos_sim_benchmark_signatures_exome, 
                  title = "PCA plot of pairwise cosine similarity of original profile VS exomes COSMIC signatures (signatures as input, benchmark data, labeled by annotation)",
                  label = "annotation")
  save_plot("benchmark/COSMIC/cosine_similarity_original_vs_cosmic/pca_cos_sim_exomes_samples_input_benchmark_annotation.png", 
            pca_plot_cos_sim_exomes_samples_input_benchmark_annotation)
  
  ### Relevant result
  pca_plot_cos_sim_exomes_sig_input_relevant <- 
    pca_with_metadata_signatures(cos_sim_samples_signatures_exomes,
                                 title = "PCA plot of pairwise cosine similarity (signatures as input)")
  save_plot("relevant_results/pca_plot_cos_sim_signatures.png", 
            pca_plot_cos_sim_exomes_sig_input_relevant, 
            height = 7, width = 12)
  pca_plot_cos_sim_exomes_samples_input_col_date_labeled_subtypes_relevant <- 
    pca_with_metadata_samples(cos_sim_samples_signatures_exomes,
                              "PCA plot of pairwise cosine similarity (samples as input)",
                              colored_by = "seq_date",
                              labels = metadata_sorted_by_samplesID$Subtype)
  save_plot("relevant_results/pca_plot_cos_sim_samples_date.png", 
            pca_plot_cos_sim_exomes_samples_input_col_date_labeled_subtypes_relevant)
  pca_plot_cos_sim_exomes_samples_input_benchmark_annotation_relevant <- 
    pca_benchmark(cos_sim_benchmark_signatures_exome, 
                  title = "PCA plot of pairwise cosine similarity (annotated samples)",
                  label = "annotation")
  save_plot("relevant_results/pca_plot_cos_sim_benchmark.png", 
            pca_plot_cos_sim_exomes_samples_input_benchmark_annotation_relevant,
            height = 6, width = 9)
  ###
}

## Relevant
if (use_case_samples_only == TRUE){
  pca_plot_cos_sim_exomes_samples_input_col_response_labeled_cancergroup_relevant <- 
    pca_with_metadata_samples(cos_sim_samples_signatures_exomes,
                              "PCA plot of pairwise cosine similarity (case subgroup)",
                              labels = metadata_sorted_by_samplesID$Subtype,
                              colored_by = "response")
  save_plot("../results_all_samples/relevant_results/pca_plot_cos_sim_case_large.png", 
            pca_plot_cos_sim_exomes_samples_input_col_response_labeled_cancergroup_relevant,
            height = 7, width = 10)
}
##


## PCA biplots
# (could use the reldist package [not installed] to use the Gini coefficient to 
# find samples specific signatures of specificity)
print("Generating PCA biplots")

# PCA biplot with all signatures

# Original plot + increase labels size and space from line
plot_pca_biplot <- function(matrix, title = "", 
                            varname_adjust = 5, 
                            point_size = 1,
                            arrows = TRUE,
                            groups = metadata_sorted_by_samplesID$Subtype){
  pca <- prcomp(matrix)
  plot <- ggbiplot(pca, obs.scale = 1, var.scale = 1, 
                   circle = F, varname.size = 4, 
                   varname.adjust = varname_adjust, ellipse = TRUE,
                   groups = groups,
                   var.axes = arrows) +
    geom_point(aes(color = groups), size = point_size) +
    scale_color_manual(values = discrete_palette) +
    theme_bw() +
    labs(title = title,
         color = "Tumor type",
         shape = "Gender") 
  return(plot)
}

biplot_cos_sim_original_vs_exomes <- plot_pca_biplot(cos_sim_samples_signatures_exomes, 
                                                     "PCA biplot of pairwise cosine similarity between original samples profiles and exomes COSMIC signatures",
                                                     varname_adjust = 13,
                                                     arrows = FALSE,
                                                     point_size = 2)
save_plot("COSMIC/cosine_similarity_original_vs_cosmic/pca_biplots/pca_biplot_cos_sim_original_vs_exomes.png", 
          biplot_cos_sim_original_vs_exomes)

# Extract signatures with variance larger than threshold for each set (for biplot)
if (use_case_samples_only == FALSE){
  variance_threshold <- 0.015
} else{
  variance_threshold <- 0.009 
}

get_large_var_signatures <- function(cos_sim_samples_signatures, threshold){
  # Transpose matrix to have signature as rows and samples as columns
  cos_sim_samples_signatures %>% t -> cos_sim_samples_signatures_t
  # Add variance to the matrix and filter for signature with variance larger than threshold
  rownames_to_column(as.data.frame(cos_sim_samples_signatures_t), 
                     var = "Signature") %>%
    as_tibble() %>%
    mutate(Variance = apply(cos_sim_samples_signatures_t[,2:ncol(cos_sim_samples_signatures_t)], 
                            1, var)) %>%
    filter(Variance > threshold) -> signatures_var_larger_threshold
  # Return the signatures vector
  return(signatures_var_larger_threshold$Signature)
}

large_var_signatures_cos_sim_original_vs_exomes <- 
  get_large_var_signatures(cos_sim_samples_signatures_exomes, variance_threshold)

# PCA biplots with signatures with variance larger than threshold
biplot_cos_sim_original_vs_exomes_large_var <- 
  plot_pca_biplot(cos_sim_samples_signatures_exomes[,large_var_signatures_cos_sim_original_vs_exomes], 
                  "PCA biplot of pairwise cosine similarity between original samples profiles and exomes COSMIC signatures with large variance",
                  varname_adjust = 7)
save_plot("COSMIC/cosine_similarity_original_vs_cosmic/pca_biplots/pca_biplot_cos_sim_original_vs_exomes_large_var.png", 
          biplot_cos_sim_original_vs_exomes_large_var)


### Heatmaps 
print("Generating heatmaps")

## Heatmap with MutationalPattern
mp_heatmap_cos_sim_original_vs_exomes_cosmic <- 
  plot_cosine_heatmap(cos_sim_samples_signatures_exomes, col_order = cosmic_order_exomes, cluster_rows = TRUE)
ggsave(paste(results_path, "COSMIC/cosine_similarity_original_vs_cosmic/heatmaps/mp_heatmap_cos_sim_original_vs_exomes_cosmic.png", sep = ""),
       mp_heatmap_cos_sim_original_vs_exomes_cosmic)

# Define rownames for samples cases
rownames_names_HRD <- concatenate_two_vectors(sample_names, 
                                              metadata_sorted_by_samplesID$Somatic_HRD_mt)
rownames_names_HRD <- gsub("F", "(F", rownames_names_HRD)
rownames_names_HRD <- gsub(" - ", ")\n", rownames_names_HRD)
rownames_names_HRD <- gsub(", ", ",\n", rownames_names_HRD)
rownames_names_HRD <- gsub(" ", "\n", rownames_names_HRD)

# Clip TMB values larger than 350 to the largest value below 300 and set minimum to 0 (only for heatmap visualization)
#tmb_index_larger350 = which(metadata_dates_TMB$TMB > 350)
#min_tmb_index = which.min(metadata_dates_TMB$TMB)
#metadata_dates_TMB$TMB[tmb_index_larger350] = max(metadata_dates_TMB$TMB[-tmb_index_larger350])
#metadata_dates_TMB$TMB[min_tmb_index] = 0

days_after_first_seq_sample <- as.numeric(metadata_dates_TMB$Date_yymmdd) -
  min(as.numeric(metadata_dates_TMB$Date_yymmdd[!is.na(metadata_dates_TMB$Date_yymmdd)])) 
# Pheatmap custom function
generate_pheatmap <- function(matrix, 
                              title = "", 
                              metadata_annotation = FALSE,
                              show_response = FALSE,
                              no_clusters = NA,
                              fontsize = 7,
                              cluster_cols = TRUE){
  # Select what to show on the annotation row (column)
  if (metadata_annotation == FALSE){
    if (use_case_samples_only == TRUE){
      if (show_response == TRUE){
        pheatmap_annotation <- data.frame("Tumor.type" = metadata_sorted_by_samplesID$Subtype,
                                          "Response" = metadata_sorted_by_samplesID$Response)  
      } else{
        pheatmap_annotation <- data.frame("Tumor.type" = metadata_sorted_by_samplesID$Subtype)
      }
      rownames(pheatmap_annotation) <- rownames_names_HRD
      rownames(matrix) <- rownames_names_HRD
    } else{
      pheatmap_annotation <- data.frame("Tumor.type" = metadata_sorted_by_samplesID$Subtype,
                                        "Group" = metadata_sorted_by_samplesID$Group) 
      rownames(pheatmap_annotation) <- sample_names
    }
  } 
  if (metadata_annotation == TRUE){
    if (use_case_samples_only == TRUE){
      pheatmap_annotation <- data.frame("Tumor.type" = metadata_sorted_by_samplesID$Subtype,
                                        "LOH" = metadata_sorted_by_samplesID$LOH_larger_15MB,
                                        "TMB" = metadata_dates_TMB$TMB,
                                        "Response" = metadata_sorted_by_samplesID$Response)
      rownames(pheatmap_annotation) <- rownames_names_HRD
      rownames(matrix) <- rownames_names_HRD
    } else{
      pheatmap_annotation <- data.frame("Tumor.type" = metadata_sorted_by_samplesID$Subtype,
                                        "TMB" = metadata_dates_TMB$TMB,
                                        "Seq.date" = days_after_first_seq_sample,
                                        "Replicate" = metadata_sorted_by_samplesID$Replicate,
                                        "Age" = metadata_sorted_by_samplesID$Age) 
      rownames(pheatmap_annotation) <- sample_names
    }
  }
  annotation_colors = list(Tumor.type = c("Prostate c." = discrete_palette[12],
                                          "Endometrial c." = discrete_palette[6],
                                          "Ovarian c." = discrete_palette[11],
                                          "Mesothelioma" = discrete_palette[9],
                                          "Breast c." = discrete_palette[2],
                                          "Melanoma" = discrete_palette[8],
                                          "Bile duct c." = discrete_palette[1],
                                          "Cervical c." = discrete_palette[3],
                                          "Vulvovaginal c." = discrete_palette[13],
                                          "Colon c." = discrete_palette[4],
                                          "Lung c." = discrete_palette[7],
                                          "CUP" = discrete_palette[5],
                                          "Multiple c." = discrete_palette[10]))
  
  final_annotation_colors = c(Tumor.type = c())
  for (tumor in unique(pheatmap_annotation$Tumor.type)){
    final_annotation_colors$Tumor.type = c(final_annotation_colors$Tumor.type,
                                           annotation_colors$Tumor.type[tumor])
  }
  if ("Response" %in% colnames(pheatmap_annotation)){
    final_annotation_colors$Response = c("0" = "palegoldenrod",
                                         "1" = "coral",
                                         "NA" = "darkgoldenrod1")
  }
  # Generate heatmap
  plot <- pheatmap(matrix,
                   annotation_row = pheatmap_annotation,
                   annotation_colors = final_annotation_colors,
                   color = viridis(10),
                   main = title,
                   fontsize_row = fontsize,
                   cutree_rows = no_clusters,
                   cluster_cols = cluster_cols)
  # Change the matrix row to the original sample names
  rownames(matrix) <- sample_names
  return(plot)
}
if (use_case_samples_only == TRUE){
  fontsize = 7
} else{
  fontsize = 6
}

# Heatmap using metadata as annotation 
pheatmap_cos_sim_original_vs_exomes_cosmic_metadata <- 
  generate_pheatmap(cos_sim_samples_signatures_exomes,
                    "Heatmap of pairwise cosine similarity between original samples profiles and exomes COSMIC signatures (annotated by tumor type and metadata)",
                    metadata_annotation = TRUE,
                    fontsize = fontsize)
save_plot("COSMIC/cosine_similarity_original_vs_cosmic/heatmaps/pheatmap_cos_sim_original_vs_exomes_cosmic_metadata.png", 
          pheatmap_cos_sim_original_vs_exomes_cosmic_metadata)

## RELEVANT
pheatmap_cos_sim_original_vs_exomes_cosmic_metadata_relevant <- 
  generate_pheatmap(cos_sim_samples_signatures_exomes,
                    "Heatmap of pairwise cosine similarity between original samples profiles and exomes COSMIC signatures",
                    metadata_annotation = TRUE,
                    fontsize = fontsize,
                    no_clusters = 9)
save_plot("relevant_results/pheatmap_cos_sim_replicate.png", 
          pheatmap_cos_sim_original_vs_exomes_cosmic_metadata_relevant)
##

# Heatmaps for case samples
if (use_case_samples_only == TRUE){
  # Tumor types as annotation
  pheatmap_cos_sim_original_vs_exomes_cosmic_ann_tumor <- 
    generate_pheatmap(cos_sim_samples_signatures_exomes,
                      "Heatmap of pairwise cosine similarity between original samples profiles and exomes COSMIC signatures (annotated by tumor type)")
  save_plot("COSMIC/cosine_similarity_original_vs_cosmic/heatmaps/pheatmap_cos_sim_original_vs_exomes_cosmic_ann_tumor.png", 
            pheatmap_cos_sim_original_vs_exomes_cosmic_ann_tumor)
  # Tumor types and treatment response as annotation
  pheatmap_cos_sim_original_vs_exomes_cosmic_ann_tumor_response <- 
    generate_pheatmap(cos_sim_samples_signatures_exomes,
                      "Heatmap of pairwise cosine similarity between original samples profiles and exomes COSMIC signatures (annotated by tumor type and treatment response)",
                      show_response = TRUE)
  save_plot("COSMIC/cosine_similarity_original_vs_cosmic/heatmaps/pheatmap_cos_sim_original_vs_exomes_cosmic_ann_tumor_response.png", 
            pheatmap_cos_sim_original_vs_exomes_cosmic_ann_tumor_response)
}

# Generate and save n heatmaps with increasing number of clusters
pheatmaps_iterative_clustered <- function(matrix,
                                          max_clusters = 10,
                                          real_time_visualization = TRUE,
                                          measure = "similarity",
                                          metadata_annotation = TRUE,
                                          fontsize = 6,
                                          de_novo_signatures = FALSE,
                                          cluster_cols = TRUE,
                                          no_flat = FALSE){
  if(no_flat == FALSE){
    flat = ""
    flat_title = ""
  } else{
    flat_title = ("(removed flat signatures)")
    flat = "removed_flat/"
  }
  if (de_novo_signatures == TRUE){
    signatures <- "de novo extracted signatures"
    directory <- "de_novo/"
  } else{
    signatures <- "exomes COSMIC signatures"
    directory <- "COSMIC/"
  }
  if (measure == "similarity"){
    title = paste("Heatmap of pairwise cosine similarity between original samples profiles and", signatures, flat_title)  
    filename = paste(directory, "cosine_similarity_original_vs_cosmic/heatmaps/", flat, "heatmaps_clustered/pheatmap_cos_sim_original_vs_signatures", sep = "")
  }
  if (measure == "contribution"){
    title = paste("Heatmap of", signatures, "relative contribution to recostructed profiles", flat_title)
    filename = paste(directory, "sig_contributions/heatmaps/", flat, "heatmaps_clustered/pheatmap_signatures_contribution", sep = "")
  }
  for (n in seq(2, max_clusters)){
    final_title = paste(title, " (", n, " clusters)", sep = "")
    final_filename = paste(filename, "_", n,  "_clusters.png", sep = "")
    plot <- generate_pheatmap(matrix,
                              final_title,
                              metadata_annotation = metadata_annotation,
                              no_clusters = n,
                              fontsize = fontsize,
                              cluster_cols = cluster_cols)
    save_plot(final_filename, plot)
    if (real_time_visualization == TRUE){
      readline(prompt = "Press [Enter] to continue")
    }
  }
}
if (use_case_samples_only == FALSE){
  pheatmaps_iterative_clustered(cos_sim_samples_signatures_exomes,
                                real_time_visualization = FALSE,
                                metadata_annotation = TRUE)
} else{
  pheatmaps_iterative_clustered(cos_sim_samples_signatures_exomes,
                                real_time_visualization = FALSE,
                                metadata_annotation = TRUE,
                                fontsize = 7) 
}

# Heatmaps without flats mutations
cos_sim_samples_signatures_exomes %>% 
  as_tibble() %>% 
  select(-S.5, -S.3, -S.25, -S.40) -> cos_sim_samples_signatures_no_flat
cos_sim_samples_signatures_no_flat <- as.data.frame(cos_sim_samples_signatures_no_flat)
rownames(cos_sim_samples_signatures_no_flat) <- sample_names

pheatmap_cos_sim_original_vs_exomes_cosmic_metadata_noflat <- 
  generate_pheatmap(cos_sim_samples_signatures_no_flat,
                    "Heatmap of pairwise cosine similarity between original samples profiles and exomes COSMIC signatures (annotated by tumor type and metadata, no flat signatures)",
                    metadata_annotation = TRUE,
                    fontsize = fontsize)
save_plot("COSMIC/cosine_similarity_original_vs_cosmic/heatmaps/removed_flat/pheatmap_cos_sim_original_vs_exomes_cosmic_metadata_noflat.png", 
          pheatmap_cos_sim_original_vs_exomes_cosmic_metadata_noflat)

pheatmaps_iterative_clustered(cos_sim_samples_signatures_no_flat,
                              fontsize = fontsize,
                              real_time_visualization = FALSE,
                              no_flat = TRUE)

## Heatmaps benchmark data
if (use_case_samples_only == FALSE){
  
  # MutationalPattern
  benchmark_cos_sim_mpheatmap <- 
    plot_cosine_heatmap(cos_sim_benchmark_signatures_exome, col_order = cosmic_order_exomes, cluster_rows = TRUE)
  ggsave(paste(results_path, "benchmark/COSMIC/cosine_similarity_original_vs_cosmic/mpheatmap_benchmark_cos_sim.png", sep = ""),
         benchmark_cos_sim_mpheatmap)
  
  # Pheatmap
  benchmark_samples$Annotation[benchmark_samples$Annotation == "HRD/BRCAness profile"] <- "HRD/BRCA p."
  pheatmap_annotation_benchmark <- data.frame("Tumor.type" = benchmark_samples$Subtype,
                                              "Annotation" = benchmark_samples$Annotation) 
  rownames(pheatmap_annotation_benchmark) <- benchmark_samples$SampleID
  annotation_colors_benchmark = list(Tumor.type = c("Prostate c." = discrete_palette[12],
                                                    "Cervical c." = discrete_palette[3],
                                                    "Colon c." = discrete_palette[4],
                                                    "Lung c." = discrete_palette[7]))
  pheatmap_benchmark <- function(matrix,
                                 title,
                                 cluster_cols = TRUE,
                                 no_clusters = NA){
  pheatmap(matrix,
           annotation_row = pheatmap_annotation_benchmark,
           annotation_colors = annotation_colors_benchmark,
           color = viridis(10),
           main = title,
           cluster_cols = cluster_cols,
           cutree_rows = no_clusters)
  }
  benchmark_cos_sim_pheatmap <- 
    pheatmap_benchmark(cos_sim_benchmark_signatures_exome,
                       "Heatmap of pairwise cosine similarity between original benchmark samples profiles and exomes COSMIC signatures")  
  save_plot("benchmark/COSMIC/cosine_similarity_original_vs_cosmic/pheatmap_benchmark_cos_sim.png", 
            benchmark_cos_sim_pheatmap)
  benchmark_cos_sim_pheatmap_clustered <- 
    pheatmap_benchmark(cos_sim_benchmark_signatures_exome,
                       "Heatmap of pairwise cosine similarity between original benchmark samples profiles and exomes COSMIC signatures (3 clusters)",
                       no_clusters = 3)  
  save_plot("benchmark/COSMIC/cosine_similarity_original_vs_cosmic/benchmark_cos_sim_pheatmap_clustered.png", 
            benchmark_cos_sim_pheatmap_clustered)
}  


###### 6. Reconstruct 96 mutational profiles
writeLines("\nRecostruction of 96 mutational profiles")

##### 6.1 Find optimal contribution of COSMIC signatures to reconstruct 96 mutational profiles

#### 6.1.2 Use all COSMIC signatures
print("Reconstruction using all COSMIC signatures")

# Fit mutation matrix to the COSMIC signatures:
fit_res_exomes <- fit_to_signatures(mut_mat, cancer_signatures_exomes)


### Contribution after fitting

## Make a boxplot for the signature contribution
print("Generating boxplots for signatures contribution")
exomes_signatures_short <- gsub("S.", "", colnames(cancer_signatures_exomes))

fit_res_contribution_long <- 
  rownames_to_column(as.data.frame(fit_res_exomes$contribution), 
                     var = "signature") %>%
  gather(key = sample, value = contribution,
         starts_with("F")) %>%
  mutate(signature = factor(exomes_signatures_short,
                            levels = exomes_signatures_short))

vector_tumortypes <- c()
for (element in metadata_sorted_by_samplesID$Subtype){
  vector_tumortypes <- c(vector_tumortypes, rep(element, 65))
}

fit_res_contribution_long <- 
  cbind(fit_res_contribution_long, vector_tumortypes)

# Just signatures contribution to samples
fit_res_contribution_long %>%
  mutate(signature = factor(signature,
                            levels = exomes_signatures_short)) %>%
  ggplot(aes(x = signature, y = contribution)) + 
  geom_boxplot() +
  ylab(label = "Contribution") +
  xlab(label = "Signatures") + 
  ylim(0, 500) +
  labs(title = "Boxplot of exomes signatures contribution") +
  theme_bw() + 
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=11)) -> boxpot_signatures_contribution

save_plot("COSMIC/sig_contributions/boxplots/boxpot_exomes_signatures_contribution.png", 
          boxpot_signatures_contribution)

# Signatures contribution grouped by tumor types
fit_res_contribution_long %>% 
  ggplot(aes(x = signature, y = contribution)) + 
  geom_boxplot() + 
  labs(title = "Boxplot of exomes signatures contribution grouped by tumor type") +
  ylab(label = "Contribution") +
  xlab(label = "Signatures") + 
  ylim(0, 2000) +
  facet_wrap(~ vector_tumortypes, 
             nrow = length(unique(vector_tumortypes))) + 
  theme_bw() +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=11)) -> boxplot_signatures_contribution_by_tumortype

save_plot("COSMIC/sig_contributions/boxplots/boxpot_exomes_signatures_contribution_by_tumortype.png", 
          boxplot_signatures_contribution_by_tumortype)

# Signature contribution for each tumor type in separate plots
boxplot_contribution_by_tumor <- function(tumor_type, 
                                          fit_res_contribution_long){
  fit_res_contribution_long %>% 
    filter(vector_tumortypes == tumor_type) %>%
    ggplot(aes(x = signature, y = contribution)) + 
    geom_boxplot() + 
    labs(title = paste("Boxplot of exomes signatures contribution in", 
                       gsub("_", " ", tolower(tumor_type)))) +
    ylab(label = "Contribution") +
    xlab(label = "Signatures") + 
    theme_bw() +
    theme(axis.text=element_text(size=9),
          axis.title=element_text(size=11))
}

boxplots_contribution_by_tumor_separated <- function(print_boxplot = FALSE,
                                                     save_boxplot = TRUE){
  print("Generating a boxplot for each tumor tissue:")
  for (tumor in unique(metadata_sorted_by_samplesID$Subtype)){
    print(paste(gsub("_", " ", tumor)))
    boxplot <- boxplot_contribution_by_tumor(tumor, fit_res_contribution_long)
    if (save_boxplot == TRUE){
      save_plot(paste("COSMIC/sig_contributions/boxplots/boxplots_each_tumor/boxpot_exomes_signatures_contribution_", 
                      gsub(" ", "_", tumor), ".png", sep = ""), 
                boxplot)
    }
    if (print_boxplot == TRUE){
      print(boxplot)
      readline(prompt="Press [Enter] to continue")
    }
  }
}

boxplots_contribution_by_tumor_separated(print_boxplot = FALSE)

## Look at the total contribution across the samples 

# Barplot of the total contribution zooming to low contribution signatures 
print("Generating barplots for signatures contribution")

fit_res_total_contribution <- 
  rownames_to_column(as.data.frame(fit_res_exomes$contribution), 
                     var = "signature") %>%
  as_tibble() %>%
  mutate(sum_contribution = apply(fit_res_exomes$contribution, 1, sum)) %>%
  mutate(signature = factor(gsub("S.", "", signature),
                            levels = exomes_signatures_short)) %>%
  select(signature, sum_contribution)

lowest_contribution_sig <- fit_res_total_contribution %>%
  filter(sum_contribution != 0) %>%
  filter(sum_contribution == min(sum_contribution))

fit_res_total_contribution %>%
  ggplot(aes(x = as.numeric(signature), y = sum_contribution)) +
  geom_bar(stat = "identity") +
  xlab(label = "Signatures") + 
  ylab(label = "Sum of mutational contribution") +
  labs(title = "Sum of mutational contribution by each signature, zoomed to the signature with lowest contribution (different than zero)") +
  facet_zoom(xy = signature == lowest_contribution_sig$signature[1], zoom.size = 0.02) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20)) +
  scale_x_continuous(label = exomes_signatures_short,
                     breaks = 1:length(exomes_signatures_short)) +
  theme_bw() +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=11)) -> total_contribution_barplot 

save_plot("COSMIC/sig_contributions/barplot_exomes_signatures_total_contribution.png", 
          total_contribution_barplot)

# Select signatures with some contribution
sum_contribution_threshold <- 10
selected_exomes <- which(rowSums(fit_res_exomes$contribution) > sum_contribution_threshold)

## Plot signatures contribution as barplot
barplot_contribution_no_filter_exomes <- 
  plot_contribution(fit_res_exomes$contribution[selected_exomes,], 
                    cancer_signatures_exomes[,selected_exomes], 
                    coord_flip = TRUE, 
                    mode = "absolute")
save_plot("COSMIC/sig_contributions/barplots/barplot_contribution_no_filter_exomes.png", 
          barplot_contribution_no_filter_exomes)

## Mutational pattern heatmap
heatmap_contribution_no_filter_exome <- 
  plot_contribution_heatmap(fit_res_exomes$contribution, 
                            cluster_samples = TRUE, 
                            method = "complete")
save_plot("COSMIC/sig_contributions/heatmaps/mutationalpattern_heatmap_contribution_no_filter_exome.png",
          heatmap_contribution_no_filter_exome)


## Obtain relative contribution
normalizer <- function(row_or_col){
  n <- row_or_col/(sum(row_or_col)) 
  return(n)
}

fit_res_exomes$contribution %>% t -> fit_res_exomes_contribution_t
fit_res_exomes_relative_contribution_t <- apply(fit_res_exomes_contribution_t, 1, normalizer)
fit_res_exomes_relative_contribution_t %>% t -> fit_res_exomes_relative_contribution_t

if (use_case_samples_only == FALSE){
  fit_res_exomes_relative_contribution_t_benchmark <- 
    fit_res_exomes_relative_contribution_t[benchmark_sample_index,]
}

## Pheatmaps
print("Generating heatmaps for signatures contribution")

# Heatmap using metadata as annotation 
pheatmap_contribution_no_filter_exomes_metadata <- 
  generate_pheatmap(fit_res_exomes_relative_contribution_t, 
                    title = "Heatmap of all exomes COSMIC signatures relative contribution to recostructed profiles", 
                    metadata_annotation = TRUE,
                    fontsize = fontsize)
save_plot("COSMIC/sig_contributions/heatmaps/pheatmap_contribution_no_filter_exomes_metadata.png",
          pheatmap_contribution_no_filter_exomes_metadata)

# Heatmaps for all samples
if (use_case_samples_only == FALSE){
  # Treatment group as annotation 
  pheatmap_contribution_no_filter_exomes <- 
    generate_pheatmap(fit_res_exomes_relative_contribution_t, 
                      title = "Heatmap of all exomes COSMIC signatures relative contribution to recostructed profiles (annotated by tumor type and treatment group)",
                      fontsize = fontsize)
  save_plot("COSMIC/sig_contributions/heatmaps/pheatmap_contribution_no_filter_exome.png",
            pheatmap_contribution_no_filter_exomes)
} else{
# Heatmaps for case samples
  # Tumor types as annotation
  pheatmap_contribution_no_filter_exomes_ann_tumor <- 
    generate_pheatmap(fit_res_exomes_relative_contribution_t,
                      "Heatmap of all exomes COSMIC signatures relative contribution to recostructed profiles (annotated by tumor type)")
  save_plot("COSMIC/sig_contributions/heatmaps/pheatmap_contribution_no_filter_exomes_ann_tumor.png", 
            pheatmap_contribution_no_filter_exomes_ann_tumor)
  # Tumor types and treatment response as annotation
  pheatmap_contribution_no_filter_exomes_ann_tumor_response <- 
    generate_pheatmap(fit_res_exomes_relative_contribution_t,
                      "Heatmap of all exomes COSMIC signatures relative contribution to recostructed profiles (annotated by tumor type and treatment response)",
                      show_response = TRUE)
  save_plot("COSMIC/sig_contributions/heatmaps/pheatmap_contribution_no_filter_exomes_ann_tumor_response.png", 
            pheatmap_contribution_no_filter_exomes_ann_tumor_response)
}

# Generate and save n heatmaps with increasing number of clusters
if (use_case_samples_only == FALSE){
pheatmaps_iterative_clustered(fit_res_exomes_relative_contribution_t,
                              real_time_visualization = FALSE,
                              measure = "contribution")
} else{
  pheatmaps_iterative_clustered(fit_res_exomes_relative_contribution_t,
                                real_time_visualization = FALSE,
                                metadata_annotation = FALSE,
                                measure = "contribution",
                                fontsize = 7)
}

# Heatmap benchmark data
if (use_case_samples_only == FALSE){
  benchmark_contribution_heatmap <- 
    pheatmap_benchmark(fit_res_exomes_relative_contribution_t_benchmark,
                       "Heatmap of exomes COSMIC signatures relative contribution to recostructed profiles of benchmark samples")  
  save_plot("benchmark/COSMIC/sig_contributions/heatmap_benchmark_contribution.png", 
            benchmark_contribution_heatmap)
}


## PCA
print("Generating PCA plots for signatures contribution")

pca_plot_contribution_exomes_sig_input_no_filter <- 
  pca_with_metadata_signatures(fit_res_exomes_relative_contribution_t,
                               added_col = "variance",
                               "PCA plot of COSMIC signatures relative contribution after fitting (signatures as input, colored by variance)")

pca_plot_contribution_exomes_samples_input_col_subtype_no_filter <- 
  pca_with_metadata_samples(fit_res_exomes_relative_contribution_t,
                            "PCA plot of COSMIC signatures relative contribution after fitting (samples as input, colored by tumor type)")

pca_plot_contribution_exomes_samples_input_col_age_labeled_no_filter <-
  pca_with_metadata_samples(fit_res_exomes_relative_contribution_t,
                            "PCA plot of COSMIC signatures relative contribution after fitting (samples as input, colored by patient age, labeled by tumor type)",
                            colored_by = "age",
                            labels = metadata_sorted_by_samplesID$Subtype)

if (use_case_samples_only == FALSE){
  pca_plot_contribution_exomes_samples_input_col_TMB_labeled_no_filter <-
    pca_with_metadata_samples(fit_res_exomes_relative_contribution_t,
                              "PCA plot of COSMIC signatures relative contribution after fitting (samples as input, colored by TMB, labeled by tumor type)",
                              colored_by = "TMB",
                              labels = metadata_sorted_by_samplesID$Subtype)
  
  pca_plot_contribution_exomes_samples_input_col_date_labeled_no_filter <-
    pca_with_metadata_samples(fit_res_exomes_relative_contribution_t,
                              "PCA plot of COSMIC signatures relative contribution after fitting (samples as input, colored by sequencing date, labeled by tumor type)",
                              colored_by = "seq_date",
                              labels = metadata_sorted_by_samplesID$Subtype)
  pca_plot_contribution_exomes_samples_input_col_group_labeled_subtype <- 
    pca_with_metadata_samples(fit_res_exomes_relative_contribution_t,
                              "PCA plot of COSMIC signatures relative contribution after fitting (samples as input, colored by group, labeled by tumor type)",
                              colored_by = "group",
                              labels = metadata_sorted_by_samplesID$Subtype)
} else{
  pca_plot_contribution_exomes_samples_input_col_subtype_labeled_cancergroup <- 
    pca_with_metadata_samples(fit_res_exomes_relative_contribution_t,
                              "PCA plot of COSMIC signatures relative contribution after fitting (samples as input, colored and labeled by cancer group)",
                              labels = metadata_sorted_by_samplesID$Cancer_group)
  pca_plot_contribution_exomes_samples_input_col_subtype_labeled_HRD <- 
    pca_with_metadata_samples(fit_res_exomes_relative_contribution_t,
                              "PCA plot of COSMIC signatures relative contribution after fitting (samples as input, colored by cancer group, labeled by Somatic HRD mutation)",
                              labels = metadata_sorted_by_samplesID$Somatic_HRD_mt)
  pca_plot_contribution_exomes_samples_input_col_response_labeled_cancergroup <- 
    pca_with_metadata_samples(fit_res_exomes_relative_contribution_t,
                              "PCA plot of COSMIC signatures relative contribution after fitting (samples as input, colored by response, labeled by cancer group)",
                              labels = metadata_sorted_by_samplesID$Subtype,
                              colored_by = "response")
  
}
# Save plots
save_plot("COSMIC/sig_contributions/pca/pca_plot_contribution_exomes_sig_input_no_filter.png",
          pca_plot_contribution_exomes_sig_input_no_filter)
save_plot("COSMIC/sig_contributions/pca/pca_plot_contribution_exomes_samples_input_col_subtype_no_filter.png",
          pca_plot_contribution_exomes_samples_input_col_subtype_no_filter)
save_plot("COSMIC/sig_contributions/pca/pca_plot_contribution_exomes_samples_input_col_age_labeled_no_filter.png",
          pca_plot_contribution_exomes_samples_input_col_age_labeled_no_filter)
if (use_case_samples_only == FALSE){
  save_plot("COSMIC/sig_contributions/pca/pca_plot_contribution_exomes_samples_input_col_TMB_labeled_no_filter.png",
            pca_plot_contribution_exomes_samples_input_col_TMB_labeled_no_filter)
  save_plot("COSMIC/sig_contributions/pca/pca_plot_contribution_exomes_samples_input_col_date_labeled_no_filter.png",
            pca_plot_contribution_exomes_samples_input_col_date_labeled_no_filter)
  save_plot("COSMIC/sig_contributions/pca/pca_plot_contribution_exomes_samples_input_col_group_labeled_subtype.png",
            pca_plot_contribution_exomes_samples_input_col_group_labeled_subtype)
} else{
  save_plot("COSMIC/sig_contributions/pca/pca_plot_contribution_exomes_samples_input_col_subtype_labeled_cancergroup.png",
            pca_plot_contribution_exomes_samples_input_col_subtype_labeled_cancergroup)
  save_plot("COSMIC/sig_contributions/pca/pca_plot_contribution_exomes_samples_input_col_subtype_labeled_HRD.png",
            pca_plot_contribution_exomes_samples_input_col_subtype_labeled_HRD)
  save_plot("COSMIC/sig_contributions/pca/pca_plot_contribution_exomes_samples_input_col_response_labeled_cancergroup.png",
            pca_plot_contribution_exomes_samples_input_col_response_labeled_cancergroup)
}

# Benchmark data
if (use_case_samples_only == FALSE){
  benchmark_pca_contribution <- 
    pca_benchmark(fit_res_exomes_relative_contribution_t_benchmark, 
                  title = "PCA plot of COSMIC signatures relative contribution to benchmark samples mutations (labeled by sample ID)",
                  label = "sample")
  save_plot("benchmark/COSMIC/sig_contributions/pca_benchmark_contribution.png", 
            benchmark_pca_contribution)
  
  benchmark_pca_contribution_annotation <- 
    pca_benchmark(fit_res_exomes_relative_contribution_t_benchmark, 
                  title = "PCA plot of COSMIC signatures relative contribution to benchmark samples mutations (labeled by annotation)",
                  label = "annotation")
  save_plot("benchmark/COSMIC/sig_contributions/pca_benchmark_contribution.png", 
            benchmark_pca_contribution_annotation)
}


## PCA biplots (all samples only)
if (use_case_samples_only == FALSE){
  biplot_contribution_exomes_no_filter <- plot_pca_biplot(fit_res_exomes_relative_contribution_t, 
                                                          "PCA biplot of samples by exomes COSMIC signatures relative contribution after fitting",
                                                          arrow = FALSE,
                                                          point_size = 2)
  save_plot("COSMIC/sig_contributions/pca_biplots/biplot_contribution_exomes_no_filter.png",
            biplot_contribution_exomes_no_filter)
}


#### 6.1.3 Fit only to signatures with cosine similarity larger than threshold
cos_sim_threshold = 0.2
paste("Reconstructing profile using signatures filtered by cosine similarity, threshold", 
      cos_sim_threshold)

# Filtering function by cosine similarity
get_filtered_cos_similarity <- function(cos_sim_samples_signatures_matrix, 
                                        threshold){
  # Save sample names
  sample_names <- rownames(cos_sim_samples_signatures_matrix)
  # Convert the matrix in a long format tibble
  rownames_to_column(as.data.frame(cos_sim_samples_signatures_matrix),  var="sample") %>% 
    as_tibble() %>% 
    gather(key = signature, value = cos_sim, starts_with("S.")) -> cos_sim_tibble_long
  # Function to filter by cosine similarity
  get_filter_cos_sim <- function(sample_name){
    cos_sim_tibble_long %>%
      filter(sample==sample_name, cos_sim > threshold) -> filtered_cos_sim_tibble
    return(filtered_cos_sim_tibble)    
  }
  # Get a list whose elements are the a filtered tibble for each sample
  list_cos_sim_tibble <- lapply(sample_names, get_filter_cos_sim)
  # Assign samples names to the names of the list
  names(list_cos_sim_tibble) = sample_names
  return(list_cos_sim_tibble)
}

# Function to check how many sample do not have n signatures
return_no_sample_with_less_n_signatures <- function(list, n){
  no_sample = 0
  for (element in list){
    if (dim(element)[1] < n){
      no_sample = no_sample + 1
    }
  }
  return(no_sample)
}

# Function to get filtered signatures
fit_to_filtered_signatures <- function(mutation_matrix,
                                       cancer_signatures,
                                       signatures_filtered,
                                       print_no_signatures = FALSE){
  n = 1
  for (element in signatures_filtered){
    # Select sample
    sample <- (names(signatures_filtered)[n])
    # Select signatures with filtered cos similarity 
    filtered_signatures_names <- element$signature
    if (print_no_signatures == TRUE){
      print(paste("Filtered signatures for", sample, ":", length(filtered_signatures_names)))
    }
    filtered_cancer_signatures <- cancer_signatures[,filtered_signatures_names]
    # Select mutation profile of sample
    selected_sample_mut_mat <- as.data.frame(mutation_matrix[,sample])
    # Reconstruct the mutational profile of the sample by fitting the signatures filtered by cosine similarity
    fit_res_one_sample <- fit_to_signatures(selected_sample_mut_mat, filtered_cancer_signatures)
    # Assign column name
    colnames(fit_res_one_sample$contribution) = sample
    colnames(fit_res_one_sample$reconstructed) = sample
    # Save the reconstructed profile
    if (n == 1){
      fit_res <- fit_res_one_sample$reconstructed
    } else{
      fit_res <- cbind(fit_res, fit_res_one_sample$reconstructed)
    }
    n = n + 1
  }
  return(fit_res)
}

# Function to stop a function without raising error
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

# Function to asses the quality of the reconstructed profiles
get_avg_cos_sim_ori_rec <- function(mutation_matrix, 
                                    reconstructed_mut_mat){
  cos_sim_ori_rec <- cos_sim_matrix(mutation_matrix, reconstructed_mut_mat)
  cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
  avg_cos_sim_ori_rec <- mean(cos_sim_ori_rec[,1])
  return(avg_cos_sim_ori_rec)
}

# Function to fit each sample to the filtered signatures by cosine similarity
similarity_filter_and_fit <- function(mutation_matrix, 
                                      cancer_signatures, 
                                      cos_sim_threshold,
                                      min_signatures_per_sample = 2,
                                      print_no_signatures = FALSE,
                                      assess_quality = FALSE){
  # Compute pairwise cosine similarity
  cos_sim_samples_signatures <- cos_sim_matrix(mutation_matrix, cancer_signatures)
  # Filter signatures by cosine similarity to original profile for each sample
  signatures_filtered_by_cos_sim <- get_filtered_cos_similarity(cos_sim_samples_signatures, 
                                                                cos_sim_threshold)
  # Stop if the minimum number of signatures per sample is not reached
  no_sample_with_few_sig <- 
    return_no_sample_with_less_n_signatures(signatures_filtered_by_cos_sim, 
                                            min_signatures_per_sample)
  if (no_sample_with_few_sig > 0){
    print(paste("Error: there are", no_sample_with_few_sig, 
                "samples with less than", min_signatures_per_sample ,"filtered signatures"))
    print("Try a smaller cosine similarity threshold")
    stop_quietly()
  }  
  # Fit to filtered signatures
  reconsturcted_profiles <- fit_to_filtered_signatures(mutation_matrix,
                                                       cancer_signatures,
                                                       signatures_filtered_by_cos_sim,
                                                       print_no_signatures)
  # Assess the quality of the reconstructed profiles
  if (assess_quality == TRUE){
    print(get_avg_cos_sim_ori_rec(mutation_matrix, reconsturcted_profiles))
  }
  return(reconsturcted_profiles)
}

# Filter and fit the mutation matrix to the filtered signatures by cosine similarity
# (Try assess_quality = TRUE, or print_no_signatures = TRUE)
fit_res_sim_filtered_exomes <-  similarity_filter_and_fit(mut_mat,
                                                          cancer_signatures_exomes,
                                                          cos_sim_threshold = cos_sim_threshold,
                                                          print_no_signatures = FALSE,
                                                          assess_quality = FALSE)  


#### 6.1.4 Fit only to signatures that had a contribution larger than threshold after fitting to all signatures
contribution_threshold = 0.02
paste("Reconstructing profile using signatures filtered by contribution, threshold", contribution_threshold)

# Normalize the signatures contribution and convert to tibble
get_relative_contribution_tibble <- function(fit_res, selected_signatures){
  fit_relative_contribution <- apply(fit_res$contribution[selected_signatures,], 2, normalizer)
  rownames_to_column(as.data.frame(fit_relative_contribution),  var="signature") %>% 
    as_tibble() -> fit_relative_contribution_tibble
  return(fit_relative_contribution_tibble)
}

# Filtering function by relative contribution
get_filtered_contribution <- function(fit_relative_contribution_tibble, threshold){
  # Extract sample names from the tibble
  sample_names <- colnames(fit_relative_contribution_tibble[,2:ncol(fit_relative_contribution_tibble)])
  # Function to filter by relative contribution
  get_sig_larger_than_threshold <- function(sample_name){
    fit_relative_contribution_tibble %>%
      filter(get(sample_name) > threshold) %>%
      select(signature, all_of(sample_name)) -> filtered_contribution_tibble
    return(filtered_contribution_tibble)
  }
  # Get a list where elements are filtered tibbles for each sample
  list_filtered_by_contribution <- lapply(sample_names, get_sig_larger_than_threshold)
  # Assign samples names to the names of the list
  names(list_filtered_by_contribution) = sample_names
  return(list_filtered_by_contribution)
}

# Function to fit each sample to the filtered signatures by relative contribution
contribution_filter_and_fit <- function(mutation_matrix, 
                                        cancer_signatures, 
                                        contribution_threshold,
                                        min_signatures_per_sample = 2,
                                        print_no_signatures = FALSE,
                                        assess_quality = FALSE){
  # Fit to all signatures
  fit_res <- fit_to_signatures(mutation_matrix, cancer_signatures)
  # Select signatures with some contribution
  selected <- which(rowSums(fit_res$contribution) > 10)
  # Get signatures relative contribution
  contribution_samples_signatures <- get_relative_contribution_tibble(fit_res, selected)
  # Get signatures filtered by relative contribution
  signatures_filtered_by_contribution <- 
    get_filtered_contribution(contribution_samples_signatures, contribution_threshold)
  # Stop if the minimum number of signatures per sample is not reached
  n = min_signatures_per_sample
  no_sample_with_few_sig <- 
    return_no_sample_with_less_n_signatures(signatures_filtered_by_contribution, n)
  if (no_sample_with_few_sig > 0){
    print(paste("Error: there are", no_sample_with_few_sig, 
                "samples with less than", n ,"filtered signatures"))
    print("Try a smaller contribution threshold")
    stop_quietly()
  }
  # Fit to filtered signatures
  reconsturcted_profiles <- fit_to_filtered_signatures(mutation_matrix,
                                                       cancer_signatures,
                                                       signatures_filtered_by_contribution,
                                                       print_no_signatures)
  # Assess the quality of the reconstructed profiles
  if (assess_quality == TRUE){
    print(paste("QC cosine similarity original VS reconstructed profiles:", 
                round(get_avg_cos_sim_ori_rec(mutation_matrix, reconsturcted_profiles), 7)))
  }
  return(reconsturcted_profiles)
}

# Filter and fit the mutation matrix to the filtered signatures by contribution
# (Try assess_quality = TRUE, or print_no_signatures = TRUE)
fit_res_contr_filtered_exomes <-  contribution_filter_and_fit(mut_mat,
                                                              cancer_signatures_exomes,
                                                              contribution_threshold = contribution_threshold,
                                                              print_no_signatures = FALSE,
                                                              assess_quality = FALSE)    


#### 7. De novo extraction
writeLines("\nDe novo extraction")

# Estimate the optimal rank
mut_mat <- mut_mat + 0.0001
estimate <- nmf(mut_mat, rank=2:5, method="brunet", nrun=10, seed=123456)
estimate_plot <- plot(estimate) +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, size = 25))
save_plot("de_novo/rank_estimation.png",
          estimate_plot)
save_plot("relevant_results/rank_estimation.png",
          estimate_plot,
          height = 7, width= 10)

no_rank = 4
denovo_signatures <- c("S.A", "S.B", "S.C", "S.D")
rank_test <- c(3,5,6)

# Extract signatures from the mutationa count matrix
nmf_res <- extract_signatures(mut_mat, rank = no_rank, nrun = 100)
colnames(nmf_res$signatures) <- denovo_signatures
rownames(nmf_res$contribution) <- denovo_signatures

# Plot the 96-profile of the signatures:
denovo_96_profile <- plot_96_profile(nmf_res$signatures, condensed = TRUE) +
  ggtitle("Mutation spectrum de novo extracted signatures") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))
save_plot("de_novo/de_novo_96_profile_plots.png",
          denovo_96_profile)
save_plot("relevant_results/de_novo_96_profiles.png",
          denovo_96_profile,
          width = 10, height = 5)

# Visualize the relative and absolute contribution of the signatures in a barplot
pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signature, 
                         mode = "relative") +
  theme(axis.text.x = element_text(angle = 90))
pc2 <- plot_contribution(nmf_res$contribution, nmf_res$signature, 
                         mode = "absolute") +
  theme(axis.text.x = element_text(angle = 90))

de_novo_barplot_contribution <- grid.arrange(pc1, pc2)

save_plot("de_novo/sig_contributions/barplots/de_novo_barplots_contributions.png",
          de_novo_barplot_contribution)

## Cosine similarity between COSMIC and novo extracted signatures
cos_sim_denovo_cosmic <- cos_sim_matrix(nmf_res$signatures, cancer_signatures_exomes)

# Heatmap
pheatmap_cos_sim_denovo_cosmic <- pheatmap(t(cos_sim_denovo_cosmic),
                                           cluster_cols = FALSE,
                                           color = viridis(10),
                                           main = "Cosine similarity between de novo \nextracted and exomes COSMIC signatures") 
  
save_plot("de_novo/pheatmap_cos_sim_denovo_cosmic.png",
          pheatmap_cos_sim_denovo_cosmic)
save_plot("relevant_results/pheatmap_cos_sim_denovo_cosmic.png",
          pheatmap_cos_sim_denovo_cosmic,
          height = 10, width=4)

# Obtain COSMIC signatures most similar to the de novo extracted ones
extract_most_similar_signature <- function(no_signature = 5,
                                           signatures){
  n = 1
  for (s in signatures){
    rownames_to_column(as.data.frame(t(cos_sim_denovo_cosmic)), var = "Signature") %>% 
      as_tibble() %>% 
      arrange(desc(get(s))) %>% 
      head(no_signature) %>% 
      select(Signature, all_of(s)) -> df
    if (n == 1){
      final_df <- df
    } else{
      final_df <- cbind(final_df, df)
    }
    n = n + 1
  }
  colnames(final_df) <- gsub("S\\.", "Sim_to_S_", colnames(final_df))
  return(final_df)
}

similar_COSMIC_signatures <- extract_most_similar_signature(no_signature = 10,
                                                            denovo_signatures)
write_tsv(similar_COSMIC_signatures, 
          paste(results_path, "de_novo/top10_similar_COSMIC_signatures.txt", sep = ""))
write_tsv(similar_COSMIC_signatures, 
          paste(results_path, "relevant_results/top10_similar_COSMIC_signatures.txt", sep = ""))
similar_COSMIC_signatures

### 7.1 Reconstructed profile from nmf

## Heatmaps de novo relative contribution
nmf_res$contribution %>% t -> nmf_res_contribution_t
nmf_res_relative_contribution_t <- apply(nmf_res_contribution_t, 1, normalizer)
nmf_res_relative_contribution_t %>% t -> nmf_res_relative_contribution_t

## Heatmap with MutationalPattern
mp_heatmap_contribution_denovo <- 
  plot_contribution_heatmap(nmf_res$contribution, 
                            cluster_samples = TRUE,
                            method = "complete")
save_plot("de_novo/sig_contributions/heatmaps/denovo_mpheatmap_contribution.png",
          mp_heatmap_contribution_denovo)

## Pheatmaps

# Heatmap using metadata as annotation 
pheatmap_contribution_denovo_metadata <- 
  generate_pheatmap(nmf_res_relative_contribution_t, 
                    title = "De novo extracted signatures relative \ncontribution to recostructed profiles", 
                    metadata_annotation = TRUE,
                    cluster_cols = FALSE,
                    no_clusters = 8)
save_plot("de_novo/sig_contributions/heatmaps/denovo_pheatmap_contribution_metadata.png",
          pheatmap_contribution_denovo_metadata)
save_plot("relevant_results/denovo_heatmap_contribution_replicate.png",
          pheatmap_contribution_denovo_metadata,
          width = 5, height = 12)
# Plots for all samples
if (use_case_samples_only == FALSE){
  # Heatmap using treatment group as annotation 
  pheatmap_contribution_denovo <- 
    generate_pheatmap(nmf_res_relative_contribution_t, 
                      title = "Heatmap of de novo extracted signatures relative contribution to recostructed profiles (annotated by tumor type and treatment group)",
                      cluster_cols = FALSE)
  save_plot("de_novo/sig_contributions/heatmaps/denovo_pheatmap_contribution.png",
            pheatmap_contribution_denovo)
} else{
# Plots case samples
  # Heatmap with only tumor types as annotation
  pheatmap_contribution_denovo_ann_tumor <- 
    generate_pheatmap(nmf_res_relative_contribution_t,
                      "Heatmap of de novo extracted signatures relative contribution to recostructed profiles (annotated by tumor type)",
                      cluster_cols = FALSE)
  save_plot("de_novo/sig_contributions/heatmaps/denovo_pheatmap_contribution_ann_tumor.png", 
            pheatmap_contribution_denovo_ann_tumor)
  # Heatmap annotated by tumor type and treatment response
  pheatmap_contribution_denovo_ann_tumor_response <- 
    generate_pheatmap(nmf_res_relative_contribution_t,
                      "Heatmap of de novo extracted signatures relative contribution to recostructed profiles (annotated by tumor type and treatment response)",
                      show_response = TRUE,
                      cluster_cols = FALSE)
  save_plot("de_novo/sig_contributions/heatmaps/denovo_pheatmap_contribution_ann_tumor_response.png", 
            pheatmap_contribution_denovo_ann_tumor_response)
}

# Generate and save n heatmaps with increasing number of clusters 
if (use_case_samples_only == FALSE){
  pheatmaps_iterative_clustered(nmf_res_relative_contribution_t,
                                real_time_visualization = FALSE,
                                measure = "contribution",
                                de_novo_signatures = TRUE,
                                cluster_cols = FALSE)
} else{
  pheatmaps_iterative_clustered(nmf_res_relative_contribution_t,
                                real_time_visualization = FALSE,
                                metadata_annotation = TRUE,
                                measure = "contribution",
                                fontsize = 7,
                                de_novo_signatures = TRUE,
                                cluster_cols = FALSE)
}


## PCA de novo relative contribution

# Generate plots
print("Generating PCA plots for de novo extracted signatures contribution")
pca_plot_contribution_denovo_sig_input <- 
  pca_with_metadata_signatures(nmf_res_relative_contribution_t,
                               added_col = "variance",
                               "PCA plot de novo extracted signatures relative contribution (signatures as input, colored by variance)")
pca_plot_contribution_denovo_samples_input_col_subtype <- 
  pca_with_metadata_samples(nmf_res_relative_contribution_t,
                            "PCA plot de novo extracted signatures relative contribution (samples as input, colored by tumor type)")
pca_plot_contribution_denovo_samples_input_col_age_labeled <-
  pca_with_metadata_samples(nmf_res_relative_contribution_t,
                            "PCA plot de novo extracted signatures relative contribution (samples as input, colored by patient age, labeled by tumor type)",
                            colored_by = "age",
                            labels = metadata_sorted_by_samplesID$Subtype)
if (use_case_samples_only == FALSE){
  pca_plot_contribution_denovo_samples_input_col_TMB_labeled <-
    pca_with_metadata_samples(nmf_res_relative_contribution_t,
                              "PCA plot de novo extracted signatures relative contribution (samples as input, colored by TMB, labeled by tumor type)",
                              colored_by = "TMB",
                              labels = metadata_sorted_by_samplesID$Subtype)
  pca_plot_contribution_denovo_samples_input_col_date_labeled <-
    pca_with_metadata_samples(nmf_res_relative_contribution_t,
                              "PCA plot de novo extracted signatures relative contribution (samples as input, colored by sequencing date, labeled by tumor type)",
                              colored_by = "seq_date",
                              labels = metadata_sorted_by_samplesID$Subtype)
} else{
  pca_plot_contribution_denovo_samples_input_col_subtype_labeled_cancergroup <- 
    pca_with_metadata_samples(nmf_res_relative_contribution_t,
                              "PCA plot of de novo extracted signatures relative contribution (samples as input, colored by and labeled by cancer group)",
                              labels = metadata_sorted_by_samplesID$Cancer_group)
  pca_plot_contribution_denovo_samples_input_col_subtype_labeled_HRD <- 
    pca_with_metadata_samples(nmf_res_relative_contribution_t,
                              "PCA plot of de novo extracted signatures relative contribution (samples as input, colored by cancer group, labeled by Somatic HRD mutation)",
                              labels = metadata_sorted_by_samplesID$Somatic_HRD_mt)
  pca_plot_contribution_denovo_samples_input_col_response_labeled_cancergroup <- 
    pca_with_metadata_samples(nmf_res_relative_contribution_t,
                              "PCA plot of de novo extracted signatures relative contribution (samples as input, colored by response, labeled by cancer group)",
                              labels = metadata_sorted_by_samplesID$Subtype,
                              colored_by = "response")
}
# Save plots
save_plot(paste("de_novo/sig_contributions/pca/denovo_pca_plot_contribution_sig_input.png", sep = ""),
          pca_plot_contribution_denovo_sig_input)
save_plot(paste("de_novo/sig_contributions/pca/denovo_pca_plot_contribution_samples_input_col_subtype.png", sep = ""),
          pca_plot_contribution_denovo_samples_input_col_subtype)
save_plot(paste("de_novo/sig_contributions/pca/denovo_pca_plot_contribution_samples_input_col_age_labeled.png", sep = ""),
          pca_plot_contribution_denovo_samples_input_col_age_labeled)
if (use_case_samples_only == FALSE){
  save_plot(paste("de_novo/sig_contributions/pca/denovo_pca_plot_contribution_samples_input_col_TMB_labeled.png", sep = ""),
            pca_plot_contribution_denovo_samples_input_col_TMB_labeled)
  save_plot(paste("de_novo/sig_contributions/pca/denovo_pca_plot_contribution_samples_input_col_date_labeled.png", sep = ""),
            pca_plot_contribution_denovo_samples_input_col_date_labeled)
} else{
  save_plot(paste("de_novo/sig_contributions/pca/denovo_pca_plot_contribution_samples_input_col_subtype_labeled_cancergroup.png", sep = ""),
            pca_plot_contribution_denovo_samples_input_col_subtype_labeled_cancergroup)
  save_plot(paste("de_novo/sig_contributions/pca/denovo_pca_plot_contribution_samples_input_col_subtype_labeled_HRD.png", sep = ""),
            pca_plot_contribution_denovo_samples_input_col_subtype_labeled_HRD)
  save_plot(paste("de_novo/sig_contributions/pca/pca_plot_contribution_denovo_samples_input_col_response_labeled_cancergroup.png", sep = ""),
            pca_plot_contribution_denovo_samples_input_col_response_labeled_cancergroup)
}


## PCA biplots de novo contribution
denovo_pca_biplot_contribution <- plot_pca_biplot(nmf_res_relative_contribution_t, 
                                                  "PCA biplot of pairwise cosine similarity between original samples profiles and exomes COSMIC signatures",
                                                  varname_adjust = 2,
                                                  point_size = 2)
save_plot("de_novo/sig_contributions/pca_biplots/denovo_pca_biplot_contribution.png", 
          denovo_pca_biplot_contribution)


#### 7.2 Benchmark samples de novo signature extraction
print("De novo extraction benchmark samples")
if (use_case_samples_only == FALSE){
  # Estimate the optimal rank
  mut_mat_benchmark <- t(t(mut_mat)[benchmark_sample_index,])
  
  estimate <- nmf(mut_mat_benchmark, rank = 2:5, method = "brunet", nrun = 10, seed = 123456)
  estimate_plot <- plot(estimate)
  save_plot("benchmark/de_novo/benchmark_rank_estimation.png",
            estimate_plot)
  
  benchmark_no_rank = 2
  benchmark_denovo_signatures <- c("S.A", "S.B")
  
  # Extract signatures from the mutationa count matrix
  nmf_res_benchmark <- extract_signatures(mut_mat_benchmark, rank = benchmark_no_rank, nrun = 100)
  colnames(nmf_res_benchmark$signatures) <- benchmark_denovo_signatures
  rownames(nmf_res_benchmark$contribution) <- benchmark_denovo_signatures
  
  # Plot the 96-profile of the signatures:
  denovo_96_profile_benchmark <- plot_96_profile(nmf_res_benchmark$signatures, condensed = TRUE)
  save_plot("benchmark/de_novo/de_novo_96_profile_plot_benchmark.png",
            denovo_96_profile_benchmark)
  
  # Visualize the relative and absolute contribution of the signatures in a barplot
  pc1 <- plot_contribution(nmf_res_benchmark$contribution, nmf_res_benchmark$signature, 
                           mode = "relative") +
    theme(axis.text.x = element_text(angle = 90))
  pc2 <- plot_contribution(nmf_res_benchmark$contribution, nmf_res_benchmark$signature, 
                           mode = "absolute") +
    theme(axis.text.x = element_text(angle = 90))
  
  de_novo_barplot_contribution_benchmark <- grid.arrange(pc1, pc2)
  
  save_plot("benchmark/de_novo/sig_contributions/de_novo_barplots_contributions_benchmark.png",
            de_novo_barplot_contribution_benchmark)
  
  ### Cosine similarity between COSMIC and novo extracted signatures
  benchmark_cos_sim_denovo_cosmic <- cos_sim_matrix(nmf_res_benchmark$signatures, cancer_signatures_exomes)
  
  # MutationalPattern similarity
  mpheatmap_cos_sim_denovo_cosmic_benchmark <- 
    plot_cosine_heatmap(benchmark_cos_sim_denovo_cosmic, cluster_rows = TRUE)
  save_plot("benchmark/de_novo/mpheatmap_cos_sim_denovo_cosmic_benchmark.png",
            mpheatmap_cos_sim_denovo_cosmic_benchmark)
  
  
  # Pheatmap similarity
  pheatmap_cos_sim_denovo_cosmic_benchmark <- 
    pheatmap(benchmark_cos_sim_denovo_cosmic,
    cluster_rows = FALSE,
    color = viridis(10),
    main = "Cosine similarity between benchmark de novo extracted signatures and exomes COSMIC signatures") 
  save_plot("benchmark/de_novo/pheatmap_cos_sim_denovo_cosmic.png",
            pheatmap_cos_sim_denovo_cosmic_benchmark)
  
  # Obtain COSMIC signatures most similar to the de novo extracted ones
  benchmark_similar_COSMIC_signatures <- extract_most_similar_signature(10,
                                                                        benchmark_denovo_signatures)
  write_tsv(benchmark_similar_COSMIC_signatures, 
            paste(results_path, "benchmark/de_novo/top10_similar_COSMIC_signatures_benchmark.txt", sep = ""))
  benchmark_similar_COSMIC_signatures
  
  ### Reconstructed profile from nmf
  
  ## Heatmaps de novo relative contribution
  nmf_res_benchmark$contribution %>% t -> nmf_res_benchmark_contribution_t
  nmf_res_benchmark_relative_contribution_t <- apply(nmf_res_benchmark_contribution_t, 1, normalizer)
  nmf_res_benchmark_relative_contribution_t %>% t -> nmf_res_benchmark_relative_contribution_t
  
  # Heatmap with MutationalPattern
  mp_heatmap_contribution_denovo_benchmark <- 
    plot_contribution_heatmap(nmf_res_benchmark$contribution, 
                              cluster_samples = TRUE,
                              method = "complete")
  save_plot("benchmark/de_novo/sig_contributions/denovo_mpheatmap_contribution.png",
            mp_heatmap_contribution_denovo_benchmark)
  
  # Pheatmaps
  pheatmap_benchmark_cos_sim_pheatmap_benchmark_clustered <- 
    pheatmap_benchmark(nmf_res_benchmark_relative_contribution_t,
                       "Heatmap of de novo extracted signatures relative contribution to recostructed profiles (3 clusters)",
                       no_clusters = 3)  
  save_plot("benchmark/de_novo/sig_contributions/denovo_pheatmap_contribution_benchmark_clustered.png",
            pheatmap_benchmark_cos_sim_pheatmap_benchmark_clustered)
  pheatmap_benchmark_cos_sim_pheatmap_benchmark <- 
    pheatmap_benchmark(nmf_res_benchmark_relative_contribution_t,
                       "Heatmap of de novo extracted signatures relative contribution to recostructed profiles")  
  save_plot("benchmark/de_novo/sig_contributions/denovo_pheatmap_contribution_benchmark.png",
            pheatmap_benchmark_cos_sim_pheatmap_benchmark)
  
  # Scatter plot de novo contribution
  denovo_benchmark_contribution_scatterplot <- 
    rownames_to_column(as.data.frame(nmf_res_benchmark_relative_contribution_t), var = "SampleID") %>% 
     as_tibble() %>%
     ggplot(aes(x = S.A, y = S.B,
                fill = benchmark_samples$Subtype, 
                label = benchmark_samples$Annotation)) +
     geom_point() +    
     geom_label_repel(segment.alpha = 0.5) +
     xlab(label = "Signature A") +
     ylab(label = "Signature B") + 
     labs(title = "De novo extracted signatures relative contribution to recostructed profiles (benchmark samples)",
          fill = "Tumor type") +
     scale_fill_manual(values = c(discrete_palette[3], discrete_palette[4], 
                                  discrete_palette[7], discrete_palette[12])) +
     theme_bw() 
  
  save_plot("benchmark/de_novo/sig_contributions/denovo_benchmark_contribution_scatterplot.png", 
            denovo_benchmark_contribution_scatterplot)
}


##### 7. Quality control 
writeLines("\nFinal quality control by cosine similarity between original and reconstructed profiles")

## Generate other reconstructed profiles for testing

### Calculate all pairwise cosine similarities between original and reconstructed mutational profile
avg_sim <- function(mat1, mat2){
  cos_sim <- cos_sim_matrix(mat1, mat2)
  cos_sim <- as.data.frame(diag(cos_sim))
  avg_cos_sim <- round(mean(cos_sim[,1]), 7)
  return(list(cos_sim = cos_sim,
              avg_cos_sim = avg_cos_sim))
}
# Reconstructed profiles using all signatures
sim_ori_exomes <- avg_sim(mut_mat, fit_res_exomes$reconstructed)
# Reconstructed profiles using signatures filtered by cosine similarity
sim_ori_exomes_sim_filtered <- avg_sim(mut_mat, fit_res_sim_filtered_exomes)
# Reconstructed profiles using signatures filtered by relative contribution
sim_ori_exomes_contr_filtered <- avg_sim(mut_mat, fit_res_contr_filtered_exomes)
# De novo reconstructed profiles from nmf
sim_ori_denovo <- avg_sim(mut_mat, nmf_res$reconstructed)


# Save as txt file
qc_sim_ori_rec <-  data.frame("Reconstructed_profile" = c("All COSMIC exomes sig.",
                                                          paste("Cosine similarity >", 
                                                                cos_sim_threshold),
                                                          paste("Relative contribution >",
                                                                contribution_threshold),
                                                          paste("De novo extraction rank", no_rank)),
                              "Avg_sim_to_original" = c(sim_ori_exomes$avg_cos_sim,
                                                        sim_ori_exomes_sim_filtered$avg_cos_sim,
                                                        sim_ori_exomes_contr_filtered$avg_cos_sim,
                                                        sim_ori_denovo$avg_cos_sim))
qc_sim_ori_rec %>% 
  write_tsv(paste(results_path, "qc_similarity_to_original.txt", sep = ""))
qc_sim_ori_rec

### Make cosine similarity barplots
generate_cos_sim_barplot <- function(cos_sim_ori_rec,
                                     title = "",
                                     lowest_ylim = 0.9,
                                     cut_off_line = 0.95){
  # Adjust data frame for plotting with gpplot
  colnames(cos_sim_ori_rec) = "cos_sim"
  cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)
  # Make barplot
  ggplot(cos_sim_ori_rec, aes(y = cos_sim, x = sample)) + 
    geom_bar(stat = "identity", fill = "skyblue4") + 
    coord_flip(ylim = c(lowest_ylim, 1)) +
    ylab("Cosine similarity\n original VS reconstructed") +
    xlab("") +
    labs(title = title) +
    theme_bw() +
    geom_hline(aes(yintercept = cut_off_line))
}

qc_barplots <- function(lowest_ylim,
                        use_lower_ylim = FALSE){
  print("Generating barplots of cosine similarity between reconstructed and original profiles")
  # Reconstructed profiles using all signatures
  cos_sim_barplot_exomes <- 
    generate_cos_sim_barplot(sim_ori_exomes$cos_sim,
                             lowest_ylim = lowest_ylim,
                             "Barplot of pairwise cosine similarity between original and reconstructed profiles (unfiltered exomes COSMIC set)")
  # De novo reconstructed profiles from nmf
  cos_sim_barplot_sim_ori_denovo <- 
    generate_cos_sim_barplot(sim_ori_denovo$cos_sim,
                             lowest_ylim = lowest_ylim,
                             "Barplot of pairwise cosine similarity between original and reconstructed profiles\n(de novo extraction)")
  
  # Save to server
  if (use_lower_ylim == TRUE){
    qc_directory = "barplots_low_ylim/"
  } else{
    qc_directory = ""
  }
  save_plot(paste("qc_cosine_similarity_original_reconstructed/", qc_directory,"cos_sim_barplot_exomes.png", sep = ""), 
            cos_sim_barplot_exomes)
  save_plot(paste("qc_cosine_similarity_original_reconstructed/", qc_directory,"cos_sim_barplot_sim_ori_denovo.png", sep = ""), 
            cos_sim_barplot_sim_ori_denovo)
}

qc_barplots(0.9)

## If there is a sample with cosine similarity less than 0.9, generate new barplots with lower y limit
lowest_sim_value <- min(c(min(sim_ori_exomes$cos_sim[1]),
                          min(sim_ori_exomes_sim_filtered$cos_sim[1]),
                          min(sim_ori_exomes_contr_filtered$cos_sim)[1],
                          min(sim_ori_denovo$cos_sim[1])))

if (lowest_sim_value < 0.9){
  qc_barplots(lowest_sim_value,
              use_lower_ylim = TRUE)
}


##### 8. Fit to COSMIC signatures that are most similar to de novo extracted ones

## Compute similarity between original profile and COSMIC signature most similar (larger than 0.7) 
## to de novo extracted signatures
  
# Extract COSMIC signatures that have similarity larger than threshold (to de novo extracted signatures)
threshold = 0.7
most_similar_signatures <- c()
for (n in seq(1, ncol(similar_COSMIC_signatures), by = 2)){
  similar_COSMIC_signatures[,n:(n+1)] %>% 
    filter(get(colnames(similar_COSMIC_signatures[n+1])) > threshold) %>% 
    select(Signature) -> iterative_similar_signatures
  most_similar_signatures <- c(most_similar_signatures, iterative_similar_signatures$Signature)
}

# Pairwise cosine similarity between COSMIC signatures (most similar to de novo extracted ones) and original profiles
cos_sim_samples_most_similar_signatures_exomes <- 
  cos_sim_matrix(mut_mat, cancer_signatures_exomes[,most_similar_signatures])

# PCA plots
pca_plot_cos_sim_exomes_samples_input_col_response_labeled_cancergroup <- 
  pca_with_metadata_samples(cos_sim_samples_most_similar_signatures_exomes,
                            "PCA plot of pairwise cosine similarity of sample original profile VS exomes COSMIC signatures (samples as input, colored by response, labeled by tumor type)",
                            labels = metadata_sorted_by_samplesID$Subtype,
                            colored_by = "response")


### Reconstruct profiles with most similar COSMIC signatures

# Fit mutation matrix to the COSMIC signatures:
fit_res_most_similar_exomes <- fit_to_signatures(mut_mat, cancer_signatures_exomes[,most_similar_signatures])

# Relative contribution after fitting
fit_res_most_similar_exomes$contribution %>% t -> similar_exomes_contribution_t
similar_exomes_relative_contribution_t <- apply(similar_exomes_contribution_t, 1, normalizer)
similar_exomes_relative_contribution_t %>% t -> similar_exomes_relative_contribution_t

# PCA plots
pca_plot_contribution_denovo_samples_input_col_response_labeled_cancergroup <- 
  pca_with_metadata_samples(similar_exomes_relative_contribution_t,
                            "PCA plot of de novo extracted signatures relative contribution (samples as input, colored by response, labeled by cancer group)",
                            labels = metadata_sorted_by_samplesID$Subtype,
                            colored_by = "response")

## Assess average cosine similarity between reconstructed and original profile
sim_ori_most_similar_exomes <- avg_sim(mut_mat, fit_res_most_similar_exomes$reconstructed)



##### 9. Output files for treatment response prediction
if (use_case_samples_only == TRUE){
  name = "case_samples"
} else{
  name = "all_samples"
}

output_dataframe <- function(matrix,
                             filename = ""){
  
  colnames(matrix) <- gsub("\\.", "_", colnames(matrix))
  matrix <- rownames_to_column(as.data.frame(matrix), var = "SampleID")
  
  write.table(matrix, 
              paste(results_path, filename, "_", name,".txt", sep = ""), 
              row.names = FALSE, append = FALSE, sep = ",", dec = ".", col.names = TRUE)
}

metadata_sorted_by_samplesID %>%
  mutate(TMB = metadata_dates_TMB$TMB) %>%
  select(2:6, TMB) -> metadata_output_file

write.table(metadata_output_file, 
            paste(results_path, "output_metadata_", name,".txt", sep = ""), 
            row.names = FALSE, append = FALSE, sep = ",", dec = ".", col.names = TRUE)

# Relative contribution of de novo extracted signatures
output_dataframe(nmf_res_relative_contribution_t,
                 paste("output_denovo_relative_contribution_rank2", name, sep = ""))
# Cosine similarity of original profile between all COSMIC signatures
output_dataframe(cos_sim_samples_signatures_exomes,
                 paste("output_cos_sim_samples_signatures_exomes", name, sep = ""))
# Cosine similarity of original profile between most similar COSMIC signatutes (to de novo extracted ones)
output_dataframe(cos_sim_samples_most_similar_signatures_exomes,
                 paste("output_cos_sim_samples_most_similar_signatures_exomes", name, sep = ""))