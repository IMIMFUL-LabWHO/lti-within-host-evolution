# Longitudinal evolution of SARS-CoV-2 during prolonged infection reveals heterogeneous evolutionary dynamics across the viral genome
## Alen Suljič
### 11.6.2026

# libraries
library(here)
library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(broom.mixed)
library(emmeans)
library(patchwork)
library(scales)
library(performance)
library(compositions)
library(vegan)
library(ggrepel)
#library(BiocManager)
#BiocManager::install("Biostrings")
library(Biostrings)

# create directories
dir.create(here("figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("supp_tables"), showWarnings = FALSE, recursive = TRUE)

# data loading
m <- read_csv(here("data", "metadata.csv"))
v <- read_csv(here("variants", "metadata.csv")) 

# create master tibble
mv <- m %>% 
  left_join(v, by = join_by("sample")) %>% 
  filter(PASS == TRUE) %>% 
  mutate(across(c(sample, REF, ALT, GFF_FEATURE, REF_CODON, REF_AA, ALT_CODON, ALT_AA,
                  lineage, patient, group, source, vax, th_antineo, th_cst, th_hma, th_antivir), factor))

summary(mv)

# data prep
mv1 <- mv %>% 
  mutate(mut_level = if_else(ALT_FREQ >= 0.5, "dominant", "minor"),
         alt2 = if_else(grepl("-", ALT), "-", ALT)) %>% 
  unite("ref_POS", REF:POS, sep = "", remove = FALSE) %>% 
  unite("nt_mut", ref_POS, sep = "", alt2, remove = FALSE) %>% 
  select(-c(ref_POS, GFF_FEATURE)) %>% 
  mutate(mut_type = case_when(grepl("-", ALT) ~ "Deletion",
                              grepl("\\+", ALT) ~ "Insertion",
                              !is.na(REF_AA) ~ "SNP",
                              ALT_AA == "*" ~ "SNP"),
         mut_type = if_else(is.na(mut_type), "Non-coding", mut_type),
         snp_type = if_else(!REF_AA == ALT_AA, "Non-synonymous", "Synonymous"),
         mut = if_else(is.na(snp_type), mut_type, snp_type),
         mut2 = if_else(mut_type == "SNP", mut, "INDEL&NC"),
         mut3 = if_else(mut2 == "INDEL&NC", mut2, "SNP"),
         across(c(sample, nt_mut, REF, ALT, REF_CODON, REF_AA, ALT_CODON, ALT_AA, mut_level, alt2,
                  mut_type, snp_type, mut), factor)) %>%
  distinct(nt_mut, sample, .keep_all = TRUE) %>% 
  mutate(gene = case_when(POS <= 265 ~ "5'UTR",
                           POS >= 266 & POS <= 21555 ~ "ORF1ab",
                           POS >= 21556 & POS <= 21562 ~ "ITR1",
                           POS >= 21563 & POS <= 25384 ~ "S",
                           POS >= 25385 & POS <= 25392 ~ "ITR2",
                           POS >= 25393 & POS <= 26220 ~ "ORF3a",
                           POS >= 26221 & POS <= 26244 ~ "ITR3",
                           POS >= 26245 & POS <= 26472 ~ "E",
                           POS >= 26473 & POS <= 26522 ~ "ITR4",
                           POS >= 26523 & POS <= 27191 ~ "M",
                           POS >= 27192 & POS <= 27201 ~ "ITR5",
                           POS >= 27202 & POS <= 27387 ~ "ORF6a",
                           POS >= 27388 & POS <= 27393 ~ "ITR6",
                           POS >= 27394 & POS <= 27887 ~ "ORF7ab",
                           POS >= 27888 & POS <= 27893 ~ "ITR7",
                           POS >= 27894 & POS <= 28259 ~ "ORF8",
                           POS >= 28260 & POS <= 28273 ~ "ITR8",
                           POS >= 28274 & POS <= 29533 ~ "N",
                           POS >= 29534 & POS <= 29557 ~ "ITR9",
                           POS >= 29558 & POS <= 29674 ~ "ORF10",
                           POS >= 29675 & POS <= 29903 ~ "3'UTR") ,
         gene = factor(gene, levels = c("5'UTR", "ORF1ab", "ITR1", "S", "ITR2", "ORF3a", "ITR3", "E",
                                          "ITR4", "M", "ITR5", "ORF6a", "ITR6", "ORF7ab", "ITR7", "ORF8",
                                          "ITR8", "N", "ITR9", "ORF10", "3'UTR"))) %>% 
  mutate(orf1ab = case_when(POS >= 266 & POS <= 805 ~ "nsp1",
                            POS >= 806 & POS <= 2719 ~ "nsp2",
                            POS >= 2720 & POS <= 8554 ~ "nsp3",
                            POS >= 8555 & POS <= 10054 ~ "nsp4",
                            POS >= 10055 & POS <= 10972 ~ "nsp5",
                            POS >= 10973 & POS <= 11842 ~ "nsp6",
                            POS >= 11843 & POS <= 12091 ~ "nsp7",
                            POS >= 12092 & POS <= 12685 ~ "nsp8",
                            POS >= 12686 & POS <= 13024 ~ "nsp9",
                            POS >= 13025 & POS <= 13441 ~ "nsp10",
                            POS >= 13442 & POS <= 16236 ~ "nsp12",
                            POS >= 16237 & POS <= 18039 ~ "nsp13",
                            POS >= 18040 & POS <= 19620 ~ "nsp14",
                            POS >= 19621 & POS <= 20658 ~ "nsp15",
                            POS >= 20659 & POS <= 21552 ~ "nsp16",
                            TRUE ~ NA_character_), .after = gene,
         orf1ab = factor(orf1ab, levels = c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8",
                                            "nsp9", "nsp10", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16")),
         genomic_region = case_when(gene == "ORF1ab" ~ "Replication complex",
                                    gene == "S" ~ "Spike",
                                    gene %in% c("E", "M", "N") ~ "Structural proteins",
                                    gene %in% c("ORF3a", "ORF6a", "ORF7ab", "ORF8", "ORF10") ~ "Accessory proteins",
                                    gene %in% c("5'UTR", "3'UTR", "ITR1", "ITR2", "ITR3", "ITR4", "ITR5","ITR6", "ITR7",
                                                "ITR8", "ITR9") ~ "Non-coding/intergenic",
                                    TRUE ~ NA_character_),
         genomic_region = factor(genomic_region, levels = c("Replication complex", "Spike", "Structural proteins", 
                                                           "Accessory proteins","Non-coding/intergenic"))) %>%
  mutate(voc = case_when(lineage == "B.1" ~ "Early lineages",
                         lineage == "B.1.1" ~ "Early lineages",
                         lineage == "B.1.6" ~ "Early lineages",
                         lineage == "B.1.160" ~ "Early lineages",
                         lineage == "B.1.210" ~ "Early lineages",
                         lineage == "B.1.258" ~ "Early lineages",
                         lineage == "B.1.1.70" ~ "Early lineages",
                         lineage == "B.1.258.17" ~ "Early lineages",
                         lineage == "B.1.177.83" ~ "Early lineages",
                         lineage == "B.1.617.2" ~ "Delta",
                         lineage == "AY.43" ~ "Delta",
                         lineage == "AY.46" ~ "Delta",
                         lineage == "AY.46.6" ~ "Delta",
                         lineage == "AY.98.1" ~ "Delta",
                         lineage == "BA.1" ~ "Omicron",
                         lineage == "BA.1.1" ~ "Omicron",
                         lineage == "BA.2" ~ "Omicron",
                         lineage == "BA.2.56" ~ "Omicron",
                         lineage == "BA.2.9" ~ "Omicron",
                         lineage == "BA.5.1" ~ "Omicron",
                         lineage == "BA.5.2" ~ "Omicron",
                         lineage == "BA.5.2.1" ~ "Omicron",
                         lineage == "BA.5.2.21" ~ "Omicron",
                         lineage == "BA.5.3.1" ~ "Omicron",
                         lineage == "BE.1.1" ~ "Omicron",
                         lineage == "BE.4.1.1" ~ "Omicron",
                         lineage == "BQ.1" ~ "Omicron",
                         lineage == "BQ.1.1" ~ "Omicron",
                         lineage == "BQ.1.1.67" ~ "Omicron",
                         lineage == "BQ.1.1.76" ~ "Omicron",
                         lineage == "XAZ" ~ "Recombinant",
                         lineage == "XBB.1.9" ~ "Recombinant",
                         lineage == "XBB.1.9.1" ~ "Recombinant",
                         lineage == "EG.5.1.6" ~ "Recombinant",
                         lineage == "FL.2" ~ "Recombinant",
                         lineage == "FL.12" ~ "Recombinant"),
         voc = if_else(voc == "Recombinant", "Omicron", voc)) %>% 
  mutate(class = case_when(REF == "A" & ALT == "C" ~ "A>C",
                           REF == "A" & ALT == "G" ~ "A>G",
                           REF == "A" & ALT == "T" ~ "A>T",
                           REF == "C" & ALT == "A" ~ "C>A",
                           REF == "C" & ALT == "G" ~ "C>G",
                           REF == "C" & ALT == "T" ~ "C>T",
                           REF == "G" & ALT == "A" ~ "G>A",
                           REF == "G" & ALT == "C" ~ "G>C",
                           REF == "G" & ALT == "T" ~ "G>T",
                           REF == "T" & ALT == "A" ~ "T>A",
                           REF == "T" & ALT == "C" ~ "T>C",
                           REF == "T" & ALT == "G" ~ "T>G"),
         class_type = if_else(class == "A>G" | class == "G>A" | class == "C>T" | class == "T>C",
                              "Transition", "Transversion"),
         sclass = case_when(class == "A>C" ~ "T>G",
                            class == "A>G" ~ "T>C",
                            class == "A>T" ~ "T>A",
                            class == "G>A" ~ "C>T",
                            class == "G>C" ~ "C>G",
                            class == "G>T" ~ "C>A"),
         sclass = if_else(is.na(sclass), class, sclass),
         across(c(class, class_type, sclass, voc), factor))

summary(mv1)

# de novo mutations
baseline_muts <- mv1 %>%
  filter(sample_n == 1) %>%
  distinct(patient, nt_mut) %>% 
  mutate(.in_baseline = TRUE)

dn <- mv1 %>%
  filter(voc != "Delta") %>% 
  left_join(baseline_muts, by = c("patient", "nt_mut")) %>%
  mutate(de_novo_mut = is.na(.in_baseline)) %>%
  select(-.in_baseline) %>% 
  group_by(patient, sample) %>%
  summarise(
    dn_muts_n = n_distinct(nt_mut[de_novo_mut]),
    dn_major = n_distinct(nt_mut[de_novo_mut & mut_level == "dominant"]),
    dn_minor = n_distinct(nt_mut[de_novo_mut & mut_level == "minor"]),
    dn_nonsyn = n_distinct(nt_mut[de_novo_mut & mut == "Non-synonymous"]),
    dn_syn = n_distinct(nt_mut[de_novo_mut & mut == "Synonymous"]),
    dn_snp = n_distinct(nt_mut[de_novo_mut & mut3 == "SNP"]),
    dn_indel = n_distinct(nt_mut[de_novo_mut & mut3 == "INDEL&NC"]),
    dn_tv = n_distinct(nt_mut[de_novo_mut & class_type == "Transversion"]),
    dn_ts = n_distinct(nt_mut[de_novo_mut & class_type == "Transition"]),
    .groups = "drop")
  
dn1 <- mv1 %>% 
  filter(voc != "Delta") %>% 
  distinct(patient, sample, .keep_all = TRUE) %>% 
  select(patient, sample, sample_n, sample_span, samples_total, duration_total, ct, lineage, voc,
         coverage, meandepth, group, source, vax, th_antineo, th_cst, th_hma, th_antivir) %>% 
  left_join(dn, by = c("patient", "sample")) %>%
  mutate(across(starts_with("dn_"), ~replace_na(.x, 0))) 

summary(dn1)

# figure 1: Longitudinal accumulation of de novo mutations
## figure 1a - observed de novo mutation burden
f1a <- dn1 %>%
  ggplot(aes(sample_span, dn_muts_n)) +
  geom_line(aes(group = patient), alpha = 0.08, linewidth = 0.25) +
  geom_point(aes(color = group), alpha = 0.55, size = 1) +
  scale_color_manual(values = c(ic = "#C43C39", nic = "#2C7FB8"),
                     labels = c(ic = "IC", nic = "nIC")) +
  labs(title = "A", x = "Days since first sequenced sample", y = "De novo mutations", color = "Group") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = ("bottom"),
        legend.box.margin = margin(t = -5),
        legend.margin = margin(t = 0, b = 0),
        legend.spacing.y = unit(0, "pt"),
        plot.margin = margin(t = 5, r = 15, b = 0, l = 5),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 10))
        
f1a

## figure 1b - model-predicted mutation accumulation
### model
#### ic
f1_ic <- glmmTMB(dn_muts_n ~ sample_span * voc + 
                 ct + meandepth + coverage + duration_total + (1 | patient),
                 family = nbinom2,
                 ziformula = ~1,
                 data = dn1 %>% filter(group == "ic"))
summary(f1_ic)

res1_ic <- DHARMa::simulateResiduals(f1_ic)
plot(res1_ic)
DHARMa::testDispersion(res1_ic)
DHARMa::testZeroInflation(res1_ic)

#### sensitivity
f1_group <- glmmTMB(dn_muts_n ~ sample_span * voc + group +
                  ct + meandepth + coverage + duration_total + (1 | patient),
                  family = nbinom2,
                  ziformula = ~1,
                  data = dn1)
summary(f1_group)

res1_group <- DHARMa::simulateResiduals(f1_group)
plot(res1_group)
DHARMa::testDispersion(res1_group)
DHARMa::testZeroInflation(res1_group)  

### plot
f1_terms <- c(
  "sample_span",
  "vocOmicron",
  "sample_span:vocOmicron",
  "ct",
  "meandepth",
  "coverage",
  "duration_total")

d_f1b <- ggpredict(
  f1_ic,
  terms = c("sample_span [all]", "voc"),
  bias_correction = TRUE) %>%
  as_tibble()

summary(d_f1b)

f1b <- d_f1b %>%
  ggplot(aes(x, predicted, color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("darkolivegreen3", "darkorchid3")) +
  scale_fill_manual(values = c("darkolivegreen3", "darkorchid3")) +
  labs(title = "B",
       x = "Days since first sequenced sample",
       y = "Predicted mutation burden",
       color = "VOC",
       fill = "VOC") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = ("bottom"),
        legend.box.margin = margin(t = -5),
        legend.margin = margin(t = 0, b = 0),
        legend.spacing.y = unit(0, "pt"),
        plot.margin = margin(t = 5, r = 15, b = 0, l = 5),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 10))

f1b

## figure 1c - adjusted effects as IRRs
tidy_irr <- function(model, keep_terms = NULL) {
  out <- broom.mixed::tidy(
    model,
    component = "cond",
    effects = "fixed",
    conf.int = TRUE,
    exponentiate = TRUE) %>%
    filter(term != "(Intercept)")
  
  if (!is.null(keep_terms)) {
    out <- out %>% filter(term %in% keep_terms)
  }
  
  out %>%
    mutate(
      irr = estimate,
      irr_low = conf.low,
      irr_high = conf.high,
      term_clean = case_when(
        term == "sample_span" ~ "Time since first sample",
        term == "groupnic" ~ "nIC group",
        term == "sample_span:groupnic" ~ "Time × nIC",
        term == "vocOmicron" ~ "Omicron",
        term == "sample_span:vocOmicron" ~ "Time × Omicron",
        term == "th_cst1" ~ "Corticosteroids",
        term == "th_antivir1" ~ "Antiviral therapy",
        term == "th_cst1:th_antivir1" ~ "Corticosteroids × Antivirals",
        term == "th_antineo1" ~ "Antineoplastic therapy",
        term == "th_hma1" ~ "Human monoclonal antibodies",
        term == "sourcemedication" ~ "Medication-associated\n immunocompromise",
        term == "vaxno" ~ "Unvaccinated",
        term == "vaxpart" ~ "Partially vaccinated",
        term == "ct" ~ "Ct value",
        term == "meandepth" ~ "Mean depth",
        term == "coverage" ~ "Genome coverage",
        term == "duration_total" ~ "Total infection duration",
        TRUE ~ term),
      term_clean = factor(term_clean)
    )
}

### plot
d_f1c <- tidy_irr(f1_ic, f1_terms) %>%
  filter(term != "meandepth") %>%
  mutate(
    term_clean = factor(term_clean,
      levels = c(
        "Time × Omicron",
        "Omicron",
        "Ct value",
        "Genome coverage",
        "Total infection duration",
        "Time since first sample")),
    sig = if_else(p.value < 0.05, "Significant", "Not significant"))

print(d_f1c, n = Inf, width = Inf)
summary(d_f1c)

f1c <- d_f1c %>%
  ggplot(aes(irr, term_clean, colour = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbar(aes(xmin = irr_low, xmax = irr_high), height = 0.2, linewidth = 0.5) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_colour_manual(values = c("Significant" = "#D73027", "Not significant" = "grey40"), guide = "none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = ("bottom"),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 10)) +
  labs(title = "C", x = "Incidence rate ratio (IRR)", y = NULL) 

f1c

fig1 <- (free(f1a) | f1b) / f1c +
  plot_layout(heights = c(1, 0.9))

fig1

ggsave(
  here("figures", "figure1_denovo_mutation_burden.tiff"),
  fig1,
  width = 174,
  height = 117,
  units = "mm",
  dpi = 300,
  compression = "lzw")

## supplementary tables fig 1:
### table a
supp_f1_model <- broom.mixed::tidy(
  f1_ic,
  effects = "fixed",
  component = "cond",
  conf.int = TRUE,
  exponentiate = FALSE) %>%
  mutate(
    IRR = exp(estimate),
    IRR_low = exp(conf.low),
    IRR_high = exp(conf.high)) %>%
  select(
    term,
    estimate,
    std.error,
    IRR,
    IRR_low,
    IRR_high,
    statistic,
    p.value)

write_csv(supp_f1_model, here("supp_tables", "Table_f1a_DeNovoMutationModel.csv"))

### table b
supp_f1_diag <- tibble(
  Metric = c(
    "Distribution",
    "Zero-inflation",
    "Dispersion parameter",
    "DHARMa dispersion test P",
    "DHARMa zero-inflation test P",
    "AIC",
    "BIC",
    "Log-likelihood",
    "Number of patients",
    "Number of samples"),
  Value = c(
    "Negative binomial (nbinom2)",
    "Intercept-only",
    sigma(f1_ic),
    DHARMa::testDispersion(res1_ic)$p.value,
    DHARMa::testZeroInflation(res1_ic)$p.value,
    AIC(f1_ic),
    BIC(f1_ic),
    as.numeric(logLik(f1_ic)),
    length(unique(dn1$patient[dn1$group == "ic"])),
    nobs(f1_ic)))

write_csv(supp_f1_diag, here("supp_tables", "Table_f1b_ModelDiagnostics.csv"))

# figure 2: Evolutionary composition of de novo mutations
## figure 2a - dominant fraction over time
dn_prop <- dn1 %>%
  mutate(
    total_af = dn_major + dn_minor,
    prop_dominant = dn_major / total_af) %>%
  filter(total_af > 0)

print(dn_prop, n = Inf, width = Inf)
summary(dn_prop)

### model
m_af_ic <- glmmTMB(
  cbind(dn_major, dn_minor) ~
    sample_span + voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_prop %>% filter(group == "ic"))

summary(m_af_ic)

res_af_ic <- DHARMa::simulateResiduals(m_af_ic)
plot(res_af_ic)
DHARMa::testDispersion(res_af_ic)
DHARMa::testZeroInflation(res_af_ic)

### sensitivity
m_af_group <- glmmTMB(
  cbind(dn_major, dn_minor) ~
    sample_span * group +
    voc + ct + coverage +
    (1|patient),
  family = betabinomial(),
  data = dn_prop)

summary(m_af_group)

res_af_group <- DHARMa::simulateResiduals(m_af_group)
plot(res_af_group)
DHARMa::testDispersion(res_af_group)
DHARMa::testZeroInflation(res_af_group)

d_f2a <- ggpredict(
  m_af_ic,
  terms = "sample_span [all]",
  bias_correction =  TRUE) %>%
  as_tibble()

summary(d_f2a)

f2a <- ggplot() +
  geom_point(data = dn_prop %>% filter(group == "ic"), aes(sample_span, prop_dominant), alpha = 0.35, size = 1) +
  geom_ribbon(data = d_f2a, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  geom_line(data = d_f2a, aes(x = x, y = predicted), color = "steelblue1", linewidth = 1.1) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = ("bottom"),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "A", x = "Days since first sequenced sample", y = "Dominant\n mutations (%)") 

f2a

## figure 2b - INDEL fraction over time
dn_prop_type <- dn1 %>%
  mutate(
    total_type = dn_indel + dn_snp,
    prop_indel = dn_indel / total_type) %>%
  filter(total_type > 0)

summary(dn_prop_type)

### model
#### ic
m_type_ic <- glmmTMB(
  cbind(dn_indel, dn_snp) ~ sample_span + 
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_prop_type %>% filter(group == "ic"))

summary(m_type_ic)

res_type_ic <- DHARMa::simulateResiduals(m_type_ic)
plot(res_type_ic)
DHARMa::testDispersion(res_type_ic)
DHARMa::testZeroInflation(res_type_ic)

#### sensitivity analysis
m_type_group <- glmmTMB(
  cbind(dn_indel, dn_snp) ~
    sample_span * group +
    voc + ct + coverage +
    (1|patient),
  family = betabinomial(),
  data = dn_prop_type)

summary(m_type_group)

res_type_group <- DHARMa::simulateResiduals(m_type_group)
plot(res_type_group)
DHARMa::testDispersion(res_type_group)
DHARMa::testZeroInflation(res_type_group)

### plot
d_f2b <- ggpredict(
  m_type_ic,
  terms = "sample_span [all]",
  bias_correction = TRUE) %>%
  as_tibble()

summary(d_f2b)

f2b <- ggplot() +
  geom_point(data = dn_prop_type %>% filter(group == "ic"), aes(sample_span, prop_indel), alpha = 0.35, size = 1) +
  geom_ribbon(data = d_f2b, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  geom_line(data = d_f2b, aes(x = x, y = predicted), color = "darkorange2", linewidth = 1.1) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = ("bottom"),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9),
        plot.title = element_text(size = 10)) +
  labs(title = "B", x = "Days since first sequenced sample", y = "INDELS (%)") 

f2b

## figure 2c - nonsynonymous fraction over time
dn_prop_snp <- dn1 %>%
  mutate(
    snp_total = dn_nonsyn + dn_syn,
    prop_nonsyn = dn_nonsyn / snp_total) %>%
  filter(snp_total > 0)

summary(dn_prop_snp)

### model
#### ic
m_snp_ic <- glmmTMB(
  cbind(dn_nonsyn, dn_syn) ~
    sample_span + voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_prop_snp %>% filter(group == "ic"))

summary(m_snp_ic)

res_snp_ic <- DHARMa::simulateResiduals(m_snp_ic)
plot(res_snp_ic)
DHARMa::testDispersion(res_snp_ic)
DHARMa::testZeroInflation(res_snp_ic)

#### sensitivity analysis
m_snp_group <- glmmTMB(
  cbind(dn_nonsyn, dn_syn) ~ sample_span * group +
    voc + ct + coverage +
    (1|patient),
  family = binomial(),
  data = dn_prop_snp)

summary(m_snp_group)

res_snp_group <- DHARMa::simulateResiduals(m_snp_group)
plot(res_snp_group)
DHARMa::testDispersion(res_snp_group)
DHARMa::testZeroInflation(res_snp_group)

### plot
d_f2c <- ggpredict(
  m_snp_ic,
  terms = "sample_span [all]",
  bias_correction = TRUE) %>%
  as_tibble()

summary(d_f2c)

f2c <- ggplot() +
  geom_point(data = dn_prop_snp %>% filter(group == "ic"), aes(sample_span, prop_nonsyn), alpha = 0.35, size = 1) +
  geom_ribbon(data = d_f2c, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  geom_line(data = d_f2c, aes(x = x, y = predicted), color = "palegreen", linewidth = 1.1) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = ("bottom"),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "C", x = "Days since first sequenced sample", y = "Non-synonymous\n SNPs (%)") 

f2c

## figure 2d - transversion fraction over time
dn_prop_tv <- dn1 %>%
  mutate(
    tv_ts_total = dn_ts + dn_tv,
    prop_tv = dn_tv / tv_ts_total) %>%
  filter(tv_ts_total > 0)

summary(dn_prop_tv)

### model
#### ic
m_tvts_ic <- glmmTMB(
  cbind(dn_tv, dn_ts) ~
    sample_span + voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_prop_tv %>% filter(group == "ic"))

summary(m_tvts_ic)

res_tvts_ic <- DHARMa::simulateResiduals(m_tvts_ic)
plot(res_tvts_ic)
DHARMa::testDispersion(res_tvts_ic)
DHARMa::testZeroInflation(res_tvts_ic)

#### sensitivity analysis
m_tvts_group <- glmmTMB(
  cbind(dn_tv, dn_ts) ~ sample_span * group +
    voc + ct + coverage +
    (1|patient),
  family = betabinomial(),
  data = dn_prop_tv)

summary(m_tvts_group)

res_tvts_group <- DHARMa::simulateResiduals(m_tvts_group)
plot(res_tvts_group)
DHARMa::testDispersion(res_tvts_group)
DHARMa::testZeroInflation(res_tvts_group)

### plot
d_f2d <- ggpredict(
  m_tvts_ic,
  terms = "sample_span [all]",
  bias_correction = TRUE) %>%
  as_tibble()

summary(d_f2d)

f2d <- ggplot() +
  geom_point(data = dn_prop_tv %>% filter(group == "ic"), aes(sample_span, prop_tv), alpha = 0.35, size = 1) +
  geom_ribbon(data = d_f2d, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  geom_line(data = d_f2d, aes(x = x, y = predicted), color = "sienna", linewidth = 1.1) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = ("bottom"),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9),
        plot.title = element_text(size = 10)) +
  labs(title = "D", x = "Days since first sequenced sample", y = "Transversions (%)") 

f2d

## figure 2e - effect-size summary
d_f2e <- bind_rows(
  tidy_irr(m_af_ic, "sample_span") %>%
    mutate(model = "Dominant vs minor"),
  tidy_irr(m_type_ic, "sample_span") %>%
    mutate(model = "INDEL vs SNP"),
  tidy_irr(m_snp_ic, "sample_span") %>%
    mutate(model = "Non-synonymous vs synonymous"),
  tidy_irr(m_tvts_ic, "sample_span") %>%
    mutate(model = "Transversion vs Transition")) %>%
  mutate(model = factor(model,
      levels = c(
        "Dominant vs minor",
        "INDEL vs SNP",
        "Non-synonymous vs synonymous",
        "Transversion vs Transition")),
      sig = if_else(p.value < 0.05, "Significant", "Not significant"))

print(d_f2e, n = Inf, width = Inf)
summary(d_f2e)

f2e <- d_f2e %>%
  ggplot(aes(irr, fct_rev(model), color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbar(aes(xmin = irr_low, xmax = irr_high), height = 0.2, orientation = "y") +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_colour_manual(values = c("Significant" = "#D73027", "Not significant" = "grey40"), guide = "none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = ("bottom"),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 10)) +
  labs(title = "E", x = "Odds ratio per day", y = NULL) 

f2e

fig2 <- (free(f2a, side = "l") | f2b) /
  (free(f2c, side = "l") | f2d) / f2e +
  plot_layout(heights = c(1, 1, 0.8),
              widths = c(1, 1),)

fig2

ggsave(
  here("figures", "figure2_mutation_composition.tiff"),
  fig2,
  width = 174,
  height = 117,
  units = "mm",
  dpi = 300,
  compression = "lzw")

## supplementary tables fig 2:
### table a
supp_f2_models <- bind_rows(
  
  tidy_irr(m_af_ic) %>%
    mutate(outcome = "Dominant vs minor"),
  
  tidy_irr(m_type_ic) %>%
    mutate(outcome = "INDEL vs SNP"),
  
  tidy_irr(m_snp_ic) %>%
    mutate(outcome = "Non-synonymous vs synonymous"),
  
  tidy_irr(m_tvts_ic) %>%
    mutate(outcome = "Transversion vs transition")) %>%
  relocate(outcome)

write_csv(supp_f2_models, here("supp_tables", "Table_f2a_MutationCompositionModels.csv"))

### table b
supp_f2_diag <- tibble(
  
  Outcome = c(
    "Dominant vs minor",
    "INDEL vs SNP",
    "Non-synonymous vs synonymous",
    "Transversion vs transition"),
  
  Primary_model = c(
    "Beta-binomial",
    "Beta-binomial",
    "Beta-binomial",
    "Beta-binomial"),
  
  Sensitivity_model = c(
    "Beta-binomial",
    "Beta-binomial",
    "Binomial",
    "Beta-binomial"),
  
  Interaction_p = c(
    summary(m_af_group)$coefficients$cond[
      "sample_span:groupnic","Pr(>|z|)"],
    summary(m_type_group)$coefficients$cond[
      "sample_span:groupnic","Pr(>|z|)"],
    summary(m_snp_group)$coefficients$cond[
      "sample_span:groupnic","Pr(>|z|)"],
    summary(m_tvts_group)$coefficients$cond[
      "sample_span:groupnic","Pr(>|z|)"]),
  
  AIC = c(
    AIC(m_af_ic),
    AIC(m_type_ic),
    AIC(m_snp_ic),
    AIC(m_tvts_ic)),
  
  Dispersion_p = c(
    DHARMa::testDispersion(res_af_ic)$p.value,
    DHARMa::testDispersion(res_type_ic)$p.value,
    DHARMa::testDispersion(res_snp_ic)$p.value,
    DHARMa::testDispersion(res_tvts_ic)$p.value),
  
  ZeroInflation_p = c(
    DHARMa::testZeroInflation(res_af_ic)$p.value,
    DHARMa::testZeroInflation(res_type_ic)$p.value,
    DHARMa::testZeroInflation(res_snp_ic)$p.value,
    DHARMa::testZeroInflation(res_tvts_ic)$p.value))

write_csv(supp_f2_diag, here("supp_tables", "Table_f2b_ModelDiagnostics.csv"))


# figure 3: genomic region composition
region_lengths <- tribble(
  ~genomic_region, ~gen_reg_len,
  "Replication complex", 21290,
  "Spike", 3822,
  "Structural proteins", 2157,
  "Accessory proteins", 1991,
  "Non-coding/intergenic", 643)

dnr <- mv1 %>%
  filter(voc != "Delta") %>%
  left_join(baseline_muts, by = c("patient", "nt_mut")) %>%
  mutate(de_novo_mut = is.na(.in_baseline)) %>%
  filter(de_novo_mut) %>%
  group_by(patient, sample, genomic_region) %>%
  summarise(
    dn_muts_n = n_distinct(nt_mut),
    dn_major = n_distinct(nt_mut[de_novo_mut & mut_level == "dominant"]),
    dn_minor = n_distinct(nt_mut[de_novo_mut & mut_level == "minor"]),
    dn_snp = n_distinct(nt_mut[de_novo_mut & mut3 == "SNP"]),
    dn_indel = n_distinct(nt_mut[de_novo_mut & mut3 == "INDEL&NC"]),
    dn_nonsyn = n_distinct(nt_mut[de_novo_mut & mut == "Non-synonymous"]),
    dn_syn = n_distinct(nt_mut[de_novo_mut & mut == "Synonymous"]),
    dn_tv = n_distinct(nt_mut[de_novo_mut & class_type == "Transversion"]),
    dn_ts = n_distinct(nt_mut[de_novo_mut & class_type == "Transition"]),
    .groups = "drop")

dn2 <- mv1 %>% 
  filter(voc != "Delta") %>% 
  distinct(patient, sample, .keep_all = TRUE) %>% 
  select(patient, sample, sample_n, sample_span, samples_total, duration_total, ct, lineage, voc,
         coverage, meandepth, group, source, vax, th_antineo, th_cst, th_hma, th_antivir) %>% 
  crossing(region_lengths) %>%
  left_join(dnr, by = c("patient", "sample", "genomic_region")) %>%
  mutate(across(starts_with("dn_"), ~replace_na(.x, 0)),
         genomic_region = factor(genomic_region, levels = c("Replication complex", "Spike", "Structural proteins",
                                                           "Accessory proteins", "Non-coding/intergenic"))) %>% 
  mutate(log_gen_reg_len = log(gen_reg_len))

summary(dn2)

## figure 3a - mutation density by genomic region
d_f3a <- dn2 %>%
  filter(group == "ic") %>%
  mutate(muts_per_kb = dn_muts_n / (gen_reg_len / 1000),
         genomic_region = factor(genomic_region, levels = c("Non-coding/intergenic", "Accessory proteins", "Structural proteins",
                                                            "Spike", "Replication complex")))
summary(d_f3a)

f3a <- d_f3a %>%
  ggplot(aes(genomic_region, muts_per_kb, fill = genomic_region)) +
  geom_boxplot(outliers = FALSE) +
  coord_cartesian(ylim = c(0, 5)) +
  scale_fill_manual(values = c("Replication complex" = "#3cb44b", "Spike" = "#4363d8", "Structural proteins" = "#e6194B",
                               "Accessory proteins" = "#ffd8b1", "Non-coding/intergenic" = "snow4")) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "A", x = NULL,  y= "De novo mutation density\n (mutations/kb)") 

f3a

## figure 3b - predicted mutation accumulation by region
### model
#### ic
m_region_count_ic <- glmmTMB(
  dn_muts_n ~ sample_span * genomic_region +
    voc + ct + coverage +
    offset(log_gen_reg_len) +
    (1 | patient),
  ziformula = ~1,
  family = nbinom2,
  data = dn2 %>% filter(group == "ic"))

summary(m_region_count_ic)

res_reg_ic <- DHARMa::simulateResiduals(m_region_count_ic)
plot(res_reg_ic)
DHARMa::testDispersion(res_reg_ic)
DHARMa::testZeroInflation(res_reg_ic)

#### sensitivity
m_region_count_group <- glmmTMB(
  dn_muts_n ~ sample_span * genomic_region + group +
    voc + ct + coverage +
    offset(log_gen_reg_len) +
    (1 | patient),
  ziformula = ~1,
  family = nbinom2,
  data = dn2)

summary(m_region_count_group)

res_reg_group <- DHARMa::simulateResiduals(m_region_count_group)
plot(res_reg_group)
DHARMa::testDispersion(res_reg_group)
DHARMa::testZeroInflation(res_reg_group)

#### interaction
m_region_count_ic_1 <- glmmTMB(
  dn_muts_n ~ sample_span * genomic_region +
    voc + ct + coverage +
    offset(log_gen_reg_len) +
    (1 | patient),
  ziformula = ~1,
  family = nbinom2,
  data = dn2 %>% filter(group == "ic"))

m_region_count_ic_2 <- glmmTMB(
  dn_muts_n ~ sample_span + genomic_region +
    voc + ct + coverage +
    offset(log_gen_reg_len) +
    (1 | patient),
  ziformula = ~1,
  family = nbinom2,
  data = dn2 %>% filter(group == "ic"))

anova(m_region_count_ic_1, m_region_count_ic_2)

### plot
d_f3b <- ggpredict(
  m_region_count_ic,
  terms = c("sample_span [all]", "genomic_region"),
  bias_correction = TRUE) %>%
  as_tibble()

summary(d_f3b)

f3b <- d_f3b %>%
  ggplot(aes(x, predicted, color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.5) +
  scale_fill_manual(values = c("Replication complex" = "#3cb44b", "Spike" = "#4363d8", "Structural proteins" = "#e6194B",
                               "Accessory proteins" = "#ffd8b1", "Non-coding/intergenic" = "snow4")) +
  scale_color_manual(values = c("Replication complex" = "#3cb44b", "Spike" = "#4363d8", "Structural proteins" = "#e6194B",
                               "Accessory proteins" = "#ffd8b1", "Non-coding/intergenic" = "snow4")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "B", x = "Days since first sequenced sample", y = "De novo mutation density\n (mutations/kb)",
       color = "Genomic region", fill = "Genomic region")

f3b

## figure 3c - dominant fraction by region
dn_region_af <- dn2 %>%
  mutate(total_af = dn_major + dn_minor,
         prop_dominant = dn_major / total_af) %>%
  filter(total_af > 0)

summary(dn_region_af)

### model
#### ic
m_region_af_ic <- glmmTMB(
  cbind(dn_major, dn_minor) ~ sample_span * genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_af %>% filter(group == "ic"))

summary(m_region_af_ic)

res_reg_af_ic <- DHARMa::simulateResiduals(m_region_af_ic)
plot(res_reg_af_ic)
DHARMa::testDispersion(res_reg_af_ic)
DHARMa::testZeroInflation(res_reg_af_ic)

#### sensitivity
m_region_af_group <- glmmTMB(
  cbind(dn_major, dn_minor) ~ sample_span * genomic_region + group +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_af)

summary(m_region_af_group)

res_reg_af_group <- DHARMa::simulateResiduals(m_region_af_group)
plot(res_reg_af_group)
DHARMa::testDispersion(res_reg_af_group)
DHARMa::testZeroInflation(res_reg_af_group)

#### interaction
m_region_af_ic_1 <- glmmTMB(
  cbind(dn_major, dn_minor) ~ sample_span * genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_af %>% filter(group == "ic"))

m_region_af_ic_2 <- glmmTMB(
  cbind(dn_major, dn_minor) ~ sample_span + genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_af %>% filter(group == "ic"))

anova(m_region_af_ic_1, m_region_af_ic_2)

### plot
d_f3c <- ggpredict(
  m_region_af_ic,
  terms = c("genomic_region"),
  bias_correction = TRUE) %>%
  as_tibble() %>%  
  mutate(x = factor(x, levels = c("Non-coding/intergenic", "Accessory proteins", "Structural proteins",
                                  "Spike", "Replication complex")))

print(d_f3c, n = Inf, width = Inf)
summary(d_f3c)

f3c <- d_f3c %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_point(aes(color = x), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, alpha = 0.8), width = 0.15) +
  coord_flip() +
  scale_color_manual(values = c("Replication complex" = "#3cb44b", "Spike" = "#4363d8", "Structural proteins" = "#e6194B",
                                "Accessory proteins" = "#ffd8b1", "Non-coding/intergenic" = "snow4")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "C", x = NULL, y = "Predicted dominant mutations (%)") 

f3c

## figure 3d - indel fraction by region
dn_region_type <- dn2 %>%
  filter(genomic_region != "Non-coding/intergenic") %>% 
  mutate(total_type = dn_indel + dn_snp,
         prop_indel = dn_indel / total_type) %>%
  filter(total_type > 0)

summary(dn_region_type)

### model
#### ic
m_region_type_ic <- glmmTMB(
  cbind(dn_indel, dn_snp) ~ sample_span * genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_type %>% filter(group == "ic"))

summary(m_region_type_ic)

res_reg_type_ic <- DHARMa::simulateResiduals(m_region_type_ic)
plot(res_reg_type_ic)
DHARMa::testDispersion(res_reg_type_ic)
DHARMa::testZeroInflation(res_reg_type_ic)

#### sensitivity
m_region_type_group <- glmmTMB(
  cbind(dn_indel, dn_snp) ~ sample_span * genomic_region + group +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_type)

summary(m_region_type_group)

res_reg_type_group <- DHARMa::simulateResiduals(m_region_type_group)
plot(res_reg_type_group)
DHARMa::testDispersion(res_reg_type_group)
DHARMa::testZeroInflation(res_reg_type_group)

#### interaction
m_region_type_ic_1 <- glmmTMB(
  cbind(dn_indel, dn_snp) ~ sample_span * genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_type %>% filter(group == "ic"))

m_region_type_ic_2 <- glmmTMB(
  cbind(dn_indel, dn_snp) ~ sample_span + genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_type %>% filter(group == "ic"))

anova(m_region_type_ic_1, m_region_type_ic_2)

### plot
d_f3d <- ggpredict(
  m_region_type_ic,
  terms = c("genomic_region"),
  bias_correction = TRUE) %>%
  as_tibble() %>%  
  mutate(x = factor(x, levels = c("Non-coding/intergenic", "Accessory proteins", "Structural proteins",
                                  "Spike", "Replication complex")))

print(d_f3d, n = Inf, width = Inf)
summary(d_f3d)

f3d <- d_f3d %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_point(aes(color = x), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, alpha = 0.8), width = 0.15) +
  coord_flip() +
  scale_color_manual(values = c("Replication complex" = "#3cb44b", "Spike" = "#4363d8", "Structural proteins" = "#e6194B",
                                "Accessory proteins" = "#ffd8b1", "Non-coding/intergenic" = "snow4")) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "D", x = NULL, y = "Predicted INDELs (%)") 

f3d

## figure 3e - nonsynonymous fraction by region
dn_region_snp <- dn2 %>%
  mutate(snp_total = dn_nonsyn + dn_syn,
         prop_nonsyn = dn_nonsyn / snp_total) %>%
  filter(snp_total > 0)

summary(dn_region_snp)

### model
#### ic
m_region_snp_ic <- glmmTMB(
  cbind(dn_nonsyn, dn_syn) ~ sample_span * genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_snp %>% filter(group == "ic"))

summary(m_region_snp_ic)

res_reg_snp_ic <- DHARMa::simulateResiduals(m_region_snp_ic)
plot(res_reg_snp_ic)
DHARMa::testDispersion(res_reg_snp_ic)
DHARMa::testZeroInflation(res_reg_snp_ic)

#### sensitivity
m_region_snp_group <- glmmTMB(
  cbind(dn_nonsyn, dn_syn) ~ sample_span * genomic_region + group +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_snp)

summary(m_region_snp_group)

res_reg_snp_group <- DHARMa::simulateResiduals(m_region_snp_group)
plot(res_reg_snp_group)
DHARMa::testDispersion(res_reg_snp_group)
DHARMa::testZeroInflation(res_reg_snp_group)

#### interaction
m_region_snp_ic_1 <- glmmTMB(
  cbind(dn_nonsyn, dn_syn) ~ sample_span * genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_snp %>% filter(group == "ic"))

m_region_snp_ic_2 <- glmmTMB(
  cbind(dn_nonsyn, dn_syn) ~ sample_span + genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_snp %>% filter(group == "ic"))

anova(m_region_snp_ic_1, m_region_snp_ic_2)

### plot
d_f3e <- ggpredict(
  m_region_snp_ic,
  terms = c("genomic_region"),
  bias_correction = TRUE) %>%
  as_tibble() %>%  
  mutate(x = factor(x, levels = c("Non-coding/intergenic", "Accessory proteins", "Structural proteins",
                                  "Spike", "Replication complex")))

print(d_f3e, n = Inf, width = Inf)
summary(d_f3e)

f3e <- d_f3e %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_point(aes(color = x), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, alpha = 0.8), width = 0.15) +
  coord_flip() +
  scale_color_manual(values = c("Replication complex" = "#3cb44b", "Spike" = "#4363d8", "Structural proteins" = "#e6194B",
                                "Accessory proteins" = "#ffd8b1", "Non-coding/intergenic" = "snow4")) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "E", x = NULL, y = "Predicted non-synonymous SNPs (%)")

f3e

## figure 3f - transversion fraction by region
dn_region_tvts <- dn2 %>%
  mutate(total_tvts = dn_tv + dn_ts,
         prop_tv = dn_tv / total_tvts) %>%
  filter(total_tvts > 0)

summary(dn_region_tvts)

### model
#### ic
m_region_tvts_ic <- glmmTMB(
  cbind(dn_tv, dn_ts) ~ sample_span * genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_tvts %>% filter(group == "ic"))

summary(m_region_tvts_ic)

res_reg_tvts_ic <- DHARMa::simulateResiduals(m_region_tvts_ic)
plot(res_reg_tvts_ic)
DHARMa::testDispersion(res_reg_tvts_ic)
DHARMa::testZeroInflation(res_reg_tvts_ic)

#### sensitivity
m_region_tvts_group <- glmmTMB(
  cbind(dn_tv, dn_ts) ~ sample_span * genomic_region + group +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_tvts)

summary(m_region_tvts_group)

res_reg_tvts_group <- DHARMa::simulateResiduals(m_region_tvts_group)
plot(res_reg_tvts_group)
DHARMa::testDispersion(res_reg_tvts_group)
DHARMa::testZeroInflation(res_reg_tvts_group)

#### interaction
m_region_tvts_ic_1 <- glmmTMB(
  cbind(dn_tv, dn_ts) ~ sample_span * genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_tvts %>% filter(group == "ic"))

m_region_tvts_ic_2 <- glmmTMB(
  cbind(dn_tv, dn_ts) ~ sample_span + genomic_region +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn_region_tvts %>% filter(group == "ic"))

anova(m_region_tvts_ic_1, m_region_tvts_ic_2)

### plot
d_f3f <- ggpredict(
  m_region_tvts_ic,
  terms = c("genomic_region"),
  bias_correction = TRUE) %>%
  as_tibble() %>%  
  mutate(x = factor(x, levels = c("Non-coding/intergenic", "Accessory proteins", "Structural proteins",
                                  "Spike", "Replication complex")))

print(d_f3f, n = Inf, width = Inf)
summary(d_f3f)

f3f <- d_f3f %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_point(aes(color = x), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, alpha = 0.8), width = 0.15) +
  coord_flip() +
  scale_color_manual(values = c("Replication complex" = "#3cb44b", "Spike" = "#4363d8", "Structural proteins" = "#e6194B",
                                "Accessory proteins" = "#ffd8b1", "Non-coding/intergenic" = "snow4")) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "F", x = NULL, y = "Predicted transversions (%)") 

f3f

## figure 3
fig3 <- (f3a | f3b) /
  (f3c | f3d) /
  (f3e | f3f) +
  plot_layout(heights = c(1.1, 1, 1))

fig3

ggsave(
  here("figures", "figure3_genomic_region_evolution.tiff"),
  fig3,
  width = 174,
  height = 234,
  units = "mm",
  dpi = 300,
  compression = "lzw")

## supplementary tables fig 3:
tidy_model_full <- function(model, outcome, family_label) {
  
  broom.mixed::tidy(
    model,
    effects = "fixed",
    component = "cond",
    conf.int = TRUE,
    exponentiate = FALSE) %>%
    mutate(
      outcome = outcome,
      family = family_label,
      effect_ratio = exp(estimate),
      effect_ratio_low = exp(conf.low),
      effect_ratio_high = exp(conf.high),
      effect_measure = if_else(
        family_label == "Zero-inflated negative binomial",
        "IRR",
        "OR")) %>%
    select(
      outcome,
      family,
      term,
      estimate,
      std.error,
      statistic,
      p.value,
      effect_measure,
      effect_ratio,
      effect_ratio_low,
      effect_ratio_high)
}
### table 5
supp_table5_coefficients <- bind_rows(
  tidy_model_full(
    m_region_count_ic,
    outcome = "Mutation density",
    family_label = "Zero-inflated negative binomial"),
  tidy_model_full(
    m_region_af_ic,
    outcome = "Dominant vs minor mutation fraction",
    family_label = "Beta-binomial"),
  tidy_model_full(
    m_region_type_ic,
    outcome = "INDEL vs SNP fraction",
    family_label = "Beta-binomial"),
  tidy_model_full(
    m_region_snp_ic,
    outcome = "Non-synonymous vs synonymous SNP fraction",
    family_label = "Beta-binomial"),
  tidy_model_full(
    m_region_tvts_ic,
    outcome = "Transversion vs transition fraction",
    family_label = "Beta-binomial")) %>%
  mutate(
    term_label = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "sample_span" ~ "Time since first sequenced sample",
      term == "genomic_regionSpike" ~ "Spike",
      term == "genomic_regionStructural proteins" ~ "Structural proteins",
      term == "genomic_regionAccessory proteins" ~ "Accessory proteins",
      term == "genomic_regionNon-coding/intergenic" ~ "Non-coding/intergenic",
      term == "vocOmicron" ~ "Omicron",
      term == "ct" ~ "Ct value",
      term == "coverage" ~ "Genome coverage",
      term == "sample_span:genomic_regionSpike" ~
        "Time × Spike",
      term == "sample_span:genomic_regionStructural proteins" ~
        "Time × structural proteins",
      term == "sample_span:genomic_regionAccessory proteins" ~
        "Time × accessory proteins",
      term == "sample_span:genomic_regionNon-coding/intergenic" ~
        "Time × non-coding/intergenic",
      TRUE ~ term)) %>%
  relocate(term_label, .after = term)


write_csv(supp_table5_coefficients, here("supp_tables", "Table_f3_S5_FullModelCoefficients.csv"))


### table 6
extract_lrt <- function(model_full, model_reduced, outcome) {
  
  tab <- as.data.frame(anova(model_reduced, model_full))
  
  tibble(
    outcome = outcome,
    full_model = deparse(formula(model_full)),
    reduced_model = deparse(formula(model_reduced)),
    df_full = tab[2, "Df"],
    df_reduced = tab[1, "Df"],
    AIC_full = tab[2, "AIC"],
    AIC_reduced = tab[1, "AIC"],
    chi_square = tab[2, "Chisq"],
    df_difference = tab[2, "Chi Df"],
    p_value = tab[2, "Pr(>Chisq)"])
}

supp_table6_lrt <- bind_rows(
  extract_lrt(
    m_region_count_ic_1,
    m_region_count_ic_2,
    "Mutation density"),
  extract_lrt(
    m_region_af_ic_1,
    m_region_af_ic_2,
    "Dominant vs minor mutation fraction"),
  extract_lrt(
    m_region_type_ic_1,
    m_region_type_ic_2,
    "INDEL vs SNP fraction"),
  extract_lrt(
    m_region_snp_ic_1,
    m_region_snp_ic_2,
    "Non-synonymous vs synonymous SNP fraction"),
  extract_lrt(
    m_region_tvts_ic_1,
    m_region_tvts_ic_2,
    "Transversion vs transition fraction"))


write_csv(supp_table6_lrt, here("supp_tables", "Table_f3_S6_TimeByRegionLRT.csv"))

## table 7
extract_group_effect <- function(model, outcome, primary_interaction_p,
                                 sensitivity_interaction_p) {
  
  tab <- broom.mixed::tidy(
    model,
    effects = "fixed",
    component = "cond",
    conf.int = TRUE,
    exponentiate = TRUE)
  
  group_row <- tab %>%
    filter(term == "groupnic")
  
  tibble(
    outcome = outcome,
    group_comparison = "nIC vs IC",
    effect_ratio = group_row$estimate,
    ci_low = group_row$conf.low,
    ci_high = group_row$conf.high,
    p_value = group_row$p.value,
    primary_time_by_region_p = primary_interaction_p,
    sensitivity_time_by_region_p = sensitivity_interaction_p,
    conclusion_changed = case_when(
      primary_interaction_p < 0.05 &
        sensitivity_interaction_p >= 0.05 ~ "Yes",
      primary_interaction_p >= 0.05 &
        sensitivity_interaction_p < 0.05 ~ "Yes",
      TRUE ~ "No"))
}

lrt_count_group <- anova(
  update(m_region_count_group, . ~ . - sample_span:genomic_region),
  m_region_count_group)

lrt_af_group <- anova(
  update(m_region_af_group, . ~ . - sample_span:genomic_region),
  m_region_af_group)

lrt_type_group <- anova(
  update(m_region_type_group, . ~ . - sample_span:genomic_region),
  m_region_type_group)

lrt_snp_group <- anova(
  update(m_region_snp_group, . ~ . - sample_span:genomic_region),
  m_region_snp_group)

lrt_tvts_group <- anova(
  update(m_region_tvts_group, . ~ . - sample_span:genomic_region),
  m_region_tvts_group)

lrt_p <- function(x) {
  as.data.frame(x)[2, "Pr(>Chisq)"]
}

supp_table7_sensitivity <- bind_rows(
  extract_group_effect(
    m_region_count_group,
    "Mutation density",
    primary_interaction_p = 0.0004082,
    sensitivity_interaction_p = lrt_p(lrt_count_group)),
  extract_group_effect(
    m_region_af_group,
    "Dominant vs minor mutation fraction",
    primary_interaction_p = 0.07528,
    sensitivity_interaction_p = lrt_p(lrt_af_group)),
  extract_group_effect(
    m_region_type_group,
    "INDEL vs SNP fraction",
    primary_interaction_p = 0.06177,
    sensitivity_interaction_p = lrt_p(lrt_type_group)),
  extract_group_effect(
    m_region_snp_group,
    "Non-synonymous vs synonymous SNP fraction",
    primary_interaction_p = 5.078e-07,
    sensitivity_interaction_p = lrt_p(lrt_snp_group)),
  extract_group_effect(
    m_region_tvts_group,
    "Transversion vs transition fraction",
    primary_interaction_p = 0.1203,
    sensitivity_interaction_p = lrt_p(lrt_tvts_group)))

write_csv(supp_table7_sensitivity, here("supp_tables", "Table_f3_S7_SensitivityAnalyses.csv"))


### table 8
extract_diagnostics <- function(model, residuals, outcome, family_label) {
  
  dispersion_test <- DHARMa::testDispersion(residuals)
  zero_test <- DHARMa::testZeroInflation(residuals)
  
  tibble(
    outcome = outcome,
    family = family_label,
    n_observations = nobs(model),
    n_patients = length(unique(model.frame(model)$patient)),
    AIC = AIC(model),
    BIC = BIC(model),
    log_likelihood = as.numeric(logLik(model)),
    random_intercept_variance =
      as.data.frame(VarCorr(model)$cond)$vcov[1],
    random_intercept_sd =
      as.data.frame(VarCorr(model)$cond)$sdcor[1],
    dispersion_parameter = sigma(model),
    DHARMa_dispersion_statistic =
      unname(dispersion_test$statistic),
    DHARMa_dispersion_p = dispersion_test$p.value,
    DHARMa_zero_ratio =
      unname(zero_test$statistic),
    DHARMa_zero_inflation_p = zero_test$p.value,
    positive_definite_hessian = model$sdr$pdHess
  )
}

supp_table8_diagnostics <- bind_rows(
  extract_diagnostics(
    m_region_count_ic,
    res_reg_ic,
    "Mutation density",
    "Zero-inflated negative binomial"),
  extract_diagnostics(
    m_region_af_ic,
    res_reg_af_ic,
    "Dominant vs minor mutation fraction",
    "Beta-binomial"),
  extract_diagnostics(
    m_region_type_ic,
    res_reg_type_ic,
    "INDEL vs SNP fraction",
    "Beta-binomial"),
  extract_diagnostics(
    m_region_snp_ic,
    res_reg_snp_ic,
    "Non-synonymous vs synonymous SNP fraction",
    "Beta-binomial"),
  extract_diagnostics(
    m_region_tvts_ic,
    res_reg_tvts_ic,
    "Transversion vs transition fraction",
    "Beta-binomial"))

write_csv(supp_table8_diagnostics, here("supp_tables", "Table_f3_S8_ModelDiagnostics.csv"))


# figure 4 - fine-grained mutation density in orf1ab
nsp_lengths <- tribble(
  ~orf1ab, ~nsp_start, ~nsp_end, ~nsp_len,
  "nsp1", 266, 805, 540,
  "nsp2", 806, 2719, 1914,
  "nsp3", 2720, 8554, 5835,
  "nsp4", 8555, 10054, 1500,
  "nsp5", 10055, 10972, 918,
  "nsp6", 10973, 11842, 870,
  "nsp7", 11843, 12091, 249,
  "nsp8", 12092, 12685, 594,
  "nsp9", 12686, 13024, 339,
  "nsp10", 13025, 13441, 417,
  "nsp12", 13442, 16236, 2795,
  "nsp13", 16237, 18039, 1803,
  "nsp14", 18040, 19620, 1581,
  "nsp15", 19621, 20658, 1038,
  "nsp16", 20659, 21552, 894)

dnn <- mv1 %>%
  filter(voc != "Delta") %>%
  left_join(baseline_muts, by = c("patient", "nt_mut")) %>%
  mutate(de_novo_mut = is.na(.in_baseline)) %>%
  filter(de_novo_mut) %>%
  group_by(patient, sample, orf1ab) %>%
  summarise(
    dn_muts_n = n_distinct(nt_mut),
    dn_major = n_distinct(nt_mut[de_novo_mut & mut_level == "dominant"]),
    dn_minor = n_distinct(nt_mut[de_novo_mut & mut_level == "minor"]),
    dn_snp = n_distinct(nt_mut[de_novo_mut & mut3 == "SNP"]),
    dn_indel = n_distinct(nt_mut[de_novo_mut & mut3 == "INDEL&NC"]),
    dn_nonsyn = n_distinct(nt_mut[de_novo_mut & mut == "Non-synonymous"]),
    dn_syn = n_distinct(nt_mut[de_novo_mut & mut == "Synonymous"]),
    dn_tv = n_distinct(nt_mut[de_novo_mut & class_type == "Transversion"]),
    dn_ts = n_distinct(nt_mut[de_novo_mut & class_type == "Transition"]),
    .groups = "drop")

dn3 <- mv1 %>% 
  filter(voc != "Delta") %>% 
  distinct(patient, sample, .keep_all = TRUE) %>% 
  select(patient, sample, sample_n, sample_span, samples_total, duration_total, ct, lineage, voc,
         coverage, meandepth, group, source, vax, th_antineo, th_cst, th_hma, th_antivir) %>% 
  crossing(nsp_lengths) %>%
  left_join(dnn, by = c("patient", "sample", "orf1ab")) %>%
  mutate(across(starts_with("dn_"), ~replace_na(.x, 0)),
         orf1ab = factor(orf1ab, levels = c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8",
                                            "nsp9", "nsp10", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16"))) %>% 
  mutate(orf1ab = relevel(orf1ab, ref = "nsp3")) %>% 
  mutate(log_nsp_len = log(nsp_len))

summary(dn3)

## figure 4a - observed mutation density by nsp
d_f4a <- dn3 %>%
  filter(group == "ic") %>%
  mutate(muts_per_kb = dn_muts_n / (nsp_len / 1000),
         orf1ab = factor(orf1ab, levels = c("nsp16", "nsp15", "nsp14", "nsp13", "nsp12", "nsp10", "nsp9", "nsp8", "nsp7",
                                            "nsp6", "nsp5", "nsp4", "nsp3", "nsp2", "nsp1")),
         nsp_funcmod = case_when(orf1ab %in% c("nsp2", "nsp1") ~ "Host interaction",
                                 orf1ab %in% c("nsp3", "nsp4", "nsp5", "nsp6") ~ "Replication organelle",
                                 orf1ab %in% c("nsp7", "nsp8", "nsp12") ~ "RNA synthesis",
                                 orf1ab %in% c("nsp9", "nsp10", "nsp13", "nsp14", "nsp15", "nsp16") ~ "RNA processing/proofreading"),
         nsp_funcmod = factor(nsp_funcmod, levels = c("Host interaction", "Replication organelle", "RNA synthesis", "RNA processing/proofreading")))


summary(d_f4a)

### plot
f4a <- d_f4a %>%
  ggplot(aes(orf1ab, muts_per_kb, fill = nsp_funcmod)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 6), n.breaks = 6) +
  scale_fill_manual(values = c("Host interaction" = "royalblue1", "Replication organelle" = "plum1", "RNA synthesis" = "springgreen3",
                               "RNA processing/proofreading" = "violetred4")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(hjust = 0.5),
        legend.title.position = "top",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "A", x = NULL, y = "De novo mutations/kb", fill = "NSP functional modules") 

f4a

## figure 4b - predicted mutation accumulation by NSP
### model
#### ic
m_nsp_count_ic <- glmmTMB(
  dn_muts_n ~ sample_span * orf1ab +
    voc + ct + coverage +
    offset(log_nsp_len) +
    (1 | patient),
  ziformula = ~1,
  family = nbinom2,
  data = dn3 %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

summary(m_nsp_count_ic)
m_nsp_count_ic$sdr$pdHess

res_nsp_count_ic <- DHARMa::simulateResiduals(m_nsp_count_ic)
plot(res_nsp_count_ic)
DHARMa::testDispersion(res_nsp_count_ic)
DHARMa::testZeroInflation(res_nsp_count_ic)

#### sensitivity
m_nsp_count_group <- glmmTMB(
  dn_muts_n ~ sample_span * orf1ab + group +
    voc + ct + coverage +
    offset(log(nsp_len)) +
    (1 | patient),
  ziformula = ~1,
  family = nbinom2,
  data = dn3,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

summary(m_nsp_count_group)

res_nsp_count_group <- DHARMa::simulateResiduals(m_nsp_count_group)
plot(res_nsp_count_group)
DHARMa::testDispersion(res_nsp_count_group)
DHARMa::testZeroInflation(res_nsp_count_group)

#### interaction
m_nsp_count_ic_1 <- glmmTMB(
  dn_muts_n ~ sample_span * orf1ab +
    voc + ct + coverage +
    offset(log_nsp_len) +
    (1 | patient),
  ziformula = ~1,
  family = nbinom2,
  data = dn3 %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

m_nsp_count_ic_2 <- glmmTMB(
  dn_muts_n ~ sample_span + orf1ab +
    voc + ct + coverage +
    offset(log_nsp_len) +
    (1 | patient),
  ziformula = ~1,
  family = nbinom2,
  data = dn3 %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

anova(m_nsp_count_ic_1, m_nsp_count_ic_2)

#### pairwise differences
dn3_ic <- dn3 %>%
  filter(group == "ic") %>%
  mutate(
    log_nsp_len = log(nsp_len),
    orf1ab = droplevels(orf1ab))

m_nsp_count_ic_3 <- glmmTMB(
  dn_muts_n ~ sample_span * orf1ab +
    voc + ct + coverage +
    offset(log_nsp_len) +
    (1 | patient),
  ziformula = ~1,
  family = nbinom2,
  data = dn3_ic,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

emm_nsp_count <- emmeans(
  m_nsp_count_ic_3,
  pairwise ~ orf1ab,
  type = "response",
  offset = log(1000),
  adjust = "fdr",
  at = list(
    sample_span = median(dn3_ic$sample_span, na.rm = TRUE),
    ct = median(dn3_ic$ct, na.rm = TRUE),
    coverage = median(dn3_ic$coverage, na.rm = TRUE)))

emm_nsp_count$emmeans
emm_nsp_count$contrasts

### plot
#### lineplot
d_f4bb <- ggpredict(
  m_nsp_count_ic,
  terms = c("sample_span [all]", "orf1ab"),
  condition = c(nsp_len = 1000),
  bias_correction = TRUE) %>%
  as_tibble() %>% 
  mutate(
    group = factor(
      group,
      levels = c(
        "nsp1","nsp2","nsp3","nsp4","nsp5","nsp6",
        "nsp7","nsp8","nsp9","nsp10",
        "nsp12","nsp13","nsp14","nsp15","nsp16")))

f4bb <- d_f4bb %>%
  ggplot(aes(x, predicted, color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.8) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "B", x = "Days since first sequenced sample",y = "Predicted length-adjusted de novo mutations", color = "NSP", fill = "NSP") 

f4bb

#### heatmap
d_f4b <- expand_grid(
  sample_span = seq(
    min(dn3_ic$sample_span, na.rm = TRUE),
    max(dn3_ic$sample_span, na.rm = TRUE),
    by = 1),
  orf1ab = levels(dn3_ic$orf1ab)) %>%
  mutate(
    voc = "Early lineages",
    ct = median(dn3_ic$ct, na.rm = TRUE),
    coverage = median(dn3_ic$coverage, na.rm = TRUE),
    log_nsp_len = log(1000),   # predicted mutations per kb
    orf1ab = factor(orf1ab, levels = levels(dn3_ic$orf1ab)),
    voc = factor(voc, levels = levels(dn3_ic$voc)),
    orf1ab = factor(orf1ab, levels = c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8",
                                       "nsp9", "nsp10", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16"))) 


# Predict over complete grid
pred_f4b <- predict(
  m_nsp_count_ic,
  newdata = d_f4b,
  type = "response",
  se.fit = TRUE,
  re.form = NA)

d_f4b <- d_f4b %>%
  mutate(
    predicted = pred_f4b$fit,
    se = pred_f4b$se.fit,
    conf.low = predicted - 1.96 * se,
    conf.high = predicted + 1.96 * se)

print(d_f4b, n = Inf, width = Inf)
summary(d_f4b)

f4b <- d_f4b %>%
  ggplot(aes(x = sample_span, y = fct_rev(orf1ab), fill = predicted)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", trans = "log10", name = "Predicted\nmutations/kb") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "B",
       x = "Days since first sequenced sample",y = "NSP") 

f4b

## figure 4c - dominant fraction by NSP
dn3_af <- dn3 %>%
  mutate(
    total_af = dn_major + dn_minor,
    prop_dominant = dn_major / total_af) %>%
  filter(total_af > 0)

summary(dn3_af)

### model
#### ic
m_nsp_af_ic <- glmmTMB(
  cbind(dn_major, dn_minor) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_af %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

summary(m_nsp_af_ic)

res_nsp_af_ic <- DHARMa::simulateResiduals(m_nsp_af_ic)
plot(res_nsp_af_ic)
DHARMa::testDispersion(res_nsp_af_ic)
DHARMa::testZeroInflation(res_nsp_af_ic)

#### sensitivity
m_nsp_af_group <- glmmTMB(
  cbind(dn_major, dn_minor) ~ sample_span * orf1ab + group +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_af,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

summary(m_nsp_af_group)

res_nsp_af_group <- DHARMa::simulateResiduals(m_nsp_af_group)
plot(res_nsp_af_group)
DHARMa::testDispersion(res_nsp_af_group)
DHARMa::testZeroInflation(res_nsp_af_group)

#### interaction
m_nsp_af_ic_1 <- glmmTMB(
  cbind(dn_major, dn_minor) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_af %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

m_nsp_af_ic_2 <- glmmTMB(
  cbind(dn_major, dn_minor) ~ sample_span + orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_af %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

anova(m_nsp_af_ic_1, m_nsp_af_ic_2)

#### pairwise differences
dn3_af_ic <- dn3_af %>%
  filter(group == "ic") %>%
  mutate(orf1ab = droplevels(orf1ab))

m_nsp_af_ic_3 <- glmmTMB(
  cbind(dn_major, dn_minor) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_af_ic,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

emm_nsp_af <- emmeans(
  m_nsp_af_ic_3,
  pairwise ~ orf1ab,
  type = "response",
  adjust = "fdr",
  at = list(
    sample_span = median(dn3_af_ic$sample_span, na.rm = TRUE),
    ct = median(dn3_af_ic$ct, na.rm = TRUE),
    coverage = median(dn3_af_ic$coverage, na.rm = TRUE)))

emm_nsp_af$emmeans
emm_nsp_af$contrasts

### plot
d_f4c <- ggpredict(
  m_nsp_af_ic,
  terms = "orf1ab",
  bias_correction = TRUE) %>%
  as_tibble() %>%
  mutate(x = factor(x, levels = c("nsp16", "nsp15", "nsp14", "nsp13", "nsp12", "nsp10", "nsp9", "nsp8", "nsp7",
                                            "nsp6", "nsp5", "nsp4", "nsp3", "nsp2", "nsp1")),
         nsp_funcmod = case_when(x %in% c("nsp2", "nsp1") ~ "Host interaction",
                                 x %in% c("nsp3", "nsp4", "nsp5", "nsp6") ~ "Replication organelle",
                                 x %in% c("nsp7", "nsp8", "nsp12") ~ "RNA synthesis",
                                 x %in% c("nsp9", "nsp10", "nsp13", "nsp14", "nsp15", "nsp16") ~ "RNA processing/proofreading"),
         nsp_funcmod = factor(nsp_funcmod, levels = c("Host interaction", "Replication organelle", "RNA synthesis", "RNA processing/proofreading")))
  

print(d_f4c, n = Inf, width = Inf)
summary(d_f4c)

f4c <- d_f4c %>%
  ggplot(aes(x, predicted)) +
  geom_point(aes(color = nsp_funcmod), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = c("Host interaction" = "royalblue1", "Replication organelle" = "plum1", "RNA synthesis" = "springgreen3",
                               "RNA processing/proofreading" = "violetred4")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "C", x = NULL, y = "Predicted dominant mutations (%)") 

f4c

## figure 4d - INDEL fraction by NSP
dn3_type <- dn3 %>%
  mutate(
    total_type = dn_indel + dn_snp,
    prop_indel = dn_indel / total_type) %>%
  filter(total_type > 0)

summary(dn3_type)

### model
#### ic
m_nsp_type_ic <- glmmTMB(
  cbind(dn_indel, dn_snp) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_type %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

summary(m_nsp_type_ic)

res_nsp_type_ic <- DHARMa::simulateResiduals(m_nsp_type_ic)
plot(res_nsp_type_ic)
DHARMa::testDispersion(res_nsp_type_ic)
DHARMa::testZeroInflation(res_nsp_type_ic)

#### sensitivity
m_nsp_type_group <- glmmTMB(
  cbind(dn_indel, dn_snp) ~ sample_span * orf1ab + group +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_type,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

summary(m_nsp_type_group)

res_nsp_type_group <- DHARMa::simulateResiduals(m_nsp_type_group)
plot(res_nsp_type_group)
DHARMa::testDispersion(res_nsp_type_group)
DHARMa::testZeroInflation(res_nsp_type_group)

#### interaction
m_nsp_type_ic_1 <- glmmTMB(
  cbind(dn_indel, dn_snp) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_type %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

m_nsp_type_ic_2 <- glmmTMB(
  cbind(dn_indel, dn_snp) ~ sample_span + orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_type %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

anova(m_nsp_type_ic_1, m_nsp_type_ic_2)

#### pairwise differences
dn3_type_ic <- dn3_type %>%
  filter(group == "ic") %>%
  mutate(orf1ab = droplevels(orf1ab))

m_nsp_type_ic_3 <- glmmTMB(
  cbind(dn_indel, dn_snp) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_type_ic,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

emm_nsp_type <- emmeans(
  m_nsp_type_ic_3,
  pairwise ~ orf1ab,
  type = "response",
  adjust = "fdr",
  at = list(
    sample_span = median(dn3_type_ic$sample_span, na.rm = TRUE),
    ct = median(dn3_type_ic$ct, na.rm = TRUE),
    coverage = median(dn3_type_ic$coverage, na.rm = TRUE)))

emm_nsp_type$emmeans
emm_nsp_type$contrasts

### plot
d_f4d <- ggpredict(
  m_nsp_type_ic,
  terms = "orf1ab",
  bias_correction = TRUE) %>%
  as_tibble() %>%
  mutate(x = factor(x, levels = c("nsp16", "nsp15", "nsp14", "nsp13", "nsp12", "nsp10", "nsp9", "nsp8", "nsp7",
                                  "nsp6", "nsp5", "nsp4", "nsp3", "nsp2", "nsp1")),
         nsp_funcmod = case_when(x %in% c("nsp2", "nsp1") ~ "Host interaction",
                                 x %in% c("nsp3", "nsp4", "nsp5", "nsp6") ~ "Replication organelle",
                                 x %in% c("nsp7", "nsp8", "nsp12") ~ "RNA synthesis",
                                 x %in% c("nsp9", "nsp10", "nsp13", "nsp14", "nsp15", "nsp16") ~ "RNA processing/proofreading"),
         nsp_funcmod = factor(nsp_funcmod, levels = c("Host interaction", "Replication organelle", "RNA synthesis", "RNA processing/proofreading")))

summary(d_f4d)

f4d <- d_f4d %>%
  ggplot(aes(x, predicted)) +
  geom_point(aes(color = nsp_funcmod), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = c("Host interaction" = "royalblue1", "Replication organelle" = "plum1", "RNA synthesis" = "springgreen3",
                                "RNA processing/proofreading" = "violetred4")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "D", x = NULL, y = "Predicted INDEL mutations (%)") 

f4d

## figure 4e - nonsynonymous SNP fraction by NSP
dn3_snp <- dn3 %>%
  mutate(total_syn = dn_nonsyn + dn_syn) %>%
  filter(total_syn > 0)

summary(dn3_snp)

### model
#### ic
m_nsp_snp_ic <- glmmTMB(
  cbind(dn_nonsyn, dn_syn) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_snp %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

summary(m_nsp_snp_ic)

res_nsp_snp_ic <- DHARMa::simulateResiduals(m_nsp_snp_ic)
plot(res_nsp_snp_ic)
DHARMa::testDispersion(res_nsp_snp_ic)
DHARMa::testZeroInflation(res_nsp_snp_ic)

#### sensitivity
m_nsp_snp_group <- glmmTMB(
  cbind(dn_nonsyn, dn_syn) ~ sample_span * orf1ab + group +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_snp,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

summary(m_nsp_snp_group)

res_nsp_snp_group <- DHARMa::simulateResiduals(m_nsp_snp_group)
plot(res_nsp_snp_group)
DHARMa::testDispersion(res_nsp_snp_group)
DHARMa::testZeroInflation(res_nsp_snp_group)

#### interaction
m_nsp_snp_ic_1 <- glmmTMB(
  cbind(dn_nonsyn, dn_syn) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_snp %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

m_nsp_snp_ic_2 <- glmmTMB(
  cbind(dn_nonsyn, dn_syn) ~ sample_span + orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_snp %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

anova(m_nsp_snp_ic_1, m_nsp_snp_ic_2)

#### pairwise differences
dn3_snp_ic <- dn3_snp %>%
  filter(group == "ic") %>%
  mutate(orf1ab = droplevels(orf1ab))

m_nsp_snp_ic_3 <- glmmTMB(
  cbind(dn_nonsyn, dn_syn) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = betabinomial(),
  data = dn3_snp_ic,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

emm_nsp_snp <- emmeans(
  m_nsp_snp_ic_3,
  pairwise ~ orf1ab,
  type = "response",
  adjust = "fdr",
  at = list(
    sample_span = median(dn3_snp_ic$sample_span, na.rm = TRUE),
    ct = median(dn3_snp_ic$ct, na.rm = TRUE),
    coverage = median(dn3_snp_ic$coverage, na.rm = TRUE)))

emm_nsp_snp$emmeans
emm_nsp_snp$contrasts

### plot
d_f4e <- ggpredict(
  m_nsp_snp_ic,
  terms = "orf1ab",
  bias_correction = TRUE) %>%
  as_tibble() %>%
  mutate(x = factor(x, levels = c("nsp16", "nsp15", "nsp14", "nsp13", "nsp12", "nsp10", "nsp9", "nsp8", "nsp7",
                                  "nsp6", "nsp5", "nsp4", "nsp3", "nsp2", "nsp1")),
         nsp_funcmod = case_when(x %in% c("nsp2", "nsp1") ~ "Host interaction",
                                 x %in% c("nsp3", "nsp4", "nsp5", "nsp6") ~ "Replication organelle",
                                 x %in% c("nsp7", "nsp8", "nsp12") ~ "RNA synthesis",
                                 x %in% c("nsp9", "nsp10", "nsp13", "nsp14", "nsp15", "nsp16") ~ "RNA processing/proofreading"),
         nsp_funcmod = factor(nsp_funcmod, levels = c("Host interaction", "Replication organelle", "RNA synthesis", "RNA processing/proofreading")))

summary(d_f4e)

f4e <- d_f4e %>%
  ggplot(aes(x, predicted)) +
  geom_point(aes(color = nsp_funcmod), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = c("Host interaction" = "royalblue1", "Replication organelle" = "plum1", "RNA synthesis" = "springgreen3",
                                "RNA processing/proofreading" = "violetred4")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "E", x = NULL, y = "Predicted non-synonymous SNPs (%)") 

f4e

## figure 4f - transversion fraction by NSP
dn3_tvts <- dn3 %>%
  mutate(
    total_tvts = dn_tv + dn_ts,
    prop_tv = dn_tv / total_tvts) %>%
  filter(total_tvts > 0)

summary(dn3_tvts)

### model
#### ic
m_nsp_tvts_ic <- glmmTMB(
  cbind(dn_tv, dn_ts) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = binomial(),
  data = dn3_tvts %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

summary(m_nsp_tvts_ic)

res_nsp_tvts_ic <- DHARMa::simulateResiduals(m_nsp_tvts_ic)
plot(res_nsp_tvts_ic)
DHARMa::testDispersion(res_nsp_tvts_ic)
DHARMa::testZeroInflation(res_nsp_tvts_ic)

#### sensitivity
m_nsp_tvts_group <- glmmTMB(
  cbind(dn_tv, dn_ts) ~ sample_span * orf1ab + group +
    voc + ct + coverage +
    (1 | patient),
  family = binomial(),
  data = dn3_tvts,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

summary(m_nsp_tvts_group)

res_nsp_tvts_group <- DHARMa::simulateResiduals(m_nsp_tvts_group)
plot(res_nsp_tvts_group)
DHARMa::testDispersion(res_nsp_tvts_group)
DHARMa::testZeroInflation(res_nsp_tvts_group)

#### interaction
m_nsp_tvts_ic_1 <- glmmTMB(
  cbind(dn_tv, dn_ts) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = binomial(),
  data = dn3_tvts %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

m_nsp_tvts_ic_2 <- glmmTMB(
  cbind(dn_tv, dn_ts) ~ sample_span + orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = binomial(),
  data = dn3_tvts %>% filter(group == "ic"),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

anova(m_nsp_tvts_ic_1, m_nsp_tvts_ic_2)

#### pairwise differences
dn3_tvts_ic <- dn3_tvts %>%
  filter(group == "ic") %>%
  mutate(orf1ab = droplevels(orf1ab))

m_nsp_tvts_ic_3 <- glmmTMB(
  cbind(dn_tv, dn_ts) ~ sample_span * orf1ab +
    voc + ct + coverage +
    (1 | patient),
  family = binomial(),
  data = dn3_tvts_ic,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

emm_nsp_tvts <- emmeans(
  m_nsp_tvts_ic_3,
  pairwise ~ orf1ab,
  type = "response",
  adjust = "fdr",
  at = list(
    sample_span = median(dn3_tvts_ic$sample_span, na.rm = TRUE),
    ct = median(dn3_tvts_ic$ct, na.rm = TRUE),
    coverage = median(dn3_tvts_ic$coverage, na.rm = TRUE)))

emm_nsp_tvts$emmeans
emm_nsp_tvts$contrasts

### plot
d_f4f <- ggpredict(
  m_nsp_tvts_ic,
  terms = c("orf1ab"),
  bias_correction = TRUE) %>%
  as_tibble() %>%
  mutate(x = factor(x, levels = c("nsp16", "nsp15", "nsp14", "nsp13", "nsp12", "nsp10", "nsp9", "nsp8", "nsp7",
                                  "nsp6", "nsp5", "nsp4", "nsp3", "nsp2", "nsp1")),
         nsp_funcmod = case_when(x %in% c("nsp2", "nsp1") ~ "Host interaction",
                                 x %in% c("nsp3", "nsp4", "nsp5", "nsp6") ~ "Replication organelle",
                                 x %in% c("nsp7", "nsp8", "nsp12") ~ "RNA synthesis",
                                 x %in% c("nsp9", "nsp10", "nsp13", "nsp14", "nsp15", "nsp16") ~ "RNA processing/proofreading"),
         nsp_funcmod = factor(nsp_funcmod, levels = c("Host interaction", "Replication organelle", "RNA synthesis", "RNA processing/proofreading")))

summary(d_f4f)

f4f <- d_f4f %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_point(aes(color = nsp_funcmod), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  coord_flip() +
  scale_y_continuous(labels = percent_format()) +
  scale_color_manual(values = c("Host interaction" = "royalblue1", "Replication organelle" = "plum1", "RNA synthesis" = "springgreen3",
                                "RNA processing/proofreading" = "violetred4")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(title = "F", x = NULL, y = "Predicted transversions (%)") 

f4f

## figure 4 combined
fig4 <- (f4a | f4b) /
  (f4c | f4d) /
  (f4e | f4f) +
  plot_layout(heights = c(1.1, 1, 1))

fig4

ggsave(
  here("figures", "figure4_nsp_evolution.tiff"),
  fig4,
  width = 174,
  height = 234,
  units = "mm",
  dpi = 300,
  compression = "lzw")

## accumulation rates
nsp_trends <- emtrends(
  m_nsp_count_ic_3,
  specs = ~ orf1ab,
  var = "sample_span")

nsp_trends

nsp_trend_pairs <- pairs(
  nsp_trends,
  adjust = "fdr")

nsp_trend_pairs

nsp_monthly_change <- summary(
  nsp_trends,
  infer = c(TRUE, TRUE),
  level = 0.95) %>%
  as.data.frame() %>%
  as_tibble() %>% 
  mutate(
    monthly_rate_ratio = exp(30 * sample_span.trend),
    monthly_percent_change =
      100 * (monthly_rate_ratio - 1),
    
    monthly_ci_low =
      100 * (exp(30 * asymp.LCL) - 1),
    
    monthly_ci_high =
      100 * (exp(30 * asymp.UCL) - 1)) %>%
  arrange(desc(monthly_percent_change))

print(nsp_monthly_change, n = Inf, width = Inf)

nsp_trend_pairs_table <- as.data.frame(
  summary(
    nsp_trend_pairs,
    infer = c(TRUE, TRUE),
    level = 0.95)) %>%
  as_tibble()

nsp_trend_pairs_table <- nsp_trend_pairs_table %>%
  mutate(
    monthly_ratio_of_rate_ratios =
      exp(30 * estimate),
    
    monthly_ratio_ci_low =
      exp(30 * asymp.LCL),
    
    monthly_ratio_ci_high =
      exp(30 * asymp.UCL))

print(nsp_trend_pairs_table, n = Inf, width = Inf)

## supplementary tables fig 4:
### table 9
supp_f4_9 <- broom.mixed::tidy(
  m_nsp_count_ic_3,
  effects = "fixed",
  component = "cond",
  conf.int = TRUE,
  exponentiate = FALSE) %>%
  mutate(
    term_label = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "sample_span" ~ "Time since first sequenced sample",
      str_detect(term, "^sample_span:orf1ab") ~
        str_replace(term, "^sample_span:orf1ab", "Time × "),
      str_detect(term, "^orf1ab") ~
        str_remove(term, "^orf1ab"),
      term == "vocOmicron" ~ "Omicron",
      term == "ct" ~ "Ct value",
      term == "coverage" ~ "Genome coverage",
      TRUE ~ term)) %>% 
  mutate(
    IRR = exp(estimate),
    IRR_low = exp(conf.low),
    IRR_high = exp(conf.high),
    term_label = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "sample_span" ~ "Time since first sequenced sample",
      term == "vocOmicron" ~ "Omicron",
      term == "ct" ~ "Ct value",
      term == "coverage" ~ "Genome coverage",
      str_detect(term, "^orf1ab") ~
        str_remove(term, "^orf1ab"),
      str_detect(term, "^sample_span:orf1ab") ~
        str_replace(
          term,
          "^sample_span:orf1ab",
          "Time × "),
      TRUE ~ term)) %>%
  select(
    term,
    term_label,
    estimate,
    std.error,
    statistic,
    p.value,
    IRR,
    IRR_low,
    IRR_high)

write_csv(supp_f4_9, here("supp_tables", "Table_f4_S9_NSPCountModelCoefficients.csv"))


### table 10 (acummulation rates)
supp_f4_9 <- nsp_monthly_change %>%
  transmute(
    NSP = as.character(orf1ab),
    daily_log_slope = sample_span.trend,
    SE = SE,
    z = z.ratio,
    p_value = p.value,
    monthly_rate_ratio = monthly_rate_ratio,
    monthly_percent_change = monthly_percent_change,
    monthly_percent_CI_low = monthly_ci_low,
    monthly_percent_CI_high = monthly_ci_high) %>%
  arrange(desc(monthly_percent_change))

supp_f4_9_formatted <- supp_f4_9 %>%
  mutate(
    daily_log_slope = round(daily_log_slope, 5),
    SE = round(SE, 5),
    monthly_rate_ratio = round(monthly_rate_ratio, 3),
    monthly_percent_change = round(monthly_percent_change, 1),
    monthly_percent_CI_low = round(monthly_percent_CI_low, 1),
    monthly_percent_CI_high = round(monthly_percent_CI_high, 1),
    p_value = scales::pvalue(p_value, accuracy = 0.001))

### table 10
supp_f4_10 <- nsp_trend_pairs_table %>%
  transmute(
    comparison = contrast,
    slope_difference = estimate,
    SE = SE,
    z = z.ratio,
    p_adjusted = p.value,
    monthly_ratio_of_rate_ratios =
      monthly_ratio_of_rate_ratios,
    monthly_ratio_CI_low =
      monthly_ratio_ci_low,
    monthly_ratio_CI_high =
      monthly_ratio_ci_high) %>%
  arrange(p_adjusted)


write_csv(supp_f4_10, here("supp_tables", "Table_f4_S10_NSPAccumulationPairwiseComparisons.csv"))

### table 11
extract_nsp_emmeans <- function(model, outcome) {
  
  emm <- emmeans(
    model,
    specs = ~ orf1ab,
    type = "response"
  )
  
  as.data.frame(summary(emm, infer = c(TRUE, TRUE))) %>%
    as_tibble() %>%
    mutate(outcome = outcome) %>%
    rename_with(
      ~ "predicted",
      any_of(c("prob", "response"))
    ) %>%
    rename_with(
      ~ "CI_low",
      any_of(c("asymp.LCL", "lower.CL"))
    ) %>%
    rename_with(
      ~ "CI_high",
      any_of(c("asymp.UCL", "upper.CL"))
    ) %>%
    select(
      outcome,
      orf1ab,
      predicted,
      SE,
      CI_low,
      CI_high
    )
}

supp_f4_11 <- extract_nsp_emmeans(
  m_nsp_af_ic_3,
  "Dominant mutation fraction") %>%
  mutate(
    predicted_percent = predicted * 100,
    CI_low_percent = CI_low * 100,
    CI_high_percent = CI_high * 100)

write_csv(supp_f4_11, here("supp_tables", "Table_f4_S11_NSPDominantFraction.csv"))

### table 12
supp_f4_12 <- extract_nsp_emmeans(
  m_nsp_type_ic_3,
  "INDEL fraction") %>%
  mutate(
    predicted_percent = predicted * 100,
    CI_low_percent = CI_low * 100,
    CI_high_percent = CI_high * 100)

write_csv(supp_f4_12, here("supp_tables", "Table_f4_S12_NSPINDELFraction.csv"))


### table 13
supp_f4_13_emm <- extract_nsp_emmeans(
  m_nsp_snp_ic_3,
  "Non-synonymous SNP fraction") %>%
  mutate(
    predicted_percent = predicted * 100,
    CI_low_percent = CI_low * 100,
    CI_high_percent = CI_high * 100)

nsp_nonsyn_trends <- emtrends(
  m_nsp_snp_ic_3,
  specs = ~ orf1ab,
  var = "sample_span")

supp_f4_13_trends <- summary(
  nsp_nonsyn_trends,
  infer = c(TRUE, TRUE),
  level = 0.95) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(
    monthly_OR = exp(30 * sample_span.trend),
    monthly_OR_low = exp(30 * asymp.LCL),
    monthly_OR_high = exp(30 * asymp.UCL),
    monthly_odds_percent_change =
      100 * (monthly_OR - 1),
    monthly_odds_percent_CI_low =
      100 * (monthly_OR_low - 1),
    monthly_odds_percent_CI_high =
      100 * (monthly_OR_high - 1)) %>%
  select(
    orf1ab,
    sample_span.trend,
    SE,
    z.ratio,
    p.value,
    monthly_OR,
    monthly_OR_low,
    monthly_OR_high,
    monthly_odds_percent_change,
    monthly_odds_percent_CI_low,
    monthly_odds_percent_CI_high)

supp_f4_13 <- supp_f4_13_emm %>%
  left_join(
    supp_f4_13_trends,
    by = "orf1ab")

write_csv(supp_f4_13, here("supp_tables", "Table_f4_S13_NSPNonsynonymousFractionAndTrends.csv"))


### table 14
supp_f4_14 <- extract_nsp_emmeans(
  m_nsp_tvts_ic_3,
  "Transversion fraction") %>%
  mutate(
    predicted_percent = predicted * 100,
    CI_low_percent = CI_low * 100,
    CI_high_percent = CI_high * 100)

write_csv(supp_f4_14, here("supp_tables", "Table_f4_S14_NSPTransversionFraction.csv"))


### table 15
extract_model_diagnostics <- function(
    model,
    residual_object,
    outcome,
    family_label,
    interaction_p = NA_real_) {
  
  dispersion_test <- DHARMa::testDispersion(
    residual_object,
    plot = FALSE
  )
  
  zero_test <- DHARMa::testZeroInflation(
    residual_object,
    plot = FALSE
  )
  
  vc <- as.data.frame(VarCorr(model)$cond)
  
  tibble(
    outcome = outcome,
    family = family_label,
    observations = nobs(model),
    patients = n_distinct(model.frame(model)$patient),
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = as.numeric(logLik(model)),
    random_intercept_variance = vc$vcov[1],
    random_intercept_SD = vc$sdcor[1],
    dispersion_parameter = sigma(model),
    DHARMa_dispersion_statistic =
      unname(dispersion_test$statistic),
    DHARMa_dispersion_p =
      dispersion_test$p.value,
    DHARMa_zero_ratio =
      unname(zero_test$statistic),
    DHARMa_zero_inflation_p =
      zero_test$p.value,
    interaction_LRT_p = interaction_p,
    positive_definite_hessian =
      model$sdr$pdHess
  )
}

supp_f4_15 <- bind_rows(
  extract_model_diagnostics(
    m_nsp_count_ic_3,
    res_nsp_count_ic,
    "Mutation accumulation",
    "Zero-inflated negative binomial",
    interaction_p = 0.000415),
  extract_model_diagnostics(
    m_nsp_af_ic_3,
    res_nsp_af_ic,
    "Dominant mutation fraction",
    "Beta-binomial",
    interaction_p = 0.00035),
  extract_model_diagnostics(
    m_nsp_type_ic_3,
    res_nsp_type_ic,
    "INDEL fraction",
    "Beta-binomial",
    interaction_p = 0.0012),
  extract_model_diagnostics(
    m_nsp_snp_ic_3,
    res_nsp_snp_ic,
    "Non-synonymous SNP fraction",
    "Beta-binomial",
    interaction_p = 9.9e-9),
  extract_model_diagnostics(
    m_nsp_tvts_ic_3,
    res_nsp_tvts_ic,
    "Transversion fraction",
    "Beta-binomial",
    interaction_p = 0.120))

write_csv(supp_f4_15, here("supp_tables", "Table_f4_S15_NSPModelDiagnostics.csv"))

# figure 5 - clinical drivers
f5_terms <- c(
  "sample_span",
  "groupnic",
  "th_cst1",
  "th_antivir1",
  "th_cst1:th_antivir1",
  "th_antineo1",
  "th_hma1",
  "vaxno",
  "vaxpart",
  "sourcemedication",
  "vocOmicron",
  "duration_total")

### model
m_f5 <- glmmTMB(
  dn_muts_n ~ sample_span * group +
    th_cst + th_antivir + th_cst:th_antivir +
    th_antineo + th_hma +
    vax + source +
    voc + duration_total +
    (1 | patient),
  ziformula = ~1,
  family = nbinom2,
  data = dn1)

summary(m_f5)

res_f5 <- DHARMa::simulateResiduals(m_f5)
plot(res_f5)
DHARMa::testDispersion(res_f5)
DHARMa::testZeroInflation(res_f5)

## figure 5a - full adjusted model
d_f5 <- tidy_irr(m_f5, f5_terms) %>%
  filter(!is.na(irr), !term == "sample_span") %>%
  mutate(
    domain = case_when(
      term %in% c("sample_span", "duration_total") ~ "Time",
      term %in% c("groupnic", "sourcemedication") ~ "Patient status",
      term %in% c("th_cst1", "th_antivir1", "th_cst1:th_antivir1",
                  "th_antineo1", "th_hma1", "vaxno", "vaxpart") ~ "Therapeutic exposure",
      term == "vocOmicron" ~ "VOC",
      term %in% c("ct", "meandepth", "coverage") ~ "Sequencing / sample quality",
      TRUE ~ "Other"),
    domain = factor(domain, levels = c("Therapeutic exposure", "Patient status", "Time", "VOC")))

print(d_f5, n = Inf, width = Inf)
summary(d_f5)

### plot
fig5 <- d_f5 %>%
  ggplot(aes(irr, term_clean)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbar(aes(xmin = irr_low, xmax = irr_high), height = 0.2, orientation = "y", colour = "grey40") +
  geom_point(size = 2) +
  scale_x_log10() +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0),
    legend.position = "none") + 
  labs(x = "Incidence rate ratio", y = NULL) 

fig5

ggsave(
  here("figures", "figure5_clinical_predictors.tiff"),
  fig5,
  width = 174,
  height = 117,
  units = "mm",
  dpi = 300,
  compression = "lzw")

## supplementary tables fig 5:
### table 13
supp_f5_13 <- broom.mixed::tidy(
  m_f5,
  effects = "fixed",
  component = "cond",
  conf.int = TRUE,
  exponentiate = FALSE) %>%
  mutate(
    IRR = exp(estimate),
    IRR_low = exp(conf.low),
    IRR_high = exp(conf.high),
    
    term_label = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "sample_span" ~ "Time since first sequenced sample",
      term == "groupnic" ~ "nIC group",
      term == "th_cst1" ~ "Corticosteroids",
      term == "th_antivir1" ~ "Antiviral therapy",
      term == "th_cst1:th_antivir1" ~
        "Corticosteroids × antiviral therapy",
      term == "th_antineo1" ~ "Antineoplastic therapy",
      term == "th_hma1" ~ "Human monoclonal antibodies",
      term == "vaxno" ~ "Unvaccinated",
      term == "vaxpart" ~ "Partially vaccinated",
      term == "sourcemedication" ~
        "Medication-associated immunocompromise",
      term == "vocOmicron" ~ "Omicron",
      term == "duration_total" ~ "Total infection duration",
      TRUE ~ term),
    
    domain = case_when(
      term %in% c(
        "th_cst1", "th_antivir1",
        "th_cst1:th_antivir1",
        "th_antineo1", "th_hma1",
        "vaxno", "vaxpart") ~ "Therapeutic exposure",
      
      term %in% c(
        "groupnic",
        "sourcemedication") ~ "Patient status",
      
      term %in% c(
        "sample_span",
        "duration_total") ~ "Time",
      
      term == "vocOmicron" ~ "VOC",
      
      TRUE ~ "Other")) %>%
  select(
    domain,
    term,
    term_label,
    estimate,
    std.error,
    statistic,
    p.value,
    IRR,
    IRR_low,
    IRR_high) %>%
  arrange(
    factor(
      domain,
      levels = c(
        "Time",
        "Patient status",
        "Therapeutic exposure",
        "VOC",
        "Other")))

supp_f5_13_formatted <- supp_f5_13 %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 2),
    IRR = round(IRR, 2),
    IRR_low = round(IRR_low, 2),
    IRR_high = round(IRR_high, 2),
    p_value = scales::pvalue(
      p.value,
      accuracy = 0.001))


write_csv(supp_f5_13_formatted, here("supp_tables", "Table_f5_S13_ModelCoefficients.csv"))


### table 14
supp_f5_14 <- tibble(
  Metric = c(
    "Distribution",
    "Zero-inflation",
    "Dispersion parameter",
    "DHARMa dispersion test P",
    "DHARMa zero-inflation test P",
    "AIC",
    "BIC",
    "Log-likelihood",
    "Number of patients",
    "Number of samples"),
  Value = c(
    "Negative binomial (nbinom2)",
    "Intercept-only",
    sigma(m_f5),
    DHARMa::testDispersion(res_f5)$p.value,
    DHARMa::testZeroInflation(res_f5)$p.value,
    AIC(m_f5),
    BIC(m_f5),
    as.numeric(logLik(m_f5)),
    length(unique(dn1$patient[dn1$group == "ic"])),
    nobs(m_f5)))

write_csv(supp_f5_14, here("supp_tables", "Table_f5_S14_ModelDiagnostics.csv"))


# mutational signature
extract_context_left <- function(fasta, sample, pos) {
  samp <- fasta[[sample]]
  left <- pos -1
  return(as.character(subseq(samp, start = left, end = left)))}

extract_context_right <- function(fasta, sample, pos) {
  samp <- fasta[[sample]]
  right <- pos +1
  return(as.character(subseq(samp, start = right, end = right)))}


seq <- readDNAStringSet("data_prep/consensus_seqs_lti_ms.fasta")
#seq

dn4 <- mv1 %>%
  filter(voc != "Delta") %>% 
  left_join(baseline_muts, by = c("patient", "nt_mut")) %>%
  mutate(de_novo_mut = is.na(.in_baseline)) %>%
  select(-.in_baseline) %>% 
  filter(de_novo_mut, mut_type == "SNP")

#summary(dn4)

## mutation signature data prep
sbs_context <- dn4 %>% 
  rowwise() %>% 
  mutate(context_left = extract_context_left(seq, sample, POS),
         context_right = extract_context_right(seq, sample, POS)) %>% 
  ungroup() %>% 
  mutate(sclass = as.character(sclass),
         class = as.character(class),
         context_left_true = if_else(!sclass == class,
                                     case_when(context_left == "A" ~ "T",
                                               context_left == "T" ~ "A",
                                               context_left == "C" ~ "G",
                                               context_left == "G" ~ "C"), context_left),
         middle = if_else(!sclass == class,
                          case_when(REF == "A" ~ "T",
                                    REF == "T" ~ "A",
                                    REF == "C" ~ "G",
                                    REF == "G" ~ "C"), REF),
         context_right_true = if_else(!sclass == class,
                                      case_when(context_right == "A" ~ "T",
                                                context_right == "T" ~ "A",
                                                context_right == "C" ~ "G",
                                                context_right == "G" ~ "C"), context_right)) %>% 
  mutate(component = paste0(context_left_true, middle, context_right_true),
         sccomp = paste(sclass, component, sep = "_")) %>%
  mutate(across(c(class, sclass, context_left, context_right, context_left_true, context_right_true,
                  middle, component, sccomp), factor)) %>% 
  filter(context_left %in% c("A", "C", "G", "T")) %>% 
  filter(context_right %in% c("A", "C", "G", "T")) 

#summary(sbs_context)

all_sbs96 <- expand_grid(
  sclass = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
  left = c("A", "C", "G", "T"),
  right = c("A", "C", "G", "T")) %>%
  mutate(
    middle = str_sub(sclass, 1, 1),
    component = paste0(left, middle, right),
    sccomp = paste(sclass, component, sep = "_")) %>%
  pull(sccomp)

sbs_sample <- sbs_context %>%
  count(patient, sample, group, voc, sample_span, ct, coverage, sccomp, name = "n") %>%
  group_by(patient, sample, group, voc, sample_span, ct, coverage) %>%
  complete(sccomp = all_sbs96, fill = list(n = 0)) %>%
  mutate(
    total_sbs = sum(n),
    prop = n / total_sbs,
    sccomp = factor(sccomp)) %>%
  ungroup()

#summary(sbs_sample)

### clr transformation
sbs_wide <- sbs_sample %>%
  select(sample, sccomp, prop) %>%
  pivot_wider(names_from = sccomp, values_from = prop, values_fill = 0) %>%
  column_to_rownames("sample")

sbs_meta <- sbs_sample %>%
  distinct(sample, patient, group, voc, sample_span, ct, coverage, total_sbs)

sbs_clr <- clr(as.matrix(sbs_wide) + 1e-6)

sbs_clr_df <- as_tibble(sbs_clr, rownames = "sample") %>%
  left_join(sbs_meta, by = "sample")

summary(sbs_clr_df)

## global spectrum test (PERMANOVA)
sbs_meta_adonis <- sbs_meta %>%
  filter(
    !is.na(sample_span),
    !is.na(group),
    !is.na(voc),
    !is.na(ct),
    !is.na(coverage)) %>%
  mutate(
    group = droplevels(factor(group)),
    voc = droplevels(factor(voc)))

sbs_clr_adonis <- sbs_clr_df %>%
  semi_join(sbs_meta_adonis, by = "sample") %>%
  arrange(match(sample, sbs_meta_adonis$sample))

sbs_meta_adonis <- sbs_meta_adonis %>%
  arrange(match(sample, sbs_clr_adonis$sample))

sbs_mat_adonis <- sbs_clr_adonis %>%
  select(all_of(all_sbs96)) %>%
  as.matrix()

rownames(sbs_mat_adonis) <- sbs_clr_adonis$sample

dist_sbs <- dist(sbs_mat_adonis)

adonis_sbs <- adonis2(
  dist_sbs ~ sample_span * group + voc + ct + coverage,
  data = sbs_meta_adonis,
  permutations = 10000,
  by = "margin")

adonis_sbs

adonis2(
  dist_sbs ~ sample_span * group + voc + ct + coverage,
  data = sbs_meta_adonis,
  permutations = 10000)

adonis2(
  dist_sbs ~ sample_span * group + voc + ct + coverage,
  data = sbs_meta_adonis,
  permutations = 10000,
  by = "terms")

dn_sig <- sbs_meta_adonis %>%
  select(sample_span, ct, coverage, voc, group)

cor(
  dn_sig %>%
    select(sample_span, ct, coverage),
  use = "pairwise.complete.obs")

adonis2(
  dist_sbs ~ sample_span,
  data = sbs_meta_adonis)

adonis2(
  dist_sbs ~ sample_span + voc,
  data = sbs_meta_adonis)

adonis2(
  dist_sbs ~ sample_span + voc + ct,
  data = sbs_meta_adonis)

bd <- betadisper(dist_sbs, sbs_meta_adonis$group)
anova(bd)
permutest(bd)

adonis2(
  dist_sbs ~ sample_span + group + voc + ct + coverage,
  data = sbs_meta_adonis)

perm_patient <- how(
  nperm = 10000,
  blocks = sbs_meta_adonis$patient)

adonis_sbs <- adonis2(
  dist_sbs ~ sample_span * group + voc + ct + coverage,
  data = sbs_meta_adonis,
  permutations = perm_patient,
  by = "margin")

sbs_meta_cap <- sbs_meta %>%
  filter(
    !is.na(sample_span),
    !is.na(group),
    !is.na(voc),
    !is.na(ct),
    !is.na(coverage))

idx <- match(sbs_meta_cap$sample, rownames(sbs_clr))

sbs_clr_cap <- sbs_clr[idx, ]

rda_fit <- rda(
  sbs_clr_cap ~ sample_span + group + voc + ct + coverage,
  data = sbs_meta_cap)

anova(rda_fit, permutations = perm_patient)

anova(rda_fit, by = "term", permutations = perm_patient)

anova(rda_fit, by = "axis", permutations = perm_patient)

RsquareAdj(rda_fit)

species_scores <- scores(rda_fit, display = "species")

species_df <- as.data.frame(species_scores) %>%
  rownames_to_column("context") %>%
  mutate(
    loading = sqrt(RDA1^2 + RDA2^2)) %>%
  arrange(desc(loading))
  
## visualization
rda_sites <- scores(rda_fit, display = "sites", scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(sbs_meta_cap, by = "sample") %>%
  mutate(group = if_else(group == "ic", "IC", "nIC"))

summary(rda_sites)

### biplot arrows for explanatory variables
rda_bp <- scores(rda_fit, display = "bp", scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  mutate(variable = case_when(variable == "sample_span" ~ "Days since\n first sequenced sample",
                              variable == "groupnic" ~ "Group nIC",
                              variable == "vocOmicron" ~ "VOC Omicron",
                              variable == "ct" ~ "Ct value",
                              variable == "coverage" ~ "Coverage"))

rda_var <- summary(rda_fit)$cont$importance[2, 1:2] * 100
rda_var

f_sig_a <- rda_sites %>%
  ggplot(aes(RDA1, RDA2)) +
  geom_point(aes(color = sample_span, shape = group), alpha = 0.4, size = 1.5) +
  geom_segment(data = rda_bp, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), inherit.aes = FALSE, arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.4) +
  geom_text_repel(data = rda_bp, aes(RDA1, RDA2, label = variable), inherit.aes = FALSE, size = 3.5) +
  scale_color_viridis_c(name = "Days since first sample", guide = guide_colorbar(barwidth = unit(4, "cm"), barheight = unit(0.25, "cm"), title.position = "left")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.box.margin = margin(t = -5),
        legend.margin = margin(t = 0, b = 0),
        legend.spacing.y = unit(0, "pt"),
        plot.margin = margin(t = 5, r = 15, b = 0, l = 5)) +
  labs(title = "A", x = paste0("RDA1 (", round(rda_var[1], 1), "%)"), y = paste0("RDA2 (", round(rda_var[2], 1), "%)"), shape = "Group")

f_sig_a

rda_loadings <- scores(rda_fit, display = "species", scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column("sccomp") %>%
  mutate(
    loading_RDA1_RDA2 = sqrt(RDA1^2 + RDA2^2),
    sclass = str_extract(sccomp, "^[ACTG]>[ACTG]"),
    component = str_extract(sccomp, "[ACTG]{3}$")) %>%
  arrange(desc(loading_RDA1_RDA2))

top_contexts <- rda_loadings %>%
  slice_head(n = 20)

top_contexts

f_sig_b <- top_contexts %>%
  mutate(sccomp = fct_reorder(sccomp, loading_RDA1_RDA2, .desc = TRUE),
         sclass = factor(sclass, levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))) %>%
  ggplot(aes(sccomp, loading_RDA1_RDA2, fill = sclass)) +
  geom_col() +
  scale_fill_manual(values = c("C>A" = "turquoise3", "C>G" = "black", "C>T" = "red3", 
                               "T>A" = "snow3", "T>C" = "olivedrab1", "T>G" = "pink1"),
                    breaks = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, size = 6, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10),
        legend.position = "none",
        legend.box.margin = margin(t = -5),
        legend.margin = margin(t = 0, b = 0),
        legend.spacing.y = unit(0, "pt"),
        plot.margin = margin(t = 0, r = 15, b = 0, l = 5)) +
  labs(title = "B", x = NULL, y = "Loading magnitude\n on RDA1/RDA2", fill = "Substitution")

f_sig_b

sbs_profile_plot <- sbs_sample %>%
  semi_join(sbs_meta_cap, by = "sample") %>%
  select(-any_of("voc")) %>%
  left_join(sbs_meta_cap %>%
              select(sample, voc), by = "sample") %>%
  group_by(voc, sccomp) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE) * 100, .groups = "drop") %>%
  mutate(sclass = str_extract(sccomp, "^[ACTG]>[ACTG]"),
         component = str_extract(sccomp, "[ACTG]{3}$"),
         sclass = factor(sclass, levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")))

print(sbs_profile_plot, n = Inf, width = Inf)
summary(sbs_profile_plot)

sbs_voc_tests <- sbs_clr_df %>%
  semi_join(sbs_meta_cap, by = "sample") %>%
  pivot_longer(
    cols = all_of(all_sbs96),
    names_to = "sccomp",
    values_to = "clr_value") %>%
  group_by(sccomp) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lm(
      clr_value ~ voc + sample_span + group + ct + coverage,
      data = .x)),
    tidy = map(fit, broom::tidy)) %>%
  select(sccomp, tidy) %>%
  unnest(tidy) %>%
  filter(term == "vocOmicron") %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BY"),
    direction = case_when(
      estimate > 0 ~ "enriched",
      estimate < 0 ~ "depleted"),
    label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ NA_character_))

print(sbs_voc_tests, n = Inf, width = Inf)
summary(sbs_voc_tests)

sbs_annot <- sbs_voc_tests %>%
  filter(!is.na(label)) %>% 
  left_join(sbs_profile_plot %>% filter(voc == "Omicron"), by = "sccomp") %>%
  mutate(
    sclass = str_extract(sccomp, "^[ACTG]>[ACTG]"),
    component = str_extract(sccomp, "[ACTG]{3}$"),
    sclass = factor(
      sclass,
      levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")),
    voc = "Omicron")

f_sig_c <- sbs_profile_plot %>%
  ggplot(aes(component, mean_prop, fill = sclass)) +
  geom_col(width = 0.85) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 7.8, ymax = 8.5, fill = sclass, alpha = 0.5)) +
  geom_text(data = sbs_annot, aes(x = component, y = (mean_prop + 0.3), label = label, color = direction), inherit.aes = FALSE, size = 0.8, fontface = "bold") +
  facet_grid(voc ~ sclass, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("C>A" = "turquoise3", "C>G" = "black", "C>T" = "red3", 
                               "T>A" = "snow3", "T>C" = "olivedrab1", "T>G" = "pink1")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, size = 4),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "none",
    legend.box.margin = margin(t = -5),
    legend.margin = margin(t = 0, b = 0),
    legend.spacing.y = unit(0, "pt"),
    plot.margin = margin(t = 0, r = 15, b = 0, l = 5)) +
  labs(title = "C", x = "SBS96 context", y = "Mean SBS96 proportion (%)", fill = "Substitution")

f_sig_c

sample_order <- sbs_meta_cap %>%
  arrange(sample_span) %>%
  pull(sample)

top_sbs <- top_contexts$sccomp

sbs_heat <- sbs_clr_df %>%
  mutate(sample = factor(sample, levels = sample_order)) %>%
  select(sample, all_of(top_sbs)) %>%
  pivot_longer(
    cols = all_of(top_sbs),
    names_to = "sccomp",
    values_to = "clr_value") %>%
  left_join(sbs_meta_cap, by = "sample") %>%
  mutate(sccomp = factor(sccomp, levels = c("T>G_TTC", "T>G_TTG", "T>G_TTT",
                                            "T>C_TTG", "T>C_GTG", "T>C_CTA", "T>C_ATA", "T>C_TTA",
                                            "T>A_TTA", "T>A_GTT", "T>A_TTT", "T>A_ATC", "T>A_GTA", "T>A_GTG", "T>A_TTC",
                                            "C>T_CCT", "C>T_GCC",
                                            "C>A_CCA", "C>A_GCA", "C>A_TCT")))

f_sig_d <- sbs_heat %>%
  filter(group == "ic") %>%
  drop_na(group) %>% 
  ggplot(aes(sample, sccomp, fill = clr_value)) +
  geom_tile() +
  #facet_grid(group ~ voc, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(name = "SBS96 context abundance", guide = guide_colorbar(barwidth = unit(4, "cm"), barheight = unit(0.25, "cm"), title.position = "left")) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 6, face = "bold"),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(size = 10),
    legend.position = "bottom",
    legend.box.margin = margin(t = -5),
    legend.margin = margin(t = 0, b = 0),
    legend.spacing.y = unit(0, "pt"),
    plot.margin = margin(t = 0, r = 15, b = 0, l = 5)) +
  labs(title = "D", x = "Days since first sequenced sample", y = "SBS96 context")

f_sig_d


fig6 <- f_sig_a / f_sig_b / f_sig_c / f_sig_d +
  plot_layout(heights = c(1, 1, 1.5, 1.5))

fig6

ggsave(
  here("figures", "figure6_mutational_signatures.tiff"),
  fig6,
  width = 174,
  height = 234,
  units = "mm",
  dpi = 300,
  compression = "lzw")

## supplementary tables fig 6:
### table 15
rda_overall <- as.data.frame(
  anova(rda_fit, permutations = perm_patient)) %>%
  rownames_to_column("effect") %>%
  as_tibble() %>%
  mutate(
    analysis = "RDA overall model",
    R2 = if_else(
      effect == "Model",
      RsquareAdj(rda_fit)$r.squared,
      NA_real_),
    adjusted_R2 = if_else(
      effect == "Model",
      RsquareAdj(rda_fit)$adj.r.squared,
      NA_real_))

rda_terms <- as.data.frame(
  anova(
    rda_fit,
    by = "term",
    permutations = perm_patient)) %>%
  rownames_to_column("effect") %>%
  as_tibble() %>%
  mutate(
    analysis = "RDA term test")

rda_axes <- as.data.frame(
  anova(
    rda_fit,
    by = "axis",
    permutations = perm_patient)) %>%
  rownames_to_column("effect") %>%
  as_tibble() %>%
  mutate(
    analysis = "RDA axis test")

supp_f6_15 <- bind_rows(
  rda_overall,
  rda_terms,
  rda_axes) 


write_csv(supp_f6_15, here("supp_tables", "Table_f6_S15_RDATestsSBS96.csv"))


### table 16
supp_f6_16 <- scores(
  rda_fit,
  display = "species",
  scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column("SBS96_context") %>%
  as_tibble() %>%
  mutate(
    substitution_class =
      str_extract(SBS96_context, "^[ACTG]>[ACTG]"),
    
    trinucleotide_context =
      str_extract(SBS96_context, "[ACTG]{3}$"),
    
    loading_magnitude_RDA1_RDA2 =
      sqrt(RDA1^2 + RDA2^2),
    
    loading_rank =
      min_rank(loading_magnitude_RDA1_RDA2)) %>%
  arrange(loading_rank)


write_csv(supp_f6_16, here("supp_tables", "Table_f6_S16_SBS96RDALoadings.csv"))

### table 17
supp_f6_17 <- sbs_voc_tests %>%
  left_join(
    sbs_profile_plot %>%
      select(
        voc,
        sccomp,
        mean_prop) %>%
      pivot_wider(
        names_from = voc,
        values_from = mean_prop,
        names_prefix = "mean_percent_"),
    by = "sccomp") %>%
  transmute(
    SBS96_context = sccomp,
    substitution_class =
      str_extract(sccomp, "^[ACTG]>[ACTG]"),
    trinucleotide_context =
      str_extract(sccomp, "[ACTG]{3}$"),
    
    estimate_CLR = estimate,
    standard_error = std.error,
    t_statistic = statistic,
    raw_p_value = p.value,
    BY_adjusted_p_value = p_adj,
    direction_in_Omicron = direction,
    
    mean_percent_early_lineages =
      `mean_percent_Early lineages`,
    
    mean_percent_Omicron =
      mean_percent_Omicron) %>%
  arrange(BY_adjusted_p_value)


write_csv(supp_f6_17, here("supp_tables", "Table_f6_S17_OmicronVsEarlySBS96.csv"))


### table 18
supp_f6_18 <- sbs_sample %>%
  semi_join(
    sbs_meta_cap,
    by = "sample") %>%
  group_by(
    voc,
    sccomp) %>%
  summarise(
    n_samples = n_distinct(sample),
    
    mean_proportion_percent =
      mean(prop, na.rm = TRUE) * 100,
    
    median_proportion_percent =
      median(prop, na.rm = TRUE) * 100,
    
    Q1_proportion_percent =
      quantile(
        prop,
        0.25,
        na.rm = TRUE) * 100,
    
    Q3_proportion_percent =
      quantile(
        prop,
        0.75,
        na.rm = TRUE) * 100,
    
    samples_with_context_n =
      sum(n > 0, na.rm = TRUE),
    
    samples_with_context_percent =
      100 * mean(n > 0, na.rm = TRUE),
    
    total_mutations =
      sum(n, na.rm = TRUE),
    
    .groups = "drop") %>%
  mutate(
    substitution_class =
      str_extract(sccomp, "^[ACTG]>[ACTG]"),
    
    trinucleotide_context =
      str_extract(sccomp, "[ACTG]{3}$"))


write_csv(supp_f6_18, here("supp_tables", "Table_f6_S18_DescriptiveSBS96Profiles.csv"))
