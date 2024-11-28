####
#
#
# bats analysis
#
#
####


## packages ----
library(tidyverse)
library(conflicted)
library(viridis)
library(lme4)
library(vegan)
library(phyloseq)
library(phyloseqCompanion)
library(ggordiplots)
library(metagMisc)
library(microbiome)
library(ecodist)
library(fantaxtic)
library(ALDEx2)
library(janitor)
library(ggmap)
theme_set(theme_classic())
conflicts_prefer(dplyr::filter())
conflicts_prefer(dplyr::select)

## data input ----
feature_table = read_csv("feature_table.csv")
taxonomy_table = read_csv("taxonomy_table.csv")
mcmaster_metadata = read_csv("mcmaster_metadata.csv")
webber_metadata = read_csv("webber_metadata.csv")

## _ data cleanup ----
# remove "KCGneg" column from feature_table and add rownames
feature_table = feature_table |> 
  select(-KCGneg) |> 
  column_to_rownames(var = "...1") #|> t() #to transpose (not needed anymore)

# check date format
str(webber_metadata$date) #character
str(mcmaster_metadata$date) #date, format: "2023-06-05"

# fix webber_metadata dates
webber_metadata = webber_metadata |> 
  mutate(date = as.Date(date, format = "%d-%B-%y"))
#check
str(webber_metadata$date) #date, format: "2023-06-05". Worked!

# merge metadata
bat_metadata = mcmaster_metadata |> 
  left_join(webber_metadata, by = c("mcmaster_sample_id", "date"))


## _ create phyloseq ----
otu = feature_table
otu = otu_table(otu, taxa_are_rows = TRUE)
tax = taxonomy_table |> column_to_rownames("...1") #set row names
tax = tax_table(as.matrix(tax))
samples = bat_metadata |> column_to_rownames("study_id")
samples = sample_data(samples)

webber_bats = phyloseq(otu, tax, samples)


# _ alpha div data ----
adiv = plot_richness(webber_bats, measures=c("Observed", "Chao1", "Shannon", "Simpson"), x = "site", color = "colony_size")

# store alpha diversity measures as new variables
alphadiv = data.table(adiv$data)

# create tibble with adiv measures for each sample
alphadiv = alphadiv |> 
  select(mcmaster_sample_id, samples, variable, value) |> 
  pivot_wider(names_from = variable, values_from = value) |> 
  mutate(mcmaster_sample_id = as.character(mcmaster_sample_id))

# bind to existing metadata
alphadiv$mcmaster_sample_id = as.numeric(alphadiv$mcmaster_sample_id)

# merge to metadata
bat_metadata = bat_metadata |> 
  left_join(alphadiv, by = "mcmaster_sample_id")

# _ subset BBB and LBB ----
bbb <- subset_samples(webber_bats, species == "EPFU")
lbb <- subset_samples(webber_bats, species == "MYLU")

df_otus_bbb = otu.matrix(bbb)
bbb_metadata = as.data.frame(as.matrix(sample_data(bbb)))


## descriptive stats ----

# _ sample sizes ----
# number of samples per site
bat_metadata |> group_by(site) |> summarise(n = n())
# 1 Alliston          17
# 2 Dunnville          8
# 3 Langton            2
# 4 Lion's Head (lbb)  13
# 5 Mount Pleasant     6
# 6 Wainfleet         20

# number of sites for EPFU
bat_metadata |> filter(species == "EPFU") |> distinct(site) #5 sites (MYLU = 1 site)

# number of adults/juveniles
bat_metadata |> filter(species == "EPFU") |> group_by(age) |> summarise(n = n()) #40 adults, 13 juveniles
bat_metadata |> filter(species == "MYLU") |> group_by(age) |> summarise(n = n()) #all 13 are adults

# number of samples per colony size
bat_metadata |> group_by(colony_size) |> summarise(n = n()) #large: 35; small: 31
#there are 3 large colonies and 3 small colonies

# number of samples per colony size per species
bat_metadata |> filter(species == "EPFU") |> group_by(colony_size) |> summarise(n = n()) #large: 22; small: 31

bat_metadata |> filter(species == "MYLU") |> group_by(colony_size) |> summarise(n = n()) #large: 13; small: 0

# number of reads ----
df_otus_bats |> 
  sum()

# _ ASVs ----
feature_table_over2 = feature_table[-(which(rowSums(feature_table)<3)),]
feature_table_over5 = feature_table[-(which(rowSums(feature_table)<6)),]
## number of ASVs in total: 2094
## number of ASVs >2 reads: 1978, >5 reads: 1626

# number of ASVs shared between EPFU and MYLU
# subset a smaller dataset with just the samples you want
common_otus = subset_samples(webber_bats, species == "EPFU" | species == "MYFU")

# remove all ASVs not found at least once in both samples
common_otus = filter_taxa(common_otus, function(x) sum(x >= 1) == (2), TRUE)
# number of common ASVs: 208

# create feature table subsets for BBB and LBB
samples_bbb = bat_metadata |> filter(species == "EPFU") |> pull(study_id)
samples_lbb = bat_metadata |> filter(species == "MYLU") |> pull(study_id)
feature_table_bbb = feature_table |> select(all_of(samples_bbb))
feature_table_lbb = feature_table |> select(all_of(samples_lbb))

# get total number of ASVs for species separately and remove ASVs with zeroes everywhere
feature_table_bbb = feature_table_bbb[-(which(rowSums(feature_table_bbb) == 0)),] #ASVs: 1939
feature_table_lbb = feature_table_lbb[-(which(rowSums(feature_table_lbb) == 0)),] #ASVs: 254

## number of ASVs per sample
# first, calculate the number of non-zero rows for each column:
# this is the number of observed ASVs in each sample
total_asvs_bbb = feature_table_bbb %>% 
  summarise(across(starts_with("KCG"), ~ sum(.x != 0))) |> 
  pivot_longer(cols = everything()) 
  # now name = sample, value = total ASVs

total_asvs_lbb = feature_table_lbb %>% 
  summarise(across(starts_with("KCG"), ~ sum(.x != 0))) |> 
  pivot_longer(cols = everything()) 
# now name = sample, value = total ASVs

total_asvs_bbb |> summarise(mean = mean(value), sd = sd(value), min = min(value), max = max(value))
# mean = 95, sd = 56, min = 35, max = 261

total_asvs_lbb |> summarise(mean = mean(value), sd = sd(value), min = min(value), max = max(value))
# mean = 74.5, sd = 15, min = 54, max = 104


# _ sequencing depth ----
depth_mean = mean(rowSums(df_otus_bats)) #65407.73
depth_sd = sd(rowSums(df_otus_bats)) #35057.05

hist(sample_sums(webber_bats), main = "Histogram: Read Counts", xlab = "Total Reads", las = 1, breaks = 40)

# plot the distribution of sequencing depth
sdt = data.table(as(sample_data(webber_bats), "data.frame"), TotalReads = sample_sums(webber_bats), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth

pSeqDepth + facet_wrap(~species)
# no obvious difference between species


## normalize number of reads using media sequencing depth
total = median(sample_sums(webber_bats))
standf = function(x, t=total) round(t * (x / sum(x)))
webber_bats_norm = transform_sample_counts(webber_bats, standf)
bbb_norm = transform_sample_counts(bbb, standf)
lbb_norm = transform_sample_counts(lbb, standf)


## abundance graphs ----
# general stats
# taxa counts
tax_count = data.table(tax_table(webber_bats),
                 TotalCounts = taxa_sums(webber_bats),
                 OTU = taxa_names(webber_bats))

ggplot(tax_count, aes(TotalCounts)) + 
  geom_histogram(bins = 50) + 
  ggtitle("Histogram of Total Counts")

# abundances
# per species
plot_bar(webber_bats, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  facet_grid(~ species) + theme(legend.position = "none") +
  ggtitle("Relative abundances in big and little brown bats")

# per site
plot_bar(webber_bats_norm, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  facet_grid(~ site) + theme(legend.position = "none")

# basic bat graph based on Phylum for normalized counts
plot_bar(webber_bats_norm, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") + facet_grid(~ species)

# or plot each age group separately
plot_bar(bbb_norm, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  ggtitle("Phylum-level abundances in BBB") #+ theme(legend.position = "none") + facet_wrap(~ site)

plot_bar(lbb_norm, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  ggtitle("Phylum-level abundances in LBB")

# heatmap
plot_heatmap(webber_bats, method = "NMDS", distance = "bray", title = "Heatmap of OTUs in BBB and LBB")

# heatmap #2: abundant OTUs (representing at least 20% of reads in at least 1 sample)
total = median(sample_sums(webber_bats))
webber_bats_abund = filter_taxa(webber_bats, function(x) sum(x > total*0.20) > 0, TRUE)
plot_heatmap(webber_bats_abund, method = "NMDS", distance = "bray") + theme(axis.text.y = element_blank()) +
  ggtitle("Abundant OTUs (≥ 20% of reads in ≥ 1 sample")

## dominant taxa ----
# top 3 most abundant phyla
fantaxtic::top_taxa(webber_bats, n = 3, tax_level = "Phylum")
# top 5 most abundant orders
fantaxtic::top_taxa(webber_bats, n = 5, tax_level = "Order")

# dominant phyla in BBB
fantaxtic::top_taxa(bbb, n = 3, tax_level = "Phylum")
#56% Firmicutes, 37% Proteobacteria, 2% Verrucomicrobiota
fantaxtic::top_taxa(lbb, n = 3, tax_level = "Phylum")
#60% Firmicutes, 27% Proteobacteria, 2% Actinobacteriota

# dominant classes and genera in BBB and LBB
fantaxtic::top_taxa(bbb, n = 3, tax_level = "Class")
#50% Bacilli, 35% Gammaproteobacteria, 10% Clostridia
fantaxtic::top_taxa(lbb, n = 3, tax_level = "Class")
#27% Bacilli, 27% Gammaproteobacteria, 42% Clostridia

fantaxtic::top_taxa(bbb, n = 3, tax_level = "Genus")
#25% Enterococcus, 10% Klebsiella, 10% Lactococcus
fantaxtic::top_taxa(lbb, n = 3, tax_level = "Genus")
#20% Plesiomonas, 22% Clostridium sensu stricto 1, 10% Paeniclostridium


## alpha diversity ----
# richness plots
plot_richness(webber_bats, measures = c("Observed","Chao1", "Shannon"),
              color = "site", shape = "species",
              title = "Alpha diversity") +
  geom_point(size = 2) + theme_classic()

# plot alpha diversity over the summer
bat_metadata |> 
  ggplot(aes(x = date, y = Chao1, color = site))+
  geom_point(size = 1) +
  ggtitle("Observed OTUs over time (BBB and LBB)")

# plot alpha diversity by age
bat_metadata |> 
  ggplot(aes(x = age, y = Chao1, color = species))+
  geom_point(size = 1) +
  ggtitle("Observed OTUs per age category (BBB and LBB)")

# plot alpha diversity by colony size
bat_metadata |> 
  ggplot(aes(x = colony_size, y = Observed)) +
  geom_boxplot() +
  #geom_point(size = 1) +
  ggtitle("Observed ASVs per colony size (BBB and LBB)") +
  xlab("Colony size") + ylab("Observed ASVs")

kruskal.test(Observed ~ site, data = bat_metadata) #p = 0.02

## beta diversity ----
# _ bray NMDS ----
bats_nmds = ordinate(webber_bats, method = "NMDS", distance = "bray")
plot_ordination(webber_bats, bats_nmds, color = "roost_type", shape = "site", title = "NMDS (Bray-Curtis)") + theme(aspect.ratio = 1) + geom_point(size = 1)

## _ nmds stats ----
# Running NMDS in vegan (metaMDS) on distance matrix
# 20 iterations = 1 min
bats_NMS_k2 =
  metaMDS(df_otus_bats,
          distance = "robust.aitchison",
          k = 2,
          maxit = 20, 
          trymax = 20,
          wascores = TRUE)
# stress values: 0.05 - 0.1 is good (can be confident in inferences from plot)

# plots
stressplot(bats_NMS_k2)
stressplot(bats_NMS_k3)
plot(bats_NMS_k2)

# _ robust Aitchison ----
# setup
df_otus_bats = otu.matrix(webber_bats)
df_metadata_bats = as.data.frame(as.matrix(sample_data(webber_bats)))

# generate Aitchison distance matrix
aitchison_bats = metaMDS(df_otus_bats, distance = "robust.aitchison")
autoplot(aitchison_bats)

# plot NMDS of Aitchison distance matrix
autoplot(aitchison_bats, layers = "sites", arrows = FALSE) #sites = samples; species = ASVs

# join to metadata
df_metadata_bats = df_metadata_bats |> rownames_to_column()
aitch_fort_bats <- fortify(aitchison_bats) |> 
  filter(score == "sites") |>
  rename(rowname = label) |> 
  left_join(df_metadata_bats)

# plot Aitchison distances using NMDS
ggplot(aitch_fort_bats, aes(x = NMDS1, y = NMDS2, colour = species)) +
  geom_point()

gg_ordiplot(aitchison_bats, groups = df_metadata_bats$species, hull = TRUE, label = TRUE, spiders = TRUE, ellipse = FALSE, pt.size = 3, plot = TRUE)

## _ aitchison and RDA ----
aitchison_bats = vegdist(feature_table, method = "robust.aitchison")
#I don't think I need to add a pseudocount of 0 (robust part of the clr can include zeroes)
adonis_bats = adonis2(aitchison_bats ~ (sample_data(webber_bats)$site), by = "terms") #terms tested sequentially
# ERROR:  'qr' and 'y' must have the same number of rows

## _ db-RDA ----
# setup
df_otus_bats = otu.matrix(webber_bats)
df_metadata_bats = as.data.frame(as.matrix(sample_data(webber_bats)))
df_metadata_bats = df_metadata_bats |> rownames_to_column()

# __ species ----
rda_bats_species = capscale(formula = df_otus_bats ~ species, data = df_metadata_bats,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_bats_species) #p = 0.001
RsquareAdj(rda_bats_species) #R^2 = 0.0597

plot(rda_bats_species) #blue x: centre of each group; red: ASVs; circles: samples
autoplot(rda_bats_species, layers = "sites", arrows = FALSE, groups = "species") #sites = samples; species = ASVs

# join rda output to metadata
rda_scores_species_df = rda_bats_species |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(bat_metadata, by = "study_id")

# plot
ggplot(rda_scores_species_df, aes(x = CAP1, y = MDS1, colour = species)) +
  geom_point() + ggtitle("RDA on species")

# __ site ----
rda_bats_site = capscale(formula = df_otus_bats ~ site, data = df_metadata_bats,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_bats_site) #p = 0.001
RsquareAdj(rda_bats_site) #R^2 = 0.0945
plot(rda_bats_site) #blue x: centre of each group; red: ASVs; circles: samples
autoplot(rda_bats_site, layers = "sites", arrows = FALSE, groups = "site") #sites = samples; species = ASVs

# join rda output to metadata
rda_scores_sites_df = rda_bats_site |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(bat_metadata, by = "study_id")
# plot the rda per site
ggplot(rda_scores_sites_df, aes(x = CAP1, y = CAP2, colour = site)) +
  geom_point()

# __ model selection ----
# run null and full models
mod0_bats = capscale(df_otus_bats ~ 1, data = df_metadata_bats, distance = "robust.aitchison", na.action = na.exclude)
mod1_bats = capscale(formula = df_otus_bats ~ species + site, data = df_metadata_bats, distance = "robust.aitchison", na.action = na.exclude)

anova(mod1_bats) #p = 0.001
RsquareAdj(mod1_bats) #ajd R^2 = 0.094

# ordistep
step_r2_bats = ordiR2step(mod0_bats, scope = formula(mod1_bats), perm.max = 200, na.action = na.exclude)
step_bats = ordistep(mod0_bats, scope = formula(mod1_bats), perm.max = 200, na.action = na.exclude)

# general distance matrix plot ----
euclidean_pcoa <- ecodist::pco(aitchison_bats)

# create a data frame from principal coordinates
euclidean_pcoa_df <- data.frame(pcoa1 = euclidean_pcoa$vectors[,1], 
                                pcoa2 = euclidean_pcoa$vectors[,2])

# create the plot
euclidean_plot <- ggplot(data = euclidean_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
  geom_point() +
  labs(x = "PC1",
       y = "PC2",
       title = "Euclidean PCoA with CLR transformation") +
  theme(title = element_text(size = 12)) # makes titles smaller
euclidean_plot

## differential abundance ----

# for clr: reads = df where rows are OTUs and cols are samples; conds: vector

## aldex by species
v_species = bat_metadata$species
species_clr = aldex.clr(otu_table(webber_bats), v_species, mc.samples = 200)
species_tt = aldex.ttest(species_clr)
aldex_species_effect = aldex.effect(species_clr, CI = TRUE)
species_aldex_all = data.frame(species_tt, aldex_species_effect)
par(mfrow = c(1,3))
aldex.plot(species_aldex_all, type = "MA", test = "welch", main = "MA plot")
aldex.plot(species_aldex_all, type = "MW", test = "welch", main = "effect plot")
aldex.plot(species_aldex_all, type = "volcano", test = "welch", main = "volcano plot") #Error in match.arg(type) : 'arg' should be one of “MW”, “MA”
summary(species_aldex_all)

# make df of differentially abundant OTUs
diff_otus = species_aldex_all |> 
  filter(we.eBH < 0.05) |> 
  #filter(overlap < 0.01) |> 
  rownames_to_column("...1") |> 
  # add OTU taxonomical info
  left_join(taxonomy_table, by = "...1") |>
  replace_na(list(Genus = "")) |>
  replace_na(list(Species = "sp.")) |>
  #concatenate to keep full OTU name
  mutate(otu_scientific_name = paste(Genus, Species, sep = " ")) |> 
  mutate(otu_scientific_name = str_trim(otu_scientific_name, side = "left")) |>
  mutate(otu_scientific_name = as.factor(otu_scientific_name)) |> 
  mutate(otu_scientific_name = fct_recode(otu_scientific_name, "NA" = "sp.")) |> 
  mutate(ci_product = effect.low * effect.high ) |> #column that returns TRUE if CI does not contain 0# |> 
  mutate(ci_no_zero = ci_product > 0) #returns TRUE if CI does not include zero

# group by OTU identity
diff_otus |> 
  group_by(otu_scientific_name)

write_csv(diff_otus, "diff_otus.csv")

# create df of differentially abundant OTUs from CI
diff_otus_ci = diff_otus |> 
  filter(ci_no_zero == TRUE)

# plot y = difference, x = OTU name
diff_otus |> 
  ggplot(aes(x = effect, y = reorder(otu_scientific_name, (effect*-1)), colour = Phylum)) +
  geom_point() +
  xlab("Effect size") +
  ylab("ASV") +
  theme(axis.text.y = element_text(face = "italic")) +
  ggtitle("Differentially abundant ASVs \n(BH-corrected)") + 
  geom_vline(xintercept = 0, linetype = "dotted")

diff_otus_ci |> 
  ggplot(aes(x = effect, y = otu_scientific_name, colour = Phylum)) +
  geom_point() +
  xlab("Effect size") +
  ylab("ASV") +
  theme(axis.text.y = element_text(face = "italic")) +
  ggtitle("ASVs whose CI does not include zero")


## BBB  ----
# plot alpha diversity by colony size
bbb_metadata |> 
  ggplot(aes(x = colony_size, y = Observed)) +
  geom_boxplot() +
  #geom_point(size = 1) +
  ggtitle("Observed ASVs per colony size (BBB)") +
  xlab("Colony size") + ylab("Observed ASVs")

kruskal.test(Observed ~ site, data = bat_metadata) #p = 0.02


## _ metaMDS ----
nmds_bbb = metaMDS(df_otus_bbb, distance = "robust.aitchison", autotransform = FALSE)

# add metadata
bbb_metadata$mcmaster_sample_id = as.character(bbb_metadata_sib$mcmaster_sample_id)

nmds_fort_bbb = nmds_bbb |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(bbb_metadata, by = "study_id")

# plot
ggplot(nmds_fort_bbb, aes(x = NMDS1, y = NMDS2, colour = site)) +
  geom_point() + ggtitle("NMDS on robust Aitchison distances (BBB)")

# spiderplot
gg_ordiplot(nmds_fort_bbb, groups = bbb_metadata$site, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE)


# _ descriptive stats ----
# heatmap #1
plot_heatmap(bbb, method = "NMDS", distance = "bray")

# heatmap #2: abundance OTUs (representing at least 20% of reads in at least 1 sample)
total = median(sample_sums(webber_bats))
bbb_abundant = filter_taxa(bbb, function(x) sum(x > total*0.20) > 0, TRUE)
plot_heatmap(bbb_abundant, method = "NMDS", distance = "bray") + theme(axis.text.y = element_blank()) +
  ggtitle("Big Brown Bat, Heatmap of abundant OTUs (≥ 20% of reads in ≥ 1 sample")

# relative abundances by colony size
plot_bar(bbb, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") + facet_wrap(~ colony_size) + theme(legend.position = "none")

# relative abundances (normalized) by age
plot_bar(bbb_norm, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") + facet_wrap(~ age)

# by site
plot_bar(bbb_norm, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  facet_grid(~ site) + theme(legend.position = "none")

# by roost type
plot_bar(bbb_norm, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  facet_grid(~ roost_type) + theme(legend.position = "none")

## _ dominant taxa ----
# top 3 most abundant phyla
fantaxtic::top_taxa(bbb, n = 3, tax_level = "Phylum")

# top 5 most abundant orders
fantaxtic::top_taxa(bbb, n = 5, tax_level = "Order")


# _ alpha diversity ----
# plot alpha diversity by site
bat_metadata |> 
  filter(species == "EPFU") |> 
  ggplot(aes(x = date, y = Chao1, color = site)) +
  #geom_point(size = 1) +
  geom_boxplot() +
  ggtitle("Chao1 over time (BBB)")

# plot alpha diversity by age
bat_metadata |>
  filter(species == "EPFU") |> 
  ggplot(aes(x = age, y = Chao1, color = site))+
  geom_point(size = 1) +
  ggtitle("Chao1 per age category (BBB)")

# plot alpha div by roost site + size
site_summary = bbb_metadata |> 
  group_by(site, colony_size) |> tally()

bbb_metadata |> 
  ggplot(aes(x = site, y = Observed, fill = colony_size)) +
  geom_boxplot() +
  labs(title = 'Big brown bat alpha diversity', x = 'site', y = 'observed ASVs') +
  scale_x_discrete(limits = c("Langton", "Wainfleet", "Alliston", "Dunnville", "Mount Pleasant")) +
  geom_text(data = site_summary, aes(site, Inf, label = n), vjust = 1)

# __ kruskal ----
# does alpha div differ by site?
kruskal.test(Observed ~ site, data = bbb_metadata)
#p = 0.0186

# does alpha div differ by age?
kruskal.test(Observed ~ age, data = bbb_metadata) 
# not significant, p = 0.22

# differ bw large and small roosts?
kruskal.test(Observed ~ colony_size, data = bbb_metadata)
kruskal.test(Chao1 ~ colony_size, data = bbb_metadata)
kruskal.test(Shannon ~ colony_size, data = bbb_metadata)
# all p < 0.001 but small colonies have higher alpha diversity

# __ glm ----
bbb_metadata = bbb_metadata |> mutate(julian = lubridate::yday(date))

glm_size_time_bbb = lmer(Observed ~ colony_size + julian + (1|site), data = bbb_metadata)
summary(glm_size_time_bbb)
jtools::summ(glm_size_time_bbb) #size p = 0.17; julian p = 0.38

glm_site_bbb = lmer(Observed ~ site + (1|julian), data = bbb_metadata)
summary(glm_site_bbb)
jtools::summ(glm_site_bbb) #there's a p for each site, don't know how to interpret that

# this is the hypothesis-guided test: size affects adiv
glm_size_bbb = lmer(Observed ~ colony_size + (1|site), data = bbb_metadata)
summary(glm_size_bbb)
jtools::summ(glm_size_bbb)
# Est.    S.E.   t val.   d.f.      p
#  53.47   21.43     2.49   1.53   0.17 (Observed)

## _ beta diversity ----
# bray NMDS
bbb_nmds = ordinate(bbb, method = "NMDS", distance = "bray")

# robust aitchison (robust takes care of zeroes so no need to add a pseudocount of 1)
aitchison_bbb = vegdist(feature_table_bbb, method = "robust.aitchison")

# plot by site
plot_ordination(bbb, bbb_nmds, color = "site", title = "NMDS of BBB samples (Bray-Curtis)") + theme(aspect.ratio = 1) + geom_point(size = 2)

# plot by type and site
plot_ordination(bbb, bbb_nmds, color = "roost_type", shape = "site", title = "NMDS of BBB samples (Bray-Curtis)") + theme(aspect.ratio = 1) + geom_point(size = 1)

# plot by roost type and site
plot_ordination(bbb, aitchison_bbb, color = "roost_type", shape = "site", title = "NMDS of BBB samples (Bray-Curtis)") + theme(aspect.ratio = 1) + geom_point(size = 1) ## this doesn't work!


## _ db-RDA ----
# __ setup ----
df_otus_bbb = otu.matrix(bbb)
bbb_adults <- subset_samples(bbb, age == "A")
df_otus_bbb_adults = otu.matrix(bbb_adults)
bbb_metadata_adults = as.data.frame(as.matrix(sample_data(bbb_adults)))
#bbb_metadata = as.data.frame(as.matrix(sample_data(bbb)))

# __ spatial pcnm ----
# compute euclidean distances from xy coordinates, generate pcnm, then run RDA
# make latitude and longitude numeric
bbb_metadata$X = as.numeric(bbb_metadata$X)
bbb_metadata$Y = as.numeric(bbb_metadata$Y)

# compute pcnm
pcnm_dist = bbb_metadata |> select(X, Y) |> dist(method = "euclidean") |> pcnm()
bbb_metadata_pcnm = bbb_metadata |> cbind(pcnm_dist[["vectors"]])

mod_spatial = rda(df_otus_bbb ~ PCNM1 + PCNM2 + PCNM3, data = bbb_metadata_pcnm)
anova(mod_spatial) #p = 0.061
RsquareAdj(mod_spatial) # R^2 = 0.042

# __ site ----
rda_bbb_site = capscale(formula = df_otus_bbb ~ site, data = bbb_metadata,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_bbb_site) #p = 0.002
plot(rda_bbb_site)
RsquareAdj(rda_bbb_site) #R^2 = 0.0365

# join rda output to metadata by extracting scores
rda_scores_site_df = rda_bbb_site |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(bbb_metadata, by = "study_id")
# plot the rda per site
ggplot(rda_scores_site_df, aes(x = CAP1, y = CAP2, colour = site)) +
  geom_point()+ ggtitle("RDA by site (big brown bats)")

## __ site (adults) ----
# exclude juveniles as they may be allogroomed by adults
rda_bbb_site_adults = capscale(formula = df_otus_bbb_adults ~ site, data = bbb_metadata_adults,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_bbb_site_adults) #p = 0.001
RsquareAdj(rda_bbb_site_adults) #R^2 = 0.0451

# join site_adults output to metadata
rda_scores_site_adults_df = rda_bbb_site_adults |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(bbb_metadata, by = "study_id")

# plot rda site_adults
ggplot(rda_scores_site_adults_df, aes(x = CAP1, y = CAP2, colour = site)) +
  geom_point() + ggtitle("RDA on sites, only adults")

## __ age ----
rda_bbb_age = capscale(formula = df_otus_bbb ~ age, data = bbb_metadata,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_bbb_age, strata = bbb_metadata$site) #p = 0.179

# join age output to metadata
rda_scores_age_df = rda_bbb_age |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(bbb_metadata, by = "study_id")

# plot by age
ggplot(rda_scores_age_df, aes(x = CAP1, y = MDS1, colour = age)) +
  geom_point()

# __ colony size ----
rda_bbb_size = capscale(formula = df_otus_bbb ~ colony_size, data = bbb_metadata,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_bbb_size, strata = bbb_metadata$site) #p = 1

# __ julian ----
# convert dates to julian day
bbb_metadata = bbb_metadata |> mutate(julian = lubridate::yday(date))

# linear time
rda_bbb_julian = capscale(formula = df_otus_bbb ~ julian, data = bbb_metadata,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_bbb_julian) #p = 0.001
RsquareAdj(rda_bbb_julian) #R^2 = 0.0134

# join rda output to metadata by extracting scores
rda_scores_julian_df = rda_bbb_julian |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(bbb_metadata, by = "study_id")

# plot the rda per site
ggplot(rda_scores_julian_df, aes(x = CAP1, y = MDS1, colour = julian, shape = site)) +
  geom_point() + ggtitle("RDA by Julian day (big brown bats)")

# __ time pcnm
pcnm_time = bbb_metadata |> select(julian) |> dist(method = "euclidean") |> pcnm()
bbb_metadata_time = bbb_metadata |> cbind(pcnm_time[["vectors"]])

# run rda + test significance
mod_time = rda(df_otus_bbb ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6  + PCNM7 + PCNM8, data = bbb_metadata_time)
anova(mod_time) #p = 0.037
RsquareAdj(mod_time) #R^2 = 0.0977

# __ time, aw ----
bbb_wa <- subset_samples(bbb, site == "Alliston" | site == "Wainfleet")
df_otus_bbb_wa = otu.matrix(bbb_wa)
bbb_metadata_wa = as.data.frame(as.matrix(sample_data(bbb_wa)))
bbb_metadata_wa = bbb_metadata_wa |> mutate(julian = lubridate::yday(date))

# linear time
rda_bbb_julian_wa = capscale(formula = df_otus_bbb_wa ~ julian, data = bbb_metadata_wa,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_bbb_julian_wa) #p = 0.002
RsquareAdj(rda_bbb_julian_wa) #R^2 = 0.0240

# join rda output to metadata by extracting scores
bbb_metadata_wa = bbb_metadata_wa |> rownames_to_column()
rda_scores_julian_wa_df = rda_bbb_julian_wa |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> 
  left_join(bbb_metadata_wa, by = "rowname")
# plot the rda per site
ggplot(rda_scores_julian_wa_df, aes(x = CAP1, y = MDS1, colour = julian, shape = site)) +
  geom_point() + ggtitle("RDA by Julian day (big brown bats)")

# pcnm
pcnm_time_wa = bbb_metadata_wa |> select(julian) |> dist(method = "euclidean") |> pcnm()
bbb_metadata_time_wa = bbb_metadata_wa |> cbind(pcnm_time_wa[["vectors"]])

# run rda + test significance
mod_time_wa = rda(df_otus_bbb_wa ~ PCNM1 + PCNM2 + PCNM3 + PCNM4, data = bbb_metadata_time_wa)
anova(mod_time_wa) #p = 0.01
RsquareAdj(mod_time_wa) #R^2 = 0.130


# __ site + time ----
rda_bbb_site_julian = capscale(formula = df_otus_bbb ~ site + julian, data = bbb_metadata,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_bbb_site_julian) #p = 0.001
RsquareAdj(rda_bbb_site_julian) #R^2 = 0.0433

# join rda output to metadata by extracting scores
rda_scores_site_julian_df = rda_bbb_site_julian |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(bbb_metadata, by = "study_id")
# plot the rda per site
ggplot(rda_scores_site_julian_df, aes(x = CAP1, y = CAP2, colour = site)) +
  geom_point()+ ggtitle("RDA by site (big brown bats)")


# __ model selection ----
# run null and full models
mod0_bbb = capscale(df_otus_bbb ~ 1, data = bbb_metadata, distance = "robust.aitchison", na.action = na.exclude)

mod1_bbb = capscale(formula = df_otus_bbb ~ site + julian, data = bbb_metadata, distance = "robust.aitchison", na.action = na.exclude)
anova(mod1_bbb) #p = 0.001
RsquareAdj(mod1_bbb) #ajd R^2 = 0.0433

# ordistep
step_r2_bbb = ordiR2step(mod0_bbb, scope = formula(mod1_bbb), perm.max = 200, na.action = na.exclude)



## core microbiome ----
# packages required
library(eulerr)
library(microbiome)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)

# core microbiome for each species + site (and how much do they overlap)
# convert to relative abundances
pseq.rel <- microbiome::transform(webber_bats, "compositional")
core.taxa.standard <- core_members(pseq.rel, detection = 0.0001, prevalence = 50/100)

# create phyloseq object of core microbiome
pseq.core <- core(pseq.rel, detection = 0.0001, prevalence = .5)
# retrieve taxa names
core.taxa <- taxa(pseq.core)

# get the taxonomy data
tax.mat <- tax_table(pseq.core)
tax.df <- as.data.frame(tax.mat)

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only those ASVs that are core based on thresholds
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
knitr::kable(head(core.taxa.class))

# 2 sets of detection levels
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)
prevalences <- seq(.05, 1, .05)

# plot
plot_core(pseq.rel.sp.bbb, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()


# Also define gray color palette
gray <- gray(seq(0,1,length=5))
p1 <- plot_core(pseq.rel, 
                plot.type = "heatmap", 
                prevalences = prevalences, 
                detections = detections, min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")
p1 <- p1 + theme_bw() + ylab("ASVs") + scale_fill_viridis()
p1

# need to reprocess the figure to include taxonomic information
library(RColorBrewer)
library(knitr)
# get the data used for plotting 
df <- p1$data 

# get the list of OTUs
list <- df$Taxa

# get the taxonomy data
tax <- as.data.frame(tax_table(pseq.rel))

# add the ASVs to last column
tax$ASV <- rownames(tax)

# select taxonomy of only those OTUs that are used in the plot
tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 

# merge all the columns into one
tax.unit <- tidyr::unite(tax2, Taxa_level,c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"), sep = "_;", remove = TRUE)

tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)

# add this new information into the plot data df
df$Taxa <- tax.unit$Taxa_level

# you can see now we have the taxonomic information
knitr::kable(head(df))

p1$data <- df
plot(p1 + theme(axis.text.y = element_text(face="italic")))

# plot at genus level
pseq.rel.gen <- aggregate_taxa(pseq.rel, "Genus")

prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)

p_genus <- plot_core(pseq.rel.gen, 
                plot.type = "heatmap", 
                prevalences = prevalences, 
                detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")
p_genus <- p_genus + theme_minimal() + ylab("Genus") + scale_fill_viridis()
p_genus

## _ bbb core genera ----
pseq.rel.bbb <- subset_samples(pseq.rel, species == "EPFU") 

# aggregate taxa
pseq.rel.gen.bbb <- aggregate_taxa(pseq.rel.bbb, "Genus")
pseq.rel.sp.bbb <- aggregate_taxa(pseq.rel.bbb, "Species")

# remove unknown
#pseq.rel.gen.bbb <- subset_taxa(pseq.rel.gen.bbb, Genus != "Unknown")

# plot
p_genus_bbb <- plot_core(pseq.rel.gen.bbb, 
                     plot.type = "heatmap", 
                     prevalences = prevalences, 
                     detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")
p_genus_bbb <- p_genus_bbb + theme_minimal() + ylab("Genus") + scale_fill_viridis() + ggtitle("Core microbiome (big brown bats)")
p_genus_bbb

# _ bbb core species ----
# aggregate taxa to species
pseq.rel.sp.bbb <- aggregate_taxa(pseq.rel.bbb, "Species") #aggregate_taxa merges all unclassified to "Unknown"

# check if any taxa unassigned at species
any(taxa_names(pseq.rel.sp.bbb) == "Unknown")

# remove unknown species; this makes it non-compositional!
pseq.rel.sp.bbb <- subset_taxa(pseq.rel.sp.bbb, Species != "Unknown")

# remake compositional
pseq.rel.sp.bbb.cp = microbiome::transform(pseq.rel.sp.bbb, "compositional")

# get genus names from taxonomy table
taxonomy_table |> 
  filter(Species == "faecium") |> pull(Genus) #Enterococcus
taxonomy_table |> 
  filter(Species == "alvei") |> pull(Genus) #Hafnia-Obesumbacterium
taxonomy_table |> 
  filter(Species == "lactis") |> pull(Genus) #Lactococcus
taxonomy_table |> 
  filter(Species == "fonticola") |> pull(Genus) #Serratia
taxonomy_table |> 
  filter(Species == "garvieae") |> pull(Genus) #Lactococcus
taxonomy_table |> 
  filter(Species == "gilvus") |> pull(Genus) #Enterococcus
taxonomy_table |> 
  filter(Species == "mundtii") |> pull(Genus) #Enterococcus

# set prevalence & detection levels
prevalences <- seq(.05, 1, .0)
det <- c(0, 0.1, 0.5, 1, 2)/100
det2 <- c(0, 0.01, 1, 0.5, 2)/100
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)

# plot core species
p_sp_bbb <- plot_core(pseq.rel.sp.bbb.cp,
                      plot.type = "heatmap",
                      prevalences = prevalences,
                      detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

p_sp_bbb <- p_sp_bbb + theme_minimal() + ylab("ASVs\n") + scale_fill_viridis() + ggtitle("Core microbiome (big brown bats)") + scale_y_discrete(labels=c("Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus_faecium" = "Enterococcus faecium", "alvei" = "Hafnia-Obesumbacterium alvei", "lactis" = "Lactococcus lactis", "fonticola" = "Serratia fonticola", "garvieae" = "Lactococcus garvieae", "gilvus" = "Enterococcus gilvus", "mundtii" = "Enterococcus mundtii")) +
  theme(axis.text.y = element_text(face = "italic"))

p_sp_bbb

# compare to 7 most dominant species
fantaxtic::top_taxa(bbb, n = 7, tax_level = "Species")


# _ lbb core sp ----
pseq.rel.lbb <- subset_samples(pseq.rel, species == "MYLU") 
# aggregate taxa to species
pseq.rel.sp.lbb <- aggregate_taxa(pseq.rel.lbb, "Species")

# Check if any taxa with no species classification. aggregate_taxa will merge all unclassified to Unknown
any(taxa_names(pseq.rel.sp.lbb) == "Unknown")

# remove unknown species; this makes it non-compositional!
pseq.rel.sp.lbb <- subset_taxa(pseq.rel.sp.lbb, Species != "Unknown") #TRUE

# remake compositional
pseq.rel.sp.lbb.cp = microbiome::transform(pseq.rel.sp.lbb, "compositional")

# set prevalence & detection levels
prevalences <- seq(.05, 1, .0)
det <- c(0, 0.1, 0.5, 1, 2)/100
det2 <- c(0, 0.01, 1, 0.5, 2)/100
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)

# plot core species
p_sp_lbb <- plot_core(pseq.rel.sp.lbb.cp,
                      plot.type = "heatmap",
                      prevalences = prevalences,
                      detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

# get genus names from taxonomy table
taxonomy_table |> 
  filter(Species == "perfringens") |> pull(Genus) #Clostridium sensu stricto 1 perfringens
taxonomy_table |> 
  filter(Species == "sordellii") |> pull(Genus) #Paeniclostridium sordellii
taxonomy_table |> 
  filter(Species == "colicanis") |> pull(Genus) #Clostridium sensu stricto 1 colicanis
taxonomy_table |> 
  filter(Species == "shigelloides") |> pull(Genus) #Plesiomonas shigelloides
taxonomy_table |> 
  filter(Species == "saprophyticus") |> pull(Genus) #Staphylococcus saprophyticus
taxonomy_table |> 
  filter(Species == "mundtii") |> pull(Genus) #Enterococcus mundtii
taxonomy_table |> 
  filter(Species == "urealyticum") |> pull(Genus) #Corynebacterium urealyticum
taxonomy_table |> 
  filter(Species == "vinsonii") |> pull(Genus)
#Bartonella vinsonii

p_sp_lbb <- p_sp_lbb + theme_minimal() + ylab("ASVs\n") + scale_fill_viridis() + ggtitle("Core microbiome (little brown bats)") + scale_y_discrete(labels=c("perfringens" = "Clostridium sensu stricto 1 perfringens", "Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus_faecium" = "Enterococcus faecium", "sordellii" = "Paeniclostridium sordellii", "colicanis" = "Clostridium sensu stricto 1 colicanis", "shigelloides" = "Plesiomonas shigelloides", "saprophyticus" = "Staphylococcus saprophyticus", "mundtii" = "Enterococcus mundtii", "urealyticum" = "Corynebacterium urealyticum", "Bacteria_Firmicutes_Bacilli_Staphylococcales_Staphylococcaceae_Staphylococcus_lentus" = "Staphylococcus lentus", "vinsonii" = "Bartonella vinsonii")) +
  theme(axis.text.y = element_text(face = "italic"))

p_sp_lbb

# compare to 7 most dominant species
fantaxtic::top_taxa(bbb, n = 7, tax_level = "Species")


# _ bbb sp core site ----
pseq.rel.alliston <- subset_samples(pseq.rel, site == "Alliston") 
pseq.rel.dunnville <- subset_samples(pseq.rel, site == "Dunnville") 
pseq.rel.langton <- subset_samples(pseq.rel, site == "Langton") 
pseq.rel.mtp <- subset_samples(pseq.rel, site == "Mount Pleasant") 
pseq.rel.wainfleet <- subset_samples(pseq.rel, site == "Wainfleet") 

# __ alliston ----
# aggregate taxa to species
pseq.rel.sp.alliston <- aggregate_taxa(pseq.rel.alliston, "Species")

# check if any taxa unassigned at species
any(taxa_names(pseq.rel.sp.alliston) == "Unknown")

# remove unknown species; this makes it non-compositional!
pseq.rel.sp.alliston <- subset_taxa(pseq.rel.sp.alliston, Species != "Unknown")

# remake compositional
pseq.rel.sp.alliston.cp = microbiome::transform(pseq.rel.sp.alliston, "compositional")

# plot core species
p_sp_alliston <- plot_core(pseq.rel.sp.alliston.cp,
                           plot.type = "heatmap",
                           prevalences = prevalences,
                           detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

p_sp_alliston <- p_sp_alliston + theme_minimal() + ylab("ASVs") + scale_fill_viridis() + ggtitle("Core microbiome \n(Alliston) (n = 17)") + scale_y_discrete(labels=c("Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus_faecium" = "Enterococcus faecium", "alvei" = "Hafnia-Obesumbacterium alvei", "lactis" = "Lactococcus lactis", "fonticola" = "Serratia fonticola", "gilvus" = "Enterococcus gilvus", "proteus" = "Hafnia-Obesumbacterium proteus", "asini" = "Enterococcus proteus asini", "freundii" = "Citrobacter freundii", "pestis" = "Yersinia sp.", "alginatilytica" = "Dysgonomonas alginatilytica")) +
  theme(axis.text.y = element_text(face = "italic"))

p_sp_alliston

# get genus names from taxonomy table
taxonomy_table |> 
  filter(Species == "proteus") |> pull(Genus) #Hafnia-Obesumbacterium
taxonomy_table |> 
  filter(Species == "asini") |> pull(Genus) #Enterococcus
taxonomy_table |> 
  filter(Species == "freundii") |> pull(Genus) #Citrobacter
taxonomy_table |> 
  filter(Species == "pestis") |> pull(Genus) #Yersinia
taxonomy_table |> 
  filter(Species == "alginatilytica") |> pull(Genus) #Dysgonomonas

# __ dunnville ----
# aggregate taxa to species
pseq.rel.sp.dunnville <- aggregate_taxa(pseq.rel.dunnville, "Species")

# check if any taxa unassigned at species
any(taxa_names(pseq.rel.sp.dunnville) == "Unknown")

# remove unknown species
pseq.rel.sp.dunnville <- subset_taxa(pseq.rel.sp.dunnville, Species != "Unknown")

# remake compositional
pseq.rel.sp.dunnville.cp = microbiome::transform(pseq.rel.sp.dunnville, "compositional")

# plot core species
p_sp_dunnville <- plot_core(pseq.rel.sp.dunnville.cp,
                            plot.type = "heatmap",
                            prevalences = prevalences,
                            detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

p_sp_dunnville <- p_sp_dunnville + theme_minimal() + ylab("ASVs") + scale_fill_viridis() + ggtitle("Core microbiome \n(Dunville) (n = 8)") + scale_y_discrete(labels=c("Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus_faecium" = "Enterococcus faecium", "alvei" = "Hafnia-Obesumbacterium alvei", "lactis" = "Lactococcus lactis", "fonticola" = "Serratia fonticola", "gilvus" = "Enterococcus gilvus", "proteus" = "Hafnia-Obesumbacterium proteus", "asini" = "Enterococcus proteus asini", "garvieae" = "Lactococcus garvieae", "Bacteria_Firmicutes_Bacilli_Entomoplasmatales_Spiroplasmataceae_Spiroplasma_citri" = "Spiroplasma citri", "cuneatus" = "Desulfovibrio cuneatus", "pittii" = "Acinetobacter pittii")) +
  theme(axis.text.y = element_text(face = "italic"))

p_sp_dunnville

taxonomy_table |> 
  filter(Species == "citri") |> pull(Genus) #Spiroplasma
taxonomy_table |> 
  filter(Species == "cuneatus") |> pull(Genus) #Desulfovibrio cuneatus
taxonomy_table |> 
  filter(Species == "pittii") |> pull(Genus) #Acinetobacter pittii


# __ langton ----
# aggregate taxa to species
pseq.rel.sp.langton <- aggregate_taxa(pseq.rel.langton, "Species")

# check if any taxa unassigned at species
any(taxa_names(pseq.rel.sp.langton) == "Unknown")

# remove unknown species
pseq.rel.sp.langton <- subset_taxa(pseq.rel.sp.langton, Species != "Unknown")

# remake compositional
pseq.rel.sp.langton.cp = microbiome::transform(pseq.rel.sp.langton, "compositional")

# plot core species
p_sp_langton <- plot_core(pseq.rel.sp.langton.cp,
                          plot.type = "heatmap",
                          prevalences = prevalences,
                          detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

p_sp_langton <- p_sp_langton + theme_minimal() + ylab("ASVs") + scale_fill_viridis() + ggtitle("Core microbiome \n(Langton) (n = 2)") + scale_y_discrete(labels=c("Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus_faecium" = "Enterococcus faecium", "alvei" = "Hafnia-Obesumbacterium alvei", "lactis" = "Lactococcus lactis", "fonticola" = "Serratia fonticola", "gilvus" = "Enterococcus gilvus", "proteus" = "Hafnia-Obesumbacterium proteus", "asini" = "Enterococcus proteus asini", "garvieae" = "Lactococcus garvieae", "Bacteria_Firmicutes_Bacilli_Entomoplasmatales_Spiroplasmataceae_Spiroplasma_citri" = "Spiroplasma citri", "cuneatus" = "Desulfovibrio cuneatus", "pittii" = "Acinetobacter pittii", "freundii" = "Citrobacter freundii", "pestis" = "Yersinia sp.", "alginatilytica" = "Dysgonomonas alginatilytica", "ixodetis" = "Citrobacter ixodetis", "bacterium" = "Lachnospiraceae bacterium")) +
  theme(axis.text.y = element_text(face = "italic"))

p_sp_langton

taxonomy_table |> 
  filter(Species == "ixodetis") |> pull(Genus) #Spiroplasma
taxonomy_table |> 
  filter(Species == "bacterium") |> pull(Genus)  #Lachnospiraceae NK4A136 group

# __ mt pleasant ----
# aggregate taxa to species
pseq.rel.sp.mtp <- aggregate_taxa(pseq.rel.mtp, "Species")

# check if any taxa unassigned at species
any(taxa_names(pseq.rel.sp.mtp) == "Unknown")

# remove unknown species
pseq.rel.sp.mtp <- subset_taxa(pseq.rel.sp.mtp, Species != "Unknown")

# remake compositional
pseq.rel.sp.mtp.cp = microbiome::transform(pseq.rel.sp.mtp, "compositional")

# plot core species
p_sp_mtp <- plot_core(pseq.rel.sp.mtp.cp,
                      plot.type = "heatmap",
                      prevalences = prevalences,
                      detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

p_sp_mtp <- p_sp_mtp + theme_minimal() + ylab("ASVs") + scale_fill_viridis() + ggtitle("Core microbiome \n(Mt Pleasant) (n = 6)") + scale_y_discrete(labels=c("Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus_faecium" = "Enterococcus faecium", "alvei" = "Hafnia-Obesumbacterium alvei", "lactis" = "Lactococcus lactis", "fonticola" = "Serratia fonticola", "gilvus" = "Enterococcus gilvus", "proteus" = "Hafnia-Obesumbacterium proteus", "asini" = "Enterococcus proteus asini", "alginatilytica" = "Dysgonomonas alginatilytica", "freundii" = "Citrobacter freundii", "pestis" = "Yersinia sp.", "salamandrae" = "Candidatus Amphibiichlamydia salamandrae")) +
  theme(axis.text.y = element_text(face = "italic"))

p_sp_mtp

taxonomy_table |> 
  filter(Species == "salamandrae") |> pull(Genus) #Candidatus Amphibiichlamydia

# __ wainfleet ----
# aggregate taxa to species
pseq.rel.sp.wainfleet <- aggregate_taxa(pseq.rel.wainfleet, "Species")

# check if any taxa unassigned at species
any(taxa_names(pseq.rel.sp.wainfleet) == "Unknown")

# remove unknown species
pseq.rel.sp.wainfleet <- subset_taxa(pseq.rel.sp.wainfleet, Species != "Unknown")

# remake compositional
pseq.rel.sp.wainfleet.cp = microbiome::transform(pseq.rel.sp.wainfleet, "compositional")

# plot core species
p_sp_wainfleet <- plot_core(pseq.rel.sp.wainfleet.cp,
                            plot.type = "heatmap",
                            prevalences = prevalences,
                            detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

p_sp_wainfleet <- p_sp_wainfleet + theme_minimal() + ylab("ASVs") + scale_fill_viridis() + ggtitle("Core microbiome \n(Wainfleet) (n = 20)") + scale_y_discrete(labels=c("Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus_faecium" = "Enterococcus faecium", "alvei" = "Hafnia-Obesumbacterium alvei", "lactis" = "Lactococcus lactis", "fonticola" = "Serratia fonticola", "gilvus" = "Enterococcus gilvus", "proteus" = "Hafnia-Obesumbacterium proteus", "asini" = "Enterococcus proteus asini", "alginatilytica" = "Dysgonomonas alginatilytica", "freundii" = "Citrobacter freundii", "pestis" = "Yersinia sp.", "mundtii" = "Enterococcus mundtii", "garvieae" = "Lactococcus garvieae", "morganii" = "Morganella morganii")) +
  theme(axis.text.y = element_text(face = "italic"))

p_sp_wainfleet

taxonomy_table |> 
  filter(Species == "mundtii") |> pull(Genus) #Enterococcus 
taxonomy_table |> 
  filter(Species == "garvieae") |> pull(Genus) #Lactococcus
taxonomy_table |> 
  filter(Species == "morganii") |> pull(Genus) #Morganella

## _ bbb gen core site ----
# __ alliston ----
# aggregate taxa to genus
pseq.rel.gen.alliston <- aggregate_taxa(pseq.rel.alliston, "Genus")

# check if any taxa unassigned at species
any(taxa_names(pseq.rel.gen.alliston) == "Unknown")

# remove unknown
pseq.rel.gen.alliston <- subset_taxa(pseq.rel.gen.alliston, Genus != "Unknown")

# remake compositional
pseq.rel.gen.alliston.cp = microbiome::transform(pseq.rel.gen.alliston, "compositional")

# plot core genera
p_gen_alliston <- plot_core(pseq.rel.gen.alliston.cp,
                            plot.type = "heatmap",
                            prevalences = prevalences,
                            detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

p_gen_alliston <- p_gen_alliston + theme_minimal() + ylab("Genus\n") + scale_fill_viridis() + ggtitle("Alliston (n = 17)") +
  theme(axis.text.y = element_text(face = "italic")) +
  scale_y_discrete(labels = c( "Bacteria_Firmicutes_Clostridia_Oscillospirales_Ruminococcaceae_Incertae Sedis" = "Incertae Sedis"))

p_gen_alliston

# __ dunnville ----
pseq.rel.gen.dunnville <- aggregate_taxa(pseq.rel.dunnville, "Genus")

any(taxa_names(pseq.rel.gen.dunnville) == "Unknown")

pseq.rel.gen.dunnville <- subset_taxa(pseq.rel.gen.dunnville, Genus != "Unknown")

pseq.rel.gen.dunnville.cp = microbiome::transform(pseq.rel.gen.dunnville, "compositional")

# plot core genera
p_gen_dunnville <- plot_core(pseq.rel.gen.dunnville.cp,
                             plot.type = "heatmap",
                             prevalences = prevalences,
                             detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

p_gen_dunnville <- p_gen_dunnville + theme_minimal() + ylab("Genus\n") + scale_fill_viridis() + ggtitle("Dunville (n = 8)") +
  theme(axis.text.y = element_text(face = "italic"))

p_gen_dunnville


# __ langton ----
pseq.rel.gen.langton <- aggregate_taxa(pseq.rel.langton, "Genus")

any(taxa_names(pseq.rel.gen.langton) == "Unknown")

pseq.rel.gen.langton <- subset_taxa(pseq.rel.gen.langton, Genus != "Unknown")

pseq.rel.gen.langton.cp = microbiome::transform(pseq.rel.gen.langton, "compositional")

# plot core genera
p_gen_langton <- plot_core(pseq.rel.gen.langton.cp,
                           plot.type = "heatmap",
                           prevalences = prevalences,
                           detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

p_gen_langton <- p_gen_langton + theme_minimal() + ylab("Genus\n") + scale_fill_viridis() + ggtitle("Langton (n = 2)") +
  theme(axis.text.y = element_text(face = "italic"))

p_gen_langton

# __ mt pleasant ----
pseq.rel.gen.mtp <- aggregate_taxa(pseq.rel.mtp, "Genus")
any(taxa_names(pseq.rel.gen.mtp) == "Unknown")

pseq.rel.gen.mtp <- subset_taxa(pseq.rel.gen.mtp, Genus != "Unknown")

pseq.rel.gen.mtp.cp = microbiome::transform(pseq.rel.gen.mtp, "compositional")

# plot core genera
p_gen_mtp <- plot_core(pseq.rel.gen.mtp.cp,
                       plot.type = "heatmap",
                       prevalences = prevalences,
                       detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

p_gen_mtp <- p_gen_mtp + theme_minimal() + ylab("Genus\n") + scale_fill_viridis() + ggtitle("Mt Pleasant (n = 6)") +
  theme(axis.text.y = element_text(face = "italic"))

p_gen_mtp

# __ wainfleet ----
pseq.rel.gen.wainfleet <- aggregate_taxa(pseq.rel.wainfleet, "Genus")

any(taxa_names(pseq.rel.gen.wainfleet) == "Unknown")

pseq.rel.gen.wainfleet <- subset_taxa(pseq.rel.gen.wainfleet, Genus != "Unknown")

pseq.rel.gen.wainfleet.cp = microbiome::transform(pseq.rel.gen.wainfleet, "compositional")

# plot core genera
p_gen_wainfleet <- plot_core(pseq.rel.gen.wainfleet.cp,
                             plot.type = "heatmap",
                             prevalences = prevalences,
                             detections = det, min.prevalence = .5) +
  xlab("Detection threshold (relative abundance)")

p_gen_wainfleet <- p_gen_wainfleet + theme_minimal() + ylab("Genus\n") + scale_fill_viridis() + ggtitle("Wainfleet (n = 20)") +
  theme(axis.text.y = element_text(face = "italic"))

p_gen_wainfleet


## _ diff ab sites ----

# subset sites
wa <- subset_samples(webber_bats, site == c("Wainfleet", "Alliston"))

wa  <- webber_bats |> 
  subset_samples(site %in% c("Wainfleet", "Alliston"))

# aldex by site
v_sites = bbb_metadata$site
v_wa = bbb_metadata |> filter(site == "Alliston" | site == "Wainfleet") |> pull(site)
site_clr = aldex.clr(otu_table(wa), v_wa, mc.samples = 200)
site_tt = aldex.ttest(site_clr)
aldex_site_effect = aldex.effect(site_clr, CI = TRUE)
site_aldex_all = data.frame(site_tt, aldex_site_effect)

# plots
par(mfrow = c(1,2))
aldex.plot(site_aldex_all, type = "MA", test = "welch", main = "MA plot")
aldex.plot(site_aldex_all, type = "MW", test = "welch", main = "effect plot")

# make df of differentially abundant ASVs
diff_otus_site = site_aldex_all |> 
  filter(we.eBH < 0.05) |> 
  #filter(overlap < 0.01) |> 
  rownames_to_column("...1") |> 
  # add OTU taxonomical info
  left_join(taxonomy_table, by = "...1") |>
  replace_na(list(Genus = "")) |>
  replace_na(list(Species = "sp.")) |>
  #concatenate to keep full OTU name
  mutate(otu_scientific_name = paste(Genus, Species, sep = " ")) |> 
  mutate(otu_scientific_name = str_trim(otu_scientific_name, side = "left")) |>
  mutate(otu_scientific_name = as.factor(otu_scientific_name)) |> 
  #mutate(otu_scientific_name = fct_recode(otu_scientific_name, "NA" = "sp.")) |> 
  mutate(ci_product = effect.low * effect.high ) |> #column that returns TRUE if CI does not contain 0# |> 
  mutate(ci_no_zero = ci_product > 0) #returns TRUE if CI does not include zero

# differentially abundant OTUs based on CI
diff_otus_site |> 
  filter(ci_no_zero == TRUE) #there are none

# plot y = difference, x = OTU name
diff_otus_site |> 
  ggplot(aes(x = effect, y = reorder(otu_scientific_name, (effect*-1)), colour = Phylum)) +
  geom_point() +
  xlab("Effect size") +
  ylab("ASV") +
  theme(axis.text.y = element_text(face = "italic")) +
  ggtitle("Differentially abundant ASVs \nWainfleet and Alliston\n(BH-corrected)") + 
  geom_vline(xintercept = 0, linetype = "dotted")

# by confidence interval
diff_otus_ci |> 
  ggplot(aes(x = effect, y = otu_scientific_name, colour = Phylum)) +
  geom_point() +
  xlab("Effect size") +
  ylab("ASV") +
  theme(axis.text.y = element_text(face = "italic")) +
  ggtitle("ASVs whose CI does not include zero")

