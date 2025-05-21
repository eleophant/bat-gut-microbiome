###
#
# bats figures
#
###

# packages ----
library(tidyverse)
library(viridis)
library(patchwork)
library(leaflet)
library(sf)
library(cowplot)
library(phyloseq)
library(lme4)
library(vegan)
library(ggordiplots)
library(microbiome)
library(ALDEx2)
theme_set(theme_classic())


# data ----
load(file = "bats_analysis_2024_09_05.RData")

# 2. map ----
# generate df with site coordinates & correct Mt Pleasant
sites = bat_metadata |> 
  select(site, X, Y) |> 
  mutate( 
    X = case_when(site == "Mount Pleasant" ~ 43.15694, TRUE ~ X),
    Y = case_when(site == "Mount Pleasant" ~ 80.2575, TRUE ~ Y)) |> 
  mutate(Y = -Y) |> 
  mutate(site = as.factor(site)) |> 
  rename(lat = X, lon = Y)
sites = sites |> mutate(Y = lat, X = lon)

# _ North America ----
# Toronto coordinates
toronto_lat = 43.6532
toronto_lon = -79.3832

leaflet(options = leafletOptions(zoomControl = FALSE))  |>
  addProviderTiles(providers$Esri.WorldImagery)  |>  # Satellite imagery of North America
  addCircleMarkers(lng = toronto_lon, lat = toronto_lat, 
                   color = "white",  # Set the color of the circle marker to white
                   fill = TRUE,      # Fill the circle
                   fillColor = "white",  # Fill the circle with white
                   fillOpacity = 1,  # Set full opacity for the circle
                   radius = 6,       # Set the radius of the circle
                   stroke = FALSE) |>    # Remove the border of the circle)
  setView(lng = toronto_lon, lat = toronto_lat, zoom = 2.4)  # Set zoom level for North America


# _ Ontario ----
leaflet(data = sites, options = leafletOptions(zoomControl = FALSE))  |>
  addProviderTiles(providers$Esri.WorldImagery)  |>
  
  # add points for roost sites
  addCircleMarkers(~lon, ~lat, label = ~site, color = "white", radius = 5) |> #, 
  #labelOptions = labelOptions(noHide = TRUE, 
  #direction = ~ifelse(site == "Dunnville", "bottom", "top"),
  #textOnly = TRUE, 
  #style = list("color" = "white", "font-size" = "16px")))  |>
  
  # Add landmark labels
  addLabelOnlyMarkers(lng = -79.3832, lat = 43.6532, label = "Toronto", 
                      labelOptions = labelOptions(noHide = TRUE, direction = "center", 
                                                  textOnly = TRUE, style = list("color" = "white", "font-size" = "16px")))  |> 
  addLabelOnlyMarkers(lng = -78, lat = 43.7, label = "Lake Ontario", 
                      labelOptions = labelOptions(noHide = TRUE, direction = "center", 
                                                  textOnly = TRUE, style = list("color" = "white", "font-size" = "16px")))  |>
  addLabelOnlyMarkers(lng = -81.2, lat = 42.15, label = "Lake Erie", 
                      labelOptions = labelOptions(noHide = TRUE, direction = "center", 
                                                  textOnly = TRUE, style = list("color" = "white", "font-size" = "16px")))  |>
  addLabelOnlyMarkers(lng = -82.2, lat = 44.5, label = "Lake Huron", 
                      labelOptions = labelOptions(noHide = TRUE, direction = "center", 
                                                  textOnly = TRUE, style = list("color" = "white", "font-size" = "16px")))  |>
  
  # add scale bar
  addScaleBar(position = "bottomright", 
              options = scaleBarOptions(imperial = FALSE, metric = TRUE, maxWidth = 200))  |>
  
  # adjust color of the scale bar
  htmlwidgets::onRender('
    function(el, x) {
      var scaleBar = document.querySelector(".leaflet-control-scale-line");
      scaleBar.style.backgroundColor = "white";
      scaleBar.style.borderColor = "white";
    }
  ')

# 3. microbiomes ----
## A. rel ab bbb ----
rel_ab_bbb = plot_bar(bbb_norm, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "E. fuscus", y = "Normalized abundance") +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme(axis.title.x = element_text(face = "italic")) +
  ggtitle("A")

## B. rel ab lbb ----
rel_ab_lbb = plot_bar(lbb_norm, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") + 
  #theme(legend.position = "none") +
  #theme(axis.title.y = element_blank()) +
  labs(x = "M. lucifugus", y = "Normalized abundance") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme(axis.title.x = element_text(face = "italic")) +
  ggtitle("B")

## C. rda species ----
rda_bats_species = capscale(formula = df_otus_bats ~ species, data = df_metadata_bats,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_bats_species) #F = 5.13, p = 0.001
RsquareAdj(rda_bats_species) #adj R^2 = 0.0597

# join rda output to metadata
rda_scores_species_df = rda_bats_species |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(bat_metadata, by = "study_id")
rda_scores_species_df = rda_scores_species_df |> 
  mutate(species = recode(species,
                          EPFU = "E. fuscus",
                          MYLU = "M. lucifugus"))

# plot
p_rda_species = ggplot(rda_scores_species_df, aes(x = CAP1, y = MDS1, colour = species)) +
  geom_point() +
  theme(legend.text = element_text(face = "italic", size = 11)) +
  ggtitle("C") +
  scale_colour_viridis_d(option = "viridis")

## D. rda site ----
rda_bbb_site = capscale(formula = df_otus_bbb ~ site, data = bbb_metadata,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_bbb_site) #F = 1.49, p = 0.002
plot(rda_bbb_site)
RsquareAdj(rda_bbb_site) #adj R^2 = 0.0365

# join rda output to metadata by extracting scores
rda_scores_site_df = rda_bbb_site |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(bbb_metadata, by = "study_id")

p_rda_site = ggplot(rda_scores_site_df, aes(x = CAP1, y = CAP2, colour = site)) +
  geom_point() + 
  theme(legend.text = element_text(size = 11)) +
  ggtitle("D") +
  scale_colour_viridis_d(option = "plasma")

## E. all ----
rel_ab_bbb + rel_ab_lbb + plot_layout(widths = c(4, 1))
p_rda_species + p_rda_site

# 4. diff ab ----
# aldex by species
v_species = bat_metadata$species
species_clr = aldex.clr(otu_table(webber_bats), v_species, mc.samples = 200)
species_tt = aldex.ttest(species_clr)
aldex_species_effect = aldex.effect(species_clr, CI = TRUE)
species_aldex_all = data.frame(species_tt, aldex_species_effect)
par(mfrow = c(1,2))
aldex.plot(species_aldex_all, type = "MA", test = "welch", main = "MA plot")
aldex.plot(species_aldex_all, type = "MW", test = "welch", main = "effect plot")

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

# create df of differentially abundant OTUs from CI
diff_otus_ci = diff_otus |> 
  filter(ci_no_zero == TRUE)

# plot y = difference, x = OTU name
diff_otus |> 
  filter(otu_scientific_name != "NA") |> 
  ggplot(aes(x = effect, y = reorder(otu_scientific_name, (effect*-1)), colour = Phylum)) +
  geom_point() +
  xlab("Effect size") +
  ylab("ASV") +
  theme(axis.text.y = element_text(face = "italic")) +
  #ggtitle("Differentially abundant ASVs \n(BH-corrected)") + 
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d()


# 5. adiv ----
site_summary = bbb_metadata |> 
  group_by(site, colony_size) |> tally()

bbb_metadata |> 
  ggplot(aes(x = site, y = Observed, fill = colony_size)) +
  geom_boxplot() +
  geom_point() +
  labs(x = 'Site', y = 'Observed ASVs', fill = "Colony size") +
  scale_x_discrete(limits = c("Langton", "Wainfleet", "Alliston", "Dunnville", "Mount Pleasant")) +
  scale_fill_viridis_d() #+
  #geom_text(data = site_summary, aes(site, Inf, label = n), vjust = 1)

glm_size_bbb = lmer(Observed ~ colony_size + (1|site), data = bbb_metadata)
summary(glm_size_bbb)
jtools::summ(glm_size_bbb) #Est. = 53.47, S.E. = 21.43, t = 2.49, d.f. = 1.53, p = 0.17

# 6. core genera ----
# set prevalence and detection levels
prevalences = seq(.05, 1, .05)
det <- c(0, 0.1, 0.5, 2)/100

# define a common theme for all plots
core_common_theme = theme_classic() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "none")

create_heatmap_plot <- function(pseq_data, prevalences, detections) {
  plot_core(pseq_data,
            plot.type = "heatmap",
            prevalences = prevalences,
            detections = det, 
            min.prevalence = .75)}

customize_plot <- function(plot, title) {
  plot + 
    core_common_theme +  # Apply common theme
    ggtitle(title) +  # Add title
    xlab("Detection threshold") +  # Label x-axis
    ylab("Genus\n") +  # Label y-axis
    scale_fill_viridis() +  # Color scale
    coord_fixed(ratio = 1)  # Adjust aspect ratio
}

# generate plots using functions
p_gen_langton <- create_heatmap_plot(pseq.rel.gen.langton.cp, prevalences, det)
p_gen_langton <- customize_plot(p_gen_langton, "Langton") # (n = 2)
p_gen_langton

p_gen_mtp <- create_heatmap_plot(pseq.rel.gen.mtp.cp, prevalences, det)
p_gen_mtp <- customize_plot(p_gen_mtp, "Mount Pleasant") #(n = 6)
p_gen_mtp

p_gen_dunnville <- create_heatmap_plot(pseq.rel.gen.dunnville.cp, prevalences, det)
p_gen_dunnville <- customize_plot(p_gen_dunnville, "Dunnville") #(n = 8)
p_gen_dunnville

p_gen_alliston <- create_heatmap_plot(pseq.rel.gen.alliston.cp, prevalences, det)
p_gen_alliston <- customize_plot(p_gen_alliston, "Alliston") #(n = 17)
p_gen_alliston

p_gen_wainfleet <- create_heatmap_plot(pseq.rel.gen.wainfleet.cp, prevalences, det)
p_gen_wainfleet = p_gen_wainfleet +
  theme_classic() +
  theme(axis.text.y = element_text(face = "italic")) +
  ggtitle("Wainfleet") +
  xlab("Detection threshold") +
  ylab("Genus\n") +
  scale_fill_viridis() + 
  coord_fixed(ratio = 1) #(n = 20)
p_gen_wainfleet
