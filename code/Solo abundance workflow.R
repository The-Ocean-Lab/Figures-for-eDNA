#making the phylum graph #final_opti_mcc has sample names, otus
#I am using the old data and im going to filter to pr1 because
#The data I sent myself is weird for phylum

phylatax <- final_opti_mcc %>%
  select(-label, -numOtus) %>%
  rename(sample_id = Group) %>%
  pivot_longer(-sample_id, names_to = "otu", values_to = "counts") %>%
  filter(sample_id != "ALGAE") %>% # Taking out algae
  rename_all(tolower) %>%
  filter(!(sample_id %in% c("F24NC1", "F24NC1_2S", "F24NC1_3S", "F24NNR1", 
                            "F24NNR1_2S", "F24NNR1_3S", "JP", "JP_101124", 
                            "JP_101124_12", "JP_101124_13", "JP_101124_14", 
                            "JP_101124_15")))%>% # Taking out everything that's not PR1
  mutate(sample_id = recode(sample_id, 
                            "JP_101124_16" = "F24OPR1",
                            "JP_101124_17" = "F24OPR1_2",
                            "JP_101124_18" = "F24OPR1_3"))
otu_counts <- phylatax #changing the name so i can  follow my old code better
#I then renamed things, i noticed that my old PR1 samples i subset before ASLO were not the correct ones when compared
#to the wright labs list, so i had to redo this, not i am gonna make a phylum abundance and the genera abundance again with this

#I am going to do the same filtering of the big dataset for taxonomy then combined (im only doing this cuz i messed up the phylum subset)

taxonomy <- read_tsv("C:/Users/Pinel/OneDrive/Desktop/R code/Figures/final.opti_mcc.0.03.cons.taxonomy")%>%
  select("OTU","Taxonomy")%>%
  rename_all(tolower)%>%
  mutate(taxonomy = str_replace_all(taxonomy,"\\(\\d+\\)",""),
         taxonomy = str_replace(taxonomy, ";$", ""))
#stringr remove parenthesis (100) with number inside next to taxons 
#"\\(\\d+\\)" is a regular expression argument and basically is (d+) which means a number greater than 1 and 
#replaced it with nothing which is why its just quotes
taxonomy<- separate(taxonomy,
                    col="taxonomy",
                    into = c("kingdom", "phylum", "class", "order", "family", "genus"), 
                    sep = ";")

relativeabund<- inner_join(otu_counts,taxonomy)%>%
  mutate(rel_abund=counts/sum(counts))%>%
  select(-counts)%>%
  pivot_longer(cols=c("kingdom", "phylum", "class", "order", "family", "genus","otu"), 
               names_to = "level",
               values_to = "taxon")
#When I inner joined everything, it filtered to just PR1 for me, perhaps cuz i did it earlier

fallpr1phylaabundance_data <- relativeabund %>%
  filter(level == "phylum") #making a new file with just the phyla

fallpr1genusabundance_data <- relativeabund %>%
  filter(level == "genus")#making a new file with just the genus

#Combining the data is done and it has a relative abundance column
#I now am determining if the mean is needed for my graph and what that means

fallpr1phylaabundance_data %>%
  group_by(sample_id, taxon) %>% 
  summarize(rel_abund=sum(rel_abund),.groups = "drop")%>%
  arrange(desc(rel_abund))
#Why: Sum of Relative Abundance (First code) tells you how much of each taxon
#is present in each sample, which is useful for understanding the overall 
#proportion of a taxon in each sample, but it doesn't account for variation across samples.

averagefallpr1phyla <-fallpr1phylaabundance_data %>%
  group_by(sample_id, taxon) %>%
  summarize(mean_rel_abund = mean(rel_abund),.groups = "drop")%>%
  arrange(desc(mean_rel_abund))
#mean relative abundance
#Why:gives you the average relative abundance of each taxon, which helps to 
#understand the dominance of taxa in your data, particularly when comparing 
#across samples. Sorting by mean_rel_abund allows you to identify which taxa 
#tend to be more abundant across your dataset.

###############################################################################
#Lets make the phyla graph- This is my final code

fallpr1phylaabundancegraph <- fallpr1phylaabundance_data %>%
  group_by(sample_id, taxon) %>% 
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  # Scale the relative abundance to ensure it sums to 100% for each sample
  group_by(sample_id) %>%
  mutate(rel_abund_scaled = rel_abund / sum(rel_abund)) %>%
  ggplot(aes(x = sample_id, y = rel_abund_scaled, fill = taxon)) +
  geom_col() +
  scale_y_continuous(name = "Relative abundance", 
                     labels = scales::percent, 
                     limits = c(0, 1)) + # Set y-axis limits to 0-1 (100%)
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Optional: Rotate x-axis labels for readability
##################################################
#I now want to change the colors for funsies

# Count the number of unique phyla
num_phyla <- fallpr1phylaabundance_data %>%
  summarise(num_unique_phyla = n_distinct(taxon)) %>%  # Count the number of distinct phyla
  pull(num_unique_phyla)  # Extract the result from the summary

num_phyla #33


phylum_colors <- c(
  "grey22", "darkcyan", "orchid1", "green", "orange", "blue", "tomato2", "olivedrab", "grey47",
  "cyan", "coral3", "darkgreen", "magenta", "palegoldenrod", "dodgerblue", "firebrick", "yellow", 
  "purple4", "lightblue", "grey77", "mediumpurple1", "tan4", "red", "darkblue", "yellowgreen",
  "slateblue", "seagreen", "gold", "plum", "mediumvioletred", "darkorange", "chartreuse4", 
  "indianred", "royalblue3", "skyblue3"
)
length(phylum_colors) #35 colors

fallpr1phylaabundancegraph <- fallpr1phylaabundance_data %>%
  group_by(sample_id, taxon) %>% 
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  # Scale the relative abundance to ensure it sums to 100% for each sample
  group_by(sample_id) %>%
  mutate(rel_abund_scaled = rel_abund / sum(rel_abund)) %>%
  ggplot(aes(x = sample_id, y = rel_abund_scaled, fill = taxon)) +
  geom_col() +
  scale_fill_manual(values = phylum_colors)+
  scale_y_continuous(name = "Relative abundance", 
                     labels = scales::percent, 
                     limits = c(0, 1)) + # Set y-axis limits to 0-1 (100%)
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Optional: Rotate x-axis labels for readability

ggsave("FallPR1PhylaAbundance.png", plot =  fallpr1phylaabundancegraph, 
       width = 12, height = 8, dpi = 300)  
##################################################
#Average phyla across PR1 in the Fall Graph

averagefallpr1phyla <-fallpr1phylaabundance_data %>%
  group_by(sample_id, taxon) %>%
  summarize(mean_rel_abund = mean(rel_abund),.groups = "drop")%>%
  arrange(desc(mean_rel_abund))
#mean relative abundance
#Why:gives you the average relative abundance of each taxon, which helps to 
#understand the dominance of taxa in your data, particularly when comparing 
#across samples. Sorting by mean_rel_abund allows you to identify which taxa 
#tend to be more abundant across your dataset.

averagefallpr1phylagraph <- averagefallpr1phyla %>%
  # Scale the mean relative abundance to ensure it sums to 100% for each sample
  group_by(sample_id) %>%
  mutate(rel_abund_scaled = mean_rel_abund / sum(mean_rel_abund)) %>%
  ggplot(aes(x = sample_id, y = rel_abund_scaled, fill = taxon)) +
  geom_col() + 
  scale_fill_manual(values = phylum_colors)+
  scale_y_continuous(name = "Relative abundance", 
                     labels = scales::percent, 
                     limits = c(0, 1)) + # Set y-axis limits to 0-1 (100%)
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Optional: Rotate x-axis labels for readability

ggsave("AverageFallPR1PhylaAbundance.png", plot = averagefallpr1phylagraph, 
       width = 12, height = 8, dpi = 300)  


##############################################################################      
#Now I want the percentages of each phyla  so I can write them down

Falltop10pr1 <- fallpr1phylaabundance_data %>%
  filter(level == "phylum") %>%  # Filter for phylum-level taxa
  group_by(sample_id, taxon) %>%  # Group by sample and taxon (phylum)
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%  # Sum relative abundance for each taxon in each sample
  group_by(sample_id) %>%  # Group by sample to get the total for each sample
  mutate(perc_abund = rel_abund / sum(rel_abund) * 100) %>%  # Calculate percentage for each phylum in each sample
  arrange(sample_id, desc(perc_abund))  # Optionally arrange by sample_id and percentage

view(Falltop10pr1)

#The top 10 phyla were found in F24NPR1 so I wanted to compare to F24OPR1

top_ten_F240PR1 <- fallpr1phylaabundance_data %>%
  filter(sample_id == "F24OPR1", level == "phylum") %>%  # Filter for the correct sample and level
  group_by(sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%  # Sum the relative abundance for each taxon
  arrange(desc(rel_abund)) %>%  # Order by relative abundance
  head(10) %>%  # Select top 10 phyla
  mutate(perc_abund = rel_abund * 100)  # Calculate percentage
view(top_ten_P240PR1)