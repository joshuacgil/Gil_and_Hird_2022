###R version 4.1.0###

###Required Packages###
ibrary(ggplot2); packageVersion("ggplot2")#v3.3.5
library(rlang);packageVersion("rlang")#v0.4.11
library(vegan); packageVersion("vegan")#v2.5.7
library(phyloseq); packageVersion("phyloseq")#v1.36.0
library(dplyr); packageVersion("dplyr")#v1.0.7
library(plyr); packageVersion("plyr")#v1.8.6
library(tibble); packageVersion("tibble")#v3.1.2
library(limma);packageVersion("limma")#3.48.3



###Load in the metadata and -omics data (KOs, ECs and Pathway data)###

metadata_omics <-read.csv(file="$metadata_file.csv", header=TRUE, sep = ",") #Load metadata metadata
otu_mat<- read.csv(file="$count_data.csv", header=TRUE, sep=',') #this is contains the number of KOs/ECs/Pathways
tax_mat<- read.csv(file="$descritpion_data.csv") #this holds the pathway/EC description data (Must have at least 2 columns)



###Make Venndiagrams for each samples###

venndf <- data.frame(otu_mat[0], 
                     +otu_mat$SampleX_Metagenome, 
                     +otu_mat$SampleX_Metatranscriptome,
                     +otu_mat$SampleX_SimulatedMetagenome)
head(venndf)
venn <- vennDiagram(venndf, include="both",
            names = c("MetaGene", "MetaTrans","16S-meta"),
            cex= 1, counts.col="Black")
venn 



###Create Phylsoeq Object using -omics data###

#define row names from the csv files

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("$ID")#change to first column's name (ID/KO/EC/pathway etc.)
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("$ID")#change to first column's name (ID/KO/EC/pathway etc.)

samples_df <- metadata_omics %>% 
  tibble::column_to_rownames("$ID") #change to name of the sample ID column 

#transform matrices (otu and taxa) as a matrix 

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#make the phyloseq object

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat) 
samples = sample_data(samples_df)

phyloseq_object <- phyloseq(OTU, TAX, samples)



###Ordinations###

#Jaccard#
set.seed(1)
phyloseq_jaccard<-ordinate(phyloseq_object, "PCoA", "jaccard", binary = TRUE)
ellipse.test<-plot_ordination(phyloseq_object, phyloseq_jaccard, color='$data_type', label=$NULL, shape='$data_type') +ggtitle("PCoA using Jaccard") + geom_point(size=4)
ellipse.test + stat_ellipse(geom="polygon",type="t",alpha = .2,linetype=1, aes(fill=data_type))



###PERMANOVA###
set.seed(1)

#Methods = jaccard as we only care about presence vs absence
jaccard_nova <- phyloseq::distancephyloseq_object method = "jaccard")

#Make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq))

# Adonis test
jac <- adonis(jaccard_nova ~ $data_type, data = sampledf, permutations = 999)
jac



###box plot###
boxplot($EC/KO/path_abundancedata~$data_type,
        data = metadata_omics,
        main = NULL,
        xlab = NULL,
        ylab = "Total EC/KO/Pathways per data type")



###box plot statistics###

#wilcox rank sum shows matrix of all pairwise comaprisons 
pairwise.wilcox.test((metadata_omics$EC/KO/Pathway_counts), metadata_omics$data_type, p.adjust.method='holm')


