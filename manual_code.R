# USER MANUAL
library(MutationalPatterns)
library("gridExtra")
vcf_files = list.files(system.file("extdata", package="MutationalPatterns"), full.names = T)
sample_names = c("colon1", "colon2", "colon3", "intestine1", "intestine2", "intestine3", "liver1", "liver2", "liver3")
vcfs = read_vcf(vcf_files, sample_names, genome = "hg19")

ref_genome = "BSgenome.Hsapiens.UCSC.hg19" 
# source("http://bioconductor.org/biocLite.R")
# biocLite(ref_genome)
library(ref_genome, character.only = T)
# metadata
tissue = c("colon", "colon", "colon", "intestine", "intestine", "intestine", "liver", "liver", "liver")

# ----- CHROMOSOME NAMES ------

# check if chromosome names in your vcfs(s) and reference genomes are the same
all(seqlevels(vcfs[[1]]) %in% seqlevels(get(ref_genome)))
# rename the seqlevels of all vcfs in list to UCSC standara
vcfs = lapply(vcfs, function(x) rename_chrom(x))

# only select autosomal chromosomes, mt dna length is different for vcf and ref genome, why??
auto = extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="auto")
vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))

# ------ BASIC FUNCTIONS ----

muts = get_muts(vcfs[[1]])
types = get_types(vcfs[[1]])
mut_context = get_mut_context(vcfs[[1]], ref_genome)
type_context = get_type_context(vcfs[[1]], ref_genome)

type_occurences = mut_type_occurences(vcfs[1], ref_genome)
type_occurences = mut_type_occurences(vcfs, ref_genome)

# ----- SPECTRA ------

# plot spectrum over all samples
plot1 = plot_spectrum(type_occurences)
# with distinction of C>T at CpG sites
plot2 = plot_spectrum(type_occurences, CT = T)
# without legend
plot3 = plot_spectrum(type_occurences, CT = T, legend = F)
# define own colors for plotting
my_colors = c("pink", "orange", "blue", "lightblue", "green", "red", "purple")
plot4 = plot_spectrum(type_occurences, CT = T, legend = T, colors = my_colors)
# plot spectrum per tissue
plot5 = plot_spectrum(type_occurences, by = tissue, CT = T)

grid.arrange(plot1, plot2, plot3, ncol=3, widths=c(3,3,1.8))
grid.arrange(plot5, plot4, ncol=2, widths=c(3,2))


# ------ 96 SPECTRUM -----

# make 96 count matrix
test_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

# plot 96 profile of three samples
plot_96_profile(test_matrix[,c(1,4,7)])

# ------ EXTRACT SIGNATURES ------

# estimate rank
estimate_rank(test_matrix, rank_range = 2:5, nrun = 50)
# extract 2 signatures
nmf_res = extract_signatures(test_matrix, rank = 2)
# provide signature names (optional)
colnames(nmf_res$signatures) = c("Signature A", "Signature B")
# plot signatures
plot_96_profile(nmf_res$signatures)

# provide signature names (optional)
rownames(nmf_res$contribution) = c("Signature A", "Signature B")

# plot signature contribution
p1 = plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "relative")
p2 = plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute")
p3 = plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute", index = c(1,2))
p4 = plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute", coord_flip = T)

grid.arrange(p1,p2, nrow=2)
grid.arrange(p3,p4, ncol=2)


# compare reconstructed 96 profile of sample 1 with orignal profile
plot_compare_profiles(test_matrix[,1], nmf_res$reconstructed[,1], profile_names = c("Original", "Reconstructed"))

# ------ REFIT SIGNATURES ------

# download signatures from pan-cancer study Alexandrov et al.
sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
cancer_signatures = read.table(sp_url, sep = "\t", header = T)
# reorder (to make the order of the trinucleotide changes the same)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
# only signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])
# plot signature 1
plot_96_profile(cancer_signatures[,1,drop=F], ymax = 0.2)

fit_res = fit_to_signatures(test_matrix, cancer_signatures)
# select signatures with some contribution
select = which(rowSums(fit_res$contribution) > 0)
# plot contribution
plot_contribution(fit_res$contribution[select,], cancer_signatures[,select], coord_flip = F, mode = "absolute")

# compare reconstructed from refit with original profile
plot_compare_profiles(test_matrix[,1], fit_res$reconstructed[,1], profile_names = c("Original", "Reconstructed \n cancer signatures"))

# ------ STRAND BIAS ----------

# get knowngenes for hg19
source("https://bioconductor.org/biocLite.R")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes_hg19 = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

# get transcriptional strand information

get_strand(vcfs[[1]], genes_hg19)

# make mutation count matrix with transcriptional information
mut_mat_s = mut_matrix_stranded(vcfs, ref_genome, genes_hg19)

strand_counts = strand_occurences(mut_mat_s, by=tissue)
strand_plot = plot_strand(strand_counts , mode = "relative")
strand_bias = strand_bias_test(strand_counts)
strand_bias_plot = plot_strand_bias(strand_bias)

grid.arrange(strand_plot, strand_bias_plot)

# Extract signatures with strand bias

# extract 2 signatures
nmf_res_strand = extract_signatures(mut_mat_s, rank = 2)

# provide signature names (optional)
colnames(nmf_res_strand$signatures) = c("Signature A", "Signature B")
# plot signatures with 192 features
sig_stranded = plot_192_profile(nmf_res_strand$signatures)

# provide signature names (optional)
rownames(nmf_res_strand$contribution) = c("Signature A", "Signature B")
# plot signature contribution
plot_contribution(nmf_res_strand$contribution, nmf_res_strand$signatures, coord_flip = T, mode = "absolute")

# plot strand bias per mutation type for each signature with significance test
sig_strand_bias_plot = plot_signature_strand_bias(nmf_res_strand$signatures)

grid.arrange(sig_stranded, sig_strand_bias_plot, ncol=2, widths=c(3,1))

# ------ RAINFALL PLOT ------

# define chromosomes of interest
chromosomes = seqnames(get(ref_genome))[1:22]
# along chromosome
vcf = vcfs[[1]]
plot_rainfall(vcfs[[1]], title = names(vcfs[1]), ref_genome = ref_genome, chromosomes = chromosomes, cex = 1)
# for chromosome 1
plot_rainfall(vcfs[[1]], title = names(vcfs[1]), ref_genome = ref_genome, chromosomes = chromosomes[1], cex = 2)




# ------ GENOMIC DISTRIBUTION -------




# Use regulation annotation data from ENCODE

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)

listEnsembl()
listMarts()
mart="ensembl"
# list datasets available from ensembl for hg19 = GrCh37
listDatasets(useEnsembl(biomart="regulation", GRCh = 37))
# Homo sapiens Regulatory Segments
regulation_segmentation = useEnsembl(biomart="regulation", dataset="hsapiens_segmentation_feature", GRCh = 37)
regulation_regulatory = useEnsembl(biomart="regulation", dataset="hsapiens_regulatory_feature", GRCh = 37)
regulation_annotated = useEnsembl(biomart="regulation", dataset="hsapiens_annotated_feature", GRCh = 37)
# list all possible filters
listFilters(regulation_segmentation)
listFilters(regulation_regulatory)
# list all posible output attributes
listAttributes(regulation_regulatory)
# list all filter options for a specific attribute
filterOptions("feature_type_name", regulation_segmentation)
filterOptions("regulatory_feature_type_name", regulation_regulatory)

# Multicell regulatory features

# CTCF Binding Site

CTCF = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name', 'cell_type_name'), 
              filters = "regulatory_feature_type_name", 
              values = "CTCF Binding Site", 
              mart = regulation_regulatory)

CTCF_g = reduce(GRanges(CTCF$chromosome_name, IRanges(CTCF$chromosome_start, CTCF$chromosome_end))) 


promoter = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                 filters = "regulatory_feature_type_name", 
                 values = "Promoter", 
                 mart = regulation_regulatory)
promoter_g = reduce(GRanges(promoter$chromosome_name, IRanges(promoter$chromosome_start, promoter$chromosome_end)))

open = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
              filters = "regulatory_feature_type_name", 
              values = "Open chromatin", 
              mart = regulation_regulatory)
open_g = reduce(GRanges(open$chromosome_name, IRanges(open$chromosome_start, open$chromosome_end))) 

promoter_flanking = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
              filters = "regulatory_feature_type_name", 
              values = "Promoter Flanking Region", 
              mart = regulation_regulatory)
promoter_flanking_g = reduce(GRanges(promoter_flanking$chromosome_name, IRanges(promoter_flanking$chromosome_start, promoter_flanking$chromosome_end))) 

TF_binding = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                          filters = "regulatory_feature_type_name", 
                          values = "TF binding site", 
                          mart = regulation_regulatory)
TF_binding_g = reduce(GRanges(TF_binding$chromosome_name, IRanges(TF_binding$chromosome_start, TF_binding$chromosome_end))) 



regions = list(promoter_g, promoter_flanking_g, CTCF_g, open_g, TF_binding_g)
names(regions) = c("Promoter", "Promoter flanking", "CTCF", "Open chromatin", "TF binding")
# rename chromosomes to UCSC standard
regions = lapply(regions, function(x) rename_chrom(x))

# GENOMIC DISTRIBUTION TESTING

# Provide regions that are surveyed/callable
surveyed_file = list.files(system.file("bed", package="MutationalPatterns"), full.names = T)
# read bed file as grange object
surveyed_list = bed_to_granges(surveyed_file, "surveyed_all")
# for this example we use the same surveyed file for each sample
surveyed_list= rep(surveyed_list, 9)

# for each sample calculate the number of observed and expected number of mutations in each genomic regions
distr = genomic_distribution(vcfs, surveyed_list, regions)
# test for significant enrichment or depletion in the genomic regions
# samples can be collapsed into groups, here: tissue type
distr_test = enrichment_depletion_test(distr, by = tissue)
# plot enrichment depletion test results
plot_enrichment_depletion(distr_test)
