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
#aa = proc.time()
test_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
#proc.time() - aa

# plot 96 profile of three samples
plot_96_profile(test_matrix[,c(1,4,7)])

# ------ EXTRACT SIGNATURES ------

# estimate rank
estimate_rank(test_matrix, rank_range = 2:5, nrun = 50)
# extract 3 signatures
nmf_res = extract_signatures(test_matrix, rank = 3)
# provide signature names (optional)
colnames(nmf_res$signatures) = c("Signature A", "Signature B" , "Signature C")
# plot signatures
plot_96_profile(nmf_res$signatures)

# provide signature names (optional)
rownames(nmf_res$contribution) = c("Signature A", "Signature B" , "Signature C")
# plot signature contribution relative
plot_contribution(nmf_res$contribution, coord_flip = T)
# plot signature contribution absolute

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

# ------ RAINFALL PLOT ------

# define chromosomes of interest
chromosomes = seqnames(get(ref_genome))[1:22]
# along chromosome
vcf = vcfs[[1]]
plot_rainfall(vcfs[[1]], title = names(vcfs[1]), ref_genome = ref_genome, chromosomes = chromosomes, cex = 1)
# for chromosome 1
plot_rainfall(vcfs[[1]], title = names(vcfs[1]), ref_genome = ref_genome, chromosomes = chromosomes[1], cex = 2)


# ------ GENOMIC DISTRIBUTION -------

bed_files = list.files("~/Documents/Organoids/ASC_multi_tissues/data/genomic_regions/", full.names = T)
region_names = c("H3k27ac", "H3K9me3", "early", "intermediate", "late", "exonic")
regions = bed_to_granges(bed_files, region_names = region_names)

surveyed_files = list.files("~/Documents/Organoids/ASC_multi_tissues/data/surveyed/", full.names = T)
surveyed_consensus = bed_to_granges(surveyed_files[1], "surveyed_consensus")
surveyed_list = c(surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus)


# rename chromosomes to UCSC standard
vcfs = lapply(vcfs, function(x) rename_chrom(x))
regions = lapply(regions, function(x) rename_chrom(x))
surveyed_list = lapply(surveyed_list, function(x) rename_chrom(x))

# GENOMIC DISTRIBUTION

x = genomic_distribution(vcfs[1:2], surveyed_list[1:2], regions[1:2])
df = enrichment_depletion_test(x)
plot_enrichment_depletion(df)

# ------ STRAND BIAS ----------

library("TxDb.Mmusculus.UCSC.mm10.knownGene")
genes_mm10 = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

seqlevels(vcf_ERCC1[[1]]) %in% seqlevels(get(mm10))
vcf_ERCC1 = lapply(vcf_ERCC1, function(x) rename_chrom(x, "UCSC"))

mut_mat = mut_matrix(vcf_ERCC1, mm10)
mut_mat_s = mut_matrix_stranded(vcf_ERCC1, mm10, genes_mm10)

sample_type = rep(c("mutant","wt"), each=3)
grid.arrange(plot_96_profile(mut_mat), plot_192_profile(mut_mat_s))

strand_bias = strand_bias(mut_mat_s, by=sample_type)

plot_strand_bias(strand_bias, mode = "absolute")

# statistical test for strand ratio
# poisson
df_strand = dcast(melt(strand_bias), type + group ~ strand, sum, subset = .(variable == "no_mutations"))
df_strand = ddply(df_strand, c("group", "type", "T", "U"), summarise, total = T+U, ratio = T/U,  p_poisson = poisson.test(c(U,T), c(total,total), r=1)$p.value)
df_strand$significant[df_strand$p_poisson < 0.05] = "*"

vcf = vcf_ERCC1[[1]]
genes = genes_mm10




