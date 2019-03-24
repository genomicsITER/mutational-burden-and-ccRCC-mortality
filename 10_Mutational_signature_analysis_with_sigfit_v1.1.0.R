##############################################################
#     MUTATIONAL SIGNATURE ANALYSIS WITH SIGFIT v.1.1.0      #
##############################################################

# Mutational signature analysis using sigfit (v1.1.0)
# https://github.com/kgori/sigfit

# Create output directory
OUT.DIR <- "path/to/outdir/results"
dir.create(OUT.DIR, showWarnings = FALSE)

# Read VCF
VCF.FILE <- "somatics_mutation_file.vcf"
vcf <- read.table(VCF.FILE, sep = "\t", comment.char = "", 
                  header = T, check.names = F, stringsAsFactors = F)
colnames(vcf)[1] <- "CHROM"

# Discard indels
vcf <- vcf[nchar(vcf$REF) == 1 & nchar(vcf$ALT) == 1, ]

# Create presence/absence table
presence <- 
    vcf[, -(1:9)] == "0/1" | 
    vcf[, -(1:9)] == "1/1" | 
    vcf[, -(1:9)] == "1/0"
presence <- presence[, c(2:10, 1)]
colnames(presence) <- paste0("P", 1:10)

# Create mutation data table for sigfit
vcf$Trinuc <- sapply(vcf$INFO, function(info) { 
    x <- strsplit(info, split = ";")[[1]][12]
    substr(x, 13, 15)
})
stopifnot(all.equal(substr(vcf$Trinuc, 2, 2), vcf$REF))

muts <- NULL
for (j in 1:ncol(presence)) {
    index <- presence[,j]
    muts <- rbind(muts, cbind("Sample" = colnames(presence)[j],
                              "Ref" = vcf$REF[index],
                              "Alt" = vcf$ALT[index],
                              "Trinuc" = vcf$Trinuc[index]))
}
muts <- as.data.frame(muts)
stopifnot(nrow(muts) == sum(presence))


# Build mutational catalogues (sigfit v1.1.0)
library(sigfit)
spectra <- build_catalogues(muts)

# Adapt COSMIC signatures to exome trinucleotide composition
data("cosmic_signatures")
cosmic.norm <- convert_signatures(cosmic_signatures, 
                                  ref_opportunities = "human-genome", model_to = "emu")
cosmic.exome <- convert_signatures(cosmic.norm, 
                                   ref_opportunities = "human-exome", model_to = "nmf")

# Initial fitting all COSMIC signatures
sig.fit.1 <- fit_signatures(spectra, cosmic.exome, iter = 5000)

# Plot fitted signature exposures to identify significantly supported signatures
plot_exposures(sig.fit.1, pdf_path = file.path(OUT.DIR, "Signature_Exposures_Initial.pdf"))

# Refined fitting of significant signatures only
# COSMIC signatures 1 and 12 have significantly non-zero exposures;
# we include sig. 5 on the basis of previous knowledge about its universality;
# we exclude sig. 11 on the basis of previous knowledge about its aetiology (temozolomide);
# we exclude sig. 30 on the basis of previous knowledge about its extreme rarity.
# (https://cancer.sanger.ac.uk/cosmic/signatures)
signif.sigs.idx <- c(1, 5, 12)
sig.fit.2 <- fit_signatures(spectra, cosmic.exome[signif.sigs.idx, ], iter = 10000)

# Plot catalogues, signatures, exposures and reconstructions
plot_all(sig.fit.2, out_path = file.path(OUT.DIR, "Signature_Fitting_Definitive"))
