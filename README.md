![](https://drive.google.com/uc?export=view&id=1r1BV4J16x2HGCmmvmGA56QpIkLK8Lerk)

# Mutational burden and ccRCC mortality
This is the repository with the whole code of our pipeline used in the analysis of mortality in metastatic clear cell kidney cancer carcinoma.


Clear cell renal cell carcinoma (ccRCC) is one of the most common kidney cancer forms and represents 2-3% of all cancers. This histologic subtype is characterized by malignant epithelial cells that show a clear cytoplasm and a compact-alveolar or acinar growth pattern interspersed with arborizing vasculature. It is described as one of the most aggressive histologic subtypes, with frequent metastatic staging at diagnosis. Around 60% of ccRCC patients die in the first 2–3 years from diagnosis.


Somatypus, a computational pipeline built upon the Platypus variant caller, was used to call genetic variation. After quality filtering, putative somatic variation was proposed if the variant was present in tumoral tissue but not in healthy tissue and was supported by at least 3 reads.


A gene-based burden of putative somatic mutations was calculated, and mutation enrichment assessed using Fisher’s exact test.
