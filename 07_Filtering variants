# Initial filtering for all variants

# Get a file with filtered INDELS: somatic INDELS
awk -F'[\t]' '$9 == "1" && $10 == "0" {print $0}' somatypus_output_indels.vcf > file1.vcf
head somatypus_output.vcf -n 1 > header1
cat header1 file1.vcf > somatics_Indels.vcf

# Get a file with filtered SNVs: somatic SNVs
awk -F'[\t]' '$9 == "1" && $10 == "0" {print $0}' somatypus_output_SNV.vcf > file2.vcf
head somatypus_output.vcf -n 1 > header2
cat header2 file2.vcf > somatics_SNVs.vcf

# Merge both somatic VCFs in one global VCF
cat somatics_SNVs.vcf somatics_Indels.vcf > somatics_without_header.vcf

# Somatic VCF sorting by position
sort -V -k1 -k2 somatics_without_header.vcf > somatics.ordered.vcf

# Get a list of somatic locus
cut -f1,2 somatics.ordered.vcf > somatics.ordered.list

# Extract the header of initial VCF from Somatypus
grep "^#" somatypus_output.vcf > header1

# Extract the body of initial VCF from Somatypus
grep -v "^#" somatics_without_header.vcf > body1

# Sorting again somatics variants detected
sort -V -k1 -k2 somatics_without_header.vcf > somatics_without_header_ordered.vcf

# Now place a header
cat header1 somatics_without_header.vcf > somatypus_final_somatics.vcf

# Check the final result with all somatic variants included
somatypus_final_somatics.vcf


#############################

# Filtering germline variants
# Get a file with filtered INDELS: germline INDELS
awk -F'[\t]' '$9 == "1" && $10 == "0" {print $0}' Somatypus_Indels_Presence.vcf > file1.vcf
head somatypus_output.vcf -n 1 > header1
cat header1 file1.vcf > germline_Indels.vcf

# Get a file with filtered SNVs: germline SNVs
awk -F'[\t]' '$9 == "1" && $10 == "0" {print $0}' Caregen_SNVs_Presence.vcf > file2.vcf
head somatypus_output.vcf -n 1 > header2
cat header2 file2.vcf > germline_SNVs.vcf

# Merge both somatic VCFs in one global VCF
cat somatics_SNVs.vcf germline_Indels.vcf > germline_without_header.vcf

# Somatic VCF sorting by position
sort -V -k1 -k2 germline_without_header.vcf > germline.ordered.vcf

# Get a list of somatic locus
cut -f1,2 germline.ordered.vcf > germline.ordered.list

# Extract the header of initial VCF from Somatypus
grep "^#" somatypus_output_germline.vcf > header1

# Extract the body of initial VCF from Somatypus
grep -v "^#" germline_without_header.vcf > body1

# Sorting again germline variants detected
sort -V -k1 -k2 germline_without_header.vcf > germline_without_header_ordered.vcf

# Now place a header
cat header1 germline_without_header.vcf > somatypus_final_germline.vcf

# Check the final result with all somatic variants included
somatypus_final_germline.vcf

#############################

# For ANNOVAR annotation of VCF using ExAC, avsnp150 and dbSNP:
table_annovar.pl VEP_Germinal.vcf /data/d2/jlorsal/ANNOVAR/humandb \
  -buildver hg19 \
  -out VEP_Germinal.annotated \
  -remove \
  -protocol exac03,avsnp150,dbnsfp33a \
  -operation f,f,f \
  -nastring . \
  -vcfinput

# For filtering variants without annotation in ExAc: 
grep -v "^#" VEP_Germinal.annotated.hg19_multianno.vcf | grep -vE "ExAC_ALL=\.;ExAC_AFR=\.;ExAC_AMR=\.;ExAC_EAS=\.;ExAC_FIN=\.;ExAC_NFE=\.;ExAC_OTH=\.;ExAC_SAS=\.;" > VEP_Germinal.annotated.hg19_multianno.ExAC.vcf
