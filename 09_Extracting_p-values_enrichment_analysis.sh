# Obtain the data from Enrichr ouputs
for size in {30..100}
do

indir="path/to/enrichr_output/${size}genes"
cd $indir

echo "---------------------------------------------------------------------"
echo "Workin in directory: ${indir}"
echo ""

for i in {1..1000}
do

echo ""
echo "Iteration (list number out of 1000): ${i}"
echo ""

# Define Jensen EnrichR outputfile
infile="enrichR_list_${size}genesiteration${i}"

# Grab the body
echo "Grabbing body..."
tail -n +2 $infile > body

# Mind that the infile has a header

# Grab "Kidney_cancer" and "Carcinoma" entries from Jensen plain-text EnrichR outputs
echo "Counting lines..."
wc -l body > temp3
awk '{ print $1 }' temp3 > numlines

# Grab the line with "Kidney_cancer" information (the line number and the information)
echo "Finding Kidney_cancer related genes..."
grep "Kidney_cancer" body > out1
#awk -F'[:]' '{ print $1 }' out1 > pos1

# Grab the line with "Carcinoma" information (the line number and the information)
echo "Finding Carcinoma related genes..."
grep "Carcinoma" body > out2
#awk -F'[:]' '{ print $1 }' out2 > pos2

# Save the results for "Kidney_cancer" into a single file (for APPEND)
echo "Saving Kidney_cancer related genes..:"
paste numlines out1 >> temp4

# Save the results for "Carcinoma" into a single file (for APPEND)
echo "Saving Carcinoma related genes..."
paste numlines out2 >> temp5

done

# Cat header and the results of all gene-lists
echo "Adding headers and outputting final files..."
# Grab the header
echo "Creating header..."
echo -e "Orden\tTotal" > temp1
head -n 1 $infile > temp2
paste temp1 temp2 > header
cat header temp4 > Kidney_cancer_Jensen_DISEASES_table_lista_${size}genes.txt
cat header temp5 > Carcinoma_Jensen_DISEASES_table_lista_${size}genes.txt
rm temp1 temp2 temp3 temp4 temp5 header numlines out1 out2 body

done

# Finally, cat all "Kidney_cancer" summary files into a single file (and the same for "Carcinoma")
for size in {30..32}
do

indir="path/to/enrichr_output/${size}genes"
infile1="${indir}/Kidney_cancer_Jensen_DISEASES_table_lista_${size}genes.txt"
infile2="${indir}/Carcinoma_Jensen_DISEASES_table_lista_${size}genes.txt"
outdir="/Users/alejandromendozaalvarez/Desktop/enrichr_output/summary"

# Do not include headers
tail -n +2 ${infile1} > ${indir}/body1
tail -n +2 ${infile2} > ${indir}/body2

# Append temporal outputs
cat ${indir}/body1 >> ${outdir}/temp1
cat ${indir}/body2 >> ${outdir}/temp2

done

# Build header for infile1
head -n 1 ${infile1} > ${outdir}/header

# Cat header and bodies into a single output file
cat ${outdir}/header ${outdir}/temp1 >> ${outdir}/Kidney_cancer_Jensen_DISEASES.txt
cat ${outdir}/header ${outdir}/temp2 >> >> ${outdir}/Carcinoma_Jensen_DISEASES.txt
rm ${outdir}/header ${outdir}/temp1 ${outdir}/temp2

# End of script
