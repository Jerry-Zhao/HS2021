### Analyzing the single-nuclear RNA-seq data
# Prepare meta data
#  ## picard dictionary file : /Users/jerry/Analysis/Tools/Picard
#   java -jar /Users/jerry/Analysis/Tools/Picard/picard.jar CreateSequenceDictionary 
#        -R "../Mus_musculus.GRCm38.93.mm10.fa" -O Mus_musculus.GRCm38.mm10.dict
#   java -jar /Users/jerry/Analysis/Tools/Picard/picard.jar CreateSequenceDictionary 
#        -R "../Homo_sapiens.GRCh38.93.hg38.fa" -O Homo_sapiens.GRCh38.hg38.dict
#  ## refFlat similar like GTF
#   ConvertToRefFlat ANNOTATIONS_FILE=../Mus_musculus.GRCm38.100.gtf SEQUENCE_DICTIONARY=Mus_musculus.GRCm38.mm10.dict 
#     OUTPUT=Mus_musculus.GRCm38.100.mm10.refFlat
#   ConvertToRefFlat ANNOTATIONS_FILE=../Homo_sapiens.GRCh38.100.gtf SEQUENCE_DICTIONARY=Homo_sapiens.GRCh38.hg38.dict
#     OUTPUT=Homo_sapiens.GRCh38.100.hg38.refFlat
#  ## reduced.gtf
#   ReduceGTF SEQUENCE_DICTIONARY=Mus_musculus.GRCm38.mm10.dict GTF=../Mus_musculus.GRCm38.100.gtf 
#     OUTPUT=Mus_musculus.GRCm38.100.mm10.reduced.gtf
#   ReduceGTF SEQUENCE_DICTIONARY=Homo_sapiens.GRCh38.hg38.dict GTF=../Homo_sapiens.GRCh38.100.gtf 
#     OUTPUT=Homo_sapiens.GRCh38.100.hg38.reduced.gtf



printf "\n\n\n!!!!!!! Hi Jerry, a new journey of Drop-seq analysis has began.!!!!!!!\n\n\n" 

echo -n "Please enter the input FASTQ file: donot need '_1.fastq': "
read INname

echo -n "Please enter the sample name: Stress48  "
read SMname

#### Step 1: FASTQ to Picard-queryname-sorted BAM :: picard <FastqToSam>
# Converts a FASTQ file to an unaligned BAM or SAM file. 
# This tool extracts read sequences and base qualities from the input FASTQ file 
# and writes them out to a new file in unaligned BAM (uBAM) format.
# Input FASTQ files can be in GZip format (with .gz extension).
# SORT_ORDER (SortOrder)   The sort order for the output sam/bam file. Default value: queryname. 
#                          This option can be set to 'null' to clear the default value. 
#                          Possible values: {unsorted, queryname, coordinate, duplicate, unknown}

printf  "\n\n\nStep 1 of 5: FASTQ to BAM\n\n\n"
cd "/Users/jerry/Analysis/Project/scRNA/FASTQ"

command1="java -jar /Users/jerry/Analysis/Tools/Picard/picard.jar FastqToSam FASTQ=${INname}_1.fastq.gz FASTQ2=${INname}_2.fastq.gz OUTPUT=${INname}_uBAM.bam SAMPLE_NAME=$SMname"
eval "$command1"

mv "${INname}_uBAM.bam" "../Analysis"



#### Step 2: Pre-alignment tagging and trimming 
#### <TagBamWithReadSequenceExtended> <FilterBam> <TrimStartingSequence> <PolyATrimmer> 
printf  "\n\n\nStep 2 of 5: Pre-alignment tagging and trimming\n\n"
cd "../Analysis" 

# Tag Cell Barcode:
command2="TagBamWithReadSequenceExtended INPUT=${INname}_uBAM.bam OUTPUT=${INname}_uBAM_tagged_Cell.bam SUMMARY=${INname}_uBAM_tagged_Cellular.bam_summary.txt BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1"
eval "$command2"
rm "${INname}_uBAM.bam" 
mv "${INname}_uBAM_tagged_Cellular.bam_summary.txt" "Summary"

# Tag Molecular Barcode:
command3="TagBamWithReadSequenceExtended INPUT=${INname}_uBAM_tagged_Cell.bam OUTPUT=${INname}_uBAM_tagged_CellMolecular.bam SUMMARY=${INname}_uBAM_tagged_Molecular.bam_summary.txt BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1"
eval "$command3" 
rm "${INname}_uBAM_tagged_Cell.bam"
mv "${INname}_uBAM_tagged_Molecular.bam_summary.txt" "Summary"

# Remove reads with low quality 
command4=" FilterBam TAG_REJECT=XQ INPUT=${INname}_uBAM_tagged_CellMolecular.bam OUTPUT=${INname}_uBAM_tagged_filtered.bam"
eval "$command4"
rm "${INname}_uBAM_tagged_CellMolecular.bam"

# trim 5' end SMART adapter sequences (>=5 contiguous bases match)
command5="TrimStartingSequence INPUT=${INname}_uBAM_tagged_filtered.bam OUTPUT=${INname}_uBAM_tagged_trimmed_smart.bam OUTPUT_SUMMARY=${INname}_uBAM_adapter_trimming_report.txt SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5"
eval "$command5"
rm "${INname}_uBAM_tagged_filtered.bam"
mv "${INname}_uBAM_adapter_trimming_report.txt" "Summary"

# trim 3'end ployA sequences (>=6 A; AAAAAA)
command6="PolyATrimmer INPUT=${INname}_uBAM_tagged_trimmed_smart.bam OUTPUT=${INname}_uBAM_mc_tagged_polyA_filtered.bam OUTPUT_SUMMARY=${INname}_uBAM_polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6 USE_NEW_TRIMMER=true"
eval "$command6"
rm "${INname}_uBAM_tagged_trimmed_smart.bam"
mv "${INname}_uBAM_polyA_trimming_report.txt" "Summary"





#### Step 3: STAR alignment: <BAM to FASTQ, mapping>
printf  "\n\n\nStep 3 of 5: STAR alignment\n\n"

# BAm to FASTQ (single-end) using Picard SamToFastq
command7="java -Xmx20g -jar /Users/jerry/Analysis/Tools/Picard/picard.jar SamToFastq INPUT=${INname}_uBAM_mc_tagged_polyA_filtered.bam FASTQ=${INname}_clean.fastq"
eval "$command7"

# STAR alignment
mv "${INname}_clean.fastq" "mapping"
cd "mapping"
command8="STAR --genomeDir /Users/jerry/Analysis/Genome/STAR/mm10 --readFilesIn ${INname}_clean.fastq --runThreadN 13 --outFilterMultimapNmax 1 --outFilterMismatchNmax 3"
eval "$command8"

mv "Log.final.out" "${INname}_STAR_summary.txt"
mv "${INname}_STAR_summary.txt" "../Summary"

# Sort sam to bam by queryname 
command9="java -Xmx20g -jar /Users/jerry/Analysis/Tools/Picard/picard.jar SortSam I=Aligned.out.sam O=Aligned.sorted.bam SO=queryname"
eval "$command9" 



#### Step 4: Merging and tagging aligned reads
printf  "\n\n\nStep 4 of 5: Merging and tagging aligned reads\n\n"

# Recovery the missing tags, which are lost during the alignment
command10="java -Xmx20g -jar /Users/jerry/Analysis/Tools/Picard/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=/Users/jerry/Analysis/Genome/Dropseq_tools/Mus_musculus.GRCm38.mm10.fa UNMAPPED_BAM=../${INname}_uBAM_mc_tagged_polyA_filtered.bam ALIGNED_BAM=Aligned.sorted.bam OUTPUT=${INname}_merged.bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false"
eval "$command10"

# TagReadWithGeneFunction: 3 tags for each read: gn[gene name]; gs[gene strand]; gf[gene function]
command11="TagReadWithGeneFunction I=${INname}_merged.bam O=${INname}_gene_tagged.bam ANNOTATIONS_FILE=/Users/jerry/Analysis/Genome/Dropseq_tools/Mus_musculus.GRCm38.100.mm10.refFlat"
eval "$command11"

# DetectBeadSubstitutionErrors - Detecting and repairing substitution errors in cell barcodes
command12="DetectBeadSubstitutionErrors I=${INname}_gene_tagged.bam O=${INname}_gene_tagged_subtitution.bam OUTPUT_REPORT=${INname}_gene_tagged_clean.substitution_report.txt"
eval "$command12"

mv "${INname}_gene_tagged_subtitution.bam" ".."
mv "${INname}_gene_tagged_clean.substitution_report.txt" "../Summary"
rm -rf *





#### Step 5: Digital Gene Expression (DGE)
printf  "\n\n\nStep 5 of 5: Digital Gene Expression (DGE)\n\n" 
cd ".."

# Extracting Digital Gene Expression (DGE) data 
command13="DigitalExpression I=${INname}_gene_tagged_subtitution.bam O=${INname}_DGE.txt.gz SUMMARY=${INname}_DGE.summary.txt NUM_CORE_BARCODES=5000 LOCUS_FUNCTION_LIST=INTRONIC"
eval "$command13"
#### CODING+UTR+Intron
### LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC ### Intron only

mv "${INname}_DGE.summary.txt" "Summary"
rm "${INname}_uBAM_mc_tagged_polyA_filtered.bam"


