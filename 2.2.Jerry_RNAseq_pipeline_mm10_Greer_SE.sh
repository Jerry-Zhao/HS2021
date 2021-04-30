### Analyzing the RNASeq data
printf "\n!!!!!!! \nHi Jerry, a new journey of paired-end RNA-seq analysis has began.\n!!!!!!!\n\n"

#infile=(GreerLab_2w_Cortex_totalRNA_rep1 GreerLab_2w_Cortex_totalRNA_rep2 GreerLab_2w_Cortex_totalRNA_rep3)
#infile=(GreerLab_8w_Cortex_totalRNA_rep1 GreerLab_8w_Cortex_totalRNA_rep2 GreerLab_8w_Cortex_totalRNA_rep3)
#infile=(Test)
#infile=(GreenbergLab_2w_Cortex_totalRNA_rep1 GreenbergLab_2w_Cortex_totalRNA_rep2 GreenbergLab_2w_Cortex_totalRNA_rep3)
#infile=(GreerLab_3d_Cortex_totalRNA_rep1 GreerLab_3d_Cortex_totalRNA_rep2 GreerLab_3d_Cortex_totalRNA_rep3 GreerLab_4w_Cortex_totalRNA_rep1 GreerLab_4w_Cortex_totalRNA_rep2 GreerLab_4w_Cortex_totalRNA_rep3 GreerLab_16w_Cortex_totalRNA_rep1 GreerLab_16w_Cortex_totalRNA_rep2 GreerLab_16w_Cortex_totalRNA_rep3)
#infile=(GreerLab_32w_Cortex_totalRNA_rep1 GreerLab_32w_Cortex_totalRNA_rep2 GreerLab_32w_Cortex_totalRNA_rep3 GreerLab_40w_Cortex_totalRNA_rep1 GreerLab_40w_Cortex_totalRNA_rep2 GreerLab_40w_Cortex_totalRNA_rep3 GreerLab_48w_Cortex_totalRNA_rep1 GreerLab_48w_Cortex_totalRNA_rep2 GreerLab_48w_Cortex_totalRNA_rep3 GreerLab_52w_Cortex_totalRNA_rep1 GreerLab_52w_Cortex_totalRNA_rep2 GreerLab_52w_Cortex_totalRNA_rep3)
infile=(GreenbergLab_2w_Cortex_Dnmt3aWT_totalRNA_rep1 GreenbergLab_2w_Cortex_Ezh2WT_totalRNA_rep1 GreenbergLab_2w_Cortex_Ezh2WT_totalRNA_rep2)

for inputname in "${infile[@]}"
do
    cd "/media/Data/Postnatal/FASTQ"
    #### Step 1: mapping use STAR
    echo "Working on sample ${inputname}  \n\n"
    printf "\n\n\n  Step 1 of 3: STAR mapping \n\n"
  
    mkdir "tmp1" ## map the temp directory for mapping
    mv "${inputname}.fastq" "tmp1" ## move FASTQ to the folder
    cd "tmp1"

    STAR --genomeDir "/home/jerry/Jerry/Genome/STAR/Ensembl_100_mm10" --readFilesIn "${inputname}.fastq" --runThreadN 40 --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 --outFilterScoreMinOverLread 0.25 --outFilterMatchNminOverLread 0.25

    mv "Aligned.out.sam" "${inputname}.sam" ## rename the alignment SAM file
    mv "Log.final.out" "${inputname}.Log" ## rename the mapping statistics file
    head -n 50 "${inputname}.Log" ## the mapping statistics

    mv "${inputname}.sam" "../../BAM"
    mv "${inputname}.fastq" ".."
#    rm "${inputname}.fastq"
    mv "${inputname}.Log" "/media/Project/mouse_RNAseq/Postnatal/Statistics"

    cd ".."
    rm -rf "tmp1"
#    gzip "${inputname}.fastq"



    #### Step 2: split the SAM file by chromosome 
    printf "\n\n\nStep 2 of 3: split sam by chromosome \n\n"
    cd "/media/Data/Split"
    awk -v tag=${inputname} 'NR>25 { print > tag".sam.chr"$3}' "/media/Data/Postnatal/BAM/${inputname}.sam"



    #### Step 3: Generate bam and BigWig file from sam file
    printf "\n\n\nStep 3 of 3: generate bam and BigWig from SAM file \n\n"
    cd "/media/Data/Postnatal/BAM"
    # sam to bam
    samtools view -bS "${inputname}.sam" -o  "${inputname}_raw.bam" ## convert sam to bam
    samtools sort "${inputname}_raw.bam" -o "${inputname}.bam"      ## sort ba
    samtools index "${inputname}.bam"                               ## index bam
    rm "${inputname}_raw.bam"                                       ## remove unsorted bam

    # bam to BigWig
    lines=`expr $(wc -l < "/media/Data/Postnatal/BAM/${inputname}.sam"| tr -d " ") - 25` ## uniquely mapped reads
    bw_value=`expr $lines / 1000000` ### The normalized bw y-axes will be uniquely mapped reads (Million)

    bamCoverage -b "${inputname}.bam" --filterRNAstrand forward --binSize 1 -p 14 -o "${inputname}_plus_${bw_value}.bw" 
    bamCoverage -b "${inputname}.bam" --filterRNAstrand reverse --binSize 1 -p 14 -o "${inputname}_minus_${bw_value}.bw"
    mv "${inputname}_plus_${bw_value}.bw" "/media/Project/mouse_RNAseq/Postnatal/Tracks" 
    mv "${inputname}_minus_${bw_value}.bw" "/media/Project/mouse_RNAseq/Postnatal/Tracks"

    rm "${inputname}.sam"
done

