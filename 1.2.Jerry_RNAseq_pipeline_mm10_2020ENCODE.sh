### Analyzing the RNASeq data
printf "\n!!!!!!! \nHi Jerry, a new journey of paired-end RNA-seq analysis has began.\n!!!!!!!\n\n"

#infile=(Test) 
#infile=(GingerasLab_8w_Cortex_mRNA_rep1 GingerasLab_8w_Cortex_mRNA_rep2)
#infile=(GingerasLab_8w_Frontal_Cortex_mRNA_rep1 GingerasLab_8w_Frontal_Cortex_mRNA_rep2)
#infile=(GingerasLab_8w_Cerebellum_mRNA_rep1 GingerasLab_8w_Cerebellum_mRNA_rep2 GingerasLab_8w_Urinary_Bladder_mRNA_rep1 GingerasLab_8w_Urinary_Bladder_mRNA_rep2)
#infile=(GingerasLab_8w_Placenta_mRNA_rep1 GingerasLab_8w_Placenta_mRNA_rep2 GingerasLab_8w_Subcutaneous_Adipose_mRNA_rep1 GingerasLab_8w_Subcutaneous_Adipose_mRNA_rep2) 
#infile=(GingerasLab_8w_Subcutaneous_Adipose_mRNA_rep1 GingerasLab_8w_Subcutaneous_Adipose_mRNA_rep2 GingerasLab_8w_Stomach_mRNA_rep1 GingerasLab_8w_Stomach_mRNA_rep2 GingerasLab_8w_Small_intestine_mRNA_rep1 GingerasLab_8w_Small_intestine_mRNA_rep2 GingerasLab_8w_Ovary_mRNA_rep1 GingerasLab_8w_Ovary_mRNA_rep2)
#infile=(GingerasLab_8w_Mammary_gland_mRNA_rep1 GingerasLab_8w_Mammary_gland_mRNA_rep2 GingerasLab_8w_Large_intestine_mRNA_rep1 GingerasLab_8w_Large_intestine_mRNA_rep2 GingerasLab_8w_Kidney_mRNA_rep1 GingerasLab_8w_Kidney_mRNA_rep2)
#infile=(GingerasLab_8w_Liver_mRNA_rep1 GingerasLab_8w_Liver_mRNA_rep2 GingerasLab_8w_Colon_mRNA_rep1 GingerasLab_8w_Colon_mRNA_rep2)
#infile=(GingerasLab_8w_Heart_mRNA_rep1 GingerasLab_8w_Heart_mRNA_rep2 GingerasLab_8w_Thymus_mRNA_rep1 GingerasLab_8w_Thymus_mRNA_rep2 GingerasLab_8w_Testis_mRNA_rep1 GingerasLab_8w_Testis_mRNA_rep2 GingerasLab_8w_Lung_mRNA_rep1 GingerasLab_8w_Lung_mRNA_rep2 GingerasLab_8w_Spleen_mRNA_rep1 GingerasLab_8w_Spleen_mRNA_rep2 GingerasLab_8w_Gonadal_FatPad_mRNA_rep1 GingerasLab_8w_Gonadal_FatPad_mRNA_rep2) 
infile=(GingerasLab_8w_Adrenal_gland_mRNA_rep1 GingerasLab_8w_Adrenal_gland_mRNA_rep2 GingerasLab_8w_Duodenum_mRNA_rep1 GingerasLab_8w_Duodenum_mRNA_rep2) 


for inputname in "${infile[@]}"
do
    cd "/media/Data/mENCODE/FASTQ"
    #### Step 1: mapping use STAR
    echo "Working on sample ${inputname}  \n\n"
    printf "\n\n\n  Step 1 of 3: STAR mapping \n\n"
  
    mkdir "tmp1" ## map the temp directory for mapping
    mv "${inputname}_1.fastq" "tmp1" ## move FASTQ to the folder
    mv "${inputname}_2.fastq" "tmp1" ## move FASTQ to the folder
    cd "tmp1"

    STAR --genomeDir "/media/Project/Genome/STAR/Ensembl_100_mm10" --readFilesIn "${inputname}_1.fastq" "${inputname}_2.fastq" --runThreadN 40 --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 --outFilterScoreMinOverLread 0.25 --outFilterMatchNminOverLread 0.25

    samtools view -f 0x2 "Aligned.out.sam" -o "${inputname}_pair.sam" ## remove singletons
    samtools view -H "Aligned.out.sam" -o "header.sam"  ## SAM headers 
    cat "header.sam" "${inputname}_pair.sam" > "${inputname}.sam" 
    #    mv "Aligned.out.sam" "${inputname}.sam" ## rename the alignment SAM file
    mv "Log.final.out" "${inputname}.Log" ## rename the mapping statistics file
    head -n 50 "${inputname}.Log" ## the mapping statistics

    mv "${inputname}.sam" "../../BAM"
    mv "${inputname}_1.fastq" ".."
    mv "${inputname}_2.fastq" ".."
#    rm "${inputname}_1.fastq"
#    rm "${inputname}_2.fastq"
    mv "${inputname}.Log" "/media/Project/mouse_RNAseq/mENCODE/Statistics"

    cd ".."
    rm -rf "tmp1"
#    gzip "${inputname}_1.fastq"
#    gzip "${inputname}_2.fastq"



    #### Step 2: split the SAM file by chromosome 
    printf "\n\n\nStep 2 of 3: split sam by chromosome \n\n"
    cd "/media/Data/Split"
    awk -v tag=${inputname} 'NR>25 { print > tag".sam.chr"$3}' "/media/Data/mENCODE/BAM/${inputname}.sam"



    #### Step 3: Generate bam and BigWig file from sam file
    printf "\n\n\nStep 3 of 3: generate bam and BigWig from SAM file \n\n"
    cd "/media/Data/mENCODE/BAM"
    # sam to bam
    samtools view -bS "${inputname}.sam" -o  "${inputname}_raw.bam" ## convert sam to bam
    samtools sort "${inputname}_raw.bam" -o "${inputname}.bam"      ## sort ba
    samtools index "${inputname}.bam"                               ## index bam
    rm "${inputname}_raw.bam"                                       ## remove unsorted bam

    # bam to BigWig
    lines=`expr $(wc -l < "/media/Data/mENCODE/BAM/${inputname}.sam"| tr -d " ") - 25` ## uniquely mapped reads
    bw_value=`expr $lines / 2000000` ### The normalized bw y-axes will be uniquely mapped read-pairs (Million)

    bamCoverage -b "${inputname}.bam" --filterRNAstrand forward --binSize 1 -p 14 -o "${inputname}_plus_${bw_value}.bw" 
    bamCoverage -b "${inputname}.bam" --filterRNAstrand reverse --binSize 1 -p 14 -o "${inputname}_minus_${bw_value}.bw"
    mv "${inputname}_plus_${bw_value}.bw" "/media/Project/mouse_RNAseq/mENCODE/Tracks" 
    mv "${inputname}_minus_${bw_value}.bw" "/media/Project/mouse_RNAseq/mENCODE/Tracks"

    rm "${inputname}.sam"
done

