#! /usr/local/perl -w
# Normalized to the total read/pairs that are uniquely mapped in exon regions. 
# In this way, the FPKM values of mRNA-seq and total RNA-seq are comparable. 

print "Please enter the Input file name, no need '_exonic_gene_Jerry.txt': GingerasLab_8w_Cortex_mRNA_rep1\n";
chomp($name=<STDIN>);

open(INA,"count/$name\_exonic_gene_Jerry.txt")||die("Can't open INA1:$!\n");

open(INGENOMEP,"/Users/jerry/Analysis/Genome/Ensembl_100/Mus_musculus.GRCm38.100.gtf_plus_final.exon")||die("Can't open INGENOMEP:$!\n");
open(INGENOMEM,"/Users/jerry/Analysis/Genome/Ensembl_100/Mus_musculus.GRCm38.100.gtf_minus_final.exon")||die("Can't open INGENOMEM:$!\n");
open(OUT,">FPKM/FPKM_$name\_ExonicGene.xls")||die("Can't write OUT:$!\n");


while(<INGENOMEP>)
{
 chomp;
 if(/^chr\w\w?\:(\d+)\-(\d+)\:(EN\w+\_[\w\.\-]+)\_exon\d+\t\+$/)
 {
  $exonlength=$2-$1+1;
  if(exists $exonlen{$3}){$exonlen{$3}+=$exonlength;}else{$exonlen{$3}=$exonlength;}
 }
 else{print "error exon: $_\n";}
}

while(<INGENOMEM>)
{
 chomp;
 if(/^chr\w\w?\:(\d+)\-(\d+)\:(EN\w+\_[\w\.\-]+)\_exon\d+\t\-$/)
 {
  $exonlength=$2-$1+1;
  if(exists $exonlen{$3}){$exonlen{$3}+=$exonlength;}else{$exonlen{$3}=$exonlength;}
 }
 else{print "error exon: $_\n";}
}

$unique_mapped=0;
while(<INA>)
{
 chomp;
 if(/^(EN\w+\_[\w\.\-]+)\t(\d+)$/){$count{$1}=$2;}
 elsif(/^####\s+unique-overlap\:\s+(\d+)\s*$/){$unique_mapped=$1;}
 else{print "error count: $_\n";}
}

foreach (keys %count)
{
 $fpkm=int($count{$_}*(1000000/$unique_mapped)*(1000/$exonlen{$_})*100)/100; ## FPKM. two digits 
 print OUT "$_\t$fpkm\n";
}

close INA; close OUT;

