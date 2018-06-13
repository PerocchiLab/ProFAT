#!/usr/bin/perl -w
use strict;
use Parallel::ForkManager;
# 
# Purpose: 
# This function is to run the RNASeq analysis in one script based on the samples.txt.
# This script is ONLY a reference script. To use it for productive, you could contact the ProFAT team 
# for more information. 
# The samples.txt provides all the necessary information of the RNASeq data and the reference genome.

############################################
# Initialize the data ...
############################################
print "Initialize the data ...\n";

my $sample_txt; # "samples.txt";
if ($#ARGV==0){
	$sample_txt = $ARGV[0];
} else {
	print "Please provide the sample description file \n";
	exit;
}

my $base_dir=""; 
my $genome_base_dir=""; 
my $species="";
my $ftc = ""; # fq, tc, tg fq: raw fastq file; tc: using trimmer and clipper; tg: use trim-galore program to trim the data 
my $pair_end = 0; # 0: double ended, 1: single ended.
my $accession="";
my $dataset_name="";
my $do_diff = 1;

my (@samples, @samples_syn, %group, %group_sample, @group_arr, %sample_dir, %sample_fastq_dir);
open(IN, $sample_txt);
while(<IN>){
	chomp;
	if (/^#(\w+)\s+(\S+)/){
		my $pro = $1;
		if ($pro eq "base_dir"){
			$base_dir = $2;
		}
		if ($pro eq "dataset_name"){
			$dataset_name=$2;
		}
		if ($pro eq "genome_base_dir"){
			$genome_base_dir = $2;
		}
		if ($pro eq "species"){
			$species = $2;
		}
		if ($pro eq "ftc"){
			$ftc = $2;
		} 
		if ($pro eq "pair_end"){
			$pair_end = $2;
		}
		if ($pro eq "accession"){
			$accession=$2;
		}
		if ($pro eq "do_diff"){
			$do_diff = $2;
		}
		next;
	}
	my @arr=split(/\s+/);
	my $sample_i = $arr[0];

	push @samples, $sample_i;
	push @samples_syn, $arr[1];
}
close(IN);

if (length($base_dir)<1 || length($genome_base_dir)<1){
	print "Please set up the base dir, genome base dir.\n";
	exit;
}

my $ftc_dir = $base_dir."RNASeq/".$dataset_name."/$ftc/";
if (-d $ftc_dir){
} else {
	mkdir($ftc_dir);
}

my $fastq_dir="";
my $data_type="";
if ($accession=~/^S/){
	my $acc_3=$accession;
	$acc_3=~s/\d\d\d$//;
	$fastq_dir=$base_dir."$acc_3/$accession/";
	$data_type="NCBI";
} elsif ($accession=~/^E/){
	$fastq_dir=$base_dir."$accession/";
	$data_type="EBI";
}

for (my $ii=0;$ii<=$#samples;$ii++){
	my $temp_dir = $fastq_dir.$ftc.".".$samples[$ii]."/";
	if (-d $temp_dir){
	} else {
		mkdir($temp_dir);
	}	
	$sample_dir{$samples[$ii]} = $temp_dir;
}

my $bam_file_name = "accepted_hits.bam";
my $sam_file_name = "accepted_hits.sam";

############################################
# Build the Bowtie ...
############################################
my $bowtie_index = $genome_base_dir."Sequence/Bowtie2Index/genome";
my $genome_fasta = $genome_base_dir."Sequence/Bowtie2Index/genome.fa";
my $gtf_file = $genome_base_dir."Annotation/Genes/genes.gtf";
my $f_index = $bowtie_index.".1.bt2";
if (file_es($f_index)){
} else {
	my $build_ind = "bowtie2-build $genome_fasta $bowtie_index";
	print $build_ind,"\n";
	system($build_ind);
}

############################################
# Do the FastQC to check the quality
############################################
print "FastqQC Check ...\n";

my $num_proc = 16;
my $pm=new Parallel::ForkManager($num_proc);

for (my $ii=0;$ii<=$#samples;$ii++){
	my $sample_i = $samples[$ii];

	$pm->start and next;

	my $sample_ftc_dir = $sample_dir{$sample_i};
	
	# get the fastq files for one sample.
	my @fastq_files;
	if ($pair_end==1){
		push @fastq_files, $fastq_dir.$sample_i."_1.fastq";
		push @fastq_files, $fastq_dir.$sample_i."_2.fastq";
	} else {
		push @fastq_files, $fastq_dir.$sample_i.".fastq";
	}

	# Handle different triming methods. 
	if ($ftc eq "fq"){ 
		# run fastq
		my @qc_res_files;
		if ($pair_end==1){
			push @qc_res_files, $sample_ftc_dir.$sample_i."_1_fastqc.html";
			push @qc_res_files, $sample_ftc_dir.$sample_i."_2_fastqc.html";
		} else {
			push @qc_res_files, $sample_ftc_dir.$sample_i."_fastqc.html";
		}

		# Run the FastQC.
		for (my $kk=0;$kk<=$#fastq_files;$kk++){
			if (file_es($qc_res_files[$kk])){
			} else {
				system("cp $fastq_files[$kk] $sample_ftc_dir");
				my $run_fastqc = "fastqc $fastq_files[$kk] --extract -o $sample_ftc_dir";
				print $run_fastqc,"\n";
				system($run_fastqc);
			}
		}
	} elsif ($ftc eq "tg"){
		my $tg_program = "Tools/Trim_Galore/trim_galore";
		my $run_trim_galore="";
		my $temp_fq_file = join(" ", @fastq_files);
		if ($pair_end==1){
			$run_trim_galore= "$tg_program --paired $temp_fq_file -o $sample_ftc_dir";
		} else {
			$run_trim_galore= "$tg_program $temp_fq_file -o $sample_ftc_dir";
		}

		my (@fastq_tg, @trimmed_file_names);
		if ($pair_end==1){
			push @fastq_tg, $sample_ftc_dir.$sample_i."_1.tg";
			push @fastq_tg, $sample_ftc_dir.$sample_i."_2.tg";

			push @trimmed_file_names, $sample_ftc_dir.$sample_i."_1_val_1.fq";
			push @trimmed_file_names, $sample_ftc_dir.$sample_i."_2_val_2.fq";
		} else {
			push @fastq_tg, $sample_ftc_dir.$sample_i.".tg";
			push @trimmed_file_names, $sample_ftc_dir.$sample_i."_trimmed.fq";
		}

		if (file_es($fastq_tg[0])){
		} else {
			print $run_trim_galore,"\n";
			system($run_trim_galore);

			for (my $kk=0;$kk<=$#fastq_tg;$kk++){
				system("mv $trimmed_file_names[$kk] $fastq_tg[$kk]");
				system("fastqc $fastq_tg[$kk] --extract -o $sample_ftc_dir");
			}
		}
	}

	$pm->finish;
}
$pm->wait_all_children();

############################################
# Do the TOPHAT alignment ...
############################################
print "Do the TOPHAT alignment ...\n";

for (my $ii=0;$ii<=$#samples;$ii++){
	my $sample_i = $samples[$ii];

	$pm->start and next;

	my $sample_ftc_dir = $sample_dir{$sample_i};
	my $ftc_label=$ftc;
	if ($ftc eq "fq"){
		$ftc_label="fastq";
	}
	
	my @fq_files;
	if ($pair_end==1){
		push @fq_files, $sample_ftc_dir.$sample_i."_1.".$ftc_label;
		push @fq_files, $sample_ftc_dir.$sample_i."_2.".$ftc_label;
	} else {
		push @fq_files, $sample_ftc_dir.$sample_i.".".$ftc_label;
	}

	my $bam_file = $sample_ftc_dir.$bam_file_name;
	if (file_es($bam_file)){
	} else {
		# Run the tophat: align the FASTQ file to the genome
		my $fastq_file = join(" ", @fq_files);
		my $run_tophat = "tophat -p $num_proc -G $gtf_file -o $sample_ftc_dir $bowtie_index $fastq_file";
		print $run_tophat,"\n";
		system($run_tophat);
	}
	
	my $sam_file = $sample_ftc_dir.$sam_file_name;
	if (file_es($sam_file)){
	} else {
		my $run_samtools = "samtools view -h -o $sam_file $bam_file";
		print $run_samtools,"\n";
		system($run_samtools);
	}
	$pm->finish;
}
$pm->wait_all_children();

############################################
# Do the  Cufflinks ...
############################################
print "Do the cufflinks ...\n";

for (my $ii=0;$ii<=$#samples;$ii++){
	my $sample_i = $samples[$ii];

	$pm->start and next;

	my $sample_ftc_dir = $sample_dir{$sample_i};

	my $bam_file = $sample_ftc_dir.$bam_file_name;
	my $gene_tracking_file = $sample_ftc_dir."genes.fpkm_tracking";
	if (file_es($gene_tracking_file)){
	} else {
		# Run the cufflinks: assembles transcripts, estimates their abundances
		my $run_cufflinks = "cufflinks -p $num_proc -G $gtf_file -o $sample_ftc_dir $bam_file";
		print $run_cufflinks,"\n";
		system($run_cufflinks);
	}

	$pm->finish;
}
$pm->wait_all_children();

############################################
# Do the HTSeq-count ...
############################################
print "Do the HTSeq-count ...\n";

for (my $ii=0;$ii<=$#samples;$ii++){
	my $sample_i = $samples[$ii];

	$pm->start and next;

	my $sample_ftc_dir = $sample_dir{$sample_i};
	my $bam_file = $sample_ftc_dir.$bam_file_name;
	my $sam_file = $sample_ftc_dir.$sam_file_name;
	my $htseq_out = $sample_ftc_dir."htseq_count.txt";
	if (file_es($htseq_out)){
	} else {
		my $run_htseq_count = "htseq-count -s no $sam_file $gtf_file $htseq_out"; # -s no: stranded=no
		print $run_htseq_count,"\n";
		system($run_htseq_count);
	}
	$pm->finish;
}
$pm->wait_all_children();

#################################################################################
# Retrieve the FPKM for all the samples.   
#################################################################################
my (%data_fpkm, %data_seq_count);

my $fpkm_output = $ftc_dir."fpkm.txt";
my $htseq_output = $ftc_dir."seqcount.txt";

if (file_es($fpkm_output)){
} else {
	print "Generate the fpkm_out\n";
	for (my $ii=0;$ii<=$#samples;$ii++){
		Load_FPKM($samples[$ii]);
	}
	Write_FPKM_Output($fpkm_output);
}

if (file_es($htseq_output)){
} else {
	print "Generate the seq count file\n";
	for (my $ii=0;$ii<=$#samples;$ii++){
		Load_HTSeqCount($samples[$ii]);
	}
	Write_HTSeq_Output($htseq_output);
}

#################################################################################
# #######################################
# Sub Function
# #######################################
#################################################################################
sub gen_deseqn{
	my $in_file=shift;
	my $out_file =shift;
	my $n_rlog = shift;
	open(IN, $in_file);
	my $head = <IN>;
	chomp($head);
	my @head_arr=split(/\t/,$head);
	my @new_head_arr;
	my @new_head_arr_ind;
	for (my $ii=0;$ii<=$#head_arr;$ii++){
		if ($head_arr[$ii]=~/(.*)_$n_rlog/){
			push @new_head_arr,$1;
			push @new_head_arr_ind, $ii;
		}
	}
	open(OUT, ">$out_file");
	print OUT "\t", join("\t",@new_head_arr),"\n";
	while(<IN>){
		chomp;
		my @arr=split(/\t/);
		my @new_row=();
		for (my $ii=0;$ii<=$#new_head_arr_ind;$ii++){
			push @new_row, $arr[$new_head_arr_ind[$ii]];
		}
		print OUT join("\t", $arr[0], @new_row),"\n";
	}
	close(IN);
	close(OUT);
}


sub file_es{ # check if file exists and has no-zero size
	my $file=shift;
	my $exist_size = 0;
	if (-e $file && -s $file){
		$exist_size = 1;
	}
	return($exist_size);
}

sub Write_FPKM_Output{
	my $file_name = shift;
	open(OUT,">$file_name");

	print OUT "\t", join("\t", @samples_syn),"\n";
	foreach my $temp_gene (keys %data_fpkm){
		my @sample_fpkm=();
		my $sym="";
		for (my $ii=0;$ii<=$#samples;$ii++){
			push @sample_fpkm, $data_fpkm{$temp_gene}{$samples[$ii]}{"FPKM"};
			$sym = $data_fpkm{$temp_gene}{$samples[$ii]}{"gene_short_name"};
		}
		print OUT $temp_gene,"\t", join("\t", @sample_fpkm),"\n";
	}
	close(OUT);
}

sub Write_HTSeq_Output{
	my $file_name = shift;
	open(OUT,">$file_name");

	print OUT "\t", join("\t", @samples_syn),"\n";
	foreach my $temp_gene (keys %data_seq_count){
		next if ($temp_gene=~/^_/);
		my @sample_count=();
		my $sym="";
		for (my $ii=0;$ii<=$#samples;$ii++){
			push @sample_count, $data_seq_count{$temp_gene}{$samples[$ii]};
		}
		print OUT $temp_gene,"\t", join("\t", @sample_count),"\n";
	}
	close(OUT);
}

sub Load_FPKM{
	my $sample =shift;
	my $sample_ftc_dir = $sample_dir{$sample};
	my $gene_fpkm_file = $sample_ftc_dir."genes.fpkm_tracking";
	open(FPKM, $gene_fpkm_file);
	my $head=<FPKM>;
	my @head_arr=split(/\t/, $head);
	while(<FPKM>){
		chomp;
		my @arr=split(/\t/);
		for (my $ii=1;$ii<=$#arr;$ii++){
			if ($head_arr[$ii] eq "FPKM"){
				$data_fpkm{$arr[0]}{$sample}{$head_arr[$ii]} += $arr[$ii];
			} else {
				$data_fpkm{$arr[0]}{$sample}{$head_arr[$ii]} = $arr[$ii];
			}
		}
	}
	close(FPKM);
}
 
sub Load_HTSeqCount{
	my $sample =shift;
	my $sample_ftc_dir = $sample_dir{$sample};
	my $htseq_count_file = $sample_ftc_dir."htseq_count.txt";
	open(COUNT, $htseq_count_file);
	while(<COUNT>){
		chomp;
		my @arr=split(/\s+/);
		$data_seq_count{$arr[0]}{$sample} = $arr[1];
	}
	close(COUNT);
}
