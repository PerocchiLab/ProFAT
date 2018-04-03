#!/usr/bin/perl -w
use strict;
# Purpose:
# This function is used to retrived all the QC files
#
my $pro_dir="/data/home/cheng/Adipocyte";
my $data_info="$pro_dir/DataInfo/RNASeq_Dataset.NoCore.txt";

my $type="fq";  # "fq", "tg"
my $add="";   # "", ".tg"

open(IN, $data_info);
my $head=<IN>;
while(<IN>){
	chomp;
	my @arr=split(/\t/);
	my $study_name = $arr[0];
	my $accession=$arr[1];
	my @samples=split(/;/, $arr[2]);
	my @t_names=split(/;/, $arr[3]);
	my $link_dir="";
	my $pair_end = $arr[5];
	if ($accession=~/(\w\w\w\d\d\d)\d+/){
		print $1,"\n";
		$link_dir="$pro_dir/DataInfo/ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/$1/$accession";
	}elsif ($accession=~/E-/){
		print $accession,"\n";
		$link_dir="$pro_dir/DataInfo/ftp.sra.ebi.ac.uk/vol1/fastq/$accession";
	}

	for (my $ii=0;$ii<=$#samples;$ii++){
		if ($pair_end==1){
			my $sample_qc = $link_dir."/$type.$samples[$ii]/$samples[$ii]"."_1$add"."_fastqc/Images/per_base_quality.png";
			my $des_dir = "$pro_dir/RNASeq/$study_name/$type/$t_names[$ii].1.png";
			system("cp $sample_qc '$des_dir'");

			$sample_qc = $link_dir."/$type.$samples[$ii]/$samples[$ii]"."_2$add"."_fastqc/Images/per_base_quality.png";
			$des_dir = "$pro_dir/RNASeq/$study_name/$type/$t_names[$ii].2.png";
			system("cp $sample_qc '$des_dir'");
			
			# QC data	
			$sample_qc = $link_dir."/$type.$samples[$ii]/$samples[$ii]"."_1$add"."_fastqc/fastqc_data.txt";
			$des_dir = "$pro_dir/RNASeq/$study_name/$type/$t_names[$ii].1.txt";
			system("cp $sample_qc '$des_dir'");

			$sample_qc = $link_dir."/$type.$samples[$ii]/$samples[$ii]"."_2$add"."_fastqc/fastqc_data.txt";
			$des_dir = "$pro_dir/RNASeq/$study_name/$type/$t_names[$ii].2.txt";
			system("cp $sample_qc '$des_dir'");
		} else {
			my $sample_qc = $link_dir."/$type.$samples[$ii]/$samples[$ii]".$add."_fastqc/Images/per_base_quality.png";
			my $des_dir = "$pro_dir/RNASeq/$study_name/$type/$t_names[$ii].png";
			system("cp $sample_qc '$des_dir'");

			$sample_qc = $link_dir."/$type.$samples[$ii]/$samples[$ii]".$add."_fastqc/fastqc_data.txt";
			$des_dir = "$pro_dir/RNASeq/$study_name/$type/$t_names[$ii].txt";
			system("cp $sample_qc '$des_dir'");
		}
	}
}
close(IN);
	
