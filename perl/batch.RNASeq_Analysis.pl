#!/usr/bin/perl -w
use strict;
use Parallel::ForkManager;
# Purpose:
# This file is the only place to run the RNASeq_Analysis.pl
# All the datasets and samples will be defined here.
#
# Author: Yiming Cheng
# Created: 2015/07/13
# Modified: 2015/11/06 

my $base_dir = "/data/home/share/Projects/Adipocyte/RNASeq/";
my $project_dir = "/data/home/share/Projects/Adipocyte/";
my $info_dir = "/data/home/share/Projects/Adipocyte/DataInfo";
my $rnaseq_set = $ARGV[0];
if (!defined($rnaseq_set)){
	$rnaseq_set="RNASeq_Dataset.txt";
}
print $rnaseq_set,"\n";

my $info_file = $info_dir."/$rnaseq_set";
my $max_proc = 30;

open(IN, $info_file);
while(<IN>){
	chomp;
	if (/^([mho]).*_(\w+)\s+/){
		my $one_line=$_;
		my @arr=split(/\s+/);
		my $abbr_s=$1;
		my $two_t = $2;

		# Only process Mouse-Core Dataset and Human Dataset.
		if ($abbr_s eq "m"){
			if (($two_t ne "BeBr") && ($two_t ne "BeW") && ($two_t ne "BrW")){
				#print "Name not standard. $two_t\n";
				#next;
			}
		}

		my $ds_name=$arr[0];
		my $accession=$arr[1];
		my @samples=split(/;/,$arr[2]);
		my @tissues=split(/;/,$arr[3]);

		if ($#samples != $#tissues){
			print "Dataset error: $ds_name\n";
			next;
		}
		create_samples($one_line);

		my $t_dir=$base_dir.$ds_name."/samples/";
		my @s_files;
		opendir(DIR,$t_dir);
		my @text_files=readdir(DIR);
		for (my $kk=0;$kk<=$#text_files;$kk++){
			if ($text_files[$kk]=~/\.txt/){
				my $sample_file = $t_dir.$text_files[$kk];
				print "Running ./RNASeq_Analysis.pl $sample_file\n";
				system("./RNASeq_Analysis.pl $sample_file");
			}
		}

	}
}
close(IN);

######################################################
# Sub Function.
######################################################
sub create_samples{
	my @arr=split(/\s+/, shift);

	my $ds_name=$arr[0];
	my $acc=$arr[1];
	my $sample=$arr[2];
	my $tissue=$arr[3];
	my $do_diff = $arr[4];
	my $pair_end=$arr[5];

	my $mh_num="";
	my $species="";
	my $genome_base_dir="";
	if ($ds_name=~/^m(\d+)/){
		$genome_base_dir="/data/home/share/Data/iGenomes/Mus_musculus/Ensembl/release_81/";
		$species="Mouse";
		$mh_num=$1;
	} elsif ($ds_name=~/^h(\d+)/) {
		$genome_base_dir="/data/home/share/Data/iGenomes/Homo_sapiens/Ensembl/release_81/";
		$species="Human";
		$mh_num=$1;
	} elsif ($ds_name=~/^o(\d+)/){
		$genome_base_dir="/data/home/share/Data/iGenomes/Opossum/Ensembl/release_81/";
		$species="Opossum";
		$mh_num=$1;
	}
	my @fq_tg=("fq", "tg");
	my @samples=split(/;/, $sample);
	my @tissues=split(/;/,$tissue);

	my %tissue_sample;
	for (my $ii=0;$ii<=$#tissues;$ii++){
		push @{$tissue_sample{$tissues[$ii]}}, $samples[$ii];
	}

	my $sample_file_dir = $base_dir.$ds_name."/samples";
	if (-d $sample_file_dir){}else{
		system("mkdir -p $sample_file_dir");
	}

	foreach my $type (@fq_tg){
		my $out_file = $sample_file_dir."/$type.txt";
		#next if (-e $out_file);

		open(OUT,">$sample_file_dir/$type.txt");
		print OUT "#base_dir\t", $project_dir,"\n";
		print OUT "#genome_base_dir\t", $genome_base_dir,"\n";
		print OUT "#species\t$species\n";
		print OUT "#ftc\t$type\n";
		print OUT "#do_diff\t$do_diff\n";
		print OUT "#pair_end\t$pair_end\n";
		print OUT "#accession\t$acc\n";
		print OUT "#dataset_name\t$ds_name\n";

		foreach my $temp_tissue (sort {$a cmp $b} keys %tissue_sample){
			my @temp_samples = @{$tissue_sample{$temp_tissue}};
			if ($#temp_samples==0){
				print OUT $temp_samples[0],"\t", $temp_tissue,"\n";
				next;
			}
			for (my $kk=0;$kk<=$#temp_samples;$kk++){
				print OUT $temp_samples[0],"\t", $temp_tissue,"\n";
			}
		}
		close(OUT);
	}		
}

