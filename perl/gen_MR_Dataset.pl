#!/usr/bin/perl -w
use strict;
# Purpose:
# Created on 2017/01/25 by Yiming
# This function is used to generate the MicroArray_Dataset.txt and RNASeq_Dataset.txt
# based on the curated file core_dataset.txt

my $datainfo_dir="/data/home/cheng/Adipocyte/DataInfo";
my $file="$datainfo_dir/Core_Dataset.txt";
my $type = "NoCore";
if ($type eq "NoCore"){
	$file="$datainfo_dir/NoCore_Dataset.txt";
	if (-e $file){} else {
		print "The NoCore Dataset file does not exist\n";
		exit;
	}
}

open(IN, $file);
my $head=<IN>;
my @head_arr=split(/\t/, $head);
my %out_data;
while(<IN>){
	chomp;
	my @arr=split(/\t/);
	my $study_name=$arr[0];
	my $accession=$arr[3];
	my $platform=$arr[4];
	my $sampleID=$arr[5];
	my $type=$arr[1];
	
	$out_data{$study_name}{"Access"} = $accession;
	$out_data{$study_name}{"Platform"} = $platform;
	push @{$out_data{$study_name}{"Samples"}}, $sampleID;
	push @{$out_data{$study_name}{"Type"}}, $type;
}

my $microarray_file="$datainfo_dir/MicroArray_Dataset.txt";
if ($type eq "NoCore"){
	$microarray_file="$datainfo_dir/MicroArray_Dataset.NoCore.txt";
}

open(M,">$microarray_file");
print M join("\t", "Name", "Accession","Samples", "TissueName", "DoDiff", "PlatForm"),"\n";

my $rnaseq_file="$datainfo_dir/RNASeq_Dataset.txt";
if ($type eq "NoCore"){
	$rnaseq_file="$datainfo_dir/RNASeq_Dataset.NoCore.txt";
}

open(R,">$rnaseq_file");
print R join("\t", "Name", "Accession","Samples", "TissueName", "DoDiff", "PairEnd"),"\n";

my $do_diff=0;
foreach my $s_name (sort {$a cmp $b} keys %out_data){
	if ($s_name=~/\w\d+M/){
		print M join("\t", $s_name,$out_data{$s_name}{"Access"},join(";",@{$out_data{$s_name}{"Samples"}}),join(";",@{$out_data{$s_name}{"Type"}}),$do_diff, $out_data{$s_name}{"Platform"}),"\n";
	} elsif ($s_name=~/\w\d+R/){
		my $is_pairend=0;
		if ($out_data{$s_name}{"Platform"} eq "PairEnd"){
			$is_pairend=1;
		}
		print R join("\t", $s_name,$out_data{$s_name}{"Access"},join(";",@{$out_data{$s_name}{"Samples"}}),join(";",@{$out_data{$s_name}{"Type"}}),$do_diff, $is_pairend),"\n";
	}	
}
close(M);
close(R);

