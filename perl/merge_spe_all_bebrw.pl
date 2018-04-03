#!/usr/bin/perl -w
use strict;
# Purpose: 
# This program is used to merge all the different studies into ONE. 
#
my $data_set=$ARGV[0];
my $data_spe="";
my $data_type="";
if (!defined($data_set)){
	print "$0 dataset_name\n";
	exit;
}
if ($data_set=~/^(\w)\d\d(\w)_\w+$/){
	$data_spe=$1;
	$data_type = $2;
} else {
	print "$0 dataset_name\n";
	exit;
}

my $base_dir="/data/home/cheng/Dev_Adipocyte/";
my %ortholog_MHO;
Load_Ortholog();

my @ftc=("fq", "tg");
my @count_files=("deseq.rlog", "fpkm", "seqcount");
my @all_bebrw_files=();
my @spe_files=();
my @out_files=();

if ($data_type eq "R"){
	foreach my $temp_ftc (@ftc){
	foreach my $temp_count (@count_files){
		push @spe_files, $base_dir."RNASeq/$data_set/$temp_ftc/$temp_count.txt";
		push @all_bebrw_files, $base_dir."RNASeq/All_BeBrW/$temp_ftc/All.$temp_count.txt";
		push @out_files, $base_dir."RNASeq/All_BeBrW/$temp_ftc/$data_set.$temp_count.txt";
	}}
} elsif ($data_type eq "M"){
	push @spe_files, $base_dir."MicroArray/$data_set/ND_$data_set.txt";
	push @all_bebrw_files, $base_dir."MicroArray/All_BeBrW/Data/All.ND.txt";
	push @out_files, $base_dir."MicroArray/All_BeBrW/Data/$data_set.ND.txt";
}
for (my $ii=0;$ii<=$#out_files;$ii++){
	my $all_bebrw_file=$all_bebrw_files[$ii];
	my $spe_file=$spe_files[$ii];
	my $out_file=$out_files[$ii];

	my @t_head;
	my (%data_bebrw, %data_all);

	open(IN, $all_bebrw_file);
	my $head=<IN>;
	chomp($head);
	$head=~s/^\t//;
	$head=~s/\t$//;
	push @t_head, $head;
	while(<IN>){
		chomp;
		my @arr=split(/\t/);
		@{$data_bebrw{$arr[0]}} = @arr[1..$#arr];
	}
	close(IN);

	open(IN, $spe_file);
	$head=<IN>;
	chomp($head);
	$head=~s/^\t//;
	$head=~s/\t$//;
	push @t_head, $head;
	while(<IN>){
		chomp;
		my @arr=split(/\t/);
		my $ortho_type="";
		if ($arr[0]=~/^ENSMODG/){
			$ortho_type="MO";
		} elsif ($arr[0]=~/^ENSG/){
			$ortho_type="MH";
		}
		if (length($ortho_type)==0){
			my $id=$arr[0];
			push @{$data_all{$id}}, @{$data_bebrw{$id}};
			push @{$data_all{$id}}, @arr[1..$#arr];
		} elsif (defined($ortholog_MHO{$ortho_type}{$arr[0]})){
			my $temp_orth=$ortholog_MHO{$ortho_type}{$arr[0]};
			if (defined($data_bebrw{$temp_orth})){
				my $id=$temp_orth.".".$arr[0];
				push @{$data_all{$id}}, @{$data_bebrw{$temp_orth}};
				push @{$data_all{$id}}, @arr[1..$#arr];
			}
		}

	}
	close(IN);

	open(OUT, ">$out_file");
	print OUT "\t", join("\t", @t_head),"\n";
	foreach my $temp_id (keys %data_all){
		print OUT $temp_id,"\t", join("\t", @{$data_all{$temp_id}}),"\n";
	}
	close(OUT);
}


######################################################################################################
## Functions
#######################################################################################################

sub Load_Ortholog{
	# Read Mouse-Human Ortholog pair
        my $file=$base_dir."Ortholog/Mouse_Human.txt";
        open(IN, $file);
        my $head=<IN>;
        while(<IN>){
                chomp;
                my @arr=split(/,/);
                my $mouse_id = $arr[0];
                my $human_id = $arr[2];
                my $ortho_type = $arr[5];
                next if (!defined($mouse_id));
                next if (!defined($human_id));
                next if (!defined($ortho_type));
                next if ($ortho_type ne "ortholog_one2one");
                $ortholog_MHO{"MH"}{$human_id} = $mouse_id;
                $ortholog_MHO{"MH"}{$mouse_id} = $human_id;
        }
        close(IN);

	# Read Opossum_Human.txt
        $file=$base_dir."Ortholog/Opossum_Human.txt";
        open(IN, $file);
        $head=<IN>;
        while(<IN>){
                chomp;
                my @arr=split(/,/);
                my $opossum_id = $arr[0];
                my $human_id = $arr[4];
                my $ortho_type = $arr[7];
                next if (!defined($opossum_id));
                next if (!defined($human_id));
                next if (!defined($ortho_type));
                next if ($ortho_type ne "ortholog_one2one");
                $ortholog_MHO{"HO"}{$human_id} = $opossum_id;
                $ortholog_MHO{"HO"}{$opossum_id} = $human_id;
        }
        close(IN);

	# Read Opossum_Mouse.txt
        $file=$base_dir."Ortholog/Opossum_Mouse.txt";
        open(IN, $file);
        $head=<IN>;
        while(<IN>){
                chomp;
                my @arr=split(/,/);
                my $opossum_id = $arr[0];
                my $mouse_id = $arr[4];
                my $ortho_type = $arr[7];
                next if (!defined($opossum_id));
                next if (!defined($mouse_id));
                next if (!defined($ortho_type));
                next if ($ortho_type ne "ortholog_one2one");
                $ortholog_MHO{"MO"}{$mouse_id} = $opossum_id;
                $ortholog_MHO{"MO"}{$opossum_id} = $mouse_id;
        }
        close(IN);
}

