#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Aug 24, 2010

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my $locallib = eval{
			require local::lib;
			local::lib->import();
			1;
		};
if ($locallib){
	use local::lib;
}


use Getopt::Long;
use List::Compare;
use POSIX;

=head1 NAME

 orthomcl.ValidateGroups.pl - compare a group set with two other sets, one relaxed param and one stringent param 

=head1 SYNOPSIS

  % orthomcl.ValidateGroups.pl --selgroups groups.txt --consparam groups.txt --relparam groups.txt --map map.txt --names names.txt
  
=head1 DESCRIPTION

 Finds the robustness of groups found by OrthoMCL. Reads in groups file produced by OrthoMCL. Compare the groups to groups 
 produced using most relaxed and most stringent parameters. Prints out statistics about groups where the 
 list of member proteins changed. Creates new group files with common,80+,60+,40+ conserved clusters.
 
 1. Clusters that retain 100% to 80% members across runs (highly conserved).
 2. Clusters that retain 80% to 60% members across runs (moderately conserved).
 3. Clusters that retain 60% to 40% members across runs (poorly conserved).
 
 Use 'cat *.faa| fgrep '>' > names.txt' to get all the protein names into a single file. This tool 
 only works for Refseq style FASTA headers.
 >gi|16262454|ref|NP_435247.1| FdoG formate dehydrogenase-O alpha subunit [Sinorhizobium meliloti 1021]

=head1 VERSION HISTORY

 Version   1.0  INPUT: 3 groups file, Mapping file, Fasta headers 
                OUTPUT: .XLS and group files (common,80+,60+,40+ conserved clusters) 
             
=head1 TODO
See github issues https://github.com/suryasaha/OrthoMCL_PostProcessing/issues

=head1 NOTE
 Make sure your protein FAA files are in Refseq format
 Example: >gi|16262454|ref|NP_435247.1| FdoG formate dehydrogenase-O alpha subunit [Sinorhizobium meliloti 1021]
 Make sure that all gi numbers are unique in the Fasta headers. If you use orthomcl.convert2RefseqFasta.pl to format the fasta headers, then make sure you use a seed value higher than existing values.
 
=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --selgroups <.txt>    Text file with protein clusters from OrthoMCL run to validate (required)
   --consparam <.txt>    Text file with protein clusters from OrthoMCL run with most conservative parameters (required)
   --relparam  <.txt>    Text file with protein clusters from OrthoMCL run with most relaxed parameters (required)
   --map       <.txt>)   Tab delimited file with abbreviations followed by full genome names used in Refseq (required)
   --names     <.txt>    Text file with protein FASTA headers (required)
      
=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut

sub mk_key{
	my (@group_arr,@prot_arr,@gi_arr,@gi_sorted_arr,$i,$key);
	@group_arr=split(' ',$_[0]);
	@gi_arr=();
	for $i (1..(scalar(@group_arr)-1)){
		@prot_arr=split(/\|/,$group_arr[$i]);#can be shortened skipping temp arrays
		push @gi_arr,$prot_arr[1];#record GIs
	}
	@gi_sorted_arr = sort {$a <=> $b} @gi_arr;
	foreach $i (@gi_sorted_arr){#make key of sorted GIs
		if(!defined($key)){$key=$i;}
		else{$key=$key.'_'.$i;}
	}
	return $key;
}

sub validate_protein_gi_number{
	my (@group_arr,@prot_arr,@gi_arr,$i);
	my $err_string = 0;
	@group_arr=split(' ',$_[0]);
	@gi_arr=();
	for $i (1..(scalar(@group_arr)-1)){
		@prot_arr=split(/\|/,$group_arr[$i]);#can be shortened skipping temp arrays
		if ( !(looks_like_number($prot_arr[1])) ) { 
			$err_string = "$prot_arr[1] is not a valid GI number";
			last;
		}
	}
	return $err_string;
}

sub validate_refseq_format{
	my $fasta_header = shift;
	my @fasta_header_arr = split(/\|/,$fasta_header);
	my $err_string = 0;
	if ( scalar @fasta_header_arr < 5 ) { $err_string = "less than 5 values";}
	if ( scalar @fasta_header_arr > 5 ) { $err_string = "more than 5 values";}
	if ( $fasta_header_arr[0] ne '>gi' ) { $err_string = "first value is not GI";}
	if ( !(looks_like_number($fasta_header_arr[1])) ) { $err_string = "$fasta_header_arr[1] is not a valid GI number";}
	if ( ($fasta_header_arr[4] !~ /\[/) || ($fasta_header_arr[4] !~ /\]/)   ) { $err_string = "genome name not present or not enclosed by []";}
	
	return $err_string;
}


my ($i,$sgrp,$cgrp,$rgrp,$map,$names);

GetOptions (
	'selgroups=s' => \$sgrp,
	'consgroups=s' => \$cgrp,
	'relgroups=s' => \$rgrp,
	'map=s' => \$map,
	'names=s' => \$names) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($sgrp) or (system('pod2text',$0), exit 1);
if (!(-e $sgrp)){print STDERR "$sgrp not found: $!\n"; exit 1;}
defined($cgrp) or (system('pod2text',$0), exit 1);
if (!(-e $cgrp)){print STDERR "$cgrp not found: $!\n"; exit 1;}
defined($rgrp) or (system('pod2text',$0), exit 1);
if (!(-e $rgrp)){print STDERR "$rgrp not found: $!\n"; exit 1;}
defined($map) or (system('pod2text',$0), exit 1);
if (!(-e $map)){print STDERR "$map not found: $!\n"; exit 1;}
defined($names) or (system('pod2text',$0), exit 1);
if (!(-e $names)){print STDERR "$names not found: $!\n"; exit 1;}

my ($rec,$ctr,$j,$k,$l,@temp,@temp1,@temp2,%sgrp_data,%cgrp_data,%rgrp_data,
$user_t,$system_t,$cuser_t,$csystem_t,@sgrp_common,@sgrp_100_80,@sgrp_80_60,
@sgrp_60_40,%abbv2name,%name2abbv,%temphash,%annot);

# open files
$ctr=0;

unless(open(INSGRP,"<$sgrp")){print "not able to open $sgrp\n\n";exit 1;}
unless(open(INCGRP,"<$cgrp")){print "not able to open $cgrp\n\n";exit 1;}
unless(open(INRGRP,"<$rgrp")){print "not able to open $rgrp\n\n";exit 1;}
unless(open(INMAP,"<$map")){print "not able to open $map\n\n";exit 1;}
unless(open(INNAMES,"<$names")){print "not able to open $names\n\n";exit 1;}
unless(open(OUTC,">Validated_clusters_common.xls")){print "not able to open Validated_clusters_common.xls\n";exit 1;}
unless(open(OUTCG,">Validated_clusters_common.txt")){print "not able to open Validated_clusters_common.txt\n";exit 1;}
unless(open(OUTC80,">Validated_clusters_80.xls")){print "not able to open Validated_clusters_80.xls\n";exit 1;}
unless(open(OUTC80G,">Validated_clusters_80.txt")){print "not able to open Validated_clusters_80.txt\n";exit 1;}
unless(open(OUTC60,">Validated_clusters_60.xls")){print "not able to open Validated_clusters_60.xls\n";exit 1;}
unless(open(OUTC60G,">Validated_clusters_60.txt")){print "not able to open Validated_clusters_60.txt\n";exit 1;}
unless(open(OUTC40,">Validated_clusters_40.xls")){print "not able to open Validated_clusters_40.xls\n";exit 1;}
unless(open(OUTC40G,">Validated_clusters_40.txt")){print "not able to open Validated_clusters_40.txt\n";exit 1;}

# read in files
#LAS_subclade1001: retli|86356474 retli|86361060 retli|86356879 retli|86356911 rleguminosarum|116249307 
while($rec=<INSGRP>){
	my $err_string = validate_protein_gi_number($rec);
	if ($err_string) {
		chomp $rec;
		print STDERR "$rec proteins are not in refseq format\n";
		print STDERR "$err_string\n";
		print STDERR "Exiting....\n";
		exit 1;
	}
	$sgrp_data{&mk_key($rec)}=$rec;
}
close(INSGRP);
print STDERR scalar(keys %sgrp_data)," records read from $sgrp ..\n";

while($rec=<INCGRP>){
	my $err_string = validate_protein_gi_number($rec);
	if ($err_string) {
		chomp $rec;
		print STDERR "$rec proteins are not in refseq format\n";
		print STDERR "$err_string\n";
		print STDERR "Exiting....\n";
		exit 1;
	}
	$cgrp_data{&mk_key($rec)}=$rec;
}
close(INCGRP);
print STDERR scalar(keys %cgrp_data)," records read from $cgrp ..\n";

while($rec=<INRGRP>){
	my $err_string = validate_protein_gi_number($rec);
	if ($err_string) {
		chomp $rec;
		print STDERR "$rec proteins are not in refseq format\n";
		print STDERR "$err_string\n";
		print STDERR "Exiting....\n";
		exit 1;
	}	
	$rgrp_data{&mk_key($rec)}=$rec;
}
close(INRGRP);
print STDERR scalar(keys %rgrp_data)," records read from $rgrp .. \n";

while($rec=<INMAP>){
	if($rec =~ /#/){next;}
	@temp=split("\t",$rec);
	chomp $temp[1];
	$abbv2name{$temp[0]}=$temp[1];
	$name2abbv{$temp[1]}=$temp[0];
#	$orgs{$temp[0]}=0;
}
close (INMAP);

while($rec=<INNAMES>){
	if($rec =~ /#/){next;}
	my $err_string = validate_refseq_format($rec);
	if ($err_string) {
		chomp $rec;
		print STDERR "$rec is not in refseq format\n";
		print STDERR "$err_string\n";
		print STDERR "Exiting....\n";
		exit 1;
	}
	@temp=split(/\|/,$rec);
	chomp $temp[4];
	$i="$temp[1]\t$temp[3]\t$temp[4]";
	my $gi = $temp[1];
	die "Protein with GI number $gi already exists\n" if exists $annot{$temp[1]};
	$annot{$temp[1]}=$i;
}
close (INNAMES);

while (($i,$j) = each (%sgrp_data)){
	my($flag,$key,$val,@sarr,@carr,$lc,$max_ovlap); $flag=0;
	#compare to consv grps
	if(exists $cgrp_data{$i}){$flag=1;}
	else{
		@sarr=split(/_/,$i);
		while(($key,$val) = each %cgrp_data){
			@carr=split(/_/,$key);
			#compare only if @carr is a smaller cluster, consv params will produce same or smaller clusters
			if(scalar(@carr) <= scalar(@sarr)){
				$lc=List::Compare->new ('--unsorted', '--accelerated',\@sarr,\@carr);
				if(!defined($max_ovlap)){
					$max_ovlap=scalar($lc->get_intersection());
					if($max_ovlap >= (0.8 * scalar(@sarr))){$flag=2;}
					elsif($max_ovlap >= (0.6 * scalar(@sarr))){$flag=3;}
					elsif($max_ovlap >= (0.4 * scalar(@sarr))){$flag=4;}
				}
				#compare to the larger cluster that has greater overlap
				else{
					$k=scalar($lc->get_intersection());
					if($k > $max_ovlap){
						$max_ovlap=$k;
						if($max_ovlap >= (0.8 * scalar(@sarr))){$flag=2;}
						elsif($max_ovlap >= (0.6 * scalar(@sarr))){$flag=3;}
						elsif($max_ovlap >= (0.4 * scalar(@sarr))){$flag=4;}
					}	
				}	
			}
		}
	}
	undef($max_ovlap);
	#compare to relx grps
	if(exists $rgrp_data{$i} && $flag==1){push @sgrp_common,$j;}
	else{
		@sarr=split(/_/,$i);
		while(($key,$val) = each %cgrp_data){
			@carr=split(/_/,$key);
			#compare only if @carr is a bigger cluster, relx params will produce bigger clusters
			if(scalar (@carr) => scalar(@sarr)){
				$lc=List::Compare->new ('--unsorted', '--accelerated',\@sarr,\@carr);
				if(!defined($max_ovlap)){
					$max_ovlap=scalar($lc->get_intersection());
					if(($max_ovlap >= (0.8 * scalar(@sarr))) && ($flag==2)){push @sgrp_100_80,$j;}
					elsif(($max_ovlap >= (0.6 * scalar(@sarr))) && ($flag==3)){push @sgrp_80_60,$j;}
					elsif(($max_ovlap >= (0.4 * scalar(@sarr))) && ($flag==4)){push @sgrp_60_40,$j;}
				}
				#compare to the larger cluster that has greater overlap
				else{
					$k=scalar($lc->get_intersection());
					if($k > $max_ovlap){
						$max_ovlap=$k;
						if(($max_ovlap >= (0.8 * scalar(@sarr))) && ($flag==2)){push @sgrp_100_80,$j;}
						elsif(($max_ovlap >= (0.6 * scalar(@sarr))) && ($flag==3)){push @sgrp_80_60,$j;}
						elsif(($max_ovlap >= (0.4 * scalar(@sarr))) && ($flag==4)){push @sgrp_60_40,$j;}
					}	
				}
			}
		}
	}
}

print OUTC "\tCOMMON CLUSTERS\n\n";
print OUTC "\tName\tNof clusters\n";
print OUTC "Selected clusters file:\t$sgrp\t",scalar(keys %sgrp_data),"\n";
print OUTC "Conserved clusters file:\t$cgrp\t",scalar(keys %cgrp_data),"\n";
print OUTC "Relaxed clusters file:\t$rgrp\t",scalar(keys %rgrp_data),"\n\n";
print OUTC "Count:\t",scalar(@sgrp_common),"\n\n";
foreach $i (@sgrp_common){
	print OUTCG $i;
	@temp=split(' ',$i);
	$temp[0]=~ s/\:$//;
	print OUTC "Cluster:\t$temp[0]\nOrganism\tGI\tAccession\tDescription\n";
	for $j (1..(scalar(@temp)-1)){
		@temp1=split(/\|/,$temp[$j]);
		@temp2=split("\t",$annot{$temp1[1]});
		$temp2[2]=~ s/\[[\s\S]+\]$//;#remove org name
		print OUTC $abbv2name{$temp1[0]},"\t",$temp2[0],"\t",$temp2[1],"\t",$temp2[2],"\n";
	}
	print OUTC "\n";
}
close(OUTC);
close(OUTCG);

print OUTC80 "\t80\%\+ CONSERVED CLUSTERS\n\n";
print OUTC80 "\tName\tNof clusters\n";
print OUTC80 "Selected clusters file:\t$sgrp\t",scalar(keys %sgrp_data),"\n";
print OUTC80 "Conserved clusters file:\t$cgrp\t",scalar(keys %cgrp_data),"\n";
print OUTC80 "Relaxed clusters file:\t$rgrp\t",scalar(keys %rgrp_data),"\n\n";
print OUTC80 "Count:\t",scalar(@sgrp_100_80),"\n\n";
foreach $i (@sgrp_100_80){
	print OUTC80G $i;
	@temp=split(' ',$i);
	$temp[0]=~ s/\:$//;
	print OUTC80 "Cluster:\t$temp[0]\nOrganism\tGI\tAccession\tDescription\n";
	for $j (1..(scalar(@temp)-1)){
		@temp1=split(/\|/,$temp[$j]);
		@temp2=split("\t",$annot{$temp1[1]});
		$temp2[2]=~ s/\[[\s\S]+\]$//;#remove org name
		print OUTC80 $abbv2name{$temp1[0]},"\t",$temp2[0],"\t",$temp2[1],"\t",$temp2[2],"\n";
	}
	print OUTC80 "\n";
}
close(OUTC80);
close(OUTC80G);

print OUTC60 "\t60\%\-80\% CONSERVED CLUSTERS\n\n";
print OUTC60 "\tName\tNof clusters\n";
print OUTC60 "Selected clusters file:\t$sgrp\t",scalar(keys %sgrp_data),"\n";
print OUTC60 "Conserved clusters file:\t$cgrp\t",scalar(keys %cgrp_data),"\n";
print OUTC60 "Relaxed clusters file:\t$rgrp\t",scalar(keys %rgrp_data),"\n\n";
print OUTC60 "Count:\t",scalar(@sgrp_80_60),"\n\n";
foreach $i (@sgrp_80_60){
	print OUTC60G $i;
	@temp=split(' ',$i);
	$temp[0]=~ s/\:$//;
	print OUTC60 "Cluster:\t$temp[0]\nOrganism\tGI\tAccession\tDescription\n";
	for $j (1..(scalar(@temp)-1)){
		@temp1=split(/\|/,$temp[$j]);
		@temp2=split("\t",$annot{$temp1[1]});
		$temp2[2]=~ s/\[[\s\S]+\]$//;#remove org name
		print OUTC60 $abbv2name{$temp1[0]},"\t",$temp2[0],"\t",$temp2[1],"\t",$temp2[2],"\n";
	}
	print OUTC60 "\n";
}
close(OUTC60);
close(OUTC60G);

print OUTC40 "\t40\%\-60\% CONSERVED CLUSTERS\n\n";
print OUTC40 "\tName\tNof clusters\n";
print OUTC40 "Selected clusters file:\t$sgrp\t",scalar(keys %sgrp_data),"\n";
print OUTC40 "Conserved clusters file:\t$cgrp\t",scalar(keys %cgrp_data),"\n";
print OUTC40 "Relaxed clusters file:\t$rgrp\t",scalar(keys %rgrp_data),"\n\n";
print OUTC40 "Count:\t",scalar(@sgrp_60_40),"\n\n";
foreach $i (@sgrp_60_40){
	print OUTC40G $i;
	@temp=split(' ',$i);
	$temp[0]=~ s/\:$//;
	print OUTC40 "Cluster:\t$temp[0]\nOrganism\tGI\tAccession\tDescription\n";
	for $j (1..(scalar(@temp)-1)){
		@temp1=split(/\|/,$temp[$j]);
		@temp2=split("\t",$annot{$temp1[1]});
		$temp2[2]=~ s/\[[\s\S]+\]$//;#remove org name
		print OUTC40 $abbv2name{$temp1[0]},"\t",$temp2[0],"\t",$temp2[1],"\t",$temp2[2],"\n";
	}
	print OUTC40 "\n";
}
close(OUTC40);
close(OUTC40G);

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print  STDERR "\nRuntime details \n";
print  STDERR "System time for process: $system_t\n";
print  STDERR "User time for process: $user_t\n";
exit;
