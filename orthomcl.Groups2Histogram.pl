#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Aug 19, 2010

use strict;
use warnings;
use Getopt::Long;
use POSIX;

=head1 NAME

 Groups2Histogram.pl - Read groups file produced by OrthoMCL and write out histogram png  

=head1 SYNOPSIS

  % Groups2Histogram.pl --group groups.txt --names names.txt --out histogram.R --desc groups.desc.xls
  
=head1 DESCRIPTION

 Outputs R code and calls R to create a histogram. Use 'cat *.faa| fgrep '>' > names.txt'
to get all the protein names into a single file. This tool only works for Refseq style
FASTA headers.

=head1 VERSION HISTORY

 Version   1.0  INPUT: Groups file, Fasta headers 
                OUTPUT: Histogram, stem chart and description .XLS file
                NOTES: None 
              
=head1 TODO

 1. FAA file of each cluster
 2. Annotation of cluster/proteins by COG number using PTT file from RefSeq
 3. Membership profile of each cluster like Supp file from Stanhope ppr
 4. Use Map file to report % of proteins in a organism participating in clusters

=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --group  <.txt>    Text file with protein clusters (required)
   --names  <.txt>    Text file with protein FASTA headers (required)
   --out    <.R>      R script to draw histogram
   --desc   <.xls>    Tab delimited file with descriptions of each group
      
=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut

my ($i,$group,$name,$desc,$out,$cmd);

GetOptions (
	'group=s' => \$group,
	'names=s' => \$name,
	'out:s' => \$out,
	'desc:s'=> \$desc) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($group) or (system('pod2text',$0), exit 1);
if (!(-e $group)){print STDERR "$group not found: $!\n"; exit 1;}
defined($name) or (system('pod2text',$0), exit 1);
if (!(-e $name)){print STDERR "$name not found: $!\n"; exit 1;}
$i=$group;
$i=~ s/.txt//;
$out ||= "$i.histogram.R";
$desc ||= "$i.desc.xls";

my ($rec,$ctr,@groups,%names,%tot_prot_ctr,%clus_prot_ctr,$j,$k,$l,@temp);
unless(open(INGROUP,$group)){print "not able to open $group\n\n";exit 1;}
unless(open(INNAMES,$name)){print "not able to open $name\n\n";exit 1;}
unless(open(R,">$out")){print "not able to open $out\n";exit 1;}
unless(open(DESC,">$desc")){print "not able to open $desc\n";exit 1;}

# read in files
$ctr=0;
while($rec=<INGROUP>){
	if($rec =~ /#/){next;}
	push @groups, [split(' ',$rec)]; $ctr++;
}
close (INGROUP);
print DESC "Total clusters:\t$ctr\n";

while($rec=<INNAMES>){
	if($rec =~ /#/){next;}
	@temp=split(/\|/,$rec);
	chomp $temp[4];
	$i="$temp[1]\t$temp[3]\t$temp[4]";
	$names{$temp[1]}=$i;
	$temp[4]=~ s/[\s\S]+\[//;
	$temp[4]=~ s/\]$//;
	if(exists $tot_prot_ctr{$temp[4]}){$tot_prot_ctr{$temp[4]}++}
	else {$tot_prot_ctr{$temp[4]}=1}
	@temp=();
}
close (INNAMES);
print DESC "\nTotal proteins in each genome\n";
while (($i,$j) = each %tot_prot_ctr){ print DESC "$i\t$j\n"}
print DESC "\n";

foreach $i (@groups){
	for $j (1..(scalar (@$i)-1)){#for each member
		@temp=split(/\|/,$i->[$j]);
		if(exists $clus_prot_ctr{$temp[0]}){$clus_prot_ctr{$temp[0]}++}
		else {$clus_prot_ctr{$temp[0]}=1}
	}
}
print DESC "Proteins participating in clusters from each genome\n";
while (($i,$j) = each %clus_prot_ctr){ print DESC "$i\t$j\n"}
print DESC "\n";

$cmd='OrthoMCL_clusters<-c(';

foreach $i (@groups){
#	print STDERR scalar(@$i),' ';
	$cmd=$cmd.(scalar(@$i)-1).',';
	$i->[0]=~ s/\:$//;
	print DESC 'Cluster name: ',$i->[0],"\n";
	for $j (1..(scalar (@$i)-1)){#for each member
		@temp=split(/\|/,$i->[$j]);
		print DESC $i->[$j],"\t",$names{$temp[1]},"\n";
	}
	print DESC "\n";
}
close (DESC);

$cmd=~ s/,$//; $cmd =$cmd."\)\nhist(OrthoMCL_clusters, xlab\=\"Number of proteins in a cluster\"\, main\=\"OrthoMCL cluster profile\"\)\n";
print R $cmd;print R "stem\(OrthoMCL_clusters\)\n"; close (R);

system("R CMD BATCH $out");

exit;



