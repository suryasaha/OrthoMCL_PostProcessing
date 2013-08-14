#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Sep 27, 2010

use strict;
use warnings;
use Getopt::Long;
use POSIX;
eval {
	require Bio::SearchIO;
};
use Bio::SearchIO; 

=head1 NAME

 orthomcl.blastreport.v1.pl - prepare a Excel report to summarize <LAS> protein hits to labelled <Rhizobiacae> pan-proteome 

=head1 SYNOPSIS

  % orthomcl.Blastreport.v1.pl --report blast.out --cutoff 10 --map map.txt --groups group.txt --out report.xls
  
=head1 DESCRIPTION

 Reads in a blast report presuming a default e-value cutoff of 10. Print out a tab separated file with
 description of the hits for each LAS protein. HITS WITH EVAL < CUTOFF ARE NOT RECORDED IN EITHER OF THE 
 OUTPUT FILES. Also lists the OrthoMCL groups hits by a <LAS> protein.
 
=head1 TODO
 The weakhits set overlaps with other sets right now. Will need to list only those proteins for which all
 hits are weak.
 
=head1 VERSION HISTORY
 Version   1.0  INPUT: blast report, Mapping file 
                OUTPUT: .XLS file with best hits/genome, .XLS file with rest of the hits and 
                        proteins with no hits              
 Version   1.1  INPUT: blast report, Mapping file, Groups file 
                OUTPUT: .XLS file with best hits/genome, .XLS file with rest of the hits and 
                        proteins with no hits, Hit groups listed in .XLS files
=head1 NOTE
 >>>>>>DO NOT CHANGE THE EXCEL FILE LAYOUT OR BLAST_CLASSES WILL BREAK!!<<<<<<
 >>>>>>DO NOT SAVE AFTER VIEWING OUTPUT FILES IN EXCEL/CALC!!<<<<<<
  
=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --report <.out>    Blast report in text formatText file with protein clusters from OrthoMCL run to validate (required)
   --cutoff <float>   A float value <10>
   --map    <.txt>    Tab delimited file with abbreviations followed by full genome names used in Refseq (required)
   --groups <.txt>    Text file with protein clusters from OrthoMCL run (required)
   --out    <text>    Excel report file
      
=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut

sub getGI{
	$_[0]=~ /^gi\|(\d+)/; return $1;
}

sub getOrg{
	my (@temp); @temp=split(/\|/,$_[0]);
	$temp[0]=~ s/[\s\S]+\[//;
	$temp[0]=~ s/\]$//;
	return $temp[0];
}

sub getClass{
	my (@temp); @temp=split(/\|/,$_[0]);
	return $temp[1];
}

sub getDesc{
	my (@temp); @temp=split(/\|/,$_[0]);
	$temp[0]=~ s/\[[\s\S]+\]//;
	return $temp[0];
}

my ($rep,$cutoff,$map,$out,$in,$result,$hit,@temp,%orgs_cnt,@orgs,@hits,$nohits,$weakhits,
$las_annot_ctr,$las_unannot_ctr,$i,$j,$rec,%ctr,$grp,%groups,%hitgrps);

GetOptions (
	'report=s' => \$rep,
	'cutoff:f' => \$cutoff,
	'map=s'    => \$map,
	'groups=s' => \$grp,
	'out:s'    => \$out) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($rep) or (system('pod2text',$0), exit 1);
if (!(-e $rep)){print STDERR "$rep not found: $!\n"; exit 1;}
defined($map) or (system('pod2text',$0), exit 1);
if (!(-e $map)){print STDERR "$map not found: $!\n"; exit 1;}
defined($grp) or (system('pod2text',$0), exit 1);
if (!(-e $grp)){print STDERR "$grp not found: $!\n"; exit 1;}
$cutoff ||=10;
$out ||="${rep}\.xls";

unless(open(OUT,">$out")){print "not able to open $out\n\n";exit 1;}
unless(open(OUT1,">multi\.$out")){print "not able to open multi\.$out\n\n";exit 1;}
unless(open(INMAP,"<$map")){print "not able to open $map\n\n";exit 1;}
unless(open(INGRP,"<$grp")){print "not able to open $grp\n\n";exit 1;}

while($rec=<INMAP>){
	if($rec =~ /#/){next;}
	@temp=split("\t",$rec);
	chomp $temp[1];
	push @orgs,$temp[1];
	$orgs_cnt{$temp[1]}=0;
}
close (INMAP);

while($rec=<INGRP>){
	if($rec =~ /#/){next;}
	@temp=split(' ',$rec);
	$temp[0]=~ s/://;
	#keys are GI's, values are group ID's
	for $i (1..$#temp){ 
		$temp[$i]=~ s/[\S]+\|//;
		$groups{$temp[$i]}=$temp[0];
	}
}
close (INGRP);

print OUT "\t\t\tBLAST report\n\n";
print OUT1 "\t\t\tBLAST report (multiple hits)\n\n";
print OUT "LAS protein\t\t";
foreach $i (@orgs){
	print OUT $i,"\t\t";
}
print OUT "\nGI\tDesc\tClass\tEval\tClass\tEval\tClass\tEval\tClass\tEval\tClass\tEval\tClass\tEval";
print OUT "\tClass\tEval\tClass\tEval\tGroups\n";
print OUT1 "LAS protein\t\tHit\nGI\tDesc\tOrg\tClass\tEval\tGroup\n";
$in = new Bio::SearchIO(-format => 'blast', -file   => $rep);
$weakhits=$nohits='';
$las_annot_ctr=$las_unannot_ctr=0;
while($result = $in->next_result ) {
	## $result is a Bio::Search::Result::ResultI compliant object
	#init data structures
	@hits=(); %ctr=%orgs_cnt;
	#get hit data
	if($result->num_hits>0){
		$las_annot_ctr++;
		while($hit = $result->next_hit ) {
	    ## $hit is a Bio::Search::Hit::HitI compliant object
	    	if ($hit->significance <= $cutoff){
		    	$rec={
		    		GI    => &getGI($hit->name),
		    		DESC  => &getDesc($hit->description),
		    		CLASS => &getClass($hit->description),
		    		SCORE => $hit->bits,
		    		EVAL  => $hit->significance,
		    		ORG   => &getOrg($hit->description),
		    	};
		    	push @hits,$rec;
		    	$ctr{&getOrg($hit->description)}++;
	    	}
	    	else{#for poor hits
	    		$weakhits=$weakhits.&getGI($result->query_name)."\t".&getDesc($result->query_description)."\t".
	    			&getOrg($hit->description)."\t".&getClass($hit->description)."\t".$hit->significance."\n";
	    	}
		}
		#print hit data
		%hitgrps=();
		if (@hits>0){ print OUT "\n",&getGI($result->query_name),"\t",&getDesc($result->query_description),"\t";}
		foreach $i (@orgs){
			 if($ctr{$i}==0){#no hits
			 	print OUT "\t\t";
			 }
			 elsif($ctr{$i}==1){#1 hit
			 	foreach $j (@hits){
			 		if($j->{ORG} eq $i){
			 			print OUT $j->{CLASS},"\t",$j->{EVAL},"\t";
			 			#in case the hit is a protein not present in any OMCL group
			 			if(exists $groups{$j->{GI}}){$hitgrps{$groups{$j->{GI}}}=1;} 
			 			last;
			 		}
			 	}
			 }
			 if($ctr{$i}>1){#if >1 hits for this genome, print weaker hit to aux XLS file
			 	my ($class,$eval,$gi); $eval=$cutoff;
			 	foreach $j (@hits){
			 		if (($j->{ORG} eq $i) && ($j->{EVAL} < $eval)){
			 			$class=$j->{CLASS}; $eval=$j->{EVAL}; $gi=$j->{GI};
			 		}
			 	}
			 	print OUT $class,"\t",$eval,"\t";#print best hit
			 	#record group ID
			 	if(exists $groups{$gi}){ $hitgrps{$groups{$gi}}=1;}#in case the hit is protein not present in OMCL group
			 	foreach $j (@hits){# print other hits
			 		if (($j->{ORG} eq $i) && ($j->{EVAL} != $eval)){
			 			print OUT1 &getGI($result->query_name),"\t",&getDesc($result->query_description),"\t";
			 			print OUT1 $i,"\t",$j->{CLASS},"\t",$j->{EVAL};
			 			#in case the hit is a protein not present in any OMCL group
			 			if(exists $groups{$j->{GI}}){print OUT1 "\t",$groups{$j->{GI}},"\n";}
			 			else {print OUT1 "\n";}
			 		}
			 	}
			 	
			 }
		}
#		if (@hits>0) {print OUT "\n";}
		#print group ids
		@temp= keys %hitgrps;
		foreach $i (@temp){ print OUT $i,"\t"}
	}
	else{
		$las_unannot_ctr++;
		$nohits=$nohits.&getGI($result->query_name)."\t".&getDesc($result->query_description)."\n";
	}
}

print OUT1 "\n\nLAS proteins with weak hits\n";
print OUT1 "LAS protein\t\tHit\nGI\tDesc\tOrg\tClass\tEval\n";
print OUT1 $weakhits;
print OUT1 "\n\nLAS proteins with no hits\n";
print OUT1 "GI\tDesc\n";
print OUT1 $nohits;

close(OUT);
close(OUT1);

print STDERR "$las_annot_ctr LAS proteins had hits (above and below threshold). $las_unannot_ctr had no hits\n";

exit;