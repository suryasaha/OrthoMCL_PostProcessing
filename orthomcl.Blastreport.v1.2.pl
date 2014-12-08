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

 orthomcl.blastreport.v1.2.pl - prepare a Excel report to summarize protein hits to labelled pan-proteome 

=head1 SYNOPSIS

  % orthomcl.Blastreport.v1.2.pl --report blast.out --cutoff 10 --map map.txt --groups group.txt --out report.xls
  
=head1 DESCRIPTION

 Reads in a blast report presuming a default e-value cutoff of 10. Print out a tab separated file with
 description of the hits for each source protein. HITS WITH EVAL < CUTOFF ARE NOT RECORDED IN EITHER OF
 THE OUTPUT FILES. Also lists the OrthoMCL groups hit by a protein.
 
 May not work if draft FAA headers not formatted like Refseq ">gi|1|art|XXX| undefined product [wACP3 draft]"
 
=head1 VERSION HISTORY
 Version   1.0  INPUT: blast report, Mapping file 
                OUTPUT: .XLS file with best hits/genome, .XLS file with rest of the hits and 
                        proteins with no hits              
 Version   1.1  INPUT: blast report, Mapping file, Groups file 
                OUTPUT: .XLS file with best hits/genome, .XLS file with rest of the hits and 
                        proteins with no hits, Hit groups listed in .XLS files
 Version   1.2  INPUT: blast report, Mapping file, Groups file 
                OUTPUT: .XLS file with best hits/genome, .XLS file with rest of the hits and 
                        proteins with no hits, Hit groups listed in .XLS files, Weak hit proteins
                        (cutoff>weak eval>best hit)do NOT overlap with hits with e-value below cutoff. 
                        Overlapping proteins printed separately now in multi file. 

=head1 NOTE
 >>>>>>DO NOT CHANGE THE EXCEL FILE LAYOUT OR BLAST_CLASSES WILL BREAK!!<<<<<<
 >>>>>>DO NOT SAVE AFTER VIEWING OUTPUT FILES IN EXCEL/CALC!!<<<<<<
  
=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --report <.out>    Blast report in text format with protein clusters from OrthoMCL run to validate (required)
   --cutoff <float>   A float value for blast Evalue cutoff <10>
   --map    <.txt>    Tab delimited file with abbreviations followed by full genome names used in Refseq (required)
   --groups <.txt>    Text file with protein clusters from OrthoMCL run (required)
   --out    <text>    Excel report file
   --num    <INT>     Number of organisms/genomes (required)
   --src    <text>    Abbv for assembly (required)
      
=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut

sub getGI{
	my (@temp); @temp=split(/\|/,$_[0]);
	return $temp[1];
}

sub getOrg{
	my (@temp); @temp=split(/\|/,$_[0]);
	$temp[0]=~ s/[\s\S]+\[//;
	$temp[0]=~ s/\]$//; #remove desc prefix
	return $temp[0];
}

sub getClass{
	my (@temp); @temp=split(/\|/,$_[0]);
	return $temp[1];
}

sub getDesc{
	my (@temp); @temp=split(/\|/,$_[0]);
	$temp[0]=~ s/\[[\s\S]+\]//; #remove org name suffix
#	if ( !defined $temp[0]){ print STDERR "*******\n".$_[0]; exit 1;}
	return $temp[0];
}

my ($rep,$cutoff,$map,$out,$in,$result,$hit,@temp,%orgs_cnt,@orgs,@hits,$nohits,$onlyweakhits,$src,
$weakhits,$src_annot_ctr,$src_unannot_ctr,$i,$j,$rec,%ctr,$grp,%groups,%hitgrps,%nohitgrps,$flag,$gen);

GetOptions (
	'report=s' => \$rep,
	'cutoff:f' => \$cutoff,
	'map=s'    => \$map,
	'groups=s' => \$grp,
	'out:s'    => \$out,
	'num=i' => \$gen,
	'src=s' => \$src) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($rep) or (system('pod2text',$0), exit 1);
if (!(-e $rep)){print STDERR "$rep not found: $!\n"; exit 1;}
defined($map) or (system('pod2text',$0), exit 1);
if (!(-e $map)){print STDERR "$map not found: $!\n"; exit 1;}
defined($grp) or (system('pod2text',$0), exit 1);
if (!(-e $grp)){print STDERR "$grp not found: $!\n"; exit 1;}
$cutoff ||=10;
$out ||="${rep}\.xls";
defined($gen) or (system('pod2text',$0), exit 1);
defined($src) or (system('pod2text',$0), exit 1);

unless(open(XLS,">$out")){print "not able to open $out\n\n";exit 1;}
unless(open(MXLS,">multi\.$out")){print "not able to open multi\.$out\n\n";exit 1;}
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
	$nohitgrps{$temp[0]}='';
	#keys are GI's, values are group ID's
	for $i (1..$#temp){ 
		$temp[$i]=~ s/[\S]+\|//;
		$groups{$temp[$i]}=$temp[0];
	}
}
close (INGRP);

print XLS "\t\t\tBLAST report\n\n";
print MXLS "\t\t\tBLAST report (multiple hits)\n\n";
print XLS "$src protein\t\t";
foreach $i (@orgs){
	print XLS $i,"\t\t";
}
print XLS "\nGI\tDesc\t";
for $i (1..$gen){print XLS "Class\tEval\t";}
print XLS "Groups\n";
print MXLS "$src protein\t\tHit\nGI\tDesc\tOrg\tClass\tEval\tGroup\n";
$in = new Bio::SearchIO(-format => 'blast', -file   => $rep);
$onlyweakhits=$weakhits=$nohits='';
$src_annot_ctr=$src_unannot_ctr=0;
while($result = $in->next_result ) {
	## $result is a Bio::Search::Result::ResultI compliant object
	#init data structures
	@hits=(); %ctr=%orgs_cnt;
	#get hit data
	if($result->num_hits>0){
		$src_annot_ctr++; $flag=0;
		while($hit = $result->next_hit ) {
	    	## $hit is a Bio::Search::Hit::HitI compliant object
	    	##debug
#	    	if ($hit->description eq ''){print STDERR "\n".$hit->name; exit 1;}
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
		    	$ctr{&getOrg($hit->description)}++; $flag=1;
	    	}
	    	elsif(($flag==0) && ($hit->significance > $cutoff)){#for protein with only poor hits
	    		$onlyweakhits=$onlyweakhits.&getGI($result->query_name)."\t".&getDesc($result->query_description)."\t".
	    			&getOrg($hit->description)."\t".&getClass($hit->description)."\t".$hit->significance."\n";
	    	}
	    	elsif ($hit->significance > $cutoff){#for protein with poor hits
	    		##debug
#	    		if ($result->query_description eq ''){print STDERR "11111111\n".$result->query_name."\n".$hit->name;}
	    		##to handle weird prokka protein descriptions that break parser
	    		if ($result->query_description ne ''){
	    			$weakhits=$weakhits.&getGI($result->query_name)."\t".&getDesc($result->query_description)."\t".
	    			&getOrg($hit->description)."\t".&getClass($hit->description)."\t".$hit->significance."\n";
	    		}
	    		else{
	    			$weakhits=$weakhits.&getGI($result->query_name)."\tNA via parser\t".
	    			&getOrg($hit->description)."\t".&getClass($hit->description)."\t".$hit->significance."\n";
	    		}
	    	}
		}
		#print hit data
		%hitgrps=();
		if (@hits>0){ 
			if ($result->query_description ne ''){
				print XLS "\n",&getGI($result->query_name),"\t",&getDesc($result->query_description),"\t";
			}
			else{
				print XLS "\n",&getGI($result->query_name),"\tNA via parser\t";
			}
		}
		foreach $i (@orgs){
			 if($ctr{$i}==0){#if no hits
			 	print XLS "\t\t";
			 }
			 elsif($ctr{$i}==1){#if 1 hit
			 	foreach $j (@hits){
			 		if($j->{ORG} eq $i){
			 			print XLS $j->{CLASS},"\t",$j->{EVAL},"\t";
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
			 	print XLS $class,"\t",$eval,"\t";#print best hit
			 	#record group ID
			 	if(exists $groups{$gi}){ $hitgrps{$groups{$gi}}=1;}#in case the hit is protein not present in OMCL group
			 	foreach $j (@hits){# print other hits
			 		if (($j->{ORG} eq $i) && ($j->{EVAL} != $eval)){
				    		##debug
#				    		if ($result->query_description eq ''){print STDERR "222222\n".$result->query_name;}
				    		##to handle weird prokka protein descriptions that break parser			 			
				    		if ($result->query_description ne ''){
				 			print MXLS &getGI($result->query_name),"\t",&getDesc($result->query_description),"\t";				    		
				    		}
				    		else{
				 			print MXLS &getGI($result->query_name),"\tNA via parser\t";
				    		}

			 			print MXLS $i,"\t",$j->{CLASS},"\t",$j->{EVAL};
			 			#in case the hit is a protein not present in any OMCL group
			 			if(exists $groups{$j->{GI}}){print MXLS "\t",$groups{$j->{GI}},"\n";}
			 			else {print MXLS "\n";}
			 		}
			 	}
			 }
		}
#		if (@hits>0) {print XLS "\n";}
		#print group ids
		@temp= keys %hitgrps;
		foreach $i (@temp){
			delete $nohitgrps{$i}; print XLS $i,"\t";
		}
	}
	else{
		$src_unannot_ctr++;
		$nohits=$nohits.&getGI($result->query_name)."\t".&getDesc($result->query_description)."\n";
	}
}

print MXLS "\n\n$src proteins with weak hits\n";
print MXLS "$src protein\t\tHit\nGI\tDesc\tOrg\tClass\tEval\n";
print MXLS $weakhits;

print MXLS "\n\n$src proteins with only weak hits\n";
print MXLS "$src protein\t\tHit\nGI\tDesc\tOrg\tClass\tEval\n";
print MXLS $onlyweakhits;
print MXLS "\n\n$src proteins with no hits\n";
print MXLS "GI\tDesc\n";
print MXLS $nohits;
print MXLS "\n\nGroups with no hits from $src proteins\n";
@temp= keys %nohitgrps;
foreach $i (@temp){ print MXLS "$i\n";}

close(XLS);
close(MXLS);

print STDERR "\n\n$src_annot_ctr $src proteins had hits (above and below threshold). $src_unannot_ctr had no hits\n";

exit;
