#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha June 13, 2013

use strict;
use warnings;
use Getopt::Long;

=head1 NAME

 orthomcl.convert2Refseq.pl - Convert protein fasta and GFF files to Refseq like format for OrthomclAdjustFasta

=head1 SYNOPSIS

  % orthomcl.convert2Refseq.pl --fas file.fas --seed 100
  
=head1 DESCRIPTION

 Reads in Fasta file produced by ACT or Artemis. Formats the sequence headers to look like Refseq headers.
 
=head1 VERSION HISTORY
 Version   1.0  INPUT : Fasta 
                OUTPUT: Formatted Fasta files
                                    
=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --fas      <.fas>    Fasta file of protein sequences (required)
   --organism <name>    Strain name (Pvirid, etc)
   --seed     <INT>     Number to use as staring ID
      
=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut



my($fas,$name,$seed,$i,$j,$rec,@temp,%idx);
GetOptions (
	'fas=s' => \$fas,
	'organism=s' => \$name,
	'seed:i'    => \$seed) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($fas) or (system('pod2text',$0), exit 1);
if (!(-e $fas)){print STDERR "$fas not found: $!\n"; exit 1;}
defined($name) or (system('pod2text',$0), exit 1);
$seed ||= 1;

unless(open(FAS,"<$fas")){print "not able to open $fas\n\n";exit 1;}
unless(open(OFAS,">formatted\.$fas")){print "not able to open formatted\.$fas\n\n";exit 1;}
#Seq names
#>CDS CDS FIG01208100: hypothetical protein 108:221 reverse MW:4035
#>CDS CDS hypothetical protein 205:348 forward MW:5087
#>CDS CDS Enoyl-[acyl-carrier-protein] reductase [NADH] (EC 1.3.1.9) 419:1195 reverse MW:28221
#>CDS CDS hypothetical protein 1521:1745 forward MW:8931
#>CDS CDS hypothetical protein 1802:2002 forward MW:7831
#>CDS CDS FIG01207832: hypothetical protein 2445:3341 reverse MW:33233
#
#>gi|190570479|ref|YP_001974837.1| chromosomal replication initiator protein DnaA [Wolbachia endosymbiont of Culex quinquefasciatus Pel]
#>gi|190570480|ref|YP_001974838.1| protein-export membrane protein SecF [Wolbachia endosymbiont of Culex quinquefasciatus Pel]

while ($rec=<FAS>){
	if($rec=~ /^>/){
#		$rec=~ s/^>//;
#		$rec=~ s/^>CDS CDS //;
#		@temp=split(' ',$rec);#splitting on space
#		#taking third last element
#		$idx{$temp[$#temp-2]}=$seed;
#		print OFAS "\>gi\|$seed\|RAST\|XXX\| ";
		@temp=split(/\|/,$rec);#splitting on |
		print OFAS "\>gi\|$seed\|"; 
		$seed++;
#		for $i (2..($#temp-3)){ $temp[$i]=~ s/\://g; print OFAS $temp[$i],' ';}
		$temp[0] =~ s/^>//;
		chomp $temp[1];
		print OFAS $temp[0].'|';
		print OFAS $temp[1];
		print OFAS " \[$name draft\]\n";
	}
	else{print OFAS $rec;}
}
close(OFAS);

exit;
