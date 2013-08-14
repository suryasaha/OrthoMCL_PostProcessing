#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha June 13, 2013

use strict;
use warnings;
use Getopt::Long;

=head1 NAME

 orthomcl.convert2Refseq.pl - Convert protein fasta and GFF files to Refseq like format for OrthomclAdjustFasta

=head1 SYNOPSIS

  % orthomcl.convert2Refseq.pl --fas file.fas --gff file.gff --organism Pvirid --seed 100
  
=head1 DESCRIPTION

 Reads in Fasta and GFF file produced by ACT or Artemis. Formats the sequence headers to look like Refseq headers.
 
=head1 VERSION HISTORY
 Version   1.0  INPUT: Protein GFF and Fasta 
                OUTPUT: Formatted GFF and Fasta files
                                    
=head1 NOTES
 No 1 to 1 correspondence between prots in FAS and GFF files. Both processed independently. Change??
 
=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --fas      <.fas>    Fasta file of protein sequences (required)
   --gff      <.gff>    Protein GFF file (required)
   --organism <name>    Strain name (Pvirid, etc)
   --seed     <INT>     Number to use as staring ID
      
=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut



my($fas,$gff,$name,$seed,$i,$j,$rec,@temp,%idx);
GetOptions (
	'fas=s' => \$fas,
	'gff=s' => \$gff,
	'organism=s' => \$name,
	'seed:i'    => \$seed) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($fas) or (system('pod2text',$0), exit 1);
if (!(-e $fas)){print STDERR "$fas not found: $!\n"; exit 1;}
defined($name) or (system('pod2text',$0), exit 1);
defined($gff) or (system('pod2text',$0), exit 1);
if (!(-e $gff)){print STDERR "$gff not found: $!\n"; exit 1;}
$seed ||= 1;

unless(open(FAS,"<$fas")){print "not able to open $fas\n\n";exit 1;}
unless(open(GFF,"<$gff")){print "not able to open $gff\n\n";exit 1;}
unless(open(OFAS,">formatted\.$fas")){print "not able to open formatted\.$fas\n\n";exit 1;}
unless(open(OGFF,">formatted\.$gff")){print "not able to open formatted\.$gff\n\n";exit 1;}
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
		@temp=split(' ',$rec);#splitting on space
		#taking third last element
		$idx{$temp[$#temp-2]}=$seed;
		print OFAS "\>gi\|$seed\|RAST\|XXX\| "; $seed++;
		for $i (2..($#temp-3)){ $temp[$i]=~ s/\://g; print OFAS $temp[$i],' ';}
		print OFAS "\[$name draft\]\n";
	}
	else{print OFAS $rec;}
}
close(OFAS);

#GFF file
#gff_seqname     artemis exon    108     221     .       -       .       ID=CDS:complement%28108..221%29;translation=MKLPCKGKSLVPPDLGVALQARNPIVLVSVLLRLYSY;product=FIG01208100:+hypothetical+protein
#gff_seqname     artemis exon    205     348     .       +       .       ID=CDS:205..348;translation=MQGNFKSNPHAQSLNPNGRSAAANTCKIITLWLICEDLSKKCCH+SSS;product=hypothetical+protein
#gff_seqname     artemis exon    419     1195    .       -       .       ID=CDS:complement%28419..1195%29;db_xref=GO:0004318;translation=MILQGKKGLITGIINKRSIAYGIAKTLSEHGAELAITYQNEVIK+EKLLPIANELNVELTLHCDVSNKETIDSAFGKIEKKWDNLDFLVHAIAFSDKNELNGR+YVDTSLKNFLNAMHISCYSFTALAQRAEKMMLNGGSLLTLSYYGAEKVMPNYNVMGLC+KAALEASVKYIACDLGPQNIRVNAISAGPIRTLASSGISDFHSISEWNRSNSPLRRNV+TIEDVGKAALYLLSDLSSGTTGEILHVDSGYNVVGMKIVD;product=Enoyl-[acyl-carrier-protein]+reductase+[NADH]+%28EC+1.3.1.9%29;eC_number=1.3.1.9

#NC_010981.1	RefSeq	gene	152	1534	.	+	.	ID=NC_010981.1:dnaA;locus_tag=WPa_0001;old_locus_tag=WP0001;db_xref=GeneID:6384583
#NC_010981.1	RefSeq	CDS	152	1531	.	+	0	ID=NC_010981.1:dnaA:unknown_transcript_1;Parent=NC_010981.1:dnaA;locus_tag=WPa_0001;old_locus_tag=WP0001;note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication%3B can also affect transcription of multiple genes including itself.;transl_table=11;product=chromosomal replication initiation protein;protein_id=YP_001974837.1;db_xref=GI:190570479;db_xref=GeneID:6384583;exon_number=1
#NC_010981.1	RefSeq	start_codon	152	154	.	+	0	ID=NC_010981.1:dnaA:unknown_transcript_1;Parent=NC_010981.1:dnaA;locus_tag=WPa_0001;old_locus_tag=WP0001;note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication%3B can also affect transcription of multiple genes including itself.;transl_table=11;product=chromosomal replication initiation protein;protein_id=YP_001974837.1;db_xref=GI:190570479;db_xref=GeneID:6384583;exon_number=1
#NC_010981.1	RefSeq	stop_codon	1532	1534	.	+	0	ID=NC_010981.1:dnaA:unknown_transcript_1;Parent=NC_010981.1:dnaA;locus_tag=WPa_0001;old_locus_tag=WP0001;note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication%3B can also affect transcription of multiple genes including itself.;transl_table=11;product=chromosomal replication initiation protein;protein_id=YP_001974837.1;db_xref=GI:190570479;db_xref=GeneID:6384583;exon_number=1
#NC_010981.1	RefSeq	gene	1692	2567	.	+	.	ID=NC_010981.1:secF;locus_tag=WPa_0002;old_locus_tag=WP0002;db_xref=GeneID:6385550
#NC_010981.1	RefSeq	CDS	1692	2564	.	+	0	ID=NC_010981.1:secF:unknown_transcript_1;Parent=NC_010981.1:secF;locus_tag=WPa_0002;old_locus_tag=WP0002;note=forms a complex with SecD and YajC%3B SecDFyajC stimulates the proton motive force-driven protein translocation%3B seems to modulate the cycling of SecA by stabilizing its membrane-inserted state and appears to be required for the release of mature proteins from the extracytoplasmic side of the membrane%3B in some organisms%2C such as Bacillus subtilis%2C SecD is fused to SecF;transl_table=11;product=preprotein translocase subunit SecF;protein_id=YP_001974838.1;db_xref=GI:190570480;db_xref=GeneID:6385550;exon_number=1
#NC_010981.1	RefSeq	start_codon	1692	1694	.	+	0	ID=NC_010981.1:secF:unknown_transcript_1;Parent=NC_010981.1:secF;locus_tag=WPa_0002;old_locus_tag=WP0002;note=forms a complex with SecD and YajC%3B SecDFyajC stimulates the proton motive force-driven protein translocation%3B seems to modulate the cycling of SecA by stabilizing its membrane-inserted state and appears to be required for the release of mature proteins from the extracytoplasmic side of the membrane%3B in some organisms%2C such as Bacillus subtilis%2C SecD is fused to SecF;transl_table=11;product=preprotein translocase subunit SecF;protein_id=YP_001974838.1;db_xref=GI:190570480;db_xref=GeneID:6385550;exon_number=1
#NC_010981.1	RefSeq	stop_codon	2565	2567	.	+	0	ID=NC_010981.1:secF:unknown_transcript_1;Parent=NC_010981.1:secF;locus_tag=WPa_0002;old_locus_tag=WP0002;note=forms a complex with SecD and YajC%3B SecDFyajC stimulates the proton motive force-driven protein translocation%3B seems to modulate the cycling of SecA by stabilizing its membrane-inserted state and appears to be required for the release of mature proteins from the extracytoplasmic side of the membrane%3B in some organisms%2C such as Bacillus subtilis%2C SecD is fused to SecF;transl_table=11;product=preprotein translocase subunit SecF;protein_id=YP_001974838.1;db_xref=GI:190570480;db_xref=GeneID:6385550;exon_number=1

while($rec=<GFF>){
	@temp=split("\t",$rec);
	$temp[0]=$name; $temp[2]='CDS'; $temp[7]=0; chomp $temp[8];
	$temp[8]=$temp[8].';db_xref=GI:'.$idx{"$temp[3]\:$temp[4]"}."\n";
	for $i (0..7){print OGFF $temp[$i],"\t";}
	print OGFF $temp[8];
}
close (OGFF);

exit;