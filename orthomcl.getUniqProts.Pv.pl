#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Dec 07, 2012

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use lib '/home/surya/bin/modules';
use SS;
eval {
	require Bio::DB::Fasta;
	require List::Compare;
	require Bio::SeqIO;
};

=head1 NAME

orthomcl.getUniqProts.Pv.pl - Set operations on Nohit and Weakhit Pvirid proteins 

=head1 SYNOPSIS

  % orthomcl.getUniqProts.Pv.pl --fas prot.fas --sfaa nohitweak1.faa --ffaa nohitweak1.faa --pfaa nohitweak1.faa --oname 
  
=head1 DESCRIPTION

This script reads in 3 protein faa files of Pviridiflave Nohits and Weakhits proteins with names in refseq 
format. 1 faa for each pan proteiome. Outputs set operation stats and faa subsets.

-Version
  1. Pvirid prots to Psyringae, Pfluorescens and Pputida OMCL annotated (core/shared/lin spc) pan-proteomes

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --vfaa    <.faa>    Refseq formatted protein fasta file (required)
   --sfaa   <.faa>     FAA file of Nohits and weak proteins Pv proteome blast to pan proteome PS(required)
   --ffaa   <.faa>     FAA file of Nohits and weak proteins Pv proteome blast to pan proteome PF(required)
   --pfaa   <.faa>     FAA file of Nohits and weak proteins Pv proteome blast to pan proteome PP(required)
   --oname  <>         Prefix for output FAA
      
=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut

sub mk_id {
	# Params: seq name >gi|5294|RAST|XXX| hypothetical protein [Pvirid4GB draft]
	# returns id after gi
	if ($_[0] =~ /^>gi\|(\d+)/) {return $1;}
	else {return;}
}

sub getGIs(){
	# Params: faa filename
	my (@sprot,$i,@temp,$ctr,$infaa);
	$ctr=0;
	$infaa = Bio::SeqIO->new('-file' => "<$_[0]",'-format' => 'fasta' );
	while ($i = $infaa->next_seq()) {
		#>gi|5294|RAST|XXX| hypothetical protein [Pvirid4GB draft]
		@temp=(); @temp=split(/\|/,$i->primary_id());
		push @sprot, $temp[1];
		$ctr++;
	}
	return @sprot;
}


my ($i,$faa,$sfaa,$ffaa,$pfaa,$oname);

GetOptions (
	'vfaa=s' => \$faa,
	'sfaa=s' => \$sfaa,
	'ffaa=s' => \$ffaa,
	'pfaa=s' => \$pfaa,
	'oname:s' => \$oname) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($faa) or (system('pod2text',$0), exit 1);
if (!(-e $faa)){print STDERR "$faa not found: $!\n"; exit 1;}
defined($sfaa) or (system('pod2text',$0), exit 1);
if (!(-e $sfaa)){print STDERR "$sfaa not found: $!\n"; exit 1;}
defined($ffaa) or (system('pod2text',$0), exit 1);
if (!(-e $ffaa)){print STDERR "$ffaa not found: $!\n"; exit 1;}
defined($pfaa) or (system('pod2text',$0), exit 1);
if (!(-e $pfaa)){print STDERR "$pfaa not found: $!\n"; exit 1;}
$oname ||= 'Pvirid';

my ($rec,$ctr,$j,$k,$l,$m,@temp,@temp1,$flag,$db,@srcgi,@sgi,@fgi,@pgi,$lc,@temp2,@temp3);

#get gi's from input FAA files
@srcgi=&getGIs($faa); #GIs from nohit/weak prots for PV
@sgi=&getGIs($sfaa); #GIs from nohit/weak prots for panprot PS
@fgi=&getGIs($ffaa); #GIs from nohit/weak prots for panprot PF
@pgi=&getGIs($pfaa); #GIs from nohit/weak prots for panprot PP
print STDERR $faa,': ',scalar @srcgi,"\n";
print STDERR "Pvirid protein counts for no/weak hits for each pan-proteome:\n";
print STDERR $sfaa,': ',scalar @sgi,"\n";
print STDERR $ffaa,': ',scalar @fgi,"\n";
print STDERR $pfaa,': ',scalar @pgi,"\n\n";

# making DB of the faa file, forcing reindex
$db = Bio::DB::Fasta->new($faa, '-makeid'=>\&mk_id, '-reindex'=>1);
if(!($db)){ die "These was a problem creating database from $faa\n";}

# Get subsets of Pv prots
# Writing out subsets of Pvirid prots

# All nohits and weak hits
$lc=List::Compare->new ('--unsorted', '--accelerated',\@sgi,\@fgi,\@pgi);
@temp=(); @temp=$lc->get_intersection();
print STDERR "Pvirid wtout hits or only weak hits to PSy, PFl and PPu: ",scalar @temp,"\n\n";
unless(open(OUTUF,">${oname}\.notPSy.notPFl.notPPu.faa")){print "not able to open ${oname}\.notPSy.notPFl.notPPu.faa\n\n";exit 1;}
foreach $i (@temp){
	$j=$db->header($i);
	if(!($j)){print STDERR "$i not found in $faa index\n"; exit;}
	print OUTUF "\>$j\n";
	print OUTUF $db->seq($i),"\n";
}
close(OUTUF);

# PFL comparison
$lc=List::Compare->new ('--unsorted', '--accelerated',\@srcgi,\@fgi);
@temp=(); @temp=$lc->get_unique; # get prots only in Pv, i.e. PF prots that had matches to PV
print STDERR "Pvirid prots with matches to prots in PFl: ",scalar @temp,"\n";

$lc=List::Compare->new ('--unsorted', '--accelerated',\@temp,\@sgi);
@temp1=(); @temp1=$lc->get_intersection;
print STDERR "Pvirid prots that match PFl but not PSy: ",scalar @temp1,"\n\n";
unless(open(OUTUF,">${oname}\.matchPFl.notPSy.faa")){print "not able to open ${oname}\.matchPFl.notPSy.faa\n\n";exit 1;}
foreach $i (@temp1){
	$j=$db->header($i);
	if(!($j)){print STDERR "$i not found in $faa index\n"; exit;}
	print OUTUF "\>$j\n";
	print OUTUF $db->seq($i),"\n";
}
close(OUTUF);

# PPU comparison
$lc=List::Compare->new ('--unsorted', '--accelerated',\@srcgi,\@pgi);
@temp=(); @temp=$lc->get_unique; # get prots only in Pv, i.e. PP prots that had matches to PV
print STDERR "Pvirid prots with matches to prots in Ppu: ",scalar @temp,"\n";

$lc=List::Compare->new ('--unsorted', '--accelerated',\@temp,\@sgi);
@temp2=(); @temp2=$lc->get_intersection;
print STDERR "Pvirid prots that match PPu but not PSy: ",scalar @temp2,"\n\n";
unless(open(OUTUF,">${oname}\.matchPPu.notPSy.faa")){print "not able to open ${oname}\.matchPPu.notPSy.faa\n\n";exit 1;}
foreach $i (@temp2){
	$j=$db->header($i);
	if(!($j)){print STDERR "$i not found in $faa index\n"; exit;}
	print OUTUF "\>$j\n";
	print OUTUF $db->seq($i),"\n";
}
close(OUTUF);

$lc=List::Compare->new ('--unsorted', '--accelerated',\@temp1,\@temp2);
@temp3=(); @temp3=$lc->get_intersection;
print STDERR "Pvirid prots that match both PFl and PPu but not PSy: ",scalar @temp3,"\n\n";
unless(open(OUTUF,">${oname}\.matchPFl.matchPPu.notPSy.faa")){print "not able to open ${oname}\.matchPFl.matchPPu.notPSy.faa\n\n";exit 1;}
foreach $i (@temp3){
	$j=$db->header($i);
	if(!($j)){print STDERR "$i not found in $faa index\n"; exit;}
	print OUTUF "\>$j\n";
	print OUTUF $db->seq($i),"\n";
}
close(OUTUF);

# Pv prot shared with Psyr but not present in Pfl or Pp
$lc=List::Compare->new ('--unsorted', '--accelerated',\@srcgi,\@sgi);
@temp=(); @temp=$lc->get_unique();# Pv prots found in Psy
print STDERR "Pvirid prots that match PSy: ",scalar @temp,"\n";
$lc=List::Compare->new ('--unsorted', '--accelerated',\@fgi,\@pgi);
@temp1=(); @temp1=$lc->get_union();# Pv prots not in Pfl or Ppu
$lc=List::Compare->new ('--unsorted', '--accelerated',\@temp,\@temp1); 
@temp2=(); @temp2=$lc->get_intersection();# Pv prots in Psy but not in Pfl or Ppu
print STDERR "Pvirid prots that match PSy but not PFl or PPu: ",scalar @temp2,"\n\n";
unless(open(OUTUF,">${oname}\.matchPSy.notPFl.notPPu.faa")){print "not able to open ${oname}\.matchPSy.notPFl.notPPu.faa\n\n";exit 1;}
foreach $i (@temp2){
	$j=$db->header($i);
	if(!($j)){print STDERR "$i not found in $faa index\n"; exit;}
	print OUTUF "\>$j\n";
	print OUTUF $db->seq($i),"\n";
}
close(OUTUF);

# Pv prot shared with Pf but not present in Psyr or Pp
$lc=List::Compare->new ('--unsorted', '--accelerated',\@srcgi,\@fgi);
@temp=(); @temp=$lc->get_unique();# Pv prots found in Pfl
print STDERR "Pvirid prots that match PFl: ",scalar @temp,"\n";
$lc=List::Compare->new ('--unsorted', '--accelerated',\@sgi,\@pgi);
@temp1=(); @temp1=$lc->get_union();# Pv prots not in Psy or Ppu
$lc=List::Compare->new ('--unsorted', '--accelerated',\@temp,\@temp1); 
@temp2=(); @temp2=$lc->get_intersection();# Pv prots in Pfl but not in Psy or Ppu
print STDERR "Pvirid prots that match PFl but not PSy or PPu: ",scalar @temp2,"\n\n";
unless(open(OUTUF,">${oname}\.matchPFl.notPSy.notPPu.faa")){print "not able to open ${oname}\.matchPFl.notPSy.notPPu.faa\n\n";exit 1;}
foreach $i (@temp2){
	$j=$db->header($i);
	if(!($j)){print STDERR "$i not found in $faa index\n"; exit;}
	print OUTUF "\>$j\n";
	print OUTUF $db->seq($i),"\n";
}
close(OUTUF);

# Pv prot shared with Pp but not present in Psyr or Pfl
$lc=List::Compare->new ('--unsorted', '--accelerated',\@srcgi,\@pgi);
@temp=(); @temp=$lc->get_unique();# Pv prots found in Ppu
print STDERR "Pvirid prots that match PPu: ",scalar @temp,"\n";
$lc=List::Compare->new ('--unsorted', '--accelerated',\@sgi,\@fgi);
@temp1=(); @temp1=$lc->get_union();# Pv prots not in Psy or Pfl
$lc=List::Compare->new ('--unsorted', '--accelerated',\@temp,\@temp1); 
@temp2=(); @temp2=$lc->get_intersection();# Pv prots in Ppu but not in Psy or Pfl
print STDERR "Pvirid prots that match PPu but not PSy or PFl: ",scalar @temp2,"\n\n";
unless(open(OUTUF,">${oname}\.matchPPu.notPFl.notPSy.faa")){print "not able to open ${oname}\.matchPPu.notPFl.notPSy.faa\n\n";exit 1;}
foreach $i (@temp2){
	$j=$db->header($i);
	if(!($j)){print STDERR "$i not found in $faa index\n"; exit;}
	print OUTUF "\>$j\n";
	print OUTUF $db->seq($i),"\n";
}
close(OUTUF);

#deleting index file
unlink "${faa}.index"; 

exit;