#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Aug 23, 2010

use strict;
use warnings;

my $locallib = eval{
			require local::lib;
			local::lib->import();
			1;
		};
if ($locallib){
	use local::lib;
}


use Getopt::Long;
use POSIX;
use Cwd;
eval {
	require Bio::DB::Fasta;
	require Bio::SeqIO;
	require Switch;
	require List::Compare;
};

use Bio::SeqIO;
use Bio::DB::Fasta;
use Switch;
use List::Compare;

=head1 NAME

 orthomcl.Proteome2sets.pl - Read groups file produced by OrthoMCL and write out protein sets 

=head1 SYNOPSIS

  % orthomcl.Proteome2sets.pl --selgroups groups.txt --consparam groups.txt --relparam groups.txt --faa allgenomes.faa --map map.txt --desc proteome_sets.xls
  
=head1 DESCRIPTION`

 Reads in groups file produced by OrthoMCL. Outputs protein sets for the core genome (present in all 
 genomes), shared but not core (present in some but not all genomes), and lineage specific (present 
 in only one genome). Use 'cat *.faa > allgenomes.faa' to get all the proteins into a single 
 file. Put all abbreviations and full genome names in a tab delimited text file. Verifies core 
 by requiring a core cluster to have 1 rep in the cluster from all strains under relaxed params 
 and 75%+ genomes under strict params.

=head1 VERSION HISTORY

 Version   1.0  INPUT: Groups file, FAA protein file, Mapping file 
                OUTPUT: Description .XLS file
 Version   2.0  INPUT: 3 Group files, FAA protein file, Mapping file 
                OUTPUT: Description .XLS file, FAA protein file for pan-proteome with labels,
                        FAA protein files for each core and shared cluster
             
=head1 NOTES
 Produces errors for getcwd() and line 282 since warnings are enabled. Redirect errors to error file or 
 null device (2>Err or 2>/dev/null)
	 Subroutine main::getcwd redefined at /usr/share/perl/5.10/Exporter.pm line 67.
	 at ./orthomcl.Proteome2sets.pl line 11
	 Use of uninitialized value $\ in regexp compilation at ./orthomcl.Proteome2sets.pl line 282.

=head1 TODO

=head1 NOTE
Make sure that all gi numbers are unique in the Fasta headers. If you use orthomcl.convert2RefseqFasta.pl to format the fasta headers, then make sure you use a seed value higher than existing values.

=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --selgroups <.txt>    Text file with protein clusters from selected OrthoMCL run (required)
   --consparam <.txt>    Text file with protein clusters from OrthoMCL run with most conservative parameters (required)
   --relparam  <.txt>    Text file with protein clusters from OrthoMCL run with most relaxed parameters (required)
   --faa       <.faa>    Concatenated FAA file of proteins from all genomes. Make sure its in the current directory (required)
   --map       <.txt>)   Tab delimited file with abbreviations followed by full genome names used in Refseq (required)
   --desc      <.xls>    OUTPUT tab delimited file with descriptions of proteins in each group for each genome
      
=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut

sub mk_key{
	my (@temp,@temp1,@temp2,$i,$key);
	@temp=split(' ',$_[0]);
	@temp2=();
	for $i (1..(scalar(@temp)-1)){
		@temp1=split(/\|/,$temp[$i]);#can be shortened skipping temp arrays
		push @temp2,$temp1[1];#record GIs
	}
	@temp = sort {$a <=> $b} @temp2;
	foreach $i (@temp){#make key of sorted GIs
		if(!defined($key)){$key=$i;}
		else{$key=$key.'_'.$i;}
	}
	return $key;
}

sub mk_id {
	if ($_[0] =~ /^>gi\|(\d+)/) {return $1;}
	else {return;}
}

sub mk_dir{
my ($name,$i);
$name=$_[0]; $i=localtime();
if (-e $name){rename $name,"${name}_$i"; unlink glob "${name}_/* $(name}/.*";rmdir ($name);}
mkdir ($name, 0755) or warn "Cannot make $name directory: $!\n";
}

my ($i,$sgrp,$cgrp,$rgrp,$faa,$desc,$map);

GetOptions (
	'selgroups=s' => \$sgrp,
	'consgroups=s' => \$cgrp,
	'relgroups=s' => \$rgrp,
	'faa=s' => \$faa,
	'map=s' => \$map,
	'desc:s'=> \$desc) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($sgrp) or (system('pod2text',$0), exit 1);
if (!(-e $sgrp)){print STDERR "$sgrp not found: $!\n"; exit 1;}
defined($cgrp) or (system('pod2text',$0), exit 1);
if (!(-e $cgrp)){print STDERR "$cgrp not found: $!\n"; exit 1;}
defined($rgrp) or (system('pod2text',$0), exit 1);
if (!(-e $rgrp)){print STDERR "$rgrp not found: $!\n"; exit 1;}
defined($faa) or (system('pod2text',$0), exit 1);
if (!(-e $faa)){print STDERR "$faa not found: $!\n"; exit 1;}
defined($map) or (system('pod2text',$0), exit 1);
if (!(-e $map)){print STDERR "$map not found: $!\n"; exit 1;}
$i=$sgrp;
$i=~ s/.txt//;
$desc ||= "$i.proteome_sets.xls";

my ($rec,$ctr,@groups,$infaa,%abbv2name,%name2abbv,%allprots,%orgs,$j,$k,$l,@temp,
%cgrp_data,%rgrp_data,$outfaa);
unless(open(INSGRP,"<$sgrp")){print "not able to open $sgrp\n\n";exit 1;}
unless(open(INCGRP,"<$cgrp")){print "not able to open $cgrp\n\n";exit 1;}
unless(open(INRGRP,"<$rgrp")){print "not able to open $rgrp\n\n";exit 1;}
unless(open(INMAP,$map)){print "not able to open $map\n\n";exit 1;}
unless(open(DESC,">$desc")){print "not able to open $desc\n";exit 1;}

# read in files
$ctr=0;
while($rec=<INSGRP>){
	if($rec =~ /#/){next;}
	push @groups, [split(' ',$rec)]; $ctr++;
}
close (INSGRP);

while($rec=<INCGRP>){$cgrp_data{&mk_key($rec)}=$rec;}
close(INCGRP);

while($rec=<INRGRP>){$rgrp_data{&mk_key($rec)}=$rec;}
close(INRGRP);

print DESC "Selected clusters file:\t$sgrp\t",$ctr,"\n";
print DESC "Conserved clusters file:\t$cgrp\t",scalar(keys %cgrp_data),"\n";
print DESC "Relaxed clusters file:\t$rgrp\t",scalar(keys %rgrp_data),"\n\n";

while($rec=<INMAP>){
	if($rec =~ /#/){next;}
	@temp=split("\t",$rec);
	chomp $temp[1];
	$abbv2name{$temp[0]}=$temp[1];
	$name2abbv{$temp[1]}=$temp[0];
	$orgs{$temp[0]}=0;
}
close (INMAP);

$infaa = Bio::SeqIO->new('-file' => "<$faa",'-format' => 'fasta' );
#init prot hash with entire pan-proteome
while ($i = $infaa->next_seq()) {
	#>gi|16262454|ref|NP_435247.1| FdoG formate dehydrogenase-O alpha subunit [Sinorhizobium meliloti 1021]
	@temp=split(/\|/,$i->id());
	$j=$i->desc();
	$j=~ s/[\s\S]+\[//;
	$j=~ s/\]$//;
	if (!exists $name2abbv{$j}){ print STDERR "\n\n$j not found in \%name2abbv\n\n"; exit 1;}
	my $gi = $temp[1];
	die "Protein with GI number $gi already exists\n" if exists $allprots{"$name2abbv{$j}\|$gi"};
	$allprots{"$name2abbv{$j}\|$gi"}='l';#mark all as lineage specific to start with
}
$infaa->close();

#find whether protein is core, shared or unique
foreach $i (@groups){
	#chk if all orgs are present, if yes, mark all prots as core, else all as shared
	#LAS_subclade2467: retli|86358613 smeliloti|15966232 rleguminosarum|116250056
	#LAS_subclade2595: retli|86360090 retli|86361198
	my (%temp_orgs);
	%temp_orgs=%orgs;
	for $j (1..(scalar (@$i)-1)){#for each member
		@temp=split(/\|/,$i->[$j]);
		$temp_orgs{$temp[0]}++;#mark if org is present in cluster
	}
	$ctr=0;
	@temp= values (%temp_orgs);
	foreach $j (@temp){if($j>0){$ctr++;}}#count number of orgs in cluster
	
	if($ctr == keys (%temp_orgs)){#mark all prots in cluster as core {c,cp,uc,ucp} if all orgs present
		#verify predicted core proteins using consv and relx run data
		my($flag,$key,$val,@sarr,@carr,$lc,$max_ovlap); $flag=0;
		@sarr=();
		for $j (1..(scalar(@$i)-1)){
			@temp=split(/\|/,$i->[$j]);#can be shortened skipping temp arrays
			push @sarr,$temp[1];#record GIs
		}
	
		#present with > 75% members in consv run
		while(($key,$val) = each %cgrp_data){
			@carr=split(/_/,$key);
			if(scalar(@carr) >= (0.75*scalar(@sarr))){
				$lc=List::Compare->new ('--unsorted', '--accelerated',\@sarr,\@carr);
				if(!defined($max_ovlap)){$max_ovlap=scalar($lc->get_intersection());}
				else{
					$k=scalar($lc->get_intersection());
					if($k > $max_ovlap){$max_ovlap=$k;}	
				}
			}
		}
		if($max_ovlap < (0.75*scalar(@sarr))){$flag=1}
		else{# qualifies the > 75% rule, checking in relx run
			undef($max_ovlap);
			#present with 100% members in relx run
			while(($key,$val) = each %rgrp_data){
				@carr=split(/_/,$key);
				if(scalar(@carr) >= (scalar(@sarr))){
					$lc=List::Compare->new ('--unsorted', '--accelerated',\@sarr,\@carr);
					if(!defined($max_ovlap)){$max_ovlap=scalar($lc->get_intersection());}
					else{
						$k=scalar($lc->get_intersection());
						if($k > $max_ovlap){$max_ovlap=$k;}
					}
				}
			}
		}
		if($max_ovlap < (scalar(@sarr))){$flag=1;}
		
		for $j (1..(scalar (@$i)-1)){
			@temp=split(/\|/,$i->[$j]);
			if($temp_orgs{$temp[0]}>1 && $flag==0){$allprots{$i->[$j]}='cp';}#core+paralog
			elsif($temp_orgs{$temp[0]}==1 && $flag==0){$allprots{$i->[$j]}='c';}#core
			elsif($temp_orgs{$temp[0]}>1 && $flag==1){$allprots{$i->[$j]}='ucp';}#putative but unverified core+paralog
			elsif($temp_orgs{$temp[0]}==1 && $flag==1){$allprots{$i->[$j]}='uc';}#putative but unverified core
		}
	}
	#mark all prots in cluster as shared {s,sp} if < all orgs present
	#Do we add a criterion of at least > x% of orgs required in a 
	#cluster before genes are labelled shared instead of just > 1 org
	elsif(($ctr < (keys (%temp_orgs))) && ($ctr > 1)){#mark all prots in cluster as shared if < all orgs present
		for $j (1..(scalar (@$i)-1)){
			@temp=split(/\|/,$i->[$j]);
			if($temp_orgs{$temp[0]}>1){$allprots{$i->[$j]}='sp';}#shared+paralog
			elsif($temp_orgs{$temp[0]}==1){$allprots{$i->[$j]}='s';}#shared
		}
	}
	elsif($ctr == 1){#mark all prots in cluster as lineage specific and paralogs
		for $j (1..(scalar (@$i)-1)){ $allprots{$i->[$j]}='lp';}
	}	
}


#write out labelled pan-proteome
$infaa = Bio::SeqIO->new('-file' => "<$faa",'-format' => 'fasta' );
$outfaa = Bio::SeqIO->new(-file => ">labelled\.$faa" , '-format' => 'Fasta');
while ($i = $infaa->next_seq()) {
	@temp=split(/\|/,$i->id());
	$j=$i->desc(); 	$j=~ s/[\s\S]+\[//;	$j=~ s/\]$//;
	$k=$i->desc();
	switch($allprots{"$name2abbv{$j}\|$temp[1]"}){
		case 'c' { $k=$k.'| core';}
		case 'cp' { $k=$k.'| core_paralogous';}
		case 'uc' { $k=$k.'| unverified_core';}
		case 'ucp' { $k=$k.'| unverified_core_paralogous';}
		case 's' { $k=$k.'| shared';}
		case 'sp' { $k=$k.'| shared_paralogous';}
		case 'l' { $k=$k.'| lineage_specific';}
		case 'lp' { $k=$k.'| lineage_specific_paralogous';}	
	}
	$i->desc($k); $outfaa->write_seq($i);
}
$infaa->close();
$outfaa->close();

#write out cluster FAA files with labels for individual proteins
&mk_dir('core_class');
&mk_dir('shared_class');

# making DB of the fasta file
my ($db,$cwd);  
#$db = Bio::DB::Fasta->new($faa, '-makeid' =>\&mk_id);
$db = Bio::DB::Fasta->new($faa, '-makeid'=>\&mk_id, '-reindex'=>1);
if(!($db)){ die "These was a problem creating database from $faa\n";}
$cwd=getcwd();
foreach $i (@groups){
	#LAS_subclade2467: retli|86358613 smeliloti|15966232 rleguminosarum|116250056
	#LAS_subclade2595: retli|86360090 retli|86361198
	$i->[0]=~ s/$\://;
	if($allprots{$i->[1]}=~ /c/){#core class{c,cp,uc,ucp}, if first is core class, then rest will be the same
		unless(open(TEMP,">${cwd}\/core_class\/$i->[0]\.faa")){print "not able to open ${cwd}\/core_class\/$i->[0]\.faa\n\n";exit 1;}
		for $j (1..(scalar (@$i)-1)){#for each member
			@temp=(); @temp=split(/\|/,$i->[$j]);
			#write with label
			$k=$db->header($temp[1]); 
			if(!($k)){ print STDERR "could not find $i->[$j] in core\n\n"; next;}
			chomp $k;
			switch($allprots{$i->[$j]}){
				case 'c' { $k=$k.'| core';}
				case 'cp' { $k=$k.'| core_paralogous';}
				case 'uc' { $k=$k.'| unverified_core';}
				case 'ucp' { $k=$k.'| unverified_core_paralogous';}
			}
			print TEMP "\>$k\n";
			print TEMP $db->seq($temp[1]),"\n";
		}
		close(TEMP);
	}
	elsif($allprots{$i->[1]}=~ /s/){#shared class{s,sp}, if first is shared class, then rest will be the same
		unless(open(TEMP,">${cwd}\/shared_class\/$i->[0]\.faa")){print "not able to open ${cwd}\/shared_class\/$i->[0]\.faa\n\n";exit 1;}
		for $j (1..(scalar (@$i)-1)){#for each member
			@temp=split(/\|/,$i->[$j]);
			#write with label
			$k=$db->header($temp[1]);
			if(!($k)){ print STDERR "could not find $i->[$j] in shared\n\n"; next;}
			chomp $k;
			switch($allprots{$i->[$j]}){
				case 's' { $k=$k.'| shared';}
				case 'sp' { $k=$k.'| shared_paralogous';}
			}
			print TEMP "\>$k\n";
			print TEMP $db->seq($temp[1]),"\n";
		}
		close(TEMP);
	}
}

#collect summary info for each category for each genome
my (%c,%cp,%uc,%ucp,%s,%sp,%l,%lp,%tot,%files);
%c=%cp=%uc=%ucp=%s=%sp=%l=%lp=%tot=%orgs;
while (($i,$j) = each %allprots){
	@temp=split(/\|/,$i);
	$tot{$temp[0]}++;
	switch($j){
		case 'c' { $c{$temp[0]}++;}
		case 'cp' { $cp{$temp[0]}++;}
		case 'uc' { $uc{$temp[0]}++;}
		case 'ucp' { $ucp{$temp[0]}++;}
		case 's' { $s{$temp[0]}++;}
		case 'sp' { $sp{$temp[0]}++;}
		case 'l' { $l{$temp[0]}++;}
		case 'lp' { $lp{$temp[0]}++;}
	}
}

print DESC "\n\nOrganism\tTotal proteins\tCore\tUnverified Core\tCore Paralogous\tUnverified Core Paralogous\tTotal Core\t\%\tShared\tShared Paralogous\tTotal Shared";
print DESC "\t\%\tLineage specific\tLineage specific and Paralogous\tTotal Lineage specific\t\%\n";
foreach $i (values %name2abbv){
	print DESC "$abbv2name{$i}\t$tot{$i}\t$c{$i}\t$uc{$i}\t$cp{$i}\t$ucp{$i}\t",$c{$i}+$cp{$i}+$uc{$i}+$ucp{$i},"\t",sprintf("%.4f",((($c{$i}+$cp{$i}+$uc{$i}+$ucp{$i})/$tot{$i})*100));
	print DESC "\t$s{$i}\t$sp{$i}\t",$s{$i}+$sp{$i},"\t",sprintf("%.4f",((($s{$i}+$sp{$i})/$tot{$i})*100)),"\t$l{$i}\t$lp{$i}\t";
	print DESC $l{$i}+$lp{$i},"\t",sprintf("%.4f",((($l{$i}+$lp{$i})/$tot{$i})*100)),"\n";
}

exit;

