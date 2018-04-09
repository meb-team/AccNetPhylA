#!/usr/bin/perl

=head1 NAME

	AccNetPhylA.pl

=head1 SYNOPSIS

	perl AccNetPhylA.pl -in *.fasta -seq aa
	perl AccNetPhylA.pl [-in file.fasta] [-seq aa] [-cp '-c 0.9 -n 5 -g 1 -aS 0.9 -T 6 -M 5000'] [-out file1csv] 
	   [-tblout file2.csv] [-threshold 1] [-fast no] [-clean yes]


=head1 DESCRIPTION

	This Perl script AccNetPhylA (Accessory Constellation Network with Phylogenetic Analysis)
	takes up most of the analysis of AccNET (Accessory Constellation Network) 
	(Val F. Lanza et al., 2016) by replacing the clustering tool kClust (Hauser et al., 2013) 
	with CD-HIT (Weizhong Li et al., 2001) for analysis with nucleotide and protein sequences.
	
	"AccNET is a comparative genomic tool for accessory genome analysis using
	bipartite networks. The software has been designed to be compatible with
	most of the Network Analysis software (i.e. Cytoscape, Gephi or R)."
 

=head1 OPTIONS

	-Help|help|h, produces this help file. 
	
	-in, Files to analyze as input.	
	
	-seq, Proteic or nucleic sequence (nt|aa).

	-cp, CD-HIT parameters. /!\ Output filename no required. Default('-c 0.9 -n 5 -g 1 -aS 0.9 -T 6 -M 5000').
	
	-out, Network filename. Default(Network.csv).
	
	-tblout, Table filename. Default(Table.csv).
	
	-threshold,	Percent of genomes to consider coregenome (values > 1 includes coregenome to network). Default(1).
	
	-fast, Skip the phylogenetic distance determination. Default(no).
	
	-clean, Remove and class files. Default(yes).
	
	-dir, Directory's name. Default(Analyse).


=head1 AUTHORS

	SIGURET Clea


=head1 VERSION

	V1


=head1 DATE

	Creation : 7.02.2017
	Last modification : 15.06.2017

=cut

#Libraries
use List::Util; qw(first max maxstr min minstr reduce shuffle sum);
use strict;
use Getopt::Long;
use warnings;
use Pod::Usage;

#Scalars
my $help;	# Help flag
my $seqTY;
my $outFile;
my $outTable;
my $clean;
my $dir;
my $threshold;
my $fast;
my $cdHitParameters;
my $indexStrain;
my $indexProt;
my $tax;
my $output;
my $numClust;
my $seqID;	
my $seq;
my $finalOutput;
my $compte;
my $comptetot;
my $l;
my $ar;
my $k;
my $j;
my $i;
my $n;
my $tmp;
my $actual;
my $firstLine;
my $map;
my $file;
my $strain;
my $alnIn;
my $alnOut;
my $dist;
my $aln;
my $member;
my $tmpHeader;
my $tmpCluster;
my $id;
my $actualNumCluster;
my $verif;
my $input;


#Tables
my @text;
my @inFiles;	#Files to analyze
my @FichClust;
my @FichClustaln;
my @net;
my @allFasta;
my @cls;
my @headers;
my @salida;
my @inNet;
my @c;
my @c2;
my @dists;
my @txt;
my @val;
my @uniqs;
my @cSynt;
my @tmpArray;
my @clusterFasta;
my @Options;

#Hashes
my %taxasHash;
my %referenceHash;
my %sequences;
my %clusterSeq;
my %numSeqClust;
my %hashSeq;
my %Synt;
my %fasta;
my %indexHeader;
my %clusters;
my %clustersMembers;
my %suma;
my %repeat;
my %clusterTwin;
my %hashTwinGroup;
my %maxCls;
my %hashDiff;







############### SET COMMAND LINE OPTIONS #################################

@Options = (
		{OPT=>"Help|help|h",	VAR=>\$help,	DESC=>"Help flag"},
		{OPT=>"in=s{,}",	VAR=>\@inFiles,	DESC=>"Files to analyze"},
		{OPT=>"seq=s", VAR=>\$seqTY,	DESC=>"Proteic or nucleic sequence"},
		{OPT=>"out=s",	VAR=>\$outFile,	DEFAULT => 'Network.csv', DESC=>"Network filename"},
		{OPT=>"cp=s",	VAR=>\$cdHitParameters,	DEFAULT => '-c 0.9 -n 5 -g 1 -aS 0.9 -T 6 -M 5000', DESC=>"'CD-HIT parameters'"},
		{OPT=>"tblout=s",	VAR=>\$outTable,	DEFAULT => 'Table.csv', DESC=>"Table filename"},
		{OPT=>"threshold=s",	VAR=>\$threshold,	DEFAULT => '1', DESC=>"Percent of genomes to consider coregenome (values > 1 includes coregenome to network)"},
		{OPT=>"fast=s",	VAR=>\$fast,	DEFAULT => 'no', DESC=>"Skip the phylogenetic distance determination"},
		{OPT=>"clean=s",	VAR=>\$clean,	DEFAULT => 'yes', DESC=>"Remove and class files"},
		{OPT=>"dir=s",	VAR=>\$dir,	DEFAULT => 'Analyse', DESC=>"Directory's name"}
		
	);


#Check options and set variables
GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

# Now setup default values.
foreach (@Options) {
	if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
	${$_->{VAR}} = $_->{DEFAULT};
	}
}

unless(@inFiles){
	unless($help){print STDERR "You must specified the input files\n";exit;}
}

unless($seqTY){
	unless($help){print STDERR "You must specified the sequence type (proteic or nucleic)\n";exit;}
}

if ($help) {
		pod2usage(-verbose => 2, -noperldoc => 1);
		exit;
}



		
###########################  InputFile renaming
$indexStrain=1;
$input = scalar(@inFiles);
open(O,">equivalence.txt");
foreach $ar (@inFiles)
{
	if($ar =~ m/\/*([\w\d\.]+).f[a-z]+$|\/*([\w\d\.]+).f[a-z]+$|\/*([\w\d\.]+).f[a-z]+$/){$tax =$1;}
	#@c = split('\.',$ar);
	#$tax=$c[0];
	$taxasHash{$indexStrain}=$tax;
	$indexProt=0;
	
	
	open(A,$ar);
	@text = <A>;
	close A;
	
	foreach $l (@text)
	{
		if($l =~ /^>/)
		{
			$indexProt++;
			$referenceHash{$indexStrain}{$indexProt} = $l; ### Save sequence information
			$referenceHash{$indexStrain}{$indexProt} =~ s/>//;
			$sequences{$indexStrain}{$indexProt} = ">$indexStrain|$indexProt\n";  ##New header
			print O ">$indexStrain|$indexProt\t$l";
		}else{
			$sequences{$indexStrain}{$indexProt} .= $l;
		}
	}
	$indexStrain++;
}
close O;

open(O,">all_fasta.tmp");

foreach $k (keys(%sequences))
{
	foreach $j (keys(%{$sequences{$k}}))
	{
		print O "$sequences{$k}{$j}";
	}
}
close O;
######################## End InputFile renaming

########## Clustering and post-proccess
if($seqTY eq "aa"){
	system("cd-hit -i all_fasta.tmp -o rendu.tmp $cdHitParameters");
}
if($seqTY eq "nt"){
	system("cd-hit-est -i all_fasta.tmp -o rendu.tmp $cdHitParameters");
}
print "\n";


open(F1, "rendu.tmp.clstr") or die ("Problème ouverture rendu.tmp.clstr");

##### rendu.tmp.clstr ###############
#>Cluster 0							#
#0	1230aa, >1|14540... at 99.92%	#
#1	178aa, >1|20957... at 100.00%	# 
#2	2731aa, >1|4065... *			#
#3	595aa, >1|19313... at 99.66%	#
#####################################
#
#or
#
##### rendu.tmp.clstr ###############
#>Cluster 0							#
#0	1230nt, >1|14540... at -/99.92%	#
#1	178nt, >1|20957... at +/100.00%	# 
#2	2731nt, >1|4065... *			#
#3	595nt, >1|19313... at -/99.66%	#
#####################################

while (my $li = <F1>){
	chomp($li);
	if($li =~ m/^>Cluster (.+)/){
		$numClust = $1;
		++$comptetot;
	}
	if($li =~ m/^[0-9]+\t[0-9]+.., >([0-9]+\|[0-9]+)/){
		my $recup = $1;
		push (@{$clusterSeq{$numClust}},$recup);
		++$numSeqClust{$numClust};
	}
}
close F1;	


open(F2, "all_fasta.tmp") or die ("Problème ouverture all_fasta.tmp");
while (my $li = <F2>){## hashing fasta file ##
	chomp($li);
	if($li =~ m/^>([0-9]+\|[0-9]+)/){
		$seqID = $1;
		$actual = $li;
		$actual =~ s/\s//g;
	}
	else{
		$hashSeq{$seqID}=$li;
		#print"$seqID\n$li\n";
		$fasta{$actual}.= $li;
	}
	
}
close F2;

foreach my $nc (sort keys %clusterSeq){
	
	foreach my $c (@{$clusterSeq{$nc}}){
		$seq = $hashSeq{$c};
		push(@{$clusters{$nc}},$seq);
		push(@{$clustersMembers{$nc}},$c);
		#print"$nc\t$c\n$seq\n";
	}
	
}	

print "\n";

################### END Clustering and post-proccess



################## Distances calculation

print "Distances calculation : ";
$actualNumCluster=0;
$| = 1;  # Turn off buffering on STDOUT.

foreach $k (keys(%clusters))
{
	$output = join('',@{$clusters{$k}});
	##############
	`touch Cluster_$k\_$numSeqClust{$k}.fasta`;
	foreach $member (@{$clustersMembers{$k}}) 
	{
		my $seqecri = $hashSeq{$member};
		`echo '>$member\n$seqecri' >> Cluster_$k\_$numSeqClust{$k}.fasta`;
		if(($member =~ m/^([0-9]+)\|[0-9]+/) && ($numSeqClust{$k} >= 3)){
			$hashDiff{$k}{$1}=1;
			
		}
	}
	if($numSeqClust{$k} >= 3){
		my $val = 0;
		my @idval;
		foreach my $nu(keys %{$hashDiff{$k}}){
			++$val;
			push (@idval,$nu);
		}
		my $idvals = join(" - ", @idval); 
		if($input >= 5){
			#if($val >= $input-5 ){
				`echo "Cluster $k : $val dossiers => $idvals" >> bilan_cluster.txt`;
			#}
		}
		if($input < 5){
			#if($val == $input){
				`echo "Cluster $k : $val dossiers => $idvals" >> bilan_cluster.txt`;
			#}
		}
	}
	##############
	undef @uniqs;
	undef %repeat;
	foreach $member (@{$clustersMembers{$k}}) 
	{
		if(exists($repeat{$member}))
		{
			$repeat{$member}++;
		}else{
			push(@uniqs,$member);
			$repeat{$member}++;
		}
	}
	
	if (scalar(@uniqs) < scalar (@inFiles) * $threshold)
	{
		push(@inNet,$k);
		if(scalar(@{$clusters{$k}}) <2)
		{
			#print scalar(@{$clusters{$k}})."\n";
			open(TMP,"Cluster_$k\_$numSeqClust{$k}.fasta");
			$firstLine = <TMP>;
			close TMP;
			$firstLine =~ s/>//;
			@c = split('\|',$firstLine);
			my $val = scalar(@inFiles)/2;
			push(@net,"Cluster_$k\_$numSeqClust{$k}\t$taxasHash{$c[0]}\t$val\tUndirected");
		}else{
			if($fast eq 'Y' | $fast eq 'Yes' | $fast eq 'yes' | $fast eq 'y')
			{
				open(TMP,"Cluster_$k\_$numSeqClust{$k}.fasta");
				@tmpArray = <TMP>;
				@clusterFasta = grep(/>/,@tmpArray);
				foreach $l (@clusterFasta)
				{
					$l =~ s/>//;
					@c = split('\|',$l);
					my $val = scalar(@inFiles)/2;
					push(@net,"Cluster_$k\_$numSeqClust{$k}\t$taxasHash{$c[0]}\t$val\tUndirected");
				}
				
				close TMP;
			}else{
				$verif = 1;
				push(@net,distance("Cluster_$k\_$numSeqClust{$k}"));
			}
		}
	}
	$actualNumCluster++;
	
	print("\rProccessing cluster $actualNumCluster of $comptetot\n");
}
print "done\n";
print "\n";

################# Output section



###Output Network
$finalOutput = join("\n",@net);
open(FINAL,">$outFile");
print FINAL "Source\tTarget\tWeight\tType\n";
print FINAL "$finalOutput";
close FINAL;


###Output Annotation Table
print "Output Annotation Table : ";
%Synt = synteny();

open(TABLE, ">$outTable");
print TABLE "ID\tType\tTwinGroup\tDescription\n";
foreach $k (keys(%taxasHash))
{
	print TABLE "$taxasHash{$k}\tGU\tGU\t$k\n";
}

foreach $k (@inNet)
{

	foreach $l (@{$clustersMembers{$k}}){
		$tmp = $l;
		$tmp =~ s/>//;
		@c = split(/\|/,$tmp);
		$id = "Cluster_$k\_$numSeqClust{$k}";
		if (exists($Synt{$id}))
		{
			print TABLE "$id\tCluster\t$Synt{$id}\t$referenceHash{$c[0]}{$c[1]}";
		}else{
			print TABLE "$id\tCluster\tSingle\t$referenceHash{$c[0]}{$c[1]}";
		}
	}

}
print "done\n";
close TABLE;


if($clean eq "yes")
{
	foreach my $nc (sort keys %clusterSeq){
		if($numSeqClust{$nc} < 3){
			system("rm Cluster_$nc\_$numSeqClust{$nc}.fasta");
		}
	}

	system("mkdir $dir/ $dir/cd-hit/");
	system("mv rendu.tmp* $dir/cd-hit/");
	if($verif == 1){
		system("rm seq_temp temp");
		system("mkdir $dir/phylip/ $dir/muscle/ $dir/trimal/ $dir/phyml/");
		system("mv *.aln $dir/muscle/");
		system("mv *.phylip $dir/trimal/");
		system("mv *.rendu $dir/phylip/");
		system("for f in `ls *_phyml_stats.txt`; do mv \$f $dir/phyml/ ; done");
		system("for f in `ls *_phyml_tree.txt`; do mv \$f $dir/phyml/ ; done");
	}
	system("for f in `ls Cluster*`; do mv \$f $dir/cd-hit/ ; done");
	system("mv all_fasta.tmp $dir/");
	system("mv equivalence.txt $dir/");
	system("mv bilan_cluster.txt $dir/");
	system("mv *.csv $dir/");
}



######################## Subrutines #############################
sub synteny
{
	
	
	foreach $l (@net)
	{
		chomp $l;
		@c = split('\t',$l);
		push(@{$clusterTwin{$c[0]}},$c[1]);
	
	}
	foreach $k (keys(%clusterTwin))
	{
		if (@{$clusterTwin{$k}}>1)
		{
			@tmpArray = sort(@{$clusterTwin{$k}});
			$tmpHeader = join('=',@tmpArray);
			push(@{$maxCls{$tmpHeader}},$k);
		}
	}

	$i=0;
	foreach $k (keys(%maxCls))
	{
		if (scalar(@{$maxCls{$k}}>1))
		{
			$i++;
			foreach $tmpCluster (@{$maxCls{$k}})
			{
				
				$hashTwinGroup{$tmpCluster} = "Twin_$i";
			
			}
		}
	}
	return %hashTwinGroup;
}

sub distance 
{
	if( -e "outfile")
	{
		system("rm outfile");
	}
	if( -e "seq_temp")
	{
		system("rm seq_temp");
	}
	if ( -e "temp")
	{
		system("rm temp");
	}
	
	$file = $_[0];
	
	undef @salida;
	undef $map;
	undef %suma;
	undef @dists;
	undef $alnIn;
	undef $alnOut;
	undef $aln;
	undef @txt;
	
	
	system("muscle -quiet -in $file.fasta -out $file.aln");
	
	system("trimal -strictplus -phylip -in $file.aln -out $file.aln.phylip >/dev/null 2>&1");
###### Protdist Process 


	open(TMP,">seq_temp");
	print TMP "$file.aln.phylip\nY\n";
	close TMP;
	if($seqTY eq "nt"){
		system("phylip dnadist < seq_temp > temp");
	}
	if($seqTY eq "aa"){
		system("phylip protdist < seq_temp > temp");
	}
	system("mv outfile $file.aln.phylip.rendu");
	if(open(D,"$file.aln.phylip.rendu"))
	{@txt = <D>;}
	close D;

	
	shift(@txt); ##remove header	

	for($i=0; $i<scalar(@txt); $i++)
	{
		if($txt[$i] =~ /^\d/)
		{
			@c = split(' ',$txt[$i]);
			@c2 = split('\|',$c[0]);			
			$actual = $c2[0];
			$suma{$actual} =0;
	
			for($j=1;$j<scalar(@c);$j++)
			{
				$suma{$actual} += abs($c[$j]);
			}
			
		}else{
			@c = split(' ',$txt[$i]);
			for($j=0;$j<scalar(@c);$j++)
			{
				$suma{$actual} += abs($c[$j]);
			}
		}
	}
	@val = values(%suma);
	@dists =  sort {$a <=> $b} @val;
	
	foreach $strain (keys(%suma))
	{
		$n = scalar(keys(%suma));
		if($suma{$strain} eq 0)
		{
			$dist = $n/2;
		}else{
			$dist = log(($n*$n)/$suma{$strain});
			if($dist <= 0)
			{
				$dist =0.01;
			}
		}
		
		push(@salida,"$file\t$taxasHash{$strain}\t$dist\tUndirected");  

	}
	if($file =~ m/Cluster_.+_([0-9]+)/){
		if($1 > 2){
			system("phyml -d $seqTY -i $file.aln.phylip");
		}
	}
	
	return @salida;
}	






