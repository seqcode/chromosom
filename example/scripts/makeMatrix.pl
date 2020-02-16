#!/usr/bin/perl

#Convert chromatin interaction bed files to interaction matrix
#Usage:
# perl makeMatrix.pl hg19.info [resolution] outname *bed

my $fileInfo = $ARGV[0];
my $resolution = $ARGV[1];
my $outfile = $ARGV[2];
my %chrLen=();
unless(open(INFO, $fileInfo)){
	die "Cannot open peaks file\n";}
@linesInfo = <INFO>;
for($i=0; $i<=$#linesInfo; $i++){
    @curr = split(/\s+/, $linesInfo[$i]);
    $chr = $curr[0]; 
    $chrLen{$chr}=$curr[1];
}

#Human genome
my @chrs = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22");

#Index the coordinates
my %coordIndex=();
my @coordStrs=();
my $numBins=0;
for ($c=0; $c<=$#chrs; $c++){
    $chr = $chrs[$c];
    $clen = $chrLen{$chr};
    
    for($x=1; $x<$clen; $x+=$resolution){
	$xe = $x+$resolution;
	if($xe>$clen){
	    $xe=$clen;
	}
	$cstr1 = $chr.":".$x;
	$cstr2 = $chr.":".$x."-".$xe;
	$coordIndex{$cstr1}=$numBins;
	$coordStrs[$numBins]=$cstr2;
	$numBins++;
    }
}

#Initialize matrix to zeros
print "Initializing $numBins x $numBins matrix\n";
my @matrix = ();
for($x=0; $x<$numBins; $x++){
    for($y=0; $y<$numBins; $y++){
	$matrix[$x][$y]=0;
    }
}

#Cycle through files, adding interactions to matrix
for($f=3; $f<=$#ARGV; $f++){
    $file = $ARGV[$f];
    print "$file\n";
    unless(open(F1, $file)){
	die "Cannot open file $file\n";}
    @lines = <F1>;
    for($i=0; $i<=$#lines; $i++){
	@curr = split(/\s+/, $lines[$i]);
	$chrA=$curr[0];
	$startA = $curr[1]+1;
	$endA = $curr[2];
	$A = $chrA.":".$startA;
	$chrB=$curr[3];
	$startB = $curr[4]+1;
	$endB = $curr[5];
	$B = $chrB.":".$startB;
	$val = $curr[6];

	$Aindex = $coordIndex{$A};
	$Bindex = $coordIndex{$B};

	$matrix[$Aindex][$Bindex]=$val;
	$matrix[$Bindex][$Aindex]=$val;
    }
    close(F1);
}


#Print the matrix
print "Printing output matrix\n";
unless(open(OUT, ">$outfile")){
    die "Cannot open output file\n";}
for($i=0; $i<$numBins; $i++){
    print OUT "\t$coordStrs[$i]";
}
print OUT "\n";
for($x=0; $x<$numBins; $x++){
    print OUT "$coordStrs[$x]";
    for($y=0; $y<$numBins; $y++){
	print OUT "\t$matrix[$x][$y]";
    }
    print OUT "\n";
}
close(OUT);
