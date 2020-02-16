#!/usr/bin/perl

#Format & log convert interchromosomal interaction matrix
#Usage:
# perl processMatrix.pl hg19.info resolution matrixfile.tsv coordsfile.coord outfile minval maxval

my $fileInfo = $ARGV[0];
my $resolution = $ARGV[1];
my $matrixFile = $ARGV[2];
my $coordsFile = $ARGV[3];
my $outfile = $ARGV[4];
my %chrLen=();

my $minVal = log($ARGV[5]);
my $maxVal = log($ARGV[6]);

#Load chromosome lengths
print "Loading chromosome info\n";
unless(open(INFO, $fileInfo)){
	die "Cannot open chromosome info file\n";}
@linesInfo = <INFO>;
my %chrLen=();
for($i=0; $i<=$#linesInfo; $i++){
    @curr = split(/\s+/, $linesInfo[$i]);
    $chr = $curr[0]; 
    $chrLen{$chr}=$curr[1];
}
close(INFO);

#Load coordinate list
print "Loading coordinates\n";
unless(open(COORDSF, $coordsFile)){
	die "Cannot open peaks file\n";}
@linesCoords = <COORDSF>;
my @coords=();
my $numBins=0;
for($i=0; $i<=$#linesCoords; $i++){
    @curr = split(/\s+/, $linesCoords[$i]);
    $chr = "chr".$curr[0]; 
    $start=$curr[1]+1;
    $end=$curr[1]+$resolution;
    if(exists($chrLen{$chr}) && $end > $chrLen{$chr}){
	$end = $chrLen{$chr};
    }
    $coords[$i] = $chr.":".$start."-".$end;
    #print "$coords[$i]\n";
    $numBins++;
}
close(COORDSF);



#Initialize matrix to zeros
print "Initializing $numBins x $numBins matrix\n";
my @matrix = ();
for($x=0; $x<$numBins; $x++){
    for($y=0; $y<$numBins; $y++){
 	$matrix[$x][$y]=0;
    }
}

#Load matrix from file and log transform
unless(open(MAT, $matrixFile)){
    die "Cannot open file $matrixFile\n";}
$i=0;
while (defined(my $line = <MAT>)){
#@linesMat = <MAT>;
#for($i=0; $i<=$#linesMat; $i++){
#    @curr = split(/\s+/, $linesMat[$i]);
    @curr = split(/\s+/, $line);
    for($j=0; $j<=$#curr; $j++){
	$val = $curr[$j];
	if($i==$j){
	    $lnval = 0;
	}else{
	    if($val==0){
		$lnval = $minVal;
	    }else{
		$lnval = log($val);
		if($lnval > $maxVal){$lnval = $maxVal;}
		if($lnval < $minVal){$lnval = $minVal;}
	    }
	}
	$matrix[$i][$j] = $lnval;
    }
    $i++;
}
close(MAT);

#Print the matrix
print "Printing output matrix\n";
unless(open(OUT, ">$outfile")){
    die "Cannot open output file\n";}
for($i=0; $i<$numBins; $i++){
    print OUT "\t$coords[$i]";
}
print OUT "\n";
for($x=0; $x<$numBins; $x++){
    print OUT "$coords[$x]";
    for($y=0; $y<$numBins; $y++){
	print OUT "\t$matrix[$x][$y]";
    }
    print OUT "\n";
}
close(OUT);
