#!/usr/bin/perl

# perl process.pl [input file] [output dir]
my $myFile = $ARGV[0];
my $outDir = $ARGV[1];

$cmd = "mkdir $outDir";
system($cmd);

unless(open(FILE, $myFile)){
	die "Cannot open file\n";}
@lines = <FILE>;

for($x=0; $x<=$#lines; $x++){
    chomp($lines[$x]);
    @curr = split(/\t/, $lines[$x]);
    if($curr[0] ne "File accession"){
	$acc = $curr[0];
	$assay = $curr[4];
	$target = $curr[18];
	$reps = $curr[30];
	$reps =~ s/,\ /-/g;
	
	$file = $acc.".bed";
	$newname = $assay."_".$target."_".$reps."_".$acc.".bed";
	
	$cmd = "mv $file $outDir/$newname";
	print "$cmd\n";
	system($cmd);
    }
}
close(FILE);
