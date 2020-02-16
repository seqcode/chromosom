#!/usr/bin/perl

my $myFile = $ARGV[0];
my $col = $ARGV[1];
my $pattern=$ARGV[2];
my $wordmatch = $ARGV[3]; #Should be -w

$w=0;
if($wordmatch eq "-w"){
    $w=1;
}
unless(open(FILE, $myFile)){
	die "Cannot open file\n";}

$count=0;
while($curr=<FILE>){    
    @currLine = split(/\s+/, $curr);
    $val =$currLine[$col];

    if($w==1){
	if($val eq $pattern){print "$curr";}
    }else{
	if($val =~ m/$pattern/){print "$curr";}
    }
#    if($val eq $pattern){$count++;}
#system("grep $val $targetFile");
}#print "$count\n";

close(FILE);
