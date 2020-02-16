#!/usr/bin/perl

my $myFile = $ARGV[0];
my $bins = $ARGV[1];
my $firstBin=-1; 
my $lastBin=-1;
my $column =0;
my $absolute=0;
my $firstRow = 0;
my $log = 0;
if($#ARGV>=2){
	$firstBin=$ARGV[2];
	if($#ARGV>=3){
		$lastBin=$ARGV[3];
		if($#ARGV>=4){
		    $column=$ARGV[4];
		    if($#ARGV>=5){
			$absolute=$ARGV[5];
			if($#ARGV>=6){
			    $firstRow=$ARGV[6];
			    if($#ARGV>=7){
				$log=$ARGV[7];
			    }
			}
		    }
		}
	}
}

unless(open(FILE, $myFile)){
	die "Cannot open file\n";}
	
@lines = <FILE>;
@scores=();
my $mean=0;
my $max=-1000000;
my $min = 1000000;

$count =0;
#Find mean, max, min
for($i=$firstRow; $i<=$#lines; $i++){
	@currLine = split(/\s+/, $lines[$i]);
	#if($#currLine==2 && $lines[$i] !~ m/loaded/){
		$scores[$count]=$currLine[$column];
	        
                if($absolute!=0){$scores[$count]=abs($scores[$count]);}

                if($log!=0 && $scores[$count]==0){$scores[$count]=1;}
                if($log!=0){$scores[$count]=log($scores[$count])/log(10);}
		#if($scores[$count]>0 && $scores[$count]<1){
		
		$mean+=$scores[$count];
		if($scores[$count]>$max){
			$max = $scores[$count];
		}elsif($scores[$count]<$min){
			$min = $scores[$count];
		}#print "$count\t$scores[$count]\n";
		
		$count++;
	#}
}$mean = $mean/$count;

#initialize the histogram
if($firstBin==-1){
	$firstBin=$min; $lastBin=$max;
}
$step = ($lastBin-$firstBin)/$bins;
my @histo=();
for($x=0; $x<=$bins; $x++){
	$histo[$x]=0;
}
#print "$firstBin\t$lastBin\n";
$stdDev=0; $tmp=0;
#Calculate standard deviation and populate histogram
for($i=0; $i<$count; $i++){
	$tmp+=($scores[$i]-$mean)*($scores[$i]-$mean);
	if(($scores[$i]-$firstBin)/$step<0){
		$histo[0]++;
	}elsif(($scores[$i]-$firstBin)/$step>$bins){
		$histo[$bins]++;
	}else{
		$histo[($scores[$i]-$firstBin)/$step]++;
	}
}
$stdDev = sqrt($tmp/$count);

@sortScores = sort {$a <=> $b} @scores;
$median = $sortScores[$#sortScores/2];
print "$myFile\n";
print "Mean\t$mean\nStdDev\t$stdDev\nMedian\t$median\n\n";

for($x=0; $x<$bins; $x++){
	if($log!=0){printf("%lf\t%.10lf\n", 10**($firstBin+($x*$step)), $histo[$x]/$count);
	}else{printf("%lf\t%.10lf\n", ($firstBin+($x*$step)), $histo[$x]/$count);
	}
}
if($log!=0){printf("More\t%.10lf\n", $histo[$bins]/$count);
}else{printf("More\t%.10lf\n", $histo[$bins]/$count);
}
print "\n\n";
close(FILE);
