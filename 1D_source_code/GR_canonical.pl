#!/usr/bin/perl -w

# Program that analyzes data from Wang-Landa simulations
# Performs a jackknife analysis on the density of states

if ($#ARGV != 0) {	
	print "\nIncorrect number of command line arguments!\n"; 
	exit;
}
 
$N = $ARGV[0];		#Number of monomers

$dT = 0.1;		#Temperature increment
$Tmin = 0.01;	#Minimum temperature
$Tmax = 3.0;	#Maximum temperature

@tmp_DOS = ();
@tmp_GR = ();
@data_DOS = ();
@data_GR = ();
@results_DOS = ();
@results_GR = ();

@E = ();
@Rdist = ();
@PE = ();

sub ltrim($);

#Open the DOS files - reads E and g(E)
open(INPUT_DOS, "<DOS_H_iter020.dat") || die("Error: could not open file! (INPUT DOS)\n");  # Open the file
@tmp_DOS = <INPUT_DOS>;
@tmp_DOS = grep(!/^#/, @tmp_DOS);
chomp(@tmp_DOS);
close(INPUT_DOS);
	
for( $i=0; $i<@tmp_DOS ; $i++)     
{
	@results_DOS = split(/\s+/,$tmp_DOS[$i]);
	$E[$i] = $results_DOS[0]*$N;     #Energy of each bin in the DOS
	#$DOS[$i] = $results_DOS[1]; 
		
	$data_DOS[0][$i] = $results_DOS[1];
	#print "$E[$i]\t$data_DOS[$r][$i]\n";

	$PE[$i]=0.0;
}			


#Open the DOS files - reads E and g(E)
open(INPUT_GR, "<gR_E.dat") || die("Error: could not open file! (gR_E)\n");  # Open the file
@tmp_GR = <INPUT_GR>;
@tmp_GR = grep(!/^#/, @tmp_GR);
chomp(@tmp_GR);
close(INPUT_GR);
	
#print @tmp_GR;
	
$j=0;	#Energy index
$k=0;	#R index
for( $i=0; $i<@tmp_GR ; $i++)     
{
	if ($tmp_GR[$i] =~ /^(\s)*$/) {$k=0; $j++; next; }   # skip blank lines
	
	@results_GR = split(/\s+/,$tmp_GR[$i]);
	$Rdist[$k] = $results_GR[1];
	$data_GR[$j][$k] = ($results_GR[2]);	
	#print "$j\t$k\t$data_GR[$j][$k]\n";
	$k++;
}			

#Test to make sure data loaded correctly - use diff on command line and compare
#for( $j=0; $j<@tmp_DOS ; $j++)     
#{
#	for( $k=0; $k<@Rdist ; $k++)     
#	{
# 		print "$E[$j]\t$Rdist[$k]\t$data_GR[$j][$k]\n";
#	}
#	print "\n";
#}  


for($k=0;$k<@Rdist;$k++)
{
	for($T=$Tmin;$T<=$Tmax;$T=$T+$dT)
	{
		
		$lambda = -1.0e300;
		$temp = 0.0;
		$Z = 0.0;
		$Gsum = 0.0;
		
		for($j=0;$j<@E;$j++)
		{
			$temp = $data_DOS[0][$j] - $E[$j]/$T;
			if($lambda < $temp)
			{
	    		$lambda = $temp;	
			}
		}
		
		for($j=0;$j<@E;$j++)
		{
			$PE[$j] = exp( $data_DOS[0][$j] - $E[$j]/$T -$lambda );
			$Z = $Z + $PE[$j];
			$Gsum = $Gsum + $data_GR[$j][$k]*$PE[$j];
			#print "$Rdist[$k]\t$j\t$PE[$j]\n";
		}
		
		if($Z>0.0)
		{
			$Gsum = $Gsum/$Z;
			#print "$lambda\t$Z\t$Gsum\n";
			print "$Rdist[$k]\t$T\t$Gsum\n";
		}

	}
	print "\n";
}



# Left trim function to remove leading whitespace
sub ltrim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}








