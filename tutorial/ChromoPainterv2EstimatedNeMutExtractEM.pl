## EXTRACTS E-M ESTIMATED N_E FOR EACH INDIVIDUAL, WEIGHTED-AVERAGES ACROSS CHROMOS, AND PUTS INTO SIMPLIFIED FILE (NOT AVERAGED ACROSS INDS):

## usage:   perl ChromoPainterv2EstimatedNeExtractEM.pl

######################################
## INPUT:

$infilePREFIX="output_estimateEM_Chr";
$infileSUFFIX=".EMprobs.out";

@chromovec=("21","22");
@chromolengths=(7087,6811);    ## number of SNPs

#######################################
## PROGRAM:

$numchrom=@chromovec;
$sumchromolengths=0;
for ($i=0; $i < $numchrom; $i+=1)
{
    $sumchromolengths=$sumchromolengths+$chromolengths[$i];
}

@loglikvec=();
@NeEstvec=();
@MutEstvec=();
for ($i=0; $i < $numchrom; $i+=1)
{
	$infile="${infilePREFIX}$chromovec[$i]${infileSUFFIX}";
	print "$infile\n";
	open(IN,"$infile");
	$count=0;
	$firstline=0;
	while(<IN>)
	{
	    $line=$_;
	    @linearray=split(/\s+/,$line);
	    $linelength=@linearray;
	    if (($linelength==1) && ($firstline==1))
	    {
		$loglikvec[$count]=$loglikvec[$count]+$loglikvalue;
		$NeEstvec[$count]=$NeEstvec[$count]+($chromolengths[$i]/$sumchromolengths)*$NeEst;
		$MutEstvec[$count]=$MutEstvec[$count]+($chromolengths[$i]/$sumchromolengths)*$MutEst;
		if ($count==0)
		{
		    #print "$k $numchrom $chromovec[$i] $loglikvalue $NeEst $NeEstvec[$count]\n";
		}
		$count=$count+1;
	    }
	    if ($linelength ne "IND")
	    {
		$linelength=@linearray;
		$loglikvalue=$linearray[1]+$linearray[2];
		$NeEst=$linearray[($linelength-2)];
		$MutEst=$linearray[($linelength-1)];
	    }
	    $firstline=1;
	}
	$loglikvec[$count]=$loglikvec[$count]+$loglikvalue;
	$NeEstvec[$count]=$NeEstvec[$count]+($chromolengths[$i]/$sumchromolengths)*$NeEst;
	$MutEstvec[$count]=$MutEstvec[$count]+($chromolengths[$i]/$sumchromolengths)*$MutEst;
        if ($i>0 && $previndcount != $count)
        {
             print "ERROR -- mismatched number of inds !!!!!!!!!!!!!!!!!!!!!!!!\n";
             print "$count $previndcount $popvec[$k] $chromovec[$i]\n";
             die;
        }
        $previndcount=$count;
	#print "$count $previndcount\n";
}
$NeEstMean=0.0;
$MutEstMean=0.0;
for ($i=0; $i < $count; $i+=1)
{
        $indval=$i+1;
       	$NeEstMean=$NeEstMean+$NeEstvec[$i];
       	$MutEstMean=$MutEstMean+$MutEstvec[$i];
}
$NeEstMean=sprintf("%.3f",$NeEstMean/$count);
$MutEstMean=sprintf("%.6f",$MutEstMean/$count);
print "$NeEstMean $MutEstMean\n";

