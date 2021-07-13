## SUMMARIZES 'ChromoPainter' OUTPUT BY SUMMING CHUNK COUNTS ACROSS ALL CHROMOSOMES

## usage:  perl ChromoPainterOutputSum.pl filenamePRE filenamePOST

## example:   perl ChromoPainterOutputSum.pl AllFrenchYoruba30gen30propchr _AllvALL.chunklengths.out

###############################
## INPUT:

@chromovec=("21","22");

###########################
## PROGRAM:

$filenamePRE = "$ARGV[0]";
$filenamePOST = "$ARGV[1]";

$numchrom = @chromovec;

     ## (I) GET CHUNK-COUNTS AND SUM ACROSS CHROMOSOMES:
$outfile = "${pathway}${filenamePRE}All${filenamePOST}";
print "$outfile\n";
open(OUT,">$outfile");
@totalcountsmat=();
$count = 0;
@recipientlabels=();
for ($i=0; $i < $numchrom; $i++)
{
    $infile = "${filenamePRE}$chromovec[$i]${filenamePOST}";
    open(IN,"$infile");
    $line=<IN>;
    if ($i==0)
    {
	print OUT "$line";
    }
    if ($i==0)
    {
	@labels=split(/\s+/,$line);
	shift(@labels);
	$numdonors = @labels;
	for ($j=0; $j<$numdonors; $j+=1)
	{
	    for ($k=0; $k<$numdonors; $k+=1)
	    {
		$totalcountsmat[$j][$k] = 0.0;
	    }
	}
    }
    $indcount=0;
    while(<IN>)
    {
	$line=$_;
	@linearray=split(/\s+/,$line);
	if ($i==0)
	{
		$recipientlabels[$indcount]=$linearray[0];
        }
	shift(@linearray);
	for ($k=0; $k<$numdonors; $k+=1)
	{
	    $totalcountsmat[$indcount][$k] = $totalcountsmat[$indcount][$k] + $linearray[$k];
	}
	$indcount=$indcount+1;
    }
    if ($i>0 && $previndcount != $indcount)
    {
	print "ERROR -- mismatched number of inds !!!!!!!!!!!!!!!!!!!!!!!!\n";
	print "$indcount $previndcount $chromovec[$i]\n";
	die; 
    }
    $previndcount=$indcount;
}
$numrecipients=@recipientlabels;

     ## (II) PRINT OUT:
for ($j=0; $j < $numrecipients; $j+=1)
{
    print OUT "$recipientlabels[$j]";
    for ($k=0; $k<$numdonors; $k+=1)
    {
	$printval = sprintf("%.3f",$totalcountsmat[$j][$k]);
	print OUT " $printval";
    }
    print OUT "\n";
}
