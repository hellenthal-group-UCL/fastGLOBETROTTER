#!/usr/bin/perl
## CONVERTS PHASED IMPUTE2 OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES

### AUTHOR: Daniel Lawson (dan.lawson@bristol.ac.uk)
### See "usageHelp" for details (run with no options)
### Copyright 2013 Daniel Lawson
### LICENCE: GPLv3
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


sub help {
print("CONVERTS PHASED IMPUTE2 OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES\n");
print("   copyright Daniel Lawson (dan.lawson@bristol.ac.uk) 2013, released as GPLv3\n");
print("usage:   perl impute2chromopainter.pl <options> impute_output_file.haps recommap_filein output_filename_prefix\n");

print("where:\n");
print("        (i) impute_output_file.haps = directory location and filename of IMPUTE2 output file with suffix \"haps\" that contains **PHASE OUTPUT** haplotypes\n");
print("        (ii) recommap_filein = directory location and filename of genetic map used to phase haplotypes in IMPUTE2 (though see -r option)\n");
print("        (iii) output_filename_prefix = directory location and filename prefix for chromopainter input files (two output files with suffixes \".haps\" and \".recomrates\" will be created by this program)\n\n");

print("<options>:\n");
print("-H:              The organism is HAPLOID, not DIPLOID (as otherwise assumed)\n");
print("-p:		By default, this script produces Finestructure-ChromoPainter style input, which differs \n");
print("			from PHASE output which omits the first line. Specify this option to create standard PHASE format\n\n");
print("-m <val>:        Change the number of Impute2 inds to read at a time.  (To minimize RAM requirements - makes program run a factor \n");
print("                 (N/<val>) slower than if you set <val>=N, where N is the total number of individuals in\n");
print("                 your \".haps\" IMPUTE2 output file, but also uses a factor of (N/<val>) less RAM\n");
print("-r <val>         OMIT the recombination file. You should then provide only 2 file names. You must provide the number of individuals in the file as an argument. Omit the recombination map and provide only \"impute_output_file.haps\" and \"output_filename_prefix\".  You should then use some other method of generating a recombination map, for example, makeuniformrecfile.pl.\n");
print("-q               Quiet mode.\n");
print("-h               This help file.\n\n");
print("!!! NOTE:  All 3 files are required; This program does NOT properly check that the imput files are correct and will not \n");
print("     produce helpful error messages if they are not!\n");
print("SPECIFICALLY: You MUST run this on the output of IMPUTE2's \"-phase\" _haps file, NOT the genotype file!\n");
die "\n";
}

use Switch;


###############################
## INPUT:



use Switch;
use Getopt::Long;
use strict;

###############################
## ARGUMENT PROCESSING

my $quiet=0;
my $verbose=1;
my $help=0;
my $phasemode=0;
my $fsmode=1; ## finestructure mode (i.e. start with an additional line containing 0)
my $haploid=0; # whether the organism is haploid/diploid
my $omitrecom=0; # whether to omit the recombination map
my $nfromrecom=0; # what we think "N" is from the recombination rate omission line
my $totalINDS=-1; # total number of individuals in the file
my $totalHAPS=-1; # total number of haplotypes in the file
my $Mb = 1000000.0;
my $numindsmaxsize=500;     ## ONLY READ IN THIS NUMBER OF IMPUTE2 INDS AT A TIME (TO MINIMIZE RAM REQUIREMENTS -- MAKES PROGRAM RUN A FACTOR OF (N/$numindsmaxsize) SLOWER THAN IF YOU SET $numindsmaxsize=N, where N is the total number of individuals in your ".haps" IMPUTE2 output file, but also uses a factor of (N/$numindsmaxsize) less RAM

my $IMPUTEinfile="";
my $recommapinfile="";
my $outfilePRE="";

GetOptions ('h|help' => \$help, 
	    'H|Haploid' => \$haploid,
	    'q|quiet' => \$quiet,
	    'r|recom=s' => \$nfromrecom,
	    'p|phase' => \$phasemode,
	    'm|map=s' => \$numindsmaxsize);

$verbose=1-$quiet;
$fsmode=1-$phasemode;
if($nfromrecom>0) {
    $omitrecom=1;
}
if($help) {help();}
if(scalar(@ARGV)==0) {help();}
if((scalar(@ARGV)!=3 && !$omitrecom)|| (scalar(@ARGV)!=2 && $omitrecom) ) {
    print "Incorrect number of arguments!\n";
    help();
}
my $fileon=0;
$IMPUTEinfile="$ARGV[$fileon]"; ++$fileon;
if(!$omitrecom) { $recommapinfile="$ARGV[$fileon]";++$fileon;}
$outfilePRE="$ARGV[$fileon]"; ++$fileon;

if($verbose){
    if($haploid) {print("Treating organisms as HAPLOID\n");
    }else{print("Treating organisms as DIPLOID\n");}
}

##############################
## PROGRAM:

              ## (I) GET RECOM-RATE INFO:
my @posvecFULL=(); # positions in the recombination rate file
my @posvec; # location of ALL snps
my @recomrateFULL=(); # recombination rates

my $ploidy=2-$haploid; # 2 if diploid, 1 if haploid
my $line;
my @linearray; # temporary array of the current line
my $nsitesFULL; # number of recombination rate change sites

if($omitrecom) {
    print ("OMITTING recombination rate file\n");
}else {
    if($verbose){ print("Reading recombination map file $recommapinfile...\n");}
    open IN,"$recommapinfile" or die "Could not open recombination map file $recommapinfile !\n";
    $line=<IN>;        ## header
    
    my $numfoundind=0; # 
    my $prevpos=0; # 
    
    while(<IN>)
    {
	my $line=$_;
	@linearray=split(/\s+/,$line);
	push(@posvecFULL,$linearray[0]);
	push(@recomrateFULL,$linearray[1]);
	if ($linearray[0] <= $prevpos)
	{
	    print "genetic map error -- positions out of order in genetic map input file ($linearray[0] <= ${prevpos})\n";
	    die;
	}
	$prevpos = $linearray[0];
    }
    $nsitesFULL=@posvecFULL; # number of recombination rate change sites
}

              ## (II) GET NUMBER OF SITES AND INDS: 
if($verbose){ print("Reading Impute format file $IMPUTEinfile...\n");}
open IN,"$IMPUTEinfile" or die "Could not open input file $IMPUTEinfile\n";
$line=<IN>;
@linearray=split(/\s+/,$line);
if($omitrecom) {
    $totalINDS=$nfromrecom;
    if($verbose){print "Using $totalINDS individuals as specified via \"-r\"\n";}

}else{
    $totalINDS=(@linearray-5)/$ploidy; # 
    if($verbose){print "Reading $totalINDS individuals\n";}

}
my $numsplits=int($totalINDS/$numindsmaxsize);
if (($numsplits*$numindsmaxsize)<$totalINDS)
{
    $numsplits=$numsplits+1;
}

if($verbose){print "Using $numsplits split(s) of the individuals (change with -m, see help)\n";}
if($verbose){print "Counting sites...\n";}

my $nsites=1;
while(<IN>)
{
    $line=$_;
    $nsites=$nsites+1;
}
if($verbose){print "Found $nsites sites...\n";}

my $warnings=0;
              ## (III) READ IN IMPUTE2 HAPLOTYPES AND MAKE CHROMOPAINTER HAPLOTYPE INPUT FILE:
    if($verbose) {print "Starting processing of IMPUTE2 file, to output file ${outfilePRE}.haps...\n";}
open OUT,">${outfilePRE}.haps" or die("Unable to open output file ${outfilePRE}.haps\n");
open OUT3,">${outfilePRE}.warnings" or die("Unable to open output file ${outfilePRE}.errors\n");
#if($fsmode==1) {
#	print OUT "0\n";
#}
my $snpcount;
for (my $a=0; $a < $numsplits; $a+=1)
{
    my $aa=$a+1;
    my $startIND=$a*$numindsmaxsize;
    my $endIND=($a+1)*$numindsmaxsize;
    if ($endIND > $totalINDS)
    {
	$endIND=$totalINDS;
    }

                              ## read in:
    open IN, "$IMPUTEinfile" or die("Error: cannot open impute file $IMPUTEinfile on pass $a\n");
    my @rsvec=();
    my @genomat=();
    #my $snpcount = 0;
    my $tokeepSNP = 1;
    my $posval=0;
    my $rsval=0;
    my $prevpos=0;
    $snpcount = 0;
    while(<IN>)
    {
	$line=$_;
	@linearray=split(/\s+/,$line);
	$rsval=$linearray[1];
	$posval=$linearray[2];
	shift(@linearray);
	shift(@linearray);
	shift(@linearray);
	shift(@linearray);
	shift(@linearray);

	$tokeepSNP = 1;
	if ($snpcount > 0 && $posval <= $prevpos)
	{
	    if ($posval < ${prevpos})
	    {
		print "genetic map error -- positions out of order in haplotype input file (${posval} < ${prevpos}). Exiting....\n";
		die;
	    }
	    if ($posval == ${prevpos})
	    {
		print "WARNING: genetic map error -- position duplicated in haplotype input file (${posval} == ${prevpos}). Will use only first instance of position.\n";
		$tokeepSNP=0;
	    }
	}
	$prevpos = $posval;

	if ($tokeepSNP==1)
	{
	    push(@rsvec,$rsval);
	    push(@posvec,$posval);
	    for (my $i=$startIND; $i < $endIND; $i+=1)
	    {
		if(($linearray[(2*$i)]!=0 && $linearray[(2*$i)]!=1)||($linearray[(2*$i)+1]!=0 && $linearray[(2*$i+1)]!=1)){
		    die("Error: Invalid values in IMPUTE2 phase output file ($IMPUTEinfile) - is this really the _haps output of a -phase run?\n");
		}
		$genomat[(($i-$startIND)*2)][$snpcount]=$linearray[(2*$i)];
		$genomat[(($i-$startIND)*2+1)][$snpcount]=$linearray[(2*$i+1)];
	    }

	    $snpcount=$snpcount+1;
	}
    }

                                ## print out:	
    if($verbose) {print "Processing split $aa of $numsplits... writing output...\n";}
    if ($a==0)
    {
	$totalHAPS=$ploidy*$totalINDS;
	print OUT "$totalHAPS\n";
	#print OUT "$totalINDS\n";
	print OUT "$snpcount\n";
	print OUT "P @posvec\n";
	#for (my $j=0; $j < $nsites; $j+=1)
	#{
	#    print OUT "S";
	#}
	#print OUT "\n";
    }
    for (my $i=0; $i < ($ploidy*($endIND-$startIND)); $i+=1)
    {
	for (my $j=0; $j < $snpcount; $j+=1)
	{
	    print OUT "$genomat[$i][$j]";
	}
	print OUT "\n";
    }
}

                 ## (IV) MAKE CHROMOPAINTER RECOM-MAP INPUT FILE:
my $start = 0;
my @recomrate=();

if($omitrecom) {
    if($verbose) {print "Not creating recombination map...\n";}
}else{
    if($verbose) {print "Creating chromopainter format recombination map...\n";}
## This section is over-complicated
#############################
    my $initialzeros=0;
    if ($posvec[0] < $posvecFULL[0])
    {
	#die "ERROR - first basepair of file is less than first basepair of genetic map!\n";
	print "WARNING - first basepair of file is less than first basepair of genetic map! Will 'impute' using the rate from the first included SNP in the region -- but this may not behave well!\n";
	$initialzeros=1;
    }
    for (my $i=0; $i < ($snpcount-1); $i+=1)
    {
	push(@recomrate,0);
	if ($posvec[$i] >= $posvecFULL[($nsitesFULL-1)])
	{
	    $recomrate[$i] = $recomrateFULL[($nsitesFULL-2)]/$Mb;
	}
	if ($posvec[$i] < $posvecFULL[($nsitesFULL-1)])
	{
	    for (my $j=$start; $j < ($nsitesFULL-1); $j+=1)
	    {
		if (($posvec[$i] >= $posvecFULL[$j]) && ($posvec[$i] < $posvecFULL[($j+1)]) && ($posvec[($i+1)] <= $posvecFULL[($j+1)]))
		{
		    $recomrate[$i] = $recomrateFULL[$j]/$Mb;  
		    $start = $j;
		    last;
		}
		if (($posvec[$i] >= $posvecFULL[$j]) && ($posvec[$i] < $posvecFULL[($j+1)]) && ($posvec[($i+1)] > $posvecFULL[($j+1)]))
		{
		    my $recomcurrent = $recomrateFULL[$j]*($posvecFULL[($j+1)]-$posvec[$i]);
		    my $endspot = $j+1;
		    if ($endspot == ($nsitesFULL-1))
		    {
			$recomcurrent = $recomcurrent + $recomrateFULL[$j]*($posvec[($i+1)]-$posvecFULL[($j+1)]);
			$recomrate[$i] = ($recomcurrent/($posvec[($i+1)]-$posvec[$i]))/$Mb;  
			last;
		    }
		    while($posvec[($i+1)] > $posvecFULL[($endspot+1)])
		    {
			$recomcurrent = $recomcurrent + $recomrateFULL[$endspot]*($posvecFULL[($endspot+1)]-$posvecFULL[$endspot]);
			$endspot = $endspot+1;
			if ($endspot == ($nsitesFULL-1))
			{
			    last;
			}
		    }
		    $recomcurrent = $recomcurrent + $recomrateFULL[$endspot]*($posvec[($i+1)]-$posvecFULL[$endspot]);
		    $recomrate[$i] = ($recomcurrent/($posvec[($i+1)]-$posvec[$i]))/$Mb;  
		    $start = $j;
		    last;
		}
	    }
	}
	if ($recomrate[$i] == 0)
	{
	    print OUT3 "Warning: zero recombination rate at SNP $i\n";
	    $warnings=1;
	    if(($posvec[$i] < $posvecFULL[0]) && ($posvec[($i+1)] <= $posvecFULL[1]))
	    {
		$recomrate[$i] = $recomrateFULL[0]/$Mb;
	    }
	    if(($posvec[$i] < $posvecFULL[0]) && ($posvec[($i+1)] > $posvecFULL[1]))
	    {
		my $recomcurrent = $recomrateFULL[0]*($posvecFULL[1]-$posvec[$i]);
		my $endspot = 1;
		while($posvec[($i+1)] > $posvecFULL[($endspot+1)])
		{
		    $recomcurrent = $recomcurrent + $recomrateFULL[$endspot]*($posvecFULL[($endspot+1)]-$posvecFULL[$endspot]);
		    $endspot = $endspot+1;
		    if ($endspot == ($nsitesFULL-1))
		    {
			last;
		    }
		}
		$recomcurrent = $recomcurrent + $recomrateFULL[$endspot]*($posvec[($i+1)]-$posvecFULL[$endspot]);
		$recomrate[$i] = ($recomcurrent/($posvec[($i+1)]-$posvec[$i]))/$Mb;  
	    }
	}
    }
## print:
    open(OUT2,">${outfilePRE}.recomrates");
    print OUT2 "start.pos recom.rate.perbp\n";
    if ($initialzeros==1)
    {
	my $stopbp=0;
	for (my $i=0; $i < ($snpcount-1); $i+=1)   ## if data positions start prior to genetic map positions, find first non-zero recom value
	{
	    if ($recomrate[$i]>0)
	    {
		$stopbp=$i;
		last;
	    }
	}
	for (my $k=0; $k < $stopbp; $k+=1)   ## set all zero values at left to this first non-zero recom value
	{
	    $recomrate[$k]=$recomrate[$stopbp];
	}
    }
    for (my $i=0; $i < ($snpcount-1); $i+=1)
    {
	my $recomI = sprintf("%.15f",$recomrate[$i]/100);      ## as should be a prob
	print OUT2 "$posvec[$i] $recomI\n";
    }
    print OUT2 "$posvec[($snpcount-1)] 0\n";
    
    if($warnings){ print "There were warnings; see ${outfilePRE}.warnings\n";}
    
    close OUT2;
    
} # end !$omitrecom

close OUT;
close OUT3;
