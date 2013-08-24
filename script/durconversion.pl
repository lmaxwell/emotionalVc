#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: durvonversion.pl
#
#        USAGE: ./durvonversion.pl  dumpfeat.pl
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Li Xian (), shysian@gmail.com
# ORGANIZATION: USTC
#      VERSION: 1.0
#      CREATED: 2013年01月11日 15时47分48秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my $prjDir = "..";

our %phoneType;
our %initials;
require("$prjDir/script/dumpfeat.pl");

#
#load config

our @testfiles;
require("$prjDir/config/config.pm");

my ($clean, $init, $extract, $prepareForTrain, $doTrainAndTest);

#控制流程
$clean           = 0;
$init            = 1;
$extract         = 1;
$prepareForTrain = 1;
$doTrainAndTest  = 1;    #训练&测试

#speech_tools path
my $speech_tools = "~/speech_tools/bin";
#SPTK tools
my $SPTK="/usr/local/SPTK/bin";
#全局参数初始化
my $emotion        = 'sad';
my $minClusterSize = 10;
my $dataDir        = "$prjDir/$emotion/dur";
my $testDir        = "$prjDir/test/dur/$emotion";
my $trainDir       = "$prjDir/train/dur/$emotion";

#声明函数 
sub isInitial();



if ($init)
{
    mkdir "$dataDir",          0755;
    mkdir "$prjDir/test/dur",  0755;
    mkdir "$prjDir/train/dur", 0755;
    mkdir "$trainDir",         0755;
    mkdir "$testDir",          0755;

	mkdir "$trainDir/F0",0755;
	mkdir "$testDir/F0",0755;

	mkdir "$trainDir/scale",0755;
	mkdir "$testDir/scale",0755;
}

my @pdir = ("$prjDir/neutral", "$prjDir/$emotion");
if ($extract)
{
    open SPDATA, ">$dataDir/sp" or die "can't open $!\n";
    my @labfiles = <$pdir[1]/lab/phone/*.lab>;
    for (@labfiles)
    {
        open ELAB, "<$_" or die "can't open $!\n";
        my $base = `basename $_ .lab`;
        chomp $base;
        open NLAB, "<$pdir[0]/lab/phone/$base.lab" or die "can't open $!\n";
        my @emotionlabs = <ELAB>;
        chomp @emotionlabs;

        #print "@emotionlabs\n";
        my (@start, @end, @phone);
        my $i = 0;
        for my $emotionlab (@emotionlabs)
        {
            ($start[$i], $end[$i], $phone[$i]) = split ' ', $emotionlab;
            $i++;
        }
        my $numOfPhone = $i;

        my @neutrallabs = <NLAB>;
        chomp @neutrallabs;
        my @nLen;

        for my $neutrallab (@neutrallabs)
        {
            my ($s, $e, $p) = split ' ', $neutrallab;
            next if ($p eq 'sp' or $p eq 'sil');
            my $len = ($e - $s) / 50000;

            push @nLen, $len;
        }

        open DATA, ">$dataDir/$base.dur" or die "can't open $!\n";

        my $j = 0;
        for my $i (1 .. $numOfPhone - 2)
        {

            my $eLen = ($end[$i] - $start[$i]) / 50000;
			my $position;
			if($i<=4)
			{
				$position=0;
			}
			elsif($i>=$numOfPhone-4)
			{
				$position=2;
			}
			else
			{
				$position=1;
			}
            if ($phone[$i] eq 'sp')
            {

                print SPDATA
                  "$eLen\tsp\t$phoneType{$phone[$i-1]}\t$phoneType{$phone[$i+1]}\n";
                next;
            }
            my $ratio = $eLen / $nLen[$j];
            if ($phone[$i - 1] eq 'sp' and $phone[$i + 1] eq 'sp')
            {

                print DATA
                  "$ratio\t$nLen[$j]\t$phoneType{$phone[$i]}\t$phoneType{$phone[$i-2]}\t$phoneType{$phone[$i+2]}\t".&isInitial($phone[$i])."\t$position\n";
                $j++;
                next;
            }
            if ($phone[$i - 1] eq 'sp')
            {

                print DATA
                  "$ratio\t$nLen[$j]\t$phoneType{$phone[$i]}\t$phoneType{$phone[$i-2]}\t$phoneType{$phone[$i+1]}\t".&isInitial($phone[$i])."\t$position\n";
                $j++;
                next;
            }

            if ($phone[$i + 1] eq 'sp')
            {

                print DATA
                  "$ratio\t$nLen[$j]\t$phoneType{$phone[$i]}\t$phoneType{$phone[$i-1]}\t$phoneType{$phone[$i+2]}\t".&isInitial($phone[$i])."\t$position\n";
                $j++;
                next;
            }
            print DATA
              "$ratio\t$nLen[$j]\t$phoneType{$phone[$i]}\t$phoneType{$phone[$i-1]}\t$phoneType{$phone[$i+1]}\t".&isInitial($phone[$i])."\t$position\n";
            $j++;
        }    #end for
        close DATA;

    }    #end for
    close SPDATA;

}
if ($prepareForTrain)
{
    open DUR, ">$trainDir/dur.feat" or die "can't open $!\n";
    my %testfiles;
    $testfiles{$_} = 0 for (@testfiles);
    my @datafiles = <$dataDir/*.dur>;
    for (@datafiles)
    {
        my $base = `basename $_ .dur`;
        chomp $base;
        next if (exists $testfiles{$base});
        open DATA, "<$_" or die "can't open $!\n";
        while (<DATA>)
        {
            print DUR $_;
        }
        close DATA;
    }
	system("awk '{printf(\$1\"\\t\");printf(\$2\"\\n\")}' $trainDir/dur.feat>$trainDir/F0/dur.feat");
	system("awk '{printf(\$2\"\\n\")}' $trainDir/dur.feat>$trainDir/scale/neutral.dur");
	system("awk '{printf(\$1\"\\n\")}' $trainDir/dur.feat>$trainDir/scale/$emotion.ratio");
	system("$SPTK/x2x +af $trainDir/scale/neutral.dur>$trainDir/scale/neutral.dur.f");
	system("$SPTK/x2x +af $trainDir/scale/$emotion.ratio>$trainDir/scale/$emotion.ratio.f");
	system("$SPTK/vopr -m $trainDir/scale/$emotion.ratio.f $trainDir/scale/neutral.dur.f>$trainDir/scale/$emotion.dur.f ");
	system("$SPTK/x2x +fa $trainDir/scale/$emotion.dur.f>$trainDir/scale/$emotion.dur");
    close DUR;
    open DESC, ">$trainDir/dur.desc" or die "can't open $!\n";
    print DESC
      "((ratio float)\n(inputDur float)\n(c.phonType m d t s a f n l sil)\n(p.phoneType m d t s a f n l sil)\n(n.phoneType m d t s a f n l sil)\n(isInitial 1 0)\n(positionInSentence 0 1 2)\n)";
	close DESC;
	
	open DESC,">$trainDir/F0/dur.desc" or die "can't open $!\n";
	print DESC "((ratio float)\n(inputDur float)\n)";
	close DESC;


	mkdir "$trainDir/F1",0755;
	system("awk '{printf(\$1\"\\t\");printf(\$2\"\\t\");printf(\$3\"\\t\");printf(\$4\"\\t\");printf(\$5\"\\t\");printf(\$6\"\\n\")}' $trainDir/dur.feat>$trainDir/F1/dur.feat");
	open DESC,">$trainDir/F1/dur.desc" or die "can't open $!\n";
    print DESC
      "((ratio float)\n(inputDur float)\n(c.phonType m d t s a f n l sil)\n(p.phoneType m d t s a f n l sil)\n(n.phoneType m d t s a f n l sil)\n(isInitial 1 0)\n)";
	close DESC;
}
	
if ($doTrainAndTest)
{

	mkdir "$testDir/F1",0755;
    open TEST, ">$testDir/dur.test" or die "can't open $!\n";
	for my $iFile(@testfiles)
	{
		print "testfile:$iFile\n";
		system("cp $dataDir/$iFile.dur $testDir/F1/$iFile.dur");
		
		#system("awk '{printf(\$2\"\\n\")}' $testDir/F1/$_.temp>$testDir/F1/$_.dur");
		#system("rm $testDir/F1/$_.temp");
		open DATA, "<$dataDir/$iFile.dur" or die "can't open $!\n";
		print TEST $_ while (<DATA>);
		close DATA;
	}
	close TEST;
=a
	for (@testfiles)
	{
		print "testfile:$_\n";
		#system("cp $dataDir/$_.dur $testDir/F1/$_.dur");
	}
???????????????????!
=cut 

	system("awk '{print \$2}' $testDir/dur.test>$testDir/inputDur ");
	system("awk '{print \$1}' $testDir/dur.test>$testDir/trueRatio");
	system("$SPTK/x2x +af $testDir/inputDur>$testDir/inputDur.f");
	system("$SPTK/sopr -m 5 $testDir/inputDur.f>$testDir/inputDur.ms.f");
	system("$SPTK/x2x +af $testDir/trueRatio>$testDir/trueRatio.f");
	system("$SPTK/vopr -m $testDir/inputDur.f $testDir/trueRatio.f>$testDir/trueDur.f");
	system("$SPTK/sopr -m 5 $testDir/trueDur.f>$testDir/trueDur.ms.f");
	system("$SPTK/rmse -l 0 $testDir/trueDur.ms.f $testDir/inputDur.f|$SPTK/x2x +fa>$testDir/rmse_ne");

	system("rm -f $testDir/rmse");
    while ($minClusterSize < 200)
    {
        system(
            "$speech_tools/wagon -data $trainDir/dur.feat -desc $trainDir/dur.desc -stop $minClusterSize -output $trainDir/dur.tree_$minClusterSize"
        );
        system(
            "$speech_tools/wagon_test -data $testDir/dur.test -desc $trainDir/dur.desc -tree $trainDir/dur.tree_$minClusterSize -o $testDir/predict_$minClusterSize -predict_val"
        );
		system("$SPTK/x2x +af $testDir/predict_$minClusterSize |$SPTK/vopr -m $testDir/inputDur.f>$testDir/outputDur.f");
		system("$SPTK/sopr -m 5 $testDir/outputDur.f>$testDir/outputDur.ms.f");
		system("$SPTK/rmse -l 0 $testDir/trueDur.ms.f $testDir/outputDur.ms.f|$SPTK/x2x +fa>>$testDir/rmse");
		system("rm -f $testDir/outputDur.f $testDir/outputDur.ms.f");
        $minClusterSize += 10;
    }
	


   #F1 feature group

	system("awk '{printf(\$1\"\\t\");printf(\$2\"\\t\");printf(\$3\"\\t\");printf(\$4\"\\t\");printf(\$5\"\\t\");printf(\$6\"\\n\")}' $testDir/dur.test>$testDir/F1/dur.test");
	$minClusterSize=10;

	system("rm -f $testDir/F1/rmse");
	while ($minClusterSize < 200)
    {
        system(
            "$speech_tools/wagon -data $trainDir/F1/dur.feat -desc $trainDir/F1/dur.desc -stop $minClusterSize -output $trainDir/F1/dur.tree_$minClusterSize"
        );
        system(
            "$speech_tools/wagon_test -data $testDir/F1/dur.test -desc $trainDir/F1/dur.desc -tree $trainDir/F1/dur.tree_$minClusterSize -o $testDir/F1/predict_$minClusterSize -predict_val"
        );
		system("$SPTK/x2x +af $testDir/F1/predict_$minClusterSize |$SPTK/vopr -m $testDir/inputDur.f>$testDir/F1/outputDur.f");
		system("$SPTK/sopr -m 5 $testDir/F1/outputDur.f>$testDir/F1/outputDur.ms.f");
		
		system("$SPTK/rmse -l 0 $testDir/trueDur.ms.f $testDir/F1/outputDur.ms.f|$SPTK/x2x +fa>>$testDir/F1/rmse");
		system("rm -f $testDir/F1/outputDur.f $testDir/F1/outputDur.ms.f");
        $minClusterSize += 10;
    }











	for my $iFile(@testfiles)
	{
		print $iFile."\n";
		my $ClusterSize=20;
 		system(
            "$speech_tools/wagon_test -data $testDir/F1/$iFile.dur -desc $trainDir/dur.desc -tree $trainDir/dur.tree_$ClusterSize -o $testDir/F1/$iFile.predict_$ClusterSize -predict_val"
        );

	}

	
	#0 Feature group
	system("awk '{printf(\$1\"\\t\");printf(\$2\"\\n\")}' $testDir/dur.test>$testDir/F0/dur.test");
	$minClusterSize=10;

	system("rm -f $testDir/F0/rmse");
	while ($minClusterSize < 200)
    {
        system(
            "$speech_tools/wagon -data $trainDir/F0/dur.feat -desc $trainDir/F0/dur.desc -stop $minClusterSize -output $trainDir/F0/dur.tree_$minClusterSize"
        );
        system(
            "$speech_tools/wagon_test -data $testDir/F0/dur.test -desc $trainDir/F0/dur.desc -tree $trainDir/F0/dur.tree_$minClusterSize -o $testDir/F0/predict_$minClusterSize -predict_val"
        );
		system("$SPTK/x2x +af $testDir/F0/predict_$minClusterSize |$SPTK/vopr -m $testDir/inputDur.f>$testDir/F0/outputDur.f");
		system("$SPTK/sopr -m 5 $testDir/F0/outputDur.f>$testDir/F0/outputDur.ms.f");
		
		system("$SPTK/rmse -l 0 $testDir/trueDur.ms.f $testDir/F0/outputDur.ms.f|$SPTK/x2x +fa>>$testDir/F0/rmse");
		system("rm -f $testDir/F0/outputDur.f $testDir/F0/outputDur.ms.f");
        $minClusterSize += 10;
    }

	#scale dur
	system("awk '{printf(\$2\"\\n\")}' $testDir/dur.test>$testDir/scale/dur.test");
	system("$SPTK/x2x +af $testDir/scale/dur.test >$testDir/scale/dur.test.f");
	system("$SPTK/vstat -l 1 -o 1 $trainDir/scale/neutral.dur.f>$testDir/scale/neutral.mean.f");
	system("$SPTK/x2x +fa $testDir/scale/neutral.mean.f>$testDir/scale/neutral.mean");
	system("$SPTK/vstat -l 1 -o 1 $trainDir/scale/$emotion.dur.f>$testDir/scale/$emotion.mean.f");
	system("$SPTK/x2x +fa $testDir/scale/$emotion.mean.f>$testDir/scale/$emotion.mean");
	open NMEAN,"<$testDir/scale/neutral.mean" or die "can't open";
	my $ndurMean=<NMEAN>;
	chomp $ndurMean;

	open EMEAN,"<$testDir/scale/$emotion.mean" or die "can't open";
	my $edurMean=<EMEAN>;
	chomp $edurMean;
	system("$SPTK/sopr -d $ndurMean -m $edurMean $testDir/scale/dur.test.f >$testDir/scale/dur.out.f ");
	system("$SPTK/sopr -m 5 $testDir/scale/dur.out.f > $testDir/scale/dur.out.ms.f");
	system("$SPTK/rmse -l 0 $testDir/trueDur.ms.f $testDir/scale/dur.out.ms.f|$SPTK/x2x +fa >$testDir/scale/rmse");

	


}



##################3
sub isInitial
{

	my $phone=shift;
	if (exists $initials{$phone})
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
