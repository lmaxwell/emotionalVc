#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: sp_gmm.pl
#
#        USAGE: ./sp_gmm.pl  ../config/config.pm  trgEmotion numMix mode
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
#      CREATED: 2013年06月08日 15时20分52秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;


if($#ARGV<3)
{
	print "USAGE: ./sp_gmm.pl  ../config/config.pm  trgEmotion numMix mode";
	exit;
}


my ($init,$init_mcep, $mcep, $doTrain, $doTest,$mode);
$mode=$ARGV[3];

if($mode eq 'init')
{
$init    = 1;
$init_mcep=0;
$mcep    = 0;
$doTrain = 0;
$doTest  = 0;
}
elsif($mode eq 'train')
{

$init    = 0;
$init_mcep=0;
$mcep    = 0;
$doTrain = 1;
$doTest  = 0;

}
elsif($mode eq 'test')
{

$init    = 0;
$init_mcep=0;
$mcep    = 1;
$doTrain = 1;
$doTest  = 1;

}

my $numMix=$ARGV[2];
my $src = 'neutral';
my $trg = $ARGV[1];
my $prjDir  = '..';
my $workDir = "$prjDir/vc/train/to$trg"."_$numMix"."mix";
my $testDir = "$workDir/test";

###load config
our @testfiles;
require("$prjDir/config/config.pm");
####

#SPTK path
my $SPTK = "/usr/local/SPTK/bin";
my $SPTK36="/usr/local/SPTK-3.6/bin";
if ($init)
{
	mkdir "$prjDir/vc/train",0755;
    mkdir "$workDir",                      0755;
    mkdir "$workDir/list",                 0755;
    mkdir "$workDir/wav",                  0755;
    mkdir "$workDir/wav/$src",          0755;
    mkdir "$workDir/wav/$trg",         0755;
    mkdir "$testDir",                      0755;
    mkdir "$testDir/wav",                  0755;
    mkdir "$testDir/wav/$src",          0755;
    mkdir "$testDir/wav/$trg",         0755;
    mkdir "$testDir/wav/$src-$trg"."_$numMix"."mix", 0755;
    my @srcWavs = <$prjDir/$src/wav/*.wav>;
    my @trgWavs = <$prjDir/$trg/wav/*.wav>;

    for my $wavFile (@srcWavs)
    {

        system("cp $wavFile $workDir/wav/$src");
    }
    for my $iFile (@testfiles)
    {

        #print "test file $iFile\n";
        system(
            "mv $workDir/wav/$src/$iFile.wav $testDir/wav/$src/$iFile.wav"
        );
    }

    for my $wavFile (@trgWavs)
    {
        system("cp $wavFile $workDir/wav/$trg");
    }
    for my $iFile (@testfiles)
    {
        system(
            "mv $workDir/wav/$trg/$iFile.wav $testDir/wav/$trg/$iFile.wav"
        );
    }
    open SRC_LIST, ">$workDir/list/$src"."_tr.list"
      or die "can't open $!\n";
    my @trainSrc = <$workDir/wav/$src/*.wav>;
    for my $wavFile (@trainSrc)
    {
        $wavFile =~ /([0-9]+).wav/;
        print SRC_LIST "$src/$1\n";
    }
    open TRG_LIST, ">$workDir/list/$trg"."_$numMix"."mix_tr.list"
      or die "can't open $!\n";
    my @trainTrg = <$workDir/wav/$trg/*.wav>;
    for my $wavFile (@trainTrg)
    {
        $wavFile =~ /([0-9]+).wav/;
        print TRG_LIST "$trg/$1\n";
    }

}

if ($init_mcep)
{
    mkdir "$prjDir/vc/train/mcepc0",          0755;
    mkdir "$prjDir/vc/train/mcepc0/$src",  0755;
    mkdir "$prjDir/vc/train/mcepc0/$trg", 0755;
    my @srcSPs = <$prjDir/$src/spectrum/*.sp>;
    for my $srcSP (@srcSPs)
    {
        $srcSP =~ /([0-9]+).sp/;
        print "extacting mcep from $srcSP\n";
        system(
            "$SPTK/x2x +af $srcSP|$SPTK/mcep -a 0.42 -m 24 -l 1024 -q 4 |x2x +fd>$prjDir/vc/train/mcepc0/$src/$1.mcep"
        );
    }
    my @trgSPs = <$prjDir/$trg/spectrum/*.sp>;
    for my $trgSP (@trgSPs)
    {
        print "extracting mcep from $trgSP\n";
        $trgSP =~ /([0-9]+).sp/;
        system(
            "$SPTK/x2x +af $trgSP|$SPTK/mcep -a 0.42 -m 24 -l 1024 -q 4 |x2x +fd>$prjDir/vc/train/mcepc0/$trg/$1.mcep"
        );
    }

}

if($mcep)
{
	mkdir "$workDir/mcepc0",0755;
	mkdir "$workDir/mcepc0/$src",0755;
	mkdir "$workDir/mcepc0/$trg",0755;
	my @mcepFiles=<$prjDir/vc/train/mcepc0/$src/*.mcep>;
	for my $mcepFile(@mcepFiles)
	{
		my $base=`basename $mcepFile`;
		chomp $base;
		#print "copy".$mcepFile."\n";
		system("cp $mcepFile $workDir/mcepc0/$src/$base");
	}

	@mcepFiles=<$prjDir/vc/train/mcepc0/$trg/*.mcep>;
	
	for my $mcepFile(@mcepFiles)
	{
		my $base=`basename $mcepFile`;
		chomp $base;
		#print "copy".$mcepFile."\n";
		system("cp $mcepFile $workDir/mcepc0/$trg/$base");
  	}
    
 	for my $iFile (@testfiles)
    {
        system(
            "mv $workDir/mcepc0/$src/$iFile.mcep $testDir/wav/$src-$trg"."_$numMix"."mix/$iFile.wav.org.mcep"
        );
    }
}
if ($doTrain)
{
    mkdir "$workDir/mcep",          0755;
    mkdir "$workDir/mcep/$src",  0755;
    mkdir "$workDir/mcep/$trg", 0755;
    mkdir "$workDir/npow",          0755;
    mkdir "$workDir/npow/$src",  0755;
    mkdir "$workDir/npow/$trg", 0755;
    system("$prjDir/vc/scripts/VCTrain train $src $trg"."_$numMix"."mix $numMix");
}

if ($doTest)
{
  for my $iFile (@testfiles)
    {

        #print "test file $iFile\n";
        system(
            "mv $workDir/wav/$src/$iFile.wav $testDir/wav/$src/$iFile.wav"
        );
    }
 	for my $iFile (@testfiles)
    {
        system(
            "mv $workDir/mcepc0/$src/$iFile.mcep $testDir/wav/$src-$trg"."_$numMix"."mix/$iFile.wav.org.mcep"
        );
    }

   
    system("$prjDir/vc/scripts/VCTrain test $src $trg"."_$numMix"."mix");

    for my $iFile (@testfiles)
    {
        system(
            "$SPTK/x2x +df $testDir/wav/$src-$trg"."_$numMix"."mix/$iFile.wav.conv.mcep|$SPTK/mgc2sp -a 0.42 -m 24 -l 1024 -o 3 |$SPTK/x2x +fa513> $testDir/wav/$src-$trg"."_$numMix"."mix/$iFile.sp "
        );
    }
	
	my $melCDdir="$testDir/melCD";
	my $testSrc="$melCDdir/test$src.mat";
	my $testTrg="$melCDdir/testConv.mat";
	mkdir $melCDdir,0755;
	system("rm -f $melCDdir/conv.data");
	#system("rm -f $melCDdir/trg.data");
	for my $iFile (@testfiles)
	{
		if($iFile eq '341')
		{
			next;
		}

		my $convFile="$testDir/wav/$src-$trg"."_$numMix"."mix/$iFile.wav.conv.mcep";
		my $trgFile="$workDir/../mcepc0/$trg/$iFile.mcep";
		system("$SPTK36/x2x +df  $convFile >$convFile.f");
		system("$SPTK36/x2x +df  $trgFile >$trgFile.f");
		my $convTemp="$melCDdir/convtemp";
		my $trgTemp="$melCDdir/trgtemp";
		system("$SPTK36/bcp +f -l 25 -s 1 -e 24 -L 24 -S 0 $convFile.f>$convTemp ");
		system("$SPTK36/bcp +f -l 25 -s 1 -e 24 -L 24 -S 0 $trgFile.f>$trgTemp ");
		my $nancheck=`$SPTK/nan $convTemp|wc -l`;
		chomp $nancheck;
		if($nancheck)
		{
			next;
		}

		my $nancheck=`$SPTK/nan $trgTemp|wc -l`;
		chomp $nancheck;
		if($nancheck)
		{
			next;
		}
		system("$SPTK36/dtw -l 24 -p 2 $convTemp $trgTemp>>$melCDdir/conv.data");

	}

	system("$SPTK36/bcp +f -l 48 -s 0 -e 23 -L 25 -S 1 $melCDdir/conv.data > $testSrc");
	system("$SPTK36/bcp +f -l 48 -s 24 -e 47 -L 25 -S 1 $melCDdir/conv.data > $testTrg");
	#system("$SPTK36/cdist -o 0 -m 24 $melCDdir/$src.mat $melCDdir/$trg.mat |x2x +fa >$melCDdir/melCD ");
	system("$SPTK36/cdist -o 0 -m 24 $testSrc $testTrg|x2x +fa >$melCDdir/melCD_test ");




}

