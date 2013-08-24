#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: spconversion.pl
#
#        USAGE: spconversion.pl config/config.pm
#
#  DESCRIPTION:GMM 频谱转换
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Li Xian (), shysian@gmail.com
# ORGANIZATION: USTC
#      VERSION: 1.0
#      CREATED: 2013年01月08日 12时59分09秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my ($init, $mcep, $doTrain, $doTest);
$init    = 0;
$mcep    = 0;
$doTrain = 1;
$doTest  = 0;

my $numMix=128;
my $emotion = 'sad';
my $prjDir  = '..';
my $workDir = "$prjDir/vc/to$emotion"."_$numMix"."mix";
my $testDir = "$workDir/test";

###load config
our @testfiles;
require("$prjDir/config/config.pm");
####

#SPTK path
my $SPTK = "/usr/local/SPTK/bin";

if ($init)
{
    mkdir "$workDir",                      0755;
    mkdir "$workDir/list",                 0755;
    mkdir "$workDir/wav",                  0755;
    mkdir "$workDir/wav/neutral",          0755;
    mkdir "$workDir/wav/$emotion",         0755;
    mkdir "$testDir",                      0755;
    mkdir "$testDir/wav",                  0755;
    mkdir "$testDir/wav/neutral",          0755;
    mkdir "$testDir/wav/$emotion",         0755;
    mkdir "$testDir/wav/neutral-$emotion"."_$numMix"."mix", 0755;
    my @neutralWavs = <$prjDir/neutral/wav/*.wav>;
    my @emotionWavs = <$prjDir/$emotion/wav/*.wav>;

    for my $wavFile (@neutralWavs)
    {

        system("cp $wavFile $workDir/wav/neutral");
    }
    for my $iFile (@testfiles)
    {

        #print "test file $iFile\n";
        system(
            "mv $workDir/wav/neutral/$iFile.wav $testDir/wav/neutral/$iFile.wav"
        );
    }

    for my $wavFile (@emotionWavs)
    {
        system("cp $wavFile $workDir/wav/$emotion");
    }
    for my $iFile (@testfiles)
    {
        system(
            "mv $workDir/wav/$emotion/$iFile.wav $testDir/wav/$emotion/$iFile.wav"
        );
    }
    open NEUTRAL_LIST, ">$workDir/list/neutral_tr.list"
      or die "can't open $!\n";
    my @trainNeutral = <$workDir/wav/neutral/*.wav>;
    for my $wavFile (@trainNeutral)
    {
        $wavFile =~ /([0-9]+).wav/;
        print NEUTRAL_LIST "neutral/$1\n";
    }
    open EMOTION_LIST, ">$workDir/list/$emotion"."_$numMix"."mix_tr.list"
      or die "can't open $!\n";
    my @trainEmotion = <$workDir/wav/$emotion/*.wav>;
    for my $wavFile (@trainEmotion)
    {
        $wavFile =~ /([0-9]+).wav/;
        print EMOTION_LIST "$emotion/$1\n";
    }

}

if ($mcep)
{
    mkdir "$workDir/mcepc0",          0755;
    mkdir "$workDir/mcepc0/neutral",  0755;
    mkdir "$workDir/mcepc0/$emotion", 0755;
    my @neutralSPs = <$prjDir/neutral/spectrum/*.sp>;
    for my $neutralSP (@neutralSPs)
    {
        $neutralSP =~ /([0-9]+).sp/;
        print "extacting mcep from $neutralSP\n";
        system(
            "$SPTK/x2x +af $neutralSP|$SPTK/mcep -a 0.42 -m 24 -l 1024 -q 3 |x2x +fd>$workDir/mcepc0/neutral/$1.mcep"
        );
    }
    my @emotionSPs = <$prjDir/$emotion/spectrum/*.sp>;
    for my $emotionSP (@emotionSPs)
    {
        print "extracting mcep from $emotionSP\n";
        $emotionSP =~ /([0-9]+).sp/;
        system(
            "$SPTK/x2x +af $emotionSP|$SPTK/mcep -a 0.42 -m 24 -l 1024 -q 3 |x2x +fd>$workDir/mcepc0/$emotion/$1.mcep"
        );
    }

}

 for my $iFile (@testfiles)
    {
        system(
            "cp $workDir/mcepc0/neutral/$iFile.mcep $testDir/wav/neutral-$emotion"."_$numMix"."mix/$iFile.wav.org.mcep"
        );
    }

if ($doTrain)
{
    mkdir "$workDir/mcep",          0755;
    mkdir "$workDir/mcep/neutral",  0755;
    mkdir "$workDir/mcep/$emotion", 0755;
    mkdir "$workDir/npow",          0755;
    mkdir "$workDir/npow/neutral",  0755;
    mkdir "$workDir/npow/$emotion", 0755;
    system("$prjDir/vc/scripts/VCTrain train neutral $emotion"."_$numMix"."mix $numMix");
}

if ($doTest)
{
  for my $iFile (@testfiles)
    {

        #print "test file $iFile\n";
        system(
            "mv $workDir/wav/neutral/$iFile.wav $testDir/wav/neutral/$iFile.wav"
        );
    }

   
    system("$prjDir/vc/scripts/VCTrain test neutral $emotion"."_$numMix"."mix");

    for my $iFile (@testfiles)
    {
        system(
            "$SPTK/x2x +df $testDir/wav/neutral-$emotion"."_$numMix"."mix/$iFile.wav.conv.mcep|$SPTK/mgc2sp -a 0.42 -m 24 -l 1024 -o 2 |$SPTK/x2x +fa513> $testDir/wav/neutral-$emotion"."_$numMix"."mix/$iFile.sp "
        );
    }
}
