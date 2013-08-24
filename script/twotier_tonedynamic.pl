#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: twotier_tonedynamic.pl
#
#        USAGE: ./twotier_tonedynamic.pl  ../config/config.pm 
#
#  DESCRIPTION: 加入声调曲线第一个DCT系数的速度约束
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Li Xian (), shysian@gmail.com
# ORGANIZATION: USTC
#      VERSION: 1.0
#      CREATED: 2013年01月11日 13时17分09秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my ($clean,$emotion,$init,$extractNeutral, $extract1, $extract2, $prepareForTrain,
    $doTrain,$doTest);

my $prjDir="..";
#
#load config
our @testfiles;
require("$prjDir/config/config.pm");

#控制流程
$clean=0;
$init=1;

$extractNeutral  = 1;
$extract1        = 1;    #提取phrase
$extract2        = 1;    #提取tone

$prepareForTrain = 1;    #准备训练数据
$doTrain         = 1;    #训练
$doTest 		 = 1;	 #测试
#全局参数初始化
my $method='tone_d1';
$emotion = 'sad';
my $doF0Normalization=1;
my $dctNumOfPhrase=3;
my $dctNum = 8;
my %mixNum=(
	"sad" => [9,24], 
	"angry" =>[4,4],
	"happy" => [4,4],

);


my $testDir="$prjDir/test/$method/$emotion/$mixNum{$emotion}[0]mix$mixNum{$emotion}[1]mix_$dctNumOfPhrase";
my $trainDir="$prjDir/train/$method/$emotion/$mixNum{$emotion}[0]mix$mixNum{$emotion}[1]mix_$dctNumOfPhrase";

=begin  BlockComment  # BlockCommentNo_1

my @testfiles = (
                 201, 221, 241, 261, 281, 301, 321, 341,
                 361, 381, 401, 421, 441, 461, 481
                );

=end    BlockComment  # BlockCommentNo_1

=cut


#SPTK TOOL
my $SPTK="/usr/local/SPTK/bin";
#MATLAB PATH
my $MATLAB ="/usr/local/matlab/bin/matlab -nodisplay -nosplash -nojvm";

#gmm tool
my $VC="$prjDir/vc/src";
my $lbg="$VC/vq/lbg";
my $gmmtrain="$VC/gmm/gmm_jde";
my $gmmpara="$VC/gmmpara/gmm_para";
my $gmmmap="$VC/gmmmap/gmmmap";


####f0 mean
my $f0mean;
my $f0statFile= "$prjDir/train/linear/$emotion/neutral.stat" ;
if(-r $f0statFile )
{
	open STAT,"<$f0statFile" or die "can't open $!\n";
	my @stats=<STAT>;
	chomp @stats;
	$f0mean=$stats[0];
	print "f0 mean  is $f0mean\n";
}
else
{
	print "no $f0statFile\n";
	exit;
}

#######
#read dict
open DICT, "<$prjDir/neutral/forcedAlignment/dict" or die "can't open";
my $dict;
$dict = do { local $/; <DICT> };

#全部词典
my @dict;
@dict = $dict =~ /[a-zA-Z]+/g;
my $dict_len = $#dict + 1;

#半元音词典
my %dict_hash;
my $i = 0;
while ($i < $dict_len)
{
    if ($dict[$i] eq 'NULL')
    {
        $dict_hash{$dict[$i - 1]} = $dict[$i+1];
    }
    $i++;
}

my @pdir = ("$prjDir/$emotion");
push @pdir, "$prjDir/neutral" if ($extractNeutral);

##################clean###############################
if($clean)
{
	system("rm -R $prjDir/$emotion/dct/$method  $testDir $trainDir");
}
###################initialize#########################
if($init)
{
	print "initialize...\n";
    mkdir "$prjDir/train",         0755;
    mkdir "$prjDir/train/$method", 0755;
    mkdir "$prjDir/train/$method/$emotion", 0755;
	mkdir "$trainDir",0755;
	mkdir "$prjDir/test",0755;
	mkdir "$prjDir/test/$method",0755;
	mkdir "$prjDir/test/$method/$emotion",0755;
	mkdir "$testDir",0755;
}


##################提取phrase dct########################
if ($extract1)
{
    for my $dir (@pdir)
    {
        open PB, "<$prjDir/neutral/forcedAlignment/newphrasebreak" or die "can't open";
        mkdir "$dir/dct/$method", 0755;

        while (<PB>)
        {
            /\"\*\/([0-9]+\.lab)\"/;
            my $fname = $1;
            print "current file is $fname\n";
            my $line1 = <PB>;
            my @words = split ' ', $line1;
            print "words are @words\n";
            my $base = `basename $fname .lab`;
            $base =~ s/\n//;
            open LAB, "<$dir/lab/word/$fname" or die "can't open";
            <LAB>;    #skip first line
            open F0, "<$dir/logf0/$base.f0" or die "can't open";
            my @f0data = <F0>;
            chomp @f0data;
            close F0;
			my @norm_f0data=map {$_-$f0mean} @f0data;
			@f0data=@norm_f0data if ($doF0Normalization);

            open PHRASE, ">$dir/dct/$method/$base.phrase"
              or die "can't open $!\n";
		   open NUM,">$dir/dct/$method/$base.num" or die "can't open $!\n";
            my ($pp, $cp, $np);
            $pp = 1;
            my $i = 1;
			my $numPhrase=0;
            for (@words)
            {
                if (/\/p/ or /\/w/)
                {
					$numPhrase++;
                    $cp = $i;
                    my $wordNumInPhrase = $cp - $pp + 1;
					print NUM "$wordNumInPhrase\n";
                    my $line            = <LAB>;
                    $line = <LAB> if ($line =~ /sp/);
                    my @startTime = split ' ', $line;
                    while ($wordNumInPhrase - 1 > 0)
                    {
                        $line = <LAB>;
                        next if ($line =~ /sp/);
                        $wordNumInPhrase--;
                    }
                    my @endTime = split ' ', $line;
                    print
                      "current phrase start:$startTime[0] ;end:$endTime[1]\n";
                    print "extract phrase...\n";
                    my @f0OfCurPhrase;
                    my $index = 0;
                    for (my $j = $startTime[0] ;
                         $j <= $endTime[1] - 50000 ;
                         $j += 50000)
                    {

                        $f0OfCurPhrase[$index] = $f0data[$j / 50000];
                        $index++;
                    }
                    print "length of current f0 :$#f0OfCurPhrase+1\n";
                    print "doint TONE...\n";
                    my $dctsize = $#f0OfCurPhrase + 1;
                    my $dctlen  = $dctNumOfPhrase-1;                     #不含直流分量
                    my @dctOfCurPhrase =
                      `echo @f0OfCurPhrase |x2x +af |dct -l $dctsize |sopr -d sqrt$dctsize|bcut -s 0 -e $dctlen|x2x +fa`;
                    chomp @dctOfCurPhrase;
                    print "dctofCurPhrase:@dctOfCurPhrase\n";
                    print PHRASE "$startTime[0] $endTime[1] ";
                    $wordNumInPhrase = $cp - $pp + 1;

                    for ($pp .. $cp)
                    {
                        print PHRASE "$words[$_-1]";
                    }

					print PHRASE " $numPhrase $wordNumInPhrase ";
                    print PHRASE "\n";
                    print PHRASE "$_\n" for (@dctOfCurPhrase);
                    print PHRASE "#\n";
                    my @curPhraseContour =
                      `echo @dctOfCurPhrase|x2x +af |idct -l $dctsize|sopr -m sqrt$dctsize|x2x +fa`;
                    chomp @curPhraseContour;
                    print "Current phrase contour:@curPhraseContour\n";
                    print "subtract phrase contour from f0\n";

                    for (my $j = $startTime[0] ;
                         $j <= $endTime[1] - 50000 ;
                         $j += 50000)
                    {

                        $f0data[$j / 50000] -=
                          $curPhraseContour[$j / 50000 - $startTime[0] / 50000];
                    }

                    $pp = $cp + 1;
                }
                $i++;

            }    #end for
			print NUM "$numPhrase\n";
            close LAB;
            close PHRASE;
			close NUM;
            print "f0 after subtracting phrase contour:@f0data\n";

            #减去phrase后的f0
            open LAB, "<$dir/lab/mlevel/$fname" or die "can't open $!\n";
            my $line2 = <LAB>;    #skip line 1
            open OUTF0, ">$dir/dct/$method/$base.f0" or die "can't open $!\n";
            while ($line2 = <LAB>)
            {
                my @text = split /\s+/, $line2;
                push(@text, 0);

                #	print "$text[3]\n";
                my $word = $text[3];
                my $initialType =
                  get_initialType($word, \@dict, \%dict_hash, $dict_len);

                if ($initialType == 1)
                {

                    print OUTF0 "$text[0] $text[1] $text[3]\n";
                    my $start = $text[0] / 50000;
                    my $end   = $text[1] / 50000 ;
                    my $i     = $start;
                    while ($i < $end)
                    {
                        print OUTF0 "$f0data[$i]\n";
                        $i++;
                    }
                    print OUTF0 ".\n";

                }
                elsif ($initialType == 2)
                {
                    print OUTF0 $text[0];
                    my $start = $text[0] / 50000;
                    $line2 = <LAB>;
                    @text = split /\s+/, $line2;
                    print OUTF0 " $text[1] $word\n";
                    my $end = $text[1] / 50000 ;
                    my $i   = $start;
                    while ($i < $end)
                    {
                        print OUTF0 "$f0data[$i]\n";
                        $i++;
                    }

                    print OUTF0 ".\n";

                }
                elsif ($initialType == 0)
                {
                    $line2 = <LAB>;
                    @text = split /\s+/, $line2;
                    print OUTF0 "$text[0] $text[1] $word\n";
                    my $start = $text[0] / 50000;
                    my $end   = $text[1] / 50000 ;
                    my $i     = $start;
                    while ($i < $end)
                    {
                        print OUTF0 "$f0data[$i]\n";
                        $i++;
                    }
                    print OUTF0 ".\n";
                }
                else
                {

                    print "may have ERROR $word\n"
                      unless ($text[2] eq 'sil' or $word eq 'sp');
                }
            }    #end while
        }    #end while
    }    #end for
}    #end if

############################提取tone dct###############################
if ($extract2)
{

    for my $dir (@pdir)
    {

        my @files = <$dir/dct/$method/*.f0>;
        foreach (@files)
        {

            my $name = `basename $_ .f0`;
            chomp($name);
            print "extracting tone dct :$name\n";
            my $dctfile = $name . 'dct';
            open F0,   "<$dir/dct/$method/$name.f0"   or die "can't open";
            open TONE, ">$dir/dct/$method/$name.tone" or die "can't open";
            my $line;
			my @alllines;
			my @alldct;
            while ($line = <F0>)
            {

                if ($line =~ /[a-z]+/)
                {

                    open TEMP, ">$dir/dct/$method/temp" or die "can't open";
					#print TONE $line;
					push @alllines,$line;
                    my $start;
                    my $end;

                    $line =~ /([0-9]+)\s([0-9]+)/;
                    $start = $1;

                    #print "start is $start\n";
                    $end = $2;

                    #print "end is $end\n";
                    while (($line = <F0>) !~ /^\.$/)
                    {
                        chomp $line;
                        print TEMP $line . "\n";
                    }
                    close(TEMP);
                    my $length = $end / 50000 - $start / 50000;
                    my @dct =
                      `x2x +af $dir/dct/$method/temp|dct -l $length|sopr -d sqrt$length|x2x +fa`;

                    #print "current length is ".$#dct."+1\n";
                    while ($#dct < $dctNum)
                    {
                        push @dct, "0\n";
                    }
                    @dct = splice @dct, 0, $dctNum;
					@alldct=(@alldct,@dct);
                    #print "@dct\n";

					#print TONE @dct;
					#print TONE "#\n";
                    system("rm -f $dir/$method/dct/temp");
                }    #end if
            }    #end while
			
			my $T=($#alldct+1)/$dctNum;
			my $dim=$dctNum;
			my @win=(-0.5,0,0.5);
			my @alldct_d;
			for(my $t=0;$t<$T;$t++)
			{
				for (my $j=0;$j<$dim;$j++)
				{
					$alldct_d[$t*($dim+1)+$j]=$alldct[$t*$dim+$j];
				}
				if($t>0 and $t<$T-1)
				{
					$alldct_d[$t*($dim+1)+$dim]=$win[0]*$alldct[($t-1)*$dim]+$win[2]*$alldct[($t+1)*$dim];
				}
				elsif($t==0)
				{
					$alldct_d[$t*($dim+1)+$dim]=$win[0]*$alldct[($t)*$dim]+$win[2]*$alldct[($t+1)*$dim];
				}
				elsif($t==$T-1)
				{	
					$alldct_d[$t*($dim+1)+$dim]=$win[0]*$alldct[($t-1)*$dim]+$win[2]*$alldct[($t)*$dim];
				}
			}
			$dim=$dim+1;
			my $i=0;
			while($i<$T)
			{
				my @temp=splice @alldct_d,0,$dim;
				my $temp_line=shift @alllines;
				print TONE $temp_line;
				printf TONE "%.7f\n",$_ for @temp;
				print TONE "#\n";
				$i++;
			}




            close(TONE);
            close(F0);
			


        }    #end foreach
    }    #end for
}    #end if

########################prepare data for train#######################################
if ($prepareForTrain)
{
    my %testfiles;
    for (@testfiles)
    {
        $testfiles{$_} = 1;
    }
    for my $iFile (201 .. 500)
    {
        print "$iFile is testfile\n" if (exists $testfiles{$iFile});
    }
    my @paraKind = ("phrase", "tone");
    for my $paraKind (@paraKind)
    {
        my @emotions = ("neutral", $emotion);
        for my $emotionKind (@emotions)
        {
            my $dir = "$prjDir/neutral";
            $dir = "$prjDir/$emotion" if ($emotionKind eq "$emotion");
            my @paraFiles = <$dir/dct/$method/*.$paraKind>;
            open TRAINDATA, ">$trainDir/$emotionKind.$paraKind"
              or die "can't open";
            foreach (@paraFiles)
            {
                my $base = `basename $_ .$paraKind`;
                $base =~ s/\n//;
                next if (exists $testfiles{$base});
                open PARA, "<$_" or die "can't open";
                local $/ = "#\n";
                while (<PARA>)
                {
                    chomp;
                    my @paras = split /\n/, $_;
                    print "@paras\n";
                    shift @paras;
                    print TRAINDATA "@paras\n";
                }    #end while

            }    #end foreach

        }    #end for
    }    #end for
}    #end if


#my $numMixOfPhrase=5;
#my $numMixOfTone=8;

#######################训练###################################
if ($doTrain)
{
	for my $paraKind("phrase","tone")
	{
		for my $emotionKind("neutral","$emotion")
		{
			system("$SPTK/x2x +af $trainDir/$emotionKind.$paraKind >$trainDir/$emotionKind.$paraKind.f");
		}
	}
	my $dctAugment=$dctNum+1;
	system("$SPTK/merge -s 0 -l $dctNumOfPhrase -L $dctNumOfPhrase +f $trainDir/neutral.phrase.f $trainDir/$emotion.phrase.f>$trainDir/neutral_$emotion.phrase.f");
	system("rm -f $trainDir/neutral.phrase.f $trainDir/$emotion.phrase.f");
	system("$SPTK/merge -s 0 -l $dctAugment -L $dctAugment +f $trainDir/neutral.tone.f $trainDir/$emotion.tone.f>$trainDir/neutral_$emotion.tone.f");
	system("rm -f $trainDir/neutral.tone.f $trainDir/$emotion.tone.f");
	my $vectorLength=2*$dctNumOfPhrase;
	system("$SPTK/gmm -l $vectorLength -m $mixNum{$emotion}[0] -f $trainDir/neutral_$emotion.phrase.f >$trainDir/neutral_$emotion.phrase.gmm.f");
	system("$SPTK/gmm -l 18 -m $mixNum{$emotion}[1] -f $trainDir/neutral_$emotion.tone.f >$trainDir/neutral_$emotion.tone.gmm.f");
	system("$SPTK/x2x +fa $trainDir/neutral_$emotion.phrase.gmm.f>$trainDir/neutral_$emotion.phrase.gmm");
	system("$SPTK/x2x +fa $trainDir/neutral_$emotion.tone.gmm.f>$trainDir/neutral_$emotion.tone.gmm");

}    #end if

=a
if ($doTrain)
{
	for my $paraKind("phrase","tone")
	{
		for my $emotionKind("neutral","$emotion")
		{
			system("$SPTK/x2x +ad $trainDir/$emotionKind.$paraKind >$trainDir/$emotionKind.$paraKind.d");
		}
	}
	
	system("$SPTK/merge -s 0 -l 3 -L 3 +d $trainDir/neutral.phrase.d $trainDir/$emotion.phrase.d>$trainDir/neutral_$emotion.phrase.d");
	system("rm -f $trainDir/neutral.phrase.d $trainDir/$emotion.phrase.d");
	system("$SPTK/merge -s 0 -l $dctNum -L $dctNum +d $trainDir/neutral.tone.d $trainDir/$emotion.tone.d>$trainDir/neutral_$emotion.tone.d");
	system("rm -f $trainDir/neutral.tone.d $trainDir/$emotion.tone.d");
	system("$lbg -dim 6 -cls 8 $trainDir/neutral_$emotion.phrase.d $trainDir/vqmean");
}    #end if
=cut
#########################测试###############################
if($doTest)
{

	for my $iFile(@testfiles)
	{
		system("cp $prjDir/neutral/dct/$method/$iFile.phrase $testDir");
		system("cp $prjDir/neutral/dct/$method/$iFile.tone $testDir");
	}
	print "generate Matlab file for converting...\n";
	open M,">$testDir/convert.m" or die "can't open $!\n";
	print M "addpath('$prjDir/script');\n";
	my @testPhraseFiles=<$testDir/*.phrase>;
	for my $testPhraseFile(@testPhraseFiles)
	{
		print M "fprintf(\'converting $testPhraseFile \\n \');\n";
		print M "fd=fopen(\'$testPhraseFile.predict\',\'w\');\n";
		local $/="#\n";
		open PHRASE,"<$testPhraseFile" or die "can't open $!\n";
		while(<PHRASE>)
		{
			chomp;
			my @phraseParas=split /\n/,$_;
			#print "@phraseParas\n";
			my $text=shift @phraseParas;	
			print M "fprintf(fd,\'$text\\n\');\n";
			print M "x=[@phraseParas]';\n";
			my $vectorLength=$dctNumOfPhrase*2;
			print M "predict=gmm_convert(x,$vectorLength,$mixNum{$emotion}[0],\'$trainDir/neutral_$emotion.phrase.gmm\');\n";
			print M "fprintf(fd,\'%f\\n\',predict);\n";
			print M "fprintf(fd,\'#\\n\');\n";	
		}
		close PHRASE;
	}#end for
	my @testToneFiles=<$testDir/*.tone>;
	for my $testToneFile(@testToneFiles)
	{
		print M "fprintf(\'converting $testToneFile \\n \');\n";
		print M "fd=fopen(\'$testToneFile.predict\',\'w\');\n";
		local $/="#\n";
		open TONE,"<$testToneFile" or die "can't open $!\n";
		my @alltext;
		my @alltone;
		my $T;
		while(<TONE>)
		{
			chomp;
			my @toneParas=split /\n/,$_;
			#print "@toneParas\n";
			my $text=shift @toneParas;
			@alltext=(@alltext,$text);
			@alltone=(@alltone,@toneParas);
			$T++;
		}
		
		print M "x=[@alltone]';\n";
		print M "predict=gmm_convert(x,18,$mixNum{$emotion}[1],\'$trainDir/neutral_$emotion.tone.gmm\',\'trajectory\');\n";
		for my $t(0..$T-1)
		{
			print M "fprintf(fd,\'$alltext[$t]\\n\');\n";
			print M "fprintf(fd,\'%f\\n\',predict($t*8+1:($t+1)*8));\n";
			print M "fprintf(fd,\'#\\n\');\n";	
		}
		print M "fclose(fd);\n";
		close TONE;
	}#end for
	print "doing convert...\n";
	close M;
	system("$MATLAB<$testDir/convert.m");

	print "resynth f0 from converted parameter...\n";

	my @predictToneFiles=<$testDir/*.tone.predict>;
	for my $predictToneFile(@predictToneFiles)
	{
		open TONE,"<$predictToneFile" or die "can't open $!\n";
		open TEMP,">$testDir/tonetemp" or die "can't open $!\n";
		local $/="#\n";
		while(<TONE>)
		{
			chomp;
			my @toneParas=split /\n/,$_;
			my $text=shift @toneParas;
			$text=~/([0-9]+)\s([0-9]+)\s/;
			my $lenOfSeg=($2-$1)/50000;
			my $bcutEnd=$lenOfSeg-1;
			my @f0OfSeg=`echo @toneParas|$SPTK/x2x +af|$SPTK/bcut +f -s 0 -e $bcutEnd|$SPTK/idct -l $lenOfSeg|$SPTK/sopr -m sqrt$lenOfSeg|$SPTK/x2x +fa`;
			print TEMP "$text\n";
			print TEMP @f0OfSeg;
			print TEMP "#\n";
		}
		close TONE;
		close TEMP;
		my $fname=`basename $predictToneFile .tone.predict`;
		$fname=~s/\n//;
		interpolate("$testDir/tonetemp","$testDir/$fname.tone.f0");
	}
	## add phrase tier
	my @predictPhraseFiles=<$testDir/*.phrase.predict>;
	for my $predictPhraseFile(@predictPhraseFiles)
	{
		open PHRASE,"<$predictPhraseFile" or die "can't open $!\n";
		open TEMP,">$testDir/phrasetemp" or die "can't open $!\n";
		local $/="#\n";
		while(<PHRASE>)
		{
			chomp;
			my @phraseParas=split /\n/,$_;
			my $text = shift @phraseParas;
			$text=~/([0-9]+)\s([0-9]+)/;
			my $lenOfPhrase=($2-$1)/50000;
			my $bcutEnd=$lenOfPhrase-1;

			local $/="\n";################
			my @f0OfPhrase=`echo @phraseParas|$SPTK/x2x +af|$SPTK/bcut +f -s 0 -e $bcutEnd|$SPTK/idct -l $lenOfPhrase|$SPTK/sopr -m sqrt$lenOfPhrase|$SPTK/x2x +fa1`;
			print TEMP "$text\n";
			#print "@f0OfPhrase";
			print "length:$#f0OfPhrase+1\n";
			chomp @f0OfPhrase;
			local $/="#\n";#################
			#print "befor unnorm:@f0OfPhrase\n";
			if($doF0Normalization)
			{
				my @f0_unnorm=map {$_+$f0mean} @f0OfPhrase;
				#print "after unnorm:@f0_unnorm\n";
				@f0OfPhrase=@f0_unnorm;
			}
			print TEMP "$_\n" for(@f0OfPhrase);
			print TEMP "#\n";
		}
		close PHRASE;
		close TEMP;
		my $fname=`basename $predictPhraseFile .phrase.predict`;
		$fname=~s/\n//;
		interpolate("$testDir/phrasetemp","$testDir/$fname.phrase.f0");
	}
	print "generate resynth matlab script\n";
	open M,">$testDir/resynth.m" or die "can't open $!\n";
	select M;
	print "addpath(\'$prjDir/script\');\n";
	print "emotion=\'$emotion\';\n";
	print "method=\'$method\';\n";
	for my $testFile(@testfiles)
	{
		print "iFile=$testFile;\n";
		print "f_morph(iFile,emotion,method,1,0,0,\'$testDir\');\n";
		print "f_morph(iFile,emotion,method,1,0,1,\'$testDir\');\n";
	}
	close M;
	select STDOUT;
	system("$MATLAB<$testDir/resynth.m");

}#end if


############################type 1:元音，半元音；type 2:浊辅音；type 0:清辅音;type -1:词典中无法找到;type: 3:sp(short pause)######################
sub get_initialType
{
    my ($word, $p_dict, $p_dict_hash, $dict_len) = @_;
    my $i;
    $i = 0;
    while ($i < $dict_len and $p_dict->[$i] ne $word)
    {
        $i++;
    }
    if ($word eq 'sp')
    {
        return 3;
    }
    if ($i >= $dict_len)
    {
        return (-1);
    }
    if (exists $p_dict_hash->{$word})
    {
        return (1);
    }
    elsif (   ($p_dict->[$i + 1] eq 'm')
           or ($p_dict->[$i + 1] eq 'n')
           or ($p_dict->[$i + 1] eq 'l')
           or ($p_dict->[$i + 1] eq 'r'))
    {
        return (2);
    }
    else
    {
        return (0);
    }

}

sub interpolate
{
	my ($f0temp,$fname)=@_;
	open TEMP,"<$f0temp" or die "can't oepn";
	open F0,">$fname" or die "can't open";
	my $start=0;
	my $end=0;
	my $nextStart;
	my $line;
	local $/="#\n";
	print "interpolate $fname\n";
	while($line=<TEMP>)
	{
		chomp($line);
		my @f0=split /\n/,$line;
		my $text=shift @f0;
		$text=~/([0-9]+)\s([0-9]+)/;
		$end=$1/50000;
		$nextStart=$2/50000;
		#	print $line."\n";
		my $vlnum=$end-$start;
		while($vlnum>0)
		{
			print F0 "0\n";
			$vlnum-=1;
		}
		print F0 "$_\n" for (@f0);
		$start=$nextStart;
	}
	close TEMP;
	close F0;
}
