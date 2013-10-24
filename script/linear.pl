#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: linear.pl
#
#        USAGE: ./linear.pl  ../config/config.pm
#
#  DESCRIPTION: 高斯归一F0变换
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Li Xian (), shysian@gmail.com
# ORGANIZATION: USTC
#      VERSION: 1.0
#      CREATED: 2013年01月10日 17时10分31秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my $prjDir='..';
our @testfiles;
require("$prjDir/config/config.pm");

##流程控制
my $init=1;
my $extract=1;
my $prepareForTrain=1;
my $doTrain=1;

my $doTest=0;
##全局参数初始化
my $emotion='angry';
my $trainDir="$prjDir/train/linear/$emotion";
my $testDir="$prjDir/test/linear/$emotion";
#######
#SPTK PATH
my $SPTK="/usr/local/SPTK/bin";
#MATLAB PATH
my $MATLAB ="/usr/local/matlab/bin/matlab -nodisplay -nosplash -nojvm";
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



if($init)
{
	print "initialize...\n";
	mkdir "$prjDir/neutral/dct",0755;
	mkdir "$prjDir/neutral/dct/linear",0755;
	mkdir "$prjDir/$emotion/dct/linear",0755;
	mkdir "$prjDir/$emotion/dct",0755;
    mkdir "$prjDir/train",         0755;
    mkdir "$prjDir/train/linear", 0755;
    mkdir "$prjDir/train/linear/$emotion", 0755;
	mkdir "$trainDir",0755;
	mkdir "$prjDir/test",0755;
	mkdir "$prjDir/test/linear",0755;
	mkdir "$prjDir/test/linear/$emotion",0755;
	mkdir "$testDir",0755;
}


#extract F0 of voiced segment
if($extract)
{
	my @pdir=("$prjDir/neutral","$prjDir/$emotion");
	for my $dir(@pdir)
	{
		my @labfiles=<$dir/lab/mlevel/*.lab>;
		for(@labfiles)
		{
			my $base=`basename $_ .lab`;
			chomp $base;
			my $f0file="$dir/logf0/$base.f0";
			my $outf0file="$dir/dct/linear/$base.f0";
			print "Extracting voiced segment from $f0file to $outf0file\n";
			open LAB,"<$_" or die "can't open $!\n";
			open F0,"<$f0file" or die "can't open $!\n";
			open OUTF0,">$outf0file" or die "can't open $!\n";
			my @f0data=<F0>;
			chomp @f0data;
			while(my $line=<LAB>)
			{
				 my @text = split /\s+/, $line;
                push(@text, 0);

                #	print "$text[3]\n";
                my $word = $text[3];
                my $initialType =
                  get_initialType($word, \@dict, \%dict_hash, $dict_len);

                if ($initialType == 1)
                {

					#print OUTF0 "$text[0] $text[1] $text[3]\n";
                    my $start = $text[0] / 50000;
                    my $end   = $text[1] / 50000 ;
                    my $i     = $start;
                    while ($i < $end)
                    {
                        print OUTF0 "$f0data[$i]\n";
                        $i++;

					}
			
				}
			 elsif ($initialType == 2)
                {
                    my $start = $text[0] / 50000;
                    $line = <LAB>;
                    @text = split /\s+/, $line;
                    my $end = $text[1] / 50000 ;
                    my $i   = $start;
                    while ($i < $end)
                    {
                        print OUTF0 "$f0data[$i]\n";
                        $i++;
                    }


                }
                elsif ($initialType == 0)
                {
                    $line = <LAB>;
                    @text = split /\s+/, $line;
					#print OUTF0 "$text[0] $text[1] $word\n";
                    my $start = $text[0] / 50000;
                    my $end   = $text[1] / 50000 ;
                    my $i     = $start;
                    while ($i < $end)
                    {
                        print OUTF0 "$f0data[$i]\n";
                        $i++;
                    }
                }
                else
                {

                    print "may have ERROR $word\n"
                      unless ($text[2] eq 'sil' or $word eq 'sp');
                }
	
			}
		}
	}
}

#准备训练数据
if($prepareForTrain)
{
    my %testfiles;
    for (@testfiles)
    {
        $testfiles{$_} = 1;
    }

	my @neutralF0Files=<$prjDir/neutral/dct/linear/*.f0>;
	my @emotionF0Files=<$prjDir/$emotion/dct/linear/*.f0>;
	printf "prepare data for train\n";
	open NEUTRAL,">$trainDir/neutral.f0" or die "can't open $!\n";
	for my $neutralF0File(@neutralF0Files)
	{
		my $base=`basename $neutralF0File .f0`;
		$base=~s/\n//;
		next if(exists $testfiles{$base});
		open F0,"<$neutralF0File" or die "can't open $!\n";
		while(<F0>)
		{
			print NEUTRAL $_;
		}
		close F0;
	}
	close NEUTRAL;
	open EMOTION,">$trainDir/$emotion.f0" or die "can't open $!\n";
	for my $emotionF0File(@emotionF0Files)
	{
		my $base=`basename $emotionF0File .f0`;
		$base=~s/\n//;
		next if(exists $testfiles{$base});
		
		open F0,"<$emotionF0File" or die "can't open $!\n";
		while(<F0>)
		{
			print EMOTION $_;
		}
		close F0;
	}
	close EMOTION;
}

#训练
if($doTrain)
{
	print ("doint train\n");
	system("$SPTK/x2x +af $trainDir/neutral.f0|$SPTK/vstat -l 1|x2x +fa>$trainDir/neutral.stat");
	system("$SPTK/x2x +af $trainDir/$emotion.f0|$SPTK/vstat -l 1|x2x +fa>$trainDir/$emotion.stat");
	print("done!\n");
}
if($doTest)
{
	print ("Test...\n");
	for my $iFile(@testfiles)
	{
		system("cp $prjDir/neutral/logf0/$iFile.f0 $testDir");
		system("mv $testDir/$iFile.f0 $testDir/$iFile.org.f0");
	}
	open NSTAT,"<$trainDir/neutral.stat" or die "can't open $!\n";
	my ($nMean,$nVariance)=<NSTAT>;
	chomp ($nMean,$nVariance);
	print "mean of neutral f0:$nMean variance of neutral f0:$nVariance\n";
	open ESTAT,"<$trainDir/$emotion.stat" or die "can't open $!\n";
	my ($eMean,$eVariance)=<ESTAT>;
	chomp ($eMean,$eVariance);
	print "mean of $emotion f0:$eMean variance of $emotion f0 $eVariance\n";
	
	for my $iFile(@testfiles)
	{
		system("$SPTK/x2x +af $testDir/$iFile.org.f0|$SPTK/sopr -s $nMean -d sqrt$nVariance -m sqrt$eVariance -a $eMean|$SPTK/x2x +fa > $testDir/$iFile.f0");
	}
	
	print "generate resynth matlab script\n";
	open M,">$testDir/resynth.m" or die "can't open $!\n";
	select M;
	print "addpath(\'$prjDir/script\');\n";
	print "emotion=\'$emotion\';\n";
	print "method=\'linear\';\n";
	for my $testFile(@testfiles)
	{
		print "iFile=$testFile;\n";
		print "f_morph(iFile,emotion,method,1,0,0,\'$testDir\');\n";
		print "f_morph(iFile,emotion,method,1,0,1,\'$testDir\');\n";
	}
	close M;
	select STDOUT;
	system("$MATLAB<$testDir/resynth.m");

}
#############################################################
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

