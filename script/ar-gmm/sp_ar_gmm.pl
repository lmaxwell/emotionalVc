#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: sp_ar_gmm.pl
#
#        USAGE: ./sp_ar_gmm.pl  ../config/config.pm  targetEmotion   numMix   order  mode
#
#  DESCRIPTION: convert spectrum using  ar-gmm
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Li Xian (), shysian@gmail.com
# ORGANIZATION: USTC
#      VERSION: 1.0
#      CREATED: 2013年06月07日 21时53分45秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;


my ($init,$init_mcep, $mcep,$preTrain, $doTrain, $doTest);

if($#ARGV<4)
{
	print "USAGE: ./sp_ar_gmm.pl  ../config/config.pm  targetEmotion   numMix   order  mode\n";
	exit;
}
my $mode = $ARGV[4];

if($mode eq 'init')
{
	$init    = 1;
	$init_mcep=1;
	$mcep    = 1;
	$preTrain = 0;
	$doTrain = 0;
	$doTest  = 0;
}
if($mode eq 'prepare')
{
	$init    = 1;
	$init_mcep=0;
	$mcep    = 0;
	$preTrain = 1;
	$doTrain = 0;
	$doTest  = 0;
}
if($mode eq 'train')
{
	$init    = 0;
	$init_mcep=0;
	$mcep    = 0;
	$preTrain = 0;
	$doTrain = 1;
	$doTest  = 0;
}
if($mode eq 'test')
{
	$init    = 0;
	$init_mcep=0;
	$mcep    = 0;
	$preTrain = 0;
	$doTrain = 0;
	$doTest  = 1;
}
my $numMix=$ARGV[2];
my $src = 'neutral';
my $trg = $ARGV[1];
my $prjDir  = '..';
my $workDir = "$prjDir/vc/train/to$trg"."_ar";
my $testDir = "$workDir/test/$numMix-mix";

my $order=$ARGV[3];
#train file list
my $srcList="$workDir/list/$src"."_tr.list";
my $trgList="$workDir/list/$trg"."_tr.list";

###load config
our @testfiles;
require("$prjDir/config/config.pm");
####

#SPTK path
my $SPTK = "/usr/local/SPTK/bin";
my $SPTK36="/usr/local/SPTK-3.6/bin";
#MATLAB PATH
my $MATLAB ="/usr/local/matlab/bin/matlab   -nodesktop -nosplash ";#-nodisplay -nosplash ";#-nojvm";
#toda
my $todaBin="$prjDir/vc/src";
my $todaScript="$prjDir/vc/scripts/GMMMAP";

if ($init)
{
	mkdir "$prjDir/vc/train",0755;
    mkdir "$workDir",                      0755;
    mkdir "$workDir/list",                 0755;
    mkdir "$workDir/wav",                  0755;
    mkdir "$workDir/wav/$src",          0755;
    mkdir "$workDir/wav/$trg",         0755;
	mkdir "$workDir/test",				0755;
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
	close SRC_LIST;
    open TRG_LIST, ">$workDir/list/$trg"."_tr.list"
      or die "can't open $!\n";
    my @trainTrg = <$workDir/wav/$trg/*.wav>;
    for my $wavFile (@trainTrg)
    {
        $wavFile =~ /([0-9]+).wav/;
        print TRG_LIST "$trg/$1\n";
    }
	close TRG_LIST;

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
            "$SPTK/x2x +af $srcSP|$SPTK/mcep -a 0.42 -m 24 -l 1024 -q 3 |x2x +fd>$prjDir/vc/train/mcepc0/$src/$1.mcep"
        );
    }
    my @trgSPs = <$prjDir/$trg/spectrum/*.sp>;
    for my $trgSP (@trgSPs)
    {
        print "extracting mcep from $trgSP\n";
        $trgSP =~ /([0-9]+).sp/;
        system(
            "$SPTK/x2x +af $trgSP|$SPTK/mcep -a 0.42 -m 24 -l 1024 -q 3 |x2x +fd>$prjDir/vc/train/mcepc0/$trg/$1.mcep"
        );
    }

}

if($mcep)
{
#	mkdir "$workDir/mcepc0",0755;               # 
#	mkdir "$workDir/mcepc0/$src",0755;          # 
#	mkdir "$workDir/mcepc0/$trg",0755;          # 
 	mkdir "$workDir/mcep",          0755;
    mkdir "$workDir/mcep/$src",  0755;
    mkdir "$workDir/mcep/$trg", 0755;
    mkdir "$workDir/npow",          0755;
    mkdir "$workDir/npow/$src",  0755;
    mkdir "$workDir/npow/$trg", 0755;

	my @mcepFiles=<$prjDir/vc/train/mcepc0/$src/*.mcep>;
	for my $mcepFile(@mcepFiles)
	{
		my $base=`basename $mcepFile .mcep`;
		chomp $base;
		#print "copy".$mcepFile."\n";
		#system("cp $mcepFile $workDir/mcepc0/$src/$base");

		system("$todaBin/extdim/extdim -nmsg -dim 25 -sd 1 $mcepFile  $workDir/mcep/$src/$base.mcep");
		system("$todaBin/extdim/extdim -nmsg -dim 25 -ed 0 $mcepFile  $workDir/npow/$src/$base.npow");
	}

	@mcepFiles=<$prjDir/vc/train/mcepc0/$trg/*.mcep>;
	
	for my $mcepFile(@mcepFiles)
	{
		my $base=`basename $mcepFile .mcep`;
		chomp $base;
		#print "copy".$mcepFile."\n";
		#system("cp $mcepFile $workDir/mcepc0/$trg/$base");
		
		system("$todaBin/extdim/extdim -nmsg -dim 25 -sd 1 $mcepFile  $workDir/mcep/$trg/$base.mcep");
		system("$todaBin/extdim/extdim -nmsg -dim 25 -ed 0 $mcepFile  $workDir/npow/$trg/$base.npow");

  	}
    
 	for my $iFile (@testfiles)
    {
        system(
            "cp $prjDir/vc/train/mcepc0/$src/$iFile.mcep $testDir/wav/$src-$trg"."_$numMix"."mix/$iFile.wav.org.mcep"
        );
    }
}
if ($preTrain)
{

	my $seg_mark_src_File="$workDir/dtw/$numMix-mix/seg_mark_src";
	my @seg_mark_src;
	for my $i(1..24)
	{
		push @seg_mark_src,0;
	}
	system("echo @seg_mark_src |x2x +ad >$seg_mark_src_File");


	open LIST,"<$srcList" or die "can't open";
	my $trainFile;
	system("rm -f $workDir/dtw/$numMix-mix/$src.mat");
	while($trainFile=<LIST>)
	{
	
		chomp $trainFile;
		my $npowf="$workDir/npow/$trainFile.npow";
		my $mcepf="$workDir/mcep/$trainFile.mcep";
		my $exmcepf="$workDir/mcep/$trainFile.exmcep";
		system("$todaBin/extfrm/extfrm -nmsg -dim 24 -lp -9.0 -npowfile $npowf $mcepf $exmcepf");
		system("$SPTK/x2x +df $npowf >$npowf.f ");
		system("cat $seg_mark_src_File >> $workDir/dtw/$numMix-mix/$src.mat");
		system("cat $exmcepf >> $workDir/dtw/$numMix-mix/$src.mat");
	
	}
	system("$SPTK/x2x +da24 $workDir/dtw/$numMix-mix/$src.mat > $workDir/dtw/$numMix-mix/$src.mat.txt");
	system("$SPTK/x2x +df $workDir/dtw/$numMix-mix/$src.mat > $workDir/dtw/$numMix-mix/$src.mat.f");
	close LIST;
	system("$SPTK36/vstat -l 24 -o 1  $workDir/dtw/$numMix-mix/$src.mat.f >$workDir/dtw/$numMix-mix/mcepMean ");
	my @mcepMean=`x2x +fa $workDir/dtw/$numMix-mix/mcepMean `;
	#print "@mcepMean";
	#exit;
	open LIST,"<$trgList" or die "can't open";

	system("rm -f $workDir/dtw/$numMix-mix/$trg.mat");
	while($trainFile=<LIST>)
	{
		chomp $trainFile;
		my $npowf="$workDir/npow/$trainFile.npow";
		my $mcepf="$workDir/mcep/$trainFile.mcep";
		my $exmcepf="$workDir/mcep/$trainFile.exmcep";
		system("$todaBin/extfrm/extfrm -nmsg -dim 24 -lp -100.0 -npowfile $npowf $mcepf $exmcepf");
		system("$SPTK/x2x +df $npowf >$npowf.f ");
		system("cat $seg_mark_src_File >> $workDir/dtw/$numMix-mix/$trg.mat");
		system("cat $mcepf >> $workDir/dtw/$numMix-mix/$trg.mat");
	
	}
	system("$SPTK/x2x +da24 $workDir/dtw/$numMix-mix/$trg.mat > $workDir/dtw/$numMix-mix/$trg.mat.txt");
	close LIST;


	#exit;#!!!!!!!!!!!!!!!!!!!!!!

	mkdir "$workDir/dtw",0755;
	mkdir "$workDir/dtw/$numMix-mix",0755;
	mkdir "$workDir/dtw/$numMix-mix/$src",0755;
	#system("rm -f $workDir/dtw/$numMix-mix/$src/*");
	my @trainFiles=`cat $trgList`;
	chomp @trainFiles;
	#print $_ for (@trainFiles);
	my $index=0;
	my $jointFile="$workDir/dtw/$numMix-mix/$src-$trg.mat";
    my $gvFile="$workDir/dtw/$numMix-mix/$src-$trg.gvtemp";
    my $gvModel="$workDir/ar-gmm/$numMix-mix/gvModel";
	my $seg_markFile="$workDir/dtw/$numMix-mix/seg_mark";
	my @seg_mark;
	for my $i(1..48)
	{
		push @seg_mark,0;
	}
	system("echo @seg_mark |x2x +af >$seg_markFile");
	
	system("rm -f $jointFile.melCD");
	for my $trainFile(`cat $srcList`)
	{
		chomp $trainFile;
		my $srcFile="$workDir/mcep/$trainFile.mcep";
		my $trgFile="$workDir/mcep/$trainFiles[$index].mcep";
		my $srcPow="$workDir/npow/$trainFile.npow";
		my $trgPow="$workDir/npow/$trainFiles[$index].npow";
		my $twfFile="$workDir/dtw/$numMix-mix/$trainFile.twf";
		my $twfdata="$workDir/dtw/$numMix-mix/$trainFile.twfdat";
		my $frmcdFile="$workDir/dtw/$numMix-mix/$trainFile.frmcd";
		my $viterbiFile="$workDir/dtw/$numMix-mix/$trainFile.viterbi";
		my $scoreFile="$workDir/dtw/$numMix-mix/$trainFile.score";
		print "$srcFile-$trgFile\n";

		my $dtwFile="$workDir/dtw/$numMix-mix/$trainFile.jmcep";
		system("$SPTK/x2x +df $srcFile >$srcFile.f");
		system("$SPTK/x2x +df $trgFile >$trgFile.f");
		system("$SPTK36/dtw -l 24 -p 2 $trgFile.f $srcFile.f -v $viterbiFile -s $scoreFile >>$dtwFile");
			system("$prjDir/script/ar-gmm/ar-gmmC/extract -l 24 -o $order -t -100.0  $viterbiFile $srcFile.f $trgFile.f $srcPow.f $trgPow.f $dtwFile");
		system("cat $scoreFile >>$jointFile.melCD") unless(-z $scoreFile);
#				
#		for my $i(0..$order)
#		{
#			system("$SPTK36/vopr -l 24 -s  $dtwFile-$i  $workDir/dtw/$numMix-mix/mcepMean > $dtwFile-norm");
#			system("mv $dtwFile-norm $dtwFile-$i");
#		}

		#	system("$todaBin/dtw/dtw -nmsg -fixed -dim 24 -ldim 23 -frmcdfile $frmcdFile -twffile $twfFile -twfdat $twfdata  $srcFile $trgFile");
		
		#my $dtwFile="$workDir/dtw/$numMix-mix/$trainFile.jmcep";
		#system("$todaBin/dtw/dtw -nmsg -fixed -dim 24 -ldim 23 -intwf $twfFile -dtwotfile $dtwFile $srcFile $trgFile");
		#system("rm -f $srcFile $trgFile $twfFile");

		if ($index==0)
		{
		   #system("mv $frmcdFile $workDir/dtw/$numMix-mix/frmcd.txt");
			system("cat $dtwFile-0 > $jointFile");
			for my $i (0..$order)
			{
			system("cat $seg_markFile $dtwFile-$i > $jointFile.$i-ar");
			}
			system("cat $seg_markFile $jointFile > $jointFile-ar");
			system("$SPTK36/vstat -l 48 -d -o 2 $dtwFile-0 > $gvFile");
		}
		else
		{
			
			#	system("mv $workDir/dtw/$numMix-mix/frmcd.txt $workDir/dtw/$numMix-mix/frmcd.txt.tmp");
			#system("cat $workDir/dtw/$numMix-mix/frmcd.txt.tmp $frmcdFile > $workDir/dtw/$numMix-mix/frmcd.txt");
			#system("rm -f $workDir/dtw/$numMix-mix/frmcd.txt.tmp $frmcdFile");

			
			system("mv $jointFile $jointFile.temp");
			system("cat $jointFile.temp $dtwFile-0 > $jointFile");

			system("mv $jointFile-ar $jointFile-ar.temp");
			system("cat $jointFile-ar.temp $seg_markFile > $jointFile-ar");
			system("mv $jointFile-ar $jointFile-ar.temp");
			system("cat $jointFile-ar.temp $dtwFile > $jointFile-ar");
			

			for my $i (0..$order)
			{
				system("mv $jointFile.$i-ar $jointFile.$i-ar.temp");
				system("cat $jointFile.$i-ar.temp $seg_markFile > $jointFile.$i-ar");
				system("mv $jointFile.$i-ar $jointFile.$i-ar.temp");
				system("cat $jointFile.$i-ar.temp $dtwFile-$i > $jointFile.$i-ar");
				#	system("rm -f $jointFile.tmp $jointFile-ar.tmp $dtwFile");
			}
			system("$SPTK36/vstat -l 48 -d -o 2 $dtwFile-0 >> $gvFile");
		}

		$index+=1;

	}
	system("$SPTK36/vstat $jointFile.melCD |x2x +fa >>$jointFile.melCDavg ");
	system("$SPTK/x2x +fa48 $jointFile-ar >$jointFile-ar.txt");
	system("$SPTK/x2x +fd $jointFile >$jointFile.f");
	#system("$SPTK36/vstart -l 48 -o 1  $jointFile > $workDir/dtw/$numMix-mix/mcepVMean");
	 
	for my $i (0..$order)
	{

		system("$SPTK/x2x +fa48 $jointFile.$i-ar >$jointFile-ar.txt-$i");
	}

	mkdir "$workDir/ar-gmm",0755;
	mkdir "$workDir/ar-gmm/$numMix-mix",0755;
	system("$SPTK36/vstat -d -l 48 -o 0 $gvFile > $gvModel");
	system("$SPTK/gmm -m 64 -l 48 -b 5 $jointFile >$workDir/ar-gmm/$numMix-mix/initmodel");
}
if($doTrain)
{
	my $jointFile="$workDir/dtw/$numMix-mix/$src-$trg.mat";
    system("$prjDir/script/ar-gmm/ar-gmmC/ar-gmm -m 64 -l 48 -o $order  -t 0.01 -v 0.001 -T 1 -N 20 -I $workDir/ar-gmm/$numMix-mix/initmodel $jointFile $workDir/ar-gmm/$numMix-mix/$src-$trg-order$order.ar-gmm > $workDir/ar-gmm/$numMix-mix/order$order.log");

	exit;
	mkdir "$workDir/cbook",0755;
	mkdir "$workDir/cbook/$numMix-mix",0755;
	mkdir "$workDir/gmm_init",0755;
	mkdir "$workDir/gmm_init/$numMix-mix",0755;
	mkdir "$workDir/gmm_jde",0755;
	mkdir "$workDir/gmm_jde/$numMix-mix",0755;
	mkdir "$workDir/gmm_para",0755;
	mkdir "$workDir/gmm_para/$numMix-mix",0755;
	system("$todaBin/vq/lbg -sd 0 -dim 48 -cls $numMix $jointFile.f $workDir/cbook/$numMix-mix/$src-$trg.cb>$workDir/cbook/$src-$trg-$numMix.cbres");

	system("$todaBin/vqlbl/vqlbl -nmsg -dim 48 -cbookfile $workDir/cbook/$numMix-mix/$src-$trg.cb1.mat \\
	-weightfile $workDir/gmm_init/$numMix-mix/$src-$trg.tmp.wght -covfile $workDir/gmm_init/$numMix-mix/$src-$trg"."1.cov $jointFile.f $workDir/gmm_init/$numMix-mix/$src-$trg.tmp.lbl");
	system("rm -f $workDir/gmm_init/$numMix-mix/$src-$trg.tmp.wght $workDir/gmm_init/$numMix-mix/$src-$trg.tmp.lbl");
	system("$todaBin/vqlbl/vqlbl -nmsg -dim 48 -cbookfile $workDir/cbook/$numMix-mix/$src-$trg.cb$numMix.mat \\
	   	-weightfile $workDir/gmm_init/$numMix-mix/$src-$trg"."$numMix.wght -covfile $workDir/gmm_init/$numMix-mix/$src-$trg"."$numMix.cov $jointFile.f $workDir/gmm_init/$numMix-mix/$src-$trg.tmp.lbl");
	system("rm -f  $workDir/gmm_init/$numMix-mix/$src-$trg.tmp.lbl");

	system("$todaBin/gmm/gmm_jde -flcoef 0.001 -dim 48 -wghtfile $workDir/gmm_init/$numMix-mix/$src-$trg"."$numMix.wght \\
		-meanfile $workDir/cbook/$numMix-mix/$src-$trg".".cb$numMix.mat -covfile $workDir/gmm_init/$numMix-mix/$src-$trg"."$numMix.cov \\
		-bcovfile $workDir/gmm_init/$numMix-mix/$src-$trg"."1.cov -dia  $jointFile.f \\
		$workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.wght $workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.mean $workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.cov >$workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.res");
	system("$todaBin/gmmpara/gmm_para -nmsg -dim 48 -ydim 24 -dia $workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.wght \\
	   	$workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.mean $workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.cov $workDir/gmm_para/$numMix-mix/$src-$trg"."$numMix.param");
	#my $dtw_res="$workDir/$src-$trg"."_frmcd.res";
	#system("csh $todaScript/get_twf.csh $workDir/mcep $workDir/dtw $todaBin $srcList $trgList mcep $dtw_res");
	#system("$prjDir/vc/scripts/VCTrain train $src $trg"."_$numMix"."mix $numMix");

	system("$SPTK/x2x +da  $workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.mean > $workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.mean.txt ");
	system("$SPTK/x2x +da  $workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.wght > $workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.wght.txt ");
	system("$SPTK/x2x +da  $workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.cov > $workDir/gmm_jde/$numMix-mix/$src-$trg"."$numMix.cov.txt ");


	#ar-gmm training
}

if ($doTest)
{
=a
  for my $iFile (@testfiles)
    {

        #print "test file $iFile\n";
        system(
            "mv $workDir/wav/$src/$iFile.wav $testDir/wav/$src/$iFile.wav"
        );
    }
=cut
	for my $iFile (@testfiles)
    {
        system(
            "cp $prjDir/vc/train/mcepc0/$src/$iFile.mcep $testDir/wav/$src-$trg"."_$numMix"."mix/$iFile.wav.org.mcep"
        );
    }
   mkdir "$testDir/wav/convert",0755;
   mkdir "$testDir/wav/convert-$order",0755;
   open M,">$testDir/wav/convert-$order/convert.m" or die "can't oepn";
   print M "addpath('$prjDir/script');\n";
   print M "addpath('$prjDir/script/ar-gmm');\n";

   for my $iFile (@testfiles)
   {
   		my $srcMcep= "$testDir/wav/$src-$trg"."_$numMix"."mix/$iFile.wav.org.mcep";
#		open SRCMCEP,"<$srcMcep" or die "can't oepn";
#		my @Stat=stat(SRCMCEP);	
#		my $data;
#		read(SRCMCEP,$data,$Stat[7]);
#		close SRCMCEP;
#		my $n=$Stat[7]/8;
#		my $T=$n/25;
#		my @src=unpack("d$n",$data);
#		#print "@src";
		print M "fprintf(\'converting $srcMcep \\n \');\n";
		print M "fd=fopen(\'$testDir/wav/convert-$order/$iFile.mcep.txt',\'w\');\n";
#		my $num=($#src+1)/25;
#		#	print $num;
#		my $expr="x=[";
#		for my $i(1 .. $num)
#		{
#			
#			$expr.="@src[ ($i-1)*25..$i*25-1]";
#			$expr.="...\n";
#		}
#		$expr .="]';\n";
#		#print M "x=[@src[0..100]]';\n";
#		print M $expr;
 		print M "fp=fopen(\'$srcMcep\',\'rb\');\n";
		print M "x=fread(fp,inf,\'float32\');\n";
		print M "fclose(fp);\n";
		print M "T=size(x,1)/25;\n x=reshape(x,25,T);\n";

		my $model= "$workDir/gmm_para/$numMix-mix/$src-$trg"."$numMix"."_ordery"."$order"."_fullbt.para.txt";
		print M "predict=gmm_ar_ABdiag_convert(x,48,$order,$order,$numMix,\'$model\',\'spBfull\');\n";
		print M "fprintf(fd,\'%f\\n\',predict);\n";

#		for my $t(0..$T-1)
#		{
#			#	print M "fprintf(fd,\'$alltext[$t]\\n\');\n";
#			print M "fprintf(fd,\'%f\\n\',predict($t*24+1:($t+1)*24));\n";
#			#print M "fprintf(fd,\'#\\n\');\n";	
#		}
		print M "fclose(fd);\n";
		system("$SPTK/x2x +df $srcMcep >$srcMcep.f");
		system("$prjDir/script/ar-gmm/ar-gmmC/ar-gmm_convert -l 48 -m $numMix -o $order  -g $workDir/ar-gmm/$numMix-mix/gvModel -i 40 -s 1 -M $workDir/ar-gmm/$numMix-mix/$src-$trg-order$order.ar-gmm $srcMcep.f $testDir/wav/convert-$order/$iFile.mcep");
		print "$prjDir/script/ar-gmm/ar-gmmC/ar-gmm_convert -l 48 -m $numMix -o $order  -g $workDir/ar-gmm/$numMix-mix/gvModel -i 20 -s 0.001 -M $workDir/ar-gmm/$numMix-mix/$src-$trg-order$order.ar-gmm $srcMcep.f $testDir/wav/convert-$order/$iFile.mcep";
		system("$SPTK/x2x +fa $testDir/wav/convert-$order/$iFile.mcep > $testDir/wav/convert-$order/$iFile.mcep.txt");
   }
  
   #exit;
   #system("$MATLAB<$testDir/wav/convert/convert.m");
   #exit;

    for my $iFile (@testfiles)
    {
        system(
            "$SPTK/x2x +af $testDir/wav/convert-$order/$iFile.mcep.txt|$SPTK/mgc2sp -a 0.42 -m 24 -l 1024 -o 2 |$SPTK/x2x -o +fa513> $testDir/wav/convert-$order/$iFile.sp "
        );
    }



}
