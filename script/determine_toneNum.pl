#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: determine_toneNum.pl
#
#        USAGE: ./determine_toneNum.pl  mode
#
#  DESCRIPTION: 不同音节声调DCT个数时的重建误差
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Li Xian (), shysian@gmail.com
# ORGANIZATION: USTC
#      VERSION: 1.0
#      CREATED: 2013年01月23日 21时21分06秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my $emotion='neutral';
my $SPTK = "/usr/local/SPTK/bin";
my $prjDir="..";
my $workDir="$prjDir/test/error_dctnum/$emotion";

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

my ($calculate,$plot);
my $mode=$ARGV[0];
if($mode eq 'calculate')
{
	 $calculate=1;
	 $plot=0;
}
elsif($mode eq 'plot')
{
	$calculate=0;
	$plot=1;
}

if($calculate==1)
{
####f0 mean
my $f0mean;
my $f0statFile= "$prjDir/train/linear/happy/neutral.stat" ;
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

#$f0mean=0;


mkdir "$prjDir/test/error_dctnum",0755;
mkdir $workDir,0755;
mkdir "$workDir/f0",0755;
#my @testfiles= (201, 215,224, 244,260, 269, 275,283,290, 295,302,316, 321, 330,341,344,352,361, 376,381, 395,403,414, 421,432, 441,456, 461,478, 489);
my @testfiles=(201..500);


for my $iFile(@testfiles)
{
	my @voicedF0=extractVoicedF0("$prjDir/$emotion/logf0/$iFile.f0","$prjDir/$emotion/lab/mlevel/$iFile.lab");
	open VF0,">$workDir/f0/$iFile.f0" or die "can't oepn $!\n";
	print VF0 $_."\n" for(@voicedF0);
	system("cp $prjDir/$emotion/dct/tone_dall/$iFile.phrase_f0 $workDir/f0/$iFile.phrase_f0");
	system("cp $prjDir/$emotion/dct/tone_dall/$iFile.f0 $workDir/f0/$iFile.tone_f0_voiced");
	interpolate("$workDir/f0/$iFile.tone_f0_voiced","$workDir/f0/$iFile.tone_f0");
}
for my $dct_num(4..20)
{
		mkdir "$workDir/$dct_num",0755;
		my @files = <$workDir/f0/*.tone_f0_voiced>;
        foreach (@files)
        {

            my $name = `basename $_ .tone_f0_voiced`;
            chomp($name);
            print "extracting tone dct :$name\n";
            open F0,   "<$workDir/f0/$name.tone_f0_voiced"   or die "can't open";
            open TONE, ">$workDir/$dct_num/$name.tone" or die "can't open";
            my $line;
            while ($line = <F0>)
            {

                if ($line =~ /[a-z]+/)
                {

                    open TEMP, ">$workDir/$dct_num/temp" or die "can't open";
                    print TONE $line;
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
                      `x2x +af $workDir/$dct_num/temp|dct -l $length|sopr -d sqrt$length|x2x +fa`;

                    #print "current length is ".$#dct."+1\n";
                    while ($#dct < $dct_num)
                    {
                        push @dct, "0\n";
                    }
                    @dct = splice @dct, 0, $dct_num;

                    #print "@dct\n";

                    print TONE @dct;
                    print TONE "#\n";
                    system("rm -f $workDir/$dct_num/temp");
                }    #end if
            }    #end while
            close(TONE);
            close(F0);
        }    #end foreach
	my @ToneFiles = <$workDir/$dct_num/*.tone>;
    for my $ToneFile (@ToneFiles)
    {
        my $fname = `basename $ToneFile .tone`;
        $fname =~ s/\n//;
        open TONE, "<$ToneFile"  or die "can't open $!\n";
		
        open REF0, ">$workDir/$dct_num/$fname.tone_f0_voiced" or die "can't open $!\n";
        local $/ = "#\n";
        while (<TONE>)
        {
            chomp;
            my @toneParas = split /\n/, $_;
            my $text = shift @toneParas;
            $text =~ /([0-9]+)\s([0-9]+)\s/;
            my $lenOfSeg = ($2 - $1) / 50000;
            my $bcutEnd  = $lenOfSeg - 1;
            my @f0OfSeg =
              `echo @toneParas|$SPTK/x2x +af|$SPTK/bcut +f -s 0 -e $bcutEnd|$SPTK/idct -l $lenOfSeg|$SPTK/sopr -m sqrt$lenOfSeg|$SPTK/x2x +fa`;
            print REF0 "$text\n";
            print REF0 @f0OfSeg;
            print REF0 ".\n";
        }
        close TONE;
        close REF0;

		interpolate("$workDir/$dct_num/$fname.tone_f0_voiced","$workDir/$dct_num/$fname.tone_f0");
		#interpolate("$workDir/$dct_num/tonetemp", "$workDir/$dct_num/$fname.tone.f0");
    }
	#calculate RMSE
	open OTEMP,">$workDir/$dct_num/originf0" or die "can't open $!\n";
	open RTEMP,">$workDir/$dct_num/reconstructf0" or die "can't open $!\n";
	for my $iFile(@testfiles)
	{
		open OF0,"<$workDir/f0/$iFile.f0" or die "can't open $!\n";
		open REF0,"<$workDir/$dct_num/$iFile.tone_f0_voiced" or die "can't open $!\n";
		
		open REF0_all,">$workDir/$dct_num/$iFile.ref0" or die "can't open $!\n";
		
		my $count=0;
		while(<OF0>)
		{
			chomp;
		 unless(/[a-z]+$/ or  /^\./)
		 {
			my $temp=$_;
		 	print OTEMP "$temp\n";
			$count++;
		 }
		}

		my @voicedF0=extractVoicedF0("$workDir/f0/$iFile.phrase_f0","$prjDir/$emotion/lab/mlevel/$iFile.lab");
		my $i=0;
		while(<REF0>)
		{
			chomp;
			 unless(/[a-z]+$/ or /^\./)
			 {
				my $temp=$_+$f0mean+$voicedF0[$i];
				print RTEMP "$temp\n";
				print REF0_all "$temp\n";
				$count--;
				$i++;
			 }
		}
		close OF0;
		close REF0;
		close REF0_all;
		if($count !=0 ) 
		{
			print "$iFile error\n";
		}

	}
	close OTEMP;
	close RTEMP;
	system("$SPTK/x2x +af $workDir/$dct_num/originf0  >$workDir/$dct_num/originf0.f");
	system("$SPTK/x2x +af $workDir/$dct_num/reconstructf0 >$workDir/$dct_num/reconstructf0.f");
	system("$SPTK/sopr -EXP $workDir/$dct_num/originf0.f >$workDir/$dct_num/originf0.hz.f");
	system("$SPTK/sopr -EXP $workDir/$dct_num/reconstructf0.f >$workDir/$dct_num/reconstructf0.hz.f");
	system("$SPTK/rmse -l 0 $workDir/$dct_num/originf0.hz.f $workDir/$dct_num/reconstructf0.hz.f|$SPTK/x2x +fa >$workDir/$dct_num/reconstructRmse");

}
}
if($plot==1)
{
	my $workDir="$prjDir/test/error_dctnum";
	for my $dir("neutral","sad","angry","happy")
	{

		#	open ALL,"$workDir/$dir/reconstructRmse" or die "can't open $!\n";
		system("rm -f $workDir/$dir/reconstructRmse");
		for my $dct_num(4..20)
		{
			system("cat $workDir/$dir/$dct_num/reconstructRmse >> $workDir/$dir/reconstructRmse");
		}
		#	close ALL;
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
	local $/=".\n";
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
sub extractVoicedF0
{
	my ($f0file,$labfile)=@_;
	open F0,"<$f0file" or die "can't open $!\n";
	open LAB,"<$labfile" or die "can't open $!\n";
  
	my @f0data=<F0>;
	chomp @f0data;
	
	my @voicedF0;
	my $line2=<LAB>;#skip 1 line;

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
				my $start = $text[0] / 50000;
				my $end   = $text[1] / 50000 ;
				my $i     = $start;
				while ($i < $end)
				{
					push @voicedF0,$f0data[$i];
					$i++;
				}

			}
			elsif ($initialType == 2)
			{
				my $start = $text[0] / 50000;
				$line2 = <LAB>;
				@text = split /\s+/, $line2;
				my $end = $text[1] / 50000 ;
				my $i   = $start;
				while ($i < $end)
				{
					push @voicedF0,$f0data[$i];
					$i++;
				}


			}
			elsif ($initialType == 0)
			{
				$line2 = <LAB>;
				@text = split /\s+/, $line2;
				my $start = $text[0] / 50000;
				my $end   = $text[1] / 50000 ;
				my $i     = $start;
				while ($i < $end)
				{
				
					push @voicedF0,$f0data[$i];
					$i++;
				}
			}
			else
			{

				print "may have ERROR $word\n"
				  unless ($text[2] eq 'sil' or $word eq 'sp');
			}
		}
		return @voicedF0;

}
