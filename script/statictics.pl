#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: statictics.pl
#
#        USAGE: ./statictics.pl  dumpfeat.pl
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
#      CREATED: 2013年01月29日 14时18分33秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
our %initials;
our %finals;
our %phoneType;
require("dumpfeat.pl");

my @testfiles = (201, 215,224, 244,260, 269, 275,283,290, 295,302,316, 321, 330,341,344,352,361, 376,381, 395,403,414, 421,432, 441,456, 461,478, 489);

my $SPTK="/usr/local/SPTK/bin";

my %testfiles;
for my $iFile(@testfiles)
{
	$testfiles{$iFile}=1;
}
my %trainfiles;
my @trainfiles;
for my $iFile(201..500)
{
	unless(exists $testfiles{$iFile})
	{
		$trainfiles{$iFile}=1;
		push @trainfiles,$iFile;
	}
}
print "testfiles are :\n";
my $numtest=printHash(\%testfiles);
print "totally:$numtest\ntrainfiles are:\n";
my $numtrain=printHash(\%trainfiles);
print "totally:$numtrain\n";

print "phone num in testfiles:";
my $num_phone;
for my $iFile(@testfiles)
{
	my $phonelabfile="../neutral/lab/phone/$iFile.lab";
	open LAB,"<$phonelabfile" or die "can't open $!\n";
	while(<LAB>)
	{
		unless(/sp/ or /sil/)
		{
			$num_phone++;
		}
	}
}
print "$num_phone\n";

print "word num in testfiles:";
my $num_word;
for my $iFile(@testfiles)
{
	my $wordlabfile="../neutral/lab/word/$iFile.lab";
	open LAB,"<$wordlabfile" or die "can't open $!\n";
	while(<LAB>)
	{
		unless(/sp/ or /sil/)
		{
			$num_word++;
		}
	}
}
print "$num_word\n";

print "phrase num in testfiles:";
my $num_phrase;
open PB,"<../neutral/forcedAlignment/newphrasebreak" or die "can't open $!\n";
while(<PB>)
{
	/([0-9]+)\.lab/;
	my $iFile=$1;
	my $line=<PB>;
	if (exists$testfiles{$iFile})
	{	
		 $num_phrase+=() =$line=~m/(\/p|\/w)/g;
	}
}
print "$num_phrase\n";


print "phone num in trainfiles:";
my $num_phone;
for my $iFile(@trainfiles)
{
	my $phonelabfile="../neutral/lab/phone/$iFile.lab";
	open LAB,"<$phonelabfile" or die "can't open $!\n";
	while(<LAB>)
	{
		unless(/sp/ or /sil/)
		{
			$num_phone++;
		}
	}
}
print "$num_phone\n";

print "word num in trainfiles:";
my $num_word;
for my $iFile(@trainfiles)
{
	my $wordlabfile="../neutral/lab/word/$iFile.lab";
	open LAB,"<$wordlabfile" or die "can't open $!\n";
	while(<LAB>)
	{
		unless(/sp/ or /sil/)
		{
			$num_word++;
		}
	}
}
print "$num_word\n";

print "phrase num in trainfiles:";
my $num_phrase;
open PB,"<../neutral/forcedAlignment/newphrasebreak" or die "can't open $!\n";
while(<PB>)
{
	/([0-9]+)\.lab/;
	my $iFile=$1;
	my $line=<PB>;
	if (exists$trainfiles{$iFile})
	{	
		 $num_phrase+=() =$line=~m/(\/p|\/w)/g;
	}
}
print "$num_phrase\n";

mkdir "../statistics",0755;
for my $emotionType("neutral","happy","angry","sad")
{
	my $dir="../statistics/$emotionType";
	mkdir $dir,0755;
	open INITIAL,">$dir/initial.dur" or die "can't open $!\n";
	open FINAL,">$dir/final.dur" or die "can't open $!\n";
	open VOWEL,">$dir/vowel.dur" or die "can't open $!\n";
	open DI,">$dir/diph.dur" or die "can't open $!\n";
	open STOP,">$dir/stop.dur" or die "can't open $!\n";
	open FRIC,">$dir/fric.dur" or die "can't open $!\n";
	open AFRIC,">$dir/afric.dur" or die "can't open $!\n";
	open NAS,">$dir/nasal.dur" or die "can't open $!\n";
	open LAT,">$dir/lateral.dur" or die "can't open $!\n";
	print $emotionType.":\n";
	for my $iFile(201..500)
	{
		open LAB,"<../$emotionType/lab/phone/$iFile.lab" or die "can't open $!\n";
		while (<LAB>)
		{
			chomp;
			my ($start,$end,$phone)=split ' ',$_;
			my $len=($end-$start)/10000;
			if(exists $initials{$phone})
			{
				print INITIAL "$len\n";
			}
			elsif(exists $finals{$phone})
			{
				print FINAL "$len\n";
			}

			if($phoneType{$phone} eq 'm')
			{
				print VOWEL "$len\n";
			}
			elsif($phoneType{$phone} eq 'd' or $phoneType{$phone} eq 't')
			{
				print DI "$len\n";
			}
			elsif($phoneType{$phone} eq 's')
			{
				print STOP "$len\n";
			}
			elsif($phoneType{$phone} eq 'f')
			{
				print FRIC "$len\n";
			}
			elsif($phoneType{$phone} eq 'a')
			{
				print AFRIC "$len\n";
			}
			elsif($phoneType{$phone} eq 'n')
			{
				print NAS "$len\n";
			}

			elsif($phoneType{$phone} eq 'l')
			{
				print LAT "$len\n";
			}
		}
	}
	system("$SPTK/x2x +af $dir/initial.dur >$dir/initial.dur.f");
	system("$SPTK/vstat -l 1 -o 1 $dir/initial.dur.f |$SPTK/x2x +fa>$dir/initialmean");
	system("$SPTK/x2x +af $dir/final.dur >$dir/final.dur.f");
	system("$SPTK/vstat -l 1 -o 1 $dir/final.dur.f |$SPTK/x2x +fa>$dir/finalmean");


	system("$SPTK/x2x +af $dir/vowel.dur >$dir/vowel.dur.f");
	system("$SPTK/vstat -l 1 -o 1 $dir/vowel.dur.f |$SPTK/x2x +fa>$dir/vowelmean");
	system("cat $dir/vowelmean");

	system("$SPTK/x2x +af $dir/diph.dur >$dir/diph.dur.f");
	system("$SPTK/vstat -l 1 -o 1 $dir/diph.dur.f |$SPTK/x2x +fa>$dir/diphmean");
	system("cat $dir/diphmean");

	system("$SPTK/x2x +af $dir/stop.dur >$dir/stop.dur.f");
	system("$SPTK/vstat -l 1 -o 1 $dir/stop.dur.f |$SPTK/x2x +fa>$dir/stopmean");
	system("cat $dir/stopmean");

	system("$SPTK/x2x +af $dir/fric.dur >$dir/fric.dur.f");
	system("$SPTK/vstat -l 1 -o 1 $dir/fric.dur.f |$SPTK/x2x +fa>$dir/fricmean");
	system("cat $dir/fricmean");

	system("$SPTK/x2x +af $dir/afric.dur >$dir/afric.dur.f");
	system("$SPTK/vstat -l 1 -o 1 $dir/afric.dur.f |$SPTK/x2x +fa>$dir/africmean");
	system("cat $dir/africmean");

	system("$SPTK/x2x +af $dir/nasal.dur >$dir/nasal.dur.f");
	system("$SPTK/vstat -l 1 -o 1 $dir/nasal.dur.f |$SPTK/x2x +fa>$dir/nasalmean");
	system("cat $dir/nasalmean");

	system("$SPTK/x2x +af $dir/lateral.dur >$dir/lateral.dur.f");
	system("$SPTK/vstat -l 1 -o 1 $dir/lateral.dur.f |$SPTK/x2x +fa>$dir/lateralmean");
	system("cat $dir/lateralmean");
}


sub printHash
{
	my $hash=shift;
	my $num;
	while(my($key,$value)=each %$hash)
	{
		print "$key=>$value ";
		$num++;
	}
	print "\n";
	return $num;
}

