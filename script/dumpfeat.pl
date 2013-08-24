#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: dumpfeat.pl
#
#        USAGE: ./dumpfeat.pl  
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
#      CREATED: 2013年01月09日 14时59分59秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;




my $prjDir="..";


my %dict;
our %initials;
our %finals;

open DICT,"<$prjDir/general/dict" or die "can't open $!\n";
while(<DICT>)
{
	my ($syllable,$inital,$final)=split /\s+/;
	$dict{$syllable}=[$inital,$final];
}
while( my($key,$value)=each %dict)
{
	#print "$key @$value\n";
	$initials{$value->[0]}=0 unless ($value->[0] eq 'NULL');
	$finals{$value->[1]}=0 ;
}
print "initials:\n";
my $numInitial=printHash(\%initials);
print "num of Initial :$numInitial\n";
print "finals\n";
my $numFinal=printHash(\%finals);
print "num of Final :$numFinal\n";


my @monophthong=("a","o","e","i","u","v","ii","iii");
my @diphthong=("ai","ei","ao","ou","ie","ve","an","en","in","ia","ua","uo","er","vn");
my @triphthong=("ang","eng","ong","ing","iao","iou","ian","iang","iong","uang","uai","van","uan","uen","uei");
my @allphthong=(@monophthong,@diphthong,@triphthong);
print "total $#allphthong+1\n";
for(@allphthong)
{
	print "$_ is wrong\n" unless(exists$finals{$_})
}

our %phoneType;
for(@monophthong)
{
	$phoneType{$_}="m";
}

for(@diphthong)
{
	$phoneType{$_}="d";
}
for(@triphthong)
{
	$phoneType{$_}="t";
}
my @stop=("b","p","d","g","k","t");
my @affricate=("z","zh","j","c","ch","q");
my @fricative=("f","s","sh","x","h","r");
my @nasal=("m","n");
my @lateral=("l");
our @allconsonant=(@stop,@affricate,@fricative,@nasal,@lateral);
print "total $#allconsonant+1\n";
for (@allconsonant)
{
	print "$_ is wrong\n" unless(exists$initials{$_});
}


for(@stop)
{
	$phoneType{$_}="s";
}

for(@affricate)
{
	$phoneType{$_}="a";
}
for(@fricative)
{
	$phoneType{$_}="f";
}

for(@nasal)
{
	$phoneType{$_}="n";
}
for(@lateral)
{
	$phoneType{$_}="l";
}
$phoneType{'sil'}='sil';
$phoneType{'sp'}='sp';

printHash(\%phoneType);
print "total ". keys %phoneType;
##################sub #########################
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

