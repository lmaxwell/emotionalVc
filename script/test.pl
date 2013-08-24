#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: test.pl
#
#        USAGE: ./test.pl  
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
#      CREATED: 2013年01月07日 21时01分29秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
my %mixNum=(
	"sad" => [5,8],
	"angry" =>3,

);
print $mixNum{"sad"}[0];

my @array=('test');
my $ref=\@array;
print  " $ref->[0]";
my $ref2=['test'];
print " $$ref2[0]";
my $ref3=\('test');
print " $$ref3\n";

my @array1=('test2');

local $/="#\n";
#my @data=("1\n","2\n","3\n","4\n");
my @data=`echo 1 2 3 4|x2x +af|bcut +f -s 0 -e 3 |idct -l 50|sopr -m sqrt4|x2x +fa`;
print "length is $#data\n";
my @temp=map {s/\n//;$_} @data;
#print "@temp\n";
my @fdata=map {$_-10} @temp;
#print "@fdata\n";
chomp @data;
#print "@data\n";
