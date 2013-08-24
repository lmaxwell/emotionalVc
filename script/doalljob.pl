#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: doalljob.pl
#
#        USAGE: ./doalljob.pl  
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
#      CREATED: 2013年08月21日 21时12分40秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my $perl="/usr/bin/perl";


for my $trg ('happy','angry')
{
	system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 6 prepare");
}
for my $order (0..5)
{
	for my $trg ('sad','happy','angry')
	{
		print "train $trg-order$order\n";
		system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 $order train");
		print "test $trg-order$order\n";
		system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 $order test");
		print "convert $trg-order$order\n";
		system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 $order resynth");
	}
}
#	system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 4 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 4 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 4 resynth");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 2 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 2 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 2 resynth");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 1 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 1 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 1 resynth");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 0 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 0 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 0 resynth");
exit;
#$trg='happy';
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 5 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 5 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 5 resynth");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 4 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 4 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 4 resynth");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 2 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 2 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 2 resynth");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 1 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 1 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 1 resynth");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 0 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 0 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 0 resynth");
#
#
#$trg='angry';
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 5 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 5 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 5 resynth");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 4 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 4 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 4 resynth");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 2 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 2 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 2 resynth");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 1 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 1 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 1 resynth");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 0 train");
#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 0 test");
#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 0 resynth");


