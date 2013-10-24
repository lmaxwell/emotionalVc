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
my $gv=0;
=a
for my $trg ('sad','happy','angry')
{
	system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 6 $gv init");
	system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 6 $gv prepare");
}
=cut
=a
for my $order (3..5)
{
	for my $trg ('sad')
	{
		print "train $trg-order$order\n";
		#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 $order $gv train");
		print "test $trg-order$order\n";
		system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 $order $gv test");
		system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 $order 1 test");
		#print "convert $trg-order$order\n";
		#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 $order $gv resynth");
		#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 $order 1 resynth");
	}
}
for my $order (0..5)
{
	for my $trg ('happy','angry')
	{
		print "train $trg-order$order\n";
		#system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 $order $gv train");
		print "test $trg-order$order\n";
		system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 $order $gv test");
		system("$perl ar-gmm/sp_ar_gmm.pl ../config/config.pm $trg 64 $order 1 test");
		#print "convert $trg-order$order\n";
		#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 $order $gv resynth");
		#system("$perl ar-gmm/twotier_nodynamic_ar_ABfull.pl ../config.pm $trg 6 12 3 64 $order 1 resynth");
	}
}
=cut
for my $trg('sad','happy','angry')
{
	open MELCD,">../vc/train/$trg.melcd" or die "can't open\n";
	for my $order(0..5)
	{
		my $temp=`cat "../vc/train/to$trg""_ar/test/64-mix/mel-CD/$order-g0/melCD_test"`;
		print MELCD $temp;
	}
	close MELCD;
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


