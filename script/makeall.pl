#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: makeall.pl
#
#        USAGE: ./makeall.pl  
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
#      CREATED: 2013年01月16日 21时57分29秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;


my $doMdall=1;

my $prjDir="..";
if($doMdall)
{
	for my $emotion("sad","happy","angry")
	{
		my $testDir="$prjDir/test/tone_dall/$emotion";
		system("mv  $testDir/mrmse.phrase  $testDir/mrmse.phrase.back ;mv $testDir/mrmse.tone $testDir/mrmse.tone.back");
		system("./twotier_tonedynamicall.pl ../config/config.pm $emotion 1 1 init");
		for (my $mix1=1;$mix1<20;$mix1+=1)
		{
			{
				
				system("./twotier_tonedynamicall.pl ../config/config.pm $emotion $mix1 $mix1 train");
			}
		}
		
		#system("./twotier_tonedynamicall.pl ../config/config.pm $emotion 7 8 resynth");
	}
}
