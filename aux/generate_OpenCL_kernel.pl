#!/usr/bin/env perl
use warnings;
use strict;
use v5.20;
my $macros_kernel_src = './macros_kernel.h';

if (not -e $macros_kernel_src  && !@ARGV) {
    die "Please specify the source for the macros NTH and NUNITS\n";
}
system($ENV{HOME}.'/Git/RefactorF4Acc/bin/refactorF4acc.pl -P translate_to_OpenCL -c rf4a_to_C.cfg adam_bondv1_feedbf_les_press_v_etc_superkernel');
system('cp  module_adam_bondv1_feedbf_les_press_v_etc_superkernel.cl module_adam_bondv1_feedbf_les_press_v_etc_superkernel_ORIG.cl');
open my $MK, '<', $macros_kernel_src or die $!;
my @ls=<$MK>;
close $MK;
my $macros_str=join(" ",map {$_=~s/\n//;s/^\s*//;s/\s*$//;s/.define\s*/-D/;s/.undef\s*/-U/;s/\s+/=/;$_} @ls);
say("cpp $macros_str -I. -P module_adam_bondv1_feedbf_les_press_v_etc_superkernel.cl > module_adam_bondv1_feedbf_les_press_v_etc_superkernel_after_CPP_for_debugging.cl");
system("cpp $macros_str -I. -P module_adam_bondv1_feedbf_les_press_v_etc_superkernel.cl > module_adam_bondv1_feedbf_les_press_v_etc_superkernel_after_CPP_for_debugging.cl");
