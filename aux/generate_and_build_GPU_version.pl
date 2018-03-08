#!/usr/bin/env perl
use v5.10;
use strict;
use warnings;

BEGIN {
	push @INC, './aux';
};

use Cwd;
use RunCpp qw( run_cpp );
use RestoreStashedLines qw( restore_stashed_lines );
use MacroFileToCmdLine qw( macro_file_to_cmd_line_str );

my $wd = cwd();

my $VV=1;

my @kernel_sources=qw(./adam.f95 ./bondv1.f95 ./feedbf.f95 ./les.f95 ./press.f95 ./velFG.f95 ./velnw.f95);
my $kernel_sources_str = join(' ',@kernel_sources);

my $main_src = './main.f95';
my $iowrite_subs = 'anime';

my $plat = 'GPU';

my $TRUST_THE_COMPILER = 1;
my $skip_step_1 = 0;
my $skip_step_2 = 0;
my $skip_step_3 = 0;
my $skip_step_4 = 1;
if (@ARGV) {
	$TRUST_THE_COMPILER = $ARGV[0];
}

my $gen_dir = 'GeneratedCode';
if ($TRUST_THE_COMPILER==1) {
	$gen_dir = 'GeneratedCodeV2';
}

# The compiler fails if this directory does not exists
if (not -d $gen_dir) {
	mkdir $gen_dir;
}
 
if (not $skip_step_1) {
if ($TRUST_THE_COMPILER) {
	chdir 'src';
	##
	say'*NOTE 2018-03-07* 
	The `AutoParallel-Fortran` compiler has built-in handling of macros via the -D and -X flags. 
	This should generate the same code as when using the `run_cpp.pl` and `restore_stashed_lines.pl` scripts. 
	However, this is _TO BE TESTED_!
	' if 0;
	
	(my $defined_macros_str, my $undef_macros_str) = macro_file_to_cmd_line_str( './macros.h','-D');
	(my $macros_to_skip_str, my $empty_str) = macro_file_to_cmd_line_str('./macros_to_skip.h','-X');
	
	say("AutoParallel-Fortran-exe $kernel_sources_str -out ../$gen_dir/ -iowrite $iowrite_subs -main $main_src -plat $plat  $defined_macros_str $macros_to_skip_str -v");
	#system('which AutoParallel-Fortran-exe');die;
	system("AutoParallel-Fortran-exe $kernel_sources_str -out ../$gen_dir/ -iowrite $iowrite_subs -main $main_src  -plat $plat  $defined_macros_str $macros_to_skip_str -v");	
	
} else {	

	say '* First, in `src`, run CPP on the code using the macros in `macros.h` and stash lines guarded with macros from `macros_to_skip.h`. This generates the file `stash.pl`' if $VV;
	
	chdir 'src';
	
	run_cpp();
	  
	##
	say '* Then, in `PostCPP`, run the OpenCL compiler `AutoParallel-Fortran-exe`. This will take a while and produce a lot of output, which you can ignore.' if $VV;
	
	chdir $wd;
	if (not -d 'PostCPP') {
		mkdir 'PostCPP';
	}

	 
	chdir 'PostCPP';
#	system('which AutoParallel-Fortran-exe');
	say("AutoParallel-Fortran-exe $kernel_sources_str -out ../$gen_dir/ -iowrite $iowrite_subs -main $main_src -v -plat $plat");
#	die cwd();
	system("AutoParallel-Fortran-exe $kernel_sources_str -out ../$gen_dir/ -iowrite $iowrite_subs -main $main_src -v -plat $plat");
	
	##
	say "* In '$gen_dir', we restore code segments that were stashed in the previous step" if $VV;
	chdir $wd;
	chdir $gen_dir;
	
	restore_stashed_lines("$wd/src/stash.pl"); 
	system('cp ./PostGen/* .');

}
}

if (not $skip_step_2) { 
	
##
say "* In `$gen_dir`, we copy the non-modified source files into the current folder, as well as some scripts and config files needed to build the OpenCL kernel." if $VV;
chdir $wd;
chdir $gen_dir;

if (not -d 'data') {
	system('cp -r ../data .');
}
if (not -d 'GIS') {
	system('cp -r ../GIS .');
}

my @sources=qw(anime.f95
aveflow.f95
bondFG.f95
bondv1.f95
boundp.f95
boundsm.f95
common_sn.f95
feedbfm.f95
grid.f95
ifdata.f95
init.f95
params_common_sn.f95
set.f95
timdata.f95
timseris.f95
macros.h
macros_to_skip.h
);

my $ref_dir = $TRUST_THE_COMPILER ? "$wd/src" : "$wd/PostCPP";
for my $src (@sources) {
   system("cp $ref_dir/$src ."); 
}

my @sources2=qw(macros_kernel.h
SConstruct.auto
generate_OpenCL_kernel.sh
rf4a_to_C.cfg
array_index_f2c1d.h);

my $ref_dir_2 = "$wd/aux";
for my $src (@sources2) {
    system("cp $ref_dir_2/$src . ");
}

  


}


if (not $skip_step_3) {
##
say '* Then we generate the actual OpenCL kernel code using `RefactorF4Acc`' if $VV;
chdir $wd;
chdir $gen_dir;

my $macros_kernel_src = './macros_kernel.h';

if (not -e $macros_kernel_src ) {
    die "Please specify the source for the macros NTH and NUNITS\n";
}
system($ENV{HOME}.'/Git/RefactorF4Acc/bin/refactorF4acc.pl -P translate_to_OpenCL -c rf4a_to_C.cfg adam_bondv1_feedbf_les_press_v_etc_superkernel');
system('cp  module_adam_bondv1_feedbf_les_press_v_etc_superkernel.cl module_adam_bondv1_feedbf_les_press_v_etc_superkernel_ORIG.cl');

# Unused, for debugging
#open my $MK, '<', $macros_kernel_src or die $!;
#my @ls=<$MK>;
#close $MK;
#my $macros_str=join(" ",map {
#	$_=~s/\n//;
#	s/^\s*//;
#	s/\s*$//;
#	s/.define\s*/-D/;
#	s/.undef\s*/-U/;
#	s/\s+/=/;
#	$_
#} @ls);
#say("cpp $macros_str -I. -P module_adam_bondv1_feedbf_les_press_v_etc_superkernel.cl > module_adam_bondv1_feedbf_les_press_v_etc_superkernel_after_CPP_for_debugging.cl");
#system("cpp $macros_str -I. -P module_adam_bondv1_feedbf_les_press_v_etc_superkernel.cl > module_adam_bondv1_feedbf_les_press_v_etc_superkernel_after_CPP_for_debugging.cl");

}
if (not $skip_step_4) {
##
chdir $wd;
chdir $gen_dir;
say  "Now we can build the OpenCL Fortran host code, setting the number of threads and compute units depending on the GPU";
say 'Note that the Scons build runs cpp on the kernel for the macros NTH, NUNITS and BARRIER_OK';
say 'Normally these are set via the `nth`, `nunits` and `dev` flags on the scons command line';
say 'But you can also put them in `macros_kernel.h`';

say("scons -f SConstruct.auto -s mcm=m dev=$plat nth=256 nunits=16");
system("scons -f SConstruct.auto -s mcm=m dev=$plat nth=256 nunits=16");
}