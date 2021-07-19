#!/usr/bin/perl 
use v5.30;
use strict;
use warnings;

# - Find the file with the main program. It is the file starting with 'gen_'
    my @main_files =  glob('gen_*.f95');
    if (scalar @main_files>1) {
        die 'Too many main files!';
    }
    my $main_file = shift @main_files;
    my $superkernel_file = $main_file;$superkernel_file=~s/^gen_//;

# - From the superkernel file, get the module name, subroutine name and stage kernel name(s)
    my $module_name = `grep -E '^module' $superkernel_file`;
    $module_name=~s/module\s+//;
    my $sub_name = `grep -E '^subroutine' $superkernel_file`;
    $sub_name=~s/subroutine\s+//;
    $sub_name=~s/\s*\(.+$//;
    my @stage_kernel_names = `grep -E '^call .+_kernel' $superkernel_file`;
    for my $stage_kernel_name (@stage_kernel_names) {
        $stage_kernel_name=~s/call\s+//;
        $stage_kernel_name=~s/\s*\(.+$//;
    }

# - In the superkernel file, in the superkernel subroutine, add the use declarations for the stage kernels
    open my $SKMF, '<', $superkernel_file or die $!;
    my @superkernel_file_lines = <$SKMF>;
    close $SKMF;
    open $SKMF, '>', $superkernel_file or die "$!";
    for my $line (@superkernel_file_lines) {
        say $SKMF $line;
        if ($line=~/^subroutine $sub_name/) {
            for my $stage_kernel_name (@stage_kernel_names) {
                say $SKMF "use singleton_module_${stage_kernel_name}, only: ${stage_kernel_name}";
            }
        }
    }
    close $SKMF;
    
# - In the main file:
    # - get the module line and correct the name
    open my $MF, '<', $main_file or die $!;
    my @main_file_lines = <$MF>;
    close $MF;
    open $MF, '>', $main_file or die $!;
    for my $line (@main_file_lines) {
        say $MF $line;
        if ($line=~/^use.+only:\s+$sub_name/) {
            say $MF "use $module_name, only : $sub_name";
        }
    }
    close $MF;
# - In the stage kernel file
for my $stage_kernel_name (@stage_kernel_names) {
    my $stage_kernel_file=$stage_kernel_name.'.f95';
    open my $SKF, '<', $stage_kernel_file or die $!;
    my @stage_kernel_file_lines = <$SKF>;
    close $SKF;
    my %scalar_functions=();
    # - find all unique calls     
    for my $line (@stage_kernel_file_lines) {
        if ($line=~/call\s+(\w+?_scal)/) {
            $scalar_functions{$1}=1;
        }
    }
    open $SKF, '>', $stage_kernel_file or die $!;
    for my $line (@stage_kernel_file_lines) {
        # Rather ad-hoc
        # replace `call get_global_id(idx,0,global_id)` by `idx=global_id`
        if ($line=~/call\s+get_global_id/) {
            say $SKF "idx = global_id"; 
        } else {
            say $SKF $line;
        }
    #     if ($line=~/subroutine\s+$stage_kernel_name/) {
    # # - add the use declarations for the called functions:    
    #         for my $f (sort keys %scalar_functions) {
    #             say $SKF "use singleton_module_${f}, only: ${f}_scal";
    #         }
    #     }
    }
    close $SKF;
}