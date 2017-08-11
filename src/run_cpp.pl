#!/usr/bin/env perl
use v5.22;
use warnings;
use strict;

use Cwd;
my $wd=cwd();

my $VV=1;

my @macros_to_keep=qw( I_ANIME )
# To support this we need extra work; to do it right quite a bit of extra work.
# I can't find anything decent on CPAN either so we'll have to do it from scratch
# Basically I need to keep a stack for the #if/#endif blocks
my @blocks_stack=();
my @stored_lines=();
# Then go through all lines, if we find a line with a matching macro, no matter if it is #if, #ifdef or #ifndef, we push all lines in the block onto @stored_lines
# We replace them by a "7188 continue" as in the other code


my @macros=qw( NO_GLOBAL_SOR TWINNED_BUFFER INLINE_BOUND_CALCS WV_TEST NO_IO NO_FILE_IO ICAL=0 IFBF=1  IADAM=0 WV_NEW );

my @srcdirs=qw(.);

my @includes = map {"-I$wd/".$_}  @srcdirs;
my $includestr = join(' ',@includes);

my @includes_l1 = map {"-I$wd/../".$_}  @srcdirs;
my $includestr_l1 = join(' ',@includes);

my @defines = map {"-D".$_}  @macros;
my $definestr = join(' ',@defines);

if (not -d "$wd/../PostCPP") {
    mkdir "$wd/../PostCPP";
}

my @srcs = glob("*.f95");
for my $src (@srcs) {
    my $srcf =$src;
#    $srcf=~s/F$/f/;
    my $cmd =  "cpp -Wno-invalid-pp-token -P  $includestr $definestr -Wno-extra-tokens $src | grep -v -E '^\\s*\\!|^\\s*\$' | perl -p -e 's/\\s*\\!.+\$//;s/\£/\#/' > $wd/../PostCPP/$srcf";
    say $cmd if $VV;
    system($cmd);
}

for my $srcdir (@srcdirs) {
    if (not -d "$wd/../PostCPP/$srcdir") {
        system("mkdir -p $wd/../PostCPP/$srcdir");
    }
    chdir "$wd/$srcdir";
    my @srcs = glob("*.f95");
    for my $src (@srcs) {
        my $srcf =$src;
        my $cmd = "cpp -Wno-invalid-pp-token -P  $includestr_l1 $definestr -Wno-extra-tokens $src  | grep -v -E '^\\s*\\!|^\\s*\$' | perl -p -e 's/\\s*\\!.+\$//;s/\£/\#/' > $wd/../PostCPP/$srcdir/$srcf";
        say $cmd if $VV;
        system( $cmd );  
    }
    chdir $wd;
}
