#!/usr/bin/env perl
use v5.22;
use warnings;
use strict;
use Data::Dumper;
use Cwd;
my $wd=cwd();

my $VV=0;

my $single=0;
my $single_src='';
if (@ARGV) {
    $single=1;
    $single_src=$ARGV[0];
}

my @macros_to_skip=();
my @defined_macros=();
my @undef_macros=();
my $no_macros=0;
if (-e 'macros.h') {
    open my $MACROS, '<', 'macros.h' or die $!;
    while(my $line=<$MACROS>) {
        chomp $line;
        $line=~s/^\s*#(\w+)\s+//;
        my $cmd = $1;
        $line=~s/\s+//;
        if ($cmd eq 'define') {
            push @defined_macros, $line;
        } elsif ($cmd eq 'undef') {
            push @undef_macros, $line;
        }
    }
    close $MACROS;
} else {
    $no_macros=1;
}

# NOTE: macros with #undef are ignored
if (-e 'macros_to_skip.h') {
    open my $MACROS, '<', 'macros_to_skip.h' or die $!;
    while(my $line=<$MACROS>) {
        chomp $line;
        $line=~s/^\s*#(\w+)\s+//;
        my $cmd = $1;        
        $line=~s/\s+//;
        if ($cmd eq 'define') {
            push @macros_to_skip, $line;
        }
    }
    close $MACROS;
} elsif ($no_macros) {
    die 'Could not find macros.h or macros_to_skip.h';
}




#die Dumper(@macros,@macros_to_skip);

if (not -d '../PostCPP') {
    mkdir '../PostCPP';
}
if (not -d '../PostCPP/PrePostCPP') {
    mkdir '../PostCPP/PrePostCPP';
}
my $stash_ref = {};


my @srcdirs=qw(.);

my @includes = map {"-I$wd/".$_}  @srcdirs;
my $includestr = join(' ',@includes);

my @includes_l1 = map {"-I$wd/../".$_}  @srcdirs;
my $includestr_l1 = join(' ',@includes);

my @defines = map {"-D".$_}  @defined_macros;
my $definestr = join(' ',@defines);
my @undefs = map {"-U".$_} @undef_macros;
my $undefstr = join(' ',@undefs);

if (not -d "$wd/../PostCPP") {
    mkdir "$wd/../PostCPP";
}

#my @srcs = $single ? ( $single_src ) : glob("*.f95");

#for my $src (@srcs) {
##    my $srcf =$src;
#    my $src_path = "$wd/../PostCPP/PrePostCPP/$src";
#    $stash_ref = stash_lines_to_skip($src,$stash_ref,\@macros_to_skip,$src_path);
#    #my $cmd_cpp =  $no_macros ? "cat ../PostCPP/PrePostCPP/$src " : "cpp -Wno-invalid-pp-token -P  $includestr $definestr -Wno-extra-tokens ../PostCPP/PrePostCPP/$src ";
#    #my $cmd_clean_up = "| grep -v -E '^\\s*\\!|^\\s*\$' | perl -p -e 's/\\s*\\!.+\$//;s/\£/\#/' > $wd/../PostCPP/$srcf";
#    #my $cmd = $cmd_cpp.$cmd_clean_up
#    #say $cmd if $VV;
#    #system($cmd);
#    my $out_path = "$wd/../PostCPP/$src";
#    run_cpp_and_clean_up( $no_macros, $includestr,  $definestr, $src_path, $out_path );
#}

for my $srcdir (@srcdirs) {
    if (not -d "$wd/../PostCPP/$srcdir") {
        system("mkdir -p $wd/../PostCPP/$srcdir");
    }
    if (not -d "$wd/../PostCPP/PrePostCPP/$srcdir") {
        system("mkdir -p $wd/../PostCPP/PrePostCPP/$srcdir");
    }    
    chdir "$wd/$srcdir";
#    my @srcs = glob("*.f95");
    my @srcs = $single ? ( $single_src ) : glob("*.f95");
    for my $src (@srcs) {
        my $src_path = "$wd/../PostCPP/PrePostCPP/$srcdir/$src";
        $stash_ref = stash_lines_to_skip($src,$stash_ref,\@macros_to_skip,$src_path);
        my $out_path = $single ? '' : "$wd/../PostCPP/$srcdir/$src";
        run_cpp_and_clean_up( $no_macros, $includestr,  $definestr, $src_path, $out_path );

#        my $cmd_cpp = "cpp -Wno-invalid-pp-token -P  $includestr_l1 $definestr -Wno-extra-tokens $wd/../PostCPP/PrePostCPP/$srcdir/$srcf ";
#        my $cmd_clean_up = " | grep -v -E '^\\s*\\!|^\\s*\$' | perl -p -e 's/\\s*\\!.+\$//;s/\£/\#/' > $wd/../PostCPP/$srcdir/$srcf";
    }
    chdir $wd;
}

open my $STASH, '>', 'stash.pl';
print $STASH Dumper($stash_ref);
close $STASH;

#$src_path = $wd/../PostCPP/PrePostCPP/$srcdir/$srcf
sub run_cpp_and_clean_up { (my $no_macros, my $includestr, my $definestr, my $src_path, my $out_path) = @_;
    my $cmd_cpp =  $no_macros ? "cat $src_path " : "cpp -Wno-invalid-pp-token -P  $includestr $definestr -Wno-extra-tokens $src_path ";
    my $redir = $out_path eq '' ? '' : '>';
    my $cmd_clean_up = "| grep -v -E '^\\s*\\!|^\\s*\$' | perl -p -e 's/\\s*\\!.+\$//;s/\£/\#/' $redir $out_path";#$wd/../PostCPP/$src_path";
#    my $cmd_clean_up = "| grep -v -E '^\\s*\\!|^\\s*\$' | perl -p -e 's/\\s*\\!.+\$//;s/\£/\#/' > $wd/../PostCPP/$srcdir/$srcf";
#    my $cmd_cpp = "cpp -Wno-invalid-pp-token -P  $includestr_l1 $definestr -Wno-extra-tokens $wd/../PostCPP/PrePostCPP/$srcdir/$srcf ";
        my $cmd = $cmd_cpp.$cmd_clean_up;
        say $cmd if $VV;
        system( $cmd );  
}

#die Dumper($stash_ref );
# To support this we need extra work; to do it right quite a bit of extra work.
# I can't find anything decent on CPAN either so we'll have to do it from scratch
# Basically I need to keep a stack for the #if/#endif blocks
# Then go through all lines, if we find a line with a matching macro, no matter if it is #if, #ifdef or #ifndef, we push all lines in the block onto @stored_lines
# We replace them by a "7188 continue" as in the other code
sub stash_lines_to_skip { (my $srcfile, my $stash_ref, my $macros_list, my $output_path) =@_;
    # $output_path='../PostCPP/PrePostCPP/'.$srcdir
    my %macros_to_skip = map {$_ => 1} @{ $macros_list };
    my @blocks_stack=();
    my @stored_lines=();
    my $placeholder = 7188;
    my @out_lines=();
    my $in_block='';
    $stash_ref->{$srcfile}={};
    open my $IN, '<', $srcfile or die $!;
    while (my $line = <$IN>) {
        chomp $line;
        if ($line=~/^\s*\#if\w*\s+(.+)\s*$/) {
            # put on the stack
            my $macro = $1;
            # for simple MACRO == ... and similar we attempt to detect them
            # FIXME: we should actually detect all macros on the line and test if any of them is in the skip list
            $macro=~s/\s*[\!\=\>\<]+.+//;
            push @blocks_stack, $macro;
            if (exists $macros_to_skip{$macro}) { 
                # put the placeholder line
                $in_block=$macro;
                push @out_lines, '  '.$placeholder.' continue';
                $stash_ref->{$srcfile}{$placeholder}=[];
                # push the line onto the list of lines for this placeholder
                push @{ $stash_ref->{$srcfile}{$placeholder}}, $line;
            } else {
            if ($in_block ne '') {
                push @{ $stash_ref->{$srcfile}{$placeholder}}, $line;
            } else {
                # do nothing
                push @out_lines, $line;
            }
            }

        } elsif ($line=~/^\s*\#endif/) { 
            # pop off the stack
            my $macro = pop @blocks_stack;
            # if this line is the end of a block, resume normal operation
            if (scalar @blocks_stack == 0 and $in_block eq $macro and $macro ne '') {
                push @{ $stash_ref->{$srcfile}{$placeholder}}, $line;
                $in_block='';
                $placeholder++;
            } else {
                if ($in_block ne '') {
                    push @{ $stash_ref->{$srcfile}{$placeholder}}, $line;
                } else {
                    # do nothing
                    push @out_lines, $line;
                }
            }
        } else {
            if ($in_block ne '') {
            # if we're in a block, stash
             push @{ $stash_ref->{$srcfile}{$placeholder}}, $line;
            } else {
            # otherwise, do nothing
             push @out_lines, $line;

            }
        }
    }
    close $IN;
    open my $OUT, '>', $output_path or die $!;
    map {say $OUT $_} @out_lines;
    close $OUT;
    if (not exists  $stash_ref->{$srcfile}{--$placeholder}) {
        delete $stash_ref->{$srcfile};
    }
    return $stash_ref;
}




