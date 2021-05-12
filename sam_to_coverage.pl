#!/usr/bin/perl

my $version = "0.2";

use strict;
use warnings;

$|++;

my $usage = "
sam_to_coverage.pl [options] -f <fasta.fai> <alignment.sam>

Output coordinates of reference bases covered in each alignment, as well as
some coverage statistics.

Required:
  -f    fasta.fai file produced by \"samtools faidx\" command
  
Options:
  -o    if selected, will simply output a file of comma separated values of
        coverages at each position in the alignment, regardless of the setting
        for -d. This file can then be ater processed by this script to output
        coverage coordinates and statistics at various depths
  -c    if selected, indicates the file given to <alignment.sam> is a
        coverage file produced earlier using option -o above
  -r    file with coordinate ranges. If given, will output coverage statistics
        over coordinates in the given range.
        File should have the following format:
          contig_id<tab>start<tab>end<tab>segment_id
          \"segment_id\" can be any string that can be used to identify the
          region in the output
  -d    minimum coverage depth (default 1)
  -v    verbose output (default is off)

";

## command line processing.
use Getopt::Std;
use vars qw( $opt_f $opt_d $opt_v $opt_o $opt_c $opt_r);
getopts('f:d:vocr:');
die $usage unless ($opt_f and @ARGV);

my $fai     = $opt_f if $opt_f;
my $minc    = $opt_d ? $opt_d : 1;
my $r_file  = $opt_r if $opt_r;

#Read in the SAM file and store coverages in a hash
my %covs;
my $total_mapped_reads = 0;
my $id;
open (my $in, "<$ARGV[0]") or die "Can't open $ARGV[0]: $!\n";
while (my $line = <$in>){
    chomp $line;
    if ($opt_c){
        if ($line =~ m/^>/){
            $id = substr($line, 1);
            next;
        }
        my @tmp = split(",", $line);
        for my $i (0 .. $#tmp){
            my $val = $tmp[$i];
            next if $val == 0;
            $covs{$id}{$i+1} = $val;
        }
        next;
    }
    next if $line =~ m/^@/;
    my @tmp = split("\t", $line);
    my ($rid, $start, $cigar) = ($tmp[2], $tmp[3], $tmp[5]);
    next if $cigar eq "*"; #unmapped reads not counted
    $total_mapped_reads++;
    print STDERR "\rchecking reads: $total_mapped_reads" if ($total_mapped_reads % 5000 == 0 and $opt_v);
    if ($cigar =~ m/[NP=X]/){
        die "Died at $cigar: contains unusual CIGAR operation.\n";
    }
    my @cigars = split(/(?<=\D)/, $cigar);
    #if ($flag & 16){ #if the read was aligned in reverse...
    #    @cigars = reverse @cigars;
    #}
    my $pos = $start;
    for my $i (0 .. $#cigars){
        my $tmp = $cigars[$i];
        if ($tmp =~ m/(\d+)([MIDSH])/){
            my ($num, $op) = ($1, $2);
            next if $op eq "S"; #soft-clipped bases don't change the position at all
            next if $op eq "H"; #hard-clipped bases don't change the position at all
            next if $op eq "I"; #insertions relative to the reference don't change the position at all
            my $end = ($pos + $num) - 1;
            unless ($op eq "D"){ #for deletions relative to the reference, advance position without adding coverage to the deleted bases
                for my $j ($pos .. $end){ #add matching ("M") bases to coverage
                    $covs{$rid}{$j}++;
                }
            }
            $pos = $end + 1;
        }
    }
}
close ($in);
print STDERR "\n" if $opt_v;

#read in the ranges file, if given
my %ranges;
if ($r_file){
    open (my $rin, "<$r_file") or die "ERROR: Can't open $r_file: $!\n";
    while (my $line = <$rin>){
        chomp $line;
        next if $line =~ m/^\s*$/;
        my ($contig, $start, $stop, $segid) = split("\t", $line);
        push @{$ranges{$contig}}, ([$start, $stop, $segid]);
    }
    close ($rin);
}

#Read in the index file and calculate coverages
open (my $fin, "<$fai") or die "Can't open $fai: $!\n";
print "Segment ID\tlength\t#of bases with < $minc coverage\t% of bases covered\taverage coverage\tmedian coverage\tmax coverage\tcovered_coords\n" if (!$opt_o);
while (my $line = <$fin>){
    chomp $line;
    my ($id, $leng) = split("\t", $line);
    print ">$id\n" if $opt_o;
    my $seg_start;
    my $max_cov = 0;
    my $num_cov = 0;
    my @coords;
    my $tot_depths = 0;
    my @cvg;
    my @per_base_covs;
    for my $i (1 .. $leng){
        if ($covs{$id}{$i}){
            my $depth = $covs{$id}{$i};
            if ($opt_o){
                push @cvg, $depth;
                next;
            }
            $max_cov = $depth if $depth > $max_cov;
            if ($depth >= $minc){
                $tot_depths += $depth;
                push @per_base_covs, $depth;
                if (!$seg_start){
                    $seg_start = $i;
                }
                next;
            } else {
                push @per_base_covs, 0;
            }
            if ($seg_start){
                my $last = $i - 1;
                push @coords, "$seg_start - $last";
                $num_cov += (($last - $seg_start) + 1);
                $seg_start = "";
                next;
            }
        } else {
            if ($opt_o){
                push @cvg, 0;
                next;
            }
            push @per_base_covs, 0;
            if ($seg_start){
                my $last = $i - 1;
                push @coords, "$seg_start - $last";
                $num_cov += (($last - $seg_start) + 1);
                $seg_start = "";
                next;
            }
        }
    }
    if ($opt_o){
        print join(",",@cvg), "\n";
        next;
    }
    if ($seg_start){
        my $last = $leng;
        push @coords, "$seg_start - $last";
        $num_cov += (($last - $seg_start) + 1);
        $seg_start = "";
    }
    my $num_less = $leng - $num_cov;
    my $pct_cov = sprintf("%.2f", 100*($num_cov / $leng));
    my $avg_depth = sprintf("%.2f", 0);
    #$avg_depth = sprintf("%.2f", ($tot_depths / $num_cov)) if $num_cov > 0;
    $avg_depth = sprintf("%.2f", ($tot_depths / $leng));
    my $median = 0;
    if (@per_base_covs > 0){
        @per_base_covs = sort{$a <=> $b}@per_base_covs;
        my $cov_num = scalar @per_base_covs;
        if ($cov_num % 2 == 0){
            $median = ($per_base_covs[$cov_num/2 - 1] + $per_base_covs[$cov_num/2]) / 2;
        } else {
            $median = ($per_base_covs[$cov_num/2 - 0.5]);
        }
    }
    print "$id\t$leng\t$num_less\t$pct_cov\t$avg_depth\t$median\t$max_cov";
    if (@coords){
        print "\t", join(", ", @coords);
    }
    print "\n";
    #determine and output range statistics, if range file was given
    if ($ranges{$id}){
        foreach my $slice (@{$ranges{$id}}){
            my ($start, $stop, $segid) = @{$slice};
            $leng = $stop - $start + 1;
            my $seg_start;
            my $max_cov = 0;
            my $num_cov = 0;
            my @coords;
            my $tot_depths = 0;
            my @per_base_covs;
            for my $i ($start .. $stop){
                if ($covs{$id}{$i}){
                    my $depth = $covs{$id}{$i};
                    $max_cov = $depth if $depth > $max_cov;
                    if ($depth >= $minc){
                        $tot_depths += $depth;
                        push @per_base_covs, $depth;
                        if (!$seg_start){
                            $seg_start = $i;
                        }
                        next;
                    } else {
                        push @per_base_covs, 0;
                    }
                    if ($seg_start){
                        my $last = $i - 1;
                        push @coords, "$seg_start - $last";
                        $num_cov += (($last - $seg_start) + 1);
                        $seg_start = "";
                        next;
                    }
                } else {
                    push @per_base_covs, 0;
                    if ($seg_start){
                        my $last = $i - 1;
                        push @coords, "$seg_start - $last";
                        $num_cov += (($last - $seg_start) + 1);
                        $seg_start = "";
                        next;
                    }
                }
            }
            if ($seg_start){
                my $last = $stop;
                push @coords, "$seg_start - $last";
                $num_cov += (($last - $seg_start) + 1);
                $seg_start = "";
            }
            my $num_less = $leng - $num_cov;
            my $pct_cov = sprintf("%.2f", 100*($num_cov / $leng));
            my $avg_depth = sprintf("%.2f", 0);
            #$avg_depth = sprintf("%.2f", ($tot_depths / $num_cov)) if $num_cov > 0;
            $avg_depth = sprintf("%.2f", ($tot_depths / $leng)); 
            my $median = 0;
            if (@per_base_covs > 0){
                @per_base_covs = sort{$a <=> $b}@per_base_covs;
                my $cov_num = scalar @per_base_covs;
                if ($cov_num % 2 == 0){
                    $median = ($per_base_covs[$cov_num/2 - 1] + $per_base_covs[$cov_num/2]) / 2;
                } else {
                    $median = ($per_base_covs[$cov_num/2 - 0.5]);
                }
            }
            print "$segid\t$leng\t$num_less\t$pct_cov\t$avg_depth\t$median\t$max_cov";
            if (@coords){
                print "\t", join(", ", @coords);
            }
            print "\n";
        }
    }
}
close ($fin);
