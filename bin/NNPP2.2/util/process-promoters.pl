# process-promoters.pl
# by Nomi Harris
# Copyright (c) 1996-8 The Regents of the University of California.
#
# Required by process-signals
# Returns @promoter structure
# Required argument:  array of promoter hits
# Optional argument:  If $highest_only is nonzero, save only highest-scoring
# hit for each overlapping interval (+/- $highest_only)
sub ProcessPromoters {
    local ( *hits, $highest_only, $min_score ) = @_;
    local ( $i, $h, @promoter, $_ );

    # Clear fields
    @promoter = ();

    $h = -1;    # We're going to increment $h before saving the first hit
    for ( $i = 0 ; $i <= $#hits ; $i++ ) {
        $_ = $hits[$i];

        if (/window size = ([0-9]+)/i) {
            $window_size = $1;
            if ($DEBUG) { printf( stderr "Window size = $window_size\n" ); }
        }
        elsif (/threshold = ([0-9\.]+)/i) {
            $threshold = $1;
        }
        elsif (/Prediction for ([^ ]*) \(([0-9]+) bases/) {
            $seq_name   = $1;
            $seq_length = $2;
        }

        # Look for hit
        elsif (
            /^Hit ([0-9]+): position.* ([0-9]+) \.\. ([0-9]+), (.*) = ([0-9]+)/i
        ) {
            ++$h;
            $promoter[$h]{start}  = $2;
            $promoter[$h]{end}    = $3;
            $promoter[$h]{signal} = $5;
            if ( $1 == 1 ) {   # Identify which kind of signal we're looking for
                if ( $4 =~ /exon end/ ) {
                    $IS_DONOR = 1;
                }
                elsif ( $4 =~ /exon start/ ) {
                    $IS_ACCEPTOR = 1;
                }
                else {
                    $IS_PROMOTER = 1;
                }
            }
        }
        elsif (/Prediction: ([0-9\.]+)/) {
            $promoter[$h]{score} = $1;
        }
    }

    if ($highest_only) {
        return (
            &RemoveOverlappingHits( *promoter, $seq_length, $highest_only ) );
    }
    else {
        return (@promoter);
    }
}

# 9/8 New method
# Sort promoters by score.
# Take best-scoring promoter (or promoters, if tie).
# Get rid of all promoters that are +/-7 away from best promoter.
# Take next-best promoter(s), repeat process.
# At the end, we should have the best promoter in each region.
sub RemoveOverlappingHits {
    local ( *promoter, $seq_length, $neighborhood ) = @_;
    local ( @forward, @reverse, $p, @sorted_by_score, @best, @best_reverse );

    undef @forward;
    undef @reverse;

    # 7 is default neighborhood (promoters can't be within +/-7 of e/o)
    $neighborhood = ( $neighborhood > 1 ) ? $neighborhood : 7;

    # Separate promoters into forward and reverse
    @forward = @promoter;
    for ( $p = 0 ; $p <= $#promoter ; $p++ ) {

        # Look for the first reverse hit
        if ( $promoter[$p]{start} > $promoter[$p]{end} ) {
            @reverse = @promoter;
            splice( @reverse, 0, $p );    # ?
            splice( @forward, $p );
            last;
        }
    }

    if ($DEBUG) {
        print("Forward strand:\nStart\tEnd\tScore\n");
        for ( $p = 0 ; $p <= $#forward ; $p++ ) {
            print(
                "$forward[$p]{start}\t$forward[$p]{end}\t$forward[$p]{score}\n"
            );
        }

        print("Reverse strand:\nStart\tEnd\tScore\n");
        for ( $p = 0 ; $p <= $#reverse ; $p++ ) {
            print(
                "$reverse[$p]{start}\t$reverse[$p]{end}\t$reverse[$p]{score}\n"
            );
        }
    }

    undef @sorted_by_score;
    @sorted_by_score = sort byscore @forward;
    if ($DEBUG) {
        print("Forward strand promoters sorted by score:\nStart\tScore\n");
        for ( $p = 0 ; $p <= $#sorted_by_score ; $p++ ) {
            print("$sorted_by_score[$p]{start}\t$sorted_by_score[$p]{score}\n");
        }
    }

    while ( $#sorted_by_score >= 0 ) {

        # pick_best_promoters destructively modifies both arrays!
        &pick_best_promoters( *sorted_by_score, *best, $neighborhood );
    }

    # Resort promoters by start position
    @best = sort bystart @best;

    # Now do reverse strand, if any
    if ( $#reverse >= 0 ) {
        @sorted_by_score = sort byscore @reverse;

        if ($DEBUG) {
            print("Reverse strand promoters sorted by score:\nStart\tScore\n");
            for ( $p = 0 ; $p <= $#sorted_by_score ; $p++ ) {
                print(
                    "$sorted_by_score[$p]{start}\t$sorted_by_score[$p]{score}\n"
                );
            }
        }

        while ( $#sorted_by_score >= 0 ) {

            # pick_best_promoters destructively modifies both arrays!
            &pick_best_promoters( *sorted_by_score, *best_reverse,
                $neighborhood );
        }

        # Resort promoters by start position
        @best_reverse = sort bystart @best_reverse;

        push( @best, @best_reverse );
    }

    return (@best);
}

# Sort in descending order by score
sub byscore {
    $b->{score} <=> $a->{score};
}

# Sort in ascending order by start
sub bystart {
    $a->{start} <=> $b->{start};
}

# Destructively modifies both arrays!
# @prom is already sorted by score
sub pick_best_promoters {
    local ( *prom, *keep, $neighborhood ) = @_;
    local ( $p, $j, $k );

    # Pick best promoter (first on list) and any that tie with it
    for (
        $p = 0 ;
        ( $p <= $#prom ) && ( $prom[$p]{score} == $prom[0]{score} ) ;
        $p++
    ) {
        if ($DEBUG) {
            print(
"\nSaving current best promoter (start = $prom[$p]{start}, score = $prom[$p]{score})\n"
            );
        }
        push( @keep, $prom[$p] );
        splice( @prom, $p, 1 );    # Get rid of this one now that it's on "keep" list

        # Remove from consideration all remaining promoters in
        # neighborhood of best one.
        for ( $j = $p ; $j <= $#prom ; $j++ ) {
            if (
                &in_neighborhood(
                    $prom[$j]{start}, $keep[-1]{start}, $neighborhood
                )
            ) {
                if ($DEBUG) {
                    print(
"Promoter $j ($prom[$j]{start}) is in neighborhood of current best ($keep[-1]{start})\n"
                    );
                    print(
"Remaining $#prom+1 promoters after removing $prom[$j]{start}: "
                    );
                }
                splice( @prom, $j, 1 );
                if ($DEBUG) {
                    for ( $k = 0 ; $k <= $#prom ; $k++ ) {
                        print("$prom[$k]{start} ");
                    }
                    print("\n");
                }
                $j--; # It's bad to tamper with loop variables inside a loop, but too bad.
            }
        }
    }
}

sub in_neighborhood {
    local ( $pos1, $pos2, $neighborhood ) = @_;

    if ( abs( $pos1 - $pos2 ) <= $neighborhood ) {
        return (1);
    }
    else {
        return (0);
    }
}

# Return the complement of a sequence
sub complement {
    local (*seq) = @_;

    for ( $c = 0 ; $c <= $#seq ; $c++ ) {
        if ( $seq[$c] =~ /[Aa]/ ) {
            $seq[$c] = "T";
        }
        elsif ( $seq[$c] =~ /[Cc]/ ) {
            $seq[$c] = "G";
        }
        elsif ( $seq[$c] =~ /[Gg]/ ) {
            $seq[$c] = "C";
        }
        elsif ( $seq[$c] =~ /[Tt]/ ) {
            $seq[$c] = "A";
        }
    }

    return (@seq);
}

1;    # So this file can be included in other programs
