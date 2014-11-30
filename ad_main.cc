//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
// For used_time()
#include <limits.h>      // For CLK_TCK.
#include <sys/times.h>   // For times().

#include "ad_graph.h"
#include "ad_util.h"

int ( *dist_func )( int, int );

int
main( int argc, char *argv[] )
{
    // Read the initial graph into g.
    ad_graph< ninfo > g;
    ginfo             gi;
    args_t            args;

    parse_args( argc, argv, args );

    float total_time;

    // Read or generate input graph.
    {
        switch ( args.which_dist ) {
        case 0:
            dist_func = uniform_dist;
            break;
        case 1:
            dist_func = normal_dist;
            break;
        case 2:
            dist_func = exp_dist;
            break;
        }

        switch ( args.mode ) {
        case 0:
        case 1:
            {
                total_time = used_time();
                g.read( gi, args );
                total_time = used_time() - total_time;
                printf( "time to read input graph= %10.2f\n", total_time );
            }
            break;
        case 2:
            {
                args.seed = randomize( args.seed );
                printf( "seed= %ld\n", args.seed );
                total_time = used_time();
                g.generate_all( gi, args );
                total_time = used_time() - total_time;
                printf( "time to generate input graph= %10.2f\n", total_time );
                g.fprint( args );
            }
            break;
        }
    }

#ifdef PRINT_GRAPH
    g.print();
#endif

    ad_graph< cninfo >  cg;

    // If g is acyclic, the SCCs of the cg are not created, so no need
    // to call clear_components().
    bool is_acyclic;

    // Find components.
    {
        total_time = used_time();
        if ( 2 != args.mode )
            is_acyclic = find_components( cg, g, gi.has_self_loop );
        else
            is_acyclic = find_components( cg, g, gi.has_self_loop, true );
        total_time = used_time() - total_time;
        printf( "time to find components= %10.2f\n", total_time );

#ifdef PRINT_GRAPH
        print_components( cg );
        printf( "\nThe component graph: \n" );
        cg.print();
#endif
    }

    // If g is acyclic, we already know its optimum cycle mean.
    if ( is_acyclic ) {

        if ( args.min_version )
            printf( "final min_lambda= infinity time= 0.00\n" );
        else
            printf( "final max_lambda= -infinity time= 0.00\n" );

    } else { 

        // If g is cyclic, we have to compute its optimum cycle mean via cg.

        float lambda;

        if ( args.min_version ) {
            for ( int run_no = 0; run_no < args.nruns; ++run_no ) {
                printf( "run_no= %d\n", run_no );

                if ( run_no && ( 0 != args.mode ) ) {
                    args.seed += run_no;
                    args.seed = randomize( args.seed );
                    generate_part_for_all_components( cg, gi, args );
                }

#ifdef REP_COUNT
                begin_count();
#endif
                total_time = used_time();
                lambda = find_min_cycle_ratio_for_components( cg, gi.total_edge_weight );
                total_time = used_time() - total_time;
#ifdef REP_COUNT
                print_count();
#endif
                printf( "final min_lambda= %10.2f time= %10.2f\n", lambda, total_time );
            }
        } else {
            for ( int run_no = 0; run_no < args.nruns; ++run_no ) {
                printf( "run_no= %d\n", run_no );

                if ( run_no && ( 0 != args.mode ) ) {
                    args.seed += run_no;
                    args.seed = randomize( args.seed );
                    generate_part_for_all_components( cg, gi, args );
                }

#ifdef REP_COUNT
                begin_count();
#endif
                total_time = used_time();
                lambda = find_max_cycle_ratio_for_components( cg, gi.total_edge_weight );
                total_time = used_time() - total_time;
#ifdef REP_COUNT
                print_count();
#endif
                printf( "final max_lambda= %10.2f time= %10.2f\n", lambda, total_time );
            }
        }

#ifdef REP_COUNT
        end_count();
#endif

        clear_components( cg );

    } // if cyclic

    return 0;
}  // main

// End of file
