//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
//
// An implementation of Lawler's minimum cycle ratio algorithm using
// Szymanski's improvement.
//
//@inproceedings{Sz92,
//author = {T.~G. Szymanski},
//title = {Computing Optimal Clock Schedules}, 
//booktitle = {Proc. 29th Design Automation Conf.},
//year = 1992,
//editor = {},
//pages = {399--404},
//organization = {ACM/IEEE},
//publisher = {},
//address = "",
//month = {},
//}
//

// Parameters: NUPDATES, NCYCLES, CYCLELEN.

// count[0] = number of iterations
// count[1] = number of passes

#include "ad_graph.h"

// More node info for Szymanski' algorithm. einfo and einfo2 fields
// are redundant but needed for efficiency.
struct ninfo_szymanski {
    float dist;    // node distance or potential.
    int   pred;    // predecessor node
    int   einfo;   // weight of the edge from this node to its pred.
#ifndef CYCLE_MEAN_VERSION
    int   einfo2;  // transit time of the edge from this node to its pred.
#endif
    int   visited; // set if visited.
    bool  changed; // set if dist is changed.
};

float
find_min_cycle_ratio_for_scc( const ad_graph< ninfo > *g, 
                              int plus_infinity,
                              float lambda_so_far )
{
    const int INTERVAL = 10; // Szymanski's interval for cycle check.

    int n = g->num_nodes();
    int m = g->num_edges();

    float f_plus_infinity = ( float ) plus_infinity;

#ifdef CYCLE_MEAN_VERSION
    float lower = f_plus_infinity;
    float upper = -f_plus_infinity;

    // STEP: Determine lower and upper bounds on lambda:
    for ( int e = 0; e < m; ++e ) {
        min2( lower, ( float ) g->edge_info( e ) );
        max2( upper, ( float ) g->edge_info( e ) );
    }

#else
#ifdef IMPROVE_LAMBDA_BOUNDS
    float upper = find_min_lambda( g, plus_infinity );
    float lower = f_plus_infinity;

    // Determine lower and upper bounds on lambda. 
    for ( int e = 0; e < m; ++e ) {
        if ( 0 < g->edge_info2( e ) )
            min2( lower, ( float ) g->edge_info( e ) / g->edge_info2( e ) );
    }
#else
    float upper = 0.0;
    float lower = f_plus_infinity;

    // Determine lower and upper bounds on lambda. Note that this step
    // is different from the cycle mean case because it is difficult to
    // define a tight upper bound. We take upper to be the sum of the
    // arc weights.
    for ( int e = 0; e < m; ++e ) {
        upper += g->edge_info( e );
        if ( 0 < g->edge_info2( e ) )
            min2( lower, ( float ) g->edge_info( e ) / g->edge_info2( e ) );
    }
#endif
#endif

    if ( lambda_so_far <= lower )
        return lambda_so_far;

    min2( upper, ( float ) 2.0 * lambda_so_far - lower );

    float lambda = upper;

    ninfo_szymanski   *more_ninfo = new ninfo_szymanski[ n ];

    // STEP: Run Lawler's algorithm.
    while ( ( upper - lower ) > EPSILON ) {

        // STEP: Determine the new lambda in the middle.
        lambda = ( upper + lower ) / 2;

#ifdef REP_COUNT
        count[ 0 ]++;
#endif
#ifdef REP_COUNT_PRINT
        printf( "REP_COUNT Iteration number= %d lambda= %10.2f\n", count[ 0 ], lambda );
#endif

        // STEP: Run Szymanski's algorithm to see if the graph has a
        // negative cycle.

#ifdef PROGRESS
        int NPASS = 0;
        int NUPDATES = 0;
        int NCYCLES = 0;
        int CYCLELEN = 0;
#endif

        // neg_cycle_found is set when a negative cycle is found.
        bool neg_cycle_found;   
        {
            // STEP2: Initialize the info of each node.
            {
                ninfo_szymanski *ptr = &more_ninfo[ SOURCE ];

                ptr->dist = 0.0;
                ptr->pred = -1;  // -1 = not yet set.
                ptr->changed = true;
                for ( int v = 1; v < n; ++v ) {
                    ++ptr;
                    ptr->dist = f_plus_infinity;
                    ptr->pred = -1;
                    ptr->changed = true;
                }
            }

            neg_cycle_found = false;

            for ( int npass = 0; npass < n; ++npass ) {

#ifdef PROGRESS
                ++NPASS;
#endif

#ifdef REP_COUNT
                ++count[ 1 ];
#endif

                bool check_cycle = ( npass == n - 1 ) || ( npass % INTERVAL == 0 );

                // STEP2: Update distances. one_changed is set when at least
                // one distance changes. Also initialize visited for 
                bool one_changed = false;
                for ( int u = 0; u < n; ++u ) {
                    if ( check_cycle )
                        more_ninfo[ u ].visited = -1; // use -1 since 0 is also a node.
                    if ( more_ninfo[ u ].changed ) {
                        more_ninfo[ u ].changed = false;

#ifdef CYCLE_MEAN_VERSION
                        float reduced_dist = more_ninfo[ u ].dist - lambda;
#else
                        float udist = more_ninfo[ u ].dist;
#endif

                        for ( int i = 0; i < g->outdegree( u ); ++i ) {
                            int e_uv = g->ith_target_edge( u, i );
                            int v = g->target( e_uv );
                            int uv_info = g->edge_info( e_uv );
#ifdef CYCLE_MEAN_VERSION
                            float new_dist = reduced_dist + uv_info;
#else
                            int uv_info2 = g->edge_info2( e_uv );
                            float new_dist = udist + uv_info - lambda * uv_info2;
#endif
                            if ( new_dist < more_ninfo[ v ].dist ) {
                                more_ninfo[ v ].dist = new_dist;
                                more_ninfo[ v ].pred = u;
                                more_ninfo[ v ].einfo = uv_info;
#ifndef CYCLE_MEAN_VERSION
                                more_ninfo[ v ].einfo2 = uv_info2;
#endif
                                more_ninfo[ v ].changed = true;
                                one_changed = true;
#ifdef PROGRESS
                                ++NUPDATES;
#endif
                            }
                        }  // for i
                    }  // if changed
                }  // for u

                // STEP2: Exit if there is a negative cycle thru SOURCE.
                if ( more_ninfo[ SOURCE ].dist < 0.0 ) {
#ifdef PROGRESS
                    printf( "COUNTERS reason to exit: dist less than 0\n" );
#endif
                    neg_cycle_found = true;
                    goto update;
                }

                // STEP2: Exit if no distance has changed.
                if ( ! one_changed ) {
#ifdef PROGRESS
                    printf( "COUNTERS reason to exit: no dist changed\n" );
#endif
                    // No negative cycle.
                    neg_cycle_found = false;
                    goto update;
                }

                // STEP2: Every so often (every INTERVAL), check for a
                // negative cycle using the pred edges. Visit them in reverse
                // order.
                if ( check_cycle ) {

                    for ( int v = 0; v < n; ++v ) {

                        if ( 0 <= more_ninfo[ v ].visited )
                            continue;

                        // Search for a new cycle. visited[u] shows from which
                        // node the search started. not_okay is set when the
                        // search leads to the source via a dummy edge
                        // (non-existent in the graph).
                        int u = v;
                        do {
                            more_ninfo[ u ].visited = v;
                            u = more_ninfo[ u ].pred;
                        } while ( ( -1 != u ) && ( -1 == more_ninfo[ u ].visited ) );

                        if ( ( -1 == u ) || ( v != more_ninfo[ u ].visited ) )
                            continue;

#ifdef PROGRESS
                        NCYCLES++;
#endif

                        // Compute the mean of the cycle found. The node u is in
                        // this cycle.
                        int w = u;
                        int total_length = 0;
                        int total_weight = 0;
#ifdef DEBUG
                        int cycle_len = 0;
#endif
                        do {
#ifdef DEBUG
                            assert( -1 != u );
                            assert( cycle_len++ < n );
                            assert( more_ninfo[ u ].pred >= 0 );
#endif
#ifdef CYCLE_MEAN_VERSION
                            ++total_length;
#else
                            total_length += more_ninfo[ u ].einfo2;
#endif
                            total_weight += more_ninfo[ u ].einfo;
                            u = more_ninfo[ u ].pred;
                        } while ( u != w );
            
                        float new_lambda = ( float ) total_weight / total_length;
                        if ( new_lambda < lambda ) {
                            // There is a negative cycle, and by definition of
                            // lambda, we know that min_lambda <= new_lambda. Thus,
                            // we set lambda to new_lambda, which will become upper
                            // later.

#ifdef PROGRESS
                            printf( "COUNTERS reason to exit: negative cycle\n" );
                            CYCLELEN += total_length;    
#endif

#ifdef IMPROVE_UPPER_BOUND
                            lambda = new_lambda;
#endif
                            neg_cycle_found = true;
                            goto update;
                        }
                    }  // for v
                }  // if every so often

            }  // for npass
        }  // End of Szymanski's algorithm.

    update:

#ifdef PROGRESS
        printf( "COUNTERS NPASS= %d NUPDATES= %d NCYCLES= %d CYCLELEN= %d\n", 
                NPASS, NUPDATES, NCYCLES, CYCLELEN );
#endif

        // STEP: Update lambda depending on whether or not a negative
        // cycle is found.
        if ( neg_cycle_found ) {  // a negative cycle
            if ( ( upper - lambda ) < EPSILON2 )
                break;
            upper = lambda;
        }
        else {  // no negative cycle
            if ( ( lambda - lower ) < EPSILON2 )
                break;
            lower = lambda;
        }

    }  // while 

#undef END_PHASE

    delete [] more_ninfo;

    return lambda;
}  // find_min_cycle_ratio_for_scc

// End of file
