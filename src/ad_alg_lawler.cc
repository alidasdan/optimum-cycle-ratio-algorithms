//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
//
// An implementation of Lawler's minimum cycle ratio algorithm using
// Bellman-Ford's shortest path algorithm.
//

#include "ad_graph.h"
#include "ad_cqueue.h"

// More node info for Lawler's algorithm.
struct ninfo_lawler {
    float dist;
    int   not_included;
};

float
find_min_cycle_ratio_for_scc( const ad_graph< ninfo > *g, 
                              int plus_infinity,
                              float lambda_so_far )
{
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

    ninfo_lawler      *more_ninfo = new ninfo_lawler[ n ];
    ad_cqueue< int >  nodeq( n + 1 );  // +1 for insertion of END_PHASE node.

#define END_PHASE -1

    while ( ( upper - lower ) > EPSILON ) {

        // Determine the new lambda in the middle.
        lambda = ( upper + lower ) / 2;

#ifdef PROGRESS
        ++NITER;
#endif

#ifdef PROGRESS
        printf( "PROGRESS Iteration number= %d lambda= %10.2f\n", NITER, lambda );
#endif

        // Subtract lambda from each edge weight.
        // Check to see if the resulting graph has a negative cycle.

        more_ninfo[ SOURCE ].dist = 0;
        more_ninfo[ SOURCE ].not_included = 0;
        for ( int v = 1; v < n; ++v ) {
            more_ninfo[ v ].dist = f_plus_infinity;
            more_ninfo[ v ].not_included = 1;
        }

        nodeq.init();
        nodeq.put( SOURCE );
        nodeq.put( END_PHASE );

        bool found = true;
        int nphase = 0;

        while ( nphase < n ) {
            int u = nodeq.get();

            if ( END_PHASE == u ) {
                nphase++;

                if ( nodeq.is_empty() ) {
                    found = false;
                    break;
                }

                nodeq.put( END_PHASE );
                continue;
            }
            else 
                more_ninfo[u].not_included = 1;

            float udist = more_ninfo[ u ].dist;

            for ( int i = 0; i < g->outdegree( u ); ++i ) {
                int v = g->ith_target_node( u, i );

#ifdef CYCLE_MEAN_VERSION
                float new_dist = udist + g->ith_target_edge_info( u, i ) - lambda;
#else
                float new_dist = udist + g->ith_target_edge_info( u, i ) 
                    - lambda * g->ith_target_edge_info2( u , i );
#endif
                if ( new_dist < more_ninfo[ v ].dist ) {
                    more_ninfo[ v ].dist = new_dist;
                    if ( more_ninfo[ v ].not_included ) {
                        more_ninfo[ v ].not_included = 0;
                        nodeq.put( v );
                    }
                }
            } // for

        } // while nphase > n

        if ( found ) {
            if ( ( upper - lambda ) < EPSILON2 )
                break;
            upper = lambda;
        }
        else {
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
