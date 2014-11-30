//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
//
// An implementation of Howard's value iteration algorithm for the
// minimum mean cycle problem.
//

// Parameters of interest: NITER, NCYCLES, CYCLELEN, CHECK_LIMIT,
// NUPDATES, not_improved.

// count[0] = number of iterations to finish
// count[1] = number of iterations to find final lambda

#include "ad_graph.h"
#include "ad_queue.h"

#define BELLMAN_FORD_LIKE
//#define MAKE_POLICY_CONNECTED

// More node info for Howard's algorithm. For a node u,
// more_ninfo[u].policy = (u, more_ninfo[u].target). Carrying target,
// einfo (weight of policy edge), and einfo2 (transit time) is
// redundant but it makes running time faster by eliminating access to
// the edge list.
struct ninfo_how {
    float dist;    // node potential.
    int   visited; // set if node is visited for some purpose.
    int   policy;  // successor edge.
    int   target;  // successor node
    int   einfo;   // weight of policy edge.
#ifndef CYCLE_MEAN_VERSION
    int   einfo2;  // transit time of policy edge.
#endif
};

float
find_min_cycle_ratio_for_scc( const ad_graph< ninfo > *g, 
                              int plus_infinity,
                              float lambda_so_far )
{
    int n = g->num_nodes();
    int m = g->num_edges();

    ninfo_how        *more_ninfo = new ninfo_how[ n ];

    float f_plus_infinity = ( float ) plus_infinity;

    // STEP: Find the initial policy graph:
    for ( int v = 0; v < n; ++v )
        more_ninfo[ v ].dist = f_plus_infinity;

    for ( int e = 0; e < m; ++e ) {
        int u = g->source( e );
        int d = g->edge_info( e );

        if ( d < more_ninfo[ u ].dist ) {
            more_ninfo[ u ].dist = ( float ) d;
            more_ninfo[ u ].policy = e;
            more_ninfo[ u ].target = g->target( e );
            more_ninfo[ u ].einfo = d;
#ifndef CYCLE_MEAN_VERSION
            more_ninfo[ u ].einfo2 = g->edge_info2( e );
#endif
        }
    }

    float lambda = lambda_so_far;

    int CHECK_LIMIT = n;
    int CHECK_COUNT = 0;

    while ( true ) {

#ifdef REP_COUNT
        count[ 0 ]++;
#endif
#ifdef REP_COUNT_PRINT
        printf( "REP_COUNT Iteration number= %d lambda= %10.2f\n", count[ 0 ], lambda );
#endif

        // STEP: Find the min mean cycle in the policy graph. Note that
        // each connected component in the policy graph has exactly one
        // cycle.
        for ( int v = 0; v < n; ++v )
            more_ninfo[ v ].visited = -1;

#ifdef PROGRESS
        int NCYCLES = 0;
        int CYCLELEN = 0;
#endif

        // At the exit of this loop, visited field of every node must be >
        // nonnegatice (-1 to be exact).
        for ( int v = 0; v < n; ++v ) {

            if ( 0 <= more_ninfo[ v ].visited )
                continue;

            // Search for a new cycle:
            int u = v;
            do {
                more_ninfo[ u ].visited = v;
                u = more_ninfo[ u ].target;
            } while ( -1 == more_ninfo[ u ].visited );

            if ( v != more_ninfo[ u ].visited )
                continue;

#ifdef PROGRESS
            NCYCLES++;
#endif

            // Compute the mean of the cycle found. Note that u is a node on
            // this cycle.
            int w = u;
            int total_weight = 0;
            int total_length = 0;
            do {
#ifdef CYCLE_MEAN_VERSION
                ++total_length;
#else
                total_length += more_ninfo[ u ].einfo2;
#endif
                total_weight += more_ninfo[ u ].einfo;
                u = more_ninfo[ u ].target;
            } while ( u != w );

            // Update lambda only if it decreases.
            float new_lambda = ( float ) total_weight / total_length;
            if ( new_lambda < lambda ) {
                lambda = new_lambda;
#ifdef PROGRESS
                CHECK_COUNT = 0;
#endif
#ifdef REP_COUNT
                count[ 1 ] = count[ 0 ];
#endif
            }

#ifdef PROGRESS
            CYCLELEN += total_length;
#endif
        } // for v

#ifdef PROGRESS
        printf( "COUNTERS NCYCLES= %d CYCLELEN= %d\n", NCYCLES, CYCLELEN);
        int NUPDATES = 0;
#endif

        if ( CHECK_COUNT++ > CHECK_LIMIT ) {
#ifdef PROGRESS
            printf( "COUNTERS reason to exit: CHECK_LIMIT\n" );
#endif
            break;
        }

        // STEP: Update the dist of the other nodes:
        bool not_improved = true;

#ifdef BELLMAN_FORD_LIKE
        for ( int e = 0; e < m; ++e ) {
            int u = g->source( e );
            int v = g->target( e );

#ifdef CYCLE_MEAN_VERSION
            float new_dist = more_ninfo[ v ].dist + 
                g->edge_info( e ) - lambda;
#else
            float new_dist = more_ninfo[ v ].dist + 
                g->edge_info( e ) - lambda * g->edge_info2( e );
#endif
            if ( EPSILON < ( more_ninfo[ u ].dist - new_dist ) ) {
                not_improved = false;
                more_ninfo[ u ].dist = new_dist;
                more_ninfo[ u ].policy = e;
                more_ninfo[ u ].target = v;
                more_ninfo[ u ].einfo = g->edge_info( e );
#ifndef CYCLE_MEAN_VERSION
                more_ninfo[ u ].einfo2 = g->edge_info2( e );
#endif
#ifdef PROGRESS
                ++NUPDATES;
#endif
            }
        }
#else
        for ( int u = 0; u < n; ++u ) {
            float d = more_ninfo[ u ].dist;
            int which = -1;
            for ( int i = 0; i < g->outdegree( u ); ++i ) {
                int v = g->ith_target_node( u, i );

#ifdef CYCLE_MEAN_VERSION
                float new_dist = more_ninfo[ v ].dist + 
                    g->ith_target_edge_info( u, i ) - lambda;
#else
                float new_dist = more_ninfo[ v ].dist + 
                    g->ith_target_edge_info( u, i ) - lambda * g->ith_target_edge_info2( u, i ); 
#endif
                if ( EPSILON < ( d - new_dist ) ) {
                    not_improved = false;
                    d = new_dist;
                    which = i;
#ifdef PROGRESS
                    ++NUPDATES;
#endif
                }
            }
            if ( -1 != which ) {
                more_ninfo[ u ].dist = d;
                more_ninfo[ u ].policy = g->ith_target_edge( u, which );
                more_ninfo[ u ].target = g->ith_target_node( u, which );
                more_ninfo[ u ].einfo = g->ith_target_edge_info( u, which );
#ifndef CYCLE_MEAN_VERSION
                more_ninfo[ u ].einfo2 = g->ith_target_edge_info2( u, which );
#endif
            }
        }
#endif

#ifdef PROGRESS
        printf( "COUNTERS NUPDATES= %d\n", NUPDATES );
#endif

        if ( not_improved ) {
#ifdef PROGRESS
            printf( "COUNTERS reason to exit: not_improved\n" );
#endif
            break;
        }
    }  // main while loop

#ifdef REP_COUNT
    count[ 0 ]++;
#endif
#ifdef REP_COUNT_PRINT
    printf( "REP_COUNT Iteration number= %d lambda= %10.2f\n", count[ 0 ], lambda );
#endif

    delete [] more_ninfo;

    return lambda;
}  // find_min_cycle_ratio_for_scc

#undef BELLMAN_FORD_LIKE
#undef MAKE_POLICY_CONNECTED

// End of file
