//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
// Using a successor graph, find the max and min realizable lambda in
// the graph g. This code is taken from ad_alg_howard.cc.

#include "ad_graph.h"

struct ninfo_init {
    int dist;
    int visited;
    int target;
    int einfo;
#ifndef CYCLE_MEAN_VERSION
    int einfo2;
#endif
};

float 
find_lambda_bound( const ad_graph< ninfo > *g, 
                   int plus_infinity, 
                   bool which )
{
    int n = g->num_nodes();
    int m = g->num_edges();

    ninfo_init *more_ninfo = new ninfo_init[ n ];

    float f_plus_infinity = ( float ) plus_infinity;
    float lambda = f_plus_infinity;

    if ( which ) {
        // Find min arc weights.
        for ( int v = 0; v < n; ++v ) {
            more_ninfo[ v ].visited = -1;
            more_ninfo[ v ].dist = plus_infinity;
        }

        for ( int e = 0; e < m; ++e ) {
            int u = g->source( e );
            int d = g->edge_info( e );
    
            if ( d < more_ninfo[ u ].dist ) {
                more_ninfo[ u ].dist = d;
                more_ninfo[ u ].target = g->target( e );
                more_ninfo[ u ].einfo = d; 
#ifndef CYCLE_MEAN_VERSION
                more_ninfo[ u ].einfo2 = g->edge_info2( e );
#endif
            }
        }
    } else {
        // Find max arc weights.
        for ( int v = 0; v < n; ++v ) {
            more_ninfo[ v ].visited = -1;
            more_ninfo[ v ].dist = -plus_infinity;
        }

        for ( int e = 0; e < m; ++e ) {
            int u = g->source( e );
            int d = g->edge_info( e );
    
            if ( d > more_ninfo[ u ].dist ) {
                more_ninfo[ u ].dist = d;
                more_ninfo[ u ].target = g->target( e );
                more_ninfo[ u ].einfo = d;      
#ifndef CYCLE_MEAN_VERSION
                more_ninfo[ u ].einfo2 = g->edge_info2( e );
#endif
            }
        }
    }

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
    
        float new_lambda = ( float ) total_weight / total_length;
        if ( new_lambda < lambda )
            lambda = new_lambda;
    } // for v

    delete [] more_ninfo;

    return lambda;
}  // find_lambda_bound

// End of file
