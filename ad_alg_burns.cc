//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
//
// An implementation of Burns's minimum cycle ratio algorithm. 
//
// @techreport{Bu91,
// author = {S.~M. Burns},
// title = {Performance Analysis and Optimization of Asynchronous Circuits},
// institution = {California Institute of Technology},
// year = 1991,
// type = {{PhD} Thesis},
// address = "",
// number = "",
// month = "",
// }
//

#include "ad_graph.h"
#include "ad_queue.h"

// More node info for Burns's algorithm.
struct ninfo_burns {
    float dist;    // node distance or potential.
    int   length;  // path length from source in topological order.
    int   indeg;   // indegree.
};

#if 0
bool search( const ad_graph< ninfo > *g, int u, bool *visited, bool *critical )
{
    visited[ u ] = true;
    for ( int i = 0; i < g->outdegree( u ); ++i ) {
        int e = g->ith_target_edge( u, i );
        if ( critical[ e ] ) {
            int v = g->ith_target_node( u, i );
            if ( visited[ v ] )
                return true;
            else
                return search( g, v, visited, critical );
        }
    }
    return false;
}
#endif

float
find_min_cycle_ratio_for_scc( const ad_graph< ninfo > *g, 
                              int plus_infinity,
                              float lambda_so_far )
{
    int n = g->num_nodes();
    int m = g->num_edges();

    // critical[e] is true if edge e is critical.

    ninfo_burns      *more_ninfo = new ninfo_burns[ n ];
    bool             *critical = new bool[ m ];
    ad_queue< int >  nodeq( n );

    float f_plus_infinity = ( float ) plus_infinity;
    float f_minus_infinity = -f_plus_infinity;

#ifdef CYCLE_MEAN_VERSION
    // STEP: Initialize lambda to the minimum of the min edge weight and
    // the previous lambda:
    float lambda = lambda_so_far;
    for ( int e = 0; e < m; ++e )
        min2( lambda, ( float ) g->edge_info( e ) );

    // STEP: Initialize the distance of each node to zero:
    for ( int v = 0; v < n; ++v )
        more_ninfo[ v ].dist = 0.0;

#else
    // STEP: Initialize node distances and lambda. critical array is
    // used here to differentiate between the edges with zero or
    // positive transit times.
    {
        // Temporarily remove all edges with positive transit time from
        // g. The remaining edges can have zero or negative transit times
        // but for our testcases, the transit time is never negative.
        for ( int e = 0; e < m; ++e ) {
            if ( 0 >= g->edge_info2( e ) )
                critical[ e ] = true;
            else
                critical[ e ] = false;
        }

        // Topologically sort the resulting graph. Note that if the
        // resulting graph is cyclic, there are cycles with zero transit
        // time and this is an error. Using the topological order, set
        // d(v) of node v to the min distance from the source.

        // Find the indegree and distance of each node and put them into
        // the nodeq.
        nodeq.init();

        for ( int v = 0; v < n; ++v ) {
            more_ninfo[ v ].indeg = 0;

            for ( int i = 0; i < g->indegree( v ); ++i ) {
                if ( critical[ g->ith_source_edge( v, i ) ] )
                    more_ninfo[ v ].indeg++;
            }

            if ( 0 == more_ninfo[ v ].indeg ) {
                more_ninfo[ v ].dist = 0.0;
                nodeq.put( v );
            } else {
                more_ninfo[ v ].dist = f_plus_infinity;
            }
        }  // for

        // Do the actual topological sorting using the previously found
        // indegrees and distances for nodes:

        int count_visited = 0;  // Number of visited nodes.
        while ( nodeq.is_not_empty() ) {
            ++count_visited;

            int u = nodeq.get();
            for ( int i = 0; i < g->outdegree( u ); ++i ) {
                if ( critical[ g->ith_target_edge( u, i ) ] ) {
                    // edge = u->v
                    int v = g->ith_target_node( u, i );

                    // dist[ v ] = min( dist[ v ], dist[ u ] + w(u,v))
                    min2( more_ninfo[ v ].dist, 
                          more_ninfo[ u ].dist + g->ith_target_edge_info( u, i ) );
                    more_ninfo[ v ].indeg--;

                    if ( 0 == more_ninfo[ v ].indeg )
                        nodeq.put( v );
                }
            }
        }  // while

        // The graph of zero transit time edges must not have a cycle. A
        // cycle shows a deadlock in the systems that this graph models.
        if ( count_visited != n ) {
            printf("ERROR: Cycles with zero transit time are not allowed.\n");
            abort();
        }
    }

    // Set lambda to the minimum of ((d(v)+w(u,v)-d(u)) / tt(u,v)) if t(u,v) > 0.
    float lambda = lambda_so_far;
    for ( int e = 0; e < m; ++e ) {
        int u = g->source( e );
        int v = g->target( e );
        if ( 0 < g->edge_info2( e ) ) {
            float delta = more_ninfo[ u ].dist + g->edge_info( e ) - more_ninfo[ v ].dist;
            min2( lambda, delta / g->edge_info2( e ) );
        }
    }
#endif

    // STEP: Iterate until the critical graph is cyclic.
    while ( true ) {

#ifdef PROGRESS
        printf( "PROGRESS Iteration number= %d lambda= %10.2f\n", NITER, lambda );
#endif

#ifdef PROGRESS
        ++NITER;
#endif

        // STEP: Find all critical edges:
        int ncrit = 0;
        for ( int e = 0; e < m; ++e ) {
            int u = g->source( e );
            int v = g->target( e );

            float delta1 = more_ninfo[ u ].dist + g->edge_info( e ) - more_ninfo[ v ].dist;
#ifdef CYCLE_MEAN_VERSION
            if ( fabs_val( lambda - delta1 ) < SMALL_EPSILON ) {
                critical[ e ] = true;
                ++ncrit;
            } else {
                critical[ e ] = false;
            }
#else 
            if ( fabs_val( lambda * g->edge_info2( e ) - delta1 ) < SMALL_EPSILON ) {
                critical[ e ] = true;
                ++ncrit;
            } else {
                critical[ e ] = false;
            }
#endif
        }  // for each edge

        // STEP: Topologically sort the critical graph:
    
        // Find the indegree and ( negative ) length of each node and put
        // them into the nodeq.
        nodeq.init();

        for ( int v = 0; v < n; ++v ) {
            more_ninfo[ v ].indeg = 0;

            for ( int i = 0; i < g->indegree( v ); ++i ) {
                if ( critical[ g->ith_source_edge( v, i ) ] )
                    more_ninfo[ v ].indeg++;
            }

            if ( 0 == more_ninfo[ v ].indeg ) {
                more_ninfo[ v ].length = 0;
                nodeq.put( v );
            } else {
                more_ninfo[ v ].length = plus_infinity;
            }
        }  // for

        // Do the actual topological sorting using the previously found
        // indegrees and lengths for nodes:

        int count_visited = 0;  // Number of visited nodes.
        while ( nodeq.is_not_empty() ) {
            ++count_visited;

            int u = nodeq.get();
            for ( int i = 0; i < g->outdegree( u ); ++i ) {
                if ( critical[ g->ith_target_edge( u, i ) ] ) {
                    // edge = u->v
                    int v = g->ith_target_node( u, i );

#ifdef CYCLE_MEAN_VERSION
                    // length[ v ] = min( length[ v ], length[ u ] - 1 )
                    min2( more_ninfo[ v ].length, more_ninfo[ u ].length - 1 );
#else
                    // length[ v ] = min( length[ v ], length[ u ] - t(u,v) )
                    min2( more_ninfo[ v ].length, more_ninfo[ u ].length - g->ith_target_edge_info2( u, i ) );
#endif

                    more_ninfo[ v ].indeg--;

                    if ( 0 == more_ninfo[ v ].indeg )
                        nodeq.put( v );
                }
            }
        }  // while

        // STEP: If the critical graph is cyclic, then the optimum lambda
        // is found, so exit.
        if ( count_visited != n )
            break;

#if 0
        /***********************/
        if ( count_visited != n ) {
            bool *visited = new bool [ n ];
            for ( int u = 0; u < n; ++u )
                visited[ u ] = false;
            bool found = false;
            for ( int u = 0; u < n; ++u ) {
                if ( ! visited[ u ] )
                    found = search( g, u, visited, critical );
            }
            if ( found )
                printf( "cycle found\n" );
            else
                printf( "no cycle found\n" );
            delete [] visited;
            break;
        }
        /***********************/
#endif

        // STEP: Find theta to update lambda as well as the distance of
        // every node:
        float theta = f_minus_infinity;
        for ( int e = 0; e < m; ++e ) {
            // e = u->v
            int u = g->source( e );
            int v = g->target( e );

#ifdef CYCLE_MEAN_VERSION
            int delta2 = more_ninfo[ v ].length + 1 - more_ninfo[ u ].length;
            if ( delta2 > 0 ) {
                float delta1 = more_ninfo[ u ].dist + g->edge_info( e ) - more_ninfo[ v ].dist;
                max2( theta, ( lambda - delta1 ) / delta2 );
            }
#else
            int delta2 = more_ninfo[ v ].length + g->edge_info2( e ) - more_ninfo[ u ].length;
            if ( delta2 > 0 ) {
                float delta1 = more_ninfo[ u ].dist + g->edge_info( e ) - more_ninfo[ v ].dist;
                max2( theta, ( lambda * g->edge_info2( e ) - delta1 ) / delta2 );
            }
#endif
        }  // for each edge

        // STEP: Using theta, update lambda as well as the distance of
        // every node:
        lambda -= theta;
        for ( int v = 0; v < n; ++v )
            more_ninfo[ v ].dist -= theta * more_ninfo[ v ].length;

    }  // main while loop

    delete [] more_ninfo;
    delete [] critical;

    return lambda;
} // find_min_cycle_ratio_for_scc

// End of file
