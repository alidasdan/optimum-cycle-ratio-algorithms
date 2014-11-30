//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
//
// An implementation of Lawler's minimum cycle ratio algorithm using
// Tarjan's shortest path algorithm. A. Goldberg's paper on Negative
// Cycle Detection describes Tarjan's algorithm.
//
//@techreport{Ta81, 
//author = {R.~E. Tarjan},
//title = {Shortest Paths},
//institution = {AT\&T Bell Laboratories, Murray Hill, NJ},
//year = 1981,
//type = {Tech. Rep.},
//address = "",
//number = {},
//month = {},
//}
//

// Parameters: NUPDATES, NCYCLES, CYCLELEN.

// count[0] = number of iterations
// count[1] = number of nodes in in Q
// count[2] = number of arcs visited (out of nodes in Q)
// count[3] = number of times T(v) is accessed
// count[4] = number of nodes in T(v)
// until a cycle is found

#include "ad_graph.h"
#include "ad_cqueue.h"

// Status of a node.
enum STATUS { 
    OUT_OF_Q = 0,  // Out of node queue.
    INACTIVE = 1,  // In the successor list but also in a deleted subtree.
    IN_Q = 2,      // Inside node queue.
    ACTIVE = IN_Q  // In the successor list and not in a deleted subtree.
};

// More node info for Tarjan's algorithm. parent is redundant since
// edge2parent's source node is parent.
struct ninfo_tarjan {
    float  dist;         // the min distance from the source.
    int    degree;       // the real degree in the tree - 1.
    int    parent;       // the parent node in the tree.
    int    edge2parent;  // the edge to the parent node.
    int    prev;         // the previous node in the successor list.
    int    next;         // the next node in the successor list.
    STATUS status;       // status as in STATUS.
};

// Note: The successor list in a tree corresponds to a preorder
// traversal of its nodes.

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

    ninfo_tarjan      *more_ninfo = new ninfo_tarjan[ n ];
    ad_cqueue< int >  nodeq( n );

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

        // STEP: Run Tarjan's algorithm to see if the graph has a negative
        // cycle.

#ifdef PROGRESS
        int QUEUELEN = 0;
        int NUPDATES = 0;
        int NEDGES = 0;
        int CYCLELEN = 0;
#endif
    
        // neg_cycle_found is set when a negative cycle is found.
        bool neg_cycle_found;
        {
            // STEP2: Initialize the info of each node.
            {
                ninfo_tarjan *ptr = &more_ninfo[ SOURCE ];

                ptr->dist = 0.0;
                ptr->degree = -1;
                ptr->prev = SOURCE;
                ptr->next = SOURCE;
                ptr->parent = SOURCE;
                ptr->edge2parent = -1;
                ptr->status = IN_Q; 

                for ( int v = 1; v < n; ++v ) {
                    ++ptr;

                    ptr->dist = f_plus_infinity;
                    ptr->degree = -1;
                    ptr->prev = -1;
                    ptr->next = -1;
                    ptr->parent = -1;
                    ptr->edge2parent = -1;
                    ptr->status = OUT_OF_Q;
                }
            }

            neg_cycle_found = false;  

            // STEP2: Enqueue the source. This explains why the status of
            // the source is IN_Q during the initialization above.
            nodeq.init();
            nodeq.put( SOURCE );

            while ( nodeq.is_not_empty() ) {

#ifdef REP_COUNT
                ++count[ 1 ];
#endif

#ifdef PROGRESS
                ++QUEUELEN;
#endif

                // STEP2: Dequeue u and determine its status.
                int u = nodeq.get();
                STATUS u_stat = more_ninfo[ u ].status;
                more_ninfo[ u ].status = OUT_OF_Q;
                if ( INACTIVE == u_stat )
                    continue;

                // STEP2: For each node v adjacent to node u, do:
#ifdef CYCLE_MEAN_VERSION
                float reduced_dist = more_ninfo[ u ].dist - lambda;
#else
                float udist = more_ninfo[ u ].dist;
#endif

                for ( int i = 0; i < g->outdegree( u ); ++i ) {
                    int v = g->ith_target_node( u, i );
                    int e_uv = g->ith_target_edge( u, i );

#ifdef PROGRESS
                    ++NEDGES;
#endif

#ifdef REP_COUNT
                    ++count[ 2 ];
#endif

                    // STEP2: If v's dist decreases, update it.
#ifdef CYCLE_MEAN_VERSION
                    float new_dist = reduced_dist + g->edge_info( e_uv );
#else
                    float new_dist = udist + g->edge_info( e_uv ) -
                        lambda * g->edge_info2( e_uv );
#endif

                    if ( new_dist < more_ninfo[ v ].dist ) {
                        more_ninfo[ v ].dist = new_dist;

#ifdef PROGRESS
                        ++NUPDATES;
#endif

                        if ( -1 != more_ninfo[ v ].prev ) {

#ifdef PROGRESS
                            int NTREENODES = 0;
#endif

#ifdef REP_COUNT
                            ++count[ 3 ];
#endif

                            // STEP2: Delete the subtree T(v) rooted at v. Also
                            // check if u is inside this subtree. If so, a negative
                            // cycle is found.
                            int before_v = more_ninfo[ v ].prev;
                            int total_degree = 0;
                            int w;
                            for ( w = v; total_degree >= 0; w = more_ninfo[ w ].next ) {

#ifdef REP_COUNT
                                ++count[ 4 ];
#endif

                                if ( w == u ) {
                                    // u is in the subtree rooted at v, so a negative
                                    // cycle is found.
                                    neg_cycle_found = true;
                                    {
#ifdef IMPROVE_UPPER_BOUND
                                        more_ninfo[ v ].parent = u;
                                        more_ninfo[ v ].edge2parent = e_uv;
                                        int x = u;
                                        int total_weight = 0;
                                        int total_length = 0;
#ifdef DEBUG
                                        int cycle_len = 0;
#endif
                                        do {
#ifdef DEBUG
                                            assert( cycle_len++ < n );
                                            assert( more_ninfo[ x ].parent >= 0 );
                                            assert( more_ninfo[ x ].edge2parent >= 0 );
#endif

#ifdef REP_COUNT
                                            ++count[ 5 ];
#endif

#ifdef CYCLE_MEAN_VERSION
                                            ++total_lengt;
#else
                                            total_length += g->edge_info2( more_ninfo[ x ].edge2parent );
#endif
                                            total_weight += g->edge_info( more_ninfo[ x ].edge2parent );
                                            x = more_ninfo[ x ].parent;
                                        } while ( x != u );

                                        float new_lambda = ( float ) total_weight / total_length;
                                        if ( new_lambda < lambda )
                                            lambda = new_lambda;

#ifdef PROGRESS
                                        CYCLELEN += total_length;
#endif

#endif
                                    }
                                    goto done;
                                }  // if w == u

                                total_degree += more_ninfo[ w ].degree;
                                more_ninfo[ w ].degree = -1;
                                more_ninfo[ w ].prev = -1;
                                if ( ACTIVE == more_ninfo[ w ].status )
                                    more_ninfo[ w ].status = INACTIVE;
#ifdef PROGRESS
                                ++NTREENODES;
#endif
                            }  // for w

                            more_ninfo[ more_ninfo[ v ].parent ].degree--;
                            more_ninfo[ before_v ].next = w;
                            more_ninfo[ w ].prev = before_v;

#ifdef PROGRESS
                            printf( "COUNTERS NTREENODES= %d\n", NTREENODES );
#endif

                        }  // if v has a prev

                        {
                            // STEP2: Make v a child of u.
                            more_ninfo[ v ].parent = u;
                            more_ninfo[ v ].edge2parent = e_uv;
                            more_ninfo[ u ].degree++;

                            // STEP2: Insert v into the successor list.
                            int after_u = more_ninfo[ u ].next;
                            more_ninfo[ u ].next = v;
                            more_ninfo[ v ].prev = u;
                            more_ninfo[ v ].next = after_u;
                            more_ninfo[ after_u ].prev = v;
                        }

                        // STEP2: Enqueue v if it is not already in the
                        // queue. Activate it if it was inactive. Note that ACTIVE
                        // and IN_Q are actually equal, so the status assignment
                        // can be moved out of these if-statements.
                        if ( OUT_OF_Q == more_ninfo[ v ].status ) {
                            more_ninfo[ v ].status = IN_Q;
                            nodeq.put( v );
                        } else {
                            more_ninfo[ v ].status = ACTIVE;
                        }
                    }  // if dist change occurs 
                }  // for each target node
            }  // while nodeq is not empty

        }  // End of Tarjan's algorithm

    done:

        // STEP: Update lambda depending on whether or not a negative
        // cycle is found.
        if ( neg_cycle_found ) {
            if ( ( upper - lambda ) < EPSILON2 )
                break;
            upper = lambda;
        }
        else {
            if ( ( lambda - lower ) < EPSILON2 )
                break;
            lower = lambda;
        }

    } // while 

    delete [] more_ninfo;

    return lambda;
}  // find_min_cycle_mean_for_scc

// End of file
