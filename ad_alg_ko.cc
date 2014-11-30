//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
//
// An implementation of Karp and Orlin's minimum cycle ratio algorithm.
//
// @article{KaOr81,
// author = {R.~M. Karp and J.~B. Orlin},
// title = {Parametric Shortest Path Algorithms with An Application 
//          to Cyclic Staffing},
// journal = {Discrete Applied Mathematics},
// year = 1981,
// volume = 3,
// number = "",
// pages = {37--45},
// month = "",
// }
//

// Parameters of interest: NTREENODES, NEDGECROSS.

// count[0] = number of iterations
// count[1] = number of nodes in T(v)
// count[2] = number of arcs entering or leaving T(v)

#include "ad_graph.h"
#include "ad_pq.h"

#ifndef ADD_SOURCE_NODE
#error "ADD_SOURCE_NODE has to be defined for Karp-Orlin's algorithm."
#endif

typedef ad_pq< float, int >::ad_pq_node_ptr pq_ptr;

// A heap node represents (edge_key(e), e) where e is an edge, and
// edge_key is defined in the reference above.

// More node info for Karp and Orlin's algorithm.
struct ninfo_ko {
    int   dist;   // the weight of path from source.
    int   length; // the length of path from source.

    // fields for subtree management
    int   degree; // the real degree in the tree - 1.
    int   parent; // the parent node in the tree.
    int   prev;   // the previous node in the inorder tree.
    int   next;   // the next node in the inorder tree.
    bool  visited; 
};

/* ARGSUSED2 */
float
find_min_cycle_ratio_for_scc( const ad_graph< ninfo > *g, 
                              int plus_infinity,
                              float lambda_so_far )
{
    int n = g->num_nodes();
    int m = g->num_edges();

    ad_pq< float, int > pq( m + 1 );

    ninfo_ko  *more_ninfo = new ninfo_ko[ n ];
    pq_ptr    *edge2heap = new pq_ptr[ m ];

    float f_plus_infinity = ( float ) plus_infinity;
    float f_minus_infinity = -f_plus_infinity;

    // STEP: Set the initial values of the fields of more_ninfo
    // array. Also construct the initial subtree rooted at the source
    // node SOURCE. Denote this tree as T(0). This tree has imaginary
    // edges to every node in the graph. These edges are not included in
    // the graph g, so the number m of edges does not include them. The
    // initial subtree is stored in more_ninfo in preorder. The preorder
    // is the circular list of the node numbers in ascending order: 0,
    // 1, 2, ..., (n-1). For a node in this list, prev and next point to
    // the previous and next nodes in the list.
    {
        ninfo_ko *ptr = &more_ninfo[ SOURCE ];

        ptr->dist    = 0;
        ptr->length  = 0;
        ptr->degree  = n - 2;  // There is an edge to every other node in g.
        ptr->prev    = n - 1;  // For circularity.
        ptr->next    = 1;
        ptr->parent  = -1;
        ptr->visited = false;

        for ( int v = 1; v < n; ++v ) {
            ++ptr;

            ptr->dist    = 0;
            ptr->length  = 1;
            ptr->degree  = -1;
            ptr->prev    = v - 1;
            ptr->next    = v + 1;
            ptr->parent  = SOURCE;
            ptr->visited = false;
        }
        ptr->next = SOURCE;  // For circularity.
    }

    // STEP: Initialize each edge key and insert them into the heap.
    {
        // Plus infinity in the heap is a special mark. Once found, it
        // indicates that g is acyclic.
        pq.put( f_plus_infinity, -1 );

        // Insert edge keys into the heap.
#ifdef CYCLE_MEAN_VERSION
        for ( int e = 0; e < m; ++e )
            edge2heap[ e ] = pq.put( ( float ) g->edge_info( e ), e );
#else
        for ( int e = 0; e < m; ++e ) {
            if ( 0 < g->edge_info2( e ) )
                edge2heap[ e ] = pq.put( ( float ) g->edge_info( e ) / g->edge_info2( e ), e );
            else
                edge2heap[ e ] = pq.put( f_plus_infinity, e );
        }
#endif
    }

    // STEP: Iterate until a cycle is found or +infinity in the heap is
    // seen.
    float lambda;
#ifdef DEBUG
    float prev_lambda = f_minus_infinity;
#endif
    while ( true ) {

        // STEP: Find the edge e_min whose key is the min key in the
        // heap. lambda is equal to its key.
        int e_min = pq.getinfo();
        lambda = pq.getkey();
        if ( lambda >= f_plus_infinity )
            goto done;

#ifdef REP_COUNT
        count[ 0 ]++;
#endif
#ifdef REP_COUNT_PRINT
        printf( "REP_COUNT Iteration number= %d lambda= %10.2f\n", count[ 0 ], lambda );
#endif

#ifdef DEBUG
        if ( lambda >= prev_lambda ) {
            prev_lambda = lambda;
        } else {
            printf( "ERROR: lambda must be nondecreasing.\n" );
            abort();
        }
#endif

#ifdef PROGRESS
        int NTREENODES = 0;
        int NEDGECROSS = 0;
        int NUPDATES = 0;
#endif

        // e_min = u -> v.
        int u = g->source( e_min );
        int v = g->target( e_min );

        // STEP: Find delta values for u and v of e_min. Update the dist
        // and length field of each node in the subtree rooted at
        // v. Denote this tree as T(v). Also construct the new tree if v
        // is not a parent of u.
        {
            // STEP: Find delta values for u->v.
#ifdef CYCLE_MEAN_VERSION
            int delta1 = more_ninfo[ u ].dist + g->edge_info( e_min ) - more_ninfo[ v ].dist;
            int delta2 = more_ninfo[ u ].length + 1 - more_ninfo[ v ].length;
#else
            int delta1 = more_ninfo[ u ].dist + g->edge_info( e_min ) - more_ninfo[ v ].dist;
            int delta2 = more_ninfo[ u ].length + g->edge_info2( e_min ) - more_ninfo[ v ].length;
#endif

            int w;
            {
                ninfo_ko *ptr_w;

                w = v;
                for ( int total_degree = 0; total_degree >= 0; w = ptr_w->next ) {
                    if ( w == u ) // Check for a cycle.
                        goto done;
                    ptr_w = &more_ninfo[ w ];
                    total_degree += ptr_w->degree;
                    ptr_w->dist += delta1;
                    ptr_w->length += delta2;
                    ptr_w->visited = true;
#ifdef PROGRESS
                    ++NTREENODES;
#endif
                }
            }

            // STEP2: Find the new tree. To get it, delete the unique edge
            // directed into v in the shortest paths tree T(0), and add the
            // edge e_min.  Note that node w comes from the above for loop,
            // and it is the next node of the last node in T(v).

            /* Here is the picture:
             *    v's (old)parent  u (v's new parent)
             *          \__      __/ 
             *             \    /
             *  before_v --> v       (Tv = T(v))
             *              / \
             *             /   \
             *            /____last_of_Tv ---> w
             * 
             * Before picture as list: 
             *    (v's old parent ... before_v (v ... last_of_Tv) w ...)
             * After pictures as lists: 
             *    (v's old parent ... before_v w ...)
             *    (u (v ... last_of_Tv) after_u)
             */

            // Remove T(v) between before_v and w.
            more_ninfo[ more_ninfo[ v ].parent ].degree--;
            int before_v = more_ninfo[ v ].prev;
            int last_of_Tv = more_ninfo[ w ].prev;
            more_ninfo[ before_v ].next = w;
            more_ninfo[ w ].prev = before_v;

            // Insert T(v) between u and its next node. Since u can be equal
            // to before_v, after_u must be assigned after the above two
            // lines.
            more_ninfo[ u ].degree++;
            more_ninfo[ v ].parent = u;
            int after_u = more_ninfo[ u ].next;
            more_ninfo[ u ].next = v;
            more_ninfo[ v ].prev = u;
            more_ninfo[ last_of_Tv ].next = after_u;
            more_ninfo[ after_u ].prev = last_of_Tv ;
        }

        // STEP (part 1): Update the edge keys. The edges whose keys
        // change are those of which exactly one end point is in T(v).
        {
            // STEP2: Go over every edge entering T(v).
            ninfo_ko *ptr_y;
            for ( int total_degree = 0, y = v; total_degree >= 0; y = ptr_y->next ) {
                ptr_y = &more_ninfo[ y ];
                total_degree += ptr_y->degree;

#ifdef REP_COUNT
                ++count[ 1 ];
#endif

                for ( int i = 0; i < g->indegree( y ); ++i ) {
                    int e = g->ith_source_edge( y, i ); // e = x->y
                    ninfo_ko *ptr_x = &more_ninfo[ g->ith_source_node( y, i ) ];

#ifdef REP_COUNT
                    ++count[ 2 ];
#endif
                    // For an entering edge e=x->y where y in T(v), x is not in
                    // T(v), i.e., x is not visited but y is visited. Only the
                    // key of an entering edge can change.
                    if ( ptr_x->visited != ptr_y->visited ) {
#ifdef CYCLE_MEAN_VERSION
                        int delta2 = ptr_x->length + 1 - ptr_y->length;
#else            
                        int delta2 = ptr_x->length + g->edge_info2( e ) - ptr_y->length;
#endif
                        if ( delta2 > 0 ) {
                            int delta1 = ptr_x->dist + g->edge_info( e ) - ptr_y->dist;
                            pq.update_key( ( float ) delta1 / delta2, edge2heap[ e ] );
                        } else {
                            pq.update_key( f_plus_infinity, edge2heap[ e ] );
                        }
#ifdef PROGRESS
                        ++NEDGECROSS;
                        ++NUPDATES;
#endif
                    }
                } // for each incoming edge
            } // for each y in the subtree
        }

        // STEP (part 2): Update the edge keys. The edges whose keys
        // change are those of which exactly one end point is in T(v).
        {
            // STEP2: Go over every edge leaving T(v).
            ninfo_ko *ptr_x;
            for ( int total_degree = 0, x = v; total_degree >= 0; x = ptr_x->next ) {
                ptr_x = &more_ninfo[ x ];
                total_degree += ptr_x->degree;

                // If e's key changes, it cannot increase. Thus, y's key
                // should be found among only the edges whose keys can
                // change. These are exactly the leaving edges.
                ninfo_ko *ptr_y;
                for ( int i = 0; i < g->outdegree( x ); ++i ) {
                    int e = g->ith_target_edge( x, i ); // e = x->y
                    ptr_y = &more_ninfo[ g->ith_target_node( x, i ) ];

#ifdef REP_COUNT
                    ++count[ 2 ];
#endif
                    // For a leaving edge e=x->y where x in T(v), y is not in
                    // T(v), i.e., x is visited but y is not visited. Only the
                    // key of a leaving edge can change.
                    if ( ptr_x->visited != ptr_y->visited ) {
#ifdef CYCLE_MEAN_VERSION
                        int delta2 = ptr_x->length + 1 - ptr_y->length;
#else            
                        int delta2 = ptr_x->length + g->edge_info2( e ) - ptr_y->length;
#endif
                        if ( delta2 > 0 ) {
                            int delta1 = ptr_x->dist + g->edge_info( e ) - ptr_y->dist;
                            pq.update_key( ( float ) delta1 / delta2, edge2heap[ e ] );
                        } else {
                            pq.update_key( f_plus_infinity, edge2heap[ e ] );
                        }
#ifdef PROGRESS
                        ++NEDGECROSS;
                        ++NUPDATES;
#endif
                    }
                } // for each incoming edge

                // Reset visited(x) for every x in T(v) here for the next tree.
                ptr_x->visited = false;
            } // for each x in the subtree
        }

#ifdef PROGRESS
        printf( "COUNTERS NTREENODES= %d NEDGECROSS= %d NUPDATES= %d\n", NTREENODES, NEDGECROSS, NUPDATES );
#endif

    }  // main while loop

 done:

    delete [] more_ninfo;
    delete [] edge2heap;

    return lambda;
}  // find_min_cycle_ratio_for_scc

// End of file
