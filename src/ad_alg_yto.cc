//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
//
// An implementation of Young, Tarjan, and Orlin's minimum cycle ratio
// algorithm.
//
// @article{YoTaOr91,
// author = {N.~E. Young and R.~E. Tarjan and J.~B. Orlin},
// title = {Faster Parametric Shortest Path and Minimum-Balance Algorithms},
// journal = {Networks},
// year = 1991,
// volume = 21,
// number = "",
// pages = {205--221},
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
#error "ADD_SOURCE_NODE has to be defined for Young-Tarjan-Orlin's algorithm."
#endif

typedef ad_pq< float, int >::ad_pq_node_ptr pq_ptr;

// A heap node represents (ekey(v), key(v)) where v is a node, and
// key(v) is a predecessor edge of v that has the smallest edge key,
// which is ekey(v), and the edge key is defined in the reference
// above. edge_key array holds edge keys.

// More node info for Young-Tarjan-Orlin's algorithm.
struct ninfo_yto {
    int    dist;      // the weight of path from source.
    int    length;    // the length of path from source.
    int    key;       // node key (the edge w/ edge key = ekey)
    float  ekey;      // the inedge w/ the min key
    pq_ptr node2heap; // ptr to the heap node w/ key = ekey

    // fields for subtree management
    int    degree;    // the real degree in the tree - 1.
    int    parent;    // the parent node in the tree.
    int    prev;      // the previous node in the inorder tree.
    int    next;      // the next node in the inorder tree.
    bool   visited;   // set if this node is in the tree
};

/* ARGSUSED2 */
float
find_min_cycle_ratio_for_scc( const ad_graph< ninfo > *g, 
                              int plus_infinity,
                              float lambda_so_far )
{
    int n = g->num_nodes();
    int m = g->num_edges();

    ad_pq< float, int > pq( n + 1 );

    ninfo_yto *more_ninfo = new ninfo_yto[ n ];
    float      *edge_key = new float[ m ];

    // Plus and minus infinity: upper and lower bounds on lambda.
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
        ninfo_yto *ptr = &more_ninfo[ SOURCE ];

        // ptr->key is not used at all.
        // ptr->ekey is not used at all.
        // ptr->node2heap is not used at all.
        ptr->dist    = 0;      // SOURCE is 0 distance to itself.
        ptr->length  = 0;      
        ptr->degree  = n - 2;  // There is an edge to every other node in g.
        ptr->prev    = n - 1;  // For circularity.
        ptr->next    = 1;
        ptr->parent  = -1;
        ptr->visited = false;  

        for ( int v = 1; v < n; ++v ) {
            ++ptr;
            // ptr->key is set later.
            // ptr->node2heap is not used at all.
#ifdef DEBUG
            ptr->key     = -1;
#endif
            ptr->dist    = 0;
            ptr->length  = 1;
            ptr->ekey    = f_plus_infinity;
            ptr->degree  = -1;
            ptr->prev    = v - 1;
            ptr->next    = v + 1;
            ptr->parent  = SOURCE;
            ptr->visited = false;
        }
        ptr->next = SOURCE; // For circularity.
    }

    // STEP: Initialize the key of each node and edge, and insert the
    // node keys into the heap. The edge key of an edge is initially its
    // weight in g. The key of a node is the inedge with the min
    // key. The ekey of a node is the key of the node's key edge.
    {
        // Find the keys.
        ninfo_yto *ptr;
        for ( int e = 0; e < m; ++e ) {
#ifdef CYCLE_MEAN_VERSION
            edge_key[ e ] = ( float ) g->edge_info( e );
#else
            if ( 0 < g->edge_info2( e ) )
                edge_key[ e ] = ( float ) g->edge_info( e ) / g->edge_info2( e );
            else
                edge_key[ e ] = f_plus_infinity;
#endif
            ptr = &more_ninfo[ g->target( e ) ];
            if ( edge_key[ e ] <= ptr->ekey ) {
                // <= is needed because f_plus_infinity <= f_plus_infinity.
                ptr->key = e;
                ptr->ekey = edge_key[ e ];
            }
        }

        // Plus infinity in the heap is a special mark. Once found, it
        // indicates that g is acyclic.
        pq.put( f_plus_infinity, -1 );

        // Insert node keys into the heap.
        ptr = &more_ninfo[ 1 ];
        for ( int v = 1; v < n; ++v, ++ptr ) {
#ifdef DEBUG
            assert( -1 != ptr->key );
#endif
            ptr->node2heap = pq.put( ptr->ekey, ptr->key );
        }
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

#if DEBUG
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
            // STEP2: Find delta values for u->v.  
#ifdef CYCLE_MEAN_VERSION
            int delta1 = more_ninfo[ u ].dist + g->edge_info( e_min ) - more_ninfo[ v ].dist;
            int delta2 = more_ninfo[ u ].length + 1 - more_ninfo[ v ].length;
#else
            int delta1 = more_ninfo[ u ].dist + g->edge_info( e_min ) - more_ninfo[ v ].dist;
            int delta2 = more_ninfo[ u ].length + g->edge_info2( e_min ) - more_ninfo[ v ].length;
#endif

            // STEP2: Using the delta values, update the dist and length of
            // each node in T(v). Also check if u is in
            // T(v). If so, a cycle is found, so
            // exit. The last check is done here for performance. If the
            // node distances are also desired, doing the check here
            // corrupts them. We can do the check in a separate loop before
            // computing the delta values.
            int w;
            {
                ninfo_yto *ptr_w;
                w = v;
                for ( int total_degree = 0; total_degree >= 0; w = ptr_w->next ) {
                    if ( w == u )  // Check for a cycle.
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

        // STEP (part 1): Update node and edge keys. The edges whose keys
        // change are those of which exactly one end point is in T(v).
        {
            // STEP2: Go over every edge entering T(v).
            ninfo_yto *ptr_y;
            for ( int total_degree = 0, y = v; total_degree >= 0; y = ptr_y->next ) {
                ptr_y = &more_ninfo[ y ];
                total_degree += ptr_y->degree;

#ifdef REP_COUNT
                ++count[ 1 ];
#endif

                // y's key must be found among ALL of its inedges from
                // STRATCH. My original implementation was wrong; the current
                // approach below is not only correct but is cleaner and
                // simpler too.
                ptr_y->ekey = f_plus_infinity;
                for ( int i = 0; i < g->indegree( y ); ++i ) {
                    int e = g->ith_source_edge( y, i ); // e = x->y
                    ninfo_yto *ptr_x = &more_ninfo[ g->ith_source_node( y, i ) ];

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
                        if ( delta2 > 0 )
                            edge_key[ e ] = ( float ) ( ptr_x->dist + g->edge_info( e ) - ptr_y->dist ) / delta2;
                        else
                            edge_key[ e ] = f_plus_infinity;
#ifdef PROGRESS
                        ++NEDGECROSS;
#endif
                    }

                    if ( edge_key[ e ] <= ptr_y->ekey ) {
                        ptr_y->key = e;
                        ptr_y->ekey = edge_key[ e ];
                    }
                } // for each incoming edge

                pq.update_node( ptr_y->ekey, ptr_y->key, ptr_y->node2heap );
#ifdef PROGRESS
                ++NUPDATES;
#endif
            } // for each y in the subtree
        }

        // STEP (part 2): Update node and edge keys. The edges whose keys
        // change are those of which exactly one end point is in T(v).
        {
            // STEP2: Go over every edge leaving T(v).
            ninfo_yto *ptr_x;
            for ( int total_degree = 0, x = v; total_degree >= 0; x = ptr_x->next ) {
                ptr_x = &more_ninfo[ x ];
                total_degree += ptr_x->degree;

                // If e's key changes, it cannot increase. Thus, y's key
                // should be found among only the edges whose keys can
                // change. These are exactly the leaving edges.
                ninfo_yto *ptr_y;
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
                        if ( delta2 > 0 )
                            edge_key[ e ] = ( float ) ( ptr_x->dist + g->edge_info( e ) - ptr_y->dist ) / delta2;
                        else
                            edge_key[ e ] = f_plus_infinity;

#ifdef PROGRESS
                        ++NEDGECROSS;
#endif

                        // Even if y's key is e, the assignments below don't hurt.
                        if ( edge_key[ e ] < ptr_y->ekey ) {
                            ptr_y->key = e;
                            ptr_y->ekey = edge_key[ e ];
                            pq.update_node( ptr_y->ekey, ptr_y->key, ptr_y->node2heap );
#ifdef PROGRESS
                            ++NUPDATES;
#endif
                        }
                    }
                } // for each incoming edge

                // Reset visited(x) for every x in T(v) here for the next tree.
                ptr_x->visited = false;
            } // for each x in the subtree
        }

#ifdef PROGRESS
        printf( "COUNTERS NTREENODES=%d NEDGECROSS= %d NUPDATES= %d\n", NTREENODES, NEDGECROSS, NUPDATES );
#endif
    }  // main while loop

 done:

    delete [] more_ninfo;
    delete [] edge_key;

    return lambda;
}  // find_min_cycle_ratio_for_scc

// End of file
