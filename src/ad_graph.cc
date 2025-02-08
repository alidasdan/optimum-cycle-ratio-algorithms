//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
//#include <new.h>
#include <new>
#include "ad_graph.h"
#include "ad_util.h"

///////////////////////////////////////////////////////////////////////
// FUNCTION PROTOTYPES:
///////////////////////////////////////////////////////////////////////

// Find the number of SCCs as well as a map comp_num of which node
// belongs to which SCC.
int 
num_sccs( const ad_graph< ninfo >& g, int *comp_num );

// Traverse g in DFS order using outedges.
void 
traverse_via_outedges( const ad_graph< ninfo >& g, 
                       int v, bool *not_visited, 
                       int& n, int *node_list );

// Traverse g in DFS order using outedges. Iterative DFS for scalability.
void 
traverse_via_outedges_iter( const ad_graph< ninfo >& g, 
                            int v, bool *not_visited, 
                            int& n, int *node_list );

// Traverse g in DFS order using outedges. This recursive version runs
// into stack size limitations for some large random graphs in the
// graph-benchmarks repo so the non-recursive version above.
void 
traverse_via_outedges_recur( const ad_graph< ninfo >& g, 
                             int v, bool *not_visited, 
                             int& n, int *node_list );

// Traverse g in DFS order using inedges.
void 
traverse_via_inedges( const ad_graph< ninfo >& g, 
                      int v, int num_comps, int *comp_num );


// Traverse g in DFS order using inedges. Iterative DFS for scalability.
void 
traverse_via_inedges_iter( const ad_graph< ninfo >& g, 
                           int v, int num_comps, int *comp_num );

// Traverse g in DFS order using inedges. This recursive version runs
// into stack size limitations for some large random graphs in the
// graph-benchmarks repo so the non-recursive version above.
void 
traverse_via_inedges_recur( const ad_graph< ninfo >& g, 
                            int v, int num_comps, int *comp_num );

// Find the properties of SCCs before their dynamic creation.
extern 
void
find_component_props( ad_graph< cninfo >& cg, 
                      const ad_graph< ninfo >& g, int *comp_num );

// Check for errors after fscanf.
void 
check_for_error( int status, int val );

// Remove comments and other junk from input files in DIMACS format.
char
remove_junk( FILE *fp );

///////////////////////////////////////////////////////////////////////
// External I/O functions:
void 
print( int x )
{
    printf( "%d", x );
}

void 
print( const ninfo& ni ) 
{
    ni.print();
}

void 
print( const cninfo& cni )
{
    cni.print();
}

///////////////////////////////////////////////////////////////////////
// Build the adjacency lists of the graph.
// pre: The lists are already allocated.
// pre: The nodes and edges are already inserted into the lists.
template< class ninfo_t >
void 
ad_graph< ninfo_t >::
build_adj()
{
#ifdef DEBUG
    assert( nlist != NULL );
    assert( elist != NULL );
#endif

    // Correct num_nodes and num_edges, and prevent further updates.
    update_nums();

    // Set the list heads based on the node degrees.
    int pre_in_sum = 0;
    int pre_out_sum = 0;
    for ( int v = 0; v < num_nodes(); ++v ) {
        first_inedge( v, pre_in_sum );
        first_outedge( v, pre_out_sum );
        pre_in_sum += indegree( v );
        pre_out_sum += outdegree( v );
    }

    int *in_inx = new int[ num_nodes() ];
    int *out_inx = new int[ num_nodes() ];

    // Initialize the indices into the lists.
    for( int v = 0; v < num_nodes(); ++v ) {
        in_inx[ v ] = 0;
        out_inx[ v ] = 0;
    }

    // Insert edge e = s->t into the adjacency lists of nodes s and t.
    for ( int e = 0; e < num_edges(); ++e )  {
        int s = source( e );
        int t = target( e );
        int ei = edge_info( e );
        int ei2 = edge_info2( e );

#ifdef CYCLE_MEAN_VERSION
        outlist[ nlist[ s ].first_outedge() + out_inx[ s ] ].set_edge( e, t, ei );
        inlist[ nlist[ t ].first_inedge() + in_inx[ t ] ].set_edge( e, s, ei );

#else
        outlist[ nlist[ s ].first_outedge() + out_inx[ s ] ].set_edge( e, t, ei, ei2 );
        inlist[ nlist[ t ].first_inedge() + in_inx[ t ] ].set_edge( e, s, ei, ei2 );
#endif

        out_inx[ s ]++;
        in_inx[ t ]++;
    }  

    delete [] out_inx;
    delete [] in_inx;
}  // build_adj

///////////////////////////////////////////////////////////////////////
// Check for errors after fscanf.
void 
check_for_error( int status, int val )
{
    if ( EOF == status ) {
        printf( "ERROR: EOF is reached before all edges are read.\n" );
        abort();
    } else if ( status != val ) {
        printf( "ERROR: Cannot read %d items.\n", val );
        abort();
    }
}  // check_for_error

// Remove comments and other junk from input files in DIMACS format.
char
remove_junk( FILE *fp ) 
{
    char buf[ MAX_LINE_SIZE ];
    int status;

    status = fscanf( fp, "%s", buf );
#ifdef DEBUG
    check_for_error( status, 1 );
#endif
    while ( ( 'c' == buf[ 0 ] ) || ( 't' == buf[ 0 ] ) || ( 'n' == buf[ 0 ] ) ) {
#ifdef DEBUG
        char *ptr = fgets(buf, sizeof( buf ), fp);
        if ( NULL == ptr ) {
            printf( "ERROR: EOF is reached before all edges are read.\n" );
            abort();
        }
#else
        ( void ) fgets(buf, sizeof( buf ), fp);
#endif
        status = fscanf( fp, "%s", buf );
#ifdef DEBUG
        check_for_error( status, 1 );
#endif
    }  // while

    return buf[ 0 ];
}

// Read the graph from file input_file. Duplicate edges are removed.
template< class ninfo_t >
void
ad_graph< ninfo_t >::
read( ginfo& gi, const args_t& args )
{
    FILE *fp;

    if ( ( fp = fopen( args.input_file, "r" ) ) == NULL ) {
        printf( "ERROR: Cannot open the input file %s.\n", args.input_file );
        abort();
    }

    char buf[ MAX_LINE_SIZE ];
#ifdef CYCLE_MEAN_VERSION
    int  status, total_weight = 0; // of the edges
#else
    int  status, total_weight = 0, total_ttime = 0; // of the edges
#endif
    bool has_self_loop = false;

    buf[ 0 ] = remove_junk( fp );

    // Read the problem line, the line starting with 'p'.
    if ( 'p' != buf[ 0 ] ) {
        printf( "ERROR: Input file %s is not in DIMACS format.\n", args.input_file );
        abort();
    }

    // Read the rest of the line after the line descriptor. This line
    // assumes that the problem name will not overflow buf.
    status = fscanf( fp, "%s%d%d", buf, &nnodes, &nedges );
#ifdef DEBUG
    check_for_error( status, 3 );
#endif

    if ( ( nnodes <= 0 ) || ( nedges < 0 ) ) {
        printf( "ERROR: Require 'nnodes > 0' and 'nedges >= 0'. \n" );
        abort();
    }
  
    if ( nedges ) {
        buf[ 0 ] = remove_junk( fp );

        if ( 'a' != buf[ 0 ] ) {
            printf( "ERROR: Input file %s is not in DIMACS format.\n", args.input_file );
            abort();
        }
    }

    // Create node, edge, and adj lists.
    create( nnodes, nedges );

    // Edge (u, v) with weight w and transit time t. The transit time is
    // ignored for the cycle mean problems.
    int u, v, w, t; 

    if ( nedges ) {
        // Read the rest of the line for edge if buf[0] is 'a'.
        status = fscanf( fp, "%d%d%d%d", &u, &v, &w, &t );
#ifdef DEBUG
        check_for_error( status, 4 );
#endif
    }

    // Read the edges. If the file doesn't contain nedges lines, we will
    // have an error after the following fscanf.
    for ( int e = 0; e < nedges; ++e ) {
#ifdef DEBUG
        if ( ( u < 1 || u > nnodes ) || ( v < 1 || v > nnodes ) )  {
            printf( "ERROR: Invalid node number.\n" );
            abort();
        }
#endif

        if ( 0 == args.mode ) {
            w -= args.offset;
        } else {
            w = ( *dist_func )( args.w1, args.w2 ) - args.offset;
        }

        total_weight += abs_val( w );
#ifdef DEBUG
        assert( 0 <= total_weight );
#endif

#ifndef CYCLE_MEAN_VERSION
        total_ttime += t;
#ifdef DEBUG
        assert( 0 <= total_ttime );

        if ( t < 1 ) {
            printf(" ERROR: Transit time must be a positive integer.\n");
            abort();
        }
#endif
#endif

        --u;
        --v;

        if ( u == v )
            has_self_loop = true;

#ifdef CYCLE_MEAN_VERSION
        if ( args.min_version )
            ins_edge( u, v, w );
        else 
            ins_edge( u, v, -w );

#else
        if ( args.min_version )
            ins_edge( u, v, w, t );
        else 
            ins_edge( u, v, -w, t );
#endif


        if ( e == nedges - 1 )
            break;

        status = fscanf( fp, "%s%d%d%d%d", buf, &u, &v, &w, &t );
#ifdef DEBUG
        check_for_error( status, 5 );
#endif

        if ( 'a' != buf[ 0 ] ) {
#ifdef DEBUG
            printf( "ERROR: An edge is expected.\n" );
            abort();
#endif
        }
    }  // for each e, read e.

    fclose( fp );

    build_adj();

#ifdef CYCLE_MEAN_VERSION
    // total_weight = 2 + |w( e )| for all e.
    gi.total_edge_weight = 2 + total_weight;
    gi.has_self_loop = has_self_loop;

#else
    // total_weight = 2 + |w( e )| for all e.
    gi.total_edge_weight = 2 + total_weight;
    gi.total_trans_time = total_ttime;
    gi.has_self_loop = has_self_loop;
#endif

#ifdef DEBUG
    assert( total_weight < gi.total_edge_weight );
    assert( 0 <= gi.total_edge_weight );
#endif
}  // read

// Generate the given graph's weights with a given distribution.
template< class ninfo_t >
void    
ad_graph< ninfo_t >::
generate_part( ginfo& gi, const args_t& args )
{
    for ( int e = 0; e < nedges; ++e ) {
        int w = ( *dist_func )( args.w1, args.w2 ) - args.offset;
        int t = ( *dist_func )( args.t1, args.t2 );

#ifdef CYCLE_MEAN_VERSION
        if ( args.min_version )
            edge_info( e, w );
        else 
            edge_info( e, -w );

#else
        if ( args.min_version )
            edge_info( e, w, t );
        else 
            edge_info( e, -w, t );
#endif
    }
}  // generate_part

// Generate an entire graph with a given distribution.
template< class ninfo_t >
void    
ad_graph< ninfo_t >::
generate_all( ginfo& gi, const args_t& args )
{
    struct tnode_t {
        int out_deg, out_inx;
    } *tnodes;

    struct tedge_t {
        int src, tar;
        int next_inx;
    } *tedges;

#ifdef CYCLE_MEAN_VERSION
    int  total_weight = 0;  // of the edges.
#else
    int  total_weight = 0, total_ttime = 0;  // of the edges.
#endif
    bool has_self_loop = false;
  
    nnodes = args.nnodes;
    nedges = args.nedges;

    if ( ( nnodes <= 0 ) || ( nedges < 0 ) ) {
        printf( "ERROR: Require 'nnodes > 0' and 'nedges >= 0'. \n" );
        abort();
    }
  
    // Create node, edge, and adj lists.
    create( nnodes, nedges );

    tnodes = new tnode_t[ nnodes ];
    tedges = new tedge_t[ nedges ];

    for ( int v = 0; v < nnodes; ++v ) {
        tnodes[ v ].out_deg = 0;
        tnodes[ v ].out_inx = -1;
    }

    for ( int e = 0; e < nedges; ++e ) {
        tedges[ e ].src = tedges[ e ].tar = -1;
        tedges[ e ].next_inx = -1;
    }

    int e = 0;

    {
        // Add a cycle around all the nodes.

        for ( int u = 0; u < ( nnodes - 1 ); ++u ) {
            int v = u + 1;

#ifdef CYCLE_MEAN_VERSION
            int w = ( *dist_func )( args.w1, args.w2 ) - args.offset;

            total_weight += abs_val( w );
#ifdef DEBUG
            assert( 0 <= total_weight );
#endif

            if ( args.min_version )
                ins_edge( u, v, w );
            else 
                ins_edge( u, v, -w );

#else
            int w = ( *dist_func )( args.w1, args.w2 ) - args.offset;
            int t = ( *dist_func )( args.t1, args.t2 );

            total_weight += abs_val( w );
            total_ttime += t;
#ifdef DEBUG
            assert( 0 <= total_weight );
            assert( 0 <= total_ttime );
#endif

            if ( args.min_version )
                ins_edge( u, v, w, t );
            else 
                ins_edge( u, v, -w, t );
#endif

            tnodes[ u ].out_deg++;
            tedges[ e ].src = u;
            tedges[ e ].tar = v;
            tedges[ e ].next_inx = tnodes[ u ].out_inx;
            tnodes[ u ].out_inx = e;
            ++e;
        }

#ifdef CYCLE_MEAN_VERSION
        int u = nnodes - 1;
        int v = 0;
        int w = ( *dist_func )( args.w1, args.w2 ) - args.offset;
        total_weight += abs_val( w );
#ifdef DEBUG
        assert( 0 <= total_weight );
#endif

        if ( args.min_version )
            ins_edge( u, v, w );
        else 
            ins_edge( u, v, -w );

#else
        int u = nnodes - 1;
        int v = 0;
        int w = ( *dist_func )( args.w1, args.w2 ) - args.offset;
        int t = ( *dist_func )( args.t1, args.t2 );

        total_weight += abs_val( w );
        total_ttime += t;
#ifdef DEBUG
        assert( 0 <= total_weight );
        assert( 0 <= total_ttime );
#endif

        if ( args.min_version )
            ins_edge( u, v, w, t );
        else 
            ins_edge( u, v, -w, t );
#endif
    
        tnodes[ u ].out_deg++;
        tedges[ e ].src = u;
        tedges[ e ].tar = v;
        tedges[ e ].next_inx = tnodes[ u ].out_inx;
        tnodes[ u ].out_inx = e;
        ++e;
    }

    for ( ; e < nedges; ++e ) {
        int u, v;

        int curr_e;
        do {
            do {
                u = uniform_dist( 0, nnodes - 1 );
                v = uniform_dist( 0, nnodes - 1 );
            } while ( u == v );

            curr_e = tnodes[ u ].out_inx;
            while ( -1 != curr_e ) {
                if ( v == tedges[ curr_e ].tar )
                    break;
                curr_e = tedges[ curr_e ].next_inx;
            }
        } while ( -1 != curr_e );

        tedges[ e ].src = u;
        tedges[ e ].tar = v;
        tedges[ e ].next_inx = tnodes[ u ].out_inx;
        tnodes[ u ].out_inx = e;

#ifdef CYCLE_MEAN_VERSION
        int w = ( *dist_func )( args.w1, args.w2 ) - args.offset;

        total_weight += abs_val( w );
#ifdef DEBUG
        assert( 0 <= total_weight );
#endif

        if ( args.min_version )
            ins_edge( u, v, w );
        else 
            ins_edge( u, v, -w );

#else
        int w = ( *dist_func )( args.w1, args.w2 ) - args.offset;
        int t = ( *dist_func )( args.t1, args.t2 );

        total_weight += abs_val( w );
        total_ttime += t;
#ifdef DEBUG
        assert( 0 <= total_weight );
        assert( 0 <= total_ttime );
#endif

        if ( args.min_version )
            ins_edge( u, v, w, t );
        else 
            ins_edge( u, v, -w, t );
#endif
    }

    delete [] tnodes;
    delete [] tedges;

    build_adj();

#ifdef CYCLE_MEAN_VERSION
    // total_weight = 2 + |w( e )| for all e.
    gi.total_edge_weight = 2 + total_weight;
    gi.has_self_loop = has_self_loop;

#else
    // total_weight = 2 + |w( e )| for all e.
    gi.total_edge_weight = 2 + total_weight;
    gi.total_trans_time = total_ttime;
    gi.has_self_loop = has_self_loop;
#endif

#ifdef DEBUG
    assert( total_weight < gi.total_edge_weight );
    assert( 0 <= gi.total_edge_weight );
#endif
}  // generate_all

///////////////////////////////////////////////////////////////////////

// Print the graph.
template< class ninfo_t >
void
ad_graph< ninfo_t >::
print( bool all_out ) const
{
    printf( "n = %d m = %d\n", num_nodes(), num_edges() );

    printf( "Edges:\n" );
    for ( int e = 0; e < num_edges(); ++e ) {
        printf( "EDGE %d : %d -> %d : w = ", e, source( e ), target( e )  );
        elist[ e ].print_info();
        printf( "\n" );
    }

    if ( !all_out )
        return;

    printf( "Out Adj:\n" );
    for ( int v = 0; v < num_nodes(); ++v ) {
        printf( "NODE %d WITH %d OUT-EDGES( s ) : ", v, outdegree( v ) );
        for ( int i = 0; i < outdegree( v ); ++i )
            printf( " e%d:%d", ith_target_edge( v, i ), ith_target_node( v, i ) );
        printf( "\n" );
    }

    printf( "In Adj:\n" );
    for ( int v = 0; v < num_nodes(); ++v ) {
        printf( "NODE %d WITH %d IN-EDGES( s ) : ", v, indegree( v ) );
        for ( int i = 0; i < indegree( v ); ++i )
            printf( " e%d:%d", ith_source_edge( v, i ), ith_source_node( v, i ) );
        printf( "\n" );
    }
}  // print

// Print the graph to a file.
template< class ninfo_t >
void
ad_graph< ninfo_t >::
fprint( const args_t& args ) const
{
    if ( !strcmp( args.dump_file, "" ) )
        return;

    FILE *fp;

    if ( ( fp = fopen( args.dump_file, "w") ) == NULL ) {
        printf( "ERROR: Cannot open the output file. Abort. \n" );
        abort();
    }

    fprintf( fp, "p generated-%ld %d %d\n", args.seed, num_nodes(), num_edges() );
    for ( int e = 0; e < num_edges(); ++e ) {
#ifdef CYCLE_MEAN_VERSION
        fprintf( fp, "a %d %d %d 1\n", 1 + source( e ), 1 + target( e ), edge_info( e ) );

#else
        fprintf( fp, "a %d %d %d %d\n", 1 + source( e ), 1 + target( e ), 
                 edge_info( e ), edge_info2( e ) );
#endif
    }

    fclose( fp );
}  // fprint

///////////////////////////////////////////////////////////////////////

// Remove duplicate edges and print graph.
template< class ninfo_t >
void
ad_graph< ninfo_t >::
print_no_duplicates( bool min_version )
{
    int n = num_nodes();
    int m = num_edges();
    int m_orig = m;

    bool *marked = new bool[ m ];
  
    for ( int e = 0; e < m; ++e )
        marked[ e ] = true;
  
    for ( int v = 0; v < n; ++v ) {
        for ( int i = 0; i < outdegree( v ); ++i ) {
            // e_i = (v, i, t, w_i)
            int e_i = ith_target_edge( v, i );  
            int t = target( e_i );
            int w_i = edge_info( e_i );

            for ( int j = i + 1; j < outdegree( v ); ++j ) {
                if ( t == ith_target_node( v, j ) ) {
                    // e_j = (v, j, t, w_j)
                    int e_j = ith_target_edge( v, j );
                    marked[ e_j ] = false;
                    --m;
                    int w_j = edge_info( e_j );
                    if ( w_i > w_j )
                        edge_info( e_i , w_j );
                }  // if t
            }  // for j 
        }  // for i 
    }  // for v

    if ( m == m_orig ) {
        printf( "No change in the input file. No output is generated.\n" );
    } else {
        printf( "p filtered %d %d\n", n, m  );
        for ( int e = 0; e < num_edges(); ++e ) {
            if ( marked[ e ] ) {
                int w = edge_info( e );
                if ( ! min_version )
                    w = -w;
                printf( "a %d %d %d\n", 1 + source( e ), 1 + target( e ), w  );
            }
        }
    }

    delete [] marked;
}  // print_no_duplicates

///////////////////////////////////////////////////////////////////////

void
print_components( ad_graph< cninfo >& cg )
{
    for ( int v = 0; v < cg.num_nodes(); ++v ) {
        printf( "\nSCC # %d\n", v  );
        cg.node_info( v ).comp->print();
    }
}

void
clear_components( ad_graph< cninfo >& cg )
{
    for ( int v = 0; v < cg.num_nodes(); ++v )
        delete cg.node_info( v ).comp;
}

// Find the component properties (number of nodes, number of edges)
// before dynamic allocation of lists (or arrays) to hold nodes,
// edges, etc.
void
find_component_props( ad_graph< cninfo >& cg, 
                      const ad_graph< ninfo >& g, int *comp_num )
{
    // Allocate only the node list for the component graph.
    cg.alloc_lists( true, false );

#ifdef DEBUG
    try {
#endif
        // Create nodes in the component graph.
        for ( int v = 0; v < cg.num_nodes(); ++v )
            cg.ins_node( cninfo( new ad_graph< ninfo > ) );
#ifdef DEBUG
    } catch ( std::bad_alloc &except ) {
        printf( "ERROR: Memory exhausted.\n" );
        abort();
    }
#endif
  
    // Compute the number of nodes of each SCC.
    for ( int v = 0; v < g.num_nodes(); ++v )
        cg.node_info( comp_num[ v ] ).comp->inc_num_nodes();

    // Compute the number of edges of the component graph and each SCC.
    for ( int e = 0; e < g.num_edges(); ++e ) {
        int s = g.source( e );
        int t = g.target( e );

        if ( comp_num[ s ] == comp_num[ t ] )
            cg.node_info( comp_num[ s ] ).comp->inc_num_edges();
        else  // Possible more edges than necessary.
            cg.inc_num_edges();
    }

#ifdef ADD_SOURCE_NODE
    // Add 1 for each additional source node.
    for ( int v = 0; v < cg.num_nodes(); ++v )
        cg.node_info( v ).comp->inc_num_nodes();
#endif

    // Allocate lists for each SCC.
    for ( int v = 0; v < cg.num_nodes(); ++v )
        cg.node_info( v ).comp->alloc_lists();

    // Allocate lists ( excluding the node list ) for the component graph.
    cg.alloc_lists( false );
}  // find_component_props

// Create the component graph cg, create each component, and attach
// them to the nodes of the component graph.
bool
find_components( ad_graph< cninfo >& cg, 
                 const ad_graph< ninfo >& g, 
                 bool has_self_loop, bool already_sc )
{
    // Construct the component graph cg from g. 

    if ( already_sc ) {
        // This is a shortcut. If g is already known to be strongly
        // connected, this shortcut copies g directly to the single node
        // of cg. There is no error checking here. That is, g must be as
        // it is known.

        cg.set_num_nodes( 1 );
        cg.set_num_edges( 0 );

#ifdef PRINT_SCC
        printf( "The component graph has %d SCCs.\n", cg.num_nodes()  );
#endif

        // Check to see if g is acyclic. If so, exit.
        if ( 0 == g.num_edges() )
            return true; 

        cg.alloc_lists( true, true );

        // Create the only node in cg.
        cg.ins_node( cninfo( new ad_graph< ninfo > ) );

        // Set parameters of the only SCC.
        cg.node_info( 0 ).comp->set_num_nodes( g.num_nodes() );
        cg.node_info( 0 ).comp->set_num_edges( g.num_edges() );

#ifdef ADD_SOURCE_NODE
        // The source node will all be numbered 0. The source node is
        // needed for Karp-Orlin and Young-Tarjan-Orlin algorithms.
        cg.node_info( 0 ).comp->inc_num_nodes();
#endif

        // Create lists.
        cg.node_info( 0 ).comp->alloc_lists();

#ifdef CYCLE_MEAN_VERSION
        for ( int e = 0; e < g.num_edges(); ++e ) {
            int s = g.source( e );
            int t = g.target( e );
#ifdef ADD_SOURCE_NODE
            cg.node_info( 0 ).comp->ins_edge( s + 1, t + 1, g.edge_info( e ) );
#else
            cg.node_info( 0 ).comp->ins_edge( s, t, g.edge_info( e ) );
#endif
        }

#else
        for ( int e = 0; e < g.num_edges(); ++e ) {
            int s = g.source( e );
            int t = g.target( e );
#ifdef ADD_SOURCE_NODE
            cg.node_info( 0 ).comp->ins_edge( s + 1, t + 1, g.edge_info( e ), g.edge_info2( e ) );
#else
            cg.node_info( 0 ).comp->ins_edge( s, t, g.edge_info( e ), g.edge_info2( e ) );
#endif
        }
#endif

    
        cg.node_info( 0 ).comp->build_adj();
        cg.build_adj();
    
    } else {

        // comp_num[v] holds the SCC where v belongs to.
        int *comp_num = new int[ g.num_nodes() ];

        // Set comp_num[v] for each node v.
        int nsccs = num_sccs( g, comp_num );

        // Set the number of nodes of the component graph.
        cg.set_num_nodes( nsccs );
#ifdef PRINT_SCC
        printf( "The component graph has %d SCCs.\n", cg.num_nodes()  );
#endif

        // Check to see if g is acyclic. If so, exit.
        if ( ( cg.num_nodes() == g.num_nodes() ) && ( !has_self_loop ) )
            return true;  

        // Determine the properties of SCCs.
        find_component_props( cg, g, comp_num );

        // Create a map from the SCC nodes to the nodes in the original
        // input graph, i.e., graph g.
        int *scc2orig = new int[ g.num_nodes() ];

#ifdef ADD_SOURCE_NODE
        // Create the source node in each SCC. The source node will all be
        // numbered 0. The source node is needed for Karp-Orlin and
        // Young-Tarjan-Orlin algorithms.
        for ( int v = 0; v < cg.num_nodes(); ++v )
            cg.node_info( v ).comp->ins_node();
#endif

        // Create (the remaining) nodes in each SCC:
        for ( int v = 0; v < g.num_nodes(); ++v )
            scc2orig[ v ] = cg.node_info( comp_num[ v ] ).comp->ins_node();

        // Create edges in the component graph and the SCCs. If the end
        // nodes of an edge have the same comp_num, insert them into the
        // same SCC; otherwise, add an edge into the component graph between
        // the corresponding SCCs if these SCCs are not adjacent already.
        // NOTE: I have removed the adjacency code, i.e., we insert an edge
        // between two SCCs even if they are already adjacent. In fact, we
        // do not check this at all.
#ifdef CYCLE_MEAN_VERSION
        for ( int e = 0; e < g.num_edges(); ++e ) {
            int s = g.source( e );
            int t = g.target( e );

            if ( comp_num[ s ] == comp_num[ t ] ) {
                cg.node_info( comp_num[ s ] ).comp->ins_edge( scc2orig[ s ], 
                                                              scc2orig[ t ], 
                                                              g.edge_info( e ) );
            } else {
                // or else if ( !cg.adjacent_slow( comp_num[ s ], comp_num[ t ] ) )
                cg.ins_edge( comp_num[ s ], comp_num[ t ], 0 );
            }
        }

#else
        for ( int e = 0; e < g.num_edges(); ++e ) {
            int s = g.source( e );
            int t = g.target( e );

            if ( comp_num[ s ] == comp_num[ t ] ) {
                cg.node_info( comp_num[ s ] ).comp->ins_edge( scc2orig[ s ], 
                                                              scc2orig[ t ], 
                                                              g.edge_info( e ),
                                                              g.edge_info2( e ) );
            } else {
                // or else if ( !cg.adjacent_slow( comp_num[ s ], comp_num[ t ] ) )
                cg.ins_edge( comp_num[ s ], comp_num[ t ], 0, 0 );
            }
        }
#endif


        delete [] scc2orig;
        delete [] comp_num;

        // Build the adjacency lists for each SCC.
        for ( int v = 0; v < cg.num_nodes(); ++v )
            cg.node_info( v ).comp->build_adj();

        // Build the adjacency lists of the component graph cg.
        cg.build_adj();
    }

    return false;  // true if g is acyclic.
}  // find_components

// Find which SCC each node blongs and return the info in
// comp_num. Also return the number of SCCs. The algorithm used is the
// one in the book by CLR. This algorithm uses depth-first search.
int
num_sccs( const ad_graph< ninfo >& g, int *comp_num )
{
    int n = g.num_nodes();

    int  *node_list = new int[ n ];
    bool *not_visited = new bool[ n ];

    for ( int v = 0; v < n; ++v ) {
        not_visited[ v ] = true;
        comp_num[ v ] = -1;
    }

    int nvisited = n;
    for ( int v = 0; v < n; ++v ) {
        if ( not_visited[ v ] )
            traverse_via_outedges( g, v, not_visited, nvisited, node_list );
    }

    delete [] not_visited;

    int num_comps = 0;
    for ( int v = 0; v < n; ++v ) {
        int u = node_list[ v ];
        if ( -1 == comp_num[ u ] ) {
            traverse_via_inedges( g, u, num_comps, comp_num );
            ++num_comps;
        }
    }

    delete [] node_list;
  
    return num_comps;
}

// Traverse the graph g using the outgoing edges. This algorithm is a
// top-level algorithm calling an iterative or recursive depth-first
// search algorithm.
void
traverse_via_outedges( const ad_graph< ninfo >& g,
                       int start_v, bool *not_visited,
                       int& n, int *node_list )
{
#ifdef DFS_RECUR
    traverse_via_outedges_recur(g, start_v, not_visited, n, node_list);
#else
    traverse_via_outedges_iter(g, start_v, not_visited, n, node_list);
#endif
}

// Traverse the graph g using the outgoing edges. This algorithm is a
// iterative depth-first search algorithm.
void
traverse_via_outedges_iter( const ad_graph< ninfo >& g,
                            int start_v, bool *not_visited,
                            int& n, int *node_list )
{
    // Structure to store node and edge index on stack
    struct stack_entry {
        int node;
        int edge_index;
    };
    
    // Allocate stack - maximum size needed is number of nodes
    int stack_size = g.num_nodes();
    stack_entry* stack = new stack_entry[ stack_size ];
    int stack_top = 0;  // Index of top of stack
    
    // Push initial node
    stack[ stack_top ].node = start_v;
    stack[ stack_top ].edge_index = 0;
    stack_top++;
    not_visited[ start_v ] = false;
    
    while ( stack_top > 0 ) {
        // Get reference to top stack entry
        stack_entry& current = stack[ stack_top - 1 ];
        
        // If we've processed all out-edges of current node
        if ( current.edge_index >= g.outdegree( current.node ) ) {
            node_list[ --n ] = current.node;
            stack_top--;  // Pop from stack
            continue;
        }
        
        // Get next unvisited target node
        int t = g.ith_target_node( current.node, current.edge_index );
        current.edge_index++;
        
        if ( not_visited[ t ] ) {
            not_visited[ t ] = false;
            // Push new node onto stack
            stack[ stack_top ].node = t;
            stack[ stack_top ].edge_index = 0;
            stack_top++;
        }
    }
    
    // Free allocated memory
    delete[] stack;
} // traverse_via_outedges_iter

// Traverse the graph g using the outgoing edges. This algorithm is a
// recursive depth-first search algorithm.
void
traverse_via_outedges_recur( const ad_graph< ninfo >& g, 
                             int v, bool *not_visited, 
                             int& n, int *node_list )
{
    not_visited[ v ] = false;

    for ( int i = 0; i < g.outdegree( v ); ++i ) {
#ifdef DEBUG
        assert( g.ith_target_edge( v, i ) != -1 );
#endif

        int t = g.ith_target_node( v, i );
        if ( not_visited[ t ] )
            traverse_via_outedges_recur( g, t, not_visited, n, node_list );
    }

    node_list[ --n ] = v;
}  // traverse_via_outedges_recur

// Traverse the graph g using the incoming edges. This algorithm is a
// top-level algorithm calling an iterative or recursive depth-first
// search algorithm.
void
traverse_via_inedges( const ad_graph< ninfo >& g,
                      int start_v, int num_comps, int *comp_num)
{
#ifdef DFS_RECUR
    traverse_via_inedges_recur(g, start_v, num_comps, comp_num);
#else
    traverse_via_inedges_iter(g, start_v, num_comps, comp_num);
#endif
   
} // traverse_via_inedges

// Traverse the graph g using the incoming edges. This algorithm is an
// iterative depth-first search algorithm.
void
traverse_via_inedges_iter( const ad_graph< ninfo >& g,
                           int start_v, int num_comps, int *comp_num)
{
    // Structure to store node and edge index on stack
    struct stack_entry {
        int node;
        int edge_index;
    };
    
    // Allocate stack - maximum size needed is number of nodes
    int stack_size = g.num_nodes();
    stack_entry* stack = new stack_entry[ stack_size ];
    int stack_top = 0;  // Index of top of stack
    
    // Push initial node
    stack[ stack_top ].node = start_v;
    stack[ stack_top ].edge_index = 0;
    stack_top++;
    comp_num[ start_v ] = num_comps;
    
    while ( stack_top > 0 ) {
        // Get reference to top stack entry
        stack_entry& current = stack[ stack_top - 1 ];
        
        // If we've processed all in-edges of current node
        if ( current.edge_index >= g.indegree( current.node ) ) {
            stack_top--;  // Pop from stack
            continue;
        }
        
        // Get next source node
        int s = g.ith_source_node( current.node, current.edge_index );
        current.edge_index++;
        
        if ( comp_num[ s ] == -1 ) {
            comp_num[ s ] = num_comps;
            // Push new node onto stack
            stack[ stack_top ].node = s;
            stack[ stack_top ].edge_index = 0;
            stack_top++;
        }
    }
    
    // Free allocated memory
    delete[] stack;
} // traverse_via_inedges_iter

// Traverse the graph g using the incoming edges. This algorithm is a
// recursive depth-first search algorithm.
void
traverse_via_inedges_recur( const ad_graph< ninfo >& g, 
                            int v, int num_comps, int *comp_num )
{
    comp_num[ v ] = num_comps;

    for ( int i = 0; i < g.indegree( v ); ++i ) {
#ifdef DEBUG
        assert( g.ith_source_edge( v, i ) != -1 );
#endif

        int s = g.ith_source_node( v, i );
        if ( -1 == comp_num[ s ] )
            traverse_via_inedges_recur( g, s, num_comps, comp_num );
    }
}  // traverse_via_inedges_recur

void
generate_part_for_all_components( ad_graph< cninfo >& cg,
                                  ginfo& gi,
                                  const args_t& args )
{
    for ( int v = 0; v < cg.num_nodes(); ++v ) {
        ad_graph< ninfo > *scc = cg.node_info( v ).comp;
        if ( scc->num_edges() ) {
            scc->generate_part( gi, args );
        }
    }
}  // generate_for_all_components

///////////////////////////////////////////////////////////////////////

float
find_min_cycle_ratio_for_components( const ad_graph< cninfo >& cg,
                                     int plus_infinity )
{
    float lambda = ( float ) plus_infinity;
    for ( int v = 0; v < cg.num_nodes(); ++v ) {

        const ad_graph< ninfo > *scc = cg.node_info( v ).comp;

        if ( scc->num_edges() ) {

#if PRINT_SCC
            // Only print info for non-trivial SCCs.
            printf( "Processing SCC # = %d with n= %d m= %d\n", 
                    v, scc->num_nodes(), scc->num_edges() );
#endif

            float lambda_for_scc = find_min_cycle_ratio_for_scc( scc, plus_infinity, lambda );
            min2( lambda, lambda_for_scc );

#if PRINT_SCC
            printf( "Lambda for SCC# %d is %10.2f\n", v, lambda_for_scc );
#endif

#if PRINT_SCC
            // Only print info for non-trivial SCCs.
            printf( "New lambda= %10.2f\n", lambda );
#endif

        }  // if
    }  // for

    return lambda;
}  // find_min_cycle_ratio_for_components

///////////////////////////////////////////////////////////////////////

// Instantiate instances of ad_graph with info types.
template class ad_graph< ninfo >;
template class ad_graph< cninfo >;

// End of file
