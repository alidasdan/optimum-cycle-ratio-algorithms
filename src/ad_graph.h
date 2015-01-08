//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
#ifndef AD_GRAPH_INCLUDED
#define AD_GRAPH_INCLUDED

#include "ad_globals.h"

// Forward declarations.
struct ginfo;
template< class ninfo_t > class ad_graph;

///////////////////////////////////////////////////////////////////////
// Node class:

template< class ninfo_t >
class ad_node {

    template< class ninfo_t2 >
    friend class ad_graph;

public:

    // Constructors:
    ad_node() 
    { 
        set_node( 0, 0, -1, -1 );
    }
#if 0
    ad_node( int id, int od, int ia, int oa )
    {
        set_node( id, od, ia, oa );
    }
#endif

    // Get and set functions:
    int indegree() const
    {
        return indeg;
    }
    void indegree( int d )
    {
        indeg = d;
    }
    int outdegree() const
    {
        return outdeg;
    }
    void outdegree( int d )
    {
        outdeg = d;
    }
    int degree() const
    {
        return indeg + outdeg;
    }

    ninfo_t info() const
    {
        return inf;
    }
    void info( const ninfo_t& ni )
    {
        inf = ni;
    }

    // Printing the info field.
    // pre: Existence of print( ninfo_t );
    void print_info() const
    {
        print( inf );
    }

private:
    int     indeg;  // In-degree.
    int     outdeg; // Out-degree.
    int     inadj;  // Index of the first edge in adjlist array.
    int     outadj; // Index of the first edge in adjlist array.
    ninfo_t inf;    // Info field.

    // Constructor helpers.
    void set_node( int id, int od, int ia, int oa )
    {
        indeg = id; outdeg = od; inadj = ia; outadj = oa;
    }

    // To get and set the first edge in- or out-index:
    int first_inedge() const
    {
        return inadj;
    }
    void first_inedge( int a )
    {
        inadj = a;
    }
    int first_outedge() const
    {
        return outadj;
    }
    void first_outedge( int a )
    {
        outadj = a;
    }

    // Increment in- and out-degree:
    void inc_indegree()
    {
        ++indeg;
    }
    void inc_outdegree()
    {
        ++outdeg;
    }
};  // ad_node

///////////////////////////////////////////////////////////////////////
// Edge class: ( I use the terms arc and edge interchangeably. )

class ad_edge {

    template< class ninfo_t >
    friend class ad_graph;

public:

    // Constructors:
    ad_edge()
    {
#ifdef DEBUG
        set_edge( -1, -1 );
#endif
    }
#if 0
    ad_edge( int s, int t )
    {
        set_edge( s, t );
    }
#endif

    // Get and set functions:
    int source() const
    {
        return src;
    }
    void source( int s )
    {
        src = s;
    }

    int target() const
    {
        return tar;
    }
    void target( int t )
    {
        tar = t;
    }

    int info() const
    {
        return inf;
    }
    void info( int ei )
    {
        inf = ei;
    }

#ifdef CYCLE_MEAN_VERSION
    // Printing the info fields.
    void print_info() const
    {
        printf( "%d %d", inf, inf2 );
    }
#else
    int info2() const
    {
        return inf2;
    }
    void info2( int ei2 )
    {
        inf2 = ei2;
    }
    void info( int ei, int ei2 )
    {
        inf = ei; inf2 = ei2;
    }

    // Printing the info fields.
    void print_info() const
    {
        printf( "%d %d", inf, inf2 );
    }
#endif

private:
    // Edge is src -> tar with info inf and inf2.
    int  src;  
    int  tar;
    int  inf;
#ifndef CYCLE_MEAN_VERSION
    int  inf2;
#endif

    // Constructor helper.
    void set_edge( int s, int t ) 
    {
        src = s; tar = t; 
    }
};  // ad_edge

///////////////////////////////////////////////////////////////////////
// InEdge class: ( I use the terms arc and edge interchangeably. )

class ad_inedge {

    template< class ninfo_t >
    friend class ad_graph;

public:

    // Constructors:
    ad_inedge()
    {
#ifdef DEBUG
        ine = -1;
#endif
    }

#ifdef CYCLE_MEAN_VERSION
    // Constructor helper.
    void set_edge( int e, int s, int ei ) 
    {
        ine = e; src = s; inf = ei;
    }
#else
    // Constructor helper.
    void set_edge( int e, int s, int ei, int ei2 ) 
    {
        ine = e; src = s; inf = ei; inf2 = ei2;
    }
#endif

    // Get and set functions:
    int source_edge() const
    {
        return ine;
    }
    void source_edge( int e )
    {
        ine = e;
    }

    int source() const
    {
        return src;
    }
    void source( int s )
    {
        src = s;
    }

    int info() const
    {
        return inf;
    }
    void info( int ei )
    {
        inf = ei;
    }

#ifdef CYCLE_MEAN_VERSION
    // Printing the info field.
    void print_info() const
    {
        printf( "%d", inf );
    }
#else
    int info2() const
    {
        return inf2;
    }
    void info2( int ei2 )
    {
        inf2 = ei2;
    }
    void info( int ei, int ei2 )
    {
        inf = ei; inf2 = ei2;
    }

    // Printing the info fields.
    void print_info() const
    {
        printf( "%d %d", inf, inf2 );
    }
#endif

private:
    // Except ine (src, tar), the rest of the data members below are
    // redundant. They are included for efficiency in that accessing an
    // inedge of a node tar (which is supposed to point to here) can
    // also return the attributes of the inedge itself w/o a subsequent
    // access to the edge array (or list).
    int ine;
    int src;  
    int inf;
#ifndef CYCLE_MEAN_VERSION
    int inf2;
#endif

    // Disable these functions.
    int target_edge() const;
    void target_edge( int e );
    int target() const;
    void target( int t );
};  // ad_inedge

///////////////////////////////////////////////////////////////////////
// Outedge class: ( I use the terms arc and edge interchangeably. )

class ad_outedge {

    template< class ninfo_t >
    friend class ad_graph;

public:

    // Constructors:
    ad_outedge()
    {
#ifdef DEBUG
        oute = -1;
#endif
    }

#ifdef CYCLE_MEAN_VERSION
    // Constructor helper.
    void set_edge( int e, int t, int ei ) 
    {
        oute = e; tar = t; inf = ei;
    }
#else
    // Constructor helper.
    void set_edge( int e, int t, int ei, int ei2 ) 
    {
        oute = e; tar = t; inf = ei; inf2 = ei2;
    }
#endif

    // Get and set functions:
    int target_edge() const
    {
        return oute;
    }
    void target_edge( int e )
    {
        oute = e;
    }

    int target() const
    {
        return tar;
    }
    void target( int t )
    {
        tar = t;
    }

    int info() const
    {
        return inf;
    }
    void info( int ei )
    {
        inf = ei;
    }

#ifdef CYCLE_MEAN_VERSION
    // Printing the info field.
    void print_info() const
    {
        printf( "%d", inf );
    }
#else
    int info2() const
    {
        return inf2;
    }
    void info2( int ei2 )
    {
        inf2 = ei2;
    }
    void info( int ei, int ei2 )
    {
        inf = ei; inf2 = ei2;
    }

    // Printing the info field.
    void print_info() const
    {
        printf( "%d %d", inf, inf2 );
    }
#endif

private:
    // Except oute (src, tar), the rest of the data members below are
    // redundant. They are included for efficiency in that accessing an
    // outedge of a node src (which is supposed to point to here) can
    // also return the attributes of the outedge itself w/o a subsequent
    // access to the edge array (or list).
    int oute;
    int tar;  
    int inf;
#ifndef CYCLE_MEAN_VERSION
    int inf2;
#endif

    // Disable these functions.
    int source_edge() const;
    void source_edge( int e );
    int source() const;
    void source( int s );
};  // ad_outedge

///////////////////////////////////////////////////////////////////////
// Graph class: 

template< class ninfo_t >
class ad_graph {
public:

    // Constructors:
    ad_graph()
    {
        create( 0, 0 );
    }
#if 0
    ad_graph( int n, int m )
    {
        create( n, m );
    }
#endif

    // Destructor:
    ~ad_graph()
    {
        clear();
    }

    // Get and set functions for the graph:
    int num_nodes() const
    {
        return nnodes;
    }
    int num_edges() const
    {
        return nedges;
    }

    // Get and set functions for the nodes and edges:
    int degree( int v ) const
    {
        return indegree( v ) + outdegree( v );
    }
    int indegree( int v ) const
    {
        return nlist[ v ].indegree();
    }
    void indegree( int v, int d ) 
    {
        nlist[ v ].indegree( d );
    }
    int outdegree( int v ) const
    {
        return nlist[ v ].outdegree();
    }
    void outdegree( int v, int d )
    {
        nlist[ v ].outdegree( d );
    }

    int source( int e ) const
    {
        return elist[ e ].source();
    }
    int target( int e ) const
    {
        return elist[ e ].target();
    }

    ninfo_t node_info( int v ) const
    {
        return nlist[ v ].info();
    }
    int node_info( int v, const ninfo_t& ni )
    {
        nlist[ v ].info( ni );
        return v;
    }
    int edge_info( int e ) const
    {
        return elist[ e ].info();
    }
    int edge_info( int e, int ei )
    {
        elist[ e ].info( ei );
        return e;
    }

#ifndef CYCLE_MEAN_VERSION
    int edge_info2( int e ) const
    {
        return elist[ e ].info2();
    }
    int edge_info2( int e, int ei2 )
    {
        elist[ e ].info2( ei2 );
        return e;
    }
    int edge_info( int e, int ei, int ei2 )
    {
        elist[ e ].info( ei, ei2 );
        return e;
    }
#endif

    int ith_source_node( int v, int i ) const
    {
        return inlist[ nlist[ v ].first_inedge() + i ].source();
    }
    int ith_source_edge( int v, int i ) const
    {
        return inlist[ nlist[ v ].first_inedge() + i ].source_edge();
    }
    int ith_source_edge_info( int v, int i ) const
    {
        return inlist[ nlist[ v ].first_inedge() + i ].info();
    }
#ifndef CYCLE_MEAN_VERSION
    int ith_source_edge_info2( int v, int i ) const
    {
        return inlist[ nlist[ v ].first_inedge() + i ].info2();
    }
#endif
    ad_inedge ith_source( int v, int i ) const
    {
        return inlist[ nlist[ v ].first_inedge() + i ];
    }

    int ith_target_node( int v, int i ) const
    {
        return outlist[ nlist[ v ].first_outedge() + i ].target();
    }
    int ith_target_edge( int v, int i ) const
    {
        return outlist[ nlist[ v ].first_outedge() + i ].target_edge();
    }
    int ith_target_edge_info( int v, int i ) const
    {
        return outlist[ nlist[ v ].first_outedge() + i ].info();
    }
#ifndef CYCLE_MEAN_VERSION
    int ith_target_edge_info2( int v, int i ) const
    {
        return outlist[ nlist[ v ].first_outedge() + i ].info2();
    }
#endif

    ad_outedge ith_target( int v, int i ) const
    {
        return outlist[ nlist[ v ].first_outedge() + i ];
    }

#if 0
    ninfo_t ith_source_node_info( int v, int i ) const
    {
        return node_info( ith_source_node( v, i ) );
    }
    ninfo_t ith_target_node_info( int v, int i ) const
    {
        return node_info( ith_target_node( v, i ) );
    }
#endif

    // Return true if the nodes u and v are adjacent; false
    // otherwise. The fast version can only be called after build_adj is
    // called. The slow version is okay to call any time.
    bool adjacent_fast( int u, int v ) const
    {
        for ( int i = 0; i < outdegree( u ); ++i ) {
            if ( ith_target_node( u, i ) == v )
                return true;
        }
        for ( int i = 0; i < outdegree( v ); ++i ) {
            if ( ith_target_node( v, i ) == u )
                return true;
        }
        return false;
    }
    bool adjacent_slow( int u, int v ) const
    {
        for ( int e = 0; e < num_edges(); ++e ) {
            if ( ( source( e ) == u && target( e ) == v ) ||
                 ( source( e ) == v && target( e ) == u ) )
                return true;
        }
        return false;
    }

    // Insertion functions:

    // Insert a node and return its number.
    int ins_node()
    {
        nlist[ ++cur_node ].set_node( 0, 0, -1, -1 );
        return cur_node;
    }
    int ins_node( const ninfo_t& ni )
    {
        return node_info( ins_node(), ni );
    }

    // Insert an edge and return its number.
    int ins_edge( int s, int t )
    {
        elist[ ++cur_edge ].set_edge( s, t );
        inc_outdegree( s );
        inc_indegree( t );
        return cur_edge;
    }
#ifdef CYCLE_MEAN_VERSION
    int ins_edge( int s, int t, int ei )
    {
        return edge_info( ins_edge( s, t ), ei );
    }
#else
    int ins_edge( int s, int t, int ei, int ei2 )
    {
        return edge_info( ins_edge( s, t ), ei, ei2 );
    }
#endif

    // I/O functions:

    // Return the total edge weight and determine if the graph has at
    // least one self-loop. If min_version is false, every edge weight
    // is negated in order to compute the max cycle mean. The graph is
    // either read from file_name or generated partially or completely.
    void read( ginfo& gi, const args_t& args );
    void generate_part( ginfo& gi, const args_t& args );
    void generate_all( ginfo& gi, const args_t& args );

    void print( bool all_out = true ) const;
    void fprint( const args_t& args ) const;

    // Build the adjacency info.
    void build_adj();

    // Functions to set num_nodes and num_edges during incremental
    // creation of the nodes and edges.
    void set_num_nodes( int n )
    {
        if ( not_already_built )
            nnodes = n;
    }
    void inc_num_nodes( int n = 1 )
    {
        if ( not_already_built )
            nnodes += n;
    }
    void set_num_edges( int m )
    {
        if ( not_already_built )
            nedges = m;
    }
    void inc_num_edges( int m = 1 )
    {
        if ( not_already_built )
            nedges += m;
    }

    void alloc_lists( bool must_alloc_nlist = true, bool must_alloc_elist = true )
    {
        if ( must_alloc_nlist )
            nlist = new ad_node<ninfo_t>[ nnodes ];

        if ( must_alloc_elist ) {
            elist = new ad_edge[ nedges ];
            inlist = new ad_inedge[ nedges ];
            outlist = new ad_outedge[ nedges ];
        }
    }  // alloc_lists

    void print_no_duplicates( bool min_version );

private:
    int nnodes;  // Number of nodes.
    int nedges;  // Number of edges.

    int cur_node; // Current node to be processed.
    int cur_edge; // Current edge to be processed.

    bool not_already_built; // To control allocation of lists.

    // Node, edge, and adjacency arrays ( inlist is for incoming edges
    // whereas outlist is for outcoming edges ).
    ad_node<ninfo_t>   *nlist;    
    ad_edge            *elist;    
    ad_inedge          *inlist;  
    ad_outedge         *outlist;  

private:
    // Get and set functions for the head pointer for the in- and
    // out-lists:
    int first_inedge( int v ) const
    {
        return nlist[ v ].first_inedge();
    }
    void first_inedge( int v, int ia )
    {
        nlist[ v ].first_inedge( ia );
    }

    int first_outedge( int v ) const
    {
        return nlist[ v ].first_outedge();
    }
    void first_outedge( int v, int oa )
    {
        nlist[ v ].first_outedge( oa );
    }

    // Increment in-, out-, and inout-degrees of nodes:
    void inc_inoutdegree( int e )
    {
        inc_outdegree( source( e ) );
        inc_indegree( target( e ) );
    }
    void inc_indegree( int v )
    {
        nlist[ v ].inc_indegree();
    }
    void inc_outdegree( int v )
    {
        nlist[ v ].inc_outdegree();
    }

    // Update num_nodes and num_edges, and prevent further updates.
    void update_nums()
    {
        // If cur_node and/or cur_edge are updated, then update nnodes
        // and/or nedges.
        if ( cur_node != -1 )
            nnodes = cur_node + 1;
        if ( cur_edge != -1 )
            nedges = cur_edge + 1;
        not_already_built = false;
    }  // update_nums

    // Create a graph with n nodes and m edges.
    void create( int n, int m )
    {
        nnodes = n;
        nedges = m;
        cur_node = -1;
        cur_edge = -1;

        not_already_built = true;

        if ( nnodes ) {
            alloc_lists();
        } else {
            nlist = NULL;
            elist = NULL;
            inlist = NULL;
            outlist = NULL;
        }
    }  // create

    // Clear the node, edge, and adjacency arrays.
    void clear() {
        if ( nnodes ) {
            delete [] nlist;
            delete [] elist;
            delete [] inlist;
            delete [] outlist;
        }
    }  // clear

};  // ad_graph

///////////////////////////////////////////////////////////////////////
// Information classes:

// Graph information.
struct ginfo {
    int  total_edge_weight;
#ifndef CYCLE_MEAN_VERSION
    int  total_trans_time;
#endif
    bool has_self_loop;
};

// Node information.
struct ninfo {
    void print() const
    {
        /* EMPTY */
    }
};

// Component information. The compiler doesn't inline operator= if I
// write the constructors in terms of this operator.
struct cninfo {

    // Constructors:
    cninfo() 
    {
        comp = NULL;
    }
    cninfo( const cninfo& ci )
    {
        comp = ci.comp;
    }
    cninfo( ad_graph< ninfo > *g )
    {
        comp = g;
    }

#if 0
    // Destructor.
    ~cninfo()
    {
        if ( comp ) {
            if ( comp->dec_num_refs() == 0 )
                delete comp;
        }
    }
#endif

    // Assignment operators:
    void operator=( const cninfo& ci )
    {
        comp = ci.comp;
    }
    void operator=( ad_graph< ninfo > *g )
    {
        comp = g;
    }

    // I/O functions.
    void print() const
    {
        comp->print();
    }

    // Pointer to the SCC.
    ad_graph< ninfo > *comp;
};

///////////////////////////////////////////////////////////////////////
// Strongly Connected Component (SCC) functions:

extern 
void
print_components( ad_graph< cninfo >& cg );

extern 
void 
clear_components( ad_graph< cninfo >& cg );

extern 
bool
find_components( ad_graph< cninfo >& cg, 
                 const ad_graph< ninfo >& g, 
                 bool has_self_loop, bool already_sc = false );

extern
void
generate_part_for_all_components( ad_graph< cninfo >& cg,
                                  ginfo& gi,
                                  const args_t& args );

///////////////////////////////////////////////////////////////////////
// Optimum Cycle mean (=ratio) functions:

// Find the min cycle ratio for a SCC g.
extern
float
find_min_cycle_ratio_for_scc( const ad_graph< ninfo > *g, 
                              int plus_infinity,
                              float lambda_so_far );

// Find the min cycle mean of the component graph cg by going over its
// SCCs using the previous function.
extern
float 
find_min_cycle_ratio_for_components( const ad_graph< cninfo >& cg, 
                                     int plus_infinity );

inline
float 
find_max_cycle_ratio_for_components( const ad_graph< cninfo >& cg,
                                     int plus_infinity )
{
    // Assuming that the edge weights are negated in the input graph.
    return -find_min_cycle_ratio_for_components( cg, plus_infinity );
}

float 
find_lambda_bound( const ad_graph< ninfo > *g, 
                   int plus_infinity, 
                   bool which );

inline
float
find_min_lambda( const ad_graph< ninfo > *g, 
                 int plus_infinity )
{
    return find_lambda_bound( g, plus_infinity, true );
}

inline
float
find_max_lambda( const ad_graph< ninfo > *g, 
                 int plus_infinity )
{
    return find_lambda_bound( g, plus_infinity, false );
}

#endif
