//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
#ifndef AD_GLOBALS_INCLUDED
#define AD_GLOBALS_INCLUDED

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef NULL
#define NULL 0
#endif

const float SMALL_EPSILON = 0.001;
const float EPSILON  = 0.01;
const float EPSILON2 = EPSILON / 2;

const int   MAX_STR_SIZE   =  255;
const int   MAX_LINE_SIZE  =  255;
const int   MAX_ALLOC_SIZE = ( 1 << 15 );  // Amount of memory to allocate.

const int   SOURCE = 0;

#ifdef REP_COUNT
const int COUNT_LEN = 6;
extern int *count;
#endif

extern int ( *dist_func )( int, int );

template< class T >
inline
T
max( T x, T y ) 
{
    return ( x > y ? x : y );
}

template< class T >
inline
T
min( T x, T y )
{
    return ( x < y ? x : y );
}

template< class T >
inline
void
max2( T& x, T y ) 
{
    if ( x < y )
        x = y;
}

template< class T >
inline
void
min2( T& x, T y )
{
    if ( x > y )
        x = y;
}

template< class T >
inline
T
abs_val( T x )
{
    return ( x > 0 ? x : -x );
}

template< class T >
inline
T
fabs_val( T x )
{
    return ( x > 0.0 ? x : -x );
}

template< class T >
inline
void
init_table( T *ptr, int from, int to, T val )
{
    for ( int i = from; i <= to; ++i )
        ptr[ i ] = val;
}

template< class T >
inline
T
compute_max( T *ptr, int from, int to, T min_val )
{
    T max_val = min_val;
    for ( int i = from; i <= to; ++i )
        max2( max_val, ptr[ i ] );

    return max_val;
}

template< class T >
inline
T
compute_min( T *ptr, int from, int to, T max_val )
{
    T min_val = max_val;
    for ( int i = from; i <= to; ++i )
        min2( min_val, ptr[ i ] );

    return min_val;
}

// Argument information structure (see ad_util.cc for usage).
typedef struct {
    int  mode;
    char input_file[ MAX_STR_SIZE ];
    char dump_file[ MAX_STR_SIZE ];
    bool min_version;
    int  offset;
    int  nruns;
    int  nnodes, nedges;
    int  which_dist;
    int  w1, w2;  // Parameters for weights.
    int  t1, t2;  // Parameters for transit times.
    long seed;
} args_t;

#endif
