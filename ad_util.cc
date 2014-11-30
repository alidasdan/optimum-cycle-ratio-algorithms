//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
// For used_time()
#include <limits.h>      // For CLK_TCK.
#include <time.h>        // For CLOCKS_PER_SEC, time()
#include <sys/times.h>   // For times().

// For log.
#include <math.h>

#include "ad_globals.h"
#include "ad_util.h"

#ifdef REP_COUNT
int *count;

// Initialize counters.
void 
begin_count()
{
    count = new int[COUNT_LEN];
    for ( int i = 0; i < COUNT_LEN; ++i )
        count[ i ] = 0;
}

// Print counters.
void 
print_count()
{
    printf( "REP_COUNT " );
    for ( int i = 0; i < COUNT_LEN; ++i )
        printf( "%d ", count[ i ] );
}

// Destroy counters.
void 
end_count()
{
    delete [] count;
}
#endif

// Measure time.
float 
used_time()
{
    struct tms now;

    times( &now );
#ifdef CLK_TCK
    return (float)( now.tms_utime + now.tms_stime ) / CLK_TCK;
#else
    // return (float)( now.tms_utime + now.tms_stime ) / sysconf(_SC_CLK_TCK)
    return (float)( clock() ) / CLOCKS_PER_SEC;
#endif
}

// Initialize the random number generator.
long
randomize( long seed )
{
    if ( -1 == seed )
        seed = time( 0 );
    srand48( seed );
    return seed;
}  // randomize

int 
uniform_dist( int min, int max )
{
    return ( int ) ( min + ( max - min + 1 ) * drand48() );
}

/*ARGSUSED1*/
int 
exp_dist( int mean, int dummy )
{
    double num01;

    while ( 1 ) {
        num01 = drand48();
        if ( 0.0 < num01 )
            break;
    }
  
    return ( int ) ( mean * ( -log( num01 ) ) );
}  // exp_dist

int
normal_dist( int mean, int sdev ) 
{
    double sum01 = 0.0;

    for ( int i = 0; i < 12; ++i )
        sum01 += drand48();

    return ( int ) ( mean + sdev * ( sum01 - 6.0 ) );
}  // normal_dist

// Argument parsing.
void
parse_args( int argc, char *argv[], args_t& args )
{
    // Mode 0: A file is given. It contains a graph with all the
    // attributes defined. There is no need to generate anything.

    // Mode 1: A file is given, It contains a graph with all the
    // attributes defined. However, the weights must be regenerated under
    // the required distribution.

    // Mode 2: No file is given. The graph as well as all its attributes
    // must be generated under the required distribution.

    // Format: [input_file] [-m 0/1/2] [-v 0/1] [-n nruns] [-o offset]
    // [-d 0/1/2] [-p n m] [-w w1 w2] [-t t1 t2] [-s seed] [-f dump_file]

    args.mode = 0; // 0, 1, 2
    strcpy( args.input_file, "" );
    strcpy( args.dump_file, "" );
    args.min_version = true; // true, false
    args.offset = 0; 
    args.nruns = 1; 
    args.nnodes = 0;  // Required arg under mode 2
    args.nedges = 0;  // Required arg under mode 2
    args.which_dist = 0; // 0=uniform, 1=normal, 2=exponential
    args.w1 = 1;
    args.w2 = 300;
    args.t1 = 1;
    args.t2 = 10;
    args.seed = -1;

    int i = 1;

    bool error_found = false;

    if ( i < argc ) {
        // IMPORTANT: The input file must be the 1st argument,
        // followed by other flags, if any.
        strcpy( args.input_file, argv[ i ] );
        if ( '-' == args.input_file[ 0 ] ) {
            args.mode = 2;
        } else {
            ++i;
        }
    } else {
        error_found = true;
    }

    while ( i < argc ) {
        if ( error_found )
            break;

        if ( !strcmp( argv[ i ], "-m" ) ) {
            if ( 1 >= ( argc - i ) ) 
                error_found = true;
            args.mode = atoi( argv[ i + 1 ] );
            switch ( args.mode ) {
            case 0: case 1: case 2: break;
            default:
                printf( "\nERROR: '-m' must be followed by 0, 1, or 2.\n" );
                error_found = true;
            }
            i += 2;
        } else if ( !strcmp( argv[ i ], "-v" ) ) {
            if ( 1 >= ( argc - i ) ) 
                error_found = true;
            int tmp = atoi( argv[ i + 1 ] );
            switch ( tmp ) {
            case 0: case 1: break;
            default:
                printf( "\nERROR: '-v' must be followed by 0 or 1.\n" );
                error_found = true;
            }
            args.min_version = ( bool ) tmp;
            i += 2;
        } else if ( !strcmp( argv[ i ], "-n" ) ) {
            if ( 1 >= ( argc - i ) ) 
                error_found = true;
            args.nruns = atoi( argv[ i + 1 ] );
            if ( args.nruns < 1 ) {
                printf( "\nERROR: '-n' must be followed by a positive integer.\n" );
                error_found = true;                
            }            
            i += 2;
        } else if ( !strcmp( argv[ i ], "-o" ) ) {
            if ( 1 >= ( argc - i ) ) 
                error_found = true;
            args.offset = atoi( argv[ i + 1 ] );
            if ( 0 == args.offset ) {
                printf( "\nERROR: '-o' must be followed by a non-zero integer.\n" );
                error_found = true;
            }
            i += 2;
        } else if ( !strcmp( argv[ i ], "-d" ) ) {
            if ( 1 >= ( argc - i ) ) 
                error_found = true;
            args.which_dist = atoi( argv[ i + 1 ] );
            switch ( args.which_dist ) {
            case 0: case 1: case 2: break;
            default:
                printf( "\nERROR: '-d' must be followed by 0, 1, or 2.\n" );
                error_found = true;
            }
            i += 2;
        } else if ( !strcmp( argv[ i ], "-p" ) ) {
            if ( 2 >= ( argc - i ) ) 
                error_found = true;
            args.nnodes = atoi( argv[ i + 1 ] );
            args.nedges = atoi( argv[ i + 2 ] );
            if ( 0 == args.nnodes ) {
                printf( "\nERROR: 'n > 0' must hold.\n" );
                error_found = true;
            }
            if ( ( 1 == args.nnodes ) && ( 1 == args.nedges ) ) {
                printf( "\nERROR: self-loops are not allowed.\n" );
                error_found = true;
            }
            if ( args.nedges < args.nnodes ) {
                printf( "\nERROR: 'm >= n' must hold.\n" );
                error_found = true;
            }
            i += 3;
        } else if ( !strcmp( argv[ i ], "-w" ) ) {
            if ( 2 >= ( argc - i ) ) 
                error_found = true;
            args.w1 = atoi( argv[ i + 1 ] );
            args.w2 = atoi( argv[ i + 2 ] );
            i += 3;
        } else if ( !strcmp( argv[ i ], "-t" ) ) {
            if ( 2 >= ( argc - i ) ) 
                error_found = true;
            args.t1 = atoi( argv[ i + 1 ] );
            args.t2 = atoi( argv[ i + 2 ] );
            i += 3;
        } else if ( !strcmp( argv[ i ], "-s" ) ) {
            if ( 1 >= ( argc - i ) ) 
                error_found = true;
            args.seed = atol( argv[ i + 1 ] );
            if ( args.seed <= 0 ) {
                printf( "\nERROR: '-s' must be followedb by a non-zero integer.\n" );
                error_found = true;
            }
            i += 2;
        } else if ( !strcmp( argv[ i ], "-f" ) ) {
            if ( 1 >= ( argc - i ) ) 
                error_found = true;
            strcpy( args.dump_file, argv[ i + 1 ] );
            i += 2;
        } else {
            error_found = true;
            break;
        }
    }

    if ( 2 == args.mode ) {
        if ( ( args.nnodes <= 0 ) || ( args.nedges < 0 ) ) {
            printf( "ERROR: Require 'nnodes > 0' and 'nedges >= 0'. \n" );
            abort();
        }        
    }

    if ( error_found ) {

        if ( argc > 1 ) {
            printf( "\nERROR: Cannot parse correctly. Your command line was:");
            for ( int i = 0; i < argc; ++i ) {
                printf( " %s", argv[ i ] );
            }
            printf( "\n" );
        }

        printf( "\nUsage: %s\n", argv[ 0 ] );
        printf( "   [input_file]     file to read input graph (MUST BE 1ST ARG)\n" );
        printf( "   [-m/ode 0/1/2]   read or generate -- see ad_util.cc for details\n" );
        printf( "   [-v/ersion 0/1]  min or max version\n" );
        printf( "   [-n nruns]       number of runs to perform\n" );
        printf( "   [-o offset]      subtract offset from every edge weight\n" );
        printf( "   [-p/aram n m]    num nodes and edges for graph generation\n" );
        printf( "   [-d/ist 0/1/2]   distribution to use -- see ad_util.c for details\n" );
        printf( "   [-w/eight w1 w2] min and max weight bounds\n" );
        printf( "   [-t/time t1 t2]  min and max transit time bounds\n" );
        printf( "   [-s seed]        random number generator seed\n" );
        printf( "   [-f dump_file]   file to dump output\n" );

        printf( "\nBelow are what is known from parsing:\n" );
        printf( "\tmode= %d\n", args.mode );
        printf( "\tinput file= %s\n", args.input_file );
        if ( args.min_version )
            printf( "\tversion= min\n" );
        else
            printf( "\tversion= max\n" );
        printf( "\toffset= %d\n", args.offset );
        printf( "\tnum runs= %d\n", args.nruns );
        printf( "\tn= %d\n", args.nnodes ); 
        printf( "\tm= %d\n", args.nedges ); 
        switch ( args.which_dist ) {
        case 0:
            printf( "\tdist= uniform\n" );
            break;
        case 1:
            printf( "\tdist= normal\n" );
            break;
        case 2:
            printf( "\tdist= exponential\n" );
            break;
        default:
            printf(" dist= unknown\n" );
        }
        printf( "\t[w1:w2]= [ %d : %d ]\n", args.w1, args.w2 ); 
        printf( "\t[t1:t2]= [ %d : %d ]\n", args.t1, args.t2 ); 
        printf( "\tseed= %ld\n", args.seed );
        printf( "\tdump file= %s\n", args.dump_file );

        exit( 0 );

    }
}  // parse_args

// End of file
