//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
#ifndef AD_UTIL_H
#define AD_UTIL_H

#ifdef REP_COUNT
void 
begin_count();
void 
print_count();
void 
end_count();
#endif

float 
used_time();

// Random number generation.
long
randomize( long seed );
int
uniform_dist( int min, int max );
int
exp_dist( int mean, int dummy );
int
normal_dist( int mean, int sdev );

// Argument parsing.
void
parse_args( int argc, char *argv[], args_t& args );

#endif
