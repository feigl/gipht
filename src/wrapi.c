#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
/* #include "gmt.h" */

/*--------------------------------------------------------------------
 *			GMT CONSTANTS MACRO DEFINITIONS
 *--------------------------------------------------------------------*/

#ifndef TWO_PI
#define TWO_PI        6.28318530717958647692
#endif
#ifndef M_PI
#define M_PI          3.14159265358979323846
#endif
double wrapi (double x)
{
    double r;
    extern double round();
    signed char c;
    
    /*OPERATOR: WRAP 1 1 wrap (A). (A in radians). */
    /*
     * wrap a value in radians onto [-pi,pi]
     *
     * Kurt Feigl 2010-DEC-04
     */
    
    /* scale on to a 1-byte integer on [-128, +127] */
    c = (signed char)rint(256.0 * x / TWO_PI);
    /* scale back into radians */
    r =  (double) c;
    r = r * TWO_PI/256.0;
    return(r);
}
