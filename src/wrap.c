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
double wrap (double x)
{
    double r;
    extern double round();
    
    /*OPERATOR: WRAP 1 1 wrap (A). (A in radians). */
    /*
     * wrap a value in radians onto [-pi,pi]
     *
     * Kurt Feigl 2014-01-06
     */
    
    r=2.0*M_PI*(x/M_PI/2.0 - round(x/M_PI/2.0)); 
        
    /* need this to prevent ambiguity */
    if (M_PI - r < 2.0 * M_PI/256.0) {
        r = r - 2.0 * M_PI;
    }
    
    return(r);
}
