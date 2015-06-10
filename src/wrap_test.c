#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "quadphase.h"
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

main(argc, argv)
int argc;
char **argv; {
    double a,b,r1,r2,rd;
    double wrap();
    double wrapi();
    double arc();
    int i,j,nx;
    
    nx = 100;
        
    for (i = -1*nx; i < nx; i++){
        a = i * 2.0 * M_PI/30;
        r1 = wrap(a);
        r2 = wrapi(a);
        rd = r1-r2;
        fprintf(stdout,"%5d %10.4f %10.4f %10.4f %10.4f\n",i,a,r1,r2,rd);
    }
    
    exit(0);
}


