/* testints

 *  example ./pha2qls <inputfile.pha> xsize ysize [-o <outputfile>]
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "quadphase.h"


int main()
 {

    signed char c, c1;
    signed short i2a;
    /*20121004    signed short lsthdr[6]; */ /* header record */
    signed int lsthdr[3]; /* header record */
    int iopt,ierr;

    unsigned int i,j;
    signed int ii;

    for (ii = -128; ii < 129; ii++){
        c = (signed char) ii;
        fprintf(stderr,"ii = %3d c = %3d\n",ii,c);
        }
    return(0);
}

/*help int8
 int8 Convert to signed 8-bit integer.
    I = int8(X) converts the elements of the array X into signed 8-bit
    integers. X can be any numeric object, such as a DOUBLE. The values
    of an int8 range from -128 to 127, or INTMIN('int8') to INTMAX('int8')*/