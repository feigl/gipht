/* QLS2PHA.C
 *
 * Reconstructs approximation of .pha file from quadtree partitioned
 * data (.qls) files created by pha2qls3.
 *
 * Each entry in the qls file describes a square patch with:
 * upper left coods, size, mean phase, X gradient, y gradient
 *
 * input: a QL .qls file ( GIPhT [Q]uadphase [L]ist format )
 *	QL .qls file structure: (binary 16-bit ints)
 *		header:	   'QL'	NP	NX	NY IThresh MINPIX		Signature Npatches Ncols Nrows
 *		patches:	i1	j1	nwidth	PV GX GY		Left_Upper_corner xySize Patch_value Gradient_X Gradient_Y
 * output: a signed 8-bit .pha file of wrapped phases
 *		an error is generated if a patch write exceeds nx*ny output array size
 * if  no output filename is specified  outoyt is built on input name
 *
 *
 * options:
 * -d# : to enable debugging
 *	# == 0	no debugging (default)
 *	# == 1	Warning to STDOUT if a patch write overwrites an existing patch
 *	# == 2	echos input qls data to STDOUT
 *	# == 3	combines 1 & 2
 *
 *  -F : "Flat" output - all points in each patch are set to the mean - gradients are ignored
 *  -o : <outputfilename> - specify output filename default is  qha_<inputqlsfile>.pha
 *  -V : verbose
 *
 * gcc qls2pha.c -lm -o qls2pha.mac64
 * gcc qls2pha.c -lm -o qls2pha.a64
 * gcc qls2pha.c -lm -o qls2pha.maci64
 *
 *
 * Author:	Lee Powell
 * Date:	2009-JUL-12
 * Version:	1.1
 *
 * Updated (from QLS2PHA3): Peter Sobol
 * Date:    2011-JUN-21
 * Now includes gradient values in the reconstruction of the .pha file
 * Updated input options
 * Fix file names
 * Add gradients as doubles
 *
 *
 * *  2011/07/19 Kurt documentation
 * *  2012/10/06 Kurt remove check on first line
 * *  2012/10/08 correct bug in Flat
 * *  2014/01/06 Kurt still looking for bug
 
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

static FILE *fpqls	 = NULL;
static FILE *fpout	 = NULL;
static FILE *fpaux	 = NULL;

static signed char *iwrk, *iaux;
static signed char *pout;
static long nx, ny, np, npix;

main(argc, argv)
int argc;
char **argv; {
//     unsigned short lstdat[6];
    signed short lstdat[6];
    signed int lsthdr[3]; /* header record */
    char *foutname, *flstname, *tp, *strtemp;
    long nw, nr, i, j, k;
    long debug=0;
//     unsigned short i1, j1, nwidth;
    long i1, j1, nwidth, i2, j2;
    double xmid, ymid;
//     signed short pv, tv, gx, gy;
    signed short tv;
    double pv, gx, gy;    
    short flatflag = 0, outflag =0;
    int iopt;
    double r;
    
    opterr = 0;
    strtemp= (char *) malloc( sizeof(char) * 256);
    fprintf(stderr, "qls2pha v. 1.2 of 2012-OCT-06\n\n");
    
    
    /* First ARG is the compressed phase list filename or "stdin" *
     * prog[0] arg[1] arg[2] arg[3] arg[4]
     * argc 1      2      3      4      5
     */
    if (argc > 1) {   // test & open the input filename or STDIN
        flstname = (char *) malloc(sizeof(char) * strlen(argv[0])+1);
        strcpy(flstname, argv[1]);
        
        if( strcasecmp("STDIN", flstname) == 0 ) {
            fpqls = stdin;
        }
        else{	// nope, must be then name for an input file
            if ((fpqls = fopen(flstname, "r")) == NULL) {
                fprintf(stderr, "Cannot open input phase list file %s\n", flstname);
                usage();
            }
            
        }
        /* Ripple command line arguments */
        opterr = 0;
        argc -= 1;
        argv += 1;
    }
    else{
        fprintf(stderr, "ERROR: Missing 1st argument. Should be name of QLS file to read.\n");
        usage();
    }
    
    
    /* Process optional command line switches */
    while ((iopt = getopt(argc, argv, "d:Fho:V")) != -1) {
        fprintf(stdout,"After getopt iopt is %c optarg is %s \n",iopt,optarg);
        switch (iopt) {
            case 'd':
                debug = strtol(optarg, 0, 10);
                // a non numeric arg leaves debug == 0 [OFF]
                break;
            case 'o':  //output filename provided
                foutname = (char *) malloc( sizeof(char) * strlen(optarg) + 1);
                memset(foutname, 0x00, sizeof( char ) * strlen(optarg) +1);
                strcpy(foutname, optarg);
                outflag=1;
                if( debug == 1 ) fprintf(stderr, "outfilename %s\n", optarg);
                break;
            case 'F':  //flat
                flatflag = 1;
                break;
            case 'V':
                debug = 1;
                break;
            case 'h':
            case '?':
                
            default:
                fprintf(stderr, "Command line option processing\n");
                usage();
        }
    }
    
    
    /* parse the output filename and open it */
    if ( !outflag ) {  //  if not provided it on the input file name
        sprintf(strtemp, "qls2pha_out.pha");
        foutname= (char *) malloc( sizeof(char) * strlen(strtemp));
        strcpy(foutname, strtemp);
    }
    if( debug == 1 ) fprintf(stderr, "foutname %s\n", foutname);
    
    if ((fpout = fopen(foutname, "w")) == NULL) {
        fprintf(stderr, "Cannot open output phase file %s\n", foutname);
        usage();
    }
    
    /* Read the header line from the quad phase list file */
//     nw = fread(lstdat, sizeof(unsigned short), 6, fpqls);
//     nw = fread(lstdat, sizeof(signed short), 6, fpqls);
//     if (nw != 6) {
//         fprintf(stderr, "Error reading quad phase list file %s\n", flstname);
//         usage();
//     }
//     
//     /* Verify the 'QL' signature and extract the NP (nPatches) NX (nCols) & NY (nRows)*/
//     if (lstdat[0] != 0x4c52 ) {
//         fprintf(stderr, "Invalid GIPhT QuadPhase List signature for file %s\n", flstname);
//         fprintf(stderr, "Expected %x Found %x\n", 0x4c52, lstdat[0]);
//         usage();
//     }
//     np=(long) lstdat[1];
//     nx=(long) lstdat[2];
//     ny=(long) lstdat[3];

        nw = fread(lsthdr, sizeof(signed int), 3, fpqls);
    if (nw != 3) {
        fprintf(stderr, "Error reading quad phase list file %s\n", flstname);
        usage();
    }
    
    np=(long) lsthdr[0];
    nx=(long) lsthdr[1];
    ny=(long) lsthdr[2];

    
    npix=nx*ny;
    if( debug == 1 ) printf("%2s\t%ld\t%ld\t%ld\n", lstdat, np, nx, ny);
    
    /* Allocate and zero an output array for reconstruction */
    //PES pout = calloc (npix,sizeof(signed char));
    pout = calloc(npix, sizeof(signed char));
    if (pout == NULL) {
        fprintf(stderr, "qls2pha3: failure to allocate space for pout buffer\n");
        exit(-1);
    }
    for(i=0;i<npix;i++)pout[i]=(signed char)0;
    
    
    /* Read through the list file generating the image */
    nr=0;	// Number or records counter to check against Npatches (in header)
//     while( (nw = fread(lstdat, sizeof(unsigned short), 6, fpqls) ) == 6) {
    while( (nw = fread(lstdat, sizeof(signed short), 6, fpqls) ) == 6) {
// Code from corresponding part of pha2qls
//         
//                 i2buf[0]=(short) i1;		        	              /* Column Index (X coordinate) of Upper Left Corner   */
//                 i2buf[1]=(short) j1;                               /* Row    Index (Y coordinate) of Upper Left Corner   */
//                 i2buf[2]=(short) i2-i1+1;		                  /* Number of pixels, width and height of square patch */
//                 i2buf[3]=(short) rint(256.0f*(double)cmdr);        /* phase value in 256^2 DN */
//                 i2buf[4]=(short) rint(256.0f*slopevector.slopex);  /* value of phase gradient in +X direction in 256^2 DN == 1 cycle/pixel */
//                 i2buf[5]=(short) rint(256.0f*slopevector.slopey);  /* value of phase gradient in -Y direction in 256^2 DN == 1 cycle/pixel */

        nr++;
        i1    =  (long)lstdat[0];
        j1    =  (long)lstdat[1];
        nwidth=  (long)lstdat[2];
        pv    = ((double)lstdat[3])/256.0/256.0; /* phase value in cycles */
        gx    = ((double)lstdat[4])/256.0/256.0; /* value of phase gradient in cycles per pixel in +X direction */
        gy    = ((double)lstdat[5])/256.0/256.0; /* value of phase gradient in cycles per pixel in -Y direction */
        i2 = i1 + nwidth - 1;
        j2 = j1 + nwidth - 1;
        xmid = (double)(i1 + i2)/2.0;
        ymid = (double)(j1 + j2)/2.0;
        if( debug == 1 ) printf("%5d %5d %5d %#10.4f %#10.4f %#10.4f\n", i1, j1, nwidth, pv, gx, gy);
                
        /* Write pv into an nwidth sized patch with the upper left corner at i1 & j1 */
        if (flatflag == 1 )  // create the output using only the mean data
            for(j=j1;j<=j2;j++) {	// outer loop is rows NY
            for(i=i1;i<=i2;i++) {	// inner loop is cols NX
                k=i+j*nx;		// index
                /* Development code to watch for exceeding array bounds */
                if( k >= npix ){
                    fprintf(stderr, "Error: Exceeding %d array bounds at i1=%d j1=%d nwidth=%d\n", npix, i1, j1, nwidth);
                    exit(-1);
                }
                
                /* debug code to watch for overwriting */
                if( debug == 1 ) {
                    tv = pout[k];
                    if( tv != 0 ) {
                        printf("Warning: Overwriting previous patch at i1=%d j1=%d nwidth=%d\n", i1, j1, nwidth);
                    }
                }
//   2012-OCT-20              pout[k]=(signed char)rint(pv/256.0) ;
                pout[k]=(signed char)rint(pv*256.0) ; /* 2012-OCT-20   coded so that 1 cycle = 256 DN */
                if ( debug == 1 ) fprintf(stdout, "%3d ", pout[k]);
            }
            if( debug == 1 ) fprintf(stdout, "\n");
            }
        else   // create the output with the mean and gradient in each patch
            for(j=j1;j<=j2;j++) {	// outer loop is rows NY
            for(i=i1;i<=i2;i++) {	// inner loop is cols NX
                k=i+j*nx;		// index
                /* Development code to watch for exceeding array bounds */
                if( k >= npix ){
                    fprintf(stdout, "Error: Exceeding %d array bounds at i1=%d j1=%d nwidth=%d\n", npix, i1, j1, nwidth);
                    exit(-1);
                }
                
                /* debug code to watch for overwriting */
                if( debug == 1 ) {
                    tv = pout[k];
                    if( tv != 0 ){
                        printf("Warning: Overwriting previous patch at i1=%d j1=%d nwidth=%d\n", i1, j1, nwidth);
                    }
                }
                r = pv + (double)(i-xmid)*gx + (double)(j-ymid)*gy; /* value in cycles */
                r = 256.0 * r; /* 2011-JUL-18 */
                //r = 256.0 * (r/256.0 - rint(r/256.0));
                if ( debug == 1){
                    if(r > 127 || r < -128) {
                        printf("Warning: Overflow at patch at i1=%d j1=%d nwidth=%d\n", i1, j1, nwidth);
                    }
                }
                pout[k] = (signed char)rint(r); /* value in DN such that 256 DN = 1 cycle */
                if( debug == 1 ) fprintf(stdout, "%3d ", pout[k]);
            }
            if( debug == 1 ) fprintf(stdout, "\n");
            }
    }
    if( debug == 1 ) fprintf(stdout, "\n");
    if( nr != np ) {
        fprintf(stderr, "qls2pha: WARNING: number of lines read from list %d does not match header Npatches %d\n", nr, np); exit(-1);
    }
    nw = fwrite(pout, sizeof(signed char), npix, fpout);
    
    if( nw != npix ) {
        fprintf(stderr, "qls2pha: Error writing output phase data\n"); 
        exit(-1);
    }    
    fprintf(stderr, "qls2pha: Successfully wrote phase data to file named: %s\n", foutname);
    
    /* finish up */
    free(pout);
    fclose(fpout);
    fclose(fpqls);
    exit(0);
}

usage(){
    fprintf(stderr, "usage: qls2pha [-h] [-f] [-d#] [-o output.pha] input.qls \n");
    fprintf(stderr, "\t-h   this help information \n");
    fprintf(stderr, "\t-F   optional Flat output gradients ignored\n");
    fprintf(stderr, "\t-o   optional specify outputfilename\n");
    fprintf(stderr, "\t-V   optional Verbose mode\n");
    fprintf(stderr, "\t-d#  optional debug levels #=[0123]\n");
    fprintf(stderr, "\t\t 0  no debugging (default)\n");
    fprintf(stderr, "\t\t 1  Warning to stdout if a patch write exceeds nx,ny bounds\n");
    fprintf(stderr, "\t\t\t or if a patch write overwrites an existing patch\n");
    fprintf(stderr, "\t\t 2  echos input qls data to stdout\n");
    fprintf(stderr, "\t\t 3  combines 1 & 2\n");
    fprintf(stderr, "\tquadphase list expansion from input.qls\n");
    fprintf(stderr, "\ta signed 8-bit wrapped phase array is written to output.pha\n");
    exit(-1);
}

//lap
//lap     v j increases
//lap     |            jn
//lap     ---------------------------
//lap     |            |              |
//lap     | floor(imid)|floor(imid)+1 |
//lap     |            |              |
//lap     |        floor(jmid)        |
//lap iw ----------------------------- ie
//lap     |        floor(jmid)+1      |
//lap     |            |              |
//lap     |            |              |
//lap     |            |              |
//lap     ---------------------------  x increases ->
//lap                  js
//lap
//lap matlab reconstruction code
//lap i1=ilist(:,1);  % Index to row of first pixel in patch
//lap j1=ilist(:,2);  % Index to col of first pixel in patch
//lap kw=ilist(:,3);  % width (and height) of square patch
//lap qp=ilist(:,4);  % ohase value coded -128 to +127 such that 256 DN = 1 cycle
//lap i2=i1+kw; % patches are square
//lap j2=j1+kw;
//lap figure
//lap hist(double(reshape(qp,numel(qp),1)),64);axis tight;
//lap xlabel('phase value (256 DN per cycle)');
//lap ylabel('Number of phase values');
//lap title('After quadtree partitioning');
//lap
