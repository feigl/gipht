/* PHA2QLS.C
 *
 * perform quad-tree partitioning on an array of phase values
 *
 * Authors:	Kurt Feigl & Peter Sobol
 *
 *
 * gcc pha2qls.c -lm -O3 -o pha2qls
 * gcc pha2qls.c -lm -O3 -o pha2qls.a64
 * gcc pha2qls.c -lm -O3 -o pha2qls.maci
 * gcc pha2qls.c -lm -O3 -o pha2qls.maci64
 * gcc pha2qls.c -lm -O3 -o pha2qls.w64   % 64-bit windows OS
 *
 * lap11jul09 removed all references to COH
 * lap11jul09 Add header with Signature & NX NY for reconstruction
 * lap21jul09 Change header Signature to QL & add Npatches & NX NY
 * lap24jul09 Changed padding to geenerate a square (nxpad == nypad)
 * lap14aug09 misc house cleaning
 * lap21aug09 change error exits report pname from argv[0]
 *
 *  pes 11/12/09 added slope calculation to model.  Each patch is modeled using 3 parameters: circular mean and slope in 2 dimensions
 *  2009-NOV-28 record both slopes
 *  2009-DEC-08 fix bug about file names
 *  2010-JAN-04 record gradients in 256*256 DN per pixel
 *  2010-JAN-11 check floating point math
 *  2010-JUN-21 version 1.1 criterion is misfit to mean, not misfit to ramp
 *  2011-MAR-28 version 1.2
 *  2011/04/20 Version 2.0 PES - updated to add command line options
 *  2011/06/06 Version 2.1 Kurt - error messages
 *  2011/06/15 Version 2.1 Kurt - documentation
 *  2011/06/21 Version 2.1 Kurt - fix bug with center instead of upper right corner
 *  2011/07/19 Version 2.2 Kurt - fix bug with sign of Y-gradient
 *  2011/07/23 Version 2.3 Kurt - add -Q option to specify number of pixels in quadrilateral patch
 *  2011/10/10 Version 2.3 Kurt - allow bigger patches
 *  2012/10/01 Version 2.4 Kurt - improve help
 *  2012/10/04 Version 2.5 Kurt - work on header, -T
 *  2012/10/06 Version 2.5b Kurt - still looking for bug
 *  2014/01/06 Version 2.5c Kurt - still looking for bug
 *
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

static FILE *fpin1   = NULL;
static FILE *fplst   = NULL;
static FILE *fpout   = NULL;
static FILE *fpgrx   = NULL;
static FILE *fpgry   = NULL;
static FILE *fptxt   = NULL;

static char *pname, txtoutflag=0, nanchar=0;
// 20121024 static unsigned short np;             /* number of patches : should be long integer */
/* static long npatches;                 /* number of patches */
/* static long maxpix=1000;              /* maximum number of pixels in a patch       */
/* static long minpix=4;                 /* minimum number of VALID pixels in a patch */
static unsigned short maxpix=1000;              /* maximum number of pixels in a patch       */
static unsigned short minpix=4;                 /* minimum number of VALID pixels in a patch */

static signed char i1thresh=16;       /* max Circular Mean Deviation for misfit to 3-parameter ramp model */
static signed char maxcmd=32;         /* max Circular Mean Deviation for misfit to 3-parameter ramp model */
static signed char  *i1arr, *i1out, *i1wrk, *i1mod;  /* arrays of phase values coded as 1 byte  per pixel */
static signed short *i2out, *i2grx, *i2gry;          /* arrays of phase values coded as 2 bytes per pixel*/
static long nx, ny, npix, npatches=0, nnulltotr=0, nnulltotw=0, nnullpatches=0, nxpad, nypad, npixpad;
static quadtree_node do_create(unsigned int tlx, unsigned int jn, unsigned int ie, unsigned int js, unsigned int level);
/* static double pi=3.141529; Corrected 20140106 */
static double pi=3.141592653589793; /* Corrected 20140106 */

typedef struct { double slopex;   double slopey;} slopexy;

/* this is a wrapper to make the routine recursive */
/* It must occur before the main */
quadtree_node
        quadtree_create(unsigned int iw, unsigned int jn, unsigned int ie, unsigned int js, unsigned int level) {
    return do_create(iw, jn, ie, js, level+1);
}

usage(){
    fprintf(stderr, "%s performs quadtree partitioning on a file of wrapped phase values\n", pname);
    fprintf(stderr, "usage pha2qls [options] input.pha nx ny \n");
    fprintf(stderr, "Required arguments: \n");
    fprintf(stderr, "\tinput.pha == Name of input binary file containing phase values coded as one signed byte per such that 256 DN = 1 cycle and values range from -128 to +127\n");
    fprintf(stderr, "\tnx        == Number of pixels in X-direction [= number of columns]\n");
    fprintf(stderr, "\tny        == Number of pixels in Y-direction [= number of rows]\n");
    fprintf(stderr, "Output: \n");
    fprintf(stderr, "\tFile containing Quadtree List of Samples (QLS) 6-column binary file containing I,J,W,Q,GX,GY [named output.qls by default]\n");
    fprintf(stderr, "Options: \n");
    fprintf(stderr, "\t-h   this help information \n");
    fprintf(stderr, "\t-L l == threshold Level   for circular mean deviation misfit to 1-parameter mean [256 DN = 1 cycle of phase] \n");
    fprintf(stderr, "\t-M m == threshold Maximum for circular mean deviation misfit to 3-parameter ramp [256 DN = 1 cycle of phase] \n");
    fprintf(stderr, "\t-N n == mininum number of pixels in patch \n");
    fprintf(stderr, "\t-o output_qls_filename.pha   == file name for output reconstructed phase values, written as an image with nx columns and ny rows \n");
    fprintf(stderr, "\t-P output_phase_filename.pha == file name for output reconstructed phase values, written as an image with nx columns and ny rows.\n");
    fprintf(stderr, "\t-X output_gradx_filename.i2  == file name for differential phase in positive X-direction (increasing column index) [2^16 = 65536 DN = 1 cycle of phase per pixel]\n");
    fprintf(stderr, "\t-Y output_gradx_filename.i2  == file name for differential phase in negative Y-direction (increasing line index)   [2^16 = 65536 DN = 1 cycle of phase per pixel]\n");
    fprintf(stderr, "\t-V Verbose mode\n");
    fprintf(stderr, "\t-d n  optional debug levels\n");
    fprintf(stderr, "\t\t -d 0  == No additional output debugging information to stderr (default)\n");
    fprintf(stderr, "\t\t -d 1  == Warning to stdout if a patch write exceeds nx,ny bounds or if a patch write overwrites an existing patch\n");
    fprintf(stderr, "\t\t -d 2  == echos input qls data to stdout\n");
    fprintf(stderr, "\t\t -d 3  == combines 1 & 2\n");
    exit(-1);
}

main(argc, argv)
int argc;
char **argv; {
    long nw = 0; /* number of words written */
    long nr = 0; /* number of words read    */
    long maxpatch;
    signed char c, c1;
    signed short i2a;
    /*20121004    signed short lsthdr[6]; /* header record */
    signed int lsthdr[3]; /* header record */
    int iopt,ierr;
    /* file names */
    char *fphaname, *flstname, *fgrxname, *fgryname, *fqspname , *fqlsname, *ftxtname;
    char outflag = 0, debug = 0, verboseflag = 1, xgradflag=0, ygradflag=0, qspoutflag=0, qlsapplyflag=0, fillflag=0;
    double change;
    unsigned int i,j;
    unsigned int iw, jn, ie, js;
    unsigned int level;
    quadtree_node top;
    void heapSort();
    char *strprefix();
    char *strextension();
    char *strreplace();
    
    pi = 4.0 * atan(1.0);
    pname = argv[0];
    
    fprintf(stdout, "pha2qls version 2.5c of 2014-JAN-06 \n");
    fprintf(stderr, "Copyright (c) 2011 Kurt L. Feigl & Peter E. Sobol\n");
    fprintf(stderr, "University of Wisconsin-Madison\n");
    fprintf(stderr, "All rights reserved\n\n");
    fprintf(stderr, "Part of General Inversion of Phase Technique (GIPhT)\n\n");
    fprintf(stderr, "Copyright (c) Kurt Feigl, Cliff Thurber  All rights reserved.\n");
    fprintf(stderr, "U.S. Patent No. 7,446,705.\n\n");
    fprintf(stderr, "Use of this computer program [GIPhT] is governed by a license agreement\n");
    fprintf(stderr, "between You and the Wisconsin Alumni Research Foundation (the\n");
    fprintf(stderr, "'Agreement'). Unauthorized use, reproduction, or distribution of this\n");
    fprintf(stderr, "program or any portion of it may result in severe civil and criminal\n");
    fprintf(stderr, "penalties, and will be prosecuted to the maximum extent possible under\n");
    fprintf(stderr, "the law. By using this program, you agree to be bound by the terms of the\n");
    fprintf(stderr, "Agreement.\n\n");
    
    /* check for machine dependence before going any farther */
    if (sizeof(signed int) != 4) {
        fprintf(stdout, "ERROR: sizeof(signed int) is %ld bytes on this machine. Expected 4 bytes.\n",sizeof(signed int));
        exit(-1);
    }
    if (sizeof(signed short) != 2) {
        fprintf(stdout, "ERROR: sizeof(signed short) is %ld bytes on this machine. Expected 2 bytes.\n",sizeof(signed short));
        exit(-1);
    }
    if (sizeof(signed char) != 1) {
        fprintf(stdout, "ERROR: sizeof(signed char) is %ld bytes on this machine. Expected 1 byte.\n",sizeof(signed char));
        exit(-1);
    }
    
    /* first argument is the phase file .pha */
    /* test & open the input filename */
    if (argc > 1) {
        fphaname = (char *) malloc(sizeof(char) * strlen(argv[1])+1);
        strcpy(fphaname, argv[1]);
        if ((fpin1 = fopen(fphaname, "r")) == NULL) {
            fprintf(stderr, "Cannot open input phase file %s\n", fphaname);
            usage();
        }
        else {
            if( verboseflag )  fprintf(stderr, "Succesfully opened input phase file %s\n", fphaname);
        }
    }
    else{
        fprintf(stderr, "Missing input .pha file argument\n");
        usage();
    }
    
    /* read the size arguments from the command line */
    if (argc > 2) {
        if ((nx = atoi(argv[2])) == 0) {
            fprintf(stderr, "%s: Cannot read NX %s\n", pname, argv[2]);
            exit(-1);
        }
        else {
            if( verboseflag ) fprintf(stderr, "nx=N(cols)     = %10ld\n", nx);
        }
        
        if ((ny = atoi(argv[3])) == 0) {
            fprintf(stderr, "%s: Cannot read NY %s\n", pname, argv[3]);
            exit(-1);
        }
        else {
            if( verboseflag ) fprintf(stderr, "ny=N(rows)     = %10ld\n", ny);
        }
    }
    else{
        fprintf(stderr, "Missing X or Y size argument\n");
        usage();
    }
    
    /* Proccess remaining options on command line */
    /* opterr = 0; */
    argc -= 3;
    argv += 3;
    
    /* switches that take an argument are followed by a colon : in the list below */
    while ((iopt = getopt(argc, argv, "d:o:FVT:N:L:n:M:X:Y:H:A:P:Q:")) != -1) {
/*      fprintf(stdout,"After getopt iopt is %c optarg is %s \n",iopt,optarg);*/
        switch (iopt) {
            case 'd':
                debug = strtol(optarg, 0, 10);
                /* a non numeric arg leaves debug == 0 [OFF] */
                break;
            case 'o':  /*output filename provided */
                flstname = (char *) malloc( sizeof(char) * strlen(optarg) + 1);
                memset(flstname, 0x00,      sizeof(char) * strlen(optarg) + 1);
                strcpy(flstname, optarg);
                outflag=1;
                break;
            case 'n':  /*Nan character specified */
                nanchar=atoi(optarg);
                break;
            case 'F': /* fill bad patches with original values */                
                fillflag=1; 
                break;
            case 'X': /*output X gradient */
                fgrxname = (char *) malloc( sizeof(char) * strlen(optarg) + 1);
                memset(fgrxname, 0x00,      sizeof(char) * strlen(optarg) + 1);
                strcpy(fgrxname, optarg);
                xgradflag=1;
                break;
            case 'Y': /*output Y gradient */
                fgryname = (char *) malloc( sizeof(char) * strlen(optarg) + 1);
                memset(fgryname, 0x00,      sizeof(char) * strlen(optarg) + 1);
                strcpy(fgryname, optarg);
                ygradflag=1;
                break;
            case 'L':  /*threshold specified */
                i1thresh =(signed char) atoi(optarg);
                break;
            case 'M': /*maxcmd specified */
                maxcmd=(signed char) atoi(optarg);
                break;
            case 'N':  /*minpix specified */
                minpix=atol(optarg);
                break;
            case 'T':
                ftxtname = (char *) malloc( sizeof(char) * strlen(optarg) + 1);
                memset(ftxtname, 0x00,      sizeof(char) * strlen(optarg) + 1);
                strcpy(ftxtname, optarg);
                ygradflag=1;
                txtoutflag=1;
                break ;
            case 'V':
                /* verboseflag=atoi(optarg); BAD: variable is char, function returns int */
                verboseflag=1; /* just make it a toggle */
                break;
            case 'P': /*output quadtree phase file */
                fqspname=(char*)malloc(sizeof(char)*strlen(optarg)+1);
                memset(fqspname, 0x00, sizeof(char)*strlen(optarg)+1);
                strcpy(fqspname, optarg);
                qspoutflag=1;
                break;
            case 'Q':  /* number of pixels in quadrilateral patch */
                maxpix=(unsigned short) atoi(optarg);
                break;
            case 'A': /*apply input QLS file to current pha file */
                fqlsname=(char*)malloc(sizeof(char)*strlen(optarg)+1);
                memset(fqlsname, 0x00, sizeof(char)*strlen(optarg)+1);
                strcpy(fqlsname, optarg);
                qlsapplyflag=1;
                break;
            case 'h':
            case '?':
            default:
                fprintf(stderr, "Command line option processing\n");
                usage();
        }
    }
    
    if (verboseflag) fprintf(stdout, "i1thresh, maxcmd, minpix: -L %d, -M %d, -N %d\n", i1thresh, maxcmd, minpix);
    
    /*open the output file */
    if ( !outflag ) {  /*  if not provided it on the input file name */
        flstname=strextension(fphaname, "qls"); /* replace suffix only */
    }
    
    if (qspoutflag){  /* setup for quadtree phase file output */
        if ((fpout = fopen(fqspname, "w")) == NULL) {
            fprintf(stderr, "Cannot open output phase file %s\n", fqspname);
            usage();
            exit(-1);
        }
        else {
            if (verboseflag) fprintf(stdout, "Opened: %s for output\n", fqspname);
        }
    }
    
    if (txtoutflag){  /* setup for quadtree list file text output */
        if (strlen(ftxtname) < 1){
            ftxtname=strextension(fphaname, "qls.txt");
        }
        if ((fptxt = fopen(ftxtname, "w")) == NULL) {
            fprintf(stderr, "Cannot open output text file for list %s\n", ftxtname);
            usage();
            exit(-1);
        }
        else {
            if (verboseflag) fprintf(stdout, "Opened: %s for output\n", ftxtname);
        }
    }
    
    if (xgradflag) {
        if ((fpgrx = fopen(fgrxname, "w")) == NULL) {
            fprintf(stderr, "%s: Cannot open output phase X-gradient file %s\n", pname, fgrxname);
            exit(-1);
        }
        else {
            if (verboseflag) fprintf(stdout, "Opened: %s\n", fgrxname);
        }
    }
    
    if (ygradflag) {
        if ((fpgry = fopen(fgryname, "w")) == NULL) {
            fprintf(stderr, "%s: Cannot open output phase Y-gradient file %s\n", pname, fgryname);
            exit(-1);
        }
        else {
            if (verboseflag) fprintf(stdout, "Opened: %s\n", fgryname);
        }
    }
    
    /* open output list file */
    if ((fplst = fopen(flstname, "w")) == NULL) {
        fprintf(stderr, "%s: Cannot open output data file with list of phase values%s\n", pname, flstname);
        exit(-1);
    }
    else {
        if (verboseflag) fprintf(stdout, "Opened for QLS output: %s\n", flstname);
    }
    
    fprintf(stderr, "Maximum misfit circular mean deviation to boxcar (-L) = %10d\n", i1thresh);
    fprintf(stderr, "Maximum misfit circular mean deviation to ramp   (-M) = %10d\n", maxcmd);
    fprintf(stderr, "Mininum Number of pixels in patch                (-N) = %10d\n", minpix);
    
    /* allocate memory */
    nxpad = (long)pow(2, ceil(log((float)nx)/log(2.0)));
    nypad = (long)pow(2, ceil(log((float)ny)/log(2.0)));
    if(nxpad>nypad) nypad=nxpad;
    if(nypad>nxpad) nxpad=nypad;
    npix = nx*ny;
    npixpad = nxpad*nypad;
    
    i1arr = calloc(npixpad, sizeof(signed char)); if (i1arr == NULL) {fprintf(stderr, "%s: failure to allocate space for i1arr buffer\n", pname); exit(-1);}
    i1out = calloc(npixpad, sizeof(signed char)); if (i1out == NULL) {fprintf(stderr, "%s: failure to allocate space for i1out buffer\n", pname); exit(-1);}
    i1wrk = calloc(npixpad, sizeof(signed char)); if (i1wrk == NULL) {fprintf(stderr, "%s: failure to allocate space for i1wrk buffer\n", pname); exit(-1);}
    i1mod = calloc(npixpad, sizeof(signed char)); if (i1mod == NULL) {fprintf(stderr, "%s: failure to allocate space for i1mod buffer\n", pname); exit(-1);}
    i2out = calloc(npixpad, sizeof(signed short));if (i2out == NULL) {fprintf(stderr, "%s: failure to allocate space for i2out buffer\n", pname); exit(-1);}
    i2grx = calloc(npixpad, sizeof(signed short));if (i2grx == NULL) {fprintf(stderr, "%s: failure to allocate space for i2grx buffer\n", pname); exit(-1);}
    i2gry = calloc(npixpad, sizeof(signed short));if (i2gry == NULL) {fprintf(stderr, "%s: failure to allocate space for i2gry buffer\n", pname); exit(-1);}
    
    maxpix = 1000;  /* max number of pixels in patch is large by default */
    fprintf(stderr, "Maximum number of pixels in patch  = %10d\n", maxpix);
    
    if (minpix > maxpix){
        fprintf(stderr, "ERROR: Minimum number of pixels (%d) in a patch is greater than maximum (%d)\n", minpix,maxpix);
        usage();
        exit(-1);
    }
    
    
    /* Start the output phase list with NP NX NY where
     * Where NP is # patches, NX & NY are the number of rows and number of col
     * /* 2012-OCT-04 */
    lsthdr[0]= (signed int) 0;	/* needs to be filled in later */
    lsthdr[1]= (signed int) nx;
    lsthdr[2]= (signed int) ny;
    /* write as a placeholder, rewite when Npatches is known */
    nw = fwrite(lsthdr, sizeof(signed int), 3, fplst);
    if (nw != 3) {
        fprintf(stderr, "%s: Error writing LstHdr to qls file!\n", flstname);
        exit(-1);
    }
    
    
    /* write header to text list */
    if (txtoutflag)  {
        ierr = fseek(fptxt, 0L, SEEK_SET);
        if (ierr != 0) {
            fprintf(stderr, "%s: ERROR: rewinding text list. IERR is %d", pname, ierr);
            exit(-1);
        }
        nw=fprintf(fptxt, "--I1-- --J1-- NWIDE- -CMDR- GRADX- GRADY- ! NX = %6ld NY = %6ld \n"
                ,nx,ny);
        if (nw < 3) {
            fprintf(stderr, "%s: Error writing LstHdr to txt file.\n nw is %ld", pname, nw);
            exit(-1);
        }
    }
    
    /* read the whole phase array */
    nr = fread(i1wrk, sizeof(signed char), npix, fpin1);
    fprintf(stderr, "N(pixels in)   = %10ld\n", nr);
    
    if (nr == npix) {
        /* copy from wrk to arr, counting nulls */
        nnulltotr = 0;
        for (j = 0; j < ny; j++){
            for (i = 0; i < nx; i++){
                c = i1wrk[i+j*nx];
                i1arr[i+j*nxpad] = c;
                if (c == 0) nnulltotr++;
            }
        }
        
       
        fprintf(stderr, "N(null in)     = %10ld\n", nnulltotr);
        
        /* pad right side with random values */
        fprintf(stderr, "Padding ends of rows from %ld to %ld\n", nx, nxpad);
        for (j = 0; j < nypad; j++){
            for (i = nx; i < nxpad; i++){
                i1arr[i+j*nxpad] = (signed char)rand() ;
            }
        }
        /* pad bottom side with random values */
        fprintf(stderr, "Padding ends of cols from %ld to %ld\n", ny, nypad);
        for (j = ny; j < nypad; j++){
            for (i = 0; i < nxpad; i++){
                i1arr[i+j*nxpad] =(signed char)rand();
            }
        }
        /* initialize output phase */
        for (j = 0; j < nypad; j++){
            for (i = 0; i < nxpad; i++){
                i1out[i+j*nxpad] = (signed char) 0;
            }
        }  
        
        /* pass through original values if requested */
        if (fillflag == 1) {
            for (j = 0; j < ny; j++){
                for (i = 0; i < nx; i++){
                    c = i1wrk[i+j*nx];
                    i1out[i+j*nxpad]=c;
                }
            }
        }
 
        /* initialize output gradients */
        for (j = 0; j < nypad; j++){
            for (i = 0; i < nxpad; i++){
                i2out[i+j*nxpad] = (signed short) 0;
            }
        }
        /* initialize output phase gradient in eastward (X) direction */
        for (j = 0; j < nypad; j++){
            for (i = 0; i < nxpad; i++){
                i2grx[i+j*nxpad] = (signed short) 0;
            }
        }
        /* initialize output phase gradient in northward (Y) direction */
        for (j = 0; j < nypad; j++){
            for (i = 0; i < nxpad; i++){
                i2gry[i+j*nxpad] = (signed short) 0;
            }
        }
        
        /* Count like a C programmer */
        iw = 0;
        ie = nxpad-1;
        jn = 0;
        js = nypad-1;
        
        /* level = (long) ceil(log(npix)/log(4))-3; */
        /* if (level < 3) level = 3; */
        /* level = 15; */ /* this is the maximum */
        level = 3;
        if (level < (long) ceil(log(nxpad)/log(2))) level = (long) ceil(log(nx)/log(2));
        if (level < (long) ceil(log(nypad)/log(2))) level = (long) ceil(log(ny)/log(2));
        if (level > 15) level = 15;
        fprintf(stderr, "N(levels)      = %10d\n", level);
        maxpatch = (long) pow(4, level);
        fprintf(stderr, "N(maxpatch)    = %10ld\n", maxpatch);
        
        if (maxpatch < npix/minpix) {
            fprintf(stderr, "WARNING: number of levels (%d) only allows %ld patches, but image contains %ld pixels.\n", level, maxpatch, npix);
        }
        top=do_create(iw, jn, ie, js, level);
    }
    else {
        fprintf(stderr, "\n%s: Wrong number of pixels in file %80s expected: %ld got: %ld\n", pname, fphaname, npix, nr); exit(-1);
    }
    fclose(fpin1);
    
    /* copy from grx to out, truncating lines and columns, and counting nulls */
    nnulltotw = 0;
    for (j = 0; j < ny; j++){
        for (i = 0; i < nx; i++){
            i2a = i2grx[i+j*nxpad];
            i2out[i+j*nx]=i2a;
            if (i2a == 0) nnulltotw++;  /* count nulls */
        }
    }
    fprintf(stderr, "N(null grx)    = %10ld\n", nnulltotw);
    
    /* Write the whole array */
    if (xgradflag){
        nw = fwrite(i2out, sizeof(signed short), npix, fpgrx);
        if (verboseflag) fprintf(stdout, "Writing gradient-X file %s\n", fgrxname);
        if (nw != npix) {
            fprintf(stderr, "ERROR: Nw      = %12ld\n", nw);
            exit(-1);
        }
    }
    
    /* copy from gry to out, truncating lines and columns, and counting nulls */
    nnulltotw = 0;
    for (j = 0; j < ny; j++){
        for (i = 0; i < nx; i++){
            i2a = i2gry[i+j*nxpad];
            i2out[i+j*nx]=i2a;
            if (i2a == 0) nnulltotw++;  /* count nulls */
        }
    }
    fprintf(stderr, "N(null gry)    = %10ld\n", nnulltotw);
    
    /* Write the whole array */
    if (ygradflag) {
        nw = fwrite(i2out, sizeof(signed short), npix, fpgry);
        if (verboseflag) fprintf(stdout, "writing gradient-y file %s\n", fgryname);
        
        if (nw != npix) {
            fprintf(stderr, "ERROR: Nw      = %12ld\n", nw);
            exit(-1);
        }
    }
    
    /* copy from output to work, truncating lines and columns, and counting nulls */
    nnulltotw = 0;
    for (j = 0; j < ny; j++){
        for (i = 0; i < nx; i++){
            c = i1out[i+j*nxpad];
            i1wrk[i+j*nx]=c;
            if (c == 0) nnulltotw++;  /* count nulls */
        }
    }
    fprintf(stderr, "N(null out)    = %10ld\n", nnulltotw);
    
    if (qspoutflag){
        /* Write the whole array */
        if (verboseflag) fprintf(stdout, "writing quadtree phase file %s\n", fqspname);
        nw = fwrite(i1wrk, sizeof(signed char), npix, fpout);
        if (nw != npix) {
            fprintf(stderr, "ERROR: Nw      = %12ld\n", nw);
            exit(-1);
        }
    }
    else {
        nw = nx * ny;
    }
    
    change = 6*sizeof(signed short)*(double)(npatches-nnullpatches)/(sizeof(signed char)*(double)npix);
    fprintf(stderr, "N(null patches)= %10ld\n", nnullpatches);
    fprintf(stderr, "N(patches)     = %10ld\n", npatches);
    fprintf(stderr, "N(OK patches)  = %10ld\n", npatches-nnullpatches);
    fprintf(stderr, "N(pixels out)  = %10ld\n", npix);
    fprintf(stderr, "Change         = %#10.1f percent\n", 100.0f * change);
    fprintf(stderr, "Compression    = %#10.1f percent\n", 100.0f * change);
    fprintf(stderr, "Reduction      = %#10.1f percent\n", 100.0f * (1.0f - change));
    
    /* write header to binary list */
    /* 20121004    lsthdr[1]=np;	/* rewrite Npatches */
    /* 20121004    nw = fwrite(lsthdr, sizeof(signed short), 6, fplst); */
    /* 2012-OCT-04 */
    lsthdr[0]= (signed int) npatches;
    lsthdr[1]= (signed int) nx;
    lsthdr[2]= (signed int) ny;
    ierr=fseek(fplst, 0L, SEEK_SET);
    if (ierr != 0) {
        fprintf(stderr, "%s: ERROR: rewinding QLS list. IERR is %d", pname, ierr);
        exit(-1);
    }
    nw = fwrite(lsthdr, sizeof(signed int), 3, fplst);
    if (nw != 3) {
        fprintf(stderr, "%s: Error rewriting LstHdr to qls file!\n", pname); exit(-1);
    }
    
    /* finish up */
    free(i1arr);
    free(i1out);
    free(i1wrk);
    free(i1mod);
    free(i2out);
    free(i2grx);
    free(i2gry);
    
    fclose(fplst);
    if (qspoutflag) fclose(fpout);
    if (txtoutflag) fclose(fptxt);
    if (xgradflag)  fclose(fpgrx);
    if (ygradflag)  fclose(fpgry);
    
    exit(0);
}

/* Really create the quadtree */
static quadtree_node
        do_create(unsigned int iw, unsigned int jn, unsigned int ie, unsigned int js, unsigned int level) {
    quadtree_node root;
    unsigned int ia, ib, ja, jb;
    float fimid, fjmid;
    long nok, nnull, nw=0, nrows, ncols, npixinpatch;
    unsigned int i, j, i1, i2, j1, j2;
    signed char cdev0, cmdr, c1;
    signed char cmeandir(signed char *ic1, long n), cdev1;
    signed char circmeandev(signed char *ic1, signed char alpha, long n);
    signed char circslopemeandev(signed char *ic1, signed char *i1mod, long nrows, long ncols);
    signed char *applyramp(signed char *ic1, signed char *i1mod, signed char cmean, slopexy slopevector, long nrows, long ncols);
    signed short i2a, i2b, i2c;
    
    void heapSort();
    signed short i2buf[6];
    slopexy slopevector,slopevector1,slopevector2,slope1(),slope2();
    
    double d2; /* temp */
    /* index values for columns, counting from left  */
    i1 = iw;  /* west edge */
    if (i1 < 0) i1 = 0;
    i2 = ie;  /* east edge */
    if (i2 > nxpad-1) i2 = nxpad-1;
    /* index values for rows, counting from top */
    j1 = jn;  /* north edge */
    if (j1 < 0) j1 = 0;
    j2 = js;  /* south edge */
    if (j2 > nypad-1) j2 = nypad-1;
    
    /* dimensions */
    nrows=i2-i1+1;
    ncols=j2-j1+1;
/*     npixinpatch = nrows*ncols; */
    
     /* copy phase values from this patch into working array for statistics include zero (null values), count and count nulls 11/2/09 Peter t*/
    nok = 0;nnull=0;
    npixinpatch = 0;
    for (j = j1; j <= j2; j++){
        for (i = i1; i <= i2; i++){
/*             i1wrk[nok] = i1arr[i+j*nxpad]; */
            c1 = i1arr[i+j*nxpad];
            i1wrk[npixinpatch] = c1;
            npixinpatch++;
            /* nok++; KLF 2009 NOV 14 */
            if (c1==0){
                nnull++;
            }
            else {
                nok++;
            }
        }
    }
    if (npixinpatch != nrows*ncols){
        fprintf(stderr, "ERROR miscount: %5d %5d %5d %5d %5ld\n", i1, i2, j1, j2, npixinpatch);
        exit(-1);
    }
    
    
    /*estimate slopes */
    slopevector=slope2(i1wrk, nrows, ncols); /* data,nrows,ncols
    
    /* calculate mean direction  */
    cmdr = cmeandir(i1wrk, npixinpatch);
    
    /* calculate circular mean deviation from average */
    cdev0 = circmeandev(i1wrk, cmdr, npixinpatch);
    
    /* calculate modeled value for all pixels in patch */
    i1mod = applyramp(i1wrk, i1mod, cmdr, slopevector, nrows, ncols);
    
    /* calculate circular mean deviation from ramp */
    cdev1 = circslopemeandev(i1wrk, i1mod, nrows, ncols);
    
    /* stop dividing if criteria are met */
    /*     if (cdev0 <= i1thresh && cdev1 <= maxcmd) { */
    /*     if (cdev0 <= i1thresh && cdev1 <= maxcmd && (nok+nnull) <= maxpix && nok >= minpix) { */
    /*    if (cdev0 <= i1thresh && cdev1 <= maxcmd && npixinpatch <= maxpix) { */
    /*  2012-10-01 Kurt */
/*     if (cdev0 <= i1thresh && cdev1 <= maxcmd && npixinpatch <= maxpix && nok >= minpix) { */
    if (cdev0 <= i1thresh && cdev1 <= maxcmd && npixinpatch <= maxpix && nok >= minpix) {
        level = 0;
        /* 20121004 npatches++; */
        /*         if (nok < minpix || nnull > nok || cdev1 > maxcmd) { /* 2011 JUN 06 Kurt */ 
        /*        if (nok < minpix || nnull > nok) { /* 2011 JUN 06 Kurt */ 
/*         if (nok < minpix) {                  /* 2012 OCT 01 Kurt */ 
/*             cmdr = 0; */
/*             nnullpatches++; */
/*         } */
        
/*         /* write the non-zero model values to the list  */ 
/* 2012-OCT-06         if (i1 < nx && i2 < nx && j1 < ny && j2 < ny && cmdr != 0) { */
/*       if (i1 < nx && i2 < nx && j1 < ny && j2 < ny) { */
        /* check if indices are in bounds */ 
        if (i1 >= 0 && i2 < nxpad && j1 >= 0 && j2 < nypad) {
            /* do not write indices in padded margins to list */
            if (i1 < nx && i2 < nx && j1 < ny && j2 < ny) {
                if (i2-i1 == j2-j1) {
                    i2buf[0]=(signed short) i1;                               /* Column Index (X coordinate) of Upper Left Corner   */
                    i2buf[1]=(signed short) j1;                               /* Row    Index (Y coordinate) of Upper Left Corner   */
                    i2buf[2]=(signed short) i2-i1+1;		                  /* Number of pixels, width and height of square patch */
                    i2buf[3]=(signed short) rint(256.0f*(double)cmdr);        /* mean phase value in 256^2 DN */
                    i2buf[4]=(signed short) rint(256.0f*slopevector.slopex);  /* value of phase gradient in +X direction in 256^2 DN == 1 cycle/pixel */
                    i2buf[5]=(signed short) rint(256.0f*slopevector.slopey);  /* value of phase gradient in -Y direction in 256^2 DN == 1 cycle/pixel */
                }
                else {
                    fprintf(stderr, "ERROR patch is not square! %5d %5d %5d %5d (%5d)\n", i1, i2, j1, j2, cmdr);
                    exit(-1);
                }
                
                nw = fwrite(i2buf, sizeof(unsigned short), 6, fplst);
                if (nw != 6) {
                    fprintf(stderr, "%s: Write error on qls file!\n", pname); exit(-1);
                }
                if (txtoutflag) nw = fprintf(fptxt,"%6d %6d %6d %6d %6d %6d\n",i2buf[0],i2buf[1],i2buf[2],i2buf[3],i2buf[4],i2buf[5]);
                
                /* 20121004 np++;	/* Npatches for header */
                npatches++; /* Npatches for header */
                
                /* write modeled value to all pixels in patch */
                for (j = j1; j <= j2; j++){
                    for (i = i1; i <= i2; i++){
                        i1out[i+j*nxpad] = i1mod[(i-i1) + (j-j1)*ncols];    /* modeled value from 3-parameter ramp model */
                    }
                }
                /* write X gradient value to all pixels in patch */
                for (j = j1; j <= j2; j++){
                    for (i = i1; i <= i2; i++){
                        d2 = rint(slopevector.slopex*256.0); /* coded on interval 256*[-128, +127] */
                        if (fabs(d2) > 0) {
                            i2a = (signed short) d2;
                            /*fprintf(stderr,"%5d %5d %12.4e %12.4e %12.4e %5d\n",i,j,d,d2,i2a); */
                        }
                        else {
                            i2a = 0;
                        }
                        i2grx[i+j*nxpad] = i2a;
                    }
                }
                /* write Y gradient value to all pixels in patch */
                for (j = j1; j <= j2; j++){
                    for (i = i1; i <= i2; i++){
                        d2 = rint(slopevector.slopey*256.0); /* coded on interval 256*[-128, +127] */
                        if (fabs(d2) > 0) {
                            i2a = (signed short) d2;
                            /*fprintf(stderr,"%5d %5d %12.4e %12.4e %12.4e %5d\n",i,j,d,d2,i2a); */
                        }
                        else {
                            i2a = 0;
                        }
                        i2gry[i+j*nxpad] = i2a;
                    }
                }
            }
        }
        else {
            fprintf(stderr, "ERROR index is out of bounds! %5d %5d %5d %5d (%5d)\n", i1, i2, j1, j2, cmdr);
            exit(-1);
        }
    }
    if(level == 0) {
        /*   fprintf(stdout,"level=0, return\n"); */
        return NULL;
    }
    else {  /* deviation is not below threshold, subdivide the patch */
        
        /*fprintf(stdout,"Subdividing patch\n"); */
        fimid = ((float)ie+(float)iw)/2.0f;
        fjmid = ((float)js+(float)jn)/2.0f;
        
        root = calloc(1, sizeof(struct s_quadtree_node));
        bbox2d_set(&root->box, iw, jn, ie, js);
        
        
        /*
         *
         * v j increases
         * |            jn
         * ---------------------------
         * |            |              |
         * | floor(imid)|floor(imid)+1 |
         * |            |              |
         * |        floor(jmid)        |
         * iw ----------------------------- ie
         * |        floor(jmid)+1      |
         * |            |              |
         * |            |              |
         * |            |              |
         * ---------------------------  x increases ->
         * js
         *
         */
        ia = (unsigned int) floor(fimid);
        ib = (unsigned int)  ceil(fimid);
        ja = (unsigned int) floor(fjmid);
        jb = (unsigned int)  ceil(fjmid);
        
        root->childs[0] = do_create(  iw,   jn,  ia,   ja,  level-1);  /* Top Left     Quadrant */
        root->childs[1] = do_create(  ib,   jn,  ie,   ja,  level-1);  /* Top Right    Quadrant */
        root->childs[2] = do_create(  iw,   jb,  ia,   js,  level-1);  /* Bottom Left  Quadrant */
        root->childs[3] = do_create(  ib,   jb,  ie,   js,  level-1);  /* Bottom Right Quadrant */
        return root;
    }
}


signed char cmeandir(signed char *ic1, long n)
/* mean direction, as defined by Mardia and Jupp [2000], page 15 */
/* input phase values with 256 DN = 1 cycle, ranging from -128 to +127 */
/* Kurt Feigl 2009-06-14 */
/* modified to work with slope calc - now nulls are in input, need to ignore them. */
/* inputs: ic1:patch data, n:number of points in patch, nnull:number of nulls in patch */
/* Peter Sobol 2009-11-13 */
{
    double c=0.0f, s=0.0f, r, t;
    signed char ic3;
    long i;
    
    /*knock out the zeros for the mean calculation if there are any (nnull) */
    
    for (i=0;i < n; i++){
        if (ic1[i]!=0){
            c = c + cos( pi * (float) ic1[i]/128.0f);
            s = s + sin( pi * (float) ic1[i]/128.0f);
        }
    }
    r = pow(c, 2)+pow(s, 2); /* no need to take the square root, we are only test r>0. */
    if (r > 0.0f) {
        if (c >= 0.0f) {
            t = atan(s/c);
        }
        else {
            t = atan(s/c) + pi;
        }
        return (signed char) rint(128.0f*t/pi);
    }
    else {
        return (signed char) 0;
    }
}

slopexy slope1(signed char *ic1, long nrows, long ncols)
/* mean slope for each direction for a square input matrix of nrows*ncols
 * /* If one or fewer elements return 0 for slope
 * /* Peter Sobol 2009-11-01 */
{
    double sx=0.0f, sy=0.0f, sum=0.0f;
    signed char cdif;
    long i, j, ncrtemp, kount=0;
    slopexy sxy;
    
    if (nrows >0 && ncols >0) {
        /* estimate slope in X-direction, along rows */
        sum=0;kount=0;
        for (j=0; j < nrows; j++) {
            ncrtemp=ncols*j;
            for (i=1; i < ncols; i++) {
                if (ic1[i+ncrtemp] !=0){
                    cdif = ic1[i+ncrtemp]-ic1[(i-1)+ncrtemp]; /* preserve modulo math */
                    sum = sum + (double) cdif;
                    kount++;
                }
                else {
                    i++; /*if zero skip the next slope calc too. */
                }
            }
        }
        if (kount > 0){
            sx = sum / (double) kount;
        }
        else {
            sx = 0.0;
        }
        
        /* estimate slope in Y-direction, along columns */
        sum=0;kount=0;
        for (j=0; j < ncols; j++) {
            for (i=1; i < nrows; i++) {
                if (ic1[i*ncols+j]!=0){
                    cdif = ic1[i*ncols+j]-ic1[(i-1)*ncols+j]; /* preserve modulo math */
                    sum = sum + (double) cdif;
                    kount++;
                }
                else {
                    i++; /* if zero skip the next slope calce too. */
                }
            }
        }
        
        if (kount > 0){
            sy= sum / (double) kount;
        }
        else {
            sy = 0.0;
        }
        
        /* slope must be less than half a cycle per pixel */
        /*
         * if (abs(sy) > 127.0) {
         * sy = 0.0;
         * }
         */
    }
    else {
        sx = 0.0;
        sy = 0.0;
    }
    
    sxy.slopex=sx;
    sxy.slopey=sy;
    return (slopexy) sxy;
}


slopexy slope2(signed char *ic1, long nrows, long ncols)
/* mean slope for each direction for a square input matrix of nrows*ncols
 *  If one or fewer elements return 0 for slope
 *  Peter Sobol 2009-11-01
 *  2011-06-21 Kurt Feigl modified for center. If too many nulls, return 0
 */
{
    
    double sc=0.0f, sr=0.0f, s1=0.0f, s2=0.0f;
    signed char cdif, c1, c2;
    long i, j, ncrtemp, kount=0, nnull;
    
    slopexy sxy;
    
    if (nrows > 0 && ncols > 0) {
        /* estimate slope in direction of increasing row number, i.e. decreasing Y */
        s1=0.0;s2=0.0f;kount=0;nnull=0;
        for (i=0; i < ncols; i++) {
            for (j=1; j < nrows; j++) {
/*                 2011-JUL-18: The next two lines are wrong and took only a couple of weeks to find!! */
/*                 c2 = ic1[(j-1)*ncols+i]; */
/*                 c1 = ic1[   j *ncols+i]; */
                c1 = ic1[(j-1)*ncols+i];
                c2 = ic1[   j *ncols+i];
                cdif = c2 - c1;
                s2 = s2 + (double) cdif;
                if (c1 !=0 && c2 != 0){
                    s1 = s1 + (double) cdif;
                    kount++;
                }
                else {
                    nnull = nnull+1;
                }
            }
        }
        if (kount > nnull){
            sr = s1 / (double) kount;
        }
        else {
            sr = s2 / (double) ((ncols-1)*nrows);
        }
        
        /* estimate slope in direction of increasing column number, i.e. increasing X */
        s1=0.0;s2=0.0f;kount=0;nnull=0;
        for (j=0; j < nrows; j++) {
            for (i=1; i < ncols; i++) {
                c2 = ic1[j*ncols + i  ];
                c1 = ic1[j*ncols + i-1];
                cdif = c2 - c1;
                s2 = s2 + (double) cdif;
                if (c1 !=0 && c2 != 0){
                    s1 = s1 + (double) cdif;
                    kount++;
                }
                else {
                    nnull = nnull+1;
                }
            }
        }
        if (kount > nnull){
            sc = s1 / (double) kount;
        }
        else {
            sc = s2 / (double) (ncols*(nrows-1));
        }
    }
    else {
        sr = 0.0;
        sc = 0.0;
    }
    sxy.slopey=sr;
    sxy.slopex=sc;
    return (slopexy) sxy;
}

signed char circmeandev(signed char *ic1, signed char alpha, long n)
/* circular mean deviation, as defined by Mardia and Jupp [2000], pages 19-20 */
/* input phase values with 256 DN = 1 cycle, ranging from -128 to +127 */
/* 2009-JUN-14 Kurt Feigl
 * 2011-JUN-21 Kurt Feigl modify to ignore nulls
 */
{
    double s1=0.0f, s2=0.0f;
    signed char ic3;
    long i, kount, nnull;
    
    kount = 0;
    nnull = 0;
    s1 = 0.0;
    s2 = 0.0;
    
    for (i=0;i < n; i++){
        ic3 = ic1[i];
        s2 = s2 + (128 - abs(128 - abs(ic3-alpha)));
        if (ic3 != 0) {
            kount = kount + 1;
            s1 = s1 + (128 - abs(128 - abs(ic3-alpha)));
        }
        else{
            nnull = nnull + 1;
        }
    }
    
    if (kount > nnull) {
        return (signed char) rint(s1/(double)kount) ;
    }
    else {
        if (n > 0) {
            return (signed char) rint(s2/(double)n) ;
        }
        else{
            return (signed char) 0;
        }
    }
}


signed char circslopemeandev(signed char *ic1, signed char *i1mod, long nrows, long ncols)
/* circular mean deviation, as defined by Mardia and Jupp [2000], pages 19-20 */
/* input phase values with 256 DN = 1 cycle, ranging from -128 to +127 */
/* Kurt Feigl 2011-JUL-20*/
{
    double r, sum;
    long j, j1, j2, i, i1, i2, kount;
    signed char cobs, cmod, cdif, cret;
    
    sum = 0.0;
    kount = 0;
    j1 = 0;
    j2 = nrows-1;
    i1 = 0;
    i2 = ncols-1;
    
    for(j=j1; j<=j2; j++){	/* outer loop is rows NY */
        for(i=i1; i<=i2; i++){	/* inner loop is cols NX */
            cobs =   ic1[ncols*j + i];
            cmod = i1mod[ncols*j + i];
            if (cobs != 0) {
                kount = kount + 1;
                cdif = cobs - cmod;
                sum = sum + (double)(128 - (signed char)abs(128 - (signed char)abs(cdif)));
            }
        }
    }
    if (kount > 0) {
        cret = (signed char) rint(sum/(double)kount) ;
    }
    else{
        cret = (signed char) 0;
    }
    return cret;
}
signed char *applyramp(signed char *ic1, signed char *i1mod, signed char cmean, slopexy slopevector, long nrows, long ncols)
/* given 3 parameters (mean, slopex, slopey) apply model values over square patch*/
/* Kurt Feigl 2011-JUL-19*/
{
    double r, dx, dy, grx, gry;
    long j, j1, j2, i, i1, i2, kount;
    signed char cobs, cmod;
    
    kount = 0;
    j1 = 0;
    j2 = nrows-1;
    i1 = 0;
    i2 = ncols-1;
    for(j=j1; j<=j2; j++){	/* outer loop is rows NY */
        for(i=i1; i<=i2; i++){	/* inner loop is cols NX */
            /* separation from midpoint */
            dy = (double)j - ((double)j1 + (double)j2)/2.0;
            dx = (double)i - ((double)i1 + (double)i2)/2.0;
            /* gradient values */
            gry = slopevector.slopey;
            grx = slopevector.slopex;
            /* ramp model */
            r = (double)cmean + dx*grx + dy*gry;
            cmod = (signed char)rint(r);
            i1mod[ncols*j + i] = cmod;
        }
    }
    return i1mod;
}
/**************** END OF COPYRIGHTED SOFTWARE ******************/
/*
 * http://www.koders.com/info.aspx?c=ProjectInfo&pid=32Z3491Z2DSPRQMW5WTN225EMD&s=quadtree
 *
 * Project: sfox - Summary
 *
 * A 3d engine focused (until now) on terrain rendering with OpenGL.
 * The final aim is to make a Starfox-like 3d shoot-them-up. Everything is written in C.
 *
 * Development Status: 1 - Planning
 * Environment: X11 Applications
 * Intended Audience: End Users/Desktop
 * License: GNU General Public License (GPL)
 * Natural Language: English, French
 * Operating System: Linux
 * Programming Language: C
 * Topic: Games/Entertainment
 *
 * Registered: 2004-Jan-01 20:16
 */

void bbox2d_set(bbox2d *bb, unsigned int iw, unsigned int jn, unsigned int ie, unsigned int js) {
    assert(ie!=iw && js!=jn);
    
    vector2_set(&bb->corners[0], iw, jn);
    vector2_set(&bb->corners[1], iw, js);
    vector2_set(&bb->corners[2], ie, js);
    vector2_set(&bb->corners[3], ie, jn);
}
void vector2_set(vector2 *v, unsigned int x, unsigned int y) {
    v->x = x;
    v->y = y;
}


/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * This program demonstrates the use of the heap sort algorithm.  For
 * more information about this and other sorting algorithms, see
 * http://linux.wku.edu/~lamonml/kb.html
 *
 */

void heapSort(unsigned char numbers[], long array_size) {
    long i;
    unsigned char temp;
    void siftDown();
    
    for (i = (array_size / 2)-1; i >= 0; i--)
        siftDown(numbers, i, array_size);
    
    for (i = array_size-1; i >= 1; i--) {
        temp = numbers[0];
        numbers[0] = numbers[i];
        numbers[i] = temp;
        siftDown(numbers, 0, i-1);
    }
}


void siftDown(unsigned char numbers[], long root, long bottom) {
    int done;
    long maxChild;
    unsigned char temp;
    
    done = 0;
    while ((root*2 <= bottom) && (!done)) {
        if (root*2 == bottom)
            maxChild = root * 2;
        else if (numbers[root * 2] > numbers[root * 2 + 1])
            maxChild = root * 2;
        else
            maxChild = root * 2 + 1;
        
        if (numbers[root] < numbers[maxChild]) {
            temp = numbers[root];
            numbers[root] = numbers[maxChild];
            numbers[maxChild] = temp;
            root = maxChild;
        }
        else
            done = 1;
    }
}

char *strprefix(char *strin, char *strrep, char *strnew){ /* replace the prefix string strrep with strnew in strin. */
    char *strtemp, *tp;
/*if strrep is not present then just prepend strnew */
    
    if (strncmp(strin, strrep, strlen(strrep))==0){ /*test for existing prefix */
        strin= &strin[strlen(strrep)];
    }
    strtemp=(char *) malloc(sizeof(char) * strlen(strin)+strlen(strnew)+1);
    
    if ( strtemp == NULL){  /* failed to reallocate memory for string */
        fprintf(stderr, "Memory error in strprefix\n");
        exit(-1);}
    
    strcpy(strtemp, strnew);
    strcat(strtemp, strin);
    return(strtemp);
    
}

char *strextension(char *strin, char *extension){ /*replace or postpend extension */
    char *strtemp, *tp ;
    
    tp=strrchr(strin, '.'); /* find last "." */
    if (tp != NULL)  {  /* delete extension after . */
        memset(tp+1, 0, sizeof(char));
    }
    strtemp= (char *) malloc(sizeof(char)* strlen(strin)+strlen(extension)+1);
    strcpy(strtemp, strin);
    strcat(strtemp, extension);
    return(strtemp);
    
/*     tp=strrchr(flstname,".");   /* test if it has .pha extention and delete
 *
 *         if (strcmp(tp,".pha") == 0){  /*delete .qls file string 
 *         flstname[strlen(flstname)-4] = 0; 
 * } 
 *     strcat(flstname,".qls");   /*add .qls extension */
    
}

char *strreplace(char *st, char *orig, char *repl) {
    static char buffer[4096];
    char *ch;
    if (!(ch = strstr(st, orig)))
        return st;
    strncpy(buffer, st, ch-st);
    buffer[ch-st] = 0;
    sprintf(buffer+(ch-st), "%s%s", repl, ch+strlen(orig));
    return buffer;
}

