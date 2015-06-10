/*
 * diffpha reads 2 files and writes the difference
 *
 * Author:	Kurt Feigl
 * Date:	June 2000
 * Version:	1.3	
 *
 *changed to "diffpha" from "diffbyte" Peter Sobol 4/12/2011
 * 
 *
 *usage: diffpha infile1 infile2 nx ny outfile -nullopt
 *          where nullopt = 0 or 1
 * nullopt=0 or 1 passes zeros in infile1 to output file 
 * (provides noise suppressio when comparing QLS recontstructions to original pha files.
 */

#include <stdlib.h>
#include <stdio.h>

main (argc, argv)
int argc;
char **argv; {
	long i,j, n, m, nx,ny,nr=0;
        signed char ic1,ic2,ic3;
	int passnull=0;
	
	FILE *fpi1 = NULL;
	FILE *fpi2 = NULL;
	FILE *fpo  = NULL;
		
	if (argc <= 1) {
		fprintf (stderr, "diffpha differences two files, byte by byte\n\n");
		fprintf (stderr, "usage: diffpha infile1 infile2 nx ny outfile\n");
		fprintf (stderr, "usage: diffpha infile1 infile2 nx ny outfile -0\n");
		fprintf (stderr, "usage: diffpha infile1 infile2 nx ny outfile -1\n");
			exit (-1);
	      }

        if ((fpi1 = fopen (argv[1], "r")) == NULL) {
	  fprintf (stderr, "diffpha: Cannot find input file %s\n", argv[1]);
	  exit (-1);
	}
	else {
	  fprintf (stderr, "opened input file %s\n", argv[1]);
	}

        if ((fpi2 = fopen (argv[2], "r")) == NULL) {
	  fprintf (stderr, "diffpha: Cannot find input file %s\n", argv[2]);
	  exit (-1);
	}
	else {
	  fprintf (stderr, "opened input file %s\n", argv[2]);
	}
        
	if ((nx = atoi (argv[3])) == 0) {
	  fprintf (stderr, "diffpha: Cannot read NX %s\n", argv[3]);
	  exit (-1);
	}
	else {
	  fprintf (stderr, "NX = %d\n", nx);
	}
        
	if ((ny = atoi (argv[4])) == 0) {
	  fprintf (stderr, "diffpha: Cannot read NY %s\n", argv[4]);
	  exit (-1);
	}
	else {
	  fprintf (stderr, "NY = %d\n", ny);
	}
        
        if ((fpo = fopen (argv[5], "w")) == NULL) {
	  fprintf (stderr, "diffpha: Cannot open output file %s\n", argv[5]);
	  exit (-1);
	}
	
	if (argc >= 7) {
	  if (strcmp(argv[6],"-0") == 0) {
	    fprintf (stderr, "diffpha: Passing nulls from infiles to outfile\n");
	    passnull=1;
	  }
	  if (strcmp(argv[6],"-1") == 0) {
	    fprintf (stderr, "diffpha: Passing nulls from infile1 to outfile\n");
	    passnull=2;
	  }
	}

	if (passnull == 1) {
	  for (j = 0; j < ny; j++){
	    for (i = 0; i < nx; i++){
	      ic1 = fgetc (fpi1);
	      ic2 = fgetc (fpi2);
	      if (ic1 != 0 && ic2 != 0) {
		ic3 = ic1 - ic2;
	      }
	      else {
		ic3 = 0;
	      }
	      fputc ((ic3),fpo);
	      nr++;
	    }
	  }
	}
	else if (passnull == 2) {
	  for (j = 0; j < ny; j++){
	    for (i = 0; i < nx; i++){
	      ic1 = fgetc (fpi1);
	      ic2 = fgetc (fpi2);
	      if (ic1 != 0) {
		ic3 = ic1 - ic2;
	      }
	      else {
		ic3 = 0;
	      }
	      fputc ((ic3),fpo);
	      nr++;
	    }
	  }
	}
	else {
	  for (j = 0; j < ny; j++){
	    for (i = 0; i < nx; i++){
	      ic1 = fgetc (fpi1);
	      ic2 = fgetc (fpi2);
	      ic3 = ic1 - ic2;
	      fputc ((ic3),fpo);
	      nr++;
	    }
	  }
	}
	
	if (nx * ny != nr) { 
	  fprintf (stderr, "diffpha: counting error. Read %dbytes. Expected %d\n", nr,nx*ny);
	  exit (-1);
	}
	else {
	  fclose (fpo);	 	
	} 
}


