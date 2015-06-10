#include <stdlib.h>
#include <stdio.h>

main (argc, argv)
int argc;
char **argv; {
  long i,j, n, m, nc_in, nl_in, nchg_out, nlhg_out, nlbd_out, ncbd_out, nc_out, nl_out, nl, nc, nr, nw, ipix, npixin, npixout;
  int *i4_in, *i4_out, i4, j4;
	
  FILE *fpi = NULL;
  FILE *fpo  = NULL;

  int swap_4byte();
  

  
  if (argc <= 1) {
    fprintf (stderr, "%s cuts the image described by a CI2 file in order to keep a smaller \narea of interest ",argv[0]);
    fprintf (stderr, "and writes it in a new .CI2 file\n");
   
    fprintf (stderr, "usage: %s input.ci2 output.ci2 nc_in nl_in nchg_out nlhg_out ncbd_out nlbd_out\n\n",argv[0]);
    exit (-1);
  }



  if ((fpi = fopen (argv[1], "r")) == NULL) {
    fprintf (stderr, "addbyte: Cannot find input file %s\n", argv[1]);
    exit (-1);
  }
  else {
    fprintf (stderr, "opened input file %s\n", argv[1]);
  }
  
  if ((fpo = fopen (argv[2], "w")) == NULL) {
    fprintf (stderr, "addbyte: Cannot open output file %s\n", argv[2]);
    exit (-1);
  }
  else {
    fprintf (stderr, "opened input file %s\n", argv[2]);
  }
  
  if ((nc_in = atoi (argv[3])) == 0) {
    fprintf (stderr, "addbyte: Cannot read NC_IN %s\n", argv[3]);
    exit (-1);
  }
  else {
    fprintf (stderr, "NC_IN = %d\n", nc_in);
  }
  
  if ((nl_in = atoi (argv[4])) == 0) {
    fprintf (stderr, "addbyte: Cannot read NL_IN %s\n", argv[4]);
    exit (-1);
  }
  else {
    fprintf (stderr, "NL_IN = %d\n", nl_in);
  }
  
  if ((nchg_out = atoi (argv[5])) == 0) {
    fprintf (stderr, "addbyte: Cannot read NCHG_OUT %s\n", argv[5]);
    exit (-1);
  }
  else {
    fprintf (stderr, "NCHG_OUT = %d\n", nchg_out);
  }
  
  if ((nlhg_out = atoi (argv[6])) == 0) {
    fprintf (stderr, "addbyte: Cannot read NLHG_OUT %s\n", argv[6]);
    exit (-1);
  }
  else {
    fprintf (stderr, "NLHG_OUT = %d\n", nlhg_out);
  }
  
    
  if ((ncbd_out = atoi (argv[7])) == 0) {
    fprintf (stderr, "addbyte: Cannot read NCBD_OUT %s\n", argv[5]);
    exit (-1);
  }
  else {
    fprintf (stderr, "NCBD_OUT = %d\n", ncbd_out);
  }
  
  if ((nlbd_out = atoi (argv[8])) == 0) {
    fprintf (stderr, "addbyte: Cannot read NLBD_OUT %s\n", argv[6]);
    exit (-1);
  }
  else {
    fprintf (stderr, "NLBD_OUT = %d\n", nlbd_out);
  }
  
  /* allocate memory for input array */
  npixin = nc_in * nl_in;
  i4_in = calloc (npixin,sizeof(int));
  if (i4_in == NULL) {fprintf(stderr, "failure to allocate memory for i4_in buffer\n"); exit(-1);}
  fprintf(stderr,"%d \n",npixin);
  
  /* allocate memory for output array */
  nc_out = ncbd_out-nchg_out+1;
  nl_out = nlbd_out-nlhg_out+1;
  npixout = nc_out * nl_out;
  fprintf(stderr,"%d \n",npixout);
  i4_out = calloc (npixout,sizeof(int));
  if (i4_out == NULL) {fprintf(stderr, "failure to allocate memory for i4_in buffer\n"); exit(-1);}
  
  /* read the whole array */
  nr = fread(i4_in,sizeof(int),npixin,fpi);
  
  if (nr == npixin) { 
    fprintf (stderr, "%s: read %d bytes for %d pixels\n", argv[0],4*nr, npixin);
  }
  else {
    fprintf (stderr, "%s: counting errror. Read %d bytes. Expected %d bytes.\n", argv[0],nr,npixin);
    exit(-1);
  }
  
  
  /* determine number of lines and columns to copy */
  nc = nc_in;
  if (nc_out < nc_in) nc = nc_out;
  nl = nl_in;
  if (nl_out < nl_in) nl = nl_out;
 

  nchg_out=nchg_out-1;
  nlhg_out=nlhg_out-1;

  for (j = 0; j < nl_out; j++){
   
    for (i = 0; i < nc_out; i++){
      
      ipix    = (j+nlhg_out)*nc_in + i+nchg_out;
      i4 = i4_in[ipix];
      j4 = swap_4byte(i4);
      i4_out[j*nc_out+i] = j4;
    }
  }
  
  /* write the whole array */

  nw = fwrite(i4_out,sizeof(int),npixout,fpo);

 
  if (nw == npixout) { 
    fprintf (stderr, "%s: wrote %d bytes for %d pixels.\n", argv[0],4*nw,npixout);
  }
  else {
    fprintf (stderr, "%s: counting errror. Wrote %dbytes. Expected %d\n", argv[0],4*nw,npixout);
    fclose (fpo);	 	
  }
  
  free (i4_in);
  free (i4_out);
  exit(0);
}



/*===========================================================================*/
/* SEED reader     |               swap_4byte              |    subprocedure */
/*===========================================================================*/
/*
	Name:		swap_4byte
	Purpose:	reorder a 4-byte word from 3210 to 0123 (MSB-first to MSB-last)
				or from 0123 to 3210
	Usage:		unsigned long int swap_4byte ();
				unsigned long int word4;
				unsigned long int result;
				result = swap_4byte (word4);
	Input:		a 4-byte word in order 3210
	Output:		a 4-byte word in order 0123
	Externals:	none
	Warnings:	none
	Errors:		none
	Called by:	anything
	Calls to:	none
	Algorithm:	Using a union between an unsigned long int and 4 chars,
				shuffle the bytes around to achieve the reverse word order.
	Notes:		none
	Problems:	none known
	References:	Halbert et al, 1988; see main routine
	Language:	C, hopefully ANSI standard
	Author:		Dennis O'Neill
	Revisions:	11/09/88  Dennis O'Neill  original version
			11/21/88  Dennis O'Neill  Production release 1.0
                        2005-MAR-01 Kurt Feigl make long
*/


int swap_4byte (word4) 
     int word4; {
  union
  {
    unsigned char character[4];
    long integer;
  } swap4byte;
  
  char temp0;
  char temp1;
  
  swap4byte.integer = word4;
  
  temp0 = swap4byte.character[0];
  temp1 = swap4byte.character[1];
  swap4byte.character[0] = swap4byte.character[2];
  swap4byte.character[1] = swap4byte.character[3];
  swap4byte.character[2] = temp0;
  swap4byte.character[3] = temp1;
  
  return (swap4byte.integer);
}

