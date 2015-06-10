/* Program "cutswappad.c" 
   updated 2005-APR-25 to go faster Kurt
   updated the 14 April 2005 by Leonardon Jeremie
   from "cut_swap_pad_ci2.c" version 1.0 by Kurt Feigl
   
   Compiled with:

   gcc cutswappad.c -lm -g -o cutswappad

*/

/*
  The program reads a .CI2 file describing a radar image and fit the datas in order to output another .CI2 file 
  describing a rectangular part of the image possibly padded with 0 if the requested output image dimensions exceed the initial
  dimensions of the input image.

*/

#include <stdlib.h>
#include <stdio.h>

main (argc, argv)
int argc;
char **argv; {
  long i,j, n, m, nc_in, nl_in, nlskip, ncskip, nc_out, nl_out, nl, nc, nr, nw, ipix, npixin, npixout;
  int *i4_in, *i4_out, i4, j4, swap;
	
  FILE *fpi = NULL;
  FILE *fpo  = NULL;

  int swap_4byte();
  
  if (argc <= 1) {
    fprintf (stderr, "\n\n%s cuts the image described by a CI2 file in order to keep a given \nbar of interest ",argv[0]);
    fprintf (stderr, "and writes it in a new .CI2 file,\nwhere it is padded to fit the user stated dimensions \n");
   
    fprintf (stderr, "usage: %s input.ci2 output.ci2 nc_in nl_in swap ncskip nlskip ncout nlout\n",argv[0]);
    fprintf (stderr, "OR %s input.ci2 output.ci2 nc_in nl_in swap ncout nlout\n",argv[0]);
    fprintf (stderr, "OR %s input.ci2 output.ci2 nc_in nl_in swap\n",argv[0]);
    fprintf (stderr, "OR %s input.ci2 output.ci2 nc_in nl_in\n\n",argv[0]);
    fprintf (stderr, "swap has to be 1 for launching the non-swapping algorithm and 2 for the swapping one.\n\n");
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
    fprintf (stderr, "opened output file %s\n", argv[2]);
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


  if (NULL == argv[5]){
    swap=1;
    ncskip=0;
    nlskip=0;
    nc_out=nc_in;
    nl_out=nl_in;
    fprintf (stderr, "NCSKIP = %d\n", ncskip);
    fprintf (stderr, "NLSKIP = %d\n", nlskip);
    fprintf (stderr, "NC_OUT = %d\n", nc_out);
    fprintf (stderr, "NL_OUT = %d\n", nl_out);
    fprintf(stderr,"Real and Imaginary parts will not be swapped...\n");
  }
  else{   
    if ((swap = atoi (argv[5])) == 0) {
	fprintf(stderr,"State 1 for not swapping and 2 for swapping...\n");
	exit(-1);
      }
      else{
	
	if(swap == 2){
	  fprintf(stderr,"Real and Imaginary parts will be swapped...\n");
	}
	else
	  {
	    if(swap == 1){
	      fprintf(stderr,"Real and Imaginary parts will not be swapped...\n");
	    }
	    else{
	      fprintf(stderr,"State 1 for not swapping and 2 for swapping...\n");
	      exit(-1);
	    }
	  }
      }
   
    if (NULL == argv[6]){
      ncskip=0;
      nlskip=0;
      nc_out=nc_in;
      nl_out=nl_in;
      fprintf (stderr, "NCSKIP = %d\n", ncskip);
      fprintf (stderr, "NLSKIP = %d\n", nlskip);
      fprintf (stderr, "NC_OUT = %d\n", nc_out);
      fprintf (stderr, "NL_OUT = %d\n", nl_out);
      
    }
    else{
      if (NULL == argv[7]) {
	fprintf(stderr,"Not enough arguments...\n");
	exit(-1);
      }
      else{
	if (NULL == argv[8]){
	  if ((nc_out = atoi (argv[6])) == 0)
	    {
	      fprintf (stderr, "addbyte: Cannot read NCOUT %s\n", argv[6]);
	      exit (-1);
	    }
	  else {
	    fprintf (stderr, "NCOUT = %d\n", nc_out);
	  }      
	  
	  
	  if ((nl_out = atoi (argv[7])) == 0) {
	    fprintf (stderr, "addbyte: Cannot read NLOUT %s\n", argv[7]);
	    exit (-1);
	  }
	  else {
	    fprintf (stderr, "NLOUT = %d\n", nl_out);
	  }
	  ncskip=0;
	  nlskip=0;
	  fprintf(stderr,"No lines and no columns skipped...\n");
	}
	else{
	  
	  if ((ncskip = atoi (argv[6])) == 0)
	    {
	      fprintf (stderr, "addbyte: Cannot read NCSKIP %s\n", argv[6]);
	      exit (-1);
	    }
	  else {
	    fprintf (stderr, "NCSKIP = %d\n", ncskip);
	  }      
	  
	  
	  if ((nlskip = atoi (argv[7])) == 0) {
	    fprintf (stderr, "addbyte: Cannot read NLSKIP %s\n", argv[7]);
	    exit (-1);
	  }
	  else {
	    fprintf (stderr, "NLSKIP = %d\n", nlskip);
	  }
	  
	  if ((nc_out = atoi (argv[8])) == 0)
	    {
	      fprintf (stderr, "addbyte: Cannot read NCOUT %s\n", argv[5]);
	      exit (-1);
	    }
	  else {
	    fprintf (stderr, "NCOUT = %d\n", nc_out);
	  }	 	  	  
	  
	  if ((nl_out = atoi (argv[9])) == 0) {
	    fprintf (stderr, "addbyte: Cannot read NLOUT %s\n", argv[9]);
	    exit (-1);
	  }
	  else {
	    fprintf (stderr, "NLOUT = %d\n", nl_out);
	  }
	  if (NULL==argv[10])
	    {
	      fprintf(stderr,"Every parameter has been set.../n/n");
	    }
	  else{
	    fprintf(stderr,"Wrong amount of arguments...\n");
	    exit(-1);
	  }
	}
      }
    }
  }
  
  /* allocate memory for input array */

  npixin = nc_in * nl_in;
  i4_in = calloc (npixin,sizeof(int));
  if (i4_in == NULL) {fprintf(stderr, "failure to allocate memory for i4_in buffer\n"); exit(-1);}
  fprintf(stderr,"%d \n",npixin);
  
  /* allocate memory for output array */
  
  npixout = nc_out * nl_out;
  fprintf(stderr,"%d \n",npixout);
  i4_out = calloc (npixout,sizeof(int));
  if (i4_out == NULL) {fprintf(stderr, "failure to allocate memory for i4_out buffer\n"); exit(-1);}
  
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
  if ((nc_out+ncskip) < nc_in) nc = nc_out+nlskip;
  nl = nl_in;
  if ((nl_out+nlskip) < nl_in) nl = nl_out+nlskip;

  /* SLOW WAY 
  for (j = nlskip; j < nl; j++){
   
    for (i = ncskip; i < nc; i++){
      
      ipix    = j*nc_in + i;
      i4 = i4_in[ipix];
      if (swap==1){
	j4 = swap_4byte(i4);
      }
      else{
	j4=i4;
	  }
      i4_out[(j-nlskip)*nc_out+i-ncskip] = j4;
    }
  }

  */

  if (swap==1){
    for (j = nlskip; j < nl; j++){
      for (i = ncskip; i < nc; i++){
	ipix    = j*nc_in + i;
	i4 = i4_in[ipix];
	j4 = swap_4byte(i4);
	i4_out[(j-nlskip)*nc_out+i-ncskip] = j4;
      }
    }
  }
 else{
    for (j = nlskip; j < nl; j++){
      for (i = ncskip; i < nc; i++){
	ipix    = j*nc_in + i;
	i4 = i4_in[ipix];
	j4 = i4;
	i4_out[(j-nlskip)*nc_out+i-ncskip] = j4;
      }
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
	Purpose:	reorder a 4-byte word from 0123 to 2301 
				or from 2301 to 0123
	Usage:		unsigned long int swap_4byte ();
				unsigned long int word4;
				unsigned long int result;
				result = swap_4byte (word4);
	Input:		a 4-byte word in order 0123
	Output:		a 4-byte word in order 2301
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
			2005-APR-14 Leonardon Jeremie (modification of the swapping order)
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
