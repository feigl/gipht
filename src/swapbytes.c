/*
 * swap byte order based on word size 
 *
 * Author:	Kurt Feigl
 * Date:	2005-MAY-22
 * Version:	1.0
 *
 * gcc swapbytes.c -o swapbytes	
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <netinet/in.h>

/*
#include <iostream.h>

float
Swap(float x)
{
  float swapped_x;
  for (int i = 0; i < 4; i++)
    ((char*)(&swapped_x))[i] = ((char*)(&x))[3 - i];
  return swapped_x;
}


NAME
       htonl, htons, ntohl, ntohs - convert values between host and network byte order

SYNOPSIS
       #include <netinet/in.h>

       uint32_t htonl(uint32_t hostlong);

       uint16_t htons(uint16_t hostshort);

       uint32_t ntohl(uint32_t netlong);

       uint16_t ntohs(uint16_t netshort);

DESCRIPTION
       The htonl() function converts the unsigned integer hostlong from host byte order to network byte order.

       The htons() function converts the unsigned short integer hostshort from host byte order to network byte order.

       The ntohl() function converts the unsigned integer netlong from network byte order to host byte order.

       The ntohs() function converts the unsigned short integer netshort from network byte order to host byte order.

       On  the  i80x86 the host byte order is Least Significant Byte first, whereas the network byte order, as used on the Inter-
       net, is Most Significant Byte first.
*/

main (argc, argv)
int argc;
char **argv; {
  long nr,nw;
  char *inword,*outword;
  int ios1,ios2,i;
  char fbytnam[80],foutnam[80],fsctnam[80];
  int bytes_per_word;
  FILE *fpi1   = NULL;
  FILE *fpout  = NULL;  
 
  if (argc <= 1) {
    fprintf (stderr, "swap_byte\n\n");
    fprintf (stderr, "usage: swap_bytes infile bytes_per_word outfile\n");
    exit (-1);
  }
  
  if ((fpi1 = fopen (argv[1], "r")) == NULL) {
    fprintf (stderr, "swap: Cannot find input file %s\n", argv[1]);
    exit (-1);
  }
  else {
    fprintf (stderr, "opened input file %s\n", argv[1]);
  }
 
 if ((bytes_per_word = atoi (argv[2])) == 0) {
    fprintf (stderr, "stack_pha3: Cannot read BYTES_PER_WORD %s\n", argv[4]);
    exit (-1);
  }
  else {
    if (bytes_per_word == 2 || bytes_per_word == 4 || bytes_per_word == 8){
      fprintf (stderr, "BYTES_PER_WORD = %d\n", bytes_per_word);
    }
    else {
      fprintf (stderr, "BAD BYTES_PER_WORD = %d\n", bytes_per_word);
    }
  }
  
  strcpy(foutnam,argv[3]);
  if ((fpout = fopen (foutnam, "w+")) == NULL) {
    fprintf (stderr, "swapbytes: Cannot open output file %s\n", foutnam);
    exit (-1);
  }
  else {
    fprintf (stdout, "Opened: %s\n", foutnam);
  }


  /* allocate memory for one word*/
  inword = calloc (1,bytes_per_word);
  if (inword == NULL) {fprintf(stderr, "failure to allocate space for inword buffer\n"); exit(-1);}
  outword = calloc (1,bytes_per_word);
  if (outword == NULL) {fprintf(stderr, "failure to allocate space for outword buffer\n"); exit(-1);}
 
  /* loop over all words */
  nr = 0;
  nw = 0;
  while ((ios1 = fread(inword,1,bytes_per_word,fpi1)) == bytes_per_word) {
    nr++; 
    for (i = 0; i < bytes_per_word; i++){
      outword[bytes_per_word-i-1]=inword[i];     
    }
/*
    if (nr < 10) {
      fprintf (stderr,"%2x %2x %2x %2x\n",inword[0],inword[1],outword[0],outword[1]);
    }
*/
    ios2 = fwrite(outword,1,bytes_per_word,fpout);
    if (ios2 == bytes_per_word) nw++;
  }

  
  if (nw == nr) {
    fprintf (stderr, "Read %12ld Wrote %12ld\n", nr,nw);
  }
  else {
    fprintf (stderr, "ERROR Read %12ld Wrote %12ld\n", nr,nw);
  }
   
  /* finish up */
  free (inword);
  free (outword);
  fclose (fpout);
  exit (0);
}

