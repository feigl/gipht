
/* Change History

    cc ps_filt2.c -O5 -xlibmieee -lm -o ps_filt2

    changed arguments
    changed file inputs from single complex number to phase and amplitude files
    improve documentation slightly
    changed how it deals with zeros in diapason format files - zero phase but non-zero amplitude is changed to phase=1
     							     - output==0 if nulled but only ~=0 if it is smoothed output
								=> change ~=0 cells to 1 or 255

20140531 - Kurt Feigl

Extract gradients
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#define NFFT     64
#define STEP    NFFT/4
#define PI 3.141592654
#define ROUND(x) ((((x)-(short)(x))>0.5)?((short)(x)+1):(short)(x))

typedef struct{float re,im;} fcpx;

float wf[NFFT][NFFT];
unsigned int nfft[3];

void fourn(float *, unsigned int *, int ndim, int isign);
fcpx **alloc_2d_c(int nrl, int nrh, int ncl, int nch);

int main(int argc, char **argv)
{
   fcpx *bufcz, **cmp;
   fcpx **sm, **seg_fft, *seg_fftb;

   float *dt;
   float phase,amplt;  	/* read these from oct files */
   float real,imag;  	/* calculate these from phase,amplt */
   double alpha;
   float t1, t2, tmp, peak;
   float w[NFFT], amp[NFFT][NFFT];

   int width;
   int count;
   int nlines=0;
   int nbytes1,nbytes2;
   int phase_int;
   int step;
   int xw,yh;
   int i,j,i1,j1;
   int ndim;
   int isign;
   int ic, lc;
   int imax,jmax; /* indices of max power */
   
   /*  FILE *int_file, *sm_file;*/
   FILE *pha_file, *amp_file, *sm_file, *int_file;  /* intfile is temporarily created form amp_file and pha_file */

   struct stat *file_stat = malloc(sizeof(struct stat));


   fprintf(stdout,"*** power spectrum smoothing of interferogram                ***\n");
   fprintf(stdout,"*** Original by Zhong Lu                                     ***\n");
   fprintf(stdout,"***   thanking C. Werner and R. Goldstein for their advice   ***\n");
   fprintf(stdout,"*** circa 2000 Tim Wright                                    ***\n");
   fprintf(stdout,"***      Modified to work with 8 Bit Phase/Amplitude...      ***\n");
   fprintf(stdout,"*** Version 3.0 20140530                                     ***\n");
   fprintf(stdout,"***      Kurt Feigl to record gradients                      ***\n");
   if(argc != 6){
/*    fprintf(stderr,"\nusage: %s <ifm> <sm-ifm> <width> <alpha> \n\n",argv[0]) ;*/
//      fprintf(stderr,"\nusage: %s <pha.oct> <amp.oct> <smp.oct> <ncol> <alpha> \n\n",argv[0]) ;
     fprintf(stderr,"\nusage: %s <pha.oct> <amp.oct> <smp.oct> <ncol> <alpha> \n\n",argv[0]) ;
     fprintf(stderr,"input parameters: \n");
//      fprintf(stderr,"pha.oct  (Input)  Phase     (Unsigned 8Bit Integer 0 = -pi; 255 = pi)\n");
     fprintf(stderr,"pha.oct  (Input)  Phase     (Unsigned 8Bit Integer -128 = -pi; 127 = pi)\n");
//      fprintf(stderr,"amp.oct  (Input)  Amplitude (Unsigned 8Bit Integer 0 - 255)\n");
/*    fprintf(stderr,"ifm      (Input)complex interferogram (4-byte for real and 4-byte for imaginary)\n");*/
     fprintf(stderr,"smp.oct  (Output) Smoothed Phase (Unsigned 8Bit Integer -128 = -pi; 127 = pi)\n");
     fprintf(stderr,"ncols    (Input)  Number of columns in all 3 files (pixels per line);\n");
/*     fprintf(stderr,"alpha    (Input)  Power factor in Goldstein & Werner technique (recommend: 0.0 < 0.7 < 1.0);\n"); */
     fprintf(stderr,"alpha    (Input)  Alpha exponent in Goldstein & Werner technique (recommend: 0.0 < 0.7 < 1.0);\n");
     exit(-1);
   }

   pha_file = fopen(argv[1],"rb");
   if (pha_file == NULL){fprintf(stderr,"cannot open phase file: %s\n",argv[1]); exit(-1);}

   amp_file = fopen(argv[2],"rb");
   if (amp_file == NULL){fprintf(stderr,"cannot open amplitude file: %s\n",argv[1]); exit(-1);}

   sm_file = fopen(argv[3],"w");
   if (sm_file == NULL){fprintf(stderr,"cannot create smoothed interferogram file: %s\n",argv[2]); exit(-1);}

   int_file = fopen("temp.cr4","wb");
   if (int_file == NULL){fprintf(stderr,"cannot create temp complex file: temp.cr4\n"); exit(-1);}

   sscanf(argv[4],"%d",&width);
   sscanf(argv[5],"%lf",&alpha);



   /* Work out number of lines in file and check that the files are the same size */

   stat(argv[1],file_stat);
   nbytes1=(int)file_stat->st_size;
   fprintf (stderr, "nbytes1 = %d\n",nbytes1);

   stat(argv[2],file_stat);
   nbytes2=(int)file_stat->st_size;
   fprintf (stderr, "nbytes2 = %d\n",nbytes2);

   if (nbytes1==nbytes2) {
      nlines=(int)(nbytes1/width);
   }
   else {
      fprintf(stderr, "\t Amplitude and Phase Files are different sizes\n");
      exit(-1);
   }
 
/*  fseek(int_file, 0L, 2);
   nlines=(int)ftell(int_file)/(width*2*sizeof(float));
   rewind(int_file);*/

   step = STEP;

   fprintf(stderr, "Input parameters:\n");
   fprintf(stderr, "\tInput Phase: %s\n", argv[1]);
   fprintf(stderr, "\tInput Amplitude: %s\n", argv[2]);
   fprintf(stderr, "\tSmoothed phase: %s\n", argv[3]);
   fprintf(stderr, "\tPower spectrum alpha: %f\n", alpha);
   fprintf(stderr, "\tFFT window size: %d\n", NFFT);
   fprintf(stderr, "\tImage width: %d\n", width);
   fprintf(stderr, "\tImage length: %d\n", nlines);

   xw=width;
   yh=nlines;
   fprintf(stdout,"\tarray width, height: %5d %5d\n",xw,yh);

   cmp = alloc_2d_c(0, nlines-1, 0,width-1);
   sm = alloc_2d_c(0,nlines-1, 0,width);
   if (cmp == NULL || sm == NULL) {fprintf(stderr, "failure to allocate space for cmp and sm buffers\n"); exit(-1);}

   seg_fftb = (fcpx *)malloc(sizeof(fcpx)*NFFT*NFFT);
   if(seg_fftb == NULL){fprintf(stderr,"failure to allocate space for complex data\n"); exit(-1);}

   seg_fft = (fcpx **)malloc(sizeof(fcpx *)*NFFT);
   if(seg_fft == NULL){fprintf(stderr,"failure to allocate space for complex data pointers\n"); exit(-1);}

   bufcz = (fcpx *)malloc(sizeof(fcpx)*width);
   if(bufcz == NULL){fprintf(stderr,"failure to allocate space for input line buffer\n"); exit(-1);}

   for(i=0; i < NFFT; i++)seg_fft[i] = seg_fftb  + i*NFFT;

   for(i=0; i < nlines; i++){
     for(j=0; j < width; j++){
      sm[i][j].re=0.0; sm[i][j].im=0.0;
     }
   }

   for(j=0; j < width; j++){bufcz[j].re=0.; bufcz[j].im=0.;}

   for (i=0; i <= NFFT/2-1; i++) {
     w[i] = 2.0*(i+1)/(float)NFFT;
     w[NFFT-i-1]=w[i];
   }

   for (i=0; i < NFFT; i++)
       for (j=0; j < NFFT; j++)wf[i][j] = w[i]*w[j];
   
   /* Merge amp and phase into "temp.cr4" */
   count=0;
   
   while (count < nbytes1)	{
       phase = (float)fgetc(pha_file);
       amplt = (float)fgetc(amp_file);
// 		if (phase==0&&amplt!=0) phase=1;
// 		real = (amplt/255)*cos(phase*2*PI/255);
// 		imag = (amplt/255)*sin(phase*2*PI/255);
       if (amplt==0) phase=0;
       real = (amplt/256)*cos(phase*2*PI/256);
       imag = (amplt/256)*sin(phase*2*PI/256);
       fwrite (&real, sizeof(float),1,int_file);
       fwrite (&imag, sizeof(float),1,int_file);
       count ++;
   }

   printf("\n Created temporary Complex File: temp.cr4\n");
   
   /* Close file for writing and reopen for reading */
   fclose (int_file);
   fclose (pha_file);
   fclose (amp_file);
   int_file = fopen("temp.cr4","r");
   if (int_file == NULL){fprintf(stderr,"cannot create temp complex file: temp.cr4\n"); exit(-1);}



   fseek(int_file, 0, 0);

   for (i=0; i < yh; i ++)

     fread((char *)cmp[i], sizeof(fcpx), width, int_file);

#if 1
   fprintf(stderr, "Please wait, I am filtering the noisy ifms\n");
   fprintf(stderr, "I am slow, but you will like the results\n\n");
   lc=0;
   
   /* loop over patches with overlap by step */
   for (i=-NFFT+step; i <= yh-step; i += step){
     for (j=-NFFT+step; j < width-step; j += step){

       dt=(float *)&seg_fft[0][0];
       ic = 0;
       for (i1=0; i1 < NFFT; i1++){
         for (j1=0; j1 < NFFT; j1++){
           /* pad with zeros */
           if ( (i+i1) < 0 || (i+i1) >= yh || (j+j1) < 0 || (j+j1) >= width ) {
             dt[ic++] = 0.0;
             dt[ic++] = 0.0;
           }
           else {
             /* normalize magnitude to unity */
             t1 = cmp[i+i1][j+j1].re;
             t2 = cmp[i+i1][j+j1].im;
             tmp = sqrt(t1*t1 + t2*t2);
             if (tmp != 0.0) {
               dt[ic++]=t1/tmp;
               dt[ic++]=t2/tmp;
             }
             else {
               dt[ic++]=0.0;
               dt[ic++]=0.0;
             }
           }
         }
       }

       /* Fourier transform */
       ndim=2, isign=1, nfft[1]=NFFT, nfft[2]=NFFT, nfft[0]=0;
       fourn(dt-1, nfft, ndim, isign);

       /* look for maximum spectral amplitude, a.k.a. power */
       peak = 0.0;
       for (i1=0; i1 < NFFT; i1++){
         for (j1=0; j1 < NFFT; j1++){
           /* calculate power */
           tmp = sqrt(seg_fft[i1][j1].re*seg_fft[i1][j1].re + seg_fft[i1][j1].im*seg_fft[i1][j1].im);
           /* if power is greater than peak, then update peak */
           if (peak < tmp) {
               peak = tmp;             
               printf("%3d %3d %5d %5d %12.4e %12.4e %12.4e\n",i1,j1,i,j,seg_fft[i1][j1].re,seg_fft[i1][j1].im,tmp);
           }
           /* record power */
           amp[i1][j1] = tmp;
         }
       }

       /* replace spectral value by amplified value */
       if (peak != 0.0) {
         for (i1=0; i1 < NFFT; i1++){
           for (j1=0; j1 < NFFT; j1++){
             seg_fft[i1][j1].re = (float)(seg_fft[i1][j1].re*pow((amp[i1][j1]/peak), alpha)/(float)NFFT);
             seg_fft[i1][j1].im = (float)(seg_fft[i1][j1].im*pow((amp[i1][j1]/peak), alpha)/(float)NFFT);
           }
         }
       }

       /* inverse Fourier transform */
       isign=-1;
       fourn(dt-1, nfft, ndim, isign);

       for (i1=0; i1 < NFFT; i1++){
         for (j1=0; j1 < NFFT; j1++){
           tmp = sqrt(seg_fft[i1][j1].re*seg_fft[i1][j1].re + seg_fft[i1][j1].im*seg_fft[i1][j1].im);
           seg_fft[i1][j1].re = (float)(seg_fft[i1][j1].re*wf[i1][j1]/tmp);
           seg_fft[i1][j1].im = (float)(seg_fft[i1][j1].im*wf[i1][j1]/tmp);
         }
       }

       /* replace phase value with filtered (smoothed) value */
       for (i1=0; i1 < NFFT; i1++){
         for (j1=0; j1 < NFFT; j1++){
           /* check that pixel is not in the padded area */
           if ( (i+i1) >= 0 && (i+i1) < yh && (j+j1) >= 0 && (j+j1) < width && cmp[i1+i][j1+j].re != 0.0 && cmp[i1+i][j1+j].im != 0.0 ) {
             sm[i1+i][j1+j].re += seg_fft[i1][j1].re;
             sm[i1+i][j1+j].im += seg_fft[i1][j1].im;
           }
         }
       }

     }

     lc += step;
     fprintf(stderr, "processing line: %5d\r", i);
   }

/* Close file for reading and reopen for writing !! There's probably a more elegant solution but this works !!  */
	fclose (int_file);
   	int_file = fopen("temp.cr4","w");
   	if (int_file == NULL){fprintf(stderr,"cannot rewrite temp complex file: temp.cr4\n"); exit(-1);}

/* Write temporary complex output file */
   fprintf(stderr, "writing output image \n");
   for (i=0; i < yh ; i++)
     fwrite((char *)sm[i], sizeof(fcpx), width, int_file);

/* Once more...Close for writing and reopen for reading */
	fclose (int_file);
   	int_file = fopen("temp.cr4","rb");
   	if (int_file == NULL){fprintf(stderr,"cannot create temp complex file: temp.cr4\n"); exit(-1);}

/* Convert complex array into octal and write Smoothed Phase File */
	count = 0;
	while (count < nbytes1)	{
		fread(&real, sizeof(float),1,int_file);
		fread(&imag, sizeof(float),1,int_file);

// /* Changed to deal with diapason nulls (phase==0) properly */		
// 		phase_int = (int)(ROUND((255*atan2(imag,real)/2/PI)));
// 		if (atan2(imag,real)==0) phase_int=0;
// 		else if (phase_int==0&&atan2(imag,real)>0) phase_int++;
// 		else if (phase_int==0&&atan2(imag,real)<0) phase_int--;
// 		while (phase_int > 255) phase_int -=255;
// 		while (phase_int < 0) phase_int +=255;	
/* Changed to deal with diapason nulls (phase==0) properly */	
        /* 20140530 use 256 */
		phase_int = (int)(ROUND((256*atan2(imag,real)/2/PI)));
		if (atan2(imag,real)==0) phase_int=0;
		else if (phase_int==0&&atan2(imag,real)>0) phase_int++;
		else if (phase_int==0&&atan2(imag,real)<0) phase_int--;
		while (phase_int > 127) phase_int -=256;
		while (phase_int < -128) phase_int +=256;	

		fputc(phase_int, sm_file);

		count ++;
	}

   	printf("\n finished***\n");




#endif
system("rm temp.cr4");
fclose (int_file);
fclose (sm_file);
   return(0);
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
void fourn(float data[], unsigned int nn[], int ndim, int isign)
{
   int idim;
   unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
   unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
   float tempi,tempr;
   double theta,wi,wpi,wpr,wr,wtemp;

   for (ntot=1,idim=1;idim<=ndim;idim++)
     ntot *= nn[idim];
   nprev=1;
   for (idim=ndim;idim>=1;idim--) {
     n=nn[idim];
     nrem=ntot/(n*nprev);
     ip1=nprev << 1;
     ip2=ip1*n;
     ip3=ip2*nrem;
     i2rev=1;
     for (i2=1;i2<=ip2;i2+=ip1) {
       if (i2 < i2rev) {
         for (i1=i2;i1<=i2+ip1-2;i1+=2) {
           for (i3=i1;i3<=ip3;i3+=ip2) {
             i3rev=i2rev+i3-i2;
             SWAP(data[i3],data[i3rev]);
             SWAP(data[i3+1],data[i3rev+1]);
           }
         }
       }
       ibit=ip2 >> 1;
       while (ibit >= ip1 && i2rev > ibit) {
         i2rev -= ibit;
         ibit >>= 1;
       }
       i2rev += ibit;
     }
     ifp1=ip1;
     while (ifp1 < ip2) {
       ifp2=ifp1 << 1;
       theta=isign*6.28318530717959/(ifp2/ip1);
       wtemp=sin(0.5*theta);
       wpr = -2.0*wtemp*wtemp;
       wpi=sin(theta);
       wr=1.0;
       wi=0.0;
       for (i3=1;i3<=ifp1;i3+=ip1) {
         for (i1=i3;i1<=i3+ip1-2;i1+=2) {
           for (i2=i1;i2<=ip3;i2+=ifp2) {
             k1=i2;
             k2=k1+ifp1;
             tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
             tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
             data[k2]=data[k1]-tempr;
             data[k2+1]=data[k1+1]-tempi;
             data[k1] += tempr;
             data[k1+1] += tempi;
           }
         }
         wr=(wtemp=wr)*wpr-wi*wpi+wr;
         wi=wi*wpr+wtemp*wpi+wi;
       }
       ifp1=ifp2;
     }
     nprev *= n;
   }
}

#undef SWAP

fcpx **alloc_2d_c(int nrl, int nrh, int ncl, int nch)
{
   int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   fcpx **m;

   m=(fcpx **)malloc((size_t)((nrow+1)*sizeof(fcpx*)));
   if (!m) fprintf(stderr, "allocation failure 1 in alloc_2d_c()");
   m += 1;
   m -= nrl;

   m[nrl]=(fcpx *) malloc((size_t)((nrow*ncol+1)*sizeof(fcpx)));
   if (!m[nrl]) fprintf(stderr, "allocation failure 2 in alloc_2d_c()");
   m[nrl] += 1;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   return m;
}
