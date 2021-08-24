#include <stdio.h>
#include "nifti2_io.h"
#include <math.h>
#include <stdlib.h>


/* ---------------------------------------
   -------- FUNCTION DECLARATIONS -------- 
   --------------------------------------- */

float * expfunc(float *x, int n, float s, float r);

float * linspace(float a, float b, int n, float linspaceArray[]);

float ressumsq(float *y, float *yhat, int n);

float sum_array(float *a, int n);

/* --------------------------------------- */


typedef struct exp_params {
  float s;
  float r;
} exp_params;


/* ---------------------------------------
   --------         MAIN          -------- 
   --------------------------------------- */

int main(int argc, char * argv[])
{

  nifti_image *nim_echo1=NULL, *nim_echo2=NULL, *nim_echo3=NULL, *nim_outparam1=NULL, *nim_outparam2=NULL;
  char *fin1=NULL, *fin2=NULL, *fin3=NULL, *fout1=NULL, *fout2=NULL;
  float param1_low, param1_high, param2_low, param2_high;
  float p1_best, p2_best, rss_old, rss_new;
  int idx;
  int ngrid=200; /* the size of parameter set, ngrid^num_params is the number of models fit in the brute-force */
 
 /* array of echo times */
  float tes[3];
  tes[0]=0.00814;
  tes[1]=0.02147;
  tes[2]=0.03480;
 
  /* array for voxel data from 3 echos and for fit values */
  float y[3], *yhat=NULL;
  
  /* parameter arrays */
  float *param1_array, *param2_array;

  /* parameter struct */
  exp_params params;

  fin1 = argv[1];
  fin2 = argv[2];
  fin3 = argv[3];

  param1_low = atof(argv[4]);
  param1_high = atof(argv[5]);
  param2_low = atof(argv[6]);
  param2_high = atof(argv[7]);

  printf("\nParameter constraints\n");
  printf("S0 = %6.3f to %6.3f",param1_low,param1_high);
  printf("\nR2* = %6.3f to %6.3f\n",param2_low,param2_high);

  nim_echo1 = nifti_image_read(fin1, 1);
  nim_echo2 = nifti_image_read(fin2, 1);
  nim_echo3 = nifti_image_read(fin3, 1);

  /* get dimensions of input */
  int sizeSlice = nim_echo1->nz;
  int sizePhase = nim_echo2->nx;
  int sizeRead = nim_echo1->ny;
  int nrep = nim_echo1->nt;
  int nx = nim_echo1->nx;
  int nxy = nim_echo1->nx * nim_echo1->ny;
  int nxyz = nim_echo1->nx * nim_echo1->ny * nim_echo1->nz;
  int nvox = nrep * nxyz;

  /* get access to the data */
  float *nim_echo1_data = nim_echo1->data;
  float *nim_echo2_data = nim_echo2->data;
  float *nim_echo3_data = nim_echo3->data;

  /* parameter outputs */
  nim_outparam1 = nim_echo1;
  nim_outparam2 = nim_echo2;
  float *nim_outparam1_data = nim_outparam1->data;
  float *nim_outparam2_data = nim_outparam2->data;

  /* create parameter arrays */
  float param1_linspaceArray[ngrid], param2_linspaceArray[ngrid];
  param1_array = linspace(param1_low, param1_high, ngrid, param1_linspaceArray);
  param2_array = linspace(param2_low, param2_high, ngrid, param2_linspaceArray);

  //for (int i=0; i<ngrid; i++)
  //printf("\n%6.3f\t%6.3f",*(param1_array + i),*(param2_array + i));

  /* brute-force exponential fit */
  /* lots of loops - will be slow */

  printf("\nBrute-force exponential fit\n");

  int count=1;

  //nrep=1;
  //sizeSlice=12;
  //sizePhase=21;
  //sizeRead=21;

  for (int timestep=0; timestep<nrep; timestep++){
  for (int islice=0; islice<sizeSlice; islice++){
    for (int iy=0; iy<sizePhase; iy++){
  	for (int ix=0; ix<sizeRead; ix++){
  	  idx = nxyz*timestep + nxy*islice + nx*ix +iy;
  	  y[0] = *(nim_echo1_data + idx);
  	  y[1] = *(nim_echo2_data + idx);
  	  y[2] = *(nim_echo3_data + idx);
	  //printf("\n%6.3f %6.3f",y[0],*(nim_echo1_data + idx));
	  if (sum_array(y,3) > 1){
  	  p1_best = *param1_array;
  	  p2_best = *param2_array;
  	  rss_old = NAN;
	    for (int p1=0; p1<ngrid; p1++){
  	    for (int p2=0; p2<ngrid; p2++){
  	      params.s = *(param1_array + p1);
  	      params.r = *(param2_array + p2);
	      yhat = expfunc(tes,3,params.s,params.r);
	      //printf("\n\ny[0] = %6.3f\tyhat[0] = %6.3f\n",y[0],*yhat);
	      //printf("y[1] = %6.3f\tyhat[1] = %6.3f\n",y[1],*(yhat+1));
	      //printf("y[2] = %6.3f\tyhat[2] = %6.3f\n",y[2],*(yhat+2));
	      rss_new = ressumsq(y,yhat,3);
	      if ( p1 == 0 && p2 == 0)
		rss_old=rss_new;
	      if (rss_new < rss_old){
		p1_best = *(param1_array + p1);
		p2_best = *(param2_array + p2);
		rss_old = rss_new;
	      } 
	    }
	  }
	  }
	  else {
	    p1_best = 0;
	    p2_best = 0;
	  }
  	  if ( count % 1000 == 0 ){
  	    printf("\nS0 = %6.3f\tR2* = %6.3f\t(voxel %d of %d)",p1_best,p2_best,count,nvox);
  	  }
	  *(nim_outparam1_data + idx) = p1_best;
	  *(nim_outparam2_data + idx) = p2_best;
	  count++;
	}
      }
    }
  }

  printf("\n\nFINISHED\n\n\t...writing model fit parameters\n\n");

  fout1="S0.nii";
  fout2="R2star.nii";
  if (nifti_set_filenames(nim_outparam1,fout1,1,1)) return 1;
  nifti_image_write(nim_outparam1);
  nifti_image_free(nim_outparam1);
  if (nifti_set_filenames(nim_outparam2,fout2,1,1)) return 1;
  nifti_image_write(nim_outparam2);
  nifti_image_free(nim_outparam2);

  return 0;
}

/* --------------------------------------- */



/* --------------------------------------
   -------- FUNCTION DEFINITIONS -------- 
   -------------------------------------- */

float *expfunc(float *x, int n, float s, float r)
{
  
  /* the returned pointer is to static data that is overwritten with each call */
  static float yret[3];
  for (int i=0; i<n; i++){
    //printf("\n *(x+i) = %6.3f\tr = %6.3f\t -1*%6.3f*%6.3f = %6.3f",*(x+i),r,r,*(x+i),(-1.0*r*(*(x+i))));
    yret[i] = s * exp(-1.0*r*(*(x+i)));
  }
  //printf("\n\nINSIDE expfunc: s = %6.3f\tr=%6.3f\n",s,r);
  //printf("\ttes[0,1,2] = {%6.3f,%6.3f,%6.3f}",x[0],x[1],x[2]);
  //printf("\tyret[0,1,2] = {%6.3f,%6.3f,%6.3f}\n",yret[0],yret[1],yret[2]);
  return yret;
}

float *linspace(float a, float b, int n, float linspaceArray[])
{
  float stepsize = (b - a)/(n-1);
  for (int i=0; i<n; i++){
    linspaceArray[i] = a+i*stepsize;
  }
  return linspaceArray;
}

float ressumsq(float *y, float *yhat, int n)
{
  float diffsq[n];
  for (int i=0; i<n; i++){
    // printf("\n*(y+i) = %6.3f\t*(yhat+i) = %6.3f",*(y+i),*(yhat+i));
    diffsq[i] = (*(y+i) - *(yhat+i)) * (*(y+i) - *(yhat+i));
  }
  float rss = sum_array(diffsq,n);
  return rss;
}

float sum_array(float *a, int n)
{
  float sum=0;
  for (int i=0; i<n; i++){
    sum = sum + *(a+i);
  }
  return sum;
}
