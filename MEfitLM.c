/* Adapted from code provided at https://www.gnu.org/software/gsl/manual/html_node/Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html */

/* Mono-exponential fitting of multi-echo fMRI data using non-linear least-squares based Levenberg Marquardt (LM) algorithm.
   set of predefined parameter estimates are used to determine initial parameter starting points prior to LM optimisation

   EPI echo data needs to be FLOAT (datatype =16, 4 bytes/voxel) and mask dataset is UNSIGNED_CHAR (datatype=2, 1 byte/voxel)

   Usage: MEfitLM echo1.nii echo2.nii echo3.nii mask.nii prefix
   
   Joe Whittaker (2018) */

#include <stdio.h>
#include "nifti2_io.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#define N 3 // define number of data points to fit (i.e. number of echos)
#define P 2 // define number of parameters (S0 and R2*)
#define sqr(x) ((x)*(x))

#define ERREX(str) (fprintf(stderr, "\n** %s\n",str),exit(1)) // taken from afni src code

/* ---------------------------------
   ----- FUNCTION DECLARATIONS -----
   --------------------------------- */

int exp_f(const gsl_vector *x, void *data, gsl_vector *f);

int exp_df(const gsl_vector *x, void *data, gsl_matrix *J);

float *mean_vol(float *a, int nx, int ny, int nz, int nt, float outputarray[]);

float *ME_LMfit_function(float *yin, float *xin, float x_init[P], float weights[N], float params[P]);

float *ME_BFfit_function(float *yin, float *xin, int n, int ngrid, float p1[2], float p2[2], float params[P]);

float *linspace(float a, float b, int n, float linspaceArray[]);

float ressumsq(float *y, float *yhat, int n);

float sum_array(float *a, int n);

float calc_rsqr(float *y, float *x, float params[P], int n);

/* --------------------------------- */

  struct data {
    size_t n;
    double *s;
    double *te;
  };



/* ---------------------------------
   -------------- MAIN -------------
   --------------------------------- */

int main(int argc, char *argv[])
{

  clock_t begin = clock();

  nifti_image *nim_echo1=NULL, *nim_echo2=NULL, *nim_echo3=NULL, *nim_mask=NULL;
  nifti_image  *nim_outparam1=NULL, *nim_outparam2=NULL, *nim_rsqr=NULL;
  char *fin1=NULL, *fin2=NULL, *fin3=NULL, *fmask=NULL;

  int idx;

  float tes[3]; // array of echo times (s)
  tes[0]=0.00814;
  tes[1]=0.02147;
  tes[2]=0.03480;

  fin1 = argv[1];
  fin2 = argv[2];
  fin3 = argv[3];
  fmask = argv[4];

  nim_echo1 = nifti_image_read(fin1,1);
  nim_echo2 = nifti_image_read(fin2,1);
  nim_echo3 = nifti_image_read(fin3,1);
  nim_mask = nifti_image_read(fmask,1);

  if (nim_echo1->datatype != 16 ||
      nim_echo2->datatype != 16 ||
      nim_echo3->datatype != 16 ||
      nim_mask->datatype != 2)
    {
      ERREX("Echo data needs to be FLOAT and mask data needs to be UNSIGNED_CHAR!");
    }

  /* get dimensions of input */
  int sizeSlice = nim_echo1->nz;
  int sizePhase = nim_echo1->nx;
  int sizeRead = nim_echo1->ny;
  int nrep = nim_echo1->nt;
  int nx = nim_echo1->nx;
  int nxy = nim_echo1->nx * nim_echo1->ny;
  int nxyz = nim_echo1->nx * nim_echo1->ny * nim_echo1->nz;

  /* get access to the data */
  float *nim_echo1_data = nim_echo1->data;
  float *nim_echo2_data = nim_echo2->data;
  float *nim_echo3_data = nim_echo3->data;
  unsigned char *nim_mask_data = nim_mask->data;

  /*for (int i=1000; i<3000; i++)
    {
      if(*(nim_mask_data + i) == 1)
      printf("%d",1);
    }
  if(1)
  exit(1);*/

  /* parameter outputs */
  nim_outparam1 = nim_echo1;
  nim_outparam2 = nim_echo2;
  nim_rsqr = nim_echo3;
  float *nim_outparam1_data = nim_outparam1->data;
  float *nim_outparam2_data = nim_outparam2->data;
  float *nim_rsqr_data = nim_rsqr->data;

  /* --- Use bruteforce fit on mean across volumes to get initial voxel-wise starting points for LM algorithm --- */

  printf("\n\nCalculating initial parameter values for LM algorithm...");

  float *nim_echo1_data_mean=NULL, *nim_echo2_data_mean=NULL, *nim_echo3_data_mean=NULL;
  float nim_echo1_data_mean_outarray[nxyz], nim_echo2_data_mean_outarray[nxyz], nim_echo3_data_mean_outarray[nxyz];
  nim_echo1_data_mean = mean_vol(nim_echo1_data,sizeRead,sizePhase,sizeSlice,nrep,nim_echo1_data_mean_outarray);
  nim_echo2_data_mean = mean_vol(nim_echo2_data,sizeRead,sizePhase,sizeSlice,nrep,nim_echo2_data_mean_outarray);
  nim_echo3_data_mean = mean_vol(nim_echo3_data,sizeRead,sizePhase,sizeSlice,nrep,nim_echo3_data_mean_outarray);

  float *params, sig[N];
  float paramsArray[P];
  float BFparam1[nxyz], BFparam2[nxyz];
  float p1range[2]={100.0,10000.0}, p2range[2]={10.0,200.0};

  for (int islice=0; islice<sizeSlice; islice++){
    for (int iy=0; iy<sizePhase; iy++){
      for (int ix=0; ix<sizeRead; ix++){
	idx = nxy*islice + nx*ix +iy;
	sig[0] = *(nim_echo1_data_mean + idx);
	sig[1] = *(nim_echo2_data_mean + idx);
	sig[2] = *(nim_echo3_data_mean + idx);

	if (sum_array(sig,3)>1 && *(nim_mask_data + idx)==1)
	  {
	    params = ME_BFfit_function(sig,tes,N,50,p1range,p2range,paramsArray);
	    BFparam1[idx] = params[0];
	    BFparam2[idx] = params[1];
	    //printf("\n%f %f %f\t%f %f",s[0],s[1],s[2],params[0],params[1]);
	  }
	else
	  {
	    BFparam1[idx] = 0.0;
	    BFparam2[idx] = 0.0;
	  }
      }
    }
  }

  /* ------------------------------------------------------------------------------------------------------------- */

  /* --- Do LM nonlinear fit on every voxel for every time point using previously calculated start points -------- */

  printf("\nVoxel-wise LM optimisation...\n");

  float weights[N] = {1.0, 1.0, 1.0};
  float x_init[2];

  for (int timestep=0; timestep<nrep; timestep++){

	  if (timestep%10 == 0)
	    printf("\t...volume %d of %d\n",timestep,nrep);

    for (int islice=0; islice<sizeSlice; islice++){
      for (int iy=0; iy<sizePhase; iy++){
	for (int ix=0; ix<sizeRead; ix++){
	  idx = nxyz*timestep + nxy*islice + nx*ix + iy;
	  sig[0] = *(nim_echo1_data + idx);
	  sig[1] = *(nim_echo2_data + idx);
	  sig[2] = *(nim_echo3_data + idx);
	  x_init[0] = BFparam1[nxy*islice + nx*ix + iy];
	  x_init[1] = BFparam2[nxy*islice + nx*ix + iy];

	  if (sum_array(sig,3) > 1 && *(nim_mask_data + nxy*islice + nx*ix +iy)==1)
	    {
	      //printf("\n%f ",sig[0]);
	      params = ME_LMfit_function(sig,tes,x_init,weights,paramsArray);
	      *(nim_outparam1_data + idx) = params[0];
	      *(nim_outparam2_data + idx) = params[1];
	      //printf("\n%f %f %f\t%f %f\t%f %f",sig[0],sig[1],sig[2],x_init[0],x_init[1],params[0],params[1]);
	      *(nim_rsqr_data + idx) = calc_rsqr(sig,tes,params,N);
	    }
	  else
	    {
	      *(nim_outparam1_data + idx) = 0.0;
	      *(nim_outparam2_data + idx) = 0.0;
	      *(nim_rsqr_data + idx) = 0.0;
	    }
	}
      }
    }
  }

  /* ------------------------------------------------------------------------------------------------------------- */

  printf("\nFINISHED\n\n\t...writing model fit parameters\n\n");

  size_t prefixLength = strlen(argv[5]);

  char fout1[prefixLength+10];
  sprintf(fout1,"%s.S0_LM.nii",argv[5]);
  char fout2[prefixLength+14];
  sprintf(fout2,"%s.R2star_LM.nii",argv[5]);
  char fout3[prefixLength+8];
  sprintf(fout3,"%s.rsqr.nii",argv[5]);

  if (nifti_set_filenames(nim_outparam1,fout1,1,1)) return 1;
  nifti_image_write(nim_outparam1);
  nifti_image_free(nim_outparam1);
  if (nifti_set_filenames(nim_outparam2,fout2,1,1)) return 1;
  nifti_image_write(nim_outparam2);
  nifti_image_free(nim_outparam2);
  if (nifti_set_filenames(nim_rsqr,fout3,1,1)) return 1;
  nifti_image_write(nim_rsqr);
  nifti_image_free(nim_rsqr);

  clock_t end = clock();
  int time_spent = (int)((end - begin) / CLOCKS_PER_SEC);

  printf("\ntotal time taken is approximately %d minutes and %d seconds\n",time_spent/60, time_spent%60);

  return 0;

}

// END OF MAIN!




/* --------------------------------- 
   ----- FUNCTION DEFINITIONS ------
   --------------------------------- */

int exp_f(const gsl_vector *x, void *data, gsl_vector *f)
{
  size_t n = ((struct data *)data)->n;
  double *s = ((struct data *)data)->s;
  double *te = ((struct data *)data)->te;

  double S0 = gsl_vector_get(x,0);
  double R2star = gsl_vector_get(x,1);

  size_t i;

  for (i=0; i<n; i++)
    {
      /* Model: S = S0 * exp(-R2star * t) */
      double t = te[i];
      double Si = S0 * exp(-R2star * t);
      gsl_vector_set(f,i,Si-s[i]);
    }

  return GSL_SUCCESS;
}

int exp_df(const gsl_vector *x, void *data, gsl_matrix *J)
{
  size_t n = ((struct data *)data)->n;
  double *te = ((struct data *)data)->te;

  double S0 = gsl_vector_get(x,0);
  double R2star = gsl_vector_get(x,1);

  size_t i;

  for (i=0; i<n; i++)
    {
      double t = te[i];
      double e = exp(-R2star * t);
      gsl_matrix_set(J,i,0,e);
      gsl_matrix_set(J,i,1,-t*S0*e);
    }
  return GSL_SUCCESS;
}

float *mean_vol(float *a, int nx, int ny, int nz, int nt, float outputarray[])
{

  /* function that returns the mean of all volumes in 4D dataset */

  int nvox = nx * ny * nz;
  for (int i=0; i<nvox; i++)
    {
      outputarray[i] = 0.0;
    }

  for (int vol=0; vol<nt; vol++){
    int count=0;
    for (int k=0; k<nz; k++){
      for (int j=0; j<ny; j++){
	for (int i=0; i<nx; i++){
	  outputarray[count] += *(a+count+(nvox*vol))/nt;
	  count++;
	}
      }
    }
  }

  return outputarray;
}

float *ME_LMfit_function(float *yin, float *xin, float x_init[P], float weights[N], float params[P])
{

  /* wrapper function for GSL nonlinear solver */

  //const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver *sol;
  int status, info;
  size_t i;
  //const size_t n = N;
  //const size_t p = P;
  size_t n = N;
  size_t p = P;

  gsl_matrix *J = gsl_matrix_alloc(n,p);
  gsl_matrix *covar = gsl_matrix_alloc(p,p);

  double y[N], xx[N], X[P], W[N];
  for (int ii=0; ii<N; ii++)
    {
      y[ii] = (double) *(yin + ii);
      xx[ii] = (double) *(xin + ii);
      W[ii] = (double) *(weights + ii);
    }
  for (int ii=0; ii<P; ii++)
    {
      X[ii] = (double) *(x_init + ii);
    } 

  struct data d = {n, y, xx};
  gsl_multifit_function_fdf f;
  gsl_vector_view x = gsl_vector_view_array(X, p);
  gsl_vector_view w = gsl_vector_view_array(W, n);

  gsl_vector *res_f;
  double chi, chi0;

  /* const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0; */
  double xtol = 1e-8;
  double gtol = 1e-8;
  double ftol = 0.0;

  //printf("%f %lf",*yin,d.s);

  f.f = &exp_f;
  f.df = &exp_df; // set to NULL for finite difference Jacobian
  f.n = n;
  f.p =p;
  f.params = &d;

  sol = gsl_multifit_fdfsolver_alloc(T,n,p);

  /* for (int i=0; i<N; i++)
    {
      y[i] = (double) *(yin + i);
      x[i] = (double) *(xin + i);
      } */

  /* initialise solver with starting points and weights */
  gsl_multifit_fdfsolver_wset(sol, &f, &x.vector, &w.vector);

  /* compute initial residual norm */
  res_f = gsl_multifit_fdfsolver_residual(sol);
  chi0 = gsl_blas_dnrm2(res_f);

  /* solver the system with a maximum of 20 iterations */
  status = gsl_multifit_fdfsolver_driver(sol,20,xtol,gtol,ftol,&info);

  gsl_multifit_fdfsolver_jac(sol, J);
  gsl_multifit_covar(J, 0.0, covar);

  /* compute final residual norm */
  chi = gsl_blas_dnrm2(res_f);

  //printf("\n%f %f",gsl_vector_get(sol->x,0),gsl_vector_get(sol->x,1));

  params[0] = (float) gsl_vector_get(sol->x,0);
  params[1] = (float) gsl_vector_get(sol->x,1);

  //DO NOT FORGET TO FREE MEMORY
  gsl_multifit_fdfsolver_free(sol);
  gsl_matrix_free(covar);
  gsl_matrix_free(J);

  return params;

}

 float *ME_BFfit_function(float *yin, float *xin, int n, int ngrid, float p1[2], float p2[2], float params[P])
{

  /* function to do brute-force exponential fit */

  /* create parameter arrays */
  float param1_linspaceArray[ngrid], param2_linspaceArray[ngrid];
  float *param1_array = linspace(p1[0], p1[1], ngrid, param1_linspaceArray);
  float *param2_array = linspace(p2[0], p2[1], ngrid, param2_linspaceArray);

  /* brute-force fit */

  float p1_best = *param1_array;
  float p2_best = *param2_array;
  float rss_new, rss_old = NAN;
  float yhat[n];

  for (int p1=0; p1<ngrid; p1++){
    for (int p2=0; p2<ngrid; p2++){

      for (int i=0; i<n; i++)
	{
	  yhat[i] = (*(param1_array + p1)) * exp( -1.0 * (*(param2_array + p2)) * (*(xin + i)) );
	}
      rss_new = ressumsq(yin,yhat,n);
      if ( p1 == 0 && p2 == 0)
	rss_old = rss_new;
      if (rss_new < rss_old){
	p1_best = *(param1_array + p1);
	p2_best = *(param2_array + p2);
	rss_old = rss_new;
      }
    }
  }

  params[0] = p1_best;
  params[1] = p2_best;

  return params;

}

float *linspace(float a, float b, int n, float linspaceArray[])
{
  float stepSize = (b - a)/(n - 1);
  for (int i=0; i<n; i++)
    {
      linspaceArray[i] = a+i*stepSize;
    }
  return linspaceArray;
}

float ressumsq(float *y, float *yhat, int n)
{
  float diffsq[n];
  for (int i=0; i<n; i++)
    {
      diffsq[i] = (*(y+i) - *(yhat+i)) * (*(y+i) - *(yhat+i));
    }
  float rss = sum_array(diffsq,n);
  return rss;
}

float sum_array(float *a, int n)
{
  float sum=0.0;
  for (int i=0; i<n; i++)
    {
      sum += *(a+i);
    }
  return sum;
}

float calc_rsqr(float *y, float *x, float params[P], int n)
{
  float ymean=0.0, rss=0.0, tss=0.0, rsqr=0.0;
  for (int i=0; i<n; i++)
    {
      ymean += *(y + i) / n;
    }
  for (int i=0; i<n; i++)
    {
      tss += sqr(*(y + i) - ymean);
    }
  for (int i=0; i<n; i++)
    {
      rss += sqr(*(y + i) - (params[0] * exp(-1.0 * params[1] * *(x + i))));
    }
  rsqr = 1.0 - (rss/tss);
  return rsqr;
}
