#include <math.h>
#include <stdlib.h>
#include "gmm.h"


static double data_min(int n, double *x)
{
   double m = x[0];
   int i;
   for (i = 1; i < n; i++)
     if (m > x[i])
        m = x[i];
   return m;
}

static double data_max(int n, double *x)
{
   double m = x[0];
   int i;
   for (i = 1; i < n; i++)
     if (m < x[i])
        m = x[i];
   return m;
}

static void Gaussian_recalc(Gaussian *g)
{
  g->coef = 1.0/sqrt(g->var*6.2831853);
}

void Gaussian_init(Gaussian *g, double mean, double var)
{
  g->mean = mean;
  g->var = var;
  Gaussian_recalc(g);
}

double Gaussian_prob(double x, Gaussian *g)
{
  double y = x - g->mean;
  return(g->coef*exp(-y*y/g->var));
}

void GMM_init_minmax(GMM *gmm, int nmix, double min, double max)
{
  int i;
  gmm->nmix = nmix;
  gmm->g = (Gaussian*)malloc(sizeof(Gaussian)*nmix);
  gmm->weight = (double*)malloc(sizeof(double)*nmix);
  for (i = 0; i < nmix; i++) {
    Gaussian_init(&gmm->g[i], min+(max-min)*i/nmix, 1.0);
    gmm->weight[i] = 1.0/nmix;
  }
}

void GMM_init(GMM *gmm, int nmix, int n, double *data)
{
  GMM_init_minmax(gmm,nmix,data_min(n,data),data_max(n,data));
}

void GMM_free(GMM *gmm)
{
  free(gmm->g);
  free(gmm->weight);
}

double GMM_prob(double x, GMM *gmm)
{
  double p = 0.0;
  int i;
  for (i = 0; i < gmm->nmix; i++) {
    p += gmm->weight[i]*Gaussian_prob(x,&gmm->g[i]);
  }
  return p;
}

double GMM_reest(int n, double *x, GMM *gmm) 
{
   double *gamma;
   double *newmean;
   double *newvar;
   double *totalgamma;
   double total_logprob = 0.0;
   int i,j;
   gamma = (double*)malloc(sizeof(double)*gmm->nmix);
   totalgamma = (double*)malloc(sizeof(double)*gmm->nmix);
   newmean = (double*)malloc(sizeof(double)*gmm->nmix);
   newvar = (double*)malloc(sizeof(double)*gmm->nmix);
   for (i = 0; i < gmm->nmix; i++) {
     newmean[i] = newvar[i] = totalgamma[i] = 0.0;
   }
   for (i = 0; i < n; i++) { 
     double tgamma = 0.0;
     for (j = 0; j < gmm->nmix; j++) {
        gamma[j] = gmm->weight[j]*Gaussian_prob(x[i],&gmm->g[j]);
        tgamma += gamma[j];
     }
     total_logprob += log(tgamma);
     for (j = 0; j < gmm->nmix; j++)
        gamma[j] /= tgamma;
     for (j = 0; j < gmm->nmix; j++) {
        double y;
        newmean[j] += gamma[j]*x[i];
        y = x[i]-gmm->g[j].mean;
        newvar[j] += gamma[j]*y*y;
        totalgamma[j] += gamma[j];
     }
  }
  for (j = 0; j < gmm->nmix; j++) {
    Gaussian_init(&gmm->g[j],newmean[j]/totalgamma[j],newvar[j]/totalgamma[j]);
    gmm->weight[j] = totalgamma[j]/n;
  }
  free(gamma);
  free(newmean);
  free(newvar);
  free(totalgamma);
  return total_logprob;
} 

double GMM_train(int n, double *x, GMM *gmm, int iter_max)
{
   double l;
   int i;
   for (i = 0; i < iter_max; i++) {
     l = GMM_reest(n,x,gmm);
   }
   return l;
}

