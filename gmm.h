#ifndef GMM_H
#define GMM_H

typedef struct {
  double mean;
  double var;
  double coef;
} Gaussian;

typedef struct {
  int nmix;
  Gaussian *g;
  double *weight;
} GMM;

#endif

void Gaussian_init(Gaussian *g, double mean, double var);
double Gaussian_prob(double x, Gaussian *g);
void GMM_init_minmax(GMM *gmm, int nmix, double min, double max);
void GMM_init(GMM *gmm, int nmax, int n, double *data);
void GMM_free(GMM *gmm);
double GMM_prob(double x, GMM *gmm);
double GMM_train(int n, double *x, GMM *gmm, int iter_max);



