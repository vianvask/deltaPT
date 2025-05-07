#include "basics.h"

vector<double> averageevolution(function<double(double)> Gamma, const double tmin, const int jtmax, const double dt, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt);

void rhoevolutionFG(vector<vector<double> > &F, vector<vector<double> > &deltaF, vector<vector<double> > &phiF, vector<vector<double> > &phiB, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt, vector<vector<double> > &FW, vector<vector<double> > &N, vector<vector<vector<double> > > &pd, double k, int J, int jdmax, rgen &mt);

void rhoevolutionCG(vector<vector<double> > &F, vector<vector<double> > &deltaC, vector<vector<double> > &phiC, vector<vector<double> > &phiB, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt, vector<vector<double> > &FW, vector<vector<double> > &N, vector<vector<vector<double> > > &pd, double k, int J, int jdmax, rgen &mt);

void rhoevolutionNG(vector<vector<double> > &F, vector<vector<double> > &deltaN, vector<vector<double> > &phiN, vector<vector<double> > &phiB, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt, vector<vector<double> > &FW, vector<vector<double> > &N, vector<vector<vector<double> > > &pd, double k, int J, int jdmax, rgen &mt);


double findtk(double k, double tkmax, vector<vector<double> > &at, vector<vector<double> > &Ht);

vector<vector<double> > Nbark(function<double(double)> Gamma, const double k, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at);

double Rstar(function<double(double)> Gamma, vector<vector<double> > &Ft, vector<vector<double> > &at, double tp);

vector<double> findtrange(function<double(double)> Gamma, double Nbarmin, double Fmin);

vector<vector<double> > Fk(vector<vector<double> > &Nk, vector<vector<vector<double> > > &pd, const double k, int J, vector<vector<double> > &tau);

vector<vector<vector<double> > > ddist(function<double(double)> Gamma, const double k, int jdmax, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at);

double Vfrac(double rj, double dj, double k);
