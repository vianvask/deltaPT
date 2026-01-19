#include "basics.h"

vector<double> averageevolution(function<double(double)> Gamma, const double tmin, const int jtmax, const double dt, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &ttau, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt);

void rhoevolutionFG(vector<vector<double> > &F, vector<vector<double> > &deltaF, vector<vector<double> > &phiF, vector<vector<double> > &phiB, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt, vector<vector<double> > &FW, vector<vector<double> > &N, vector<vector<vector<double> > > &pd, double k, int J, int jdmax, rgen &mt);

void rhoevolutionCG(vector<vector<double> > &F, vector<vector<double> > &deltaC, vector<vector<double> > &phiC, vector<vector<double> > &phiB, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt, vector<vector<double> > &FW, vector<vector<double> > &N, vector<vector<vector<double> > > &pd, double k, int J, int jdmax, rgen &mt);

void rhoevolutionNG(vector<vector<double> > &F, vector<vector<double> > &deltaN, vector<vector<double> > &phiN, vector<vector<double> > &phiB, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt, vector<vector<double> > &FW, vector<vector<double> > &N, vector<vector<vector<double> > > &pd, double k, int J, int jdmax, rgen &mt);

void rhoevolutionNG_new(vector<vector<double> >& F, vector<vector<double> >& deltaC, vector<vector<double> >& deltaN, vector<vector<double> >& phiN, vector<vector<double> >& phiB, vector<vector<double> >& vN, vector<vector<double> >& Boundary, vector<vector<double> >& R, vector<vector<double> >& LaplphiN, vector<vector<double> >& taut, vector<vector<double> >& at, vector<vector<double> >& Ht, vector<vector<double> >& rhoRt, vector<vector<double> >& rhoVt, vector<vector<double> >& FW, vector<vector<double> >& N, vector<vector<vector<double> > >& pd, double k, int J, int jdmax, rgen& mt, vector<vector<vector<double> > >& gamma_w_vec, vector<vector<double> >& boundaryJ);

void rhoevolutionNG_new2(vector<vector<double> >& F, vector<vector<double> > &deltaC, vector<vector<double> >& deltaN, vector<vector<double> >& phiN, vector<vector<double> >& phiB, vector<vector<double> > &vN, vector<vector<double> > &Boundary, vector<vector<double> > &R, vector<vector<double> > &LaplphiN, vector<vector<double> >& deltaN_Y, vector<vector<double> >& phiN_Y, vector<vector<double> >& phiB_Y, vector<vector<double> >& vN_Y, vector<vector<double> >& Boundary_Y, vector<vector<double> >& R_Y, vector<vector<double> >& deltaPnad_Y, vector<vector<double> >& taut, vector<vector<double> >& at, vector<vector<double> >& Ht, vector<vector<double> >& rhoRt, vector<vector<double> >& rhoVt, vector<vector<double> >& FW, vector<vector<double> >& N, vector<vector<vector<double> > >& pd, double k, int J, int jdmax, rgen& mt, vector<vector<vector<double> > >& gamma_w_vec, vector<vector<double> >& boundaryJ);


double findtk(double k, double tkmax, vector<vector<double> > &at, vector<vector<double> > &Ht);

vector<vector<double> > Nbark(function<double(double)> Gamma, const double k, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at);

double Rstar(function<double(double)> Gamma, vector<vector<double> > &Ft, vector<vector<double> > &at, double tp);

vector<double> findtrange(function<double(double)> Gamma, double Nbarmin, double Fmin);

vector<vector<double> > Fk(function<double(double)> Gamma, vector<vector<double> >& Nk,  const double k, int J, vector<vector<double> >& taut, vector<vector<double> >& at);

vector<vector<double> > Fk2(vector<vector<double> >& Nk, vector<vector<vector<double> > >& pd, const double k, int J, vector<vector<double> >& taut);

vector<vector<vector<double> > > gamma_w_fun(vector<vector<double> >& taut, vector<vector<double> >& at, vector<vector<double> >& Ht, vector<vector<double>>& R_n_vs_R_H, const double rhoV0);

vector<double> gamma_c_at_t(const int jtmax, const int jtn, vector<vector<double> >& taut, vector<vector<double> >& at, vector<vector<double> >& Ht, vector<vector<vector<double> > >& gamma_w_vec);

vector<vector<vector<Row> > > gamma_w_c_fun(vector<vector<double> >& taut, vector<vector<double> >& at, vector<vector<double> >& Ht, vector<vector<vector<double> > > &gamma_w_vec);

vector<vector<double> > boundary_term_fun(function<double(double)> Gamma, const double k, int J, vector<vector<double> >& Nk, vector<vector<double> >& at, vector<vector<double>>& Ht, vector<vector<double> >& taut, vector<vector<double> >& Ft, vector<vector<vector<double> > >& gamma_w_vec);

vector<vector<double> > boundary_term_fun2(const double k, int J, vector<vector<double> >& Nk, vector<vector<vector<double> > >& pd, vector<vector<double> >& at, vector<vector<double>>& Ht, vector<vector<double> >& taut, vector<vector<double> >& Ft, vector<vector<vector<double> > >& gamma_w_vec);

vector<vector<vector<double> > > ddist(function<double(double)> Gamma, const double k, int jdmax, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &ttau);

vector<vector<vector<double> > > ddist2(function<double(double)> Gamma, const double k, int jdmax, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at);

double Vfrac(double rj, double dj, double k);

double dArea(double rj, double dj, double k);
