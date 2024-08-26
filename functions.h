#include "basics.h"

vector<double> averageevolution(function<double(double)> Gamma, const double tmin, const int jtmax, const double dt, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht);

double findtk(double k, double tkmax, vector<vector<double> > &at, vector<vector<double> > &Ht);

vector<vector<double> > Nbark(function<double(double)> Gamma, const double k, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at);

vector<double> findtrange(function<double(double)> Gamma);

vector<vector<double> > Fk(vector<vector<double> > &Nk, vector<vector<vector<double> > > &pd, const double k, int J, vector<vector<double> > &tau);

vector<int> jtlist(vector<vector<double> > &Nk, int J, rgen &mt);

vector<vector<vector<double> > > ddist(function<double(double)> Gamma, const double k, int jdmax, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at);

double Vfrac(double tau, double tauj, double dj, double k);

vector<vector<double> > rhoevolution(vector<vector<double> > &Ft, vector<vector<double> > &Ht);
