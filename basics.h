#include <complex>
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <fstream>

using namespace std;

typedef mt19937_64 rgen;

const double PI = 3.141592653589793238463;
const double beta = 1.0;
const complex<double> I(0.0, 1.0);
const complex<double> zero(0.0, 0.0);

string to_string_prec(const double a, const int n);

class bubble {
    public:
        vector<double> x;
        double t;
        double tau;
};

int sgn(double a);
int imod(int i, int n);

double radius(double tau, double taun);
double randomreal(double x1, double x2, rgen &mt);
vector<double> sphererand(double R, rgen &mt);
vector<double> sphererandr(double r, rgen &mt);
vector<double> cubicrand(double L, rgen &mt);
double norm(const vector<double> &X);
double distance(const vector<double> &X, const vector<double> &Y);
double interpolate(double x, vector<vector<double> > &y);
double interpolate(double x, vector<double> &y);
vector<double> interpolate3(double x, vector<vector<double> > &y);
double findrootG(double y, double dx, vector<vector<double> > &list);
double findrootG(double y, double dx, vector<double> &list);
