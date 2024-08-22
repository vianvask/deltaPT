#include <complex>
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <fstream>

using namespace std;

typedef mt19937_64 rgen;

const double PI = 3.141592653589793238463;

string to_string_prec(const double a, const int n);
double radius(double tau, double taun);
double randomreal(double x1, double x2, rgen &mt);
double interpolate(double x, vector<vector<double> > &y);
double findrootG(double y, double dx, vector<vector<double> > &list);
