#include "basics.h"

string to_string_prec(const double a, const int n) {
    ostringstream out;
    out.precision(n);
    out << fixed << a;
    return out.str();
}

// bubble radius as a function of conformal time
double radius(double tau, double taun) {
    double R = -1.0;
    if (tau >= taun)
        R = tau-taun;
    return R;
}

// random real number in the range (x_1,x_2)
double randomreal(double x1, double x2, rgen &mt) {
    long double r01 = mt()/(1.0*mt.max());
    return (x1 + (x2-x1)*r01);
}

// linear interpolation
double interpolate(double x, vector<vector<double> > &y) {
    int n = y.size();
    if (x > y[n-1][0]) {
        cout << "Warning: the point lies above of the interpolation range." << endl;
        cout << x << "   " << y[n-1][0] << endl;
        return y[n-1][1];
    }
    if (x < y[0][0]) {
        cout << "Warning: the point lies below of the interpolation range." << endl;
        cout << x << "   " << y[0][0] << endl;
        return y[0][1];
    }
    double dx = y[1][0] - y[0][0];
    int jx = (int) ((x-y[0][0])/dx);
    if (jx < n-1) {
        return y[jx][1] + (y[jx+1][1] - y[jx][1])*(x - y[jx][0])/dx;
    }
    return y[n-1][1];
}

// finds the x for which y(x)=y for a growing function y(x)
double findrootG(double y, double dx, vector<vector<double> > &list) {
    int n = list.size();
    double xmin = list[0][0];
    double xmax = list[n-1][0];
    double x = (xmax+xmin)/2.0;
    while (xmax-xmin > dx) {
        if (interpolate(x, list) > y) {
            xmax = x;
        } else {
            xmin = x;
        }
        x = (xmax+xmin)/2.0;
    }
    return x;
}
