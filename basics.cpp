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

// generates one point within radius R
vector<double> sphererand(double R, rgen &mt) {
    double r, ctheta, stheta, phi;
    vector<double> x(3);
    r = R*pow(randomreal(0.0,1.0,mt),1.0/3.0);
    ctheta = randomreal(-1.0,1.0,mt);
    stheta = sqrt(1.0-ctheta*ctheta);
    phi = randomreal(0,2.0*PI,mt);
    x[0] = r*cos(phi)*stheta;
    x[1] = r*sin(phi)*stheta;
    x[2] = r*ctheta;
    return x;
}

// generates one point at radius r
vector<double> sphererandr(double r, rgen &mt) {
    double ctheta, stheta, phi;
    vector<double> x(3);
    ctheta = randomreal(-1.0,1.0,mt);
    stheta = sqrt(1.0-ctheta*ctheta);
    phi = randomreal(0,2.0*PI,mt);
    x[0] = r*cos(phi)*stheta;
    x[1] = r*sin(phi)*stheta;
    x[2] = r*ctheta;
    return x;
}

// norm of a vector
double norm(const vector<double> &X) {
    double d2 = 0.0;
    for (int j = 0; j < X.size(); j++) {
        d2 += X[j]*X[j];
    }
    return sqrt(d2);
}

// distance between two points
double distance(const vector<double> &X, const vector<double> &Y) {
    double d2 = 0.0;
    for (int j = 0; j < X.size(); j++) {
        d2 += (X[j]-Y[j])*(X[j]-Y[j]);
    }
    return sqrt(d2);
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

// linear interpolation
double interpolate(double x, vector<double> &y) {
    int n = y.size();
    
    int jx = (int) x;
    if (jx < n-1) {
        return y[jx] + (y[jx+1] - y[jx])*(x-jx);
    }
    return y[n-1];
}

// linear interpolation
vector<double> interpolate3(double x, vector<vector<double> > &y) {
    int n = y.size();
    vector<double> r(2);
    if (x > y[n-1][0]) {
        cout << "Warning: the point lies outside of the interpolation range." << endl;
        r[0] = y[n-1][1];
        r[1] = y[n-1][2];
        return r;
    }
    if (x < y[0][0]) {
        cout << "Warning: the point lies outside of the interpolation range." << endl;
        r[0] = y[0][1];
        r[1] = y[0][2];
        return r;
    }
    double dx = y[1][0] - y[0][0];
    int jx = (int) ((x-y[0][0])/dx);
    if (jx < n-1) {
        r[0] = y[jx][1] + (y[jx+1][1] - y[jx][1])*(x - y[jx][0])/dx;
        r[1] = y[jx][2] + (y[jx+1][2] - y[jx][2])*(x - y[jx][0])/dx;
        return r;
    }
    r[0] = y[n-1][1];
    r[1] = y[n-1][2];
    return r;
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

// finds the position of y in a growing list
double findrootG(double y, double dx, vector<double> &list) {
    int n = list.size();
    double xmin = 0.0;
    double xmax = n-1.0;
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
