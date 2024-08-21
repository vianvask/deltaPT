#include "functions.h"

// evolution of the Universe on average, returns kmax
vector<double> averageevolution(function<double(double)> Gamma, const double tmin, const int jtmax, const double dt, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht) {
    
    // initial state in vacuum dominance:
    double t = tmin, a = exp(tmin), tau = 1.0 - exp(-tmin);
    double rhoV = 1.0, rhoV0 = 1.0, rhoR = 0.0;
    
    double H = 1.0, F = 1.0, F0 = 1.0;
    vector<double> tmp(2);
    double Nt, kmax = 0.0, tkmax;
    
    for (int jt = 0; jt < jtmax; jt++) {
        tmp[0] = t;
        tmp[1] = F;
        Ft.push_back(tmp);
        tmp[1] = H;
        Ht.push_back(tmp);
        tmp[1] = a;
        at.push_back(tmp);
        tmp[1] = tau;
        taut.push_back(tmp);
        
        //calculate the false vacuum fraction
        Nt = 0.0;
        F0 = F;
        for (int j = 0; j < taut.size(); j++) {
            Nt += 4.0*PI/3.0*dt*Gamma(taut[j][0])*pow(at[j][1]*radius(tau,taut[j][1]), 3.0);
        }
        F = exp(-Nt);
        
        // compute the scale factor and conformal time
        H = sqrt(rhoV + rhoR);
        a += dt*H*a;
        tau += dt/a;
        
        if (a*H > kmax) {
            kmax = a*H;
            tkmax = t;
        }
        
        // update the energy densities
        rhoV = F;
        rhoR += -4.0*H*rhoR*dt - (F-F0);
        
        t += dt;
    }
    
    tmp[0] = kmax;
    tmp[1] = tkmax;
    
    return tmp;
}

// evolution of the Universe on average, returns kmax
vector<vector<double> > rhoevolution(vector<vector<double> > &Ft, vector<vector<double> > &Ht) {
    
    // initial state in vacuum dominance:
    const double dt = Ft[1][0] - Ft[0][0];
    
    vector<double> tmp(2);
    double rhoV, rhoR = 0.0;
    vector<vector<double> > rho;
    for (int jt = 0; jt < Ft.size(); jt++) {
        rhoV = Ft[jt][1];
        if (jt>0) {
            rhoR += -4.0*Ht[jt][1]*rhoR*dt - (Ft[jt][1]-Ft[jt-1][1]);
        }
        
        tmp[0] = Ft[jt][0];
        tmp[1] = rhoV + rhoR;
        
        rho.push_back(tmp);
    }
    
    return rho;
}

// evolution of the Universe on average, returns kmax
double findtk(double k, double tkmax, vector<vector<double> > &at, vector<vector<double> > &Ht) {
    
    vector<vector<double> > aH;
    vector<double> tmp(2);
    for (int jt = 0; jt < at.size(); jt++) {
        if (at[jt][0] > tkmax) {
            if (at[jt][1]*Ht[jt][1] < k) {
                return at[jt][0];
            }
        }
    }
    return at[at.size()-2][0];
}

// expected number of bubbles in sphere of radius 1/k
vector<vector<double> > Nbark(function<double(double)> Gamma, const double k, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at) {
    
    const double dt = at[1][0] - at[0][0];
    
    double t, tau, Np, N = 0.0;
    vector<double> tmp(3);
    vector<vector<double> > Nt;
    for (int jt = 0; jt < at.size(); jt++) {
        t = at[jt][0];
        tau = taut[jt][1];
        
        Np = N;
        N = 0.0;
        for (int j = 0; j < jt; j++) {
            N += 4.0*PI/3.0*dt*Gamma(taut[j][0])*pow(at[j][1]*(1.0/k + radius(tau,taut[j][1])),3.0);
        }
        
        tmp[0] = t;
        tmp[1] = N;
        tmp[2] = (N-Np)/dt;
        Nt.push_back(tmp);
    }
    
    return Nt;
}

// worse approximation of the false vacuum fraction neglecting the first J bubbles
vector<vector<double> > Fk(function<double(double)> Gamma, const double k, int J, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at) {
    
    const double dt = at[1][0] - at[0][0];
    vector<vector<double> > Nk = Nbark(Gamma, k, Ft, taut, at);
    
    double t, tau, Nt;
    vector<double> tmp(2);
    vector<vector<double> > F;
    for (int jt = 0; jt < at.size(); jt++) {
        t = at[jt][0];
        tau = taut[jt][1];
        
        // integrate the region where Nbark > J
        Nt = 0.0;
        for (int j = 0; j < taut.size(); j++) {
            if (Nk[j][1] > J && taut[j][1] < tau) {
                Nt += 4.0*PI/3.0*dt*Gamma(taut[j][0])*pow(at[j][1]*radius(tau,taut[j][1]), 3.0);
            }
        }
        tmp[0] = t;
        tmp[1] = exp(-Nt);
        F.push_back(tmp);
    }
    return F;
}

// better approximation of the false vacuum fraction neglecting the first J bubbles
vector<vector<double> > Fk(vector<vector<double> > &Nk, vector<vector<vector<double> > > &pd, const double k, int J, vector<vector<double> > &taut) {
    
    const double dt = taut[1][0] - taut[0][0];
    const double dd = pd[0][1][0] - pd[0][0][0];
    
    double t, tau, tauj, Nt;
    vector<double> tmp(2);
    vector<vector<double> > F;
    for (int jt = 0; jt < taut.size(); jt++) {
        tau = taut[jt][1];
        
        // integrate the region where Nbark > J
        Nt = 0.0;
        for (int j = 0; j < taut.size(); j++) {
            tauj = taut[j][1];
            if (Nk[j][1] > J && tau > tauj) {
                for (int jd = 0; jd < pd[0].size(); jd++) {
                    Nt += dt*dd*Nk[j][2]*pd[j][jd][2]*Vfrac(tau,tauj,pd[j][jd][0],k);
                }
            }
        }
        tmp[0] = t;
        tmp[1] = exp(-Nt);
        F.push_back(tmp);
    }
    return F;
}

// generate times tj for j<J
vector<int> jtlist(vector<vector<double> > &Nk, int J, rgen &mt) {
    const double dt = Nk[1][0] - Nk[0][0];
    vector<int> jtlist(J, Nk.size() - 1);
    int j = 0;
    for (int jt = 0; jt < Nk.size(); jt++) {
        if (dt*Nk[jt][2] > 1.0) {
            cout << "Warning: too long nucleation timestep, N = " << Nk[jt][1] << "." << endl;
        }
        if (dt*Nk[jt][2] > randomreal(0.0,1.0,mt)) {
            jtlist[j] = jt;
            j++;
        }
        if (j >= J) {
            return jtlist;
        }
    }
    cout << "Warning: nucleation finished at j = " << j << "." << endl;
    return jtlist;
}

// PDF of the nucleation times of the fitst j<J bubbles
vector<vector<double> > ptj(function<double(double)> Gamma, const double k, int J, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at) {
    
    const double dt = at[1][0] - at[0][0];
    vector<vector<double> > Nk = Nbark(Gamma, k, Ft, taut, at);
    
    double p, ptot;
    vector<double> tmp(taut.size());
    vector<vector<double> > pt;
    for (int j = 1; j <= J; j++) {
        ptot = 0.0;
        for (int jt = 0; jt < taut.size(); jt++) {
            p = Nk[jt][2]*pow(Nk[jt][1],j-1.0)/tgamma(j)*exp(-Nk[jt][1]);
            ptot += p;
            tmp[jt] = p;
        }
        pt.push_back(tmp);
    }
    return pt;
}

// CDF of nucleation distances
vector<vector<vector<double> > > pd(function<double(double)> Gamma, const double k, int jdmax, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at) {
    
    const double dt = at[1][0] - at[0][0];
    vector<vector<double> > Nk = Nbark(Gamma, k, Ft, taut, at);
    
    double dmax = 1.0/k + radius(taut[taut.size()-1][1],taut[0][1]);
    double dd = dmax/(1.0*(jdmax-1));
    
    double tau, d, p, ptot;
    vector<vector<double> > tmp(jdmax, vector<double> (3));
    vector<vector<vector<double> > > pd(taut.size(), vector<vector<double> > (jdmax, vector<double> (2)));
    for (int jt = 0; jt < taut.size(); jt++) {
        tau = taut[jt][1];
        
        ptot = 0.0;
        for (int jd = 0; jd < jdmax; jd++) {
            d = dd*jd;
            
            // integrate over time up to t
            p = 0.0;
            for (int j = 0; j < jt; j++) {
                if (1.0/k + max(0.0,radius(tau,taut[j][1])) > d) {
                    p += 4.0*PI*dt*Gamma(taut[j][0])*pow(at[j][1], 3.0)*pow(d,2.0)/Nk[jt][1];
                }
            }
            ptot += p;
            tmp[jd][0] = d;
            tmp[jd][1] = dd*ptot;
            tmp[jd][2] = p;
        }
        pd[jt] = tmp;
    }
    return pd;
}

double Vint(double d, double R, double r) {
    if (R+r > d && abs(R-r) < d) {
        return PI*pow(R+r-d,2.0)*(d*d + 2.0*d*r - 3.0*r*r + 2.0*d*R + 6.0*r*R - 3.0*R*R)/(12.0*d);
    }
    if (R+r > d && R-r >= d) {
        return 4.0*PI/3.0*pow(r,3.0);
    }
    if (R+r > d && r-R >= d) {
        return 4.0*PI/3.0*pow(R,3.0);
    }
    return 0.0;
}

double Vfrac(double tau, double tauj, double dj, double k) {
    double Vk = 4.0*PI/3.0*pow(k,-3.0);
    if (tau > tauj) {
        return Vint(dj, max(0.0,dj-1.0/k)+radius(tau,tauj), 1.0/k)/Vk;
    }
    return 0.0;
}

vector< vector<int> > getAllSubsets(vector<int> &set){
    vector<vector<int> > subset;
    vector<int> empty;
    subset.push_back(empty);

    vector< vector<int> > subsetTemp;
    for (int i = 0; i < set.size(); i++) {
        subsetTemp = subset;
        
        for (int j = 0; j < subsetTemp.size(); j++)
            subsetTemp[j].push_back( set[i] );

        for (int j = 0; j < subsetTemp.size(); j++)
            subset.push_back( subsetTemp[j] );
    }
    return subset;
}
