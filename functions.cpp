#include "functions.h"

// evolution of the Universe on average, returns {kmax,tkmax}
vector<double> averageevolution(function<double(double)> Gamma, const double tmin, const int jtmax, const double dt, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt) {
    
    // initial state in vacuum dominance:
    double H = 1.0;
    double t = tmin, a = exp(H*tmin), tau = (1.0 - exp(-H*tmin))/H;
    double rhoR = 0.0001*3.0*pow(H,2.0)/(8.0*PI);
    double rhoV0 = 3.0*pow(H,2.0)/(8.0*PI) - rhoR;
    double F = 1.0, F0 = 1.0;
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
        
        tmp[1] = rhoR;
        rhoRt.push_back(tmp);
        tmp[1] = rhoV0*F;
        rhoVt.push_back(tmp);
        
        // compute the Hubble rate, scale factor and conformal time
        H = sqrt(8.0*PI*(rhoV0*F + rhoR)/3.0);
        a += H*a*dt;
        tau += dt/a;
        
        if (a*H > kmax) {
            kmax = a*H;
            tkmax = t;
        }
        
        // update the false vacuum fraction
        Nt = 0.0;
        F0 = F;
        for (int j = 0; j < taut.size(); j++) {
            Nt += 4.0*PI/3.0*dt*Gamma(taut[j][0])*pow(at[j][1]*radius(tau,taut[j][1]), 3.0);
        }
        F = exp(-Nt);
        
        // update the radiation energy density
        rhoR += -rhoV0*(F-F0) - 4.0*H*rhoR*dt;
        
        t += dt;
    }
    
    tmp[0] = kmax;
    tmp[1] = tkmax;
    
    return tmp;
}

 // evolution of the total energy density
 void rhoevolutionFG(vector<vector<double> > &F, vector<vector<double> > &deltaF, vector<vector<double> > &phiF, vector<vector<double> > &phiB, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt, vector<vector<double> > &FW, vector<vector<double> > &N, vector<vector<vector<double> > > &pd, double k, int J, int jdmax, rgen &mt) {
     
     const double dt = at[1][0] - at[0][0];
     double rhoV0 = rhoVt[0][1];
     
     double t = Ht[0][0];
     vector<double> tmp(2);
     tmp[0] = t;
     tmp[1] = 1.0;
     F.push_back(tmp);
     
     tmp[1] = 0.0;
     deltaF.push_back(tmp);
     phiF.push_back(tmp);
     phiB.push_back(tmp);
     
     double deltarhoV = 0.0, deltarhoR = 0.0, deltarho = 0.0, deltaP = 0.0, deltaq = 0.0, Phi = 0.0, PhiB = 0.0;
     double tau, H, a, FS, rhoV, rhoVb, rhoRb, rho, P, deltarhoV0;
     vector<double> tauj, dj, rj;
     
     int jb = 0;
     for (int jt = 1; jt < at.size(); jt++) {
         t = taut[jt][0];
         tau = taut[jt][1];
         H = Ht[jt][1];
         a = at[jt][1];
         
         rhoVb = rhoVt[jt][1];
         rhoRb = rhoRt[jt][1];
         rho = rhoRb + rhoVb;
         P = rhoRb/3.0 - rhoVb;
         
         Phi = -4.0*PI*deltaq/H;
         
         // try to generate a bubble
         if (jb < J && sqrt(abs(1.0 + 2.0*Phi))*dt*N[jt][2] > randomreal(0.0,1.0,mt)) {
             tauj.push_back(taut[jt][1]);
             dj.push_back(findrootG(randomreal(0.0,pd[jt][jdmax-1][1],mt), 0.001, pd[jt]));
             rj.push_back(0.0);
             jb++;
         }
         
         // compute the false vacuum fraction
         FS = 1.0;
         for (int j = 0; j < jb; j++) {
             rj[j] += sqrt(abs(1.0 + 2.0*Phi))*dt/a;
             FS *= 1.0 - Vfrac(rj[j], dj[j], k);
         }
         rhoV = rhoV0*FS*FW[jt][1];
         
         deltarhoV0 = deltarhoV;
         deltarhoV = rhoV - rhoVb;
                  
         deltarhoR += -3.0*H*(deltarho+deltaP)*dt - (deltarhoV-deltarhoV0) + pow(k/a,2.0)*deltaq*dt + 4.0*PI*(rho+P)*(3.0*H*deltaq-deltarho)/H*dt;
         deltaq += -deltaP*dt - 3.0*H*deltaq*dt + 4.0*PI*(rho+P)*deltaq/H*dt;
         
         deltarho = deltarhoR + deltarhoV;
         deltaP = deltarhoR/3.0 - deltarhoV;
         
         PhiB = 4.0*PI*pow(a/k,2.0)*(3.0*H*deltaq-deltarho);
         
         tmp[0] = Ht[jt][0];
         tmp[1] = rhoV/rhoV0;
         F.push_back(tmp);
         
         tmp[1] = deltarho/rho;
         deltaF.push_back(tmp);
         
         tmp[1] = Phi;
         phiF.push_back(tmp);
         
         tmp[1] = PhiB;
         phiB.push_back(tmp);
     }
 }


// evolution of the total energy density
void rhoevolutionCG(vector<vector<double> > &F, vector<vector<double> > &deltaC, vector<vector<double> > &phiC, vector<vector<double> > &phiB, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt, vector<vector<double> > &FW, vector<vector<double> > &N, vector<vector<vector<double> > > &pd, double k, int J, int jdmax, rgen &mt) {
    
    const double dt = at[1][0] - at[0][0];
    double rhoV0 = rhoVt[0][1];
    
    double t = Ht[0][0];
    vector<double> tmp(2);
    tmp[0] = t;
    tmp[1] = 1.0;
    F.push_back(tmp);
    
    tmp[1] = 0.0;
    deltaC.push_back(tmp);
    phiC.push_back(tmp);
    phiB.push_back(tmp);
    
    double deltarhoV = 0.0, deltarhoR = 0.0, deltarho = 0.0, deltaP = 0.0, B = 0.0, Phi = 0.0, Psi = 0.0, PhiB = 0.0;
    double tau, H, a, FS, rhoV, rhoVb, rhoRb, rho, P, deltarhoV0, dPsi;
    vector<double> tauj, dj, rj;
    
    int jb = 0;
    for (int jt = 0; jt < at.size(); jt++) {
        t = taut[jt][0];
        tau = taut[jt][1];
        H = Ht[jt][1];
        a = at[jt][1];
        
        rhoVb = rhoVt[jt][1];
        rhoRb = rhoRt[jt][1];
        rho = rhoRb + rhoVb;
        P = rhoRb/3.0 - rhoVb;
        
        Phi = -deltaP/(P+rho);
        
        // try to generate a bubble
        if (jb < J && sqrt(abs(1.0 + 2.0*Phi))*dt*N[jt][2] > randomreal(0.0,1.0,mt)) {
            tauj.push_back(taut[jt][1]);
            dj.push_back(findrootG(randomreal(0.0,pd[jt][jdmax-1][1],mt), 0.001, pd[jt]));
            rj.push_back(0.0);
            jb++;
        }
        
        // compute the false vacuum fraction
        FS = 1.0;
        for (int j = 0; j < jb; j++) {
            rj[j] += sqrt(abs(1.0 + 2.0*Phi))*dt/a;
            FS *= 1.0 - Vfrac(rj[j], dj[j], k);
        }
        rhoV = rhoV0*FS*FW[jt][1];
        
        deltarhoV0 = deltarhoV;
        deltarhoV = rhoV - rhoVb;
                 
        B = (4.0*PI*pow(a,2.0)*deltarho + pow(k,2.0)*Psi)/(pow(k,2.0)*a*H);
        dPsi = -H*Phi;
        deltarhoR += -3.0*H*(deltarho+deltaP)*dt - (deltarhoV-deltarhoV0) + (rho+P)*(3.0*dPsi - pow(k,2.0)*B/a)*dt;
        Psi += dPsi*dt;
        
        deltarho = deltarhoR + deltarhoV;
        deltaP = deltarhoR/3.0 - deltarhoV;
        
        PhiB = Psi - a*H*B;
        
        tmp[0] = t;
        tmp[1] = rhoV/rhoV0;
        F.push_back(tmp);
        
        tmp[1] = deltarho/rho;
        deltaC.push_back(tmp);
        
        tmp[1] = Phi;
        phiC.push_back(tmp);
        
        tmp[1] = PhiB;
        phiB.push_back(tmp);
    }
}


// evolution of the total energy density
void rhoevolutionNG(vector<vector<double> > &F, vector<vector<double> > &deltaN, vector<vector<double> > &phiN, vector<vector<double> > &phiB, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt, vector<vector<double> > &FW, vector<vector<double> > &N, vector<vector<vector<double> > > &pd, double k, int J, int jdmax, rgen &mt) {
    
    const double dt = at[1][0] - at[0][0];
    double rhoV0 = rhoVt[0][1];
    
    double t = Ht[0][0];
    vector<double> tmp(2);
    tmp[0] = t;
    tmp[1] = 1.0;
    F.push_back(tmp);
    
    tmp[1] = 0.0;
    deltaN.push_back(tmp);
    phiN.push_back(tmp);
    phiB.push_back(tmp);
    
    double deltarhoV = 0.0, deltarhoR = 0.0, deltarho = 0.0, deltaP = 0.0, v = 0.0, Phi = 0.0, PhiB = 0.0;
    double tau, H, a, FS, rhoV, rhoVb, rhoRb, rho, P, deltarhoV0, dPhi;
    vector<double> tauj, dj, rj;
    
    int jb = 0;
    for (int jt = 0; jt < at.size(); jt++) {
        t = taut[jt][0];
        tau = taut[jt][1];
        H = Ht[jt][1];
        a = at[jt][1];

        rhoVb = rhoVt[jt][1];
        rhoRb = rhoRt[jt][1];
        rho = rhoRb + rhoVb;
        P = rhoRb/3.0 - rhoVb;
        
        // try to generate a bubble
        if (jb < J && sqrt(abs(1.0 + 2.0*Phi))*dt*N[jt][2] > randomreal(0.0,1.0,mt)) {
            tauj.push_back(taut[jt][1]);
            dj.push_back(findrootG(randomreal(0.0,pd[jt][jdmax-1][1],mt), 0.001, pd[jt]));
            rj.push_back(0.0);
            jb++;
        }
        
        // compute the false vacuum fraction
        FS = 1.0;
        for (int j = 0; j < jb; j++) {
            rj[j] += sqrt(abs(1.0 + 2.0*Phi))*dt/a;
            FS *= 1.0 - Vfrac(rj[j], dj[j], k);
        }
        rhoV = rhoV0*FS*FW[jt][1];
        
        deltarhoV0 = deltarhoV;
        deltarhoV = rhoV - rhoVb;
        
        v = (4.0*PI*pow(a,2.0)*deltarho + pow(k,2.0)*Phi)/(12.0*PI*pow(a,3.0)*H*(rho+P));
        dPhi = -(4.0*PI*deltarho + (pow(k/a,2.0)+3.0*pow(H,2.0))*Phi)/(3.0*H);
        
        deltarhoR += -3.0*H*(deltarho+deltaP)*dt - (deltarhoV-deltarhoV0) + (rho+P)*(3.0*dPhi + pow(k,2.0)*v/a)*dt;
        Phi += dPhi*dt;
        
        deltarho = deltarhoR + deltarhoV;
        deltaP = deltarhoR/3.0 - deltarhoV;
        
        PhiB = Phi;
        
        tmp[0] = t;
        tmp[1] = rhoV/rhoV0;
        F.push_back(tmp);
        
        tmp[1] = deltarho/rho;
        deltaN.push_back(tmp);
        
        tmp[1] = Phi;
        phiN.push_back(tmp);
        
        tmp[1] = PhiB;
        phiB.push_back(tmp);
    }
}




// time of horizon reentry of scale k
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


// computes the characteristic bubble radius
double Rstar(function<double(double)> Gamma, vector<vector<double> > &Ft, vector<vector<double> > &at, double tp) {
    const double dt = at[1][0] - at[0][0];
    
    int jt = 0;
    double Rstar = 0.0;
    while (at[jt][0] < tp) {
        Rstar += dt*Ft[jt][1]*Gamma(at[jt][0])*pow(at[jt][1],3.0);
        jt++;
    }
    Rstar += (tp - at[jt-1][0])*Ft[jt][1]*Gamma(at[jt][0])*pow(at[jt][1],3.0);
    Rstar = pow(Rstar/pow(at[jt-1][1] + (tp - at[jt-1][0])*(at[jt][1]-at[jt-1][1]),3.0), -1.0/3.0);
    return Rstar;
}


// finds the time range where the computation should be performed
vector<double> findtrange(function<double(double)> Gamma, double Nbarmin, double Fmin) {
    vector<double> trange(2);
    vector<vector<double> > Ft, taut, at, Ht, rhoRt, rhoVt;
    vector<double> tmp(2);
    
    int jtmax = 6000;
    double dt = 0.001;
    tmp = averageevolution(Gamma, -3.0, jtmax, dt, Ft, taut, at, Ht, rhoRt, rhoVt);
        
    double kmax = tmp[0];
        
    vector<vector<double> > Nk = Nbark(Gamma, kmax, Ft, taut, at);
        
    vector<vector<double> > Nt(jtmax, vector<double> (2,0.0));
    for (int jt = 0; jt < jtmax; jt++) {
        Nt[jt][0] = Nk[jt][0];
        Nt[jt][1] = Nk[jt][1];
    }
    
    vector<vector<double> > Tt(jtmax, vector<double> (2,0.0));
    for (int jt = 0; jt < jtmax; jt++) {
        Tt[jt][0] = Ft[jt][0];
        Tt[jt][1] = 1-Ft[jt][1];
    }
    
    trange[0] = findrootG(Nbarmin, dt, Nt);
    trange[1] = findrootG(1.0-Fmin, dt, Tt);
    
    return trange;
}

// the false vacuum fraction neglecting the first J bubbles
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
                    Nt += dt*dd*Nk[j][2]*pd[j][jd][2]*Vfrac(radius(tau,tauj),pd[j][jd][0],k);
                }
            }
        }
        tmp[0] = t;
        tmp[1] = exp(-Nt);
        F.push_back(tmp);
    }
    return F;
}

// CDF of nucleation distances
vector<vector<vector<double> > > ddist(function<double(double)> Gamma, const double k, int jdmax, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at) {
    
    const double dt = at[1][0] - at[0][0];
    vector<vector<double> > Nk = Nbark(Gamma, k, Ft, taut, at);
    
    double dmax = 1.0/k + radius(taut[taut.size()-1][1],taut[0][1]);
    double dd = dmax/(1.0*(jdmax-1));
    
    double tau, d, p, ptot;
    vector<vector<double> > tmp(jdmax, vector<double> (3));
    vector<vector<vector<double> > > Cd(taut.size(), vector<vector<double> > (jdmax, vector<double> (2)));
    for (int jt = 0; jt < taut.size(); jt++) {
        tau = taut[jt][1];
        
        ptot = 0.0;
        for (int jd = 0; jd < jdmax; jd++) {
            d = dd*jd;
            
            // integrate over time up to t
            p = 0.0;
            for (int j = 0; j < jt; j++) {
                if (1.0/k + max(0.0,radius(tau,taut[j][1])) > d) {
                    p += 4.0*PI*dt*Gamma(taut[j][0])*pow(at[j][1],3.0)*pow(d,2.0)/Nk[jt][1];
                }
            }
            ptot += p;
            tmp[jd][0] = d;
            tmp[jd][1] = dd*ptot;
            tmp[jd][2] = p;
        }
        Cd[jt] = tmp;
    }
    return Cd;
}

// volume of intersection of two bubbles separated by distance d
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

double Vfrac(double rj, double dj, double k) {
    double Vk = 4.0*PI/3.0*pow(k,-3.0);
    return Vint(dj, max(0.0,dj-1.0/k)+rj, 1.0/k)/Vk;
}


// surface area of intersection of two bubbles separated by distance d
double Sint(double d, double R, double r) {
    if (R+r > d && abs(R-r) < d) {
        return PI*R*(pow(r,2.0)-pow(R-d,2.0))/d;
    }
    if (R+r > d && R-r >= d) {
        return 4.0*PI*pow(r,2.0);
    }
    if (R+r > d && r-R >= d) {
        return 4.0*PI*pow(R,2.0);
    }
    return 0.0;
}

double rhowall(double rj, double dj, double k) {
    double Vj = 4.0*PI/3.0*pow(dj,3.0);
    double Vk = 4.0*PI/3.0*pow(k,-3.0);
    return Sint(dj, max(0.0,dj-1.0/k)+rj, 1.0/k)*Vj/Vk;
}
