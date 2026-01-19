#include "functions.h"

// evolution of the Universe on average, returns {kmax,tkmax}
vector<double> averageevolution(function<double(double)> Gamma, const double tmin, const int jtmax, const double dt, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &ttau, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &rhoRt, vector<vector<double> > &rhoVt) {
    
    // initial state in vacuum dominance:
    double H = 1.0;
    double t = tmin, a = exp(H*tmin), tau = (1.0 - exp(-H*tmin))/H;
    double rhoR = 0.0001*3.0*pow(H,2.0)/(8.0*PI);
    double rhoV0 = 3.0*pow(H,2.0)/(8.0*PI) - rhoR; // G -> 1
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

        ttau.emplace_back(tmp.rbegin(), tmp.rend());
        
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
        tmp[1] = rhoV / rhoV0;
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
    for (int jt = 1; jt < at.size(); jt++) { // changed starting point from 0 to 1, because vector size mismatch and 2 values at 0 otherwise
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
        
        PhiB = Phi; // Bardeen potential, it's simply equal to Phi in NG
        
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

void rhoevolutionNG_new(vector<vector<double> >& F, vector<vector<double> >& deltaC, vector<vector<double> >& deltaN, vector<vector<double> >& phiN, vector<vector<double> >& phiB, vector<vector<double> >& vN, vector<vector<double> >& Boundary, vector<vector<double> >& R, vector<vector<double> >& LaplphiN, vector<vector<double> >& taut, vector<vector<double> >& at, vector<vector<double> >& Ht, vector<vector<double> >& rhoRt, vector<vector<double> >& rhoVt, vector<vector<double> >& FW, vector<vector<double> >& N, vector<vector<vector<double> > >& pd, double k, int J, int jdmax, rgen& mt, vector<vector<vector<double> > >& gamma_w_vec, vector<vector<double> >& boundaryJ) {

    const double dt = at[1][0] - at[0][0];
    double rhoV0 = rhoVt[0][1];

    double t = Ht[0][0];
    vector<double> tmp(2);
    vector<double> tmpFnotj(2); //[F, dF]
    tmp[0] = t;
    tmp[1] = 1.0;
    F.push_back(tmp);

    tmp[1] = 0.0;
    deltaN.push_back(tmp);
    phiN.push_back(tmp);
    phiB.push_back(tmp);
    vN.push_back(tmp);
    R.push_back(tmp);
    LaplphiN.push_back(tmp);
    Boundary.push_back(tmp);


    deltaC.push_back(tmp);


    double deltarhoV = 0.0, deltarhoR = 0.0, deltarho = 0.0, deltaP = 0.0, Phi = 0.0, PhiB = 0.0, LaplPhi = 0.0;
    double deltarho_C = 0; // delta in comoving gauge, computed with new method
    double cs, deltaPnad, PhiY = 0.0, dPhidtY = 0.0, deltaY = 0.0, vY, ddPhidtY;
    double tau, H, a, FS, rhoV, rhoVb, rhoRb, rho, P, deltarhoV0, dPhi, boundary, boundary_temp, Fnotj0, v;
    double dR, R_n, dA;
    vector<double> tauj, dj, rj;
    vector<int> jtn; // iterators for nucleation times
    vector<vector<double> > gamma_w_c(J); // we store Lorentz factor for all possible collision times, for all J bubbles we simulate - gamma[jb][jc][gamma(t, tc(jc), tn(jtn(jb))]
    // we iterate jc s.t jc = 0 corresponds to t = tn(jtn) + dt (first possible collision time for bubble nucleated at tn)
    vector<vector<vector<double> > > dFnotj(J); // change of volume fraction of all the bubbles other than jth, \dot F_k^{\neq j}(t) * dt

    // initialize this guy manually, bc the loop starts at jt = 1 (not jt = 0)
    tmpFnotj[0] = 1;
    tmpFnotj[1] = 0.;
    for (int j = 0; j < J; j++) {
        dFnotj[j].reserve(at.size());
        dFnotj[j].emplace_back(tmpFnotj);
    }


    //ofstream outfile, outfileBub;
    //outfile.open("boundary_full_5.000_gammaperbeta_0.000_13.dat");
   // outfileBub.open("bubbles_5.000_gammaperbeta_0.000_13.dat");

    int jb = 0;
    for (int jt = 1; jt < at.size(); jt++) { // changed starting point from 0 to 1, because vector size mismatch and 2 values for the 0th element otherwise
        t = taut[jt][0];
        tau = taut[jt][1];
        H = Ht[jt][1];
        a = at[jt][1];

        boundary = boundaryJ[jt][1]; // add the contribution from j>j_c bubbles

        //cout << t << "   " << jb << endl;

        // background quantities
        rhoVb = rhoVt[jt][1];
        rhoRb = rhoRt[jt][1];
        rho = rhoRb + rhoVb;
        P = rhoRb / 3.0 - rhoVb;

        // compute \dot F_k^{\neq j}(t) * dt - worth checking if it's OK
        if (jb == 0) {
            tmpFnotj[0] = 1;
            tmpFnotj[1] = 0.;
            for (int j = 0; j < J; j++) {
                dFnotj[j].emplace_back(tmpFnotj);
            }
        }
        else {
            for (int j = 0; j < J; j++) {
                Fnotj0 = dFnotj[j].back()[0];
                tmpFnotj[0] = 1;
                for (int i = 0; i < jb; i++) {
                    if (i != j) {
                        tmpFnotj[0] *= 1. - Vfrac(rj[i], dj[i], k);
                    }
                }
                tmpFnotj[1] = tmpFnotj[0] - Fnotj0; //dF
                dFnotj[j].emplace_back(tmpFnotj);
            }
        }

        // compute the gamma factor for all the bubbles nucleated before jt
        for (int j = 0; j < jb; j++) {
            dR = taut[jt][1] - taut[jt - 1][1]; // we change dt to dR, solving for gamma(t(R))
            for (int jc = 0; jc < gamma_w_c[j].size(); jc++) { // evolve the Lorentz factor for the walls that has collided before jt
                R_n = taut[jt][1] - taut[jtn[j]][1]; // radius of bubble nucleated at t_n
                gamma_w_c[j][jc] -= dR * (gamma_w_c[j][jc] / R_n) * (2. + 3. * R_n * at[jt][1] * Ht[jt][1]);
            }
            gamma_w_c[j].emplace_back(gamma_w_vec[jtn[j]][jt - jtn[j]][4]); // walls colliding precisely at jt
        }

        //cout << "boundary start\n";

        // compute the boundary term - we can do it now cause bubble nucleating at jt wouldn't influence it anyway
        for (int j = 0; j < jb; j++) {
            dA = dArea(rj[j], dj[j], k);
            //cout << "dA = " <<  dA << endl;
            if (dA != 0.) {
                boundary_temp = dFnotj[j][jt][0] * gamma_w_vec[jtn[j]][jt - jtn[j]][4]; // uncollided part
                //cout << "boundary = " << boundary_temp << endl;
                for (int jc = 0; jc < gamma_w_c[j].size(); jc++) {
                    boundary_temp -= dFnotj[j][jc + jtn[j] + 1][1] * gamma_w_c[j][jc];
                }
                boundary_temp *= dA * sqrt(abs(1.0 + 2.0 * Phi)) * pow(a, -2.0);
                //cout << "boundary = " << boundary_temp << endl;
            }
            else {
                boundary_temp = 0;
            }
            boundary += boundary_temp;
        }

        //cout << "boundary end\n";

        //outfile << t << "   " << boundary;

        // try to generate a bubble
        if (jb < J && sqrt(abs(1.0 + 2.0 * Phi)) * dt * N[jt][2] > randomreal(0.0, 1.0, mt)) {
            tauj.push_back(taut[jt][1]);
            dj.push_back(findrootG(randomreal(0.0, pd[jt][jdmax - 1][1], mt), 0.001, pd[jt]));
            rj.push_back(0.0);
            if (dj[jb] <= 1. / k) { // bubble nucleated inside
                jtn.push_back(jt);
            }
            else { // bubble nucleated outside
                jtn.push_back(finditerG(tauj[jb] - dj[jb] + 1. / k, taut));
                if (jtn[jb] == jt) { // if finditerG() screws things up
                    jtn[jb] -= 1;
                }
                gamma_w_c[jb].reserve(at.size() - jtn[jb] - 1); // reserve space for gamma_w_c, -1 because we forbid collision at nucleation time
                gamma_w_c[jb] = gamma_c_at_t(jt, jtn[jb], taut, at, Ht, gamma_w_vec);
            }
            //outfileBub << jt << "   " << dj[jb] << endl;
            jb++;
        }

        // compute the false vacuum fraction
        FS = 1.0;
        for (int j = 0; j < jb; j++) {
            rj[j] += sqrt(abs(1.0 + 2.0 * Phi)) * dt / a;
            FS *= 1.0 - Vfrac(rj[j], dj[j], k);
        }
        rhoV = rhoV0 * FS * FW[jt][1];

        deltarhoV0 = deltarhoV; // in previous timestep
        deltarhoV = rhoV - rhoVb;

        // full system in NG
        dPhi = dt * (-4 * PI * deltarho - 3 * H * H * Phi + pow(a, -2.0) * LaplPhi) / (3 * H);

        deltarhoR += -(deltarhoV - deltarhoV0) + 3.0 * (rho + P) * dPhi + (boundary - 3 * H * (deltarho + deltaP)) * dt;
        Phi += dPhi;
        LaplPhi += (4 * PI * a * a * boundary - H * LaplPhi) * dt;

        // delta_C in a smart way
        deltarho_C += (boundary - 3. * H * deltarho_C) * dt;

        deltarho = deltarhoR + deltarhoV;
        deltaP = deltarhoR / 3.0 - deltarhoV;
        v = (-dPhi / dt - H * Phi) / (4 * PI * a * (rho + P));

        PhiB = Phi; // Bardeen potential, it's simply equal to Phi in NG

        tmp[0] = t;
        tmp[1] = rhoV / rhoV0;
        F.push_back(tmp);

        //outfile << "   " << rhoV / rhoV0 << endl;

        tmp[1] = deltarho / rho;
        deltaN.push_back(tmp);

        tmp[1] = Phi;
        phiN.push_back(tmp);

        tmp[1] = PhiB;
        phiB.push_back(tmp);

        tmp[1] = v;
        vN.push_back(tmp);

        tmp[1] = Phi - a * H * v;
        R.push_back(tmp);

        tmp[1] = boundary;
        Boundary.push_back(tmp);

        tmp[1] = LaplPhi;
        LaplphiN.push_back(tmp);


        tmp[1] = deltarho_C / rho;
        deltaC.push_back(tmp);

    }
    //outfile.close();
    //outfileBub.close();
}

void rhoevolutionNG_new2(vector<vector<double> >& F, vector<vector<double> >& deltaC, vector<vector<double> >& deltaN, vector<vector<double> >& phiN, vector<vector<double> >& phiB, vector<vector<double> >& vN, vector<vector<double> >& Boundary, vector<vector<double> >& R, vector<vector<double> >& LaplphiN,  vector<vector<double> >& deltaN_Y, vector<vector<double> >& phiN_Y, vector<vector<double> >& phiB_Y, vector<vector<double> >& vN_Y, vector<vector<double> >& Boundary_Y, vector<vector<double> >& R_Y, vector<vector<double> >& deltaPnad_Y, vector<vector<double> >& taut, vector<vector<double> >& at, vector<vector<double> >& Ht, vector<vector<double> >& rhoRt, vector<vector<double> >& rhoVt, vector<vector<double> >& FW, vector<vector<double> >& N, vector<vector<vector<double> > >& pd, double k, int J, int jdmax, rgen& mt, vector<vector<vector<double> > >& gamma_w_vec, vector<vector<double> >& boundaryJ) {

    const double dt = at[1][0] - at[0][0];
    double rhoV0 = rhoVt[0][1];

    double t = Ht[0][0];
    vector<double> tmp(2);
    vector<double> tmpFnotj(2); //[F, dF]
    tmp[0] = t;
    tmp[1] = 1.0;
    F.push_back(tmp);

    tmp[1] = 0.0;
    deltaN.push_back(tmp);
    phiN.push_back(tmp);
    phiB.push_back(tmp);
    vN.push_back(tmp);
    R.push_back(tmp);
    LaplphiN.push_back(tmp);
    Boundary.push_back(tmp);

    deltaN_Y.push_back(tmp);
    phiN_Y.push_back(tmp);
    phiB_Y.push_back(tmp);
    vN_Y.push_back(tmp);
    R_Y.push_back(tmp);
    deltaPnad_Y.push_back(tmp);
    Boundary_Y.push_back(tmp);

    deltaC.push_back(tmp);

    vector<double> EoS_w(2); // [w(t-dt), w(t)] - for Yann's system of eqns (to compute c_s^2)
    EoS_w[0] = (rhoRt[0][1] / 3. - rhoVt[0][1]) / (rhoRt[0][1] + rhoVt[0][1]);


    double deltarhoV = 0.0, deltarhoR = 0.0, deltarho = 0.0, deltaP = 0.0, Phi = 0.0, PhiB = 0.0, LaplPhi = 0.0;
    double deltarho_C = 0; // delta in comoving gauge, computed with new method
    double cs, deltaPnad, PhiY = 0.0, dPhidtY = 0.0, deltaY = 0.0, vY, ddPhidtY;
    double tau, H, a, FS, rhoV, rhoVb, rhoRb, rho, P, deltarhoV0, dPhi, boundary, boundary_temp, Fnotj0, v;
    double dR, R_n, gamma0, dA;
    vector<double> tauj, dj, rj;
    vector<int> jtn; // iterators for nucleation times
    vector<vector<double> > gamma_w_c(J); // we store Lorentz factor for all possible collision times, for all J bubbles we simulate - gamma[jb][jc][gamma(t, tc(jc), tn(jtn(jb))]
                                          // we iterate jc s.t jc = 0 corresponds to t = tn(jtn) + dt (first possible collision time for bubble nucleated at tn)
    vector<vector<vector<double> > > dFnotj(J); // change of volume fraction of all the bubbles other than jth, \dot F_k^{\neq j}(t) * dt

    // initialize this guy manually, bc the loop starts at jt = 1 (not jt = 0)
    tmpFnotj[0] = 1;
    tmpFnotj[1] = 0.;
    for (int j = 0; j < J; j++) {
        dFnotj[j].push_back(tmpFnotj);
    }


    //ofstream outfile, outfileBub;
    //outfile.open("boundary_full_5.000_gammaperbeta_0.000_13.dat");
   // outfileBub.open("bubbles_5.000_gammaperbeta_0.000_13.dat");

    int jb = 0;
    for (int jt = 1; jt < at.size(); jt++) { // changed starting point from 0 to 1, because vector size mismatch and 2 values for the 0th element otherwise
        t = taut[jt][0];
        tau = taut[jt][1];
        H = Ht[jt][1];
        a = at[jt][1];

        boundary = boundaryJ[jt][1]; // add the contribution from j>j_c bubbles

        //cout << t << "   " << jb << endl;

        // background quantities
        rhoVb = rhoVt[jt][1];
        rhoRb = rhoRt[jt][1];
        rho = rhoRb + rhoVb;
        P = rhoRb / 3.0 - rhoVb;

        // compute \dot F_k^{\neq j}(t) * dt - worth checking if it's OK
        if (jb == 0) {
            tmpFnotj[0] = 1;
            tmpFnotj[1] = 0.;
            for (int j = 0; j < J; j++) {
                dFnotj[j].push_back(tmpFnotj);
            }
        }
        else {
            for (int j = 0; j < J; j++) {
                Fnotj0 = dFnotj[j].back()[0];
                tmpFnotj[0] = 1;
                for (int i = 0; i < jb; i++) {
                    if (i != j) {
                        tmpFnotj[0] *= 1. - Vfrac(rj[i], dj[i], k);
                    }
                }
                tmpFnotj[1] = tmpFnotj[0] - Fnotj0; //dF
                dFnotj[j].push_back(tmpFnotj);
            }
        }

        // compute the gamma factor for all the bubbles nucleated before jt
        for (int j = 0; j < jb; j++) {
            dR = taut[jt][1] - taut[jt - 1][1]; // we change dt to dR, solving for gamma(t(R))
            for (int jc = 0; jc < gamma_w_c[j].size(); jc++) { // evolve the Lorentz factor for the walls that has collided before jt
                R_n = radius(taut[jt][1], taut[jtn[j]][1]); // radius of bubble nucleated at t_n
                gamma0 = gamma_w_c[j][jc];
                gamma_w_c[j][jc] -= dR * (2. * gamma0 / R_n) * (1 + 3 * R_n * at[jt][1] * Ht[jt][1] / 2);
            }
            gamma_w_c[j].push_back(gamma_w_vec[jtn[j]][jt - jtn[j]][4]); // walls colliding precisely at jt
        }

        //cout << "boundary start\n";

        // compute the boundary term - we can do it now cause bubble nucleating at jt wouldn't influence it anyway
        for (int j = 0; j < jb; j++) {
            dA = dArea(rj[j], dj[j], k);
            //cout << "dA = " <<  dA << endl;
            if (dA != 0.) {
                boundary_temp = dFnotj[j][jt][0] * gamma_w_vec[jtn[j]][jt - jtn[j]][4]; // uncollided part
                //cout << "boundary = " << boundary_temp << endl;
                for (int jc = 0; jc < gamma_w_c[j].size(); jc++) {
                    boundary_temp -= dFnotj[j][jc + jtn[j] + 1][1] * gamma_w_c[j][jc];
                }
                boundary_temp *= dA * sqrt(abs(1.0 + 2.0 * Phi)) * pow(a, -2.0);
                //cout << "boundary = " << boundary_temp << endl;
            }
            else {
                boundary_temp = 0;
            }
            boundary += boundary_temp;
        }

        //cout << "boundary end\n";

        //outfile << t << "   " << boundary;

        // try to generate a bubble
        if (jb < J && sqrt(abs(1.0 + 2.0 * Phi)) * dt * N[jt][2] > randomreal(0.0, 1.0, mt)) {
            tauj.push_back(taut[jt][1]);
            dj.push_back(findrootG(randomreal(0.0, pd[jt][jdmax - 1][1], mt), 0.001, pd[jt]));
            rj.push_back(0.0);
            if (dj[jb] <= 1. / k) { // bubble nucleated inside
                jtn.push_back(jt);
            }
            else { // bubble nucleated outside
                jtn.push_back(finditerG(tauj[jb] - dj[jb] + 1. / k, taut));
                if (jtn[jb] == jt) { // if finditerG() screws things up
                    jtn[jb] -= 1;
                }
                gamma_w_c[jb] = gamma_c_at_t(jt, jtn[jb], taut, at, Ht, gamma_w_vec);
            }
            //outfileBub << jt << "   " << dj[jb] << endl;
            jb++;
        }

        // compute the false vacuum fraction
        FS = 1.0;
        for (int j = 0; j < jb; j++) {
            rj[j] += sqrt(abs(1.0 + 2.0 * Phi)) * dt / a;
            FS *= 1.0 - Vfrac(rj[j], dj[j], k);
        }
        rhoV = rhoV0 * FS * FW[jt][1];

        deltarhoV0 = deltarhoV; // in previous timestep
        deltarhoV = rhoV - rhoVb;

        // full system in NG
        dPhi = dt * (-4 * PI * deltarho - 3*H*H * Phi + pow(a, -2.0) * LaplPhi) / (3 * H);

        deltarhoR += -(deltarhoV - deltarhoV0) + 3.0 * (rho + P) * dPhi + (boundary - 3 * H * (deltarho + deltaP)) * dt;
        Phi += dPhi;
        LaplPhi += (4 * PI * a * a * boundary - H * LaplPhi) * dt;

        // delta_C in a smart way
        deltarho_C += (boundary - 3. * H * deltarho_C) * dt;

        //Yann's system
        EoS_w[1] = (rhoRb / 3. - rhoVb) / (rhoRb + rhoVb); // new value of w
        cs = EoS_w[1] - (EoS_w[1] - EoS_w[0]) / (dt * 3 * H * (1 + EoS_w[1])); // sound speed (squared)
        deltaPnad = (1. / 3. - cs) * deltaY - (4. / 3.) * deltarhoV / (rhoRb + rhoVb);
        ddPhidtY = dt * (-4 * H * dPhidtY + 3 * H * H * EoS_w[1] * PhiY + (3./2.) * H * H * (cs * deltaY + deltaPnad)); // separate, to increment all values in the same way
        
        deltaY += dt * (boundary / (rhoRb + rhoVb) - 3 * H * (deltaPnad + (cs - EoS_w[1]) * deltaY) + 3 * (1 + EoS_w[1]) * dPhidtY);
        PhiY += dt * dPhidtY;
        dPhidtY += ddPhidtY;

        EoS_w[0] = EoS_w[1]; // replace old value by new one, which becomes the old for the next timestep
        vY = (-dPhidtY - H * PhiY) / (4 * PI * a * (rho + P));

        deltarho = deltarhoR + deltarhoV;
        deltaP = deltarhoR / 3.0 - deltarhoV;
        v = (-dPhi / dt - H * Phi) / (4 * PI * a * (rho + P));

        PhiB = Phi; // Bardeen potential, it's simply equal to Phi in NG

        tmp[0] = t;
        tmp[1] = rhoV / rhoV0;
        F.push_back(tmp);

        //outfile << "   " << rhoV / rhoV0 << endl;

        tmp[1] = deltarho / rho;
        deltaN.push_back(tmp);

        tmp[1] = Phi;
        phiN.push_back(tmp);

        tmp[1] = PhiB;
        phiB.push_back(tmp);

        tmp[1] = v;
        vN.push_back(tmp);

        tmp[1] = Phi - a * H * v;
        R.push_back(tmp);

        tmp[1] = boundary;
        Boundary.push_back(tmp);

        tmp[1] = LaplPhi;
        LaplphiN.push_back(tmp);


        tmp[1] = deltarho_C / rho;
        deltaC.push_back(tmp);


        tmp[1] = deltaY;
        deltaN_Y.push_back(tmp);

        tmp[1] = PhiY;
        phiN_Y.push_back(tmp);

        tmp[1] = PhiY;
        phiB_Y.push_back(tmp);

        tmp[1] = vY;
        vN_Y.push_back(tmp);

        tmp[1] = PhiY - a * H * vY;
        R_Y.push_back(tmp);

        tmp[1] = boundary;
        Boundary_Y.push_back(tmp);

        tmp[1] = deltaPnad;
        deltaPnad_Y.push_back(tmp);
    }
    //outfile.close();
    //outfileBub.close();
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
    vector<vector<double> > Ft, taut, at, Ht, rhoRt, rhoVt, ttau;
    vector<double> tmp(2);
    
    int jtmax = 6000;
    double dt = 0.001;
    tmp = averageevolution(Gamma, -3.0, jtmax, dt, Ft, taut, ttau, at, Ht, rhoRt, rhoVt);
        
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

// the false vacuum fraction neglecting the first J bubbles - integral over nucleation distance being done analytically
vector<vector<double> > Fk(function<double(double)> Gamma, vector<vector<double> >& Nk, const double k, int J, vector<vector<double> >& taut, vector<vector<double> >& at) {
    const double dt = taut[1][0] - taut[0][0];

    double t, tau, tauj, Nt, jt_jc, tau_jc;
    vector<double> tmp(2);
    vector<vector<double> > F;
    bool jc_reached = 0;
    double Nt0;
    vector<double> term(6); // auxilliary terms to simplify integration
    for (int jt = 0; jt < taut.size(); jt++) {
        tau = taut[jt][1];
        t = taut[jt][0];

        // integrate the region where Nbark > J
        Nt = 0.0;
        if (jc_reached == 0) {
            if (Nk[jt][1] >= J) {
                jc_reached = 1;
                jt_jc = jt;
                tau_jc = taut[jt_jc][1];
                // cout << taut[taut.size() - 1][1] - tau_jc << "   " << k * (taut[taut.size() - 1][1] - tau_jc) << endl; // check how much bubbles nucleated at jc grow until the end of the simulation
            }
        }
        else if (tau - taut[jt_jc][1] <= 1 / k) { // different integration regimes due to different behavior of V_int
            for (int j = 0; j < jt_jc; j++) { // bubbles nucleated before t_jc
                tauj = taut[j][1];
                term[0] = pow(k, 3.0) * pow(tau, 3.0);
                term[1] = k * tau_jc * (6 + k * tau_jc * (6 + k * tau_jc));
                term[2] = 3 * pow(k, 3.0) * tau * tau * (tau_jc - 2 * tauj);
                term[3] = 6 * k * tauj * (1 + k * tau_jc) * (4 + k * tau_jc);
                term[4] = 6 * k * k * tauj * tauj * (4 + k * tau_jc);
                term[5] = 3 * k * tau * (6 + k * (6 * tau_jc + k * tau_jc * tau_jc - 2 * tauj * (3 + k * tauj)));
                Nt += dt * Gamma(taut[j][0]) * pow(at[j][1], 4.0) * (-1 / 24) * PI * pow(tau - tau_jc, 3.0) * (term[0] - term[1] + term[2] + term[3] - term[4] - term[5]);
                if (Nt < 0) {
                    cout << "A   " << Nt << endl;
                }
            }
            for (int j = jt_jc; j < jt; j++) { // bubbles nucleated after t_jc
                tauj = taut[j][1];
                Nt0 = Nt;
                Nt += dt * Gamma(taut[j][0]) * pow(at[j][1], 4.0) * (-1 / 24) * PI * k * pow(tau - tauj, 4.0) * (-18 + k * k * (tau - tauj) * (tau - tauj)); // bubbles nucleated outside
                if (Nt < Nt0) {
                    cout << "B   " << Nt << endl;
                }
                Nt0 = Nt;
                Nt += dt * Gamma(taut[j][0]) * pow(at[j][1], 4.0) * ((4 * PI / 3) * pow(1 / k - tau + tauj, 3.0) * pow(k, 3.0) * pow(tau - tauj, 3.0) + (PI / 8) * k * pow(tau - tauj, 4.0) * (26 + k * (tau - tauj) * (-32 + 11 * k * (tau - tauj)))); // bubbles nucleated inside
                if (Nt < Nt0) {
                    cout << "C   " << Nt << endl;
                }
            }
            //cout << Nt << endl;
        }
        else {
            cout << " KURWA VITTU\n"; // if you see this message, better use Fk2
            Nt = 100;
        }
        tmp[0] = t;
        tmp[1] = exp(-Nt);
        F.push_back(tmp);
    }
    return F;

}

// the false vacuum fraction neglecting the first J bubbles - 2d integral - varied precision of integration over pd
vector<vector<double>> Fk2(vector<vector<double>>& Nk, vector<vector<vector<double>>>& pd, const double k, int J, vector<vector<double>>& taut)
{
    const int dd_prec = 10; // determines how precisely pd is probed, make sure it's a divisor of jdmax (defined in main)
    const double dt = taut[1][0] - taut[0][0];
    const double dd = dd_prec * (pd[0][1][0] - pd[0][0][0]);

    double t, tau, tauj, Nt, t_jc;
    vector<double> tmp(2);
    vector<vector<double> > F;
    bool jc_reached = 0;
    for (int jt = 0; jt < taut.size(); jt++) {
        tau = taut[jt][1];
        t = taut[jt][0];

        // integrate the region where Nbark > J
        Nt = 0.0;
        if (jc_reached == 0) {
            if (Nk[jt][1] >= J) {
                jc_reached = 1;
                t_jc = jt;
            }
        }
        else {
            for (int j = t_jc; j < jt; j++) {
                tauj = taut[j][1];
                for (int jd = 0; jd < pd[0].size(); jd += dd_prec) {
                    Nt += dt * dd * Nk[j][2] * pd[j][jd][2] * Vfrac(radius(tau, tauj), pd[j][jd][0], k);
                }
            }
        }

        tmp[0] = t;
        tmp[1] = exp(-Nt);
        F.push_back(tmp);
    }
    return F;
}

// the evolution of wall gamma factor (times sigma) before collision obtained by solving (2.5) from 2305.04924
// we approximate the \sqrt(1-1/gamma^2) term by 1, since it's only relevant for R~R_H, when gamma >> 1
// we additionally assume P_{fric} << \Delta V and we set \Delta V to 1 for the purpose of computing gamma*sigma
vector<vector<vector<double>>> gamma_w_fun(vector<vector<double>>& taut, vector<vector<double>>& at, vector<vector<double>>& Ht, vector<vector<double>>& R_n_vs_R_H, const double rhoV0)
{
    vector<vector<vector<double> > > gamma_w_tmp; // gamma_w[jt_n][jt][t, t_n, tau, tau_n, gamma]
    vector<vector<vector<double> > > gamma_w_vec; // gamma_w[jt][jt_n][t, t_n, tau, tau_n, gamma]
    vector< vector<double> > tmp2d;
    vector<double> tmp(5), tmpR(5);
    double t, t_n, tau, tau_n, gamma, gamma0, dR, R_n;
    bool bR_n_vs_R_H; // find when bubble nucleated at t_n grows to R_H;
    // we want to run first over all t > t_n keeping t_n fixed
    for (int jt_n = 0; jt_n < taut.size()-2; jt_n++) {
        tmp2d.clear();
        bR_n_vs_R_H = 0;
        t_n = taut[jt_n][0];
        tau_n = taut[jt_n][1];
        t = t_n;
        tau = tau_n;

        // we put i.c. manually for t = t_n (+dt)
        // we use gamma(R) -> a*R*DeltaV/(2*sigma) as R->0
        gamma = 0; // R(t_n;t_n)=0

        tmp[0] = t;
        tmp[1] = t_n;
        tmp[2] = tau;
        tmp[3] = tau_n;
        tmp[4] = gamma;

        tmp2d.push_back(tmp);

        t = taut[jt_n + 1][0];
        tau = taut[jt_n + 1][1];
        gamma = rhoV0 * radius(taut[jt_n + 1][1], taut[jt_n][1]) * at[jt_n + 1][1] / 2.;

        tmp[0] = t;
        tmp[2] = tau;
        tmp[4] = gamma;

        tmp2d.push_back(tmp);

        tmpR[0] = t_n;
        tmpR[2] = tau_n;

        for (int jt = jt_n + 2; jt < taut.size(); jt++) {
            dR = taut[jt][1] - taut[jt - 1][1]; // we change dt to dR, solving for gamma(t(R))
            R_n = radius(taut[jt][1], taut[jt_n][1]); // radius of bubble nucleated at t_n;

            tmp[0] = taut[jt][0];
            tmp[2] = taut[jt][1];

            gamma0 = gamma;
            gamma += dR * (at[jt][1] * rhoV0 - (2 * gamma0 / R_n) * (1 + 3 * R_n * at[jt][1] * Ht[jt][1] / 2));
            tmp[4] = gamma;

            tmp2d.push_back(tmp);

            if (R_n * at[jt][1] * Ht[jt][1] > 1 && bR_n_vs_R_H == 0) { //
                bR_n_vs_R_H = 1;
                tmpR[1] = taut[jt][0];
                tmpR[3] = taut[jt][1];
                tmpR[4] = gamma;

                R_n_vs_R_H.push_back(tmpR);
            }
        }
        gamma_w_tmp.push_back(tmp2d);

    }
    // we add last two values of jt_n manually
    tmp2d.clear();


    tmp[0] = taut[taut.size() - 2][0]; // t
    tmp[1] = tmp[0]; // t_n
    tmp[2] = taut[taut.size() - 2][1]; // tau
    tmp[3] = tmp[2]; // tau_n
    tmp[4] = 0; // gamma

    tmp2d.push_back(tmp);

    tmp[0] = taut[taut.size() - 1][0]; // t
    tmp[2] = taut[taut.size() - 1][1]; // tau
    tmp[4] = rhoV0 * radius(taut[taut.size() - 1][1], taut[taut.size() - 2][1]) * at[taut.size() - 1][1] / 2; // gamma
    tmp2d.push_back(tmp);

    gamma_w_tmp.push_back(tmp2d);

    tmp2d.clear();


    tmp[1] = tmp[0];
    tmp[3] = tmp[2];
    tmp[4] = 0;

    tmp2d.push_back(tmp);
    gamma_w_tmp.push_back(tmp2d);


    return gamma_w_tmp;
}

// gamma(jtmax, jtn, jtc), where jtc represents the collision time and jtmax gives the t at which gamma is measured
vector<double> gamma_c_at_t(const int jtmax, const int jtn, vector<vector<double>>& taut, vector<vector<double>>& at, vector<vector<double>>& Ht, vector<vector<vector<double>>>& gamma_w_vec)
{
    vector<double> gamma_c;
    double gamma, gamma0, dR, R_n;

    // gamma_c.push_back(0.); // no collision happens at t_c = t_n

    // loop over all possible collision times
    for (int jtc = jtn + 1; jtc <= jtmax; jtc++) {
        gamma = gamma_w_vec[jtn][jtc - jtn][4]; // begin with gamma at the moment of collision
        // evolve gamma past the collision according to (2.5) from 2305.04924 with 0 at the r.h.s.
        for (int jt = jtc + 1; jt <= jtmax; jt++) {

            dR = taut[jt][1] - taut[jt - 1][1]; // we change dt to dR, solving for gamma(t(R))
            R_n = radius(taut[jt][1], taut[jtn][1]); // radius of bubble nucleated at t_n
            gamma0 = gamma;
            gamma -= dR * (2. * gamma0 / R_n) * (1 + 3 * R_n * at[jt][1] * Ht[jt][1] / 2);
        }
        gamma_c.push_back(gamma);
    }

    return gamma_c;
}

// the computation of Lorentz factor after bubble collision at time t_n < t_n_c <= t
// obtained by solving (2.5) from 2305.04924 with 0 at the r.h.s.
vector<vector<vector<Row>>> gamma_w_c_fun(vector<vector<double>>& taut, vector<vector<double>>& at, vector<vector<double>>& Ht, vector<vector<vector<double>>>& gamma_w_vec)
{
    //change to gamma_w_c[jt_n][jt_c][[jt][t,t_n,t_c,tau,tau_n,tau_c,gamma*sigma]
    vector<vector<vector<Row> > > gamma_w_c_tmp; // gamma_w_c[jt_n][jt][[jt_c][t,t_n,t_c,tau,tau_n,tau_c,gamma*sigma]
    vector<vector<Row> > tmp3d;
    vector<Row> tmp2d;
    Row tmp;
    double gamma0, dR, R_n;

    gamma_w_c_tmp.reserve(gamma_w_vec.size() - 1);

    size_t jt_n = 1280;

    //for (size_t jt_n = 0; jt_n < gamma_w_vec.size() - 1; jt_n++) {
        tmp3d.clear();
        tmp3d.reserve(gamma_w_vec.size() - jt_n);

        tmp.t_n = gamma_w_vec[jt_n][0][1];
        tmp.tau_n = gamma_w_vec[jt_n][0][3];

        cout << tmp.t_n << endl;


        for (size_t jt_c = 1; jt_c < gamma_w_vec[jt_n].size(); jt_c++) { // t_c happens only after t_n, hence we start from jt_c = 1
            tmp2d.clear();
            tmp2d.reserve(gamma_w_vec[jt_n].size() - jt_c);

            tmp.t_c = gamma_w_vec[jt_n][jt_c][0];
            tmp.tau_c = gamma_w_vec[jt_n][jt_c][2];
            tmp.gamma = gamma_w_vec[jt_n][jt_c][4]; // we start with gamma at the moment of collision

            cout << "\r" << tmp.t_c << "    " << flush;

            // add 0th element here
            tmp.t = tmp.t_c;
            tmp.tau = tmp.tau_c;

            tmp2d.push_back(tmp);

            for (int jt = jt_c + 1; jt < gamma_w_vec[jt_n].size(); jt++) {
                tmp.t = gamma_w_vec[jt_n][jt][0];
                tmp.tau = gamma_w_vec[jt_n][jt][2];
                
                dR = taut[jt + jt_n][1] - taut[jt + jt_n - 1][1]; // we change dt to dR, solving for gamma(t(R))
                R_n = radius(taut[jt + jt_n][1], taut[jt_n][1]); // radius of bubble nucleated at t_n
                gamma0 = tmp.gamma;
                tmp.gamma -= dR * (2. * gamma0 / R_n) * (1 + 3 * R_n * at[jt + jt_n][1] * Ht[jt + jt_n][1] / 2);


                tmp2d.push_back(tmp);
            }


            tmp3d.push_back(tmp2d);
        }

        gamma_w_c_tmp.push_back(tmp3d);
    //}


    return gamma_w_c_tmp;
}

vector<vector<double>> boundary_term_fun(function<double(double)> Gamma, const double k, int J, vector<vector<double>>& Nk, vector<vector<double>>& at, vector<vector<double>>& Ht, vector<vector<double>>& taut, vector<vector<double>>& Ft, vector<vector<vector<double>>>& gamma_w_vec)
{
    vector<vector<double>> boundary_vec(taut.size());
    vector<double> tmp(2);
    double Rtt0, Rttn, Rt0tn; // R(t; t_*), R(t; t_n), R(t_*; t_n)
    double gamma0, dR;
    double dt = taut[1][0] - taut[0][0];
    double energy_sum;
    vector<vector<double> > dF;

    tmp[0] = Ft[0][0];
    tmp[1] = 0;
    dF.push_back(tmp);

    for (int jtF = 1; jtF < Ft.size(); jtF++) {
        tmp[0] = Ft[jtF][0];
        tmp[1] = Ft[jtF][1] - Ft[jtF - 1][1];
        dF.push_back(tmp);
    }

    // boundary term = 0 before J = j_c is reached
    tmp[1] = 0;
    int jt = 0;
    while (Nk[jt][1] < J) {
        tmp[0] = Nk[jt][0];
        boundary_vec[jt] = tmp;
        jt++;
    }
    int jt0 = jt; // point to the start of approximation

    // initialize the rest of boundary_vec with 0's
    while (jt < taut.size()) {
        tmp[0] = Nk[jt][0];
        boundary_vec[jt] = tmp;
        jt++;
    }

    // it is better to perform large loop over the nucleation time since it's easier to control the evolution of Lorentz factor
    // when the nucleation time remains fixed
    vector<double> gamma_w_c;
    for (int jt_n = 0; jt_n < jt0; jt_n++) {
        gamma_w_c = gamma_c_at_t(jt0, jt_n, taut, at, Ht, gamma_w_vec); // gamma(t_*, t_n, t_n < t_c <= t_*)[jtc] (jtc = 0 for t_c = t_n + dt or jt_n + 1)
        Rt0tn = radius(taut[jt0][1], taut[jt_n][1]);

        //jt = jt0 gives 0 bc R(t, t_*) = 0
        for (jt = jt0 + 1; jt < taut.size(); jt++) {

            dR = taut[jt][1] - taut[jt - 1][1];
            Rttn = radius(taut[jt][1], taut[jt_n][1]);
            Rtt0 = radius(taut[jt][1], taut[jt0][1]);
            energy_sum = Ft[jt][1] * gamma_w_vec[jt_n][jt - jt_n][4];

            // evolve the collided part of gamma
            for (int i = 0; i < gamma_w_c.size(); i++) {
                gamma0 = gamma_w_c[i];
                gamma_w_c[i] -= dR * (2. * gamma0 / Rttn) * (1 + 3 * Rttn * at[jt][1] * Ht[jt][1] / 2);
                energy_sum -= gamma_w_c[i] * dF[jt_n + i + 1][1]; // jt_n + i + 1 = jtc
            }
            gamma_w_c.push_back(gamma_w_vec[jt_n][jt - jt_n][4]);
            energy_sum -= gamma_w_c[gamma_w_c.size() - 1] * dF[jt][1]; // add part at t_c = t

            if (Rtt0 <= 2 / k) { // Theta function
                boundary_vec[jt][1] += 3 * PI * pow(k, 3.0) / (4 * pow(at[jt][1], 2.0)) * Rtt0 * (2 / k - Rtt0) * dt * (Rttn + Rt0tn) * (2 / k + Rttn + Rt0tn) * Gamma(taut[jt_n][0]) * pow(at[jt_n][1], 3.0) * energy_sum;
            }
            else {
                cout << "2/k - R(t, t_*) < 0\n";
            }
        }
    }



    return boundary_vec;
}

vector<vector<double>> boundary_term_fun2(const double k, int J, vector<vector<double>>& Nk, vector<vector<vector<double>>>& pd, vector<vector<double>>& at, vector<vector<double>>& Ht, vector<vector<double>>& taut, vector<vector<double>>& Ft, vector<vector<vector<double>>>& gamma_w_vec)
{
    vector<vector<double>> boundary_vec(taut.size());
    vector<double> tmp(2);
    double tau, tauj;
    int itn;
    double Rttn; // Rj =  R(t; t_n)
    double gamma0, dR;
    double dt = taut[1][0] - taut[0][0];
    double energy_sum, temp_sum;
    vector<vector<double> > dF;
    vector<double> gamma_w_c;
    const double dd = 10 * (pd[0][1][0] - pd[0][0][0]);

    tmp[0] = Ft[0][0];
    tmp[1] = 0;
    dF.push_back(tmp);

    for (int jtF = 1; jtF < Ft.size(); jtF++) {
        tmp[0] = Ft[jtF][0];
        tmp[1] = Ft[jtF][1] - Ft[jtF - 1][1];
        dF.push_back(tmp);
    }

    // boundary term = 0 before J = j_c is reached
    tmp[1] = 0;
    int jt = 0;
    while (Nk[jt][1] < J) {
        tmp[0] = Nk[jt][0];
        boundary_vec[jt] = tmp;
        jt++;
    }
    int jt0 = jt; // point to the start of approximation

    for (int jt = jt0; jt < Ft.size(); jt++) { // compute B at different times
        tau = taut[jt][1];
        tmp[0] = taut[jt][0];
        energy_sum = 0;
        for (int jtt = jt0; jtt <= jt; jtt++) { // integrate over tilde t
            tauj = taut[jtt][1]; 
            temp_sum = 0;
            for (int jd = 0; jd < pd[0].size(); jd += 10) { // integrate over d
                Rttn = max(0., pd[jtt][jd][0] - 1 / k) + tau - tauj;
                if (pd[jtt][jd][0] < Rttn + 1 / k && pd[jtt][jd][0] > abs(1 / k - Rttn)) { // theta functions
                    if (pd[jtt][jd][0] <= 1 / k) {
                        itn = jtt; // iterator of nucleation time of bubble appearing at jtt, nucleated at d
                    }
                    else {
                        itn = finditerG(tauj - pd[jtt][jd][0] + 1 / k, taut); // nucleation time of bubble appearing at jtt, nucleated at d
                    }
                    temp_sum += dd * pd[jtt][jd][2] * (1 / pd[jtt][jd][0]) * (pd[jtt][jd][0] * pd[jtt][jd][0] - Rttn * Rttn - pow(k, -2.0)) * Ft[jt][1] * gamma_w_vec[itn][jt - itn][4];
                    gamma_w_c = gamma_c_at_t(jt, itn, taut, at, Ht, gamma_w_vec); // gamma(t, t_n, t_n < t_c <= t)[jtc] (jtc = 0 for t_c = t_n + dt or jt_n + 1)
                    for (int jtc = itn+1; jtc < jt; jtc+=10) { // collided part
                        temp_sum -= 10 * dd * pd[jtt][jd][2] * (1 / pd[jtt][jd][0]) * (pd[jtt][jd][0] * pd[jtt][jd][0] - Rttn * Rttn - pow(k, -2.0)) * dF[jtc][1] * gamma_w_c[jtc - itn - 1]; // integrand
                    }
                }
            }
            energy_sum += dt * Nk[jtt][2] * temp_sum; // * dN/dt
            // cout << taut[jt][0] << "   " << taut[jtt][0] << "   " << energy_sum << "   " << Ft[jt][1] << endl;
        }
        energy_sum *= pow(at[jt][1], -2.0) * (3. / 4.) * pow(k, 3.0);
        cout << taut[jt][0] << "   " << energy_sum << endl;
        tmp[1] = energy_sum;
        boundary_vec[jt] = tmp;
    }
    return boundary_vec;
}

// CDF of nucleation distances
vector<vector<vector<double> > > ddist(function<double(double)> Gamma, const double k, int jdmax, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &ttau) {

    const double dt = at[1][0] - at[0][0];
    vector<vector<double> > Nk = Nbark(Gamma, k, Ft, taut, at);

    double dmax = 1.0 / k + radius(taut[taut.size() - 1][1], taut[0][1]);
    double dd = dmax / (1.0 * (jdmax - 1));

    double tau_n, t_n, a_n, d, p, ptot;
    vector<vector<double> > tmp(jdmax, vector<double>(3));
    vector<vector<vector<double> > > Cd(taut.size(), vector<vector<double> >(jdmax, vector<double>(3)));
    for (int jt = 0; jt < taut.size(); jt++) {

        ptot = 0.0;
        for (int jd = 0; jd < jdmax; jd++) {
            d = dd * jd;

            p = 0.;
            // distribution for bubbles nucleating at jt
            if (1.0 / k  > d) {
                p = 4.0 * PI * Gamma(taut[jt][0]) * pow(at[jt][1], 3.0) * pow(d, 2.0) / Nk[jt][2];
            }
            // distribution for bubbles entering our volume at jt from outside
            else if ((1.0 / k) + radius(taut[jt][1], taut[0][1]) > d) {
                tau_n = taut[jt][1] + 1. / k - d; // (conformal) nucleation time of a bubble nucleating at distance d and reaching V(k) at jt
                t_n = findrootG(tau_n,taut[1][0] - taut[0][0], taut); // cosmic nucleation time
                a_n = interpolate(t_n, at); 
                p = 4.0 * PI * Gamma(t_n) * (pow(a_n, 4.0) / at[jt][1]) * pow(d, 2.0) / Nk[jt][2]; // be careful with integrating \delta(R(t;t_n))
            }
            ptot += p;
            tmp[jd][0] = d;
            tmp[jd][1] = dd * ptot;
            tmp[jd][2] = p;
        }
        Cd[jt] = tmp;
    }
    return Cd;
}

// CDF of nucleation distances - old
vector<vector<vector<double> > > ddist2(function<double(double)> Gamma, const double k, int jdmax, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at) {
    
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

double dArea(double rj, double dj, double k)
{
    double R = rj + max(0., dj - 1 / k);
    if (dj >= abs(R - 1 / k)) {
        return (3. / 4.) * (pow(k, 3.0) / dj) * (dj * dj - R * R - pow(k, -2.0));
    }
    else
        return 0;
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
