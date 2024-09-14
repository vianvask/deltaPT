#include "functions.h"

int main (int argc, char *argv[]) {
    
    // parameters of the nucleation rate:
    const double beta = atof(argv[1]);
    const double gammapbeta = atof(argv[2]);
    
    const int Nsim = 1000000; // number of realizations
    const int J = 50; // average F for j>J

    clock_t time_req = clock();
    cout << setprecision(4) << fixed;
    
    rgen mt(time(NULL)*(1+beta)); // random number generator
      
    // bubble nucleation rate, units are chosen such that H0 = 1;
    function<double(double)> Gamma = [beta, gammapbeta](double t) {
        return exp(beta*t - pow(gammapbeta*beta*t,2.0)/2.0);
    };
    
    int jdmax = 100;
    
    // find that time range so that barN(tmin,kmax)=10^-7 and barF(tmax)=10^-3
    vector<double> trange = findtrange(Gamma, 0.0000001, 0.001);
    int jtmax = 3200;
    double tmin = trange[0]; double tmax = trange[1];
    double dt = (tmax-tmin)/(1.0*jtmax);
   
    // average evolution
    vector<vector<double> > Ft, taut, at, Ht;
    vector<double> tmp(2);
    tmp = averageevolution(Gamma, tmin, jtmax, dt, Ft, taut, at, Ht);
    
    double kmax = tmp[0];
    double tkmax = tmp[1];
    
    string filename, filename2, filename3;
    ofstream outfileF, outfileD, outfileT;
    
    filename = "tkR_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileF.open(filename.c_str());
    
    outfileF << tkmax << "   " << kmax << "   " << interpolate(tkmax, Ht) << "   " << Rstar(Gamma, Ft, at, tkmax) << endl;
    outfileF.close();
        
    // list of k/kmax values
    vector<double> klist {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0};
      
    // compute the times when the scales k re-enter horizon
    int jkmax = klist.size();
    vector<double> tklist(jkmax);
    filename = "klist_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileF.open(filename.c_str());
    for (int jk = 0; jk < jkmax; jk++) {
        klist[jk] = kmax*klist[jk];
        tklist[jk] = findtk(klist[jk], tkmax, at, Ht);
        outfileF << klist[jk] << "   " << tklist[jk] << "   " << interpolate(tklist[jk],Ht) << endl;
    }
    outfileF.close();
    
    double k;
    vector<vector<vector<double> > > FkW(klist.size());
    vector<vector<vector<double> > > Nk(klist.size());
    vector<vector<vector<vector<double> > > > pdk(klist.size());
    for (int jk = 0; jk < jkmax; jk++) {
        k = klist[jk];
        
        // evolution of expected number of bubbles, Nk[jk][jt][N,dN/dt]
        Nk[jk] = Nbark(Gamma, k, Ft, taut, at);
        
        // distributions of nucleation distances, pdk[jk][jt][jd][d,CDF]
        pdk[jk] = ddist(Gamma, k, jdmax, Ft, taut, at);
        
        // evolution from j>J bubbles FkW[jk][jt][t,F]
        FkW[jk] = Fk(Nk[jk], pdk[jk], k, J+1, taut);
    }
    
    // output t, F and FW for k=0.9*kmax
    filename = "FkW_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileF.open(filename.c_str());
    for (int jt = 0; jt < at.size(); jt+=1) {
        outfileF << Ft[jt][0] << "   " << Ft[jt][1]<< "   " << FkW[12][jt][1] << "   " << Nk[12][jt][1] << "   " << Nk[12][jt][2];
        outfileF  << endl;
    }
    outfileF.close();
        
    // generate Nsim realizations of tj and dj for the first j<J bubbles and compute the evolution of F and rho
    double tau, FkS;
    vector<int> jtj(J);
    vector<double> tauj(J), dj(J);
    vector<vector<double> > F(jtmax, vector<double > (2));
    vector<vector<double> > rho(jtmax, vector<double > (2));

    int nint = 10000;
    vector<double> x(3);
    vector<vector<double> > xlist(nint, vector<double> (3));
    
    filename = "Fk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileF.open(filename.c_str());
    filename2 = "deltak_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD.open(filename2.c_str());
    filename3 = "tj_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileT.open(filename3.c_str());
    
    // initialize delta binning
    int Nbins = 10000, jbin;
    double xbin = 1.0/(1.0*Nsim*2.0/(1.0*Nbins));
    double rhoavg, delta;
    vector<vector<double> > deltabins(Nbins, vector<double> (jkmax+1));
    for (int j = 0; j < Nbins; j++) {
        deltabins[j][0] = 2.0*j/(1.0*Nbins-1.0)-1.0;
        for (int jk = 0; jk < kmax; jk++) {
            deltabins[j][jk+1] = 0.0;
        }
    }
    
    // generate Nsim realizations
    for (int js = 0; js < Nsim; js++) {
        if (js%32 == 0) {
            cout << "\r" << js/(1.0*Nsim) << "    " << flush;
        }
        for (int jk = 0; jk < jkmax; jk++) {
            k = klist[jk];
            
            // generate times t_j
            jtj = jtlist(Nk[jk], J, mt);
            
            // output times t_j for j<J and k=0.9kmax
            if (jk == 12) {
                for (int j = 0; j < J; j++) {
                    outfileT << taut[jtj[j]][0] << "   ";
                }
                outfileT << endl;
            }
            
            // generate distrances d_j
            for (int j = 0; j < J; j++) {
                tauj[j] = taut[jtj[j]][1];
                dj[j] = findrootG(randomreal(0.0,pdk[jk][jtj[j]][jdmax-1][1],mt), 0.001, pdk[jk][jtj[j]]);
            }
            
            // compute the false vacuum fraction
            for (int jt = 0; jt < jtmax; jt++) {
                tau = taut[jt][1];
                FkS = 1.0;
                
                for (int j = 0; j < J; j++) {
                    FkS *= 1.0-Vfrac(tau, tauj[j], dj[j], k);
                }
                
                tmp[0] = taut[jt][0];
                tmp[1] = FkS*FkW[jk][jt][1];
                F[jt] = tmp;
            }
            
            rho = rhoevolution(F, Ht);
            rhoavg = pow(interpolate(tklist[jk], Ht),2.0);
            
            // output F and rho for k=0.9kmax from the first 10 simulations
            if (js < 10 && jk == 12) {
                for (int jt = 0; jt < at.size(); jt++) {
                    outfileF << F[jt][1] << "   ";
                    outfileD << rho[jt][1]/pow(Ht[jt][1],2.0)-1.0 << "   ";
                }
                outfileF << endl;
                outfileD << endl;
            }
            
            // bin the delta distribution
            delta = interpolate(tklist[jk], rho)/rhoavg - 1.0;
            jbin = max(0, min(Nbins-1, (int) round(Nbins*(delta+1.0)/2.0)));
            deltabins[jbin][jk+1] += xbin;
        }
    }
    
    cout << "\r" << "1.0000    " << endl;
    outfileF.close();
    outfileD.close();
    outfileT.close();
    
    // output delta distribution
    filename2 = "deltabinsk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD.open(filename2.c_str());
    for (int j = 0; j < Nbins; j++) {
        outfileD << deltabins[j][0] << "    ";
        for (int jk = 0; jk < jkmax; jk++) {
            outfileD << deltabins[j][jk+1] << "    ";
        }
        outfileD << endl;
    }
    outfileD.close();
    
    time_req = clock() - time_req;
    cout << "total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " minutes." << endl;
    
    return 0;
}
