#include "functions.h"

int main (int argc, char *argv[]) {
    
    // parameters of the nucleation rate:
    const double beta = atof(argv[1]);
    const double gammapbeta = atof(argv[2]);
    
    const int Nsim = 10000; // number of realizations
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
    vector<double> trange = findtrange(Gamma, 0.0000001, 0.00001);
    int jtmax = 3200;
    double tmin = trange[0]; double tmax = trange[1];
    double dt = (tmax-tmin)/(1.0*jtmax);
   
    // average evolution
    vector<vector<double> > Ft, taut, at, Ht, rhoRt, rhoVt;
    vector<double> tmp(2);
    tmp = averageevolution(Gamma, tmin, jtmax, dt, Ft, taut, at, Ht, rhoRt, rhoVt);
        
    double kmax = tmp[0];
    double tkmax = tmp[1];
    
    string filename, filename2, filename3, filename4;
    ofstream outfileF, outfileD1, outfileD2, outfileD3;
    
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
    
    filename = "Fk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileF.open(filename.c_str());
    filename2 = "deltak_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD1.open(filename2.c_str());
    filename3 = "phik_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD2.open(filename3.c_str());
    filename4 = "phiBk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD3.open(filename4.c_str());
    
    // initialize delta binning
    int Nbins = 10000, jbin;
    double xbin = 1.0/(1.0*Nsim*2.0/(1.0*Nbins));
    vector<vector<double> > deltabins(Nbins, vector<double> (jkmax+1));
    for (int j = 0; j < Nbins; j++) {
        deltabins[j][0] = 2.0*j/(1.0*Nbins-1.0)-1.0;
        for (int jk = 0; jk < kmax; jk++) {
            deltabins[j][jk+1] = 0.0;
        }
    }
    vector<vector<double> > phibins = deltabins;
    vector<vector<double> > phiBbins = deltabins;
    
    // generate Nsim realizations
    double x0;
    vector<vector<double> > F, delta, phi, phiB;
    for (int js = 0; js < Nsim; js++) {
        if (js%32 == 0) {
            cout << "\r" << js/(1.0*Nsim) << "    " << flush;
        }
        for (int jk = 0; jk < jkmax; jk++) {
            F.clear(); delta.clear(); phi.clear(); phiB.clear();
            rhoevolutionCG(F, delta, phi, phiB, taut, at, Ht, rhoRt, rhoVt, FkW[jk], Nk[jk], pdk[jk], klist[jk], J, jdmax, mt);
            
            // output F and rho for k=0.9kmax from the first 10 simulations
            if (js < 10 && jk == 12) {
                for (int jt = 0; jt < at.size(); jt++) {
                    outfileF << F[jt][1] << "   ";
                    outfileD1 << delta[jt][1] << "   ";
                    outfileD2 << phi[jt][1] << "   ";
                    outfileD3 << phiB[jt][1] << "   ";
                }
                outfileF << endl;
                outfileD1 << endl;
                outfileD2 << endl;
                outfileD3 << endl;
            }
            
            // bin the delta and zeta distributions
            x0 = interpolate(tklist[jk], delta);
            jbin = max(0, min(Nbins-1, (int) round(Nbins*(x0+1.0)/2.0)));
            deltabins[jbin][jk+1] += xbin;
            
            x0 = interpolate(tklist[jk], phi);
            jbin = max(0, min(Nbins-1, (int) round(Nbins*(x0+1.0)/2.0)));
            phibins[jbin][jk+1] += xbin;
            
            x0 = interpolate(tklist[jk], phiB);
            jbin = max(0, min(Nbins-1, (int) round(Nbins*(x0+1.0)/2.0)));
            phiBbins[jbin][jk+1] += xbin;
        }
    }
    
    cout << "\r" << "1.0000    " << endl;
    outfileF.close();
    outfileD1.close();
    outfileD2.close();
    outfileD3.close();
    
    // output delta distribution
    filename = "deltabinsk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD1.open(filename.c_str());
    filename2 = "phibinsk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD2.open(filename2.c_str());
    filename3 = "phiBbinsk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD3.open(filename3.c_str());
    for (int j = 0; j < Nbins; j++) {
        outfileD1 << deltabins[j][0] << "    ";
        outfileD2 << phibins[j][0] << "    ";
        outfileD3 << phiBbins[j][0] << "    ";
        for (int jk = 0; jk < jkmax; jk++) {
            outfileD1 << deltabins[j][jk+1] << "    ";
            outfileD2 << phibins[j][jk+1] << "    ";
            outfileD3 << phiBbins[j][jk+1] << "    ";
        }
        outfileD1 << endl;
        outfileD2 << endl;
        outfileD3 << endl;
    }
    outfileD1.close();
    outfileD2.close();
    outfileD3.close();

    time_req = clock() - time_req;
    cout << "total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " minutes." << endl;
    
    return 0;
}
