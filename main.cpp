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
    
    int jdmax = 2000; // precision in binning of nucleation distance generating function
    
    // find that time range so that barN(tmin,kmax)=10^-7 and barF(tmax)=10^-3
    vector<double> trange = findtrange(Gamma, 0.0000001, 0.00001);
    int jtmax = 3200; // number of timesteps
    double tmin = trange[0]; double tmax = trange[1];
    double dt = (tmax-tmin)/(1.0*jtmax);
   
    // average evolution
    vector<vector<double> > Ft, taut, at, Ht, rhoRt, rhoVt, ttau;
    vector<double> tmp(2);
    tmp = averageevolution(Gamma, tmin, jtmax, dt, Ft, taut, ttau, at, Ht, rhoRt, rhoVt);
        
    double kmax = tmp[0]; // maximum value of k that reaches the Hubble horizon
    double tkmax = tmp[1]; // time at which kmax reaches the horizon
    
    string filename, filename2, filename3, filename4;
    ofstream outfileF, outfileD1, outfileD2, outfileD3, outfileD4;

    // print evolution of the average patch
    filename = "bkg_beta_" + to_string_prec(beta, 3) + "_gammaperbeta_" + to_string_prec(gammapbeta, 3) + ".dat";
    outfileF.open(filename.c_str());

    for (int i = 0; i < at.size(); i++) {
        outfileF << at[i][0] << "   " << at[i][1] << "   " << Ht[i][1] << "   " << taut[i][1] << "   " << rhoVt[i][1] << "   " << rhoRt[i][1] << endl;
    }
    outfileF.close();
    
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
    vector<vector<vector<double> > > FkW2(klist.size()); // temp copy for comaprison between two computation methods
    vector<vector<vector<double> > > Nk(klist.size());
    vector<vector<vector<double> > > gamma_w; // Lorentz factor of uncollided bubble walls
    vector<vector<vector<double> > > boundary_term(klist.size()); // evaluate the energy flux from j>j_c bubbles from their time and distance distributions
    // vector<vector<vector<double> > > boundary_term2(klist.size()); // cross-check with different computation method

    // R_n_vs_R_H to be removed
    vector< vector<double> > R_n_vs_R_H; // check when bubble nucleated at time t_n grows to the horizon size and what is the value of gamma at that time [t_n, t, tau_n, t, gamma*sigma]
    // evolution of wall gamma factor for uncollided bubble gamma_w[jt][jt_n][t,t_n,tau,tau_n,gamma*sigma]
    cout << "Starting computation of gamma\n";
    gamma_w = gamma_w_fun(taut, at, Ht, R_n_vs_R_H, rhoVt[0][1]); // gamma_w[jt_n][jt][t,t_n,tau,tau_n,gamma*sigma] 

    // initialize delta binning
    int Nbins = 10000, jbin;
    double xbin = 1.0 / (1.0 * Nsim * 2.0 / (1.0 * Nbins));
    vector<vector<double> > deltaNbins(Nbins, vector<double>(jkmax + 1));
    for (int j = 0; j < Nbins; j++) {
        deltaNbins[j][0] = 2.0 * j / (1.0 * Nbins - 1.0) - 1.0;
        for (int jk = 0; jk < kmax; jk++) {
            deltaNbins[j][jk + 1] = 0.0;
        }
    }
    vector<vector<double> > phibins = deltaNbins;
    vector<vector<double> > deltaCbins = deltaNbins;
    vector<vector<double> > Rbins = deltaNbins;

    // this part can be computed separately to generate intermediate results independent of simulations - switched off as it messes up with RAM
    /*
    
    vector<vector<vector<vector<double> > > > pdk(klist.size());

    for (int jk = 0; jk < jkmax; jk++) {
        k = klist[jk];
        cout << "jk = " << jk << ", k = " << k << endl;

        cout << "Compute Nk\n";
        // evolution of expected number of bubbles, Nk[jk][jt][t,N,dN/dt]
        Nk[jk] = Nbark(Gamma, k, Ft, taut, at);
        if (interpolate(tklist[jk], Nk[jk]) < J) {
            cout << "WARNING: Nbark at horizon reentry is smaller than J\n";
        }

        cout << "Compute pdk\n";
        // distributions of nucleation distances, pdk[jk][jt][jd][d,CDF]
        pdk[jk] = ddist(Gamma, k, jdmax, Ft, taut, at, ttau);


        cout << "Compute FkW\n";
        // evolution from j>J bubbles FkW[jk][jt][t,F]
        FkW[jk] = Fk(Gamma, Nk[jk], k, J + 1, taut, at);

        cout << "Compute FkW2\n";
        // evolution from j>J bubbles FkW[jk][jt][t,F]
        FkW2[jk] = Fk2(Nk[jk], pdk[jk], k, J + 1, taut);

        cout << "boundary_term\n";
        // boundary_term[jk][jt][t, \mathcal B(t,k)]
        boundary_term[jk] = boundary_term_fun(Gamma, k, J+1, Nk[jk], at, Ht, taut, Ft, gamma_w);

        cout << "computed boundary_term\n";

        // numerical check, switched off
        
        // cout << "boundary_term2\n";
        // boundary_term2[jk] = boundary_term_fun2(k, J + 1, Nk[jk], pdk[jk], at, Ht, taut, Ft, gamma_w_tmp);
        

    }
    */
    
    //vector<vector<vector<Row> > > gamma_w_c; // Lorentz factor of collided bubble walls - outdated idea, keep for making plots, then remove

    //gamma_w_c = gamma_w_c_fun(taut, at, Ht, gamma_w);

    
    // bunch of outputs below
    /*
    //output the boundary term
    filename = "boundary_approx_" + to_string_prec(beta, 3) + "_gammaperbeta_" + to_string_prec(gammapbeta, 3) + ".dat";
    outfileF.open(filename.c_str());
    for (int jt = 0; jt < boundary_term[0].size(); jt++) {
        outfileF << boundary_term[0][jt][0]; // time
        // for(int jk = 12; jk < 13; jk++){
        for (int jk = 0; jk < jkmax; jk++) {
            outfileF << "   " << boundary_term[jk][jt][1]; // value of boundary_term for given k
        }
        outfileF << endl;
    }
    outfileF.close();
    
    
    //output the boundary term
    filename = "boundary_check_" + to_string_prec(beta, 3) + "_gammaperbeta_" + to_string_prec(gammapbeta, 3) + ".dat";
    outfileF.open(filename.c_str());
    for (int jt = 0; jt < boundary_term[12].size(); jt++) {
        outfileF << boundary_term[12][jt][0] << "   " << boundary_term[12][jt][1] << "   " << boundary_term2[12][jt][1] << endl;
    }
    outfileF.close();
    */
    /*
    //output pdk for k = 0.9kmax
    for (int i = 0; i < pdk[12].size(); i++) {
        if (i % 320 == 1) {
            filename = "pdk_" + to_string_prec(beta, 3) + "_gammaperbeta_" + to_string_prec(gammapbeta, 3) + "_t_" + to_string_prec(taut[i][0], 3) + ".dat";
            outfileF.open(filename.c_str());
            for (int jd = 0; jd < pdk[12][i].size(); jd++) {
                outfileF << pdk[12][i][jd][0] << "   " << pdk[12][i][jd][1] << "   " << pdk[12][i][jd][2] << endl;
            }
            outfileF.close();
        }
    }

    cout << "1.0 / k = " << 1.0 / klist[12] << endl;

    //output R_n vs R_H
    filename = "R_n_vs_R_H_" + to_string_prec(beta, 3) + "_gammaperbeta_" + to_string_prec(gammapbeta, 3) + ".dat";
    outfileF.open(filename.c_str());
    for (int jt = 0; jt < R_n_vs_R_H.size(); jt++) {
        outfileF << R_n_vs_R_H[jt][0] << "   " << R_n_vs_R_H[jt][1] << "   " << R_n_vs_R_H[jt][2] << "   " << R_n_vs_R_H[jt][3] << "   " << R_n_vs_R_H[jt][4] << endl;

    }
    outfileF.close();
    
    // output the evolution of Lorentz factor
    for (int i = 0; i < gamma_w.size(); i++) {
        if (i % 320 == 0) {
            filename = "gamma_w2_" + to_string_prec(beta, 3) + "_gammaperbeta_" + to_string_prec(gammapbeta, 3) + "_t_" + to_string_prec(gamma_w[i][0][0], 3) + ".dat";
            outfileF.open(filename.c_str());
            for (int jt = 0; jt < gamma_w[i].size(); jt++) {
                outfileF << gamma_w[i][jt][0] << "   " << gamma_w[i][jt][1] << "   " << gamma_w[i][jt][2] << "   " << gamma_w[i][jt][3] << "   " << gamma_w[i][jt][4] << endl;
            }
            outfileF.close();
        }
    }

    // we want to print out the evolution of Lorentz factor after collision for bubble nucleated at t = -0.89 (jt = 1280)
    
    int itn = 1280; 
    vector<vector<Row> > gamma_w_c_print;
    vector<Row> tmp2d;

    gamma_w_c_print = gamma_w_c[0];

    cout << endl << gamma_w_c_print.size() << endl;

    for (int i = 0; i < gamma_w_c_print.size(); i++) {
        if (i % ((gamma_w_c_print.size()+1)/10) == 1) {
            cout << i << endl;
            filename = "gamma_w_c2_" + to_string_prec(beta, 3) + "_gammaperbeta_" + to_string_prec(gammapbeta, 3) + "_tn_" + to_string_prec(gamma_w_c_print[i][0].t_n, 3) + "_tc_" + to_string_prec(gamma_w_c_print[i][0].t_c, 3) + ".dat";
            outfileF.open(filename.c_str());
            for (int jt = 0; jt < gamma_w_c_print[i].size(); jt++) {
                outfileF << gamma_w_c_print[i][jt].t << "   " << gamma_w_c_print[i][jt].t_n << "   " << gamma_w_c_print[i][jt].t_c << "   " << gamma_w_c_print[i][jt].tau << "   " << gamma_w_c_print[i][jt].tau_n << "   " << gamma_w_c_print[i][jt].tau_c << "   " << gamma_w_c_print[i][jt].gamma << endl;
            }
            outfileF.close();
        }
    }
    */

        
    filename = "Fk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileF.open(filename.c_str());
    filename2 = "deltak_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD1.open(filename2.c_str());
    filename3 = "phik_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD2.open(filename3.c_str());
    filename4 = "phiBk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD3.open(filename4.c_str());

    ofstream outfile_c; // combined results from first 10 sims - print out everything we can
    string filename_c;
    
    vector<vector<vector<double> > > pdk; // this object is the most expensive in terms of RAM memory, so we want to keep it only for the value of k we're currently working on

    // generate Nsim realizations
    double x0;
    vector<vector<double> > F, delta, phi, phiB, v, Boundary, R, Lapl_phi;
    //vector<vector<double> > delta_Y, phi_Y, phiB_Y, v_Y, Boundary_Y, R_Y, deltaPnad_Y; // Yann's proposed system
    vector<vector<double> > delta_C; // new, super smart method
    // big loop over different volumes of evolving patch
    for (int jk = 0; jk < jkmax; jk++) {
        k = klist[jk];
        cout << "jk = " << jk << ",k = " << k << endl;

        cout << "Compute Nk\n";
        // evolution of expected number of bubbles, Nk[jk][jt][t,N,dN/dt]
        Nk[jk] = Nbark(Gamma, k, Ft, taut, at);
        if (interpolate(tklist[jk], Nk[jk]) < J) {
            cout << "WARNING: Nbark at horizon reentry is smaller than J\n";
        }

        cout << "Compute pdk\n";
        // distributions of nucleation distances, pdk[jt][jd][d,CDF]
        pdk = ddist(Gamma, k, jdmax, Ft, taut, at, ttau);


        cout << "Compute FkW\n";
        // evolution from j>J bubbles FkW[jk][jt][t,F]
        FkW[jk] = Fk(Gamma, Nk[jk], k, J + 1, taut, at);

        cout << "Compute FkW2\n";
        // evolution from j>J bubbles FkW[jk][jt][t,F]
        FkW2[jk] = Fk2(Nk[jk], pdk, k, J + 1, taut);

        cout << "boundary_term\n";
        // boundary_term[jk][jt][t, \mathcal B(t,k)]
        boundary_term[jk] = boundary_term_fun(Gamma, k, J + 1, Nk[jk], at, Ht, taut, Ft, gamma_w);

        cout << "Starting simulations\n";

        // now simulations go
        for (int js = 0; js < Nsim; js++) {
            if (js % 32 == 0) {
                cout << "\r" << js / (1.0 * Nsim) << "    " << flush;
            }
            F.clear(); delta_C.clear(); delta.clear(); phi.clear(); phiB.clear(); v.clear(); Boundary.clear(); R.clear(); Lapl_phi.clear();
            rhoevolutionNG_new(F, delta_C, delta, phi, phiB, v, Boundary, R, Lapl_phi, taut, at, Ht, rhoRt, rhoVt, FkW2[jk], Nk[jk], pdk, klist[jk], J, jdmax, mt, gamma_w, boundary_term[jk]);
            
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

                filename_c = "results_combined_" + to_string_prec(beta, 3) + "_gammaperbeta_" + to_string_prec(gammapbeta, 3) + "_" + to_string(js) + ".dat";
                outfile_c.open(filename_c.c_str());
                for (int jt = 0; jt < at.size(); jt++) {
                    outfile_c << F[jt][0] << "   " << F[jt][1] << "   " << delta_C[jt][1] << "   " << delta[jt][1] << "   " << phi[jt][1] << "   " << v[jt][1] << "   " << Boundary[jt][1] << "   " << R[jt][1] << "   " << Lapl_phi[jt][1] << endl;
                }

                outfile_c.close();
            }

            // bin the delta_N, delta_C, phi and R distributions
            x0 = interpolate(tklist[jk], delta);
            jbin = max(0, min(Nbins-1, (int) round(Nbins*(x0+1.0)/2.0)));
            deltaNbins[jbin][jk+1] += xbin;

            x0 = interpolate(tklist[jk], delta_C);
            jbin = max(0, min(Nbins - 1, (int)round(Nbins * (x0 + 1.0) / 2.0)));
            deltaCbins[jbin][jk + 1] += xbin;
            
            x0 = interpolate(tklist[jk], phi);
            jbin = max(0, min(Nbins-1, (int) round(Nbins*(x0+1.0)/2.0)));
            phibins[jbin][jk+1] += xbin;
            
            x0 = interpolate(tklist[jk], R);
            jbin = max(0, min(Nbins-1, (int) round(Nbins*(x0+1.0)/2.0)));
            Rbins[jbin][jk+1] += xbin;
        }

        cout << "\r" << "1.0000    " << endl;
    }
    
    
    outfileF.close();
    outfileD1.close();
    outfileD2.close();
    outfileD3.close();

    // output t, F and FW for k=0.9*kmax
    filename = "FkW_beta_" + to_string_prec(beta, 3) + "_gammaperbeta_" + to_string_prec(gammapbeta, 3) + ".dat";
    outfileF.open(filename.c_str());
    for (int jt = 0; jt < at.size(); jt += 1) {
        outfileF << Ft[jt][0] << "   " << Ft[jt][1] << "   " << FkW[12][jt][1] << "   " << FkW2[12][jt][1] << "   " << Nk[12][jt][1] << "   " << Nk[12][jt][2];
        outfileF << endl;
    }
    outfileF.close();
    
    // output delta distribution
    filename = "deltaNbinsk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD1.open(filename.c_str());
    filename2 = "phibinsk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD2.open(filename2.c_str());
    filename3 = "deltaCbinsk_beta_" + to_string_prec(beta,3) + "_gammaperbeta_" + to_string_prec(gammapbeta,3) + ".dat";
    outfileD3.open(filename3.c_str());
    filename4 = "Rbinsk_beta_" + to_string_prec(beta, 3) + "_gammaperbeta_" + to_string_prec(gammapbeta, 3) + ".dat";
    outfileD4.open(filename4.c_str());
    for (int j = 0; j < Nbins; j++) {
        outfileD1 << deltaNbins[j][0] << "    ";
        outfileD2 << phibins[j][0] << "    ";
        outfileD3 << deltaCbins[j][0] << "    ";
        outfileD4 << Rbins[j][0] << "    ";
        for (int jk = 0; jk < jkmax; jk++) {
            outfileD1 << deltaNbins[j][jk+1] << "    ";
            outfileD2 << phibins[j][jk+1] << "    ";
            outfileD3 << deltaCbins[j][jk+1] << "    ";
            outfileD4 << Rbins[j][jk+1] << "    ";
        }
        outfileD1 << endl;
        outfileD2 << endl;
        outfileD3 << endl;
        outfileD4 << endl;
    }
    outfileD1.close();
    outfileD2.close();
    outfileD3.close();
    outfileD4.close();

    time_req = clock() - time_req;
    cout << "total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " minutes." << endl;
    
    return 0;
}
