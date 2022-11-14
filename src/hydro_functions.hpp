#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

void mat_initialize(std::vector<double> &v0, std::vector<double> &rho0, std::vector<double> &Pm, std::vector<double> &e0,
                    std::vector<double> &v, std::vector<double> &e, std::vector<double> &rho, std::vector<double> &as, std::vector<double> &Pr,
                    std::vector<double> &P, std::vector<double> &T0_m, double dx, double gamma, double xf, int p){

    //Number of cells
    int n_cells = v0.size();
    double ie;

    for(int i=0; i<n_cells; i++){

        if(p==0){ // Mach 1.2 Problem
            if((i+0.5)*dx < xf/2){
                rho0[i] = 1;
                rho[i] = rho0[i];
                e0[i] = 2.60510396E+14/rho0[i];
                e[i] = e0[i];
                v0[i] = 1.52172533E+7;
                v[i] = v0[i];
            }else{
                rho0[i] = 1.29731782;
                rho[i] = rho0[i];
                e0[i] = 3.13573034E+14/rho0[i];
                e[i] = e0[i];
                v0[i] = 1.17297805E+7;
                v[i] = v0[i];
            }
        }
        else if (p==1){ //Mach 3 Problem
            if((i+0.5)*dx < xf/2){
                rho0[i] = 1;
                rho[i] = rho0[i];
                e0[i] = 8.68367987E+14/rho0[i];
                e[i] = e0[i];
                v0[i] = 3.80431331E+7;
                v[i] = v0[i];
            }else{
                rho0[i] = 3.00185103;
                rho[i] = rho0[i];
                e0[i] = 1.83229115E+15/rho0[i];
                e[i] = e0[i];
                v0[i] = 1.26732249E+7;
                v[i] = v0[i];
            }
        }

        ie = e0[i] - (0.5*pow(v0[i],2));
        Pm[i] = (gamma - 1) * ie * rho0[i];
        P[i] = Pm[i] + Pr[i];
        as[i] = sqrt((gamma * Pm[i]) / rho0[i]);
        

    }
}

void flux(const std::vector<double> &v0, const std::vector<double> &rho0, const std::vector<double> &P, const std::vector<double> &e0, 
          const std::vector<double> &Er, std::vector<double> &as, std::vector<double> &lu_rho, std::vector<double> &lu_v, 
          std::vector<double> &lu_e, std::vector<double> &lu_Er, std::vector<double> &grad_u, double dt, double gamma) {
    
    int n_cells = v0.size();
    int n_nodes = n_cells + 1;

    // left and right of node Flux values
    std::vector<double> F_rho_l(n_nodes,0);
    std::vector<double> F_rho_r(n_nodes,0);
    std::vector<double> F_v_l(n_nodes,0);
    std::vector<double> F_v_r(n_nodes,0);
    std::vector<double> F_E_l(n_nodes,0);
    std::vector<double> F_E_r(n_nodes,0);

    //Left and right conserved values
    std::vector<double> rho_l(n_nodes,0);
    std::vector<double> rho_r(n_nodes,0);
    std::vector<double> v_l(n_nodes,0);
    std::vector<double> v_r(n_nodes,0);
    std::vector<double> E_l(n_nodes,0);
    std::vector<double> E_r(n_nodes,0);
    std::vector<double> Er_l(n_nodes,0);
    std::vector<double> Er_r(n_nodes,0);
    std::vector<double> P_l(n_nodes,0);
    std::vector<double> P_r(n_nodes,0);

    //Rusanov flux values
    std::vector<double> Frus_rho_l(n_cells,0);
    std::vector<double> Frus_rho_r(n_cells,0);
    std::vector<double> Frus_v_l(n_cells,0);
    std::vector<double> Frus_v_r(n_cells,0);
    std::vector<double> Frus_e_l(n_cells,0);
    std::vector<double> Frus_e_r(n_cells,0);

    //Flux limiter values
    std::vector<double> r_rho(n_nodes,0);
    std::vector<double> r_v(n_nodes,0);
    std::vector<double> r_e(n_nodes,0);
    std::vector<double> r_Er(n_nodes,0);
    std::vector<double> r_P(n_nodes,0);
    std::vector<double> Fl_rho(n_nodes,0);
    std::vector<double> Fl_v(n_nodes,0);
    std::vector<double> Fl_e(n_nodes,0);
    std::vector<double> Fl_Er(n_nodes,0);
    std::vector<double> Fl_P(n_nodes,0);

    //left and right Rad Flux values
    std::vector<double> F_Er_l(n_nodes,0);
    std::vector<double> F_Er_r(n_nodes,0);
    std::vector<double> Frus_Er_l(n_cells,0);
    std::vector<double> Frus_Er_r(n_cells,0);

    //Grad u term
    std::vector<double> F_u_l(n_nodes,0);
    std::vector<double> F_u_r(n_nodes,0);
    std::vector<double> Frus_u_l(n_cells,0);
    std::vector<double> Frus_u_r(n_cells,0);

    //lambda values
    double lpr, lpl;
    //values of alpha
    std::vector<double> Sp(n_nodes,0);

    //Riemann Solver
    for(int i=0; i<n_nodes; i++){

        //Values of lambda for Riemann Solver
        if(i==0){
            lpr = abs(v0[i]) + as[i];
            lpl = lpr;
        }else if(i==n_nodes-1){
            lpl = abs(v0[i-1]) + as[i-1];
            lpr = lpl;
        }else{
            lpl = abs(v0[i-1]) + as[i-1];
            lpr = abs(v0[i]) + as[i];
        }

        Sp[i] = std::max(lpr,lpl);

    }

    //Find flux limiter values
    for(int i=0; i<n_nodes; i++){
        if(i==0){
            r_rho[i] = 0;
            r_v[i] = 0;
            r_e[i] = 0;
            r_Er[i] = 0;
            r_P[i] = 0;
        }else if(i==n_nodes-1){
            r_rho[i] = 0;
            r_v[i] = 0;
            r_e[i] = 0;
            r_Er[i] = 0;
            r_P[i] = 0;
        }else{
            if(rho0[i+1] - rho0[i] == 0){
                r_rho[i] = 0;
            }else{
                r_rho[i] = (rho0[i] - rho0[i-1]) / (rho0[i+1] - rho0[i]);
            }
            if(v0[i+1] - v0[i]== 0){
                r_v[i] = 0;
            }else{
                r_v[i] = (v0[i]- v0[i-1]) / (v0[i+1] - v0[i]);
            }
           if(e0[i+1] - e0[i] == 0){
                r_e[i] = 0;
            }else{
                r_e[i] = (e0[i] - e0[i-1])/ (e0[i+1] - e0[i]);
            }
            if(Er[i+1] - Er[i] == 0){
                r_Er[i] = 0;
            }else{
                r_Er[i] = (Er[i] - Er[i-1])/ (Er[i+1] - Er[i]);
            }if(P[i+1] - P[i] == 0){
                r_P[i] = 0;
            }else{
                r_P[i] = (P[i] - P[i-1])/ (P[i+1] - P[i]);
            }
        }
        
        Fl_rho[i] = (r_rho[i] + abs(r_rho[i])) / (1 + abs(r_rho[i]));
        Fl_v[i] = (r_v[i] + abs(r_v[i])) / (1 + abs(r_v[i]));
        Fl_e[i] = (r_e[i] + abs(r_e[i])) / (1 + abs(r_e[i]));
        Fl_Er[i] = (r_Er[i] + abs(r_Er[i])) / (1 + abs(r_Er[i]));
        Fl_P[i] = (r_P[i] + abs(r_P[i])) / (1 + abs(r_P[i]));
        
    }

    //2nd order reconstruction to get L/R values @ i-1/2 for conserved values

    for(int i=0; i<n_nodes; i++){
        if(i==0){
            rho_r[i] = rho0[i] - ((0.5 * Fl_rho[i]) * (rho0[i+1] - rho0[i]));
            rho_l[i] = rho_r[i];

            v_r[i] = v0[i] - 0.5 * Fl_v[i] * (v0[i+1] - v0[i]);
            v_l[i] = v_r[i];

            E_r[i] = e0[i] - 0.5 * Fl_e[i] * (e0[i+1] - e0[i]);
            E_l[i] = E_r[i];

            Er_r[i] = Er[i] - 0.5 * Fl_Er[i] * (Er[i+1] - Er[i]);
            Er_l[i] = Er_r[i];
            
            P_r[i] = P[i] - 0.5 * Fl_P[i] * (P[i+1] - P[i]);
            P_l[i] = P_r[i];

        }else if(i==n_nodes-1){
            rho_l[i] = rho0[i-1] + ((0.5 * Fl_rho[i-1]) * (rho0[i] - rho0[i-1]));
            rho_r[i] = rho_l[i];

            v_l[i] = v0[i-1] + 0.5 * Fl_v[i-1] * (v0[i] - v0[i-1]);
            v_r[i] = v_l[i];

            E_l[i] = e0[i-1] + 0.5 * Fl_e[i-1] * (e0[i] - e0[i-1]);
            E_r[i] = E_l[i];

            Er_l[i] = Er[i-1] + 0.5 * Fl_Er[i-1] * (Er[i] - Er[i-1]);
            Er_r[i] = Er_l[i];

            P_l[i] = P[i-1] + 0.5 * Fl_P[i-1] * (P[i] - P[i-1]);
            P_r[i] = P_l[i];

        }else{
            rho_l[i] = rho0[i-1] + ((0.5 * Fl_rho[i-1]) * (rho0[i] - rho0[i-1]));
            rho_r[i] = rho0[i] - ((0.5 * Fl_rho[i]) * (rho0[i+1] - rho0[i]));

            v_l[i] = v0[i-1] + 0.5 * Fl_v[i-1] * (v0[i] - v0[i-1]);
            v_r[i] = v0[i] - 0.5 * Fl_v[i] * (v0[i+1] - v0[i]);

            E_l[i] = e0[i-1] + 0.5 * Fl_e[i-1] * (e0[i] - e0[i-1]);
            E_r[i] = e0[i] - 0.5 * Fl_e[i] * (e0[i+1] - e0[i]);
            
            Er_l[i] = Er[i-1] + 0.5 * Fl_Er[i-1] * (Er[i] - Er[i-1]);
            Er_r[i] = Er[i] - 0.5 * Fl_Er[i] * (Er[i+1] - Er[i]);

            P_l[i] = P[i-1] + 0.5 * Fl_P[i-1] * (P[i] - P[i-1]);
            P_r[i] = P[i] - 0.5 * Fl_P[i] * (P[i+1] - P[i]);
            
        }
    }

    //Calculate the Left and right Flux values from the 2nd order conserved values
    for(int i=0; i<n_nodes; i++){
        if(i==0){
            F_rho_r[i] = rho_r[i] * (v_r[i]);
            F_rho_l[i] = F_rho_r[i];

            F_v_r[i] = (rho_r[i] * v_r[i] * v_r[i]) + P_r[i];
            F_v_l[i] = F_v_r[i];

            F_E_r[i] = rho_r[i] * E_r[i]*v_r[i] + P_r[i]*v_r[i];
            F_E_l[i] = F_E_r[i];

            F_Er_r[i] = Er_r[i] * v_r[i];
            F_Er_l[i] = F_Er_r[i];

            F_u_r[i] = v_r[i];
            F_u_l[i] = F_u_r[i];

        }else if(i==n_nodes-1){
            F_rho_l[i] = rho_l[i] * (v_l[i]);
            F_rho_r[i] = F_rho_l[i];

            F_v_l[i] = (rho_l[i] * v_l[i] * v_l[i]) + P_l[i];
            F_v_r[i] = F_v_l[i];

            F_E_l[i] = rho_l[i] * E_l[i]*v_l[i] + P_l[i]*v_l[i];
            F_E_r[i] = F_E_l[i];

            F_Er_l[i] = Er_l[i] * v_l[i];
            F_Er_r[i] = F_Er_l[i];

            F_u_l[i] = v_l[i];
            F_u_r[i] = F_u_l[i];

        }else{
            F_rho_l[i] = rho_l[i] * (v_l[i]);
            F_rho_r[i] = rho_r[i] * (v_r[i]);

            F_v_l[i] =  (rho_l[i] * v_l[i] * v_l[i]) + P_l[i];
            F_v_r[i] = (rho_r[i] * v_r[i] * v_r[i]) + P_r[i];

            F_E_l[i] = rho_l[i] * E_l[i]*v_l[i] + P_l[i]*v_l[i];
            F_E_r[i] = rho_r[i] * E_r[i]*v_r[i] + P_r[i]*v_r[i];

            F_Er_l[i] = Er_l[i] * v_l[i];
            F_Er_r[i] = Er_r[i] * v_r[i];

            F_u_l[i] = v_l[i];
            F_u_r[i] = v_r[i];
        }    
    }

    for(int i=0; i<n_cells; i++){
        if(i==0){
            Frus_rho_r[i] = 0.5 * (F_rho_l[i+1] + F_rho_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1] - rho_l[i+1]));
            Frus_rho_l[i] = Frus_rho_r[i];

            Frus_v_r[i] = 0.5 * (F_v_l[i+1] + F_v_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1]*v_r[i+1] - rho_l[i+1]*v_l[i+1]));
            Frus_v_l[i] = Frus_v_r[i];

            Frus_e_r[i] =  0.5 * (F_E_l[i+1] + F_E_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1]*E_r[i+1] - rho_l[i+1]*E_l[i+1]));
            Frus_e_l[i] = Frus_e_r[i];

            Frus_Er_r[i] =  0.5 * (F_Er_l[i+1] + F_Er_r[i+1]) - (0.5 * Sp[i+1] * (Er_r[i+1] - Er_l[i+1]));
            Frus_Er_l[i] = Frus_Er_r[i];

            Frus_u_r[i] =  0.5 * (F_u_l[i+1] + F_u_r[i+1]);
            Frus_u_l[i] = Frus_u_r[i];

        }else if(i==n_cells-1){
            Frus_rho_l[i] = 0.5 * (F_rho_l[i] + F_rho_r[i]) - (0.5 * Sp[i] * (rho_r[i] - rho_l[i]));
            Frus_rho_r[i] = Frus_rho_l[i];

            Frus_v_l[i] = 0.5 * (F_v_l[i] + F_v_r[i]) - (0.5 * Sp[i] * (rho_r[i]*v_r[i] - rho_l[i]*v_l[i]));
            Frus_v_r[i] = Frus_v_l[i];

            Frus_e_l[i] = 0.5 * (F_E_l[i] + F_E_r[i]) - (0.5 * Sp[i] * (rho_r[i]*E_r[i] - rho_l[i]*E_l[i]));
            Frus_e_r[i] = Frus_e_l[i];

            Frus_Er_l[i] = 0.5 * (F_Er_l[i] + F_Er_r[i]) - (0.5 * Sp[i] * (Er_r[i] - Er_l[i]));
            Frus_Er_r[i] = Frus_Er_l[i];

            Frus_u_l[i] = 0.5 * (F_u_l[i] + F_u_r[i]);
            Frus_u_r[i] = Frus_u_l[i];

        }else{
            Frus_rho_l[i] = 0.5 * (F_rho_l[i] + F_rho_r[i]) - (0.5 * Sp[i] * (rho_r[i] - rho_l[i]));
            Frus_rho_r[i] = 0.5 * (F_rho_l[i+1] + F_rho_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1] - rho_l[i+1]));

            Frus_v_l[i] = 0.5 * (F_v_l[i] + F_v_r[i]) - (0.5 * Sp[i] * (rho_r[i]*v_r[i] - rho_l[i]*v_l[i]));
            Frus_v_r[i] = 0.5 * (F_v_l[i+1] + F_v_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1]*v_r[i+1] - rho_l[i+1]*v_l[i+1]));

            Frus_e_l[i] = 0.5 * (F_E_l[i] + F_E_r[i]) - (0.5 * Sp[i] * (rho_r[i]*E_r[i] - rho_l[i]*E_l[i]));
            Frus_e_r[i] = 0.5 * (F_E_l[i+1] + F_E_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1]*E_r[i+1] - rho_l[i+1]*E_l[i+1]));

            Frus_Er_l[i] = 0.5 * (F_Er_l[i] + F_Er_r[i]) - (0.5 * Sp[i] * (Er_r[i] - Er_l[i]));
            Frus_Er_r[i] = 0.5 * (F_Er_l[i+1] + F_Er_r[i+1]) - (0.5 * Sp[i+1] * (Er_r[i+1] - Er_l[i+1]));

            Frus_u_l[i] = 0.5 * (F_u_l[i] + F_u_r[i]);
            Frus_u_r[i] = 0.5 * (F_u_l[i+1] + F_u_r[i+1]);

        }
    }

    //Find the l(U) values
    for(int i=0; i<n_cells; i++){
        
        lu_rho[i] = -(Frus_rho_r[i] - Frus_rho_l[i]);
        lu_v[i] = -(Frus_v_r[i] - Frus_v_l[i]);
        lu_e[i] = -(Frus_e_r[i] - Frus_e_l[i]);

        //Rad flux term
        lu_Er[i] = (Frus_Er_r[i] - Frus_Er_l[i]);

        //grad(u) term
        grad_u[i] = (Frus_u_r[i] - Frus_u_l[i]);
    }

}

void eularian_calcs(const std::vector<double> &v0, const std::vector<double> &rho0, const std::vector<double> &e0, const std::vector<double> &Pr,
                    const std::vector<double> &lu_rho, const std::vector<double> &lu_v, const std::vector<double> &lu_e, const std::vector<double> &grad_u,
                    const std::vector<double> &cv, std::vector<double> &Pm, std::vector<double> &v, std::vector<double> &e, std::vector<double> &rho,   
                    std::vector<double> &P, std::vector<double> &Th, double dt, double dx, double gamma) {
    int n = v0.size();
    for(int j=0; j<n; j++){
        rho[j] = rho0[j] * dx + dt * lu_rho[j];
        rho[j] = rho[j] / dx;
        v[j] = v0[j]*rho0[j]*dx + dt*lu_v[j];
        v[j] = v[j] / (rho[j]*dx);
        e[j] = e0[j]*rho0[j]*dx + dt * (lu_e[j] + (Pr[j] * grad_u[j]));
        e[j] = e[j] / (rho[j]*dx);
        
        //Calculate internal energy
        double ie = e[j] - (0.5*pow(v[j],2));

        //Calculate new pressure
        Pm[j] = (gamma - 1) * ie * rho[j];

        P[j] = Pm[j] + Pr[j];

        //hydro temp
        Th[j] =  ie/cv[j];
    }                    
}

void reassign(const std::vector<double> &Pm, const std::vector<double> &v, const std::vector<double> &e, const std::vector<double> &rho,
              std::vector<double> &v0, std::vector<double> &rho0,  std::vector<double> &e0, std::vector<double> &as,
              double gamma){

    int n_cells = v0.size();

    //Re-assign values
    for(int i=0; i<n_cells; i++){
        v0[i] = v[i];
        rho0[i] = rho[i];
        e0[i] = e[i];

        //sound speed in material
        as[i] = sqrt((gamma * Pm[i]) / rho0[i]);
        
    }
}

