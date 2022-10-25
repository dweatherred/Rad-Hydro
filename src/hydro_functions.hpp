#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

void mat_initialize(std::vector<double> &m,std::vector<double> &v0, std::vector<double> &rho0, std::vector<double> &Pm, std::vector<double> &e0,
                    std::vector<double> &v, std::vector<double> &e, std::vector<double> &rho, std::vector<double> &as, std::vector<double> &Pr,
                    std::vector<double> &P, std::vector<double> &T0_m, std::vector<double> &Th, double dx, double gamma, double xf){

    //Number of cells
    int n_cells = v0.size();

    for(int i=0; i<n_cells; i++){
        
        
        if((i+0.5)*dx < xf/2){
            rho0[i] = 1;
            rho[i] = rho0[i];
            v0[i] = -1.52172533E+7; //1.52172533E+7  ; //3.80431331E+7;
            v[i] = v0[i];
            m[i] = rho0[i] * dx;
            Pm[i] = (rho[i] * T0_m[i])/gamma;
            P[i] = Pm[i] + Pr[i];
            e0[i] = (P[i] / (rho0[i] * (gamma - 1))) + (0.5 * (v0[i]*v0[i]));
            e[i] = e0[i];
        }else{
            rho0[i] = 1.29731782; //3.00185103;
            rho[i] = rho0[i];
            v0[i] = -1.17297805E+7; //1.17297805E+7; //1.26732249E+7;
            v[i] = v0[i];
            m[i] = rho0[i] * dx;
            Pm[i] = (rho[i] * T0_m[i])/gamma;
            P[i] = Pm[i] + Pr[i];
            e0[i] = (P[i] / (rho0[i] * (gamma - 1))) + (0.5 * (v0[i]*v0[i]));
            e[i] = e0[i];
        }
        
       
       /*
       if(i*dx > 0.2 && i*dx < 0.4 ){
            rho0[i] = 3.0;
            rho[i] = rho0[i];
            v0[i] = 0;
            v[i] = v0[i];
            m[i] = rho0[i] * dx;
            Pm[i] = 1.0;
            e0[i] = Pm[i] / (rho0[i] * (gamma - 1));
            e[i] = e0[i];
            Th[i] = T0_m[i];
        }else if(i*dx > 0.6 && i*dx < 0.8 ){
            rho0[i] = 3.0;
            rho[i] = rho0[i];
            v0[i] = 0;
            v[i] = v0[i];
            m[i] = rho0[i] * dx;
            Pm[i] = 1.0;
            e0[i] = Pm[i] / (rho0[i] * (gamma - 1));
            e[i] = e0[i];
            Th[i] = T0_m[i];
        }else{
            rho0[i] = 1.0;
            rho[i] = rho0[i];
            v0[i] = 0;
            v[i] = v0[i];
            m[i] = rho0[i] * dx;
            Pm[i] = 0.1;
            e0[i] = Pm[i] / (rho0[i] * (gamma - 1));
            e[i] = e0[i];
        }
        */
        
        as[i] = sqrt((gamma * Pm[i]) / rho0[i]);

    }
}

void flux_pre(std::vector<double> &m, std::vector<double> &v0, std::vector<double> &rho0, std::vector<double> &P, std::vector<double> &e0,
              std::vector<double> &v, std::vector<double> &e, std::vector<double> &rho, std::vector<double> &as, std::vector<double> &lu_rho,
              std::vector<double> &lu_v, std::vector<double> &lu_e, std::vector<double> &lu_Er, std::vector<double> &grad_u,
              std::vector<double> &Er, double dt, double gamma, double dx){

    int n_cells = m.size();
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

    std::cout << " PREDICTOR " << std::endl;

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

            Frus_u_l[i] =  0.5 * (F_u_l[i] + F_u_r[i]);
            Frus_u_r[i] = Frus_u_l[i];

        }else{
            Frus_rho_l[i] = 0.5 * (F_rho_l[i] + F_rho_r[i]) - (0.5 * Sp[i] * (rho_r[i] - rho_l[i]));
            Frus_rho_r[i] = 0.5 * (F_rho_l[i+1] + F_rho_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1] - rho_l[i+1]));

            Frus_v_l[i] = 0.5 * (F_v_l[i] + F_v_r[i]) - (0.5 * Sp[i] * (rho_r[i]*v_r[i] - rho_l[i]*v_l[i]));
            Frus_v_r[i] = 0.5 * (F_v_l[i+1] + F_v_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1]*v_r[i+1] - rho_l[i+1]*v_l[i+1]));

            Frus_e_l[i] = 0.5 * (F_E_l[i] + F_E_r[i]) - (0.5 * Sp[i] * (rho_r[i]*E_r[i] - rho_l[i]*E_l[i]));
            Frus_e_r[i] =  0.5 * (F_E_l[i+1] + F_E_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1]*E_r[i+1] - rho_l[i+1]*E_l[i+1]));

            Frus_Er_l[i] = 0.5 * (F_Er_l[i] + F_Er_r[i]) - (0.5 * Sp[i] * (Er_r[i] - Er_l[i]));
            Frus_Er_r[i] =  0.5 * (F_Er_l[i+1] + F_Er_r[i+1]) - (0.5 * Sp[i+1] * (Er_r[i+1] - Er_l[i+1]));

            Frus_u_l[i] =  0.5 * (F_u_l[i] + F_u_r[i]);
            Frus_u_r[i] =  0.5 * (F_u_l[i+1] + F_u_r[i+1]);

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

void eularian_calcs_pre(std::vector<double> &v0, std::vector<double> &rho0, std::vector<double> &Pm_pre, std::vector<double> &e0,
                    std::vector<double> &v_pre, std::vector<double> &e_pre, std::vector<double> &rho_pre, std::vector<double> &lu_rho,
                    std::vector<double> &lu_v, std::vector<double> &lu_e, std::vector<double> &Pr, std::vector<double> &grad_u, 
                    double dt, double dx, double gamma){

    int n = v0.size();
    for(int j=0; j<n; j++){
        rho_pre[j] = rho0[j] * dx + dt/2 * lu_rho[j];
        rho_pre[j] = rho_pre[j] / dx;
        v_pre[j] = v0[j]*rho0[j]*dx + dt/2*lu_v[j];
        v_pre[j] = v_pre[j] / (rho_pre[j]*dx);
        e_pre[j] = e0[j]*rho0[j]*dx + dt/2 * (lu_e[j] + Pr[j] * grad_u[j]);
        e_pre[j] = e_pre[j] / (rho_pre[j]*dx);

        //Calculate internal energy
        double ie = e_pre[j] - (0.5*pow(v_pre[j],2));

        //Calculate predicted material Pressure
        Pm_pre[j] =  (gamma - 1) * ie * rho_pre[j];

        std::cout << " j " << j << " Pm_pre " << Pm_pre[j] << " rho_pre " << rho_pre[j] << " ie " << ie << 
        " v_pre " << v_pre[j] << " e0 " << e0[j] << " e_pre " << e_pre[j] << " grad_u " << grad_u[j] << std::endl;

    }

}

//! Corrector Values

void flux_cor(std::vector<double> &m, std::vector<double> &P, std::vector<double> &Pm_pre, std::vector<double> &v_pre, std::vector<double> &e_pre,
              std::vector<double> &rho_pre, std::vector<double> &as, std::vector<double> &lu_rho, std::vector<double> &lu_v,
              std::vector<double> &lu_e, std::vector<double> &lu_Er, std::vector<double> &grad_u, double dt, double gamma, double dx){


    int n_cells = m.size();
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
    std::vector<double> r_P(n_nodes,0);
    std::vector<double> Fl_rho(n_nodes,0);
    std::vector<double> Fl_v(n_nodes,0);
    std::vector<double> Fl_e(n_nodes,0);
    std::vector<double> Fl_P(n_nodes,0);

    //Grad u term
    std::vector<double> F_u_l(n_nodes,0);
    std::vector<double> F_u_r(n_nodes,0);
    std::vector<double> Frus_u_l(n_cells,0);
    std::vector<double> Frus_u_r(n_cells,0);

    //lambda values
    double lpr, lpl;
    //values of alpha
    std::vector<double> Sp(n_nodes,0);

    //sound speed in material
    for(int i=0; i<n_cells; i++){
        
        as[i] = sqrt((gamma * Pm_pre[i]) / rho_pre[i]);

    }

    //Riemann Solver
    for(int i=0; i<n_nodes; i++){

        //Values of lambda for Riemann Solver
        if(i==0){
            lpr = abs(v_pre[i]) + as[i];
            lpl = lpr;
        }else if(i==n_nodes-1){
            lpl = abs(v_pre[i-1]) + as[i-1];
            lpr = lpl;
        }else{
            lpl = abs(v_pre[i-1]) + as[i-1];
            lpr = abs(v_pre[i]) + as[i];
        }
        
        Sp[i] = std::max(lpr,lpl);

    }

    //Find flux limiter values
    for(int i=0; i<n_nodes; i++){
        if(i==0){
            r_rho[i] = 0;
            r_v[i] = 0;
            r_e[i] = 0;
            r_P[i] = 0;
        }else if(i==n_nodes-1){
            r_rho[i] = 0;
            r_v[i] = 0;
            r_e[i] = 0;
            r_P[i] = 0;
        }else{
            if(rho_pre[i+1] - rho_pre[i] == 0){
                r_rho[i] = 0;
            }else{
                r_rho[i] = (rho_pre[i] - rho_pre[i-1]) / (rho_pre[i+1] - rho_pre[i]);
            }
            if(v_pre[i+1] - v_pre[i]== 0){
                r_v[i] = 0;
            }else{
                r_v[i] = (v_pre[i]- v_pre[i-1]) / (v_pre[i+1] - v_pre[i]);
            }
           if(e_pre[i+1] - e_pre[i] == 0){
                r_e[i] = 0;
            }else{
                r_e[i] = (e_pre[i] - e_pre[i-1])/ (e_pre[i+1] - e_pre[i]);
            }if(P[i+1] - P[i] == 0){
                r_P[i] = 0;
            }else{
                r_P[i] = (P[i] - P[i-1])/ (P[i+1] - P[i]);
            }
        }
        
        /*
        Fl_rho[i] = (r_rho[i] + abs(r_rho[i])) / (1 + abs(r_rho[i]));
        Fl_v[i] = (r_v[i] + abs(r_v[i])) / (1 + abs(r_v[i]));
        Fl_e[i] = (r_e[i] + abs(r_e[i])) / (1 + abs(r_e[i]));
        Fl_P[i] = (r_P[i] + abs(r_P[i])) / (1 + abs(r_P[i]));
        */
        
        Fl_rho[i] = (2 * r_rho[i]) / (pow(r_rho[i],2) + 1);
        Fl_v[i] = (2 * r_v[i]) / (pow(r_v[i],2) + 1);
        Fl_e[i] = (2 * r_e[i]) / (pow(r_e[i],2) + 1);
        Fl_P[i] = (2 * r_P[i]) / (pow(r_P[i],2) + 1);
        

    }

    //2nd order reconstruction to get L/R values @ i-1/2 for conserved values
    std::cout << " Corrector " << std::endl;
    for(int i=0; i<n_nodes; i++){
        if(i==0){
            rho_r[i] = rho_pre[i] - ((0.5 * Fl_rho[i]) * (rho_pre[i+1] - rho_pre[i]));
            rho_l[i] = rho_r[i];

            v_r[i] = v_pre[i] - 0.5 * Fl_v[i] * (v_pre[i+1] - v_pre[i]);
            v_l[i] = v_r[i];

            E_r[i] = e_pre[i] - 0.5 * Fl_e[i] * (e_pre[i+1] - e_pre[i]);
            E_l[i] = E_r[i];

            P_r[i] = P[i] - 0.5 * Fl_P[i] * (P[i+1] - P[i]);
            P_l[i] = P_r[i];

        }else if(i==n_nodes-1){
            rho_l[i] = rho_pre[i-1] + ((0.5 * Fl_rho[i-1]) * (rho_pre[i] - rho_pre[i-1]));
            rho_r[i] = rho_l[i];

            v_l[i] = v_pre[i-1] + 0.5 * Fl_v[i-1] * (v_pre[i] - v_pre[i-1]);
            v_r[i] = v_l[i];

            E_l[i] = e_pre[i-1] + 0.5 * Fl_e[i-1] * (e_pre[i] - e_pre[i-1]);
            E_r[i] = E_l[i];

            P_l[i] = P[i-1] + 0.5 * Fl_P[i-1] * (P[i] - P[i-1]);
            P_r[i] = P_l[i];

        }else{
            rho_l[i] = rho_pre[i-1] + ((0.5 * Fl_rho[i-1]) * (rho_pre[i] - rho_pre[i-1]));
            rho_r[i] = rho_pre[i] - ((0.5 * Fl_rho[i]) * (rho_pre[i+1] - rho_pre[i]));

            v_l[i] = v_pre[i-1] + 0.5 * Fl_v[i-1] * (v_pre[i] - v_pre[i-1]);
            v_r[i] = v_pre[i] - 0.5 * Fl_v[i] * (v_pre[i+1] - v_pre[i]);

            E_l[i] = e_pre[i-1] + 0.5 * Fl_e[i-1] * (e_pre[i] - e_pre[i-1]);
            E_r[i] = e_pre[i] - 0.5 * Fl_e[i] * ( e_pre[i+1] - e_pre[i]);


            P_l[i] = P[i-1] + 0.5 * Fl_P[i-1] * (P[i] - P[i-1]);
            P_r[i] = P[i] - 0.5 * Fl_P[i] * (P[i+1] - P[i]);
        }
    
        std::cout << " i " << i << " rho_l " << rho_l[i] << " rho_r " << rho_r[i] << " v_l " << v_l[i] << " v_r " << v_r[i] << 
        " as " << as[i] << " Pm_pre " << Pm_pre[i] << std::endl; 
    }

    //Calculate the Left and right Flux values from the 2nd order conserved values
    for(int i=0; i<n_nodes; i++){
        if(i==0){
            F_rho_r[i] = rho_r[i] * v_r[i];
            F_rho_l[i] = F_rho_r[i];

            F_v_r[i] = (rho_r[i] * v_r[i] * v_r[i]) + P_r[i];
            F_v_l[i] = F_v_r[i];

            F_E_r[i] = rho_r[i] * E_r[i]*v_r[i] + P_r[i]*v_r[i];
            F_E_l[i] = F_E_r[i];

            F_u_r[i] = v_r[i];
            F_u_l[i] = F_u_r[i];

        }else if(i==n_nodes-1){
            F_rho_l[i] = rho_l[i] * v_l[i];
            F_rho_r[i] = F_rho_l[i];

            F_v_l[i] = (rho_l[i] * v_l[i] * v_l[i]) + P_l[i];
            F_v_r[i] = F_v_l[i];

            F_E_l[i] = rho_l[i] * E_l[i]*v_l[i] + P_l[i]*v_l[i];
            F_E_r[i] = F_E_l[i];

            F_u_l[i] = v_l[i];
            F_u_r[i] = F_u_l[i];

        }else{
            F_rho_l[i] = rho_l[i] * v_l[i];
            F_rho_r[i] = rho_r[i] * v_r[i];

            F_v_l[i] =  (rho_l[i] * v_l[i] * v_l[i]) + P_l[i];
            F_v_r[i] = (rho_r[i] * v_r[i] * v_r[i]) + P_r[i];

            F_E_l[i] = rho_l[i] * E_l[i]*v_l[i] + P_l[i]*v_l[i];
            F_E_r[i] = rho_r[i] * E_r[i]*v_r[i] + P_r[i]*v_r[i];

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

            Frus_u_r[i] =  0.5 * (F_u_l[i+1] + F_u_r[i+1]);
            Frus_u_l[i] = Frus_u_r[i];

        }else if(i==n_cells-1){
            Frus_rho_l[i] = 0.5 * (F_rho_l[i] + F_rho_r[i]) - (0.5 * Sp[i] * (rho_r[i] - rho_l[i]));
            Frus_rho_r[i] = Frus_rho_l[i];

            Frus_v_l[i] = 0.5 * (F_v_l[i] + F_v_r[i]) - (0.5 * Sp[i] * (rho_r[i]*v_r[i] - rho_l[i]*v_l[i]));
            Frus_v_r[i] = Frus_v_l[i];

            Frus_e_l[i] = 0.5 * (F_E_l[i] + F_E_r[i]) - (0.5 * Sp[i] * (rho_r[i]*E_r[i] - rho_l[i]*E_l[i]));
            Frus_e_r[i] = Frus_e_l[i];

            Frus_u_l[i] =  0.5 * (F_u_l[i] + F_u_r[i]);
            Frus_u_r[i] = Frus_u_l[i];

        }else{
            Frus_rho_l[i] = 0.5 * (F_rho_l[i] + F_rho_r[i]) - (0.5 * Sp[i] * (rho_r[i] - rho_l[i]));
            Frus_rho_r[i] = 0.5 * (F_rho_l[i+1] + F_rho_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1] - rho_l[i+1]));

            Frus_v_l[i] = 0.5 * (F_v_l[i] + F_v_r[i]) - (0.5 * Sp[i] * (rho_r[i]*v_r[i] - rho_l[i]*v_l[i]));
            Frus_v_r[i] = 0.5 * (F_v_l[i+1] + F_v_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1]*v_r[i+1] - rho_l[i+1]*v_l[i+1]));

            Frus_e_l[i] = 0.5 * (F_E_l[i] + F_E_r[i]) - (0.5 * Sp[i] * (rho_r[i]*E_r[i] - rho_l[i]*E_l[i]));
            Frus_e_r[i] =  0.5 * (F_E_l[i+1] + F_E_r[i+1]) - (0.5 * Sp[i+1] * (rho_r[i+1]*E_r[i+1] - rho_l[i+1]*E_l[i+1]));

            Frus_u_l[i] =  0.5 * (F_u_l[i] + F_u_r[i]);
            Frus_u_r[i] =  0.5 * (F_u_l[i+1] + F_u_r[i+1]);

        }

    }

    //Find the l(U) values
    for(int i=0; i<n_cells; i++){

        lu_rho[i] = -(Frus_rho_r[i] - Frus_rho_l[i]);
        lu_v[i] = -(Frus_v_r[i] - Frus_v_l[i]);
        lu_e[i] = -(Frus_e_r[i] - Frus_e_l[i]);

        //grad(u) term
        grad_u[i] = (Frus_u_r[i] - Frus_u_l[i]);
    }
}

void eularian_calcs_cor(std::vector<double> &v0, std::vector<double> &rho0, std::vector<double> &Pm, std::vector<double> &e0,
                        std::vector<double> &v, std::vector<double> &e, std::vector<double> &rho, std::vector<double> &lu_rho,
                        std::vector<double> &lu_v, std::vector<double> &lu_e, std::vector<double> &Pr,  std::vector<double> &grad_u,
                        std::vector<double> &Tm, std::vector<double> &Th, std::vector<double> &cv, double dt, double dx, double gamma){

    int n = v0.size();
    for(int j=0; j<n; j++){
        rho[j] = rho0[j] * dx + dt * lu_rho[j];
        rho[j] = rho[j] / dx;
        v[j] = v0[j]*rho0[j]*dx + dt*lu_v[j];
        v[j] = v[j] / (rho[j]*dx);
        e[j] = e0[j]*rho0[j]*dx + dt * (lu_e[j] + (Pr[j] * grad_u[j]));
        e[j] = e[j] / (rho[j]*dx);
        
        //Calculate internal energy
        double ie = e[j] - (0.5*(v[j] * v[j]));
        double ie_initial = e0[j] - (0.5*pow(v0[j],2));
        double delta = ie - ie_initial;

        //Calculate new pressure
        Pm[j] =  (gamma - 1) * ie * rho[j];

        //hydro temp
        Th[j] = (delta / cv[j]) + Tm[j];

        std::vector<double> as(n);
        as[j] = sqrt((gamma * Pm[j]) / rho[j]);

        std::cout.precision(17);
        std::cout << " j " << j << " rho " << rho[j] << " Pm " << Pm[j] << " as " << as[j] << " ie " << ie <<
        " v " << v[j] << " e " << e[j] <<  " lu_rho " << lu_rho[j] << std::endl;
        
    }

    std::cin.get();

}

void eularian_reassign_so(std::vector<double> &m, std::vector<double> &v0, std::vector<double> &rho0, std::vector<double> &Pm, std::vector<double> &e0,
                          std::vector<double> &v, std::vector<double> &e, std::vector<double> &rho, std::vector<double> &as, std::vector<double> &lu_rho,
                          std::vector<double> &lu_v, std::vector<double> &lu_e, std::vector<double> &grad_u, double gamma, double dx){

    int n_cells = m.size();

    //Re-assign values
    for(int i=0; i<n_cells; i++){
        v0[i] = v[i];
        rho0[i] = rho[i];
        e0[i] = e[i];
        m[i] = rho0[i] * dx;

        //sound speed in material
        as[i] = sqrt((gamma * Pm[i]) / rho0[i]);

        lu_rho[i] = 0;
        lu_v[i] = 0;
        lu_e[i] = 0; 
        grad_u[i] = 0;
        
    }


}

