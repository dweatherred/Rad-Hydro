#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

//!Initialize Radiaition and Material Properties
double evaluate_opc(const double T){
    
    return  423; //788; // 1.0E+6/ pow(T, 3); 5.7735E+6;
}


double evaluate_cv(const double T, double a){
    
    return  1.4467E+12;  //4 * a * pow(T, 3); 1.3487E+11; 1.4467E+12;
}

void rad_initialize(std::vector<double> &T0_r, std::vector<double> &E0_r, double E_total,
                    std::vector<double> &T0_m, std::vector<double> &E0_m, std::vector<double> &Ek,
                    std::vector<double> &Tk, std::vector<double> &opc, std::vector<double> &abs,
                    std::vector<double> &emis, std::vector<double> &rho, std::vector<double> &cv,
                    std::vector<double> &Pr, double a, double c, double dx, double xf, double gamma){

    int n = abs.size();


    for(int j=0; j<n; j++){

        //Temperature
        if((j + 0.5)*dx < xf/2){
            T0_r[j] = 100;
            T0_m[j] = 100;
        }else{
            T0_r[j] = 207.8;
            T0_m[j] = 207.8;

        }
        //Initialize specific heat capacity in each cell
        cv[j] = evaluate_cv(T0_m[j], a);

        //Energy
        E0_r[j] = (a * pow(T0_r[j],4)); // erg/cc
        E0_m[j] = rho[j] * cv[j] * T0_m[j]; // erg/cc         
        E_total = E0_r[j] + (rho[j] * cv[j] * T0_m[j]);

        //Inital opacity in each cell
        opc[j] = evaluate_opc(T0_m[j]);

        //Abs and Emis
        abs[j] = opc[j] * c * E0_r[j];
        emis[j] = opc[j] * a * c * pow(T0_m[j],4);

        //Ek and Tk for Newton's Method
        Ek[j] = E0_r[j];
        Tk[j] = T0_m[j];

        //Radiation Pressure
        Pr[j] = E0_r[j]/3;
    }
}


void mat_initialize(std::vector<double> &m,std::vector<double> &v0, std::vector<double> &rho0, std::vector<double> &Pm, std::vector<double> &e0,
                    std::vector<double> &v, std::vector<double> &e, std::vector<double> &rho, std::vector<double> &as, std::vector<double> &Pr,
                    std::vector<double> &P, std::vector<double> &T0_m, double dx, double gamma, double xf){  

    //Number of cells
    int n_cells = v0.size();

    for(int i=0; i<n_cells; i++){
        
        if((i + 0.5)*dx < xf/2){
            rho0[i] = 1;
            rho[i] = rho0[i];
            v0[i] = 2.536E+7;
            v[i] = v0[i];
            m[i] = rho0[i] * dx;
            Pm[i] = (rho[i] * T0_m[i])/gamma;
            e0[i] = Pm[i] / (rho0[i] * (gamma - 1));
            e[i] = e0[i];
        }else{
            rho0[i] = 2.286;
            rho[i] = rho0[i];
            v0[i] = 1.11E+7; 
            v[i] = v0[i];
            m[i] = rho0[i] * dx;
            Pm[i] = (rho[i] * T0_m[i])/gamma;
            e0[i] = Pm[i] / (rho0[i] * (gamma - 1));
            e[i] = e0[i];
        }

        as[i] = sqrt((gamma * Pm[i]) / rho0[i]);

        //Total Pressure
        P[i] = Pr[i] + Pm[i];

    }
}

//!Flux Calculations
void calc_flux(std::vector<double> &m, std::vector<double> &v0, std::vector<double> &rho0, std::vector<double> &P, std::vector<double> &e0,
               std::vector<double> &as, std::vector<double> &E0_r, std::vector<double> &lu_rho, std::vector<double> &lu_v, std::vector<double> &lu_e, 
               std::vector<double> &lu_Er, std::vector<double> &grad_u){

    int n_cells = m.size();
    int n_nodes = n_cells + 1;

    // left and right Material Flux values
    std::vector<double> F_rho_l(n_nodes,0);
    std::vector<double> F_rho_r(n_nodes,0);
    std::vector<double> F_v_l(n_nodes,0);
    std::vector<double> F_v_r(n_nodes,0);
    std::vector<double> F_E_l(n_nodes,0);
    std::vector<double> F_E_r(n_nodes,0);

    //Rusanov Material flux values
    std::vector<double> Frus_rho_l(n_cells,0);
    std::vector<double> Frus_rho_r(n_cells,0);
    std::vector<double> Frus_v_l(n_cells,0);
    std::vector<double> Frus_v_r(n_cells,0);
    std::vector<double> Frus_e_l(n_cells,0);
    std::vector<double> Frus_e_r(n_cells,0);

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

    for(int i=0; i<n_nodes; i++){
        //Find +- flux values
        if(i==0){
            F_rho_r[i] = rho0[i] * (v0[i]);
            F_rho_l[i] = F_rho_r[i];

            F_v_r[i] = (rho0[i] * v0[i] * v0[i]) + P[i];
            F_v_l[i] = F_v_r[i];

            F_E_r[i] = (rho0[i] * e0[i] * v0[i]) + (P[i]*v0[i]);
            F_E_l[i] = F_E_r[i];

            F_Er_r[i] = E0_r[i] * v0[i];
            F_Er_l[i] = F_Er_r[i];

            F_u_r[i] = v0[i];
            F_u_l[i] = F_u_r[i];

        }else if(i==n_nodes-1){
            F_rho_l[i] = rho0[i-1] * (v0[i-1]);
            F_rho_r[i] = F_rho_l[i];

            F_v_l[i] = (rho0[i-1] * v0[i-1] * (v0[i-1])) + P[i-1];
            F_v_r[i] = F_v_l[i];

            F_E_l[i] = (rho0[i-1] * e0[i-1]*(v0[i-1]) + P[i-1]*v0[i-1]);
            F_E_r[i] = F_E_l[i];

            F_Er_l[i] = E0_r[i-1] * v0[i-1];
            F_Er_r[i] = F_Er_l[i];

            F_u_l[i] = v0[i-1];
            F_u_r[i] = F_u_l[i];            

        }else{
            F_rho_l[i] = rho0[i-1] * (v0[i-1]);
            F_rho_r[i] = rho0[i] * (v0[i]);

            F_v_l[i] = (rho0[i-1] * v0[i-1] * (v0[i-1])) + P[i-1];
            F_v_r[i] = (rho0[i] * v0[i] * (v0[i])) + P[i];

            F_E_l[i] = (rho0[i-1]*e0[i-1]*v0[i-1] + P[i-1]*v0[i-1]);
            F_E_r[i] = (rho0[i] * e0[i]*(v0[i]) + P[i]*v0[i]);

            F_Er_l[i] = E0_r[i-1] * v0[i-1];
            F_Er_r[i] = E0_r[i] * v0[i];

            F_u_l[i] = v0[i-1];
            F_u_r[i] = v0[i];            
        }
    }

    for(int i=0; i<n_cells; i++){
        if(i==0){
            Frus_rho_r[i] = 0.5 * (F_rho_l[i+1] + F_rho_r[i+1]) - (0.5 * Sp[i+1] * (rho0[i+1] - rho0[i]));
            Frus_rho_l[i] = Frus_rho_r[i];

            Frus_v_r[i] = 0.5 * (F_v_l[i+1] + F_v_r[i+1]) - (0.5 * Sp[i+1] * (rho0[i+1]*v0[i+1] - rho0[i]*v0[i]));
            Frus_v_l[i] = Frus_v_r[i];

            Frus_e_r[i] =  0.5 * (F_E_l[i+1] + F_E_r[i+1]) - (0.5 * Sp[i+1] * (rho0[i+1]*e0[i+1] - rho0[i]*e0[i]));
            Frus_e_l[i] = Frus_e_r[i];

            Frus_Er_r[i] =  0.5 * (F_Er_l[i+1] + F_Er_r[i+1]) - (0.5 * Sp[i+1] * (E0_r[i+1] - E0_r[i]));
            Frus_Er_l[i] = Frus_Er_r[i];

            Frus_u_r[i] =  0.5 * (F_u_l[i+1] + F_u_r[i+1]);
            Frus_u_l[i] = Frus_u_r[i];

        }else if(i==n_cells-1){
            Frus_rho_l[i] = 0.5 * (F_rho_l[i] + F_rho_r[i]) - (0.5 * Sp[i] * (rho0[i] - rho0[i-1]));
            Frus_rho_r[i] = Frus_rho_l[i];

            Frus_v_l[i] = 0.5 * (F_v_l[i] + F_v_r[i]) - (0.5 * Sp[i] * (rho0[i] * v0[i] - rho0[i-1] * v0[i-1]));
            Frus_v_r[i] = Frus_v_l[i];

            Frus_e_l[i] = 0.5 * (F_E_l[i] + F_E_r[i]) - (0.5 * Sp[i] * (rho0[i] * e0[i] - rho0[i-1]*e0[i-1]));
            Frus_e_r[i] = Frus_e_l[i];

            Frus_Er_l[i] = 0.5 * (F_Er_l[i] + F_Er_r[i]) - (0.5 * Sp[i] * (E0_r[i] - E0_r[i-1]));
            Frus_Er_r[i] = Frus_Er_l[i];

            Frus_u_l[i] =  0.5 * (F_u_l[i] + F_u_r[i]);
            Frus_u_r[i] = Frus_u_l[i];
        
        }else{
            Frus_rho_l[i] = 0.5 * (F_rho_l[i] + F_rho_r[i]) - (0.5 * Sp[i] * (rho0[i] - rho0[i-1]));
            Frus_rho_r[i] = 0.5 * (F_rho_l[i+1] + F_rho_r[i+1]) - (0.5 * Sp[i+1] * (rho0[i+1] - rho0[i]));

            Frus_v_l[i] = 0.5 * (F_v_l[i] + F_v_r[i]) - (0.5 * Sp[i] * (rho0[i]*v0[i] - rho0[i-1]*v0[i-1]));
            Frus_v_r[i] = 0.5 * (F_v_l[i+1] + F_v_r[i+1]) - (0.5 * Sp[i+1] * (rho0[i+1]*v0[i+1] - rho0[i]*v0[i]));

            Frus_e_l[i] = 0.5 * (F_E_l[i] + F_E_r[i]) - (0.5 * Sp[i] * (rho0[i] * e0[i] - rho0[i-1]*e0[i-1]));
            Frus_e_r[i] =  0.5 * (F_E_l[i+1] + F_E_r[i+1]) - (0.5 * Sp[i+1] * (rho0[i+1]*e0[i+1] - rho0[i]*e0[i]));

            Frus_Er_l[i] = 0.5 * (F_Er_l[i] + F_Er_r[i]) - (0.5 * Sp[i] * (E0_r[i] - E0_r[i-1]));
            Frus_Er_r[i] =  0.5 * (F_Er_l[i+1] + F_Er_r[i+1]) - (0.5 * Sp[i+1] * (E0_r[i+1] - E0_r[i]));

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

//!Eularian hydro calcs and reassign

void eularian_calcs(std::vector<double> &v0, std::vector<double> &rho0, std::vector<double> &Pm, std::vector<double> &e0,
                    std::vector<double> &v, std::vector<double> &e, std::vector<double> &rho, std::vector<double> &lu_rho, 
                    std::vector<double> &lu_v, std::vector<double> &lu_e, std::vector<double> &Pr,  std::vector<double> &grad_u, 
                    double dt, double dx, double gamma){

    int n = v0.size();
    for(int j=0; j<n; j++){
        rho[j] = rho0[j] * dx + dt * lu_rho[j];
        rho[j] = rho[j] / dx;
        v[j] = v0[j]*rho0[j]*dx + dt*lu_v[j];
        v[j] = v[j] / (rho[j]*dx);
        e[j] = e0[j]*rho0[j]*dx + dt * (lu_e[j] + Pr[j] * grad_u[j]);
        e[j] = e[j] / (rho[j]*dx);

        //Calculate internal energy
        double ie = e[j] - (0.5*pow(v[j],2));

        //Calculate new pressure
        Pm[j] =  (gamma - 1) * ie * rho[j];
        //std::cout << " j " << j << " ie " << ie << " e " << e[j] << " v " << v[j] << " Pm " << Pm[j] << std::endl;

    }
}

void eularian_reassign(std::vector<double> &m, std::vector<double> &v0, std::vector<double> &rho0, std::vector<double> &Pm, std::vector<double> &e0,
                       std::vector<double> &v, std::vector<double> &e, std::vector<double> &rho, std::vector<double> &as, double gamma, double dx){

    int n_cells = m.size();

    //Re-assign values
    for(int i=0; i<n_cells; i++){
        v0[i] = v[i];
        rho0[i] = rho[i];
        e0[i] = e[i];
        m[i] = rho0[i] * dx;

        //sound speed in material
        as[i] = sqrt((gamma * Pm[i]) / rho0[i]);
    }
   
}

//! Radiation MMC

void rad_mmc(std::vector<double> &Es_r,std::vector<double> &E0_r, std::vector<double> &Pr, 
             std::vector<double> &grad_u, std::vector<double> &lu_Er, double dt, double dx){

    int n_cells = E0_r.size();

    for(int i=0; i<n_cells; i++){
        Es_r[i] = E0_r[i] - ((dt/dx) * (lu_Er[i] + Pr[i]*grad_u[i]));
    }
    
}

//!Radiation Solve

// Calculate diffusion coefficients
void setup(std::vector<double> &opc, std::vector<double> &T0_m, std::vector<double> &D, 
           std::vector<double> &Dp, std::vector<double> &Dm, std::vector<double> &cv,
           double c, double a, int p=0){

    int n=opc.size();

    //Setup Loop (Diffusion Coefficients and Opacity), Constant over timestep
    for(int i=0; i<n; i++){
        opc[i] = evaluate_opc(T0_m[i]);
        cv[i] = evaluate_cv(T0_m[i], a);
        D[i] = c /(3 * opc[i]);
    }

    if(p ==0){

        for(int l=0; l<n; l++){
            //+1/2 and -1/2 Diffusion Terms
            if(l==n-1){
                Dp[l] = D[l];
            }else{
                Dp[l] =  2 * D[l] * D[l+1] / (D[l] + D[l+1]);
            }

            if(l==0){
                Dm[l] = D[l];
            }else{
                Dm[l] = 2 * D[l] * D[l-1] / (D[l] + D[l-1]);
            }
        }

    }else{ //For marshak problem

        for(int l=0; l<n; l++){
            //+1/2 and -1/2 Diffusion Terms
            if(l==n-1){
                Dp[l] = D[l];
            }else{
                double Tf = (T0_m[l] + T0_m[l+1]) * 0.5; 
                Dp[l] = c/(3 * evaluate_opc(Tf));
            }

            if(l==0){
                Dm[l] = D[l];
            }else{
                double Tf = (T0_m[l] + T0_m[l-1]) * 0.5;
                Dm[l] = c/(3 * evaluate_opc(Tf));
            }
        }
    }
}

//Newtons Method
void newtons(std::vector<double> &Ek, std::vector<double> &Tk, std::vector<double> &E0_r, 
             std::vector<double> &T0_m, std::vector<double> &opc, const double dt, 
             const double c, const double a, std::vector<double> &rho, std::vector<double> &cv, 
             const int iter, int p=0){

    int n=Ek.size();
    std::vector<double> det(n);std::vector<double> dEk(n);std::vector<double> dTk(n);
    std::vector<double> Ek1(n);std::vector<double> Tk1(n);std::vector<double> FE(n);
    std::vector<double> FT(n);std::vector<double> dEE(n);std::vector<double> dET(n);
    std::vector<double> dTE(n);std::vector<double> dTT(n);
    double TE_mat_inv[2][2];
    double FE_m[2][2];
    double FT_m[2][2];

    if(p==0){ //This is the case where Ek is not constant
        
        //Loop over cells
        for(int i=0; i<n; i++){
            //Newtons
            for(int p = 0; p < iter; p++){
                //Calculate FE(Ek,Tk) FT(Ek,Tk)
                FE[i] = ((Ek[i] - E0_r[i])/dt) + (opc[i] * c * Ek[i]) - (opc[i] * a * c * pow(Tk[i],4));  
                FT[i] = ((rho[i] * cv[i]) * ((Tk[i] - T0_m[i])/dt)) - (opc[i] * c * Ek[i]) + (opc[i] * a * c * pow(Tk[i],4));

                //Calculate matrix values and fill matrix
                dEE[i] = (1/dt) + (opc[i] * c);
                dET[i] = -4 * opc[i] * a * c * pow(Tk[i], 3);
                dTE[i] = -opc[i] * c;
                dTT[i] = ((rho[i] * cv[i]) / dt) + (4 * opc[i] * c * a * pow(Tk[i], 3));
                double TE_mat[2][2] = {{dEE[i], dET[i]}, {dTE[i], dTT[i]}};

                //Invert Matrix
                det[i] = (dEE[i] * dTT[i]) - (dET[i] * dTE[i]);
                double TE_mat_adj[2][2] = {{TE_mat[1][1], -TE_mat[0][1]}, {-TE_mat[1][0], TE_mat[0][0]}};

                for(int x=0; x<2; x++){
                    for(int y=0; y<2; y++){
                        TE_mat_inv[x][y] =  TE_mat_adj[x][y] * (1/det[i]);
                    }
                }

                //Get dEk and dTk
                dEk[i] = -TE_mat_inv[0][0]*FE[i]-TE_mat_inv[0][1]*FT[i];
                dTk[i] = -TE_mat_inv[1][0]*FE[i]-TE_mat_inv[1][1]*FT[i];

                //solve for Ek+1 and Tk+1
                Ek1[i] = Ek[i] + dEk[i]; 
                Tk1[i] = Tk[i] + dTk[i];

                //Reassign
                Ek[i] = Ek1[i];
                Tk[i] = Tk1[i];

                if(fabs(dEk[i]) && fabs(dTk[i]) < 1.0E-3){
                    break;
                }
            }
        }
    
    }else{ //This is the case where Ek is constant (implicit method)
        
        //loop over cells
        for(int i=0; i<n; i++){
            //newton step
            for(int k=0; k<iter; k++){

                //Calculate FT with "known" Ek and guess "Tk"
                FT[i] = ((rho[i] * cv[i]) * ((Tk[i] - T0_m[i])/dt)) - (opc[i] * c * Ek[i]) + (opc[i] * a * c * pow(Tk[i],4));

                //Derivative of FT with respect to Tk
                dTT[i] = ((rho[i] * cv[i]) / dt) + (4 * opc[i] * c * a * pow(Tk[i], 3));

                //Solve for dTk
                dTk[i] = -FT[i] / dTT[i];

                //Solve for Tk+1
                Tk1[i] = Tk[i] + dTk[i];

                //Reassign Tk as Tk+1
                Tk[i] = Tk1[i];

                if(fabs(dTk[i]) < 1.0E-3){
                    break;
                }
            }
        }
    }

}

//Matrix setup for Implicit 1-D problem with diffusion

void matrix(std::vector<double> &plus, std::vector<double> &mid, std::vector<double> &minus, std::vector<double> &rs, 
            std::vector<double> &Tk, std::vector<double> &Es_r, std::vector<double> &opc, 
            std::vector<double> &Dp, std::vector<double> &Dm, double dt, const double dx, const double c, const double a, 
            int p=0){

    int n=plus.size();

    if(p==0){ //zero boundary condition
        
        for(int i=0; i<n; i++){

            //upper diagonal
            plus[i] = -Dp[i] / pow(dx,2); 

            //main Diagonal
            mid[i] = 1/dt + Dp[i]/pow(dx,2) + Dm[i]/pow(dx,2) + (opc[i] * c);

            //lower diagonal
            minus[i] = -Dm[i] / pow(dx,2);

            //Right Side of Equation
            rs[i] = (opc[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);
            
        }

    }else if (p==1) {//Marshak wave boundary conditon
        for(int j=0; j<n; j++){
            //upper diagonal
            plus[j] = -Dp[j] / pow(dx,2); 

            //main Diagonal
            if(j==0){
                mid[j] = 1/dt + c/(2*dx) + Dp[j]/pow(dx,2) + (opc[j] * c);
            }else{
                mid[j] = 1/dt + Dp[j]/pow(dx,2) + Dm[j]/pow(dx,2) + (opc[j] * c);
            }

            //lower diagonal
            minus[j] = -Dm[j] / pow(dx,2);

            //Right Side of Equation
            if(j==0){
                rs[j] = opc[j] * a * c * pow(Tk[j],4) + ( c * a * pow(150,4) * 0.5 )/dx  + Es_r[j]/dt;
            }else{
                rs[j] = (opc[j] * a * c * pow(Tk[j],4)) + (Es_r[j]/dt);
            }
        }
    }else{
        for(int i=0; i<n; i++){
            if(i==0){
                //upper diagonal
                plus[i] = -Dp[i] / pow(dx,2); 

                //main Diagonal
                mid[i] = 1/dt + Dp[i]/pow(dx,2) + (opc[i] * c);

                //lower diagonal
                minus[i] = 0;

                //Right Side of Equation
                rs[i] = (opc[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);
            }else if (i==n-1){

                //upper diagonal
                plus[i] = 0; 

                //main Diagonal
                mid[i] =  mid[i] = 1/dt + Dm[i]/pow(dx,2) + (opc[i] * c);

                //lower diagonal
                minus[i] =  -Dm[i] / pow(dx,2);

                //Right Side of Equation
                rs[i] = (opc[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);
            }else{
                //upper diagonal
                plus[i] = -Dp[i] / pow(dx,2); 

                //main Diagonal
                mid[i] = 1/dt + Dp[i]/pow(dx,2) + Dm[i]/pow(dx,2) + (opc[i] * c);

                //lower diagonal
                minus[i] = -Dm[i] / pow(dx,2);

                //Right Side of Equation
                rs[i] = (opc[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);
            }
        }
    }
}

void residual(std::vector<double> &plus, std::vector<double> &mid, std::vector<double> &minus,
              std::vector<double> &Er, std::vector<double> &rs, std::vector<double> &res0)
{

    int n=plus.size();
    for(int j=0; j<n; j++){
                
      if( j==0)
      {
        double res(0);
        res = mid[j]*Er[j] + plus[j]*Er[j+1] - rs[j];
        res0[j] = -res;
        
//        std::cout <<j << " res: " << res<< " " << Er[j]
//                  << " rs: " << rs[j]<<std::endl;
        
      }
      else if(j==n-1)
      {
        double res(0);
        res = minus[j]*Er[j-1]+ mid[j]*Er[j]  - rs[j];
        res0[j] = -res;
        
//        std::cout <<j << " res: " << res<< " " << Er[j]<<std::endl;
      }
      else
      {
        double res(0);
        res = minus[j]*Er[j-1]+ mid[j]*Er[j] + plus[j]*Er[j+1] - rs[j];
        res0[j] = -res;
//        std::cout <<j << " res: " << res<< " " << Er[j]<< " " << rs[j]<<std::endl;
      }
    }//j
    
}

//Thomas Algorithm
void thomas_alg(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, 
                std::vector<double> &r, std::vector<double> &sol){
    
    int nx = a.size();
    std::vector<double> tmp(nx,0);
    std::vector<double> tmp2(nx,0);

    tmp[0] = c[0]/b[0];
    for(int i=1;i<nx-1;++i){
        tmp[i] = c[i]/(b[i]-a[i]*tmp[i-1]);
    }

    tmp2[0] = r[0]/b[0];
    for(int i=1;i<nx;++i)
    {
        tmp2[i] = (r[i]-a[i]*tmp2[i-1])/(b[i]-a[i]*tmp[i-1]);
    }
    
    sol[nx-1] = tmp2[nx-1];
    for(int i=nx-2;i>=0;--i){
        sol[i] = tmp2[i]-tmp[i]*sol[i+1];
    }
}

//Energy deposition step
void e_dep(std::vector<double> &e, std::vector<double> &cv, std::vector<double> &T0_m, std::vector<double> &Tk){

    int n=e.size();

    for(int i=0; i<n; i++){

        if(Tk[i]- T0_m[i] < 0){
            e[i] = e[i];
        }else{
            e[i] = e[i] + (cv[i] * (Tk[i]- T0_m[i]));
        }
    }

}

//Post Processing
void reassign(std::vector<double> &E0_r, std::vector<double> &T0_r, std::vector<double> &E0_m, 
              std::vector<double> &T0_m, std::vector<double> &Tk, std::vector<double> &Ek,
              std::vector<double> &abs, std::vector<double> &emis, std::vector<double> &opc, 
              std::vector<double> &rho, std::vector<double> &cv, std::vector<double> &Pr, 
              std::vector<double> &Pm, std::vector<double> &P, double a, double c, double gamma){

    int n=E0_r.size(); 

    for(int q=0; q<n; q++){
    //Radiation Terms
        E0_r[q] = Ek[q];
        T0_r[q] = pow((E0_r[q]/a), 0.25);
        //Material Terms
        T0_m[q] = Tk[q];
        E0_m[q] = rho[q] * cv[q] * T0_m[q];

        //Absorption and Emission Terms
        abs[q] = opc[q] * c * E0_r[q];
        emis[q] = opc[q] * a * c * pow(T0_m[q], 4); 
        
        //Radiation Pressure
        Pr[q] = E0_r[q] / 3;

        //Total Pressure
        P[q] = Pr[q] + Pm[q];
    }
}