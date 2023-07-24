#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

//!Initialize Radiaition and Material Properties
double evaluate_sigma_a(const double T, const double rho){
    
    return pow(rho,2) * 4.494E+8/pow(T,3.5);  //pow(T, 3); //577.35; //788; // 1.0E+6/ pow(T, 3); pow(rho,2) * 4.494E+8/pow(T,3.5);
}

double evaluate_sigma_s(const double T, const double rho){
    
    return rho * 0.4006; //pow(T, 3); //577.35; //788; // 1.0E+6/ pow(T, 3);
}

double evaluate_cv(const double T, double a){
    
    return 1.4467E+12; //4 * a * pow(T, 3); 1.3487E+11; 1.4467E+12;
}


void rad_initialize(std::vector<double> &T0_r, std::vector<double> &E0_r, double E_total,
                    std::vector<double> &T0_m, std::vector<double> &E0_m, std::vector<double> &Ek,
                    std::vector<double> &Tk, std::vector<double> &sigma_a, std::vector<double> &sigma_s, 
                    std::vector<double> &sigma_t, std::vector<double> &abs, std::vector<double> &emis,
                    std::vector<double> &rho, std::vector<double> &cv, std::vector<double> &Pr, 
                    std::vector<double> &Er_in, double a, double c, double dx, double xf, int p){

    int n = abs.size();

    for(int j=0; j<n; j++){
        //Temperature
        if(p==0){//Mach 1.2
            if((j+0.5)*dx < xf/2){
                T0_r[j] = 100;
                T0_m[j] = 100;
            }else{
                T0_r[j] = 119.476;
                T0_m[j] = 119.476;
            }
        }else if(p==1){//Mach 3
            if((j+0.5)*dx < xf/2){
                T0_r[j] = 100;
                T0_m[j] = 100;
            }else{
                T0_r[j] = 366.260705;
                T0_m[j] = 366.260705;
            }
        }else if(p==2){//Mach 45
            if((j+0.5)*dx < xf/2){
                T0_r[j] = 100;
                T0_m[j] = 100;
            }else{
                T0_r[j] = 8.357785413124984E+03;
                T0_m[j] = 8.357785413124984E+03;
            }
        }else if(p==3){//Sod Shock Tube
            T0_r[j] = 0;
            T0_m[j] = 0;
        }else if(p==4){//Marshak
            T0_r[j] = 0.025;
            T0_m[j] = 0.025;
        }
        //Initialize specific heat capacity in each cell
        cv[j] = evaluate_cv(T0_m[j], a);

        //Energy
        E0_r[j] = (a * pow(T0_r[j],4)); // erg/cc
        E0_m[j] = rho[j] * cv[j] * T0_m[j]; // erg/cc
        E_total = E0_r[j] + (rho[j] * cv[j] * T0_m[j]);
        Er_in[j] = E0_r[j];

        //Inital opacity in each cell
        sigma_a[j] = evaluate_sigma_a(T0_m[j], rho[j]);
        sigma_s[j] = evaluate_sigma_a(T0_m[j], rho[j]);
        sigma_t[j] = sigma_s[j] + sigma_a[j];

        //Ek and Tk for Newton's Method
        Ek[j] = E0_r[j];
        Tk[j] = T0_m[j];

        //Radiation Pressure
        if (p==3){
            Pr[j] = 0;
        }
        else{
            Pr[j] = E0_r[j]/3.0; 
        }
    }
}
void moving_shock_rad_init(std::vector<double> &T0_r, std::vector<double> &E0_r, double E_total,
                           std::vector<double> &T0_m, std::vector<double> &E0_m, std::vector<double> &Ek,
                           std::vector<double> &Tk, std::vector<double> &sigma_a, std::vector<double> &sigma_s, 
                           std::vector<double> &sigma_t, std::vector<double> &abs, std::vector<double> &emis, 
                           std::vector<double> &rho, std::vector<double> &cv, std::vector<double> &Pr, 
                           std::vector<double> &Er_in, double a, double c, double dx, double xf){
    int n = abs.size();
    int q=0;

    for(int j=0; j<n; j++){
        if((j+0.5)*dx < xf/2){
                T0_r[j] = 100;
                T0_m[j] = 100;
            }else{
                T0_r[j] = 120.5;
                T0_m[j] = 120.5;
        }

        //Initialize specific heat capacity and opacitu in each cell
        cv[j] = evaluate_cv(T0_m[j], a);
        sigma_a[j] = evaluate_sigma_a(T0_m[j], rho[j]);
        sigma_s[j] = evaluate_sigma_a(T0_m[j], rho[j]);
        sigma_t[j] = sigma_s[j] + sigma_a[j];
        //Energy
        E0_r[j] = (a * pow(T0_r[j],4)); // erg/cc
        E0_m[j] = rho[j] * cv[j] * T0_m[j]; // erg/cc
        E_total = E0_r[j] + (rho[j] * cv[j] * T0_m[j]);
        Er_in[j] = E0_r[j];

        //Ek and Tk for Newton's Method
        Ek[j] = E0_r[j];
        Tk[j] = T0_m[j];
        //Radiation Pressure
        Pr[j] = E0_r[j]/3.0; 
    }
}

//! Radiation MMC
void rad_mmc(const std::vector<double> &E0_r, const std::vector<double> &Pr, const std::vector<double> &grad_u,
             const std::vector<double> &lu_Er, std::vector<double> &Es_r, double dt, double dx){

    int n_cells = E0_r.size();

    for(int i=0; i<n_cells; i++){
        Es_r[i] = E0_r[i] - ((dt/dx) * (lu_Er[i] + Pr[i]*grad_u[i]));
    }

}

//!Radiation Energy Density and Material Temperature Solve
// Calculate diffusion coefficients
void setup(std::vector<double> &sigma_a, std::vector<double> &sigma_s, std::vector<double> &sigma_t, 
           const std::vector<double> &T0_m, std::vector<double> &D, std::vector<double> &Dp, 
           std::vector<double> &Dm, std::vector<double> &cv, std::vector<double> &rho, const double c, 
           const double a, int p){

    int n=sigma_a.size();

    //Setup Loop (Diffusion Coefficients and Opacity), Constant over timestep
    for(int i=0; i<n; i++){
        sigma_a[i] = evaluate_sigma_a(T0_m[i], rho[i]);
        sigma_s[i] = evaluate_sigma_a(T0_m[i], rho[i]);
        sigma_t[i] = sigma_s[i] + sigma_a[i];
        cv[i] = evaluate_cv(T0_m[i], a);
        D[i] = c /(3 * sigma_t[i]);
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
                sigma_a[l] = evaluate_sigma_a(Tf, rho[l]);
                sigma_s[l] = evaluate_sigma_a(Tf, rho[l]);
                sigma_t[l] = sigma_s[l] + sigma_a[l];
                Dp[l] = c/(3 * sigma_t[l]);
            }

            if(l==0){
                Dm[l] = D[l];
            }else{
                double Tf = (T0_m[l] + T0_m[l-1]) * 0.5;
                sigma_a[l] = evaluate_sigma_a(Tf, rho[l]);
                sigma_s[l] = evaluate_sigma_a(Tf, rho[l]);
                sigma_t[l] = sigma_s[l] + sigma_a[l];
                Dm[l] = c/(3 * sigma_t[l]);
            }
        }
    }
}

//Newtons Method
void newtons(const std::vector<double> &Ek, const std::vector<double> &Th, const std::vector<double> &sigma_a, const std::vector<double> &rho, 
             const std::vector<double> &cv, std::vector<double> &Tk, const double dt, const double c, const double a, const int iter){

    int n=Ek.size();
    std::vector<double> dTk(n); std::vector<double> Tk1(n);
    std::vector<double> FT(n); std::vector<double> dTT(n);

    for(int i=0; i<n; i++){
        //newton step
        for(int k=0; k<iter; k++){

            //Calculate FT with "known" Ek and guess "Tk"
            FT[i] = ((rho[i] * cv[i]) * ((Tk[i] - Th[i])/dt)) - (sigma_a[i] * c * Ek[i]) + (sigma_a[i] * a * c * pow(Tk[i],4));

            //Derivative of FT with respect to Tk
            dTT[i] = ((rho[i] * cv[i]) / dt) + (4 * sigma_a[i] * c * a * pow(Tk[i], 3));

            //Solve for dTk
            dTk[i] = -FT[i] / dTT[i];

            //Solve for Tk+1
            Tk1[i] = Tk[i] + dTk[i];

            //Reassign Tk as Tk+1
            Tk[i] = Tk1[i];

            if(fabs(dTk[i]) < 1.0E-5){
                break;
            }
        }
    }
}

//Matrix setup for Implicit 1-D problem with diffusion

void matrix(const std::vector<double> &Tk, const std::vector<double> &Es_r, const std::vector<double> &sigma_a, 
            const std::vector<double> &Dp, const std::vector<double> &Dm, std::vector<double> &plus, 
            std::vector<double> &mid, std::vector<double> &minus, std::vector<double> &rs, 
            double dt, const double dx, const double c, const double a, int p=0){

    int n=plus.size();

    if(p==0){ //zero boundary condition
        for(int i=0; i<n; i++){
            //upper diagonal
            plus[i] = -Dp[i] / pow(dx,2); 
            //main Diagonal
            mid[i] = 1/dt + Dp[i]/pow(dx,2) + Dm[i]/pow(dx,2) + (sigma_a[i] * c);
            //lower diagonal
            minus[i] = -Dm[i] / pow(dx,2);
            //Right Side of Equation
            rs[i] = (sigma_a[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);            
        }
    }else if (p==1) {//Marshak wave boundary conditon
        for(int j=0; j<n; j++){
            //upper diagonal
            plus[j] = -Dp[j] / pow(dx,2); 
            //main Diagonal
            if(j==0){
                mid[j] = 1/dt + c/(2*dx) + Dp[j]/pow(dx,2) + (sigma_a[j] * c);
            }else{
                mid[j] = 1/dt + Dp[j]/pow(dx,2) + Dm[j]/pow(dx,2) + (sigma_a[j] * c);
            }
            //lower diagonal
            minus[j] = -Dm[j] / pow(dx,2);
            //Right Side of Equation
            if(j==0){
                rs[j] = sigma_a[j] * a * c * pow(Tk[j],4) + ( c * a * pow(150,4) * 0.5 )/dx  + Es_r[j]/dt;
            }else{
                rs[j] = (sigma_a[j] * a * c * pow(Tk[j],4)) + (Es_r[j]/dt);
            }
        }
    }else if(p==2){ //Reflected Boundary Condition
        for(int i=0; i<n; i++){
            if(i==0){
                //upper diagonal
                plus[i] = -Dp[i] / pow(dx,2); 
                //main Diagonal
                mid[i] = 1/dt + Dp[i]/pow(dx,2) + (sigma_a[i] * c);
                //lower diagonal
                minus[i] = 0;
                //Right Side of Equation
                rs[i] = (sigma_a[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);
            }else if (i==n-1){
                //upper diagonal
                plus[i] = 0; 
                //main Diagonal
                mid[i] =  1/dt + Dm[i]/pow(dx,2) + (sigma_a[i] * c);
                //lower diagonal
                minus[i] =  -Dm[i] / pow(dx,2);
                //Right Side of Equation
                rs[i] = (sigma_a[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);
            }else{
                //upper diagonal
                plus[i] = -Dp[i] / pow(dx,2); 
                //main Diagonal
                mid[i] = 1/dt + Dp[i]/pow(dx,2) + Dm[i]/pow(dx,2) + (sigma_a[i] * c);
                //lower diagonal
                minus[i] = -Dm[i] / pow(dx,2);
                //Right Side of Equation
                rs[i] = (sigma_a[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);
            }
        }
    }else if(p==3){ //inflow outflow
        double T_l_incident = 100;
        for(int i=0; i<n; i++){
            //inflow to left boundary
            if(i==0){
                //upper diagonal
                plus[i] = -Dp[i] / pow(dx,2); 
                //middle
                mid[i] = 1/dt + c/(2*dx) + Dp[i]/pow(dx,2) + (sigma_a[i] * c);
                //lower diagonal
                minus[i] = -Dm[i] / pow(dx,2);
                //right hand side
                rs[i] = sigma_a[i] * a * c * pow(Tk[i],4) + (c * a * pow(T_l_incident,4) * 0.5 )/dx  + Es_r[i]/dt;
            }else if(i==n-1){
                 //upper diagonal
                plus[i] = 0; 
                //main Diagonal
                mid[i] =  1/dt + Dm[i]/pow(dx,2) + (sigma_a[i] * c);
                //lower diagonal
                minus[i] =  -Dm[i] / pow(dx,2);
                //Right Side of Equation
                rs[i] = (sigma_a[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);
            }else{
                //upper diagonal
                plus[i] = -Dp[i] / pow(dx,2); 
                //main Diagonal
                mid[i] = 1/dt + Dp[i]/pow(dx,2) + Dm[i]/pow(dx,2) + (sigma_a[i] * c);
                //lower diagonal
                minus[i] = -Dm[i] / pow(dx,2);
                //Right Side of Equation
                rs[i] = (sigma_a[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);
            }
        }

    }
}

void residual(const std::vector<double> &plus, const std::vector<double> &mid, const std::vector<double> &minus,
              const std::vector<double> &Er, const std::vector<double> &rs, std::vector<double> &res0)
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
void e_dep(const std::vector<double> &cv, const std::vector<double> &Th, const std::vector<double> &Tk, std::vector<double> &e){

    int n=e.size();
    for(int i=0; i<n; i++){
        e[i] = e[i] + (cv[i] * (Tk[i]- Th[i]));
    }
}

//Post Processing
void reassign(const std::vector<double> &Tk, const std::vector<double> &Ek, const std::vector<double> &rho,
              const std::vector<double> &cv, std::vector<double> &E0_r, std::vector<double> &T0_r, 
              std::vector<double> &E0_m, std::vector<double> &T0_m, std::vector<double> &Pr, 
              std::vector<double> &Pm, std::vector<double> &P, const double a, const double c, int p){

    int n=E0_r.size(); 

    for(int q=0; q<n; q++){
    //Radiation Terms
        E0_r[q] = Ek[q];
        T0_r[q] = pow((E0_r[q]/a), 0.25);
        //Material Terms
        T0_m[q] = Tk[q];
        E0_m[q] = rho[q] * cv[q] * T0_m[q];
        
        //Radiation Pressure
        if (p==2){
            Pr[q] = 0;
        }
        else{
            Pr[q] = E0_r[q]/3.0; 
        }

        //Total Pressure
        P[q] = Pr[q] + Pm[q];
    }
}