#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "hydro_functions.hpp"
#include "rad_functions.hpp"

void so_eularian_rh(double dt, double ti, double tf, double dx, double xf, int n_cells){

    //Initialize timing
    const int n_time = tf/dt;
    double time = ti;
    double dt_half = dt/2;

    //value of gamma
    double mat_gamma = 5.0/3.0;
    double rad_gamma = 1.4;

    //Constants
    const int a = 137;
    double c = 3E+10; // cm/s
    const int n_iter = 1000000000;

    //constant to switch inputs and boundary conditions
    int p;

    //for printing to files
    std::ofstream myfile;

    //Flux terms
    std::vector<double> lu_rho(n_cells,0);
    std::vector<double> lu_v(n_cells,0);
    std::vector<double> lu_e(n_cells,0);
    std::vector<double> lu_w(n_cells,0);
    std::vector<double> lu_Er(n_cells);
    std::vector<double> grad_u(n_cells);

    //!Predictor Values
    //Material and Total values
    std::vector<double> v0(n_cells,0); //initial velocity in cell
    std::vector<double> v_pre(n_cells,0); //updated velocity in cell
    std::vector<double> rho0(n_cells,0); //initial density in cell
    std::vector<double> rho_pre(n_cells,0); //mass in cell
    std::vector<double> P_pre(n_cells,0); //Total pressure in cell
    std::vector<double> Pm_pre(n_cells,0); //Material Pressure in cell
    std::vector<double> e0(n_cells,0); //inital energy in cell
    std::vector<double> e_pre(n_cells,0); //updated energy in cell
    std::vector<double> as(n_cells,0); //speed of sound
    std::vector<double> Th(n_cells, 0); //Hydro Pressure


    //!Corrector Values
    //Material Values
    std::vector<double> v(n_cells,0); //velocity in cell
    std::vector<double> rho(n_cells,0); // density in cell
    std::vector<double> P(n_cells,0); //pressure in cell
    std::vector<double> Pm(n_cells,0); //material pressure in cell
    std::vector<double> e(n_cells,0); // energy in cell

    //!Radiation values
    std::vector<double> T0_r(n_cells);
    std::vector<double> T0_m(n_cells);
    std::vector<double> Es_r(n_cells);
    std::vector<double> E0_r(n_cells);
    std::vector<double> E0_m(n_cells);
    std::vector<double> Pr(n_cells,0); //radiation pressure in cell
    double E_total;

    //! Non predictor or corrector values
    std::vector<double> abs(n_cells);
    std::vector<double> emis(n_cells);
    std::vector<double> opc(n_cells);
    std::vector<double> cv(n_cells);
    std::vector<double> D(n_cells);
    std::vector<double> Dp(n_cells);
    std::vector<double> Dm(n_cells);
    std::vector<double> Ek(n_cells);
    std::vector<double> Tk(n_cells);
    std::vector<double> plus(n_cells);
    std::vector<double> mid(n_cells);
    std::vector<double> minus(n_cells); 
    std::vector<double> rs(n_cells);
    std::vector<double> Er(n_cells);
    std::vector<double> res0(n_cells,0);
    std::vector<double> res1(n_cells,0);    


    //!Initialize Values
    rad_initialize(T0_r, E0_r, E_total, T0_m, E0_m, Ek, Tk, opc, abs, emis, rho, cv, Pr, a, c, dx, xf, p=1); //p=0 mach 1.2 p=1 mach 3
    mat_initialize(v0, rho0, Pm, e0, v, e, rho, as, Pr, P, T0_m, dx, mat_gamma, xf, p=1); // p=0 mach 1.2 p=1 mach 3

    for(int i=0; i<n_time; i++){

        //Hydro Step
        reassign(Pm, v, e, rho, v0, rho0,  e0, as, mat_gamma);
        
        //Predictor Step
        flux(v0, rho0, P, e0, Er, as, lu_rho, lu_v, lu_e, lu_Er, grad_u, dt, mat_gamma);

        eularian_calcs(v0, rho0, e0, Pr, lu_rho, lu_v, lu_e, grad_u, cv, Pm_pre, v_pre, e_pre, rho_pre, P_pre, Th, dt_half, dx, mat_gamma);

        //Corrector Step
        flux(v_pre, rho_pre, P_pre, e_pre, Er, as, lu_rho, lu_v, lu_e, lu_Er, grad_u, dt, mat_gamma);

        eularian_calcs(v0, rho0, e0, Pr, lu_rho, lu_v, lu_e, grad_u, cv, Pm, v, e, rho, P, Th, dt, dx, mat_gamma);
        
        //Radiation MMC step
        rad_mmc(E0_r, Pr, grad_u, lu_Er, Es_r, dt, dx);

        //Radiation Setup
        setup(opc, T0_m, D, Dp, Dm, cv, c, a);

        //Iteration Loop for implicit temp and energy
        for(int k=0; k<n_iter; k++){

            //Newtons method to solve for temp
            newtons(Ek, Th, opc, rho, cv, Tk, dt, c, a, n_iter);

            //set up matrix
            matrix(Tk, Es_r, opc, Dp, Dm, plus, mid, minus, rs, dt, dx, c, a, p=2); //p=0 zero boundary, p=1 marshak, p=2 reflective

            //calculate residuals          
            residual(plus, mid, minus, Ek, rs, res0);

            //Thomas Algorithm to solve tridiagonal matrix
            thomas_alg(minus, mid, plus, res0, Er);

            //Check to break out of iteration loop
            double delta_max = -1;
            for(int y=0; y<n_cells; y++){
                double delta = fabs(Er[y]) / Ek[y];
                delta_max = std::max(delta_max, delta);
                Ek[y] += Er[y];
            }

            if(delta_max < 1.0E-10){
                break;
            }
    
        }

        //Energy Deposition Step
        e_dep(cv, Th, Tk, e);

        //Reassign rad values
        reassign(Tk, Ek, opc, rho, cv, E0_r, T0_r, E0_m, T0_m, abs, emis, Pr, Pm, P, a, c);

        double time = dt * i;

        std::cout << " time " << time << std::endl;
    }
    
    myfile.open("../results/rad_hydro_shock_so.dat");
    for(int i=0; i<n_cells; i++){
        myfile << i*dx << " " << T0_r[i] << " " << T0_m[i] << " " << rho[i] << " " << v[i] << " " << e[i] << " " << P[i] <<  "\n";
        
    } 
    myfile.close();
    
}