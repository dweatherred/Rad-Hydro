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

    //value of gamma
    double mat_gamma = 5.0/3.0;
    double rad_gamma = 1.4;

    //Constants
    const int a = 137;
    double c = 3E+10; // cm/s
    const int n_iter = 10000000;

    //constant to switch zero to marshak boundary
    int p = 1;

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
    std::vector<double> m(n_cells,0); //mass in cell
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
    rad_initialize(T0_r, E0_r, E_total, T0_m, E0_m, Ek, Tk, opc, abs, emis, rho, cv, Pr, a, c, dx, xf);
    mat_initialize(m,v0, rho0, Pm, e0, v, e, rho, as, Pr, P, T0_m, Th, dx, mat_gamma, xf);

    for(int i=0; i<n_time; i++){

        //Hydro Step
        eularian_reassign_so(m, v0, rho0, Pm, e0, v, e, rho, as, lu_rho, lu_v, lu_e, grad_u,mat_gamma, dx);
        
        //Calculte flux values prior to hydro and rad calcs
        flux_pre(m, v0, rho0, P, e0, v, e, rho, as, lu_rho, lu_v, lu_e, lu_Er, grad_u, E0_r, dt, mat_gamma, dx);
        
        eularian_calcs_pre(v0, rho0, Pm_pre, e0, v_pre, e_pre, rho_pre, lu_rho, lu_v, lu_e, Pr, grad_u, dt, dx, mat_gamma);

        flux_cor(m, P, Pm_pre, v_pre, e_pre, rho_pre, as, lu_rho, lu_v, lu_e, lu_Er, grad_u, dt, mat_gamma, dx);

        eularian_calcs_cor(v0, rho0, Pm, e0, v, e, rho, lu_rho, lu_v, lu_e, Pr, grad_u, T0_m, Th, cv, dt, dx, mat_gamma);
        
        //Radiation MMC step
        rad_mmc(Es_r, E0_r, Pr, grad_u, lu_Er, dt, dx);

        //Radiation Solve
        setup(opc, Th, D, Dp, Dm, cv, c, a);

        //Iteration Loop for implicit temp and energy
        for(int k=0; k<n_iter; k++){

            //Newtons method to solve for temp
            newtons(Ek, Tk, E0_r, Th, opc, dt, rho, cv, c, a, n_iter);

            //set up matrix
            matrix(plus, mid, minus, rs, Tk, Es_r, opc, Dp, Dm, dt, dx, c, a, p=2); //p=0 zero boundary, p=1 marshak, p=2 reflective

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

            if(delta_max < 1.0E-3){
                break;
            }
    
        }

        //Energy Deposition Step
        e_dep(e, cv, Th, Tk);

        //Reassign rad values
        reassign(E0_r, T0_r, E0_m, T0_m, Tk, Ek, abs, emis, opc, rho, cv, Pr, Pm, P, a, c);

        double time = dt * i;

        std::cout << " time " << time << std::endl;
    }
    
    myfile.open("../results/rad_hydro_shock_so.dat");
    for(int i=0; i<n_cells; i++){
        myfile << i*dx << " " << T0_r[i] << " " << T0_m[i] << " " << rho[i] << " " << v[i] << " " << e[i] << " " << P[i] <<  "\n";
        
    } 
    myfile.close();
    
    /*

    myfile.open("../results/hydro_test.dat");
    for(int i=0; i<n_cells; i++){
        myfile << (i*dx) << " " << e[i] << " " << rho[i] << " " << P[i] << " " << v[i] << "\n";
    }
    myfile.close();
    */
}