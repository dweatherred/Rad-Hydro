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
    double S;

    //Constants
    const int a = 137;
    double c = 3E+10; // cm/s
    const int n_iter = 1000000;

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

    //!Initial Values That dont get edited
    std::vector<double> v_in(n_cells,0); //velocity in cell
    std::vector<double> rho_in(n_cells,0); // density in cell
    std::vector<double> P_in(n_cells,0); // total pressure in cell
    std::vector<double> e_in(n_cells,0); // energy in cell
    std::vector<double> Er_in(n_cells,0); // Radiation energy in cell

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
    std::vector<double> sigma_a(n_cells);
    std::vector<double> sigma_s(n_cells);
    std::vector<double> sigma_t(n_cells);
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
    rad_initialize(T0_r, E0_r, E_total, T0_m, E0_m, Ek, Tk, sigma_a, sigma_s, sigma_t, abs, emis, rho, cv, Pr, Er_in, a, c, dx, xf, p=2); //p=0 mach 1.2, p=1 mach 3, p=2 mach 45, p=3 sod shock tube, p=4 marshak
    mat_initialize(v0, rho0, Pm, e0, v, e, rho, as, Pr, P, T0_m,  v_in, e_in, rho_in, P_in, dx, mat_gamma, xf, p=2); //p=0 mach 1.2, p=1 mach 3, p=2 mach 45, p=3 sod shock tube, p=4 marshak

    //!For Moving Shock
    //moving_shock_rad_init(T0_r, E0_r, E_total, T0_m, E0_m, Ek, Tk, sigma_a, sigma_s, sigma_t, abs, emis, rho, cv, Pr, Er_in, a, c, dx, xf);
    //moving_shock_mat_init(v0, rho0, Pm, e0, v, e, rho, as, Pr, P, T0_m, cv, v_in, e_in, rho_in, P_in, S, dx, mat_gamma, xf, dt);
    for(int i=0; i<n_time; i++){

        //Hydro Step
        reassign(Pm, v, e, rho, v0, rho0,  e0, as, mat_gamma);
        
        //Predictor Step
        flux(v0, rho0, P, e0, Er, v_in, e_in, rho_in, P_in, Er_in, as, lu_rho, lu_v, lu_e, lu_Er, grad_u);

        eularian_calcs(v0, rho0, e0, Pr, lu_rho, lu_v, lu_e, grad_u, cv, Pm_pre, v_pre, e_pre, rho_pre, P_pre, Th, as, dt_half, dx, mat_gamma, xf);

        //Corrector Step
        flux(v_pre, rho_pre, P_pre, e_pre, Er, v_in, e_in, rho_in, P_in, Er_in, as, lu_rho, lu_v, lu_e, lu_Er, grad_u);

        eularian_calcs(v0, rho0, e0, Pr, lu_rho, lu_v, lu_e, grad_u, cv, Pm, v, e, rho, P, Th, as, dt, dx, mat_gamma, xf);
        
        //Radiation MMC step
        rad_mmc(E0_r, Pr, grad_u, lu_Er, Es_r, dt, dx);

        //Radiation Setup
        setup(sigma_a, sigma_s, sigma_t, T0_m, D, Dp, Dm, cv, rho, c, a, p=0); //p=0 sod and fixed shock, p=2 marshak

        //Iteration Loop for implicit temp and energy
        for(int k=0; k<n_iter; k++){

            //Newtons method to solve for temp
            newtons(Ek, Th, sigma_a, rho, cv, Tk, dt, c, a, n_iter);

            //set up matrix
            matrix(Tk, Es_r, sigma_a, Dp, Dm, plus, mid, minus, rs, dt, dx, c, a, p=2); //p=0 zero boundary, p=1 marshak, p=2 reflective, p=3 inflow/reflective

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

            if(delta_max < 1.0E-5){
                break;
            }
        }

        //Energy Deposition Step
        e_dep(cv, Th, Tk, e);

        //Reassign rad values
        reassign(Tk, Ek, rho, cv, E0_r, T0_r, E0_m, T0_m, Pr, Pm, P, a, c, p=0); //p=2 sod shock, p=0 else

        double time = dt * i;

        std::cout << " time " << time << std::endl;
    }
    
    
    myfile.open("../results/mach_45_diff_shock_70_ns.dat");
    for(int i=0; i<n_cells; i++){
        myfile << i*dx << " " << T0_r[i] << " " << T0_m[i] << " " << rho[i] << " " << v[i] << " " << e[i] << " " << P[i] <<  "\n";
    } 
    myfile.close();
    

    /*
    myfile.open("../results/FM_moving_shock_15000.dat");
    for(int i=0; i<n_cells; i++){
        myfile << i*dx << " " << T0_r[i] << " " << T0_m[i] << " " << "\n";
    } 
    myfile.close();
    */
}