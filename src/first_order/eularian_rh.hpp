#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

//include main function file
#include "eularian_rh_func.hpp"


void eularian_rh(double dt, double ti, double tf, double dx, double xf, int n_cells){

    //Needed Constants
    const int n_time = tf/dt;

    //!Hydro Variables

    //constants
    double gamma = 1.66666667;

    //Conserved values
    std::vector<double> m(n_cells,0); //mass in cell
    std::vector<double> v0(n_cells,0); //initial velocity in cell
    std::vector<double> rho0(n_cells,0); //initial density in cell
    std::vector<double> e0(n_cells,0); //inital energy in cell
    std::vector<double> as(n_cells,0); //speed of sound
    std::vector<double> v(n_cells,0); //velocity in cell
    std::vector<double> rho(n_cells,0); // density in cell
    std::vector<double> P(n_cells,0); //total pressure in cell
    std::vector<double> Pm(n_cells,0); //Material Pressure in cell
    std::vector<double> e(n_cells,0); // total energy in cell

    //Flux terms
    std::vector<double> lu_rho(n_cells,0);
    std::vector<double> lu_v(n_cells,0);
    std::vector<double> lu_e(n_cells,0);
    std::vector<double> lu_w(n_cells,0);

    //!Rad Variables

    //Constants
    const int a = 137;
    double c = 3E+10; // cm/s
    const int n_iter = 10000000;

    //Temperature
    std::vector<double> T0_r(n_cells);
    std::vector<double> T0_m(n_cells);
    
    //Energy
    std::vector<double> Es_r(n_cells);
    std::vector<double> E0_r(n_cells);
    std::vector<double> E0_m(n_cells);
    double E_total;

    //Radiation Pressure
    std::vector<double> Pr(n_cells);

    //Absorption and Emission
    std::vector<double> abs(n_cells);
    std::vector<double> emis(n_cells);

    //Non-constant Opacity
    std::vector<double> opc(n_cells);

    //cell specific density and opacity
    std::vector<double> cv(n_cells);

    //Diffusion Term
    std::vector<double> D(n_cells);
    std::vector<double> Dp(n_cells);
    std::vector<double> Dm(n_cells);

    //Newton Iteration Guesses
    std::vector<double> Ek(n_cells);
    std::vector<double> Tk(n_cells);

    //Flux Values
    std::vector<double> lu_Er(n_cells);
    std::vector<double> grad_u(n_cells);

    //matrix for thomas algorithm
    std::vector<double> plus(n_cells);
    std::vector<double> mid(n_cells);
    std::vector<double> minus(n_cells); 
    std::vector<double> rs(n_cells);
    std::vector<double> Er(n_cells);

    //residual values
    std::vector<double> res0(n_cells,0);
    std::vector<double> res1(n_cells,0);    

    //constant to switch zero to marshak boundary
    int p = 1;

    //!Initialize Values
    rad_initialize(T0_r, E0_r, E_total, T0_m, E0_m, Ek, Tk, opc, abs, emis, rho, cv, Pr, a, c, dx, xf, gamma);
    mat_initialize(m,v0, rho0, Pm, e0, v, e, rho, as, Pr, P, T0_m, dx, gamma, xf);

    for(int i=0; i<n_time; i++){

        //Hydro Step
        eularian_reassign(m, v0, rho0, Pm, e0, v, e, rho, as, gamma, dx);
        
        //Calculte flux values prior to hydro and rad calcs
        calc_flux(m, v0, rho0, P, e0, as, E0_r, lu_rho, lu_v, lu_e, lu_Er, grad_u);
        
        eularian_calcs(v0, rho0, P, e0, v, e, rho, lu_rho, lu_v, lu_e, Pr, grad_u, dt, dx, gamma);
        
        //Radiation MMC step
        rad_mmc(Es_r, E0_r, Pr, grad_u, lu_Er, dt, dx);

        //Radiation Solve
        setup(opc, T0_m, D, Dp, Dm, cv, c, a);

        //Iteration Loop for implicit temp and energy
        for(int k=0; k<n_iter; k++){

            //Newtons method to solve for temp
            newtons(Ek, Tk, E0_r, T0_m, opc, dt, c, a, rho, cv, n_iter, p);

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
        e_dep(e, cv, T0_m, Tk);

        //Reassign rad values
        reassign(E0_r, T0_r, E0_m, T0_m, Tk, Ek, abs, emis, opc, rho, cv, Pr, Pm, P, a, c, gamma);

        //std::cout.precision(17);
        //for(int j=0; j<n_cells; j++){
        //    int i = j;
        //    std::cout << " i " << i <<  " E0_r " << E0_r[i] << " E0_m " << E0_m[i] <<  " T0_r " << T0_r[i] <<  " T0_m " << T0_m[i] << " rho " << rho[i] 
        //    << " v " << v[i] << " P " << P[i] <<  " e " << e[i] << std::endl;
       //}

        double time = dt * i;
        std::cout << " time " << time << std::endl;

        //std::cin.get();

    }

    std::ofstream myfile;

    /*
    myfile.open("../results/rad_hydro_shock_tube.dat");
    for(int i=0; i<n_cells; i++){
         myfile << i*dx << " " << E0_r[i] << " " << rho[i] << " " << v[i] << " " <<P[i] << "\n";
    }
    myfile.close();
    */
   
    /*
    myfile.open("../results/rad_hydro_marshak.dat");
    for(int i=0; i<n_cells; i++){
        myfile << i*dx << " " << T0_r[i] << " " << T0_m[i] << "\n";
        
    } 
    myfile.close();
    */

   /*
    myfile.open("../results/rad_hydro_suolson.dat");
    for(int i=0; i<n_cells; i++){
        myfile << i*dx << " " << T0_r[i] << " " << T0_m[i] << "\n";
        
    } 
    myfile.close();
    */

    
    myfile.open("../results/rad_hydro_shock_trying.dat");
    for(int i=0; i<n_cells; i++){
        myfile << i*dx << " " << T0_r[i] << " " << T0_m[i] << " " << rho[i] << " " << v[i] << " " << e[i] << " " << P[i] <<  "\n";
        
    } 
    myfile.close();
    
    
}
