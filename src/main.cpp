//
// Main function file for the Eularian Radiation-Hydrodynamics Problem. This file reads the output from Lua_Script.lua
// and passes them into the appropriate functions.
//

#include <stdio.h>
#include <iostream>
#include <fstream>

//Include lua hpp file
#include "../Lua542/include/lua.hpp"

//Different Rad-Hydro Function Files
#include "eularian_rh_so.hpp"

double linear_interpolate(const double xp, std::vector<double> &x_tab, std::vector<double> &y_tab){

    //First binary search data to see if point is in table
    double first = 0.0;
    double last = x_tab.size();
    double middle = round((first + last)/2.0);

    while(first <= last){
        if(x_tab[middle] < xp){
            first = middle + 1.0;
        }else if(x_tab[middle] > xp){
            last = middle - 1.0;
        }else{
            return y_tab[middle];
        }
        middle = round((first + last)/2.0);
    }

    //If point does not exist in the table, linear interpolate
    if(first > last){
        return y_tab[last] + ((y_tab[first] - y_tab[last])/(x_tab[first] - x_tab[last])) * (xp - x_tab[last]);
    }
}

//This is an error check for lua
bool CheckLua(lua_State *L, int r){

	if (r != LUA_OK){
		std::string errormsg = lua_tostring(L, -1);
		std::cout << errormsg << std::endl; 
		return false;
	}
	return true;
}

int main (int argc, char *argv[]){

    //Timing Variables
    double dt, ti, tf, n_time;

    //Spacial Variables 
    double xi, xf, n_cells, dx;
    
    //Opacity Values
    double alpha_a, alpha_s, nTa, nTs, nrho_a, nrho_s;

    //Material Constants
    double gamma, cv;

    //Boundary Temperatures
    double Tbc_l, Tbc_r;

    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    for(int i=1; i<argc; i++){
        
        if(std::string(argv[1]) == "Mach_45"){
            if (CheckLua(L, luaL_dofile(L, "../inputs/Mach_45/Mach_45.lua"))){
                lua_getglobal(L, "timing");
                if (lua_istable(L, -1)){
                    lua_pushstring(L, "dt");
                    lua_gettable(L, -2);
                    dt = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                    lua_pushstring(L, "time_initial");
                    lua_gettable(L, -2);
                    ti = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                    lua_pushstring(L, "time_final");
                    lua_gettable(L, -2);
                    tf = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                    n_time = (ti + tf) / dt;
                }
                lua_getglobal(L, "spatial_cells");
                if (lua_istable(L, -1)){
                    lua_pushstring(L, "xi");
                    lua_gettable(L, -2);
                    xi = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                    lua_pushstring(L, "xf");
                    lua_gettable(L, -2);
                    xf = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                    lua_pushstring(L, "num_cells");
                    lua_gettable(L, -2);
                    n_cells = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                    dx = (xf- xi)/n_cells;
                }
                lua_getglobal(L, "opacity");
                if (lua_istable(L, -1)){
                    lua_pushstring(L, "alpha_a");
                    lua_gettable(L, -2);
                    alpha_a = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                    lua_pushstring(L, "alpha_s");
                    lua_gettable(L, -2);
                    alpha_s = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                    lua_pushstring(L, "nrho_s");
                    lua_gettable(L, -2);
                    nrho_s = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                    lua_pushstring(L, "nrho_a");
                    lua_gettable(L, -2);
                    nrho_s = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                    lua_pushstring(L, "nT_a");
                    lua_gettable(L, -2);
                    nrho_s = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                    lua_pushstring(L, "nT_s");
                    lua_gettable(L, -2);
                    nrho_s = (double)lua_tonumber(L, -1);
                    lua_pop(L, 1);
                }
                //Get constants
                lua_getglobal(L, "gamma");
                gamma = (double)lua_tonumber(L, -1);
                lua_getglobal(L, "cv");
                cv = (double)lua_tonumber(L, -1);
                
                //Boundary Temperatures
                lua_getglobal(L, "Tbc_l");
                Tbc_l = (double)lua_tonumber(L, -1);
                lua_getglobal(L, "Tbc_r");
                Tbc_r = (double)lua_tonumber(L, -1);
            }
        }
    }

    if (CheckLua(L, luaL_dofile(L, "lua_radhydro.lua"))){

        lua_getglobal(L, "timing");
        if (lua_istable(L, -1)){
            lua_pushstring(L, "dt");
            lua_gettable(L, -2);
            dt = (double)lua_tonumber(L, -1);
            lua_pop(L, 1);
            lua_pushstring(L, "time_initial");
            lua_gettable(L, -2);
            ti = (double)lua_tonumber(L, -1);
            lua_pop(L, 1);
            lua_pushstring(L, "time_final");
            lua_gettable(L, -2);
            tf = (double)lua_tonumber(L, -1);
            lua_pop(L, 1);
        }

        lua_getglobal(L, "spatial_cells");
        if (lua_istable(L, -1)){

            lua_pushstring(L, "xi");
            lua_gettable(L, -2);
            xi = (double)lua_tonumber(L, -1);
            lua_pop(L, 1);

            lua_pushstring(L, "xf");
            lua_gettable(L, -2);
            xf = (double)lua_tonumber(L, -1);
            lua_pop(L, 1);

            lua_pushstring(L, "num_cells");
            lua_gettable(L, -2);
            n_cells = (double)lua_tonumber(L, -1);
            lua_pop(L, 1);

            dx = (xf- xi)/n_cells;

	    }
    }

    //call main functions for different rad-hydro methods
    so_eularian_rh(dt, ti, tf, dx, xf, n_cells);


    return 0;
}