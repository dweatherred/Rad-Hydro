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

//This is an error check for lua
bool CheckLua(lua_State *L, int r){

	if (r != LUA_OK){
		std::string errormsg = lua_tostring(L, -1);
		std::cout << errormsg << std::endl; 
		return false;
	}
	return true;
}

int main (){

    //Function Parameters
    double dt, ti, tf, xi, xf, dx;
    int n_cells;

    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

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
            n_cells = (int)lua_tonumber(L, -1);
            lua_pop(L, 1);

            dx = (xf- xi)/n_cells;

	    }
    }

    //call main functions for different rad-hydro methods
    so_eularian_rh(dt, ti, tf, dx, xf, n_cells);


    return 0;
}