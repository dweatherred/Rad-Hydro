
timing = {}
timing["dt"] = 1.0E-13
timing["time_initial"] = 0
timing["time_final"] = 1.0E-7;

spatial_cells = {}
spatial_cells["xi"] = 0 --cm
spatial_cells["xf"] = 0.2 --cm
spatial_cells["num_cells"] = 15000

--[[
Sod Shock
timing = {}
timing["dt"] = 9.95586E-06
timing["time_initial"] = 0
timing["time_final"] = 0.2 --seconds

spatial_cells = {}
spatial_cells["xi"] = 0 --cm
spatial_cells["xf"] = 1 --cm
spatial_cells["num_cells"] = 5000
]]--

--[[
Radiative
timing = {}
timing["dt"] = 1.0E-14 --.84062e-05
timing["time_initial"] = 0
timing["time_final"] = 5.0E-9

spatial_cells = {}
spatial_cells["xi"] = 0 --cm
spatial_cells["xf"] = 0.4--cm
spatial_cells["num_cells"] = 15000
]]--