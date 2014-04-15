
f90 genesis.f90 -o genesis -I/opt/local/netcdf/include -L/opt/local/netcdf/lib -lnetcdf ngtopt.o global.o ncread_era.o ncread_jra.o ncread_ncep.o ncread_um.o establ.o find_latlon.o interp_linear.o search_grid.o geostr.o surfdist.o radian.o charney_phillips.o create_namelist.o
