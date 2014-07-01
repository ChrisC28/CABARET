ifort -c NetCDFModule.f90 -I/home/chris/Downloads/netcdf-4.0.1/include
ifort -c testCABARET2D.f90  -I/home/chris/Downloads/netcdf-4.0.1/include

ifort -o cabaret2d.exe testCABARET2D.o NetCDFModule.o -L/home/chris/Downloads/netcdf-4.0.1/lib/  -lnetcdf -lcurl -lz -lblas -llapack
