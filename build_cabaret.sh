ifort -c NetCDFModule.f90 -I/home/chris/Downloads/netcdf-4.0.1/include
ifort -c testCABARET.f90  -I/home/chris/Downloads/netcdf-4.0.1/include

ifort -o out.exe testCABARET.o NetCDFModule.o -L/home/chris/Downloads/netcdf-4.0.1/lib/  -lnetcdf -lcurl -lz -lblas -llapack
