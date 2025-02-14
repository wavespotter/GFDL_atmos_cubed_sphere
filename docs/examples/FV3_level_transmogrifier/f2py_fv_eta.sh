rm -rf *.o *.mod *.so
f2py -c -m fv_eta fv_eta.F90 -I`pwd`
