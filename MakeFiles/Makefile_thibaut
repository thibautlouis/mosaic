python/mcm_code_full.so: python/mcm_code_full.f90 python/wigner3j_sub.f
	cd python && F90FLAGS="-fopenmp -fPIC -Ofast -ffree-line-length-none" f2py-2.7 -c -m mcm_code_full mcm_code_full.f90 wigner3j_sub.f -lgomp
clean:
	rm -f python/*.so
	rm -rf python/*.so.dSYM
	rm -rf  python/*.pyc
