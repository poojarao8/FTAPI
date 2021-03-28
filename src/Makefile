CC=/nfs/scratch/jmelvin/lib/openmpi/bin/mpicc
F77=/nfs/scratch/jmelvin/lib/openmpi/bin/mpif77
all:
	if [ 0 -eq 1 ]; then cd pak/cgnslib_3.0 && $(MAKE); fi
	# if [ 0 -eq 1 ]; then cd pak/pcgnslib_0.2.0 && $(MAKE); fi
	cd front  && $(MAKE)
	  
	if [ 0 -eq 1 ]; then cd pak/linpak && $(MAKE); fi
	if [ 0 -eq 1 ]; then cd pak/blas && $(MAKE); fi

clean:
	cd util && $(MAKE) clean
	cd front && $(MAKE) clean
	cd intfc && $(MAKE) clean
	 
	if [ 0 -eq 1 ]; then cd pak/linpak && $(MAKE) clean; fi
	if [ 0 -eq 1 ]; then cd pak/blas && $(MAKE) clean; fi
	if [ 0 -eq 1 ]; then cd pak/cgnslib_3.0 && $(MAKE) clean; fi
	# if [ 0 -eq 1 ]; then cd pak/pcgnslib_0.2.0 && $(MAKE) clean; fi

diff:
	cd driver && $(MAKE) diff
	cd util && $(MAKE) diff 
	cd front && $(MAKE) diff
	cd hyp && $(MAKE) diff
	cd intfc && $(MAKE) diff
	cd tri  &&  $(MAKE) diff
	cd gas  &&  $(MAKE) diff
