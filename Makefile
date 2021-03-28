build_cpu = x86_64
header_files = include/FronTier/util/plotdecs.h include/FronTier/util/cdecs.h include/FronTier/util/uprotos.h include/FronTier/util/navdecs.h include/FronTier/util/vmalloc.h include/FronTier/util/fnamedebug.h include/FronTier/intfc/ilocprotos.h include/FronTier/intfc/iprotos.h include/FronTier/intfc/iMeshP.h include/FronTier/intfc/int.h include/FronTier/intfc/iGeom.h include/FronTier/intfc/iBase.h include/FronTier/intfc/geom.h include/FronTier/intfc/iTaps.h include/FronTier/intfc/userint.h include/FronTier/intfc/iMesh.h include/FronTier/intfc/triangledefs.h include/FronTier/intfc/int_amr.h include/FronTier/intfc/iRel.h include/FronTier/intfc/iloc.h include/FronTier/intfc/table.h include/FronTier/intfc/array.h include/FronTier/front/fapi.h include/FronTier/front/fprotos.h include/FronTier/front/fpatrecon.h include/FronTier/front/fdecs.h include/FronTier/front/fuserint.h include/FronTier/front/ftapi.h include/FronTier/front/frp.h include/FronTier/front/fvelo.h 
#include_dirs=include include/FronTier include/FronTier/front include/FronTier/util include/FronTier/intfc
include_dirs=include/FronTier/intfc include/FronTier/front include/FronTier/util
lib_only:
	cd src && $(MAKE)
	$(MAKE) incs 

all: 
	cd src && $(MAKE)
	$(MAKE) incs 
	cd airfoil && $(MAKE)
	cd cell && $(MAKE)
	cd cFluid && $(MAKE)
	cd cim && $(MAKE)
	cd inverse && $(MAKE)
	cd crystal && $(MAKE)
	cd curvature && $(MAKE)
	cd finance && $(MAKE)
	cd frgb && $(MAKE)
	cd iFluid && $(MAKE)
	cd melting && $(MAKE)
	cd pde && $(MAKE)
	cd MonteCarlo && $(MAKE)
	cd poisson && $(MAKE)
	cd subsurf && $(MAKE)
clean:
	cd src && $(MAKE) clean
	cd example && $(MAKE) clean
	cd example3d && $(MAKE) clean
	cd iFluid && $(MAKE) clean
	cd cFluid && $(MAKE) clean
	cd cim && $(MAKE) clean
	cd inverse && $(MAKE) clean
	cd airfoil && $(MAKE) clean
	cd cell && $(MAKE) clean
	cd crystal && $(MAKE) clean
	cd curvature && $(MAKE) clean
	cd finance && $(MAKE) clean
	cd frgb && $(MAKE) clean
	cd melting && $(MAKE) clean
	cd pde && $(MAKE) clean
	cd MonteCarlo && $(MAKE) clean
	cd poisson && $(MAKE) clean
	cd subsurf && $(MAKE) clean
	cd iTaps && $(MAKE) clean
	-cd lib/$(build_cpu)/; rm -rf gas *.a
	-rm -rf incs
tar:
	cd src/gas && $(MAKE) tar	
export:
	cd src/gas && $(MAKE) export

diff:
	cd src && $(MAKE) diff
	-diff -r . $(diffdir) > DIFF

incs: $(include_dirs) include/FronTier.h $(header_files)

include/FronTier:
	mkdir include include/FronTier
$(include_dirs): include/FronTier
	mkdir $@

include/FronTier.h: $(include_dirs) 
	@echo "#include <FronTier/front/fdecs.h>" > include/FronTier.h
$(header_files): include/FronTier/%: src/% 
	@echo "updating: "$@; \
	sed -e "1,$$ s,include <,include <FronTier/,g" $^ | \
	sed -e "1,$$ s,include <FronTier/hdf,include <hdf,g" | \
	sed -e "1,$$ s,include <FronTier/mpi,include <mpi,g" | \
	sed -e "1,$$ s,include <FronTier/mfhdf,include <mfhdf,g" | \
	sed -e "1,$$ s,FronTier/cdecs.h,FronTier/util/cdecs.h,g" | \
	sed -e "1,$$ s,FronTier/vmalloc.h,FronTier/util/vmalloc.h,g" | \
	sed -e "1,$$ s,FronTier/fnamedebug.h,FronTier/util/fnamedebug.h,g"|  \
	sed -e "1,$$ s,FronTier/uprotos.h,FronTier/util/uprotos.h,g"| \
	sed -e "1,$$ s,FronTier/navdecs.h,FronTier/util/navdecs.h,g" | \
	sed -e "1,$$ s,FronTier/gd.h,gd.h,g" | \
	sed -e "1,$$ s,FronTier/gdfonts.h,gdfonts.h,g" | \
	sed -e "1,$$ s,FronTier/gdfontl.h,gdfontl.h,g" | \
	sed -e "1,$$ s,FronTier/gdfontt.h,gdfontt.h,g" | \
	sed -e "1,$$ s,FronTier/stdlib.h,stdlib.h,g" | \
	sed -e "1,$$ s,FronTier/stdio.h,stdio.h,g" | \
	sed -e "1,$$ s,FronTier/stdint.h,stdint.h,g" | \
	sed -e "1,$$ s,FronTier/string.h,string.h,g" | \
	sed -e "1,$$ s,FronTier/unistd.h,unistd.h,g" | \
	sed -e "1,$$ s,FronTier/ctype.h,ctype.h,g" | \
	sed -e "1,$$ s,FronTier/math.h,math.h,g" | \
	sed -e "1,$$ s,FronTier/limits.h,limits.h,g" | \
	sed -e "1,$$ s,FronTier/float.h,float.h,g" | \
	sed -e "1,$$ s,FronTier/errno.h,errno.h,g" | \
	sed -e "1,$$ s,FronTier/libgen.h,libgen.h,g" | \
	sed -e "1,$$ s,FronTier/algorithm,algorithm,g" | \
	sed -e "1,$$ s,FronTier/assert,assert,g" | \
	sed -e "1,$$ s,FronTier/string,string,g" > $@
$(header_files): $(include_dirs)









test:	
	cd example/; $(MAKE)

TARFLIST = $$ME/README \
$$ME/build.sh \
$$ME/advbuild.sh \
$$ME/install-sh \
$$ME/config.status \
$$ME/config.sub \
$$ME/config.guess \
$$ME/devel-deps.inc.in \
$$ME/iMesh-Defs.inc.in \
$$ME/configure.in \
$$ME/src/util/mkdep.pl \
$$ME/src/util/mkfiles.pl \
$$ME/src/util/ftrules \
$$ME/src/util/gasrules \
$$ME/src/pak/Makefile.in \
$$ME/src/pak/odepack/Makefile.in \
$$ME/src/pak/odepack/*.f \
$$ME/src/pak/odepack/*.c \
$$ME/src/pak/dierckx/Makefile.in \
$$ME/src/pak/dierckx/*.f \
$$ME/src/pak/linpak/Makefile.in \
$$ME/src/pak/linpak/*.f \
$$ME/src/pak/blas/Makefile.in \
$$ME/src/pak/blas/*.f \
$$ME/configure.in \
$$ME/Makefile.in \
$$ME/src/Makefile.in \
$$ME/example/*.cpp \
$$ME/example/Makefile.in \
$$ME/example/README \
$$ME/example3d/*.cpp \
$$ME/example3d/Makefile.in \
$$ME/curvature/*.cpp \
$$ME/curvature/*.h \
$$ME/curvature/Makefile.in \
$$ME/example3d/README \
$$ME/subsurf/*.cpp \
$$ME/subsurf/*.h \
$$ME/subsurf/Makefile.in \
$$ME/subsurf/in-* \
$$ME/airfoil/*.cpp \
$$ME/airfoil/*.h \
$$ME/airfoil/Makefile.in \
$$ME/airfoil/in-* \
$$ME/airfoil/mrun* \
$$ME/crystal/*.cpp \
$$ME/crystal/*.h \
$$ME/crystal/Makefile.in \
$$ME/crystal/README \
$$ME/crystal/in-* \
$$ME/crystal/mrun* \
$$ME/melting/*.cpp \
$$ME/melting/*.h \
$$ME/melting/Makefile.in \
$$ME/melting/in-* \
$$ME/melting/mrun* \
$$ME/iFluid/*.cpp \
$$ME/iFluid/*.h \
$$ME/iFluid/Makefile.in \
$$ME/iFluid/in-* \
$$ME/iFluid/mrun* \
$$ME/cFluid/*.cpp \
$$ME/cFluid/*.h \
$$ME/cFluid/Makefile.in \
$$ME/cFluid/in-* \
$$ME/cFluid/mrun* \
$$ME/cim/*.cpp \
$$ME/cim/*.h \
$$ME/cim/Makefile.in \
$$ME/cim/in-* \
$$ME/cim/mrun* \
$$ME/inverse/*.cpp \
$$ME/inverse/*.h \
$$ME/inverse/Makefile.in \
$$ME/inverse/in-* \
$$ME/inverse/mrun* \
$$ME/solver/*.cpp \
$$ME/solver/*.h \
$$ME/solver/Makefile.in \
$$ME/weno/*.cpp \
$$ME/weno/*.h \
$$ME/weno/Makefile.in \
$$ME/weno/in-* \
$$ME/frgb/*.cpp \
$$ME/frgb/*.h \
$$ME/frgb/Makefile.in \
$$ME/frgb/in-* \
$$ME/frgb/mrun* \
$$ME/README \
$$ME/finance/*.cpp \
$$ME/finance/*.h \
$$ME/finance/Makefile.in \
$$ME/finance/in-* \
$$ME/finance/mrun* \
$$ME/cell/*.cpp \
$$ME/cell/*.h \
$$ME/cell/Makefile.in \
$$ME/pde/*.cpp \
$$ME/pde/in-* \
$$ME/pde/Makefile.in \
$$ME/MonteCarlo/*.h \
$$ME/MonteCarlo/*.cpp \
$$ME/MonteCarlo/in-* \
$$ME/MonteCarlo/Makefile.in \
$$ME/poisson/*.cpp \
$$ME/poisson/*.h \
$$ME/poisson/Makefile.in \
$$ME/iTaps/*.cpp \
$$ME/iTaps/Makefile.in \
$$ME/src/intfc/*.c \
$$ME/src/intfc/*.h  \
$$ME/src/intfc/Makefile.in \
$$ME/src/util/*.h \
$$ME/src/util/*.c \
$$ME/src/util/Makefile.in \
$$ME/src/front/*.c \
$$ME/src/front/*.h \
$$ME/src/front/Makefile.in \
$$ME/src/front/README \
$$ME/src/pak/hypre-1.6.0.tar.gz \
$$ME/src/pak/sn_ellip/*.C \
$$ME/src/pak/sn_ellip/*.h \
$$ME/src/pak/sn_ellip/Makefile.in \

tar:
	cd ./ && ME=`basename $$PWD` && cd .. && tar -cf "`date +FronTier.Lite.%m_%d_%y.tar`" ${TARFLIST} && \
	gzip "`date +FronTier.Lite.%m_%d_%y.tar`" && mv `date +FronTier.Lite.%m_%d_%y.tar`.gz $$ME/.
