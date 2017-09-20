#
# Elemental: A framework for distributed memory dense linear algebra.
#
# Copyright 2009-2010 Jack Poulson
#

srcdir = src
incdir = include
testdir = test
libdir = lib
bindir = bin

library = libelemental.a

#WITHOUT_COMPLEX = 
# Compile flags:
#   WITHOUT_COMPLEX: if defined, no complex datatypes are implemented
#   FUNDERSCORE: if defined, all BLAS/LAPACK wrappers assume underscores
#   POOL_MEMORY: if defined, Memory class only accumulates until destruction
#   RELEASE: if defined, callstack is not maintained and debug checks are off
CXX = mpicxx
#CXXFLAGS = -DFUNDERSCORE -I$(incdir)
CXXFLAGS = -DWITHOUT_COMPLEX -DFUNDERSCORE -I$(incdir)
CXXFLAGS_DEBUG = -g -Wall $(CXXFLAGS)
CXXFLAGS_RELEASE = -O3 -Wall -DRELEASE $(CXXFLAGS)
#LDFLAGS = -llapack -lblas -L/usr/lib
LDFLAGS = -mkl -L/usr/lib
AR = ar
ARFLAGS = rc

################################################################################
# Only developers should edit past this point.                                 #
################################################################################

# Source/object organization
coredir = core
corefiles = Environment.cpp \
            Grid.cpp \
            Matrix.cpp \
            DistMatrix/MC_MR.cpp \
            DistMatrix/MC_Star.cpp \
            DistMatrix/MR_MC.cpp \
            DistMatrix/MR_Star.cpp \
            DistMatrix/Star_MC.cpp \
            DistMatrix/Star_MR.cpp \
            DistMatrix/Star_Star.cpp \
            DistMatrix/Star_VC.cpp \
            DistMatrix/Star_VR.cpp \
            DistMatrix/VC_Star.cpp \
            DistMatrix/VR_Star.cpp \

coresrc = $(addprefix $(coredir)/,$(corefiles))

blasdir = BLAS
blasfiles = Level2/Gemv/Gemv.cpp \
            Level2/Gemv/GemvN.cpp \
            Level2/Ger/Ger.cpp \
            Level3/Gemm/Gemm.cpp \
            Level3/Gemm/GemmNN.cpp \

blassrc = $(addprefix $(blasdir)/,$(blasfiles))

lapackdir = LAPACK
lapackfiles = #Chol/Chol.cpp \
              #Chol/CholL.cpp \
              #Chol/CholU.cpp \
              
lapacksrc = $(addprefix $(lapackdir)/,$(lapackfiles))

# The entire list of source files relative to $(srcdir)
src = $(coresrc) $(blassrc) $(lapacksrc)

includefiles = Elemental.h \
               ElementalBLAS.h \
               ElementalBLAS_Internal.h \
               ElementalDistMatrix.h \
               ElementalDistMatrix_MC_MR.h \
               ElementalDistMatrix_MC_Star.h \
               ElementalDistMatrix_MR_MC.h \
               ElementalDistMatrix_MR_Star.h \
               ElementalDistMatrix_Star_MC.h \
               ElementalDistMatrix_Star_MR.h \
               ElementalDistMatrix_Star_Star.h \
               ElementalDistMatrix_Star_VC.h \
               ElementalDistMatrix_Star_VR.h \
               ElementalDistMatrix_VC_Star.h \
               ElementalDistMatrix_VR_Star.h \
               ElementalEnvironment.h \
               ElementalGrid.h \
               ElementalMemory.h \
               ElementalMatrix.h \
               ElementalPartitioning.h \
               ElementalRandom.h \
               ElementalTypes.h \
               wrappers/BLAS.h \
               wrappers/MPI.h \

includes = $(addprefix $(incdir)/,$(includefiles)) \
		   $(srcdir)/$(coredir)/DistMatrix/DistMatrixMacros.h

################################################################################
# make                                                                         #
################################################################################
libdir_debug = $(libdir)/debug
libdir_release = $(libdir)/release
obj_debug = $(addprefix $(libdir_debug)/,$(src:.cpp=.o))
obj_release = $(addprefix $(libdir_release)/,$(src:.cpp=.o))
library_debug = $(libdir_debug)/$(library)
library_release = $(libdir_release)/$(library)

# This is the default target
.PHONY : lib
lib: release debug

.PHONY : debug
debug: $(library_debug) 

.PHONY : release
release: $(library_release)

$(library_debug): $(obj_debug)
	@echo "[ debug ] Creating $@"
	@$(AR) $(ARFLAGS) $@ $^

$(library_release): $(obj_release)
	@echo "[release] Creating $@"
	@$(AR) $(ARFLAGS) $@ $^

# Object files must depend upon headers because we inline functions
$(libdir_debug)/%.o: $(srcdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[ debug ] Compiling $<"
	@$(CXX) $(CXXFLAGS_DEBUG) -c -o $@ $<

$(libdir_release)/%.o: $(srcdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[release] Compiling $<"
	@$(CXX) $(CXXFLAGS_RELEASE) -c -o $@ $<

################################################################################
# make test                                                                    #
################################################################################
bindir_debug = $(bindir)/debug
bindir_release = $(bindir)/release

tests = DistMatrix \
        BLAS/Gemm \
        BLAS/Gemv \
        BLAS/Ger \
        #BLAS/Hemm \
        #BLAS/Her2k \
        #BLAS/Herk \
        #BLAS/Symm \
        #BLAS/Symv \
        #BLAS/Syr2k \
        #BLAS/Syrk \
        #BLAS/Trmm \
        #BLAS/Trsm \
        #LAPACK/Chol \
        #LAPACK/LU \
        #LAPACK/Tridiag \
        #LAPACK/Trinv 

testobjs = $(addsuffix .o, $(tests))

tests_debug = $(addprefix $(bindir_debug)/, $(tests))
testobjs_debug = $(addprefix $(bindir_debug)/, $(testobjs))
tests_release = $(addprefix $(bindir_release)/, $(tests))

.PHONY : test
test: test-release test-debug

.PHONY : test-debug
test-debug: $(tests_debug) $(testobjs_debug)

.PHONY : test-release
test-release: $(tests_release)

$(bindir_debug)/%: $(bindir_debug)/%.o $(library_debug)
	@echo "[ debug ] Creating $@"
	@$(CXX) -o $@ $^ $(LDFLAGS)

$(bindir_release)/%: $(bindir_release)/%.o $(library_release)
	@echo "[release] Creating $@"
	@$(CXX) -o $@ $^ $(LDFLAGS)

$(bindir_debug)/%.o: $(testdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[ debug ] Compiling $<"
	@$(CXX) $(CXXFLAGS_DEBUG) -c -o $@ $<

$(bindir_release)/%.o: $(testdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[release] Compiling $<"
	@$(CXX) $(CXXFLAGS_RELEASE) -c -o $@ $<

################################################################################
# make clean                                                                   #
################################################################################
.PHONY : clean
clean: 
	@rm -Rf lib/
	@rm -Rf bin/

.PHONY : clean-debug
clean-debug:
	@rm -Rf lib/debug
	@rm -Rf bin/debug

.PHONY : clean-release
clean-release:
	@rm -Rf lib/release
	@rm -Rf bin/release

