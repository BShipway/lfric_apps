##############################################################################
# (c) Crown copyright 2017-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

PROJECT_NAME = gungho_model
export PROFILE ?= fast-debug

# Select option for PSyclone transformation. Select "none" for no transformation.
export PSYCLONE_TRANSFORMATION ?= minimum

# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

export EXTRA_ROSE_META = jules
export META_DIR = $(WORKING_DIR)/../rose-meta
export ROSE_META_DIRS = $(APPS_ROOT_DIR)/rose-meta $(META_DIR)

export INTERNAL_DEPENDENCIES = $(CORE_ROOT_DIR)/infrastructure \
                               $(CORE_ROOT_DIR)/components/driver \
                               $(CORE_ROOT_DIR)/components/science \
                               $(CORE_ROOT_DIR)/components/inventory \
                               $(CORE_ROOT_DIR)/components/lfric-xios \
                               $(APPS_ROOT_DIR)/science/gungho

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

# Check PSyclone transformation option is valid
ifneq ("$(PSYCLONE_TRANSFORMATION)","none")
  export OPTIMISATION_PATH ?= $(wildcard optimisation/$(PSYCLONE_TRANSFORMATION))
  ifeq ("$(wildcard $(OPTIMISATION_PATH)/psykal/global.py)","")
    $(error PSyclone transformation dir "$(PSYCLONE_TRANSFORMATION)" does not have a global.py script)
  endif
endif

.PHONY: default
default: build

.PHONY: documentation doc docs
documentation doc docs: document-api

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)

##############################################################################
# Documentation
#
.PHONY: document-api
document-api: PROJECT       = gungho
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation

##############################################################################
# Build
#
.PHONY: build
build: export BIN_DIR     ?= $(PROJECT_DIR)/bin
build: export CXX_LINK     = TRUE
build: export PROGRAMS    := $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -print)))
build: export WORKING_DIR := $(WORKING_DIR)

build: export LDFLAGS_GROUPS = OPENMP

# A 'kokkos' transformation emits C++ alongside the PSy layer. Only those
# options need a C++ toolchain configured, so every other transformation keeps
# building exactly as before.
#
# The test is on the first hyphen-separated word rather than the whole name, so
# that a variant of the kokkos transformation is configured like it. 'kokkos'
# alone was the condition until 'kokkos-timed' was added and silently got no
# CXXFLAGS: mpic++ was then invoked bare, and the generated translation units
# failed on '#include <Kokkos_Core.hpp>: No such file or directory' at the
# first compile after a PSyclone step that had gone perfectly. Matching the
# word keeps 'minimum-timed' and 'production-timed' out while admitting
# 'kokkos-timed' and any later 'kokkos-<variant>'.
ifeq ($(firstword $(subst -, ,$(PSYCLONE_TRANSFORMATION))),kokkos)
  ifndef KOKKOS_ROOT
    $(error PSYCLONE_TRANSFORMATION=kokkos requires KOKKOS_ROOT to point at a Kokkos installation)
  endif
  # USE_KOKKOS is set for both languages: it selects the live body of
  # driver_kokkos_mod.F90 and of kokkos_runtime.cpp, which between them are the
  # only place the Kokkos runtime is started and stopped. The C++ rule in
  # compile.mk passes CXXFLAGS but not the Fortran macro list, so the same
  # macro has to be named twice.
  build: export PRE_PROCESS_MACROS += USE_KOKKOS
  # KOKKOS_ACCEL names an accelerator variant of the image -- cuda90, cuda80 --
  # and switches every generated .cpp from a host compile to a device one. The
  # caller sets it, points OMPI_CXX at nvcc_wrapper, and nothing below has to
  # know which card: the architecture is read out of the Kokkos install a few
  # lines down. Unset, this file behaves exactly as it did.
  #
  # -fopenmp is the host compiler's flag and nvcc does not have it. Handed to
  # nvcc bare it is 'Unknown option', so it travels through -Xcompiler on an
  # accelerator build. Kokkos here is built with the OpenMP host backend as
  # well as CUDA, so it is still needed.
  KOKKOS_OPENMP_FLAG := -fopenmp
  ifdef KOKKOS_ACCEL
    KOKKOS_OPENMP_FLAG := -Xcompiler -fopenmp
  endif
  # -I$(WORKING_DIR) is for "lfric_kokkos_staging.hpp", the header every
  # generated region includes to obtain its Views. A capture writes it beside
  # the regions it generates; the build driver writes it into the working
  # directory as well, so that a region written anywhere under the tree finds
  # one copy of it.
  build: export CXXFLAGS += -std=c++17 $(KOKKOS_OPENMP_FLAG) -DUSE_KOKKOS -I$(KOKKOS_ROOT)/include -I$(WORKING_DIR)
  ifdef KOKKOS_ACCEL
    # The card is READ, not written. KokkosTargets.cmake records what a project
    # compiling against this Kokkos must pass, so the cuda80 image gets sm_80
    # without this file naming an image and an architecture cannot be misspelt
    # here into a build that compiles, links and then fails at its first launch
    # with 'no kernel image is available for execution on the device'.
    #
    # It cannot be left out: nvcc's own default is the toolkit's, sm_80 for
    # CUDA 13, so a region built without this line is a cubin for the wrong
    # card rather than no cubin at all -- a failure no compile-time error
    # reports. tests/test_cuda_model_build.sh in psy-ir-aidev is the assertion.
    KOKKOS_ARCH_FLAGS := $(shell grep -hoE -- '-arch=sm_[0-9]+' \
        $(wildcard $(KOKKOS_ROOT)/lib*/cmake/Kokkos/KokkosTargets.cmake) \
        2>/dev/null | sort -u)
    ifeq ($(KOKKOS_ARCH_FLAGS),)
      $(error KOKKOS_ACCEL=$(KOKKOS_ACCEL) but no -arch=sm_NN in $(KOKKOS_ROOT)/lib*/cmake/Kokkos/KokkosTargets.cmake: this Kokkos does not say which card it was built for, and a default written here would be a guess)
    endif
    # --expt-extended-lambda admits the __device__ lambda KOKKOS_LAMBDA expands
    # to; --expt-relaxed-constexpr lets device code call the constexpr host
    # functions Kokkos' own headers are full of. Neither is optional for a
    # translation unit that instantiates a Kokkos::parallel_for.
    build: export CXXFLAGS += $(KOKKOS_ARCH_FLAGS) --expt-extended-lambda --expt-relaxed-constexpr
    # KOKKOS_FMAD=off turns off nvcc's fused multiply-add contraction, so that
    # a device build's arithmetic can be compared with the host build's. It is
    # nvcc's flag and is rejected by a host compiler, which is why it is inside
    # this block and not beside the profile flags.
    ifeq "$(KOKKOS_FMAD)" "off"
      build: export CXXFLAGS += --fmad=false
    endif
  endif
  # The C++ optimisation level, matched to the Fortran side's for the same
  # profile. There is nothing to inherit: lfric_core's compiler tables define
  # flag groups for Fortran only, and compile.mk's .cpp rule passes CXXFLAGS
  # and nothing else. gfortran.mk maps NO_OPTIMISATION to -O0,
  # SAFE_OPTIMISATION to -Og and RISKY_OPTIMISATION to -Ofast; the three lines
  # below say the same thing to g++.
  #
  # WITHOUT THEM THE GENERATED C++ IS COMPILED AT -O0, because g++ optimises
  # nothing by default where the Fortran rule always passes a group. That is
  # not a small difference for this code. A generated kernel is Kokkos views,
  # lambdas and index arithmetic, nearly all of it in headers instantiated in
  # our own translation unit, so -O0 leaves every view accessor and every
  # functor call out of line -- while libkokkoscore itself, coming from Spack,
  # is optimised, which is why the omission does not show up as an obviously
  # broken build.
  #
  # It was found by measurement rather than by reading: stage 7 timed the
  # Kokkos regions against the Fortran loops they replaced and got ratios of
  # 6-32x, which sent someone looking for a compile line.
  # psy-ir-aidev/validation/performance/README.md records both readings.
  #
  # -Ofast carries -ffast-math, which is as much of a hazard in C++ as it is in
  # Fortran. It is mirrored rather than softened because a C++ side compiled
  # differently from the Fortran side is the defect being fixed here, and a
  # profile that wants IEEE semantics should say so for both languages at once.
  #
  # -g is on all three because the Fortran side's DEBUG group is: every
  # FFLAG_GROUPS line below names DEBUG, production included.
  #
  # There is deliberately no else-error here. The FFLAG_GROUPS block below
  # covers the same three profiles and errors on anything else, so an
  # unrecognised PROFILE stops the build there rather than reaching a link.
  # A second copy of that error would only be a second thing to keep in step.
  #
  # nvcc takes a NUMBER after -O and nothing else: -Og is 'nvcc fatal : "g":
  # expected a number' and -Ofast the same, before a single line of the region
  # is read. So an accelerator build maps the three profiles onto the levels
  # nvcc has. -Og becomes -O2 rather than -O1 because the Fortran side's
  # SAFE_OPTIMISATION is -Og, whose point is optimised-but-debuggable, and -O2
  # is the nearest nvcc offers; -Ofast becomes -O3 WITHOUT -ffast-math, which
  # nvcc would apply to the host pass only. The device-side counterpart of
  # fast-math is fused multiply-add contraction, which nvcc does by default and
  # KOKKOS_FMAD=off above turns off -- so the numerics knob for device code is
  # that one, not this.
  ifdef KOKKOS_ACCEL
    ifeq "$(PROFILE)" "full-debug"
      build: export CXXFLAGS += -O0 -g
    else ifeq "$(PROFILE)" "fast-debug"
      build: export CXXFLAGS += -O2 -g
    else ifeq "$(PROFILE)" "production"
      build: export CXXFLAGS += -O3 -g
    endif
  else ifeq "$(PROFILE)" "full-debug"
    build: export CXXFLAGS += -O0 -g
  else ifeq "$(PROFILE)" "fast-debug"
    build: export CXXFLAGS += -Og -g
  else ifeq "$(PROFILE)" "production"
    build: export CXXFLAGS += -Ofast -g
  endif
  build: export LDFLAGS += -L$(KOKKOS_ROOT)/lib -Wl,-rpath=$(KOKKOS_ROOT)/lib
  ifdef KOKKOS_ACCEL
    # THE CUDA RUNTIME HAS TO BE ON THE LINK LINE, and it is not enough that
    # libkokkoscore.so already needs it. An nvcc-compiled translation unit calls
    # cudaMemcpyToSymbol and friends ITSELF -- that is how a __device__ variable
    # in a region is initialised -- and ld's default
    # --no-copy-dt-needed-entries refuses to satisfy a reference from an object
    # out of a library that is only reachable through another library's NEEDED:
    #
    #   ld: lib/libgungho_model.a(smooth_orog_kokkos.o): undefined reference to
    #       symbol 'cudaMemcpyToSymbol@@libcudart.so.13'
    #   ld: /opt/view/lib64/libcudart.so.13: DSO missing from command line
    #
    # The Fortran executable is linked by mpif90, which knows nothing of CUDA,
    # so the toolkit's library directory is found the way the extracted-region
    # harness finds it: from nvcc's own position, not from a path written here.
    KOKKOS_CUDA_LIB := $(shell cd "$$(dirname "$$(command -v nvcc)")/../lib64" 2>/dev/null && pwd)
    ifeq ($(KOKKOS_CUDA_LIB),)
      $(error KOKKOS_ACCEL=$(KOKKOS_ACCEL) but no nvcc on PATH, so the CUDA runtime cannot be named on the link line)
    endif
    build: export LDFLAGS += -L$(KOKKOS_CUDA_LIB) -Wl,-rpath=$(KOKKOS_CUDA_LIB)
    build: export EXTERNAL_DYNAMIC_LIBRARIES += cudart
    # AND THE DRIVER LIBRARY, from the toolkit's STUB, which is linked against
    # and deliberately NOT put on an rpath. libkokkoscore.so needs libcuda.so.1
    # -- cuStreamGetCtx, cuCtxPushCurrent_v2, cuCtxGetDevice -- and that library
    # comes from a driver installation, not from a toolkit, so it is absent on
    # a host with no card and the executable cannot be linked there at all:
    #
    #   ld: libkokkoscore.so: undefined reference to `cuStreamGetCtx'
    #
    # The toolkit ships a stub carrying the same soname and no implementation,
    # which is exactly what the linker needs and nothing more. At run time the
    # loader takes libcuda.so.1 from the real driver that `apptainer --nv` binds
    # in on a GPU node. On the rpath instead, the stub would satisfy the loader
    # too, and the model would run against a library whose every entry point
    # returns an error.
    ifneq ($(wildcard $(KOKKOS_CUDA_LIB)/stubs),)
      build: export LDFLAGS += -L$(KOKKOS_CUDA_LIB)/stubs
      build: export EXTERNAL_DYNAMIC_LIBRARIES += cuda
    endif
  endif
  # libxios.a is built against the pre-C++11 libstdc++ string ABI; Kokkos is
  # built against the current one. The two disagree about the layout of
  # std::string, and so about any class holding one -- but the mangled names of
  # affected templates such as std::__detail::_Scanner<char> are identical in
  # both, because those names carry no ABI tag. Linking generated Kokkos code
  # into the executable makes the linker export the executable's copies, which
  # are XIOS's; ELF then gives them precedence over libkokkoscore's own, so
  # Kokkos calls into code compiled for the other ABI and segfaults before main.
  # Keeping XIOS's copies out of .dynsym leaves each library binding to its own,
  # which is what an unmixed build would have done.
  build: export LDFLAGS += -Wl,--exclude-libs,libxios.a
  build: export EXTERNAL_DYNAMIC_LIBRARIES += kokkoscore
endif

ifeq "$(PROFILE)" "full-debug"
build: export FFLAG_GROUPS = OPENMP DEBUG WARNINGS INIT RUNTIME NO_OPTIMISATION FORTRAN_STANDARD
else ifeq "$(PROFILE)" "fast-debug"
build: export FFLAG_GROUPS = OPENMP DEBUG WARNINGS SAFE_OPTIMISATION FORTRAN_STANDARD
else ifeq "$(PROFILE)" "production"
build: export FFLAG_GROUPS = OPENMP DEBUG WARNINGS RISKY_OPTIMISATION
else
  $(error Unrecognised profile "$(PROFILE)". Must be one of full-debug, fast-debug or production)
endif
# The namelist loaders are generated before the internal dependencies are
# imported, not after them, because importing a subproject PSyclones it, and a
# transformation may have to resolve a module that only exists once these have
# run. planet_config_mod is the case that forced the order: LFRicKokkosTrans
# refuses moist_dyn_factors_alg_mod in science/gungho with "cannot type
# 'recip_epsilon' without the source of 'planet_config_mod'", because
# rose_picker has not yet written that module into WORKING_DIR/configuration.
# psyclone_psykal.mk passes -d WORKING_DIR, which PSyclone searches
# recursively for modules, so generating them first is all that is needed.
#
# Neither step reads anything the import loop produces: extract_meta.mk works
# from APPS_ROOT_DIR and dependencies.yaml, and configuration.mk from
# META_FILE_DIR/rose-meta.conf and ROSE_META_DIRS. configuration.mk is passed
# SOURCE_DIR only because its header mentions one; no rule in it uses it.
#
# Extracting the application's own source stays where it was, after the
# import loop, so that the precedence between an application file and a
# subproject file of the same name is unchanged.
build: ALWAYS
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting External Rose Metadata)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(APPS_ROOT_DIR)/build/extract/extract_meta.mk
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Generating $(PROJECT) namelist loaders)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
	          PROJECT=gungho \
	          SOURCE_DIR=source \
	          META_FILE_DIR=$(META_FILE_DIR)
	$(call MESSAGE,========================================)
	$(call MESSAGE,Importing internal dependencies...)
	$(call MESSAGE,========================================)
	$(Q)for SUBPROJECT in $(INTERNAL_DEPENDENCIES) ; do \
	    $(MAKE) $(QUIET_ARG) -f $$SUBPROJECT/build/import.mk ; done
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting $(PROJECT_NAME))
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=source
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning $(PROJECT_NAME))
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=source
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Analysing $(PROJECT) build dependencies)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Compiling $(PROJECT))
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk

##############################################################################
# Pat_build step - This runs after the build step as relies on exe

# Pass $(PROJECT_NAME) which is exe name gungho_model which pat_report needs
pat_build: build
	$(call MESSAGE,=========================================================)
	$(call MESSAGE, Begin Pat Build of $(PROJECT_NAME) executable)
	$(call MESSAGE,=========================================================)
	bash $(PAT_BUILD_PATH) $(PROJECT_NAME)
	$(call MESSAGE,=========================================================)
	$(call MESSAGE, Pat Build of $(PROJECT_NAME) complete as $(PROJECT_NAME)+pat)
	$(call MESSAGE,=========================================================)

##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"$(PROJECT_NAME) work space")
	$(Q)-if [ $(WORKING_DIR) != *[\*]* ] && [ -d $(WORKING_DIR) ] ; then rm -r $(WORKING_DIR) ; fi
	$(call MESSAGE,Removing,"$(PROJECT_NAME) documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"$(PROJECT_NAME) binaries")
	$(Q)-if [ -d bin ] ; then rm -r bin ; fi
	$(Q)-if [ -d test ] ; then rm -r test ; fi
