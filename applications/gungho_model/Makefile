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
  build: export CXXFLAGS += -std=c++17 -fopenmp -DUSE_KOKKOS -I$(KOKKOS_ROOT)/include
  build: export LDFLAGS += -L$(KOKKOS_ROOT)/lib -Wl,-rpath=$(KOKKOS_ROOT)/lib
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
