PROCESSOR := $(shell uname -m)
F90=gfortran
FFLAGS=-g -C -O2 -ffree-form -I/opt/local/include #-fbounds-check -Wtabs -fcheck=all

OBJFLAGS ?=$(FFLAGS)
LDFLAGS= -L/opt/local/lib -lnetcdf -lnetcdff -framework vecLib

SRCDIR = _src

SOURCES = $(SRCDIR)/tfcn.f90 \
					$(SRCDIR)/driver.f90 \
					$(SRCDIR)/qinit.f90 \
					$(SRCDIR)/init2d.f90 \
					$(SRCDIR)/output2d.f90\
					$(SRCDIR)/strangSplit.f90\
					$(SRCDIR)/updateVelocities.f90\
					$(SRCDIR)/updateSoln1d.f90\
					$(SRCDIR)/projectAverages.f90\
					$(SRCDIR)/numFlux.f90\
					$(SRCDIR)/fluxFunction.f90\
					$(SRCDIR)/evaluateExpansion.f90\

MODULES = $(SRCDIR)/mDGmod.f90 \
					$(SRCDIR)/commonTestParameters.f90 \

INPUTS  = inputs.nl

.PHONY= 2d_test clean

OBJECTS :=$(notdir $(SOURCES:.f90=.o))
MODOBJ  :=$(notdir $(MODULES:.f90=.o))

# Set the order in which directories are searched when looking for targets
VPATH = $(SRCDIR)

all: $(SOURCES) $(INPUTS) test_2d_modal

2d_test: test_2d_modal
	./test_2d_modal

test_2d_modal: $(MODOBJ) $(OBJECTS) $(INPUTS) split_2d_modal.f90
	$(F90) $(FFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(MODOBJ) *.mod test_2d_modal

%.o: %.f90
		$(F90) $(OBJFLAGS) -c -o $@ $^
