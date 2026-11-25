# =================================================================
# Intel ifx Makefile for Lattice Boltzmann (Fortran + Coarrays)
# =================================================================

FC      = ifx
EXE     = lbm
SRCDIR  = .
OBJDIR  = objects
MODDIR  = $(OBJDIR)/mod
CONFIG  = lbm.conf

# Base compilation flags
FFLAGS  = -free -O3 -coarray=distributed -coarray-config-file=$(CONFIG) -module $(MODDIR)

# Sources and Objects
SRC_MOD = $(SRCDIR)/LatticeBoltzmann.f90
SRC_APP = $(SRCDIR)/LatticeBoltzmannMain.f90

OBJ_MOD = $(OBJDIR)/LatticeBoltzmann.o
OBJ_APP = $(OBJDIR)/LatticeBoltzmannMain.o

OBJS    = $(OBJ_MOD) $(OBJ_APP)

# -------------------------------------------------
# Build Rules
# -------------------------------------------------
all: directories $(EXE)

directories:
	@mkdir -p $(OBJDIR) $(MODDIR)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

# Compile the module first
$(OBJ_MOD): $(SRC_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

# Compile the main program
$(OBJ_APP): $(SRC_APP) | $(OBJ_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR) $(EXE)

.PHONY: all directories clean