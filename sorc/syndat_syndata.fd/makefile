SHELL=/bin/sh
#
# If this makefile is invoked by the command line
#       make -f makefile
# it will compile the fortran modules indicated by SRCS into the object
# modules indicated by OBJS and produce an executable named a.out.
#
# If it is invoked by the command line
#       make -f makefile a.out.prof
# it will compile the fortran modules indicated by SRCS into the object
# modules indicated by OBJS and produce an executable which profiles
# named a.out.prof.
#
# To remove all the objects but leave the executables use the command line
#       make -f makefile clean
#
# To remove everything but the source files use the command line
#       make -f makefile clobber
#
# To remove the source files created by /usr/bin/fmgen and this makefile
# use the command line
#       make -f makefile void
#
# The parameters SRCS and OBJS should not need to be changed.  If, however,
# you need to add a new module add the name of the source module to the
# SRCS parameter and add the name of the resulting object file to the OBJS
# parameter.  The new modules are not limited to fortran, but may be C, YACC,
# LEX, or CAL.  An explicit rule will need to be added for PASCAL modules.
#
SRCS = 	syndata.f

OBJS =  syndata.o


# Tunable parameters
#
# FC		Name of the fortran compiling system to use
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
#
FC = ifort
#LDFLAGS = -Xlinker -t -mkl
#LDFLAGS = -Xlinker -M -mkl
LDFLAGS = -mkl
LIBS = -L/nwprod/lib -lw3nco_4 -lw3emc_4 -lbufr_4_64 -lsigio_4 -lsp_4 -lbacio_4
CMD = syndat_syndata

# For now DO NOT comment out DEBUG options - want to run in ops
#  with debugger on to make sure no problems
DEBUG =		 -ftrapuv -check all -check nooutput_conversion -fp-stack-check -fstack-protector
FFLAGS = -O2 -convert big_endian -list -assume noold_ldout_format -auto -g -traceback $(DEBUG)

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.
#
all:		$(CMD)

$(CMD):		$(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)


clean:
	-rm -f $(OBJS)

clobber:	clean
	-rm -f $(CMD) $(CMD).prof

void:	clobber
	-rm -f $(SRCS) makefile

install:
	mv $(CMD) ../../exec
