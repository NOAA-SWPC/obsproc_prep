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
SRCS = 	prepacqc.f acftobs_qc.f tranQCflags.f output_acqc_prof.f output_acqc_noprof.f input_acqc.f indexc40.f sub2mem_mer.f sub2mem_um.f

OBJS =   prepacqc.o acftobs_qc.o tranQCflags.o output_acqc_prof.o output_acqc_noprof.o input_acqc.o indexc40.o sub2mem_mer.o sub2mem_um.o


# Tunable parameters
#
# FC		Name of the fortran compiling system to use
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
#
FC = ifort
#LDFLAGS = -Xlinker -t
#LDFLAGS = -Xlinker -M
LIBS = -L/nwprod/lib -lw3nco_4 -lw3emc_4 -lbufr_4_64
CMD = prepobs_prepacqc

##DEBUG =		 -ftrapuv -check all -check nooutput_conversion -fp-stack-check -fstack-protector
FFLAGS = -O2 -convert big_endian -list -assume noold_ldout_format -fixed -g -traceback $(DEBUG)

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
