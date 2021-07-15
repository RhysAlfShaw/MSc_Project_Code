##############################################################
###                                                        ###
###           Make File for the NBD Code                   ###
###                 pcc - 18/10/2020                       ###
###                                                        ###
##############################################################

.KEEP_STATE:

.SUFFIXES: .f90 .o
FC = gfortran 
FFLAGS = -O3 -ffixed-line-length-132

FOBS = global_parameters.o global_variables.o \
	nbd.o allocate_memory.o \
	get_accel_info_faster.o \
	get_accel_for_active.o \
	get_ics.o \
	hermite_4thorder_PECn_fixed_dt.o \
        hermite_4thorder_PECn_block_dt.o \
	hermite_4thorder_PECn_global_dt.o \
	initialise_run.o \
	read_param_file.o write_a_snapshot.o \
	diagnostics.o get_energies_and_momentum.o \
	synchronise_posvel_for_output.o

EXEC = nbd

.f90.o :
	$(FC) $(FLAGS) -c $*.f90

$(EXEC): $(FOBS)
	$(FC) $(FOBS) $(FFLAGS) -o $(EXEC)

clean:
	\rm -f *.o *.mod $(EXEC)
