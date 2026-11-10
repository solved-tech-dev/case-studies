include Makefile.inc
BIN = ../bin
all  :  fwt

%.o : %.mod


.SUFFIXES: .f90 .o
.f90.o:
	$(FC) $(OPTFF) $(INC) -c $*.f90

.SUFFIXES: .f .o
.f.o:
	$(FC) $(OPTF) $(INC) -c $*.f
.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) $(INC) -c $*.c


clean:
	$(RM) *.o  *~ *.mod  \#* *_genmod.f90


SUBINV = MPI_TIME.o\
flushc.o \
flush.o \
func_test.o\
sub_MPI_ABORT.o\
sub_read_freqdom_freq_list.o\
sub_damp.o\
sub_modext.o\
sub_average.o\
sub_read_acqui.o\
sub_read_models.o \
sub_PML.o\
sub_evalmatrix.o\
sub_evalmatrix_fs.o\
sub_evalmatrix_TTI.o\
sub_evalmatrix_TTI_fs.o\
sub_rho2b.o\
sub_v2mu.o\
sub_vcmplx.o\
sub_Hicks_constant.o\
sub_Hicks_extraction.o\
sub_Hicks.o\
sub_sample_Hicks.o\
sub_dAdmi.o\
sub_gradient_for_base_param.o\
sub_drho5.o\
sub_d2rho5.o\
sub_FDFD.o\
sub_fill_RHS_SOURCES.o\
sub_fill_RHS_RECEIVERS.o\
sub_fill_RHS_SOURCES_HICKS.o\
sub_fill_RHS_RECEIVERS_HICKS.o\
sub_fill_RHS_Jv.o\
sub_fill_RHS_lambda.o\
sub_restriction.o\
sub_Tikhonov.o\
sub_precond.o\
sub_prior.o\
sub_modeling.o\
sub_FWI.o\
sub_converter.o\
sub_freeze.o\
sub_Shin_preco.o\
sub_init_optim.o\
sub_FWI_OPTIM.o
#sub_FWI_PSTD.o\
#sub_FWI_PNLCG.o\
#sub_FWI_LBFGS.o\
#sub_FWI_PLBFGS.o\
#sub_FWI_TRN.o\
#sub_FWI_PTRN.o\

#DEPENDANCIES WITH MAKEFILE AND .h include files
$(SUBINV) toy2dac.o: ../include/*.h Makefile

# INVERSION CODE
fwt : $(SUBINV) toy2dac.o
	$(FL) -o $(BIN)/toy2dac$(EXT) $(OPTL) $(SUBINV) toy2dac.o $(LIB) 

