#
#  This file is part of MUMPS 5.0.1, released
#  on Thu Jul 23 17:08:29 UTC 2015
#
#
#  This file is part of MUMPS 5.0beta, built on Tue Jan 10 16:43:31 UTC 2012
#
topdir = ..
libdir = $(topdir)/lib

default:	d

.PHONY: default s d c z mumps_lib clean

s:
	$(MAKE) ARITH=s mumps_lib
d:
	$(MAKE) ARITH=d mumps_lib
c:
	$(MAKE) ARITH=c mumps_lib
z:
	$(MAKE) ARITH=z mumps_lib

include $(topdir)/Makefile.inc

mumps_lib:    $(libdir)/libmumps_common$(PLAT)$(LIBEXT) \
              $(libdir)/lib$(ARITH)mumps$(PLAT)$(LIBEXT)

OBJS_COMMON_MOD = \
        ana_omp_m.o\
	tools_common_mod.o\
        mumps_static_mapping.o\
        mumps_sol_es.o\
	fac_future_niv2_mod.o\
        mumps_comm_ibcast.o\
	mumps_ooc_common.o\
        double_linked_list.o\
        fac_asm_build_sort_index_m.o\
        fac_asm_build_sort_index_ELT_m.o\
        omp_tps_common_m.o\
        mumps_l0_omp_m.o\
        front_data_mgt_m.o\
        fac_maprow_data_m.o\
        fac_descband_data_m.o\
        fac_ibct_data_m.o

OBJS_COMMON_OTHER = \
        ana_orderings.o\
        ana_AMDMF.o\
        bcast_errors.o\
        estim_flops.o\
        mumps_type_size.o \
        mumps_type2_blocking.o \
        mumps_version.o \
        tools_common.o \
        mumps_print_defined.o \
        mumps_common.o\
        mumps_orderings.o\
        mumps_size.o\
        mumps_io.o\
        mumps_io_basic.o\
        mumps_io_thread.o\
        mumps_io_err.o\
        ana_set_ordering.o\
        mumps_numa.o

OBJS_MOD =   \
        $(ARITH)ana_aux_par.o \
        $(ARITH)ana_lr.o\
        $(ARITH)fac_asm_master_m.o\
        $(ARITH)fac_asm_master_ELT_m.o\
        $(ARITH)omp_tps_m.o\
        $(ARITH)mumps_comm_buffer.o\
        $(ARITH)mumps_load.o\
        $(ARITH)mumps_lr_data_m.o\
        $(ARITH)mumps_ooc_buffer.o\
        $(ARITH)mumps_ooc.o\
        $(ARITH)mumps_struc_def.o\
        $(ARITH)static_ptr_m.o\
        $(ARITH)lr_core.o\
        $(ARITH)fac_lr.o\
        $(ARITH)fac_omp_m.o\
        $(ARITH)fac_par_m.o\
        $(ARITH)fac_front_LU_type1.o\
        $(ARITH)fac_front_LU_type2.o\
        $(ARITH)fac_front_LDLT_type1.o\
        $(ARITH)fac_front_LDLT_type2.o\
        $(ARITH)fac_front_aux.o\
        $(ARITH)fac_front_type2_aux.o

OBJS_OTHER = \
        $(ARITH)ini_driver.o\
        $(ARITH)ana_driver.o\
        $(ARITH)fac_driver.o\
        $(ARITH)sol_driver.o\
        $(ARITH)end_driver.o\
        $(ARITH)ana_aux_ELT.o\
        $(ARITH)ana_aux.o\
        $(ARITH)ana_dist_m.o\
        $(ARITH)ana_LDLT_preprocess.o\
        $(ARITH)ana_reordertree.o\
        $(ARITH)arrowheads.o\
        $(ARITH)bcast_int.o\
        $(ARITH)fac_asm_ELT.o\
        $(ARITH)fac_asm.o\
        $(ARITH)fac_b.o\
        $(ARITH)fac_distrib_distentry.o\
        $(ARITH)fac_distrib_ELT.o\
        $(ARITH)fac_lastrtnelind.o\
        $(ARITH)fac_mem_alloc_cb.o\
        $(ARITH)fac_mem_compress_cb.o\
        $(ARITH)fac_mem_free_block_cb.o\
        $(ARITH)fac_mem_stack_aux.o\
        $(ARITH)fac_mem_stack.o\
        $(ARITH)fac_process_band.o\
        $(ARITH)fac_process_blfac_slave.o\
        $(ARITH)fac_process_blocfacto_LDLT.o\
        $(ARITH)fac_process_blocfacto.o\
        $(ARITH)fac_process_bf.o\
        $(ARITH)fac_process_end_facto_slave.o\
        $(ARITH)fac_process_contrib_type1.o\
        $(ARITH)fac_process_contrib_type2.o\
        $(ARITH)fac_process_contrib_type3.o\
        $(ARITH)fac_process_maprow.o\
        $(ARITH)fac_process_master2.o\
        $(ARITH)fac_process_message.o\
        $(ARITH)fac_process_root2slave.o\
        $(ARITH)fac_process_root2son.o\
        $(ARITH)fac_process_rtnelind.o\
        $(ARITH)fac_root_parallel.o\
        $(ARITH)fac_scalings.o\
        $(ARITH)fac_determinant.o\
        $(ARITH)fac_scalings_simScaleAbs.o\
        $(ARITH)fac_scalings_simScale_util.o\
        $(ARITH)fac_sol_pool.o\
        $(ARITH)fac_type3_symmetrize.o\
        $(ARITH)ini_defaults.o\
        $(ARITH)mumps_c.o\
        $(ARITH)mumps_driver.o\
        $(ARITH)mumps_f77.o\
        $(ARITH)mumps_iXamax.o\
        $(ARITH)ana_mtrans.o\
        $(ARITH)ooc_panel_piv.o\
        $(ARITH)rank_revealing.o\
        $(ARITH)sol_aux.o\
        $(ARITH)sol_bwd_aux.o\
        $(ARITH)sol_bwd.o\
        $(ARITH)sol_c.o\
        $(ARITH)sol_fwd_aux.o\
        $(ARITH)sol_fwd.o\
        $(ARITH)sol_matvec.o\
        $(ARITH)sol_root_parallel.o\
        $(ARITH)tools.o\
        $(ARITH)type3_root.o\
        $(ARITH)mumps_restart.o

$(libdir)/libmumps_common$(PLAT)$(LIBEXT):      $(OBJS_COMMON_MOD) $(OBJS_COMMON_OTHER)
	$(AR)$@ $?
	$(RANLIB) $@

$(libdir)/lib$(ARITH)mumps$(PLAT)$(LIBEXT):    $(OBJS_MOD) $(OBJS_OTHER)
	$(AR)$@ $?
	$(RANLIB) $@

# Dependencies between modules:
$(ARITH)mumps_load.o:		$(ARITH)mumps_comm_buffer.o \
				$(ARITH)mumps_struc_def.o \
				fac_future_niv2_mod.o
$(ARITH)mumps_ooc.o: 		$(ARITH)mumps_struc_def.o \
				$(ARITH)mumps_ooc_buffer.o \
				mumps_ooc_common.o
$(ARITH)mumps_driver.o: 	$(ARITH)static_ptr_m.o
$(ARITH)mumps_ooc_buffer.o: 	mumps_ooc_common.o
$(ARITH)ana_aux_par.o:          $(ARITH)mumps_struc_def.o \
				tools_common_mod.o
$(ARITH)mumps_comm_buffer.o:	mumps_comm_ibcast.o
$(ARITH)fac_asm_master_m.o:	omp_tps_common_m.o \
                                fac_ibct_data_m.o \
				$(ARITH)omp_tps_m.o \
				fac_asm_build_sort_index_m.o \
				$(ARITH)mumps_comm_buffer.o \
				$(ARITH)mumps_load.o
$(ARITH)fac_lastrtnelind.o:     $(ARITH)mumps_load.o
$(ARITH)fac_asm_master_ELT_m.o:	omp_tps_common_m.o \
                                fac_ibct_data_m.o \
				$(ARITH)omp_tps_m.o \
				fac_asm_build_sort_index_ELT_m.o \
				$(ARITH)mumps_comm_buffer.o \
				$(ARITH)mumps_load.o
$(ARITH)fac_omp_m.o:		$(ARITH)fac_asm_master_m.o\
				$(ARITH)fac_asm_master_ELT_m.o\
				$(ARITH)fac_front_LU_type1.o\
				$(ARITH)fac_front_LDLT_type1.o\
				$(ARITH)mumps_load.o\
				$(ARITH)omp_tps_m.o\
				omp_tps_common_m.o\
				mumps_l0_omp_m.o
$(ARITH)ana_lr.o:		$(ARITH)lr_core.o
$(ARITH)mumps_lr_data_m.o:	$(ARITH)lr_core.o\
                                $(ARITH)mumps_struc_def.o
$(ARITH)fac_lr.o:		$(ARITH)lr_core.o\
                                $(ARITH)ana_lr.o
$(ARITH)fac_par_m.o:            $(ARITH)mumps_load.o\
				$(ARITH)mumps_ooc.o\
				$(ARITH)fac_lr.o\
				$(ARITH)fac_asm_master_m.o\
				$(ARITH)fac_asm_master_ELT_m.o\
                                $(ARITH)omp_tps_m.o\
                                $(ARITH)fac_front_LU_type1.o\
                                $(ARITH)fac_front_LU_type2.o\
                                $(ARITH)fac_front_LDLT_type1.o\
                                $(ARITH)fac_front_LDLT_type2.o\
                                omp_tps_common_m.o\
				mumps_l0_omp_m.o

$(ARITH)fac_front_aux.o     : $(ARITH)mumps_comm_buffer.o\
                                 $(ARITH)mumps_load.o\
                                 $(ARITH)mumps_ooc.o\
                                 mumps_ooc_common.o\
                                 mumps_l0_omp_m.o\
                                 mumps_comm_ibcast.o

$(ARITH)fac_front_type2_aux.o : $(ARITH)mumps_ooc.o\
                                $(ARITH)fac_front_aux.o

$(ARITH)fac_front_LU_type1.o   : $(ARITH)fac_front_aux.o\
                                 $(ARITH)mumps_ooc.o\
                                 $(ARITH)lr_core.o\
                                 $(ARITH)ana_lr.o\
                                 $(ARITH)fac_lr.o
$(ARITH)fac_front_LU_type2.o   : $(ARITH)fac_front_aux.o\
                                 $(ARITH)fac_front_type2_aux.o\
                                 $(ARITH)mumps_ooc.o\
                                 $(ARITH)lr_core.o\
                                 $(ARITH)ana_lr.o\
                                 $(ARITH)fac_lr.o\
                                 fac_ibct_data_m.o

$(ARITH)fac_front_LDLT_type1.o : $(ARITH)fac_front_aux.o\
                                 $(ARITH)mumps_ooc.o

$(ARITH)fac_front_LDLT_type2.o : $(ARITH)fac_front_aux.o\
                                 $(ARITH)fac_front_type2_aux.o\
                                 $(ARITH)mumps_ooc.o\
                                 $(ARITH)mumps_load.o\
                                 fac_ibct_data_m.o

fac_maprow_data_m.o:            front_data_mgt_m.o
fac_descband_data_m.o:          front_data_mgt_m.o
fac_ibct_data_m.o  :            front_data_mgt_m.o
ana_omp_m.o:              	double_linked_list.o
mumps_comm_ibcast.o:		fac_future_niv2_mod.o
fac_asm_build_sort_index_m.o:	omp_tps_common_m.o
fac_asm_build_sort_index_ELT_m.o:omp_tps_common_m.o

# Compile modules before the rest
$(OBJS_COMMON_OTHER):$(OBJS_COMMON_MOD)
$(OBJS_OTHER):$(OBJS_COMMON_MOD) $(OBJS_MOD)

.SUFFIXES: .c .F .o
.F.o:
	$(FC) $(OPTF) $(INCS) $(IORDERINGSF) $(ORDERINGSF) -I. -I../include -c $*.F $(OUTF)$*.o
.c.o:
	$(CC) $(OPTC) $(INCS) -I../include $(CDEFS) $(IORDERINGSC) $(ORDERINGSC) -c $*.c $(OUTC)$*.o

$(ARITH)mumps_c.o:	mumps_c.c
	$(CC) $(OPTC) $(INCS) $(CDEFS) -DMUMPS_ARITH=MUMPS_ARITH_$(ARITH) \
	      $(IORDERINGSC) $(ORDERINGSC) -I../include -c mumps_c.c $(OUTC)$@


clean:
	$(RM) *.o *.mod
