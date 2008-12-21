# this makefile is an abbreviated form of src/snes/examples/tutorials/makefile

CFLAGS           =
CPPFLAGS         =
MANSEC           = SNES
EXAMPLESC	 = ex5.c poisson.c bugVecMax.c

include ${PETSC_DIR}/bmake/common/base


poisson: poisson.o  chkopts
	-${CLINKER} -o poisson poisson.o ${PETSC_SNES_LIB}
	${RM} poisson.o


bugVecMax: bugVecMax.o  chkopts
	-${CLINKER} -o bugVecMax bugVecMax.o ${PETSC_SNES_LIB}
	${RM} bugVecMax.o


ex5: ex5.o chkopts
	-${CLINKER} -o ex5 ex5.o ${PETSC_SNES_LIB}
	${RM} ex5.o


runex5:
	-@${MPIEXEC} -np 4 ./ex5 -snes_mf -da_processors_x 4 -da_processors_y 1 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5_1.tmp 2>&1; \
	   if (${DIFF} output/ex5_1.out ex5_1.tmp) then true; \
	   else echo "Possible problem with ex5_1, diffs above"; fi; \
	   ${RM} -f ex5_1.tmp
runex5_2:
	-@${MPIEXEC} -np 4 ./ex5 -da_processors_x 2 -da_processors_y 2 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5_2.tmp 2>&1; \
	   if (${DIFF} output/ex5_2.out ex5_2.tmp) then true; \
	   else echo "Possible problem with ex5_2, diffs above"; fi; \
	   ${RM} -f ex5_2.tmp
runex5_3:
	-@${MPIEXEC} -np 4 ./ex5 -snes_mf -da_processors_x 2 -da_processors_y 2 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5_3.tmp 2>&1;\
	   if (${DIFF} output/ex5_1.out ex5_3.tmp) then true; \
	   else echo "Possible problem with ex5_3, diffs above"; fi; \
	   ${RM} -f ex5_3.tmp
runex5_4:
	-@${MPIEXEC} -np 4 ./ex5 -da_processors_x 2 -da_processors_y 2 -snes_monitor_short -pc_type asm \
	   -pc_asm_blocks 4 -pc_asm_overlap 1 -ksp_gmres_cgs_refinement_type refine_always > ex5_4.tmp 2>&1; \
	   if (${DIFF} output/ex5_4.out ex5_4.tmp) then true; \
	   else echo "Possible problem with ex5_4, diffs above"; fi; \
	   ${RM} -f ex5_4.tmp
runex5_5:
	-@${MPIEXEC} -np 4 ./ex5 -da_processors_x 2 -da_processors_y 2 -snes_monitor_short -pc_type bjacobi \
	   -snes_view -ksp_gmres_cgs_refinement_type refine_always > ex5_5.tmp 2>&1; \
	   if (${DIFF} output/ex5_5.out ex5_5.tmp) then true; \
	   else echo "Possible problem with ex5_5, diffs above"; fi; \
	   ${RM} -f ex5_5.tmp
runex5_6:
	-@${MPIEXEC} -np 4 ./ex5 -snes_fd -da_processors_x 2 -da_processors_y 2 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5_6.tmp 2>&1; \
	   if (${DIFF} output/ex5_6.out ex5_6.tmp) then true; \
	   else echo "Possible problem with ex5_6, diffs above"; fi; \
	   ${RM} -f ex5_6.tmp
runex5_7:
	-@${MPIEXEC} -np 3 ./ex5 -snes_fd  -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5_7.tmp 2>&1;\
	   if (${DIFF} output/ex5_7.out ex5_7.tmp) then true; \
	   else echo "Possible problem with ex5_7, diffs above"; fi; \
	   ${RM} -f ex5_7.tmp
runex5_8:
	-@${MPIEXEC} -np 2 ./ex5 -snes_monitor_solution -ksp_gmres_cgs_refinement_type refine_always > /dev/null 2>&1
runex5_9:
	-@${MPIEXEC} -np 4 ./ex5 -snes_mf -snes_monitor_short>ex5_9.tmp 2>&1;\
	   if (${DIFF} output/ex5_9.out ex5_9.tmp) then true; \
	   else echo "Possible problem with ex5_9, diffs above"; fi; \
	   ${RM} -f ex5_9.tmp
runex5_10:
	-@${MPIEXEC} -np 2 ./ex5 -ksp_diagonal_scale -ksp_monitor_short -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5_10.tmp 2>&1;\
	   if (${DIFF} output/ex5_10.out ex5_10.tmp) then true; \
	   else echo "Possible problem with ex5_10, diffs above"; fi; \
	   ${RM} -f ex5_10.tmp
runex5_11:
	-@${MPIEXEC} -np 4 ./ex5 -snes_mf -da_processors_x 4 -da_processors_y 1 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always -matlab_function > ex5_11.tmp 2>&1; \
	   if (${DIFF} output/ex5_1.out ex5_11.tmp) then true; \
	   else echo "Possible problem with ex5_11, diffs above"; fi; \
	   ${RM} -f ex5_11.tmp
runex5_13:
	-@for F in " " -matlab_function ; do \
            for J in " " -fd_jacobian -adic_jacobian; do \
              ${MPIEXEC} -np 2 ./ex5 -snes_monitor_short $$F $$J > ex5_13.tmp 2>&1; \
              if (${DIFF} output/ex5_13.out ex5_13.tmp) then true; \
	      else echo "Possible problem with ex5_13, diffs above " $F $J; fi; \
              ${RM} -f ex5_13.tmp; \
              ${MPIEXEC} -np 2 ./ex5 -snes_monitor_short $$F $$J -adicmf_jacobian > ex5_13.tmp 2>&1; \
              if (${DIFF} output/ex5_13.out ex5_13.tmp) then true; \
	      else echo "Possible problem with ex5_13, diffs above " $F $J; fi; \
              ${RM} -f ex5_13.tmp; \
            done; \
          done
runex5_14:
	-@${MPIEXEC} -np 4 ./ex5 -snes_monitor_short -fd_jacobian_ghosted > ex5_14.tmp 2>&1; \
	   if (${DIFF} output/ex5_13.out ex5_14.tmp) then true; \
	   else echo "Possible problem with ex5_14, diffs above "; fi; \
           ${RM} -f ex5_14.tmp; 

TESTEXAMPLES_C		       = runex5 runex5_2 runex5_3 runex5_4 \
                                 runex5_5 runex5_6 runex5_7 runex5_9 runex5_13 runex5_14
TESTEXAMPLES_C_X11	       = runex5_8 runex5_10
TESTEXAMPLES_MATLAB_ENGINE     = runex5_11 runex5_13

include ${PETSC_DIR}/bmake/common/test

