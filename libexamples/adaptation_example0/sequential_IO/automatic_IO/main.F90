!>
!> Example of use of the parmmg library (basic use of mesh adaptation)
!>
!> \author Algiane Froehly (InriaSoft)
!> \version 1
!> \copyright GNU Lesser General Public License.

PROGRAM main

  IMPLICIT NONE

  !> Include the parmmg library hader file
  ! if the header file is in the "include" directory
  ! #include "libparmmgf.h"
  ! if the header file is in "include/parmmg"
#include "parmmg/libparmmgf.h"

  MMG5_DATA_PTR_T    :: parmesh
  INTEGER            :: ier,rank,i,nsols,argc,hasMet,hasSol
  CHARACTER(len=300) :: exec_name,filename,metname,solname,fileout,tmp
  CHARACTER(len=312) :: solout

  INTEGER            :: typSol(MMG5_NSOLS_MAX)

  PRINT*,"  -- FORTRAN TEST FOR PARMMGLIB"

  CALL MPI_Init( ier )
  CALL MPI_Comm_rank( MPI_COMM_WORLD, rank, ier )

  IF ( rank==0 ) THEN
     PRINT*,"  -- TEST PARMMGLIB"
  ENDIF

  argc =  COMMAND_ARGUMENT_COUNT();
  CALL get_command_ARGUMENT(0, exec_name)


  IF ( (argc<2) .AND. rank == 0 ) THEN
     PRINT*, " Usage: ",TRIM(ADJUSTL(exec_name)),&
          " filein fileout [[-sol metfile]/[-met metfile]] [-field solfile]"
     CALL EXIT(1);
  ENDIF

  CALL get_command_ARGUMENT(1, filename)
  CALL get_command_ARGUMENT(2, fileout)

  hasSol  = 0
  hasMet  = 0
  IF ( argc>2 ) THEN
     i = 3
     DO WHILE (i<=argc)
        CALL get_command_ARGUMENT(i, tmp)

        IF ( TRIM(tmp) ==  "-met" .OR. TRIM(tmp) ==  "-sol" ) THEN
           i = i + 1
           hasMet = 1
           CALL get_command_ARGUMENT(i, metname)
        ELSEIF ( TRIM(tmp) ==  "-field" ) THEN
           i = i + 1
           hasSol = 1
           CALL get_command_ARGUMENT(i, solname)
           solout = trim(fileout) // "-solphys.sol"
        ELSE
           PRINT *, "Unexpected argument: ", tmp
           CALL MPI_Abort(MPI_COMM_WORLD,1,ier)
        ENDIF
        i = i + 1
     ENDDO
  ENDIF

  !> ------------------------------ STEP   I --------------------------
  !> 1) Initialisation of th parmesh structures
  !! args of InitMesh:
  ! PMMG_ARG_start: we start to give the args of a variadic func
  ! PMMG_ARG_ppParMesh: next arg will be a pointer over a PMMG_pParMesh
  ! &parmesh: pointer toward your PMMG_pParMesh
  ! PMMG_ARG_pMesh: initialization of a mesh inside the parmesh.
  ! PMMG_ARG_pMet: init a metric inside the parmesh
  ! PMMG_ARG_dim: next arg will be the mesh dimension
  ! 3: mesh dimension
  ! PMMG_MPIComm: next arg will be the MPI COmmunicator
  ! MPI_COMM_WORLD: MPI communicator
  !

  parmesh = 0

  CALL PMMG_Init_parMesh(PMMG_ARG_start,                 &
       PMMG_ARG_ppParMesh,parmesh,                       &
       PMMG_ARG_pMesh,PMMG_ARG_pMet,                     &
       PMMG_ARG_dim,%val(3),PMMG_ARG_MPIComm,%val(MPI_COMM_WORLD), &
       PMMG_ARG_end);

  !> 2) Build mesh in PMMG format
  !> Two solutions: just use the PMMG_loadMesh_centralized function that will
  !> read a .mesh(b) file formatted or manually set your mesh using the
  !> PMMG_Set* functions

  !> with PMMG_loadMesh_centralized function
  CALL PMMG_loadMesh_centralized(parmesh,trim(filename),len(trim(filename)),ier)
  IF (ier  .NE.  1) THEN
     CALL MPI_Abort(MPI_COMM_WORLD,1,ier);
  ENDIF

  !> 3) Try to load a metric in PMMG format
  !> Two solutions: just use the PMMG_loadMet_centralized function that will
  !>    read a .sol(b) file formatted or manually set your metric using the PMMG_Set*
  !>    functions

  !> With PMMG_loadMet_centralized function
  IF ( hasMet == 1  ) THEN
     CALL PMMG_loadMet_centralized(parmesh,trim(metname),len(trim(metname)),ier)
     IF ( ier .NE. 1 ) THEN
        CALL MPI_Abort(MPI_COMM_WORLD,1,ier);
     ENDIF
  ENDIF

  !> 4) Build solutions in PMMG format
  !> Two solutions: just use the PMMG_loadAllSols_centralized function that
  !! will read a .sol(b) file formatted or manually set your solutions using
  !! the PMMG_Set* functions

  !> With PMMG_loadAllSols_centralized function
  IF ( hasSol == 1 ) THEN
     CALL PMMG_loadAllSols_centralized(parmesh,trim(solname),len(trim(solname)),ier)
     IF ( ier /= 1 )THEN
        CALL MPI_Abort(MPI_COMM_WORLD,1,ier)
     ENDIF
  ENDIF

  !> ------------------------------ STEP  II --------------------------
  !>  No surface adaptation */
  CALL PMMG_Set_iparameter( parmesh, PMMG_IPARAM_nosurf, 1, ier )
  IF ( ier /= 1 )THEN
     CALL MPI_Abort(MPI_COMM_WORLD,1,ier)
  ENDIF

  !> remesh function
  CALL PMMG_parmmglib_centralized(parmesh,ier)

  IF ( ier /= PMMG_STRONGFAILURE ) THEN
     !> ------------------------------ STEP III --------------------------
     !> get results
     !> Two solutions: just use the PMMG_saveMesh/PMMG_saveSol functions
     !! that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
     !! using the PMMG_getMesh/PMMG_getSol functions

     !> 1) Automatically save the mesh
     CALL PMMG_saveMesh_centralized(parmesh,trim(fileout),len(trim(fileout)),ier)
     IF ( ier /= 1 ) THEN
        PRINT*,"UNABLE TO SAVE MESH"
        ier = PMMG_STRONGFAILURE
     ENDIF

     !> 2) Automatically save the metric
     CALL PMMG_saveMet_centralized(parmesh,trim(fileout),len(trim(fileout)),ier)
     IF ( ier==0 ) THEN
        PRINT*,"UNABLE TO SAVE METRIC"
        ier = PMMG_LOWFAILURE
     ENDIF

     !> 3) Automatically save the solutions if needed
     CALL PMMG_Get_solsAtVerticesSize(parmesh,nsols,%val(0),typsol,ier)
     IF ( nsols /= 0 ) THEN
        CALL PMMG_saveAllSols_centralized(parmesh,trim(solout),len(trim(solout)),ier)
        IF ( ier /= 1 ) THEN
           PRINT*,"UNABLE TO SAVE SOLUTIONS"
           ier = PMMG_LOWFAILURE
        ENDIF
     ENDIF
  ELSE
     PRINT*,"BAD ENDING OF PARMMGLIB: UNABLE TO SAVE MESH"
  ENDIF

  !> 4) Free the PMMG structures
  CALL PMMG_Free_all ( PMMG_ARG_start,     &
       PMMG_ARG_ppParMesh,parmesh,         &
       PMMG_ARG_end);


  CALL MPI_Finalize(ier);

END PROGRAM main
