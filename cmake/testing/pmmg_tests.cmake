IF( BUILD_TESTING )
  include( CTest )

  set( CI_DIR  ${CMAKE_BINARY_DIR}/testparmmg CACHE PATH "path to test meshes repository" )
  set( CI_DIR_RESULTS  ${CMAKE_BINARY_DIR}/TEST_OUTPUTS )
  file( MAKE_DIRECTORY ${CI_DIR_RESULTS} )
  get_filename_component(PARENT_DIR ${CI_DIR} DIRECTORY)


  IF ( NOT ONLY_LIBRARY_TESTS )

    IF ( NOT EXISTS ${CI_DIR} )
      EXECUTE_PROCESS(
        COMMAND ${GIT_EXECUTABLE} clone https://gitlab.inria.fr/ParMmg/testparmmg.git
        WORKING_DIRECTORY ${PARENT_DIR}
        )
      EXECUTE_PROCESS(
        COMMAND ${GIT_EXECUTABLE} checkout 2182d86771506d99971226a4ae2a9eab137e127f
        WORKING_DIRECTORY ${CI_DIR}
        )
    ENDIF()

    set ( mesh_size 16384 )
    set ( myargs -niter 2 -metis-ratio 82 -v 5 )

    # remesh 2 sets of matching mesh/sol files (which are the output of mmg3d)
    # on 1,2,4,6,8 processors
    foreach( MESH cube-unit-dual_density cube-unit-int_sphere )
      foreach( NP 1 2 4 6 8 )
        add_test( NAME ${MESH}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${CI_DIR}/Cube/${MESH}.mesh
          -out ${CI_DIR_RESULTS}/${MESH}-${NP}-out.mesh
          -m 11000 -mesh-size ${mesh_size} ${myargs})
      endforeach()
    endforeach()

    # remesh a unit cube with two different solution files on 1,2,4,6,8 processors
    foreach( MESH dual_density int_sphere )
      foreach( NP 1 2 4 6 8 )
        add_test( NAME cube-unit-coarse-${MESH}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${CI_DIR}/Cube/cube-unit-coarse.mesh
          -sol ${CI_DIR}/Cube/cube-unit-coarse-${MESH}.sol
          -out ${CI_DIR_RESULTS}/${MESH}-${NP}-out.mesh
          -mesh-size ${mesh_size} ${myargs} )
      endforeach()
    endforeach()

    # remesh a non constant anisotropic test case: a torus with a planar shock
    # on 1,2,4,6,8 processors
    foreach( TYPE anisotropic-test )
      foreach( NP 1 2 4 6 8 )
        add_test( NAME ${TYPE}-torus-with-planar-shock-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${CI_DIR}/Torus/torusholes.mesh
          -sol ${CI_DIR}/Torus/torusholes.sol
          -out ${CI_DIR_RESULTS}/${TYPE}-torus-with-planar-shock-${NP}-out.mesh
          -mesh-size ${mesh_size} ${myargs} -nosurf )
      endforeach()
    endforeach()

    ###############################################################################
    #####
    #####        Tests options (on 1, 6 and 8 procs)
    #####
    ###############################################################################

    # Default option: no metric
    foreach( NP 1 6 8 )
      add_test( NAME Sphere-${NP}
        COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
        ${CI_DIR}/Sphere/sphere.mesh
        -out ${CI_DIR_RESULTS}/sphere-${NP}-out.mesh
        -mesh-size ${mesh_size} ${myargs} )
    endforeach()

    # Option without arguments
    foreach( OPTION "optim" "optimLES" "nosurf" "noinsert" "noswap"  )
      foreach( NP 1 6 8 )
        add_test( NAME Sphere-optim-${OPTION}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          -${OPTION}
          ${CI_DIR}/Sphere/sphere.mesh
          -out ${CI_DIR_RESULTS}/sphere-${OPTION}-${NP}-out.mesh
          -mesh-size ${mesh_size} ${myargs} )
      endforeach()
    endforeach()

    # Option with arguments
    SET ( OPTION
      "-v"
      "-hsiz"
      "-hausd"
      "-hgrad"
      "-hgrad"
      "-hmax"
      "-nr"
      "-ar" )

    SET ( VAL
      "5"
      "0.05"
      "0.005"
      "1.1"
      "-1"
      "0.05"
      ""
      "10" )

    SET ( NAME
      "v5"
      "hsiz0.05"
      "hausd0.005"
      "hgrad1.1"
      "nohgrad"
      "hmax0.05"
      "nr"
      "ar10" )

    SET ( MESH_SIZE
      "16384"
      "163840"
      "16384"
      "16384"
      "16384"
      "16384"
      "16384"
      "16384" )

    LIST(LENGTH OPTION nbTests_tmp)
    MATH(EXPR nbTests "${nbTests_tmp} - 1")

    FOREACH ( test_idx RANGE ${nbTests} )
      LIST ( GET OPTION    ${test_idx} test_option )
      LIST ( GET VAL       ${test_idx} test_val )
      LIST ( GET NAME      ${test_idx} test_name )
      LIST ( GET MESH_SIZE ${test_idx} test_mesh_size )

      FOREACH( NP 1 6 8 )
        add_test( NAME Sphere-optim-${test_name}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${test_option} ${test_val}
          ${CI_DIR}/Sphere/sphere.mesh
          -out ${CI_DIR_RESULTS}/sphere-${test_name}-${NP}-out.mesh
          -m 11000 -mesh-size ${test_mesh_size} ${myargs} )
      ENDFOREACH()
    ENDFOREACH ( )

    ###############################################################################
    #####
    #####        Tests fields interpolation with or without metric
    #####
    ###############################################################################
    add_test( NAME InterpolationFields-withMet-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Interpolation/coarse.mesh
      -out ${CI_DIR_RESULTS}/InterpolationFields-withMet-withFields-4-out.mesh
      -field ${CI_DIR}/Interpolation/sol-fields-coarse.sol
      -sol field3_iso-coarse.sol )

    add_test( NAME InterpolationFields-hsiz-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Interpolation/coarse.mesh
      -out ${CI_DIR_RESULTS}/InterpolationFields-hsiz-withFields-4-out.mesh
      -field ${CI_DIR}/Interpolation/sol-fields-coarse.sol
      -hsiz 0.2 )

    add_test( NAME InterpolationFields-noMet-withFields-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Interpolation/coarse.mesh
      -out ${CI_DIR_RESULTS}/InterpolationFields-noMet-withFields-4-out.mesh
      -field ${CI_DIR}/Interpolation/sol-fields-coarse.sol )

    add_test( NAME InterpolationFields-refinement-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Cube/cube-unit-coarse
      -out ${CI_DIR_RESULTS}/InterpolationFields-refinement-4-out.mesh
      -field ${CI_DIR}/Interpolation/cube-unit-coarse-field.sol )

  ENDIF()


  ###############################################################################
  #####
  #####        Tests that needs the PARMMG LIBRARY
  #####
  ###############################################################################

  SET ( PMMG_LIB_TESTS
    libparmmg_centralized_auto_example0
    libparmmg_centralized_auto_cpp_example0
    libparmmg_centralized_manual_example0_io_0
    libparmmg_centralized_manual_example0_io_1
    #libparmmg_distributed_manual_example0
    )

  SET ( PMMG_LIB_TESTS_MAIN_PATH
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/sequential_IO/automatic_IO/main.c
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/sequential_IO/automatic_IO/main.cpp
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/sequential_IO/manual_IO/main.c
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/sequential_IO/manual_IO/main.c
    #${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/parallel_IO/manual_IO/main.c
    )

  SET ( PMMG_LIB_TESTS_INPUTMESH
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube.mesh
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube.mesh
    ""
    ""
    #${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube.mesh
    )

  SET ( PMMG_LIB_TESTS_INPUTMET
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube-met.sol
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube-met.sol
    ""
    ""
    #${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube-met.sol
    )

  SET ( PMMG_LIB_TESTS_INPUTSOL
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube-solphys.sol
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube-solphys.sol
    ""
    ""
    #""
    )

  SET ( PMMG_LIB_TESTS_OUTPUTMESH
    ${CI_DIR_RESULTS}/io-seq-auto-cube.o.mesh
    ${CI_DIR_RESULTS}/io-seq-auto-cpp-cube.o.mesh
    ${CI_DIR_RESULTS}/io-seq-manual-cube_io_0.o
    ${CI_DIR_RESULTS}/io-seq-manual-cube_io_1.o
    #${CI_DIR_RESULTS}/io-seq-par-cube.o
    )

  SET ( PMMG_LIB_TESTS_OPTIONS
    "-met"
    "-met"
    "0"
    "1"
    #"-met"
    )

  SET ( PMMG_LIB_TESTS_FIELDOPT
    "-field"
    "-field"
    ""
    ""
    #"-met"
    )

  IF ( LIBPARMMG_STATIC )
    SET ( lib_name lib${PROJECT_NAME}_a )
  ELSEIF ( LIBPARMMG_SHARED )
    SET ( lib_name lib${PROJECT_NAME}_so )
  ELSE ()
    MESSAGE(WARNING "You must activate the compilation of the static or"
      " shared ${PROJECT_NAME} library to compile this tests." )
  ENDIF ( )

  #####         Fortran Tests

  IF ( MPI_Fortran_FOUND )
    SET( CMAKE_Fortran_COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}" )
    SET( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} ${MPI_LINK_FLAGS}" )
    SET( FORTRAN_LIBRARIES ${MPI_Fortran_LIBRARIES} )

    LIST ( APPEND PMMG_LIB_TESTS libparmmg_fortran_centralized_auto_example0
      # libparmmg_centralized_manual_example0_io_0
      # libparmmg_centralized_manual_example0_io_1
      # libparmmg_distributed_manual_example0
      )

    LIST ( APPEND PMMG_LIB_TESTS_MAIN_PATH
      ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/sequential_IO/automatic_IO/main.F90
      # ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/sequential_IO/manual_IO/main.c
      # ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/sequential_IO/manual_IO/main.c
      # ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/parallel_IO/manual_IO/main.c
      )

    LIST ( APPEND PMMG_LIB_TESTS_INPUTMESH
      ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube.mesh
      #""
      #""
      #${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube.mesh
      )

    LIST ( APPEND PMMG_LIB_TESTS_INPUTMET
      ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube-met.sol
     # ""
     # ""
     # ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube-met.sol
     )

    LIST ( APPEND PMMG_LIB_TESTS_INPUTSOL
      ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube-solphys.sol
      #""
      #""
      #""
      )

    LIST ( APPEND PMMG_LIB_TESTS_OUTPUTMESH
      ${CI_DIR_RESULTS}/io-seq-auto-cube.o.mesh
      #${CI_DIR_RESULTS}/io-seq-manual-cube_io_0.o
      #${CI_DIR_RESULTS}/io-seq-manual-cube_io_1.o
      #${CI_DIR_RESULTS}/io-seq-par-cube.o
       )

     LIST ( APPEND PMMG_LIB_TESTS_OPTIONS
       "-met"
       #"0"
       #"1"
       #"-met"
       )

     LIST ( APPEND PMMG_LIB_TESTS_FIELDOPT
       "-field"
       #"0"
       #"1"
       #"-met"
       )
  ENDIF ( )

  LIST(LENGTH PMMG_LIB_TESTS nbTests_tmp)
  MATH(EXPR nbTests "${nbTests_tmp} - 1")

  LIST ( APPEND lib_name ${FORTRAN_LIBRARIES})
  LIST ( APPEND lib_name ${MPI_CXX_LIBRARIES})

  FOREACH ( test_idx RANGE ${nbTests} )
    LIST ( GET PMMG_LIB_TESTS            ${test_idx} test_name )
    LIST ( GET PMMG_LIB_TESTS_MAIN_PATH  ${test_idx} main_path )
    LIST ( GET PMMG_LIB_TESTS_INPUTMESH  ${test_idx} input_mesh )
    LIST ( GET PMMG_LIB_TESTS_INPUTMET   ${test_idx} input_met )
    LIST ( GET PMMG_LIB_TESTS_INPUTSOL   ${test_idx} input_sol )
    LIST ( GET PMMG_LIB_TESTS_OUTPUTMESH ${test_idx} output_mesh )
    LIST ( GET PMMG_LIB_TESTS_OPTIONS    ${test_idx} options )
    LIST ( GET PMMG_LIB_TESTS_FIELDOPT   ${test_idx} field )

    ADD_LIBRARY_TEST ( ${test_name} ${main_path} "copy_pmmg_headers" "${lib_name}" )

    FOREACH( NP 1 2 6 )
      ADD_TEST ( NAME ${test_name}-${NP} COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP}
        $<TARGET_FILE:${test_name}>
        ${input_mesh} ${output_mesh} ${options} ${input_met} ${field} ${input_sol})
    ENDFOREACH()

  ENDFOREACH ( )

  IF ( NOT ONLY_LIBRARY_TESTS )

    # Localization test
    SET ( main_path  ${CI_DIR}/WaveSurface/locate.c )
    SET ( input_mesh ${CI_DIR}/WaveSurface/wave.mesh )

    SET ( test_name  WaveSurface_locate_saddle )
    SET ( output_mesh ${CI_DIR_RESULTS}/out_locate_wave_saddle.mesh )
    ADD_LIBRARY_TEST ( ${test_name} ${main_path} "copy_pmmg_headers" "${lib_name}" )
    ADD_TEST ( NAME ${test_name} COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${test_name}> ${input_mesh} ${output_mesh} 0.375 0.375 0.51 217 )

    SET ( test_name  WaveSurface_locate_concave )
    SET ( output_mesh ${CI_DIR_RESULTS}/out_locate_wave_concave.mesh )
    ADD_LIBRARY_TEST ( ${test_name} ${main_path} "copy_pmmg_headers" "${lib_name}" )
    ADD_TEST ( NAME ${test_name} COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${test_name}> ${input_mesh} ${output_mesh} 0.5 0.5 0.755 5934 )

    SET ( test_name  WaveSurface_locate_convex )
    SET ( output_mesh ${CI_DIR_RESULTS}/out_locate_wave_convex.mesh )
    ADD_LIBRARY_TEST ( ${test_name} ${main_path} "copy_pmmg_headers" "${lib_name}" )
    ADD_TEST ( NAME ${test_name} COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${test_name}> ${input_mesh} ${output_mesh} 0.5 0.5 0.74 3888 )

    SET ( test_name  WaveSurface_locate_exhaustive )
    SET ( output_mesh ${CI_DIR_RESULTS}/out_locate_wave_exhaustive.mesh )
    ADD_LIBRARY_TEST ( ${test_name} ${main_path} "copy_pmmg_headers" "${lib_name}" )
    ADD_TEST ( NAME ${test_name} COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${test_name}> ${input_mesh} ${output_mesh} 0.5 0.5 0.74 2494 )

    SET ( test_name  WaveSurface_locate_inexistent )
    SET ( output_mesh ${CI_DIR_RESULTS}/out_locate_wave_inexistent.mesh )
    ADD_LIBRARY_TEST ( ${test_name} ${main_path} "copy_pmmg_headers" "${lib_name}" )
    ADD_TEST ( NAME ${test_name} COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${test_name}> ${input_mesh} ${output_mesh} 0.5 0.5 0.25 3888 )
    SET_PROPERTY( TEST ${test_name} PROPERTY WILL_FAIL ON )

    # Surface interpolation tests
    SET ( input_mesh ${CI_DIR}/WaveSurface/wave.mesh )
    SET ( input_met  ${CI_DIR}/WaveSurface/wave-met.sol )
    SET ( test_name  WaveSurface_interp )

    FOREACH( NP 1 4 )
      add_test( NAME ${test_name}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${input_mesh} -sol ${input_met}
          -out ${CI_DIR_RESULTS}/${test_name}-${NP}-out.mesh
          -niter 1 -nobalance -v 10 )
    ENDFOREACH()


    # Sequential test
    SET ( test_name  LnkdList_unitTest )
    SET ( main_path  ${CI_DIR}/LnkdList_unitTest/main.c )

    ADD_LIBRARY_TEST ( ${test_name} ${main_path} "copy_pmmg_headers" "${lib_name}" )
    ADD_TEST ( NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}> )

    SET ( test_name  API_set_XName )
    SET ( main_path  ${CI_DIR}/API/PMMG_set_XName/main.c )

    ADD_LIBRARY_TEST ( ${test_name} ${main_path} "copy_pmmg_headers" "${lib_name}" )
    ADD_TEST ( NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}> )

    # Distributed API test
    SET ( PMMG_DISTR_LIB_TESTS
      libparmmg_distributed_manual_example0
      libparmmg_distributed_automatic_example0
      )
    SET ( PMMG_DISTR_LIB_TESTS_MAIN_PATH
      ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/parallel_IO/manual_IO/main.c
      ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/parallel_IO/automatic_IO/main.c
      )
    SET ( PMMG_DISTR_LIB_TESTS_INPUTMESH
      ""
      ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube.mesh
      )
    SET ( PMMG_DISTR_LIB_TESTS_INPUTMET
      ""
      ""
      )
    SET ( PMMG_DISTR_LIB_TESTS_INPUTSOL
      ""
      ""
      )
    SET ( PMMG_DISTR_LIB_TESTS_OUTPUTMESH
      ${CI_DIR_RESULTS}/io-par-manual-cube.o
      ${CI_DIR_RESULTS}/io-par-automatic-cube.o
      )

    LIST(LENGTH PMMG_DISTR_LIB_TESTS nbTests_tmp)
    MATH(EXPR nbTests "${nbTests_tmp} - 1")

    FOREACH ( test_idx RANGE ${nbTests} )
      LIST ( GET PMMG_DISTR_LIB_TESTS            ${test_idx} test_name )
      LIST ( GET PMMG_DISTR_LIB_TESTS_MAIN_PATH  ${test_idx} main_path )
      LIST ( GET PMMG_DISTR_LIB_TESTS_INPUTMESH  ${test_idx} input_mesh )
      LIST ( GET PMMG_DISTR_LIB_TESTS_INPUTMET   ${test_idx} input_met )
      LIST ( GET PMMG_DISTR_LIB_TESTS_INPUTSOL   ${test_idx} input_sol )
      LIST ( GET PMMG_DISTR_LIB_TESTS_OUTPUTMESH ${test_idx} output_mesh )

      ADD_LIBRARY_TEST ( ${test_name} ${main_path} "copy_pmmg_headers" "${lib_name}" )

      FOREACH( niter 0 3 )
        FOREACH( API_mode 0 1 )
          FOREACH( NP 1 2 4 )
            ADD_TEST ( NAME ${test_name}_niter_${niter}-API_${API_mode}-${NP} COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP}
              $<TARGET_FILE:${test_name}>
              ${input_mesh} ${output_mesh}_niter_${niter}-API_${API_mode}-${NP} ${niter} ${API_mode} ${input_met} )
          ENDFOREACH()
        ENDFOREACH()
      ENDFOREACH()

    ENDFOREACH()
  ENDIF()

  # Distributed lib test
  SET ( PMMG_DISTR_LIB_TESTS
    libparmmg_distributed_external_example0
    libparmmg_distributed_external_gen_mesh
    )
  SET ( PMMG_DISTR_LIB_TESTS_MAIN_PATH
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/parallel_IO/external_IO/main.c
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/parallel_IO/external_IO/gen_distributedMesh.c
    )
  SET ( PMMG_DISTR_LIB_TESTS_INPUTMESH
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube_in
    ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube
    )
  SET ( PMMG_DISTR_LIB_TESTS_OUTPUTMESH
    ${CI_DIR_RESULTS}/io-par-external-cube.o
    ""
    )

  LIST(LENGTH PMMG_DISTR_LIB_TESTS nbTests_tmp)
  MATH(EXPR nbTests "${nbTests_tmp} - 1")

  FOREACH ( test_idx RANGE ${nbTests} )
    LIST ( GET PMMG_DISTR_LIB_TESTS            ${test_idx} test_name )
    LIST ( GET PMMG_DISTR_LIB_TESTS_MAIN_PATH  ${test_idx} main_path )
    LIST ( GET PMMG_DISTR_LIB_TESTS_INPUTMESH  ${test_idx} input_mesh )
    LIST ( GET PMMG_DISTR_LIB_TESTS_OUTPUTMESH ${test_idx} output_mesh )

    ADD_LIBRARY_TEST ( ${test_name} ${main_path} "copy_pmmg_headers" "${lib_name}" )

    FOREACH( API_mode 0 1 )
      FOREACH( NP 4 )
        ADD_TEST ( NAME ${test_name}_API_${API_mode}-${NP} COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP}
          $<TARGET_FILE:${test_name}>
          ${input_mesh} ${output_mesh}_API_${API_mode}-${NP} ${API_mode} )
      ENDFOREACH()
    ENDFOREACH()

  ENDFOREACH()


ENDIF()
