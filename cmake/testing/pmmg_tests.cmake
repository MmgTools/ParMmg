IF( BUILD_TESTING )
  include( CTest )

  set( CI_DIR  ${CMAKE_BINARY_DIR}/Tests )
  file( MAKE_DIRECTORY ${CI_DIR} )
  set( CI_DIR_RESULTS  ${CI_DIR}/output )
  file( MAKE_DIRECTORY ${CI_DIR_RESULTS} )
  set( CI_DIR_INPUTS  "../../testparmmg" CACHE PATH "path to test meshes repository" )

  # remesh 2 sets of matching mesh/sol files (which are the output of mmg3d)
  # on 1,2,4,6,8 processors
  foreach( MESH cube-unit-dual_density cube-unit-int_sphere )
    foreach( NP 1 2 4 6 8 )
      add_test( NAME ${MESH}-${NP}
                COMMAND ${MPIEXEC} -np ${NP} $<TARGET_FILE:${PROJECT_NAME}>
                ${CI_DIR_INPUTS}/Cube/${MESH}.mesh
                -out ${CI_DIR_RESULTS}/${MESH}-${NP}-out.mesh
                -m 11000 )
    endforeach()
  endforeach()

  # remesh a unit cube with two different solution files on 1,2,4,6,8 processors
  foreach( MESH dual_density int_sphere )
    foreach( NP 1 2 4 6 8 )
      add_test( NAME cube-unit-coarse-${MESH}-${NP}
                COMMAND ${MPIEXEC} -np ${NP} $<TARGET_FILE:${PROJECT_NAME}>
                ${CI_DIR_INPUTS}/Cube/cube-unit-coarse.mesh
                -sol ${CI_DIR_INPUTS}/Cube/cube-unit-coarse-${MESH}.sol
                -out ${CI_DIR_RESULTS}/${MESH}-${NP}-out.mesh )
    endforeach()
  endforeach()

  # remesh a non constant anisotropic test case: a torus with a planar shock
  # on 1,2,4,6,8 processors
  foreach( TYPE anisotropic-test )
    foreach( NP 1 2 4 6 8 )
	    add_test( NAME ${TYPE}-torus-with-planar-shock-${NP}
                COMMAND ${MPIEXEC} -np ${NP} $<TARGET_FILE:${PROJECT_NAME}>
                ${CI_DIR_INPUTS}/Torus/torusholes.mesh
                -sol ${CI_DIR_INPUTS}/Torus/torusholes.sol
                -out ${CI_DIR_RESULTS}/${TYPE}-torus-with-planar-shock-${NP}-out.mesh )
    endforeach()
  endforeach()

  ###############################################################################
  #####
  #####        Tests that needs the PARMMG LIBRARY
  #####
  ###############################################################################

  SET ( PMMG_LIB_TESTS LnkdList_unitTest )
  SET ( PMMG_LIB_TESTS_MAIN_PATH ${CI_DIR_INPUTS}/LnkdList_unitTest/main.c )

  IF ( LIBPARMMG_STATIC )
    SET ( lib_name lib${PROJECT_NAME}_a )
  ELSEIF ( LIBPARMMG_SHARED )
    SET ( lib_name lib${PROJECT_NAME}_so )
  ELSE ()
    MESSAGE(WARNING "You must activate the compilation of the static or"
      " shared ${PROJECT_NAME} library to compile this tests." )
  ENDIF ( )

  LIST(LENGTH PMMG_LIB_TESTS nbTests_tmp)
  MATH(EXPR nbTests "${nbTests_tmp} - 1")

  FOREACH ( test_idx RANGE ${nbTests} )
    LIST ( GET PMMG_LIB_TESTS           ${test_idx} test_name )
    LIST ( GET PMMG_LIB_TESTS_MAIN_PATH ${test_idx} main_path )

    ADD_LIBRARY_TEST ( ${test_name} ${main_path} copy_pmmg_headers ${lib_name} )

    ADD_TEST ( NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}> )
  ENDFOREACH ( )

ENDIF()
