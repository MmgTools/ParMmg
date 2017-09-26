IF( BUILD_TESTING )
  include( CTest )

  set( CI_DIR  ${CMAKE_BINARY_DIR}/Tests )
  file( MAKE_DIRECTORY ${CI_DIR} )
  set( CI_DIR_RESULTS  ${CI_DIR}/output )
  file( MAKE_DIRECTORY ${CI_DIR_RESULTS} )
  set( CI_DIR_INPUTS  "../../TestParMmg" CACHE PATH "path to test meshes repository" )

  foreach( MESH cube_unit-dual_density cube_unit-int_sphere )
    foreach( NP 1 2 4 6 8 )
      add_test( NAME ${MESH}-${NP}
                COMMAND ${MPIEXEC} -np ${NP} $<TARGET_FILE:${PROJECT_NAME}>
                ${CI_DIR_INPUTS}/Cube/${MESH}.mesh
                -out ${CI_DIR_RESULTS}/${MESH}-${NP}-out.mesh )
    endforeach()
  endforeach()
ENDIF()
