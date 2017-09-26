IF( BUILD_TESTING )
  enable_testing()

  set( CI_DIR  ${CMAKE_BINARY_DIR}/Tests )
  file( MAKE_DIRECTORY ${CI_DIR} )
  set( CI_DIR_RESULTS  ${CI_DIR}/output )
  file( MAKE_DIRECTORY ${CI_DIR_RESULTS} )

  #NIKOS #include( CTest )
  #  add_test(NAME <name> [CONFIGURATIONS [Debug|Release|...]]
  #            [WORKING_DIRECTORY dir]
  #            COMMAND <command> [arg1 [arg2 ...]])
    #              ${CMAKE_BINARY_DIR dir}
    add_test( NAME Cube_HalfSmall_halfBig-2P
      COMMAND ${MPIEXEC} -np 2 $<TARGET_FILE:${PROJECT_NAME}>
                       ../../TestParMmg/Cube/cube_unit-dual_density.mesh
                       -out ${CI_DIR_RESULTS}/cube_unit-dual_density-out2P.mesh )
    add_test( NAME Cube_HalfSmall_halfBig-4P
      COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:${PROJECT_NAME}>
                       ../../TestParMmg/Cube/cube_unit-dual_density.mesh
                       -out ${CI_DIR_RESULTS}/cube_unit-dual_density-out4P.mesh )
    add_test( NAME Cube_HalfSmall_halfBig-8P
      COMMAND ${MPIEXEC} -np 8 $<TARGET_FILE:${PROJECT_NAME}>
                       ../../TestParMmg/Cube/cube_unit-dual_density.mesh
                       -out ${CI_DIR_RESULTS}/cube_unit-dual_density-out8P.mesh )
ENDIF()
