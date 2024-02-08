IF( BUILD_TESTING )
  include( CTest )

  set( CI_DIR  ${CMAKE_BINARY_DIR}/testparmmg CACHE PATH "path to test meshes repository" )
  set( CI_DIR_RESULTS  ${CMAKE_BINARY_DIR}/TEST_OUTPUTS )
  file( MAKE_DIRECTORY ${CI_DIR_RESULTS} )
  get_filename_component(PARENT_DIR ${CI_DIR} DIRECTORY)


  IF ( NOT ONLY_LIBRARY_TESTS )

    FIND_PACKAGE ( Git )

    IF ( Git_FOUND )

      IF ( NOT EXISTS ${CI_DIR} )
        EXECUTE_PROCESS(
          COMMAND ${GIT_EXECUTABLE} clone https://gitlab.inria.fr/ParMmg/testparmmg.git --filter=blob:none
          WORKING_DIRECTORY ${PARENT_DIR}
          )
      ENDIF()
      EXECUTE_PROCESS(
        COMMAND ${GIT_EXECUTABLE} -C ${CI_DIR} fetch
        COMMAND ${GIT_EXECUTABLE} -C ${CI_DIR} checkout b3fece6cb6afbcd73962c7586aafa211af396e4c
        TIMEOUT 20
        WORKING_DIRECTORY ${CI_DIR}
        #COMMAND_ECHO STDOUT
        )
    ENDIF ( )

    set ( mesh_size 16384 )
    set ( myargs -niter 2 -metis-ratio 82 -v 5 )

    # remesh 2 sets of matching mesh/sol files (which are the output of mmg3d)
    # on 1,2,4,6,8 processors
    foreach( MESH cube-unit-dual_density cube-unit-int_sphere )
      foreach( NP 1 2 4 6 8 )
        add_test( NAME ${MESH}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${CI_DIR}/Cube/${MESH}.meshb
          -out ${CI_DIR_RESULTS}/${MESH}-${NP}-out.mesh
          -m 11000 -mesh-size ${mesh_size} ${myargs})
      endforeach()
    endforeach()

    # remesh a unit cube with two different solution files on 1,2,4,6,8 processors
    foreach( MESH dual_density int_sphere )
      foreach( NP 1 2 4 6 8 )
        add_test( NAME cube-unit-coarse-${MESH}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${CI_DIR}/Cube/cube-unit-coarse.meshb
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
          -mesh-size ${mesh_size} ${myargs} )
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
        ${CI_DIR}/Sphere/sphere.meshb
        -out ${CI_DIR_RESULTS}/sphere-${NP}-out.mesh
        -mesh-size ${mesh_size} ${myargs} )
    endforeach()

    # Option without arguments
    foreach( OPTION "optim" "optimLES" "noinsert" "noswap"  )
      foreach( NP 1 6 8 )
        add_test( NAME Sphere-optim-${OPTION}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          -${OPTION}
          ${CI_DIR}/Sphere/sphere.meshb
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
          ${CI_DIR}/Sphere/sphere.meshb
          -out ${CI_DIR_RESULTS}/sphere-${test_name}-${NP}-out.mesh
          -m 11000 -mesh-size ${test_mesh_size} ${myargs} )
      ENDFOREACH()
    ENDFOREACH ( )

    ### test openbdy mode on 6 procs
    add_test( NAME opnbdy_peninsula-6
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 6 $<TARGET_FILE:${PROJECT_NAME}>
      -opnbdy
      ${CI_DIR}/OpnBdy_peninsula/peninsula.mesh
      -out ${CI_DIR_RESULTS}/opnbdy-peninsula.o.mesh
      )

    add_test( NAME opnbdy_island-6
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 6 $<TARGET_FILE:${PROJECT_NAME}>
      -opnbdy
      ${CI_DIR}/OpnBdy_island/island.mesh
      -out ${CI_DIR_RESULTS}/opnbdy-island.o.mesh
      )

    ###############################################################################
    #####
    #####        Test centralized/distributed I/O (on multidomain and openbdy tests)
    #####
    ###############################################################################

    add_test( NAME opnbdy_island-8
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 8 $<TARGET_FILE:${PROJECT_NAME}>
      -opnbdy -distributed-output
      ${CI_DIR}/OpnBdy_island/island.mesh
      -out ${CI_DIR_RESULTS}/opnbdy-island-distrib.o.mesh
      )
    add_test( NAME opnbdy_island-8-rerun
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 8 $<TARGET_FILE:${PROJECT_NAME}>
      -opnbdy -centralized-output
      ${CI_DIR_RESULTS}/opnbdy-island-distrib.o.mesh
      )
    set_tests_properties(opnbdy_island-8-rerun PROPERTIES DEPENDS opnbdy_island-8 )

    add_test( NAME multidom_wave-8
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 8 $<TARGET_FILE:${PROJECT_NAME}>
      -distributed-output -ar 89 -nobalance
      ${CI_DIR}/WaveSurface/wave.mesh
      -out ${CI_DIR_RESULTS}/multidom-wave-distrib.o.mesh
      ${myargs}
      )
    add_test( NAME multidom_wave-8-rerun
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 8 $<TARGET_FILE:${PROJECT_NAME}>
      -centralized-output -ar 89
      ${CI_DIR_RESULTS}/multidom-wave-distrib.o.mesh
      ${myargs}
      )

    set_tests_properties(multidom_wave-8-rerun PROPERTIES DEPENDS multidom_wave-8 )

    add_test( NAME multidom_wave-8-distrib_parRidge
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 8 $<TARGET_FILE:${PROJECT_NAME}>
      -centralized-output -ar 89
      ${CI_DIR}/WaveSurface_distrib/multidom-wave-distrib.o.mesh
      -out ${CI_DIR_RESULTS}/multidom-wave-distrib_parRidge-out.mesh
      ${myargs}
      )

    # Tests for distributed pvtu output with dots in filename.
    # Replacement of dots by dashes.
    IF ( (NOT VTK_FOUND) OR USE_VTK MATCHES OFF )
      set(OutputVtkErr "VTK library not found.")
    ENDIF ( )

    set(OutputVtkRenameFilename "3D-cube-PvtuOut-2-a-o.pvtu")
    set(OutputVtkRenameWarning  "## WARNING: Filename has been changed.")

    add_test( NAME PvtuOut-RenameOut-2
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/LevelSet/1p_centralized/3D-cube.mesh
      -out ${CI_DIR_RESULTS}/3D-cube-PvtuOut-2.a.o.pvtu)

    set_property(TEST PvtuOut-RenameOut-2
      PROPERTY PASS_REGULAR_EXPRESSION
      "${OutputVtkRenameFilename}.*${OutputVtkRenameWarning};
       ${OutputVtkRenameWarning}.*${OutputVtkRenameFilename}")

    ###############################################################################
    #####
    #####        Tests fields interpolation with or without metric
    #####
    ###############################################################################
    add_test( NAME InterpolationFields-withMet-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Interpolation/coarse.meshb
      -out ${CI_DIR_RESULTS}/InterpolationFields-withMet-withFields-4-out.mesh
      -field ${CI_DIR}/Interpolation/sol-fields-coarse.sol
      -sol field3_iso-coarse.sol
      -mesh-size 60000 ${myargs} )

    add_test( NAME InterpolationFields-hsiz-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Interpolation/coarse.meshb
      -out ${CI_DIR_RESULTS}/InterpolationFields-hsiz-withFields-4-out.mesh
      -field ${CI_DIR}/Interpolation/sol-fields-coarse.sol
      -mesh-size 60000 -hsiz 0.2 ${myargs} )

    add_test( NAME InterpolationFields-noMet-withFields-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Interpolation/coarse.meshb
      -out ${CI_DIR_RESULTS}/InterpolationFields-noMet-withFields-4-out.mesh
      -field ${CI_DIR}/Interpolation/sol-fields-coarse.sol
      -mesh-size 60000 ${myargs} )

    add_test( NAME InterpolationFields-refinement-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Cube/cube-unit-coarse
      -out ${CI_DIR_RESULTS}/InterpolationFields-refinement-4-out.mesh
      -field ${CI_DIR}/Interpolation/cube-unit-coarse-field.sol ${myargs} )

    ###############################################################################
    #####
    #####        Tests distributed surface adaptation
    #####
    ###############################################################################

    # Run the test only if the mesh distribution has succeed
    FOREACH( API_mode 0 1 )
      FOREACH( NP 2 4 8 )

        ADD_TEST( NAME DistribSurf-A319-gen-${API_mode}-${NP}
          COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP}
          $<TARGET_FILE:libparmmg_distributed_external_gen_mesh>
          ${CI_DIR}/A319_gmsh/A319_in_a_box.mesh
          ${CI_DIR_RESULTS}/A319_in_a_box_${API_mode}-${NP}.mesh ${API_mode} )

        ADD_TEST( NAME DistribSurf-A319-adp-${API_mode}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${CI_DIR_RESULTS}/A319_in_a_box_${API_mode}-${NP}.mesh
          -out ${CI_DIR_RESULTS}/DistribSurf-A319-${API_mode}-${NP}-out.mesh
          -optim -v 6 -hmin 20 -hausd 5 -centralized-output
          ${myargs} )

        SET_TESTS_PROPERTIES(DistribSurf-A319-adp-${API_mode}-${NP}
          PROPERTIES DEPENDS DistribSurf-A319-gen-${API_mode}-${NP})


        ADD_TEST( NAME DistribSphere_NOM-gen-${API_mode}-${NP}
          COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP}
          $<TARGET_FILE:libparmmg_distributed_external_gen_mesh>
          ${CI_DIR}/Sphere_NOM/sphere_nom.meshb
          ${CI_DIR_RESULTS}/sphere_nom_${API_mode}-${NP}.mesh ${API_mode} )

        ADD_TEST( NAME DistribSphere_NOM-adp-${API_mode}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${CI_DIR_RESULTS}/sphere_nom_${API_mode}-${NP}.mesh
          -out ${CI_DIR_RESULTS}/DistribSphere_NOM-${API_mode}-${NP}-out.mesh
          -optim -v 6 -hmin 0.5 -hausd 2 -centralized-output
          ${myargs} -niter 3 ) #override previous value of -niter

        SET_TESTS_PROPERTIES(DistribSphere_NOM-adp-${API_mode}-${NP}
          PROPERTIES DEPENDS DistribSphere_NOM-gen-${API_mode}-${NP})


        ADD_TEST( NAME DistribTorus-gen-${API_mode}-${NP}
          COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP}
          $<TARGET_FILE:libparmmg_distributed_external_gen_mesh>
          ${CI_DIR}/Torus/torusholes.mesh
          ${CI_DIR_RESULTS}/torusholes_${API_mode}-${NP}.mesh ${API_mode} )

        ADD_TEST( NAME DistribTorus-adp-${API_mode}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${CI_DIR_RESULTS}/torusholes_${API_mode}-${NP}.mesh
          -out ${CI_DIR_RESULTS}/torusholes-${API_mode}-${NP}-out.mesh
          -v 6 -centralized-output
          ${myargs} -niter 3 ) #override previous value of -niter

        SET_TESTS_PROPERTIES(DistribTorus-adp-${API_mode}-${NP}
          PROPERTIES DEPENDS DistribTorus-gen-${API_mode}-${NP})

      ENDFOREACH()
    ENDFOREACH()

    # Test to verify the patch on update MG_REF tag.
    # This test fail if the tag MG_REF is not updated by PMMG_updateTagRef_node in PMMG_update_analys.
    # See ParMmg PR#103
    add_test( NAME update-ref-tag
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/LevelSet/2p_toygeom/cube-distributed-faces-nomat-1edge.mesh -v 10 -hsiz 0.1
      -out ${CI_DIR_RESULTS}/update-ref-tag.o.mesh)

    ###############################################################################
    #####
    #####        Test isovalue mode - ls discretization
    #####
    ###############################################################################
    #--------------------------------
    #--- CENTRALIZED INPUT (CenIn)
    #--------------------------------
    # Tests of ls discretization for centralized mesh input
    foreach( NP 1 2 4 8 )
      add_test( NAME ls-CenIn-${NP}
        COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
        ${CI_DIR}/LevelSet/1p_centralized/3D-cube.mesh
        -ls 0.0
        -sol ${CI_DIR}/LevelSet/1p_centralized/3D-cube-ls.sol
        -out ${CI_DIR_RESULTS}/3D-cube-ls-CenIn-${NP}.o.mesh)
    endforeach()

    # Check that the ls file is correctly opened with or without the ls value given
    set(lsOpenFile "3D-cube-ls.sol OPENED")
    set(lsOpenFileDefault "3D-cube.sol  NOT FOUND. USE DEFAULT METRIC.")

    # Test of opening ls file when ls val is given
    foreach( NP 1 2)
      add_test( NAME ls-arg-option-openlsfile-lsval-${NP}
        COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
        ${CI_DIR}/LevelSet/1p_centralized/3D-cube.mesh
        -ls 0.0
        -sol ${CI_DIR}/LevelSet/1p_centralized/3D-cube-ls.sol
        -out ${CI_DIR_RESULTS}/ls-arg-option-openlsfile-lsval-${NP}.o.mesh)
      set_property(TEST ls-arg-option-openlsfile-lsval-${NP}
        PROPERTY PASS_REGULAR_EXPRESSION "${lsOpenFile}")
    endforeach()

    # Test of opening ls file when ls val is not given
    foreach( NP 1 2)
      add_test( NAME ls-arg-option-openlsfile-nolsval-${NP}
        COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
        ${CI_DIR}/LevelSet/1p_centralized/3D-cube.mesh
        -ls
        -sol ${CI_DIR}/LevelSet/1p_centralized/3D-cube-ls.sol
        -out ${CI_DIR_RESULTS}/ls-arg-option-openlsfile-nolsval-${NP}.o.mesh)
      set_property(TEST ls-arg-option-openlsfile-nolsval-${NP}
        PROPERTY PASS_REGULAR_EXPRESSION "${lsOpenFile}")
    endforeach()

    # Test of opening ls file with a default name when ls val is given
    foreach( NP 1 2)
      add_test( NAME ls-arg-option-openlsfiledefault-lsval-${NP}
        COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
        ${CI_DIR}/LevelSet/1p_centralized/3D-cube.mesh
        -ls 0.0
        -out ${CI_DIR_RESULTS}/ls-arg-option-openlsfiledefault-lsval-${NP}.o.mesh)
      set_property(TEST ls-arg-option-openlsfiledefault-lsval-${NP}
        PROPERTY PASS_REGULAR_EXPRESSION "${lsOpenFileDefault}")
    endforeach()

    # Test of opening ls file with a default name when ls val is not given
    foreach( NP 1 2)
      add_test( NAME ls-arg-option-openlsfiledefault-nolsval-${NP}
        COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
        ${CI_DIR}/LevelSet/1p_centralized/3D-cube.mesh
        -ls
        -out ${CI_DIR_RESULTS}/ls-arg-option-openlsfiledefault-nolsval-${NP}.o.mesh)
      set_property(TEST ls-arg-option-openlsfiledefault-nolsval-${NP}
        PROPERTY PASS_REGULAR_EXPRESSION "${lsOpenFileDefault}")
    endforeach()

    # Tests for ls + met for centralized mesh input
    foreach( NP 1 2 4 8 )
    add_test( NAME ls-CenIn-met-${NP}
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/LevelSet/1p_centralized/3D-cube.mesh
      -ls 0.0
      -sol ${CI_DIR}/LevelSet/1p_centralized/3D-cube-ls.sol
      -met ${CI_DIR}/LevelSet/1p_centralized/3D-cube-metric.sol
      -out ${CI_DIR_RESULTS}/3D-cube-ls-CenIn-met-${NP}.o.mesh)
    endforeach()

    # Tests of distributed pvtu output when ls mode
    foreach( NP 1 2 4 8 )
      add_test( NAME ls-CenIn-DisOut-${NP}
        COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
        ${CI_DIR}/LevelSet/1p_centralized/3D-cube.mesh
        -ls 0.0
        -sol ${CI_DIR}/LevelSet/1p_centralized/3D-cube-ls.sol
        -out ${CI_DIR_RESULTS}/3D-cube-ls-CenIn-DisOut-${NP}-out.pvtu)

      IF ( (NOT VTK_FOUND) OR USE_VTK MATCHES OFF )
        set_property(TEST ls-CenIn-DisOut-${NP}
          PROPERTY PASS_REGULAR_EXPRESSION "${OutputVtkErr}")
      ENDIF ( )

    endforeach()

    #--------------------------------
    #--- DISTRIBUTED INPUT (DisIn)
    #--------------------------------
    add_test( NAME ls-DisIn-ReadLs-2
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/LevelSet/2p_distributed/3D-cube.mesh -v 10
      -ls 0.01
      -sol ${CI_DIR}/LevelSet/2p_distributed/3D-cube-ls.sol
      -out ${CI_DIR_RESULTS}/ls-DisIn-ReadLs-2.o.mesh)
    set(lsReadFile "3D-cube-ls.0.sol OPENED")
    set_property(TEST ls-DisIn-ReadLs-2
      PROPERTY PASS_REGULAR_EXPRESSION "${lsReadFile}")

    # Test Medit and hdf5 distributed inputs, with npartin < npart or npartin ==
    # npart with mesh only or mesh+metric.

    ## Medit distributed with npart = 2 and  npartin = 1, only mesh and hdf5 output using .h5 ext
    add_test( NAME Medit-DisIn-MeshOnly-2
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Parallel_IO/Medit/1p/cube-unit-coarse.mesh -v 5
      -out ${CI_DIR_RESULTS}/Medit-DisIn-MeshOnly-2.o.h5)

    ## Medit distributed with npart = 2 and  npartin = 1, mesh+met and hdf5 output using .xdmf ext
    add_test( NAME Medit-DisIn-MeshAndMet-2
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:${PROJECT_NAME}>
      -in ${CI_DIR}/Parallel_IO/Medit/1p/cube-unit-coarse-with-met -v 5
      -out ${CI_DIR_RESULTS}/Medit-DisIn-MeshAndMet-2.o.xdmf)

    ## Medit distributed with npart = 4 and  npartin = 4, only mesh .h5 ext
    add_test( NAME Medit-DisIn-MeshOnly-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      -in ${CI_DIR}/Parallel_IO/Medit/4p/cube-unit-coarse.mesh -v 5
      ${CI_DIR_RESULTS}/Medit-DisIn-MeshOnly-4.o.h5)

    ## Medit distributed with npart = 6 and  npartin = 4, only mesh .xdmf ext
    add_test( NAME Medit-DisIn-MeshOnly-6
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 6 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Parallel_IO/Medit/4p/cube-unit-coarse -v 5
      ${CI_DIR_RESULTS}/Medit-DisIn-MeshOnly-6.o.xdmf)

    ## hdf5 distributed with npart = 2 and  npartin = 1, only mesh and h5 output
    add_test( NAME hdf5-DisIn-MeshOnly-2
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Parallel_IO/hdf5/1p/cube-unit-coarse.h5 -v 5
      -out ${CI_DIR_RESULTS}/hdf5-DisIn-MeshOnly-2.o.h5)

    ## hdf5 distributed with npart = 2 and  npartin = 1, mesh+met and xdmf (h5) output
    add_test( NAME hdf5-DisIn-MeshAndMet-2
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Parallel_IO/hdf5/1p/cube-unit-coarse-with-met.h5 -v 5
      -out ${CI_DIR_RESULTS}/hdf5-DisIn-MeshAndMet-2.o.xdmf)

    ## hdf5 distributed with npart = 8 and  npartin = 4, mesh+met and h5 output
    add_test( NAME hdf5-DisIn-MeshAndMet-8
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 8 $<TARGET_FILE:${PROJECT_NAME}>
      -in ${CI_DIR}/Parallel_IO/hdf5/4p/cube-unit-coarse-with-met.h5 -v 5
      ${CI_DIR_RESULTS}/hdf5-DisIn-MeshAndMet-8.o.h5)

    ## hdf5 distributed with npart = 8 and  npartin = 4, mesh only and medit centralized output
    add_test( NAME hdf5-DisIn-MeshOnly-8
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 8 $<TARGET_FILE:${PROJECT_NAME}>
      -in ${CI_DIR}/Parallel_IO/hdf5/4p/cube-unit-coarse.h5 -v 5 -centralized-output
      -out ${CI_DIR_RESULTS}/hdf5-DisIn-MeshOnly-8.o.mesh)

    ## hdf5 distributed with npart = 4 and  npartin = 4, mesh+met and h5 output
    add_test( NAME hdf5-DisIn-MeshAndMet-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      -in ${CI_DIR}/Parallel_IO/hdf5/4p/cube-unit-coarse-with-met.h5 -v 5
      ${CI_DIR_RESULTS}/hdf5-DisIn-MeshAndMet-8.o.h5)

    ## hdf5 distributed with npart = 4 and  npartin = 4, mesh only and medit centralized output
    add_test( NAME hdf5-DisIn-MeshOnly-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      -in ${CI_DIR}/Parallel_IO/hdf5/4p/cube-unit-coarse.h5 -v 5 -centralized-output
      -out ${CI_DIR_RESULTS}/hdf5-DisIn-MeshOnly-8.o.mesh)


    IF ( (NOT HDF5_FOUND) OR USE_HDF5 MATCHES OFF )
      SET(expr "HDF5 library not found")
      SET_PROPERTY(
        TEST Medit-DisIn-MeshOnly-2 Medit-DisIn-MeshAndMet-2 Medit-DisIn-MeshOnly-4
        Medit-DisIn-MeshOnly-6 hdf5-DisIn-MeshOnly-2 hdf5-DisIn-MeshAndMet-2
        hdf5-DisIn-MeshAndMet-8  hdf5-DisIn-MeshOnly-8
        hdf5-DisIn-MeshAndMet-4  hdf5-DisIn-MeshOnly-4
        PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
    ENDIF ( )


    ###############################################################################
    #####
    #####        Test with fields input and output
    #####
    ###############################################################################
    #--------------------------------
    #--- DISTRIBUTED INPUT (DisIn)
    #--------------------------------
    # Test to read  distributed input  fields in Medit format
    # and  to write distributed output fields in VTK   format
    add_test( NAME fields-DisIn-DisOutVTK-2
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/LevelSet/2p_distributed/3D-cube.mesh -v 10
      -field ${CI_DIR}/LevelSet/2p_distributed/3D-cube-fields.sol
      -out ${CI_DIR_RESULTS}/3D-cube-fields-DisIn-DisOutVTK-2-out.pvtu)

    set(InputDistributedFields "3D-cube-fields.0.sol OPENED")
    set(OutputVtkFields "Writing mesh, metric and fields.")

    set_property(TEST fields-DisIn-DisOutVTK-2
      PROPERTY PASS_REGULAR_EXPRESSION
      "${InputDistributedFields}.*${OutputVtkFields};
      ${OutputVtkFields}.*${InputDistributedFields}")

    # Test to write distributed output fields and metric in Medit format
    add_test( NAME fields-DisIn-DisOutMesh-2
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/LevelSet/2p_distributed/3D-cube.mesh
      -field ${CI_DIR}/LevelSet/2p_distributed/3D-cube-fields.sol
      -out ${CI_DIR_RESULTS}/3D-cube-fields-DisIn-DisOutMesh-2.o.mesh)

    set(OutputFieldsName "3D-cube-fields.o.0.sol OPENED.")
    set(OutputMetricName "3D-cube-fields-DisIn-DisOutMesh-2.o.0.sol OPENED.")
    set_property(TEST fields-DisIn-DisOutMesh-2
      PROPERTY PASS_REGULAR_EXPRESSION
      "${OutputFieldsName}.*${OutputMetricName};${OutputMetricName}.*${OutputFieldsName}")

    # Test saving of solution fields on 4 procs at hdf5 format
    add_test( NAME hdf5-CenIn-DisOutHdf5-4
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${PROJECT_NAME}>
      ${CI_DIR}/Interpolation/coarse.meshb -v 5
      -out ${CI_DIR_RESULTS}/hdf5-CenIn-DisOutHdf5-4.o.h5)

    IF ( (NOT HDF5_FOUND) OR USE_HDF5 MATCHES OFF )
      SET(expr "HDF5 library not found")
      SET_PROPERTY(
        TEST hdf5-CenIn-DisOutHdf5-4
        PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
    ENDIF ( )

  ENDIF()

  ###############################################################################
  #####
  #####        Tests that needs the PARMMG LIBRARY
  #####
  ###############################################################################
  SET ( LIB_TESTS OFF )

  IF ( LIBPARMMG_STATIC )
    SET ( lib_name lib${PROJECT_NAME}_a )
    SET ( LIB_TESTS ON )
  ELSEIF ( LIBPARMMG_SHARED )
    SET ( LIB_TESTS ON )
    SET ( lib_name lib${PROJECT_NAME}_so )
  ENDIF ( )

  if ( LIB_TESTS )
    #----------------- library examples in the ParMmg repo
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

    # Usable on 2 procs only
    ADD_LIBRARY_TEST ( libparmmg_distributed_manual_opnbdy
      ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/parallel_IO/manual_IO/opnbdy.c
      "copy_pmmg_headers" "${lib_name}"
      )

    ADD_TEST ( NAME libparmmg_distributed_manual_opnbdy-2
      COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2
      $<TARGET_FILE:libparmmg_distributed_manual_opnbdy>
      ${CI_DIR_RESULTS}/io-par-manual-opnbdy.o.mesh )

    #####         Fortran Tests
    IF ( MPI_Fortran_FOUND )
      SET( CMAKE_Fortran_COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}" )
      SET( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} ${MPI_LINK_FLAGS}" )
      SET( FORTRAN_LIBRARIES ${MPI_Fortran_LIBRARIES} )

      LIST ( APPEND PMMG_LIB_TESTS libparmmg_fortran_centralized_auto_example0
        # libparmmg_fortran_centralized_manual_example0_io_0
        # libparmmg_fortran_centralized_manual_example0_io_1
        # libparmmg_fortran_distributed_manual_example0
        )

      LIST ( APPEND PMMG_LIB_TESTS_MAIN_PATH
        ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/sequential_IO/automatic_IO/main.F90
        # ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/sequential_IO/manual_IO/main.F90
        # ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/sequential_IO/manual_IO/main.F90
        # ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/parallel_IO/manual_IO/main.F90
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

    # Distributed lib test
    ADD_LIBRARY_TEST ( libparmmg_distributed_external_gen_mesh
      ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/parallel_IO/external_IO/gen_distributedMesh.c
      "copy_pmmg_headers" "${lib_name}" )
    ADD_LIBRARY_TEST ( libparmmg_distributed_external_example0
      ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/parallel_IO/external_IO/main.c
      "copy_pmmg_headers" "${lib_name}" )

    # Run the test only if the mesh distribution has succeed
    FOREACH( API_mode 0 1 )
      FOREACH( NP 4 )
        ADD_TEST ( NAME  libparmmg_distributed_external_gen_mesh_API_${API_mode}-${NP}
          COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP}
          $<TARGET_FILE:libparmmg_distributed_external_gen_mesh>
          ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example0/cube
          ${CI_DIR_RESULTS}/cube-distrib_API_${API_mode}-${NP}.mesh ${API_mode} )

        ADD_TEST ( NAME  libparmmg_distributed_external_example0_API_${API_mode}-${NP}
          COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP}
          $<TARGET_FILE:libparmmg_distributed_external_example0>
          ${CI_DIR_RESULTS}/cube-distrib_API_${API_mode}-${NP}
          ${CI_DIR_RESULTS}/io-par-external-cube_API_${API_mode}-${NP} ${API_mode} )

        set_tests_properties(libparmmg_distributed_external_example0_API_${API_mode}-${NP}
          PROPERTIES DEPENDS libparmmg_distributed_external_gen_mesh_API_${API_mode}-${NP})
      ENDFOREACH()
    ENDFOREACH()

    # Distributed analysis
    SET( test_name libparmmg_distributed_example1 )
    ADD_LIBRARY_TEST ( ${test_name}
      ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example1/main.c
      "copy_pmmg_headers" "${lib_name}" )

    FOREACH( API_mode 0 )
      FOREACH( NP 4 )
        ADD_TEST ( NAME  libparmmg_distributed_example1_wave${API_mode}-${NP}
          COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP}
          $<TARGET_FILE:${test_name}>
          ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example1/wave
          ${CI_DIR_RESULTS}/io-par_wave_${API_mode}-${NP} ${API_mode} )
      ENDFOREACH()
    ENDFOREACH()

    FOREACH( API_mode 0 )
      FOREACH( NP 4 )
        ADD_TEST ( NAME  libparmmg_distributed_example1_brickLS${API_mode}-${NP}
          COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP}
          $<TARGET_FILE:${test_name}>
          ${PROJECT_SOURCE_DIR}/libexamples/adaptation_example1/brickLS
          ${CI_DIR_RESULTS}/io-par_wave_${API_mode}-${NP} ${API_mode} )
      ENDFOREACH()
    ENDFOREACH()


    #----------------- Tests using the library in the testparmmg repos
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

      FOREACH( NP 1 )
        add_test( NAME ${test_name}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${input_mesh} -sol ${input_met}
          -out ${CI_DIR_RESULTS}/${test_name}-${NP}-out.mesh
          -niter 3 -nobalance -v 10 )
      ENDFOREACH()

      SET ( input_mesh ${CI_DIR}/Tennis/tennis.meshb )
      SET ( input_met  ${CI_DIR}/Tennis/tennis.sol )
      SET ( test_name  TennisSurf_interp )

      FOREACH( NP 1 4 )
        add_test( NAME ${test_name}-${NP}
          COMMAND ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${PROJECT_NAME}>
          ${input_mesh} -sol ${input_met}
          -out ${CI_DIR_RESULTS}/${test_name}-${NP}-out.mesh
          -niter 3 -nobalance -v 10 )
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

      # 2 procs tests
      ADD_LIBRARY_TEST ( opnbdy-along-interface
        ${CI_DIR}/Parallel_IO/manual_IO/opnbdy-along-interface.c "copy_pmmg_headers"
        "${lib_name}" )

      ADD_TEST ( NAME opnbdy-along-interface-FaceComm
        COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2
        $<TARGET_FILE:opnbdy-along-interface>
        ${CI_DIR_RESULTS}/opnbdy-along-interface-FaceComm 0)

      ADD_TEST ( NAME opnbdy-along-interface-NodeComm
        COMMAND  ${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} 2
        $<TARGET_FILE:opnbdy-along-interface>
        ${CI_DIR_RESULTS}/opnbdy-along-interface-FaceComm 1)

    ENDIF(  NOT ONLY_LIBRARY_TESTS )
  ENDIF ( LIB_TESTS )
ENDIF()
