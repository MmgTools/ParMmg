###############################################################################
#####
#####         Generation of a fortran header file
#####
###############################################################################

MACRO ( GENERATE_FORTRAN_HEADER name
    in_dir in_file include_dir out_dir out_file
    )
  # Wrap add_custom_command into add_custom target to remove dpendencies from
  # the custom command and thus allow parallel build.
  ADD_CUSTOM_COMMAND (
    OUTPUT ${out_dir}/${out_file}
    COMMAND genheader ${out_dir}/${out_file} ${in_dir}/${in_file} ${include_dir}
    ${PROJECT_SOURCE_DIR}/scripts/genfort.pl
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS genheader ${in_dir}/${in_file}
    ${PROJECT_SOURCE_DIR}/scripts/genfort.pl
    COMMENT "Generating Fortran header for ${name}"
    )

  ADD_CUSTOM_TARGET (
    ${name}_fortran_header
    ALL
    DEPENDS ${out_dir}/${out_file}
    )

ENDMACRO ( )

###############################################################################
#####
#####         Copy an automatically generated header file to another place
#####
###############################################################################
MACRO ( COPY_FORTRAN_HEADER
    in_dir in_file out_dir out_file
    file_dependencies
    target_name
    )
  # Wrap add_custom_command into add_custom target to remove dpendencies from
  # the custom command and thus allow parallel build.
  ADD_CUSTOM_COMMAND (
    OUTPUT  ${out_dir}/${out_file}
    COMMAND ${CMAKE_COMMAND} -E copy  ${in_dir}/${in_file} ${out_dir}/${out_file}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS ${file_dependencies} ${in_dir}/${in_file}
    COMMENT "Copying ${in_dir}/${in_file} in ${out_dir}/${out_file}"
    )

  ADD_CUSTOM_TARGET ( ${target_name} ALL
    DEPENDS ${out_dir}/${out_file} )
  ADD_DEPENDENCIES ( ${target_name} ${file_dependencies} )

ENDMACRO ( )

###############################################################################
#####
#####         Copy an automatically generated header file to another place
#####         and create the associated target
#####
###############################################################################
MACRO ( COPY_FORTRAN_HEADER_AND_CREATE_TARGET
    binary_dir include_dir )

  COPY_FORTRAN_HEADER (
    ${binary_dir} libparmmgtypesf.h
    ${include_dir} libparmmgtypesf.h
    pmmg_fortran_header copy_libpmmgtypesf )

  COPY_FORTRAN_HEADER (
    ${binary_dir}
    libparmmgf.h ${include_dir}
    libparmmgf.h
    pmmg_fortran_header copy_libpmmgf
    )

  ADD_CUSTOM_TARGET(copy_pmmg_headers ALL
    DEPENDS
    copy_libpmmgf copy_libpmmgtypesf
    ${include_dir}/libparmmg.h
    ${include_dir}/libparmmgtypes.h )

ENDMACRO ( )

###############################################################################
#####
#####         Add a library to build and needed include dir, set its
#####         properties, add link dependencies and the install rule
#####
###############################################################################

MACRO ( ADD_AND_INSTALL_LIBRARY
    target_name target_type sources output_name )

  ADD_LIBRARY ( ${target_name} ${target_type} ${sources} )

  IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
    INCLUDE_DIRECTORIES ( ${target_name} PRIVATE
      ${COMMON_BINARY_DIR} ${COMMON_SOURCE_DIR} ${CMAKE_BINARY_DIR}/include )
  ELSE ( )
    TARGET_INCLUDE_DIRECTORIES ( ${target_name} PRIVATE
      ${COMMON_BINARY_DIR} ${COMMON_SOURCE_DIR} ${CMAKE_BINARY_DIR}/include )
  ENDIF ( )

  SET_TARGET_PROPERTIES ( ${target_name}
    PROPERTIES OUTPUT_NAME ${output_name} )

  SET_PROPERTY(TARGET ${target_name} PROPERTY C_STANDARD 99)

  TARGET_LINK_LIBRARIES ( ${target_name} ${LIBRARIES} )

  INSTALL ( TARGETS ${target_name}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    COMPONENT lib)

ENDMACRO ( )


###############################################################################
#####
#####         Add a library test
#####
###############################################################################

MACRO ( ADD_LIBRARY_TEST target_name main_path target_dependency lib_name )
  ADD_EXECUTABLE ( ${target_name} ${main_path} )
  ADD_DEPENDENCIES( ${target_name} ${target_dependency} )

  IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
    INCLUDE_DIRECTORIES ( ${target_name} PUBLIC ${CMAKE_BINARY_DIR}/include )
  ELSE ( )
    TARGET_INCLUDE_DIRECTORIES ( ${target_name} PUBLIC ${CMAKE_BINARY_DIR}/include )
  ENDIF ( )

  IF ( WIN32 AND ((NOT MINGW) AND USE_SCOTCH) )
    MY_ADD_LINK_FLAGS ( ${target_name} "/SAFESEH:NO" )
  ENDIF ( )

  TARGET_LINK_LIBRARIES ( ${target_name}  ${lib_name} )
  INSTALL(TARGETS ${target_name} RUNTIME DESTINATION bin COMPONENT appli )

ENDMACRO ( )
