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
    COMMAND genheader_pmmg ${out_dir}/${out_file} ${in_dir}/${in_file} ${include_dir}
    ${PROJECT_SOURCE_DIR}/scripts/genfort.pl
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS genheader_pmmg ${in_dir}/${in_file}
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
MACRO ( COPY_HEADER
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
#####         Copy source and automatically generated header files to another
#####         place and create the associated target
#####
###############################################################################
MACRO ( COPY_HEADERS_AND_CREATE_TARGET
    source_dir binary_dir include_dir )

  ADD_CUSTOM_TARGET(pmmgtypes_header ALL
    DEPENDS
    ${source_dir}/libparmmgtypes.h )

  ADD_CUSTOM_TARGET(pmmgversion_header ALL
    DEPENDS
    ${binary_dir}/pmmgversion.h )

  ADD_CUSTOM_TARGET(pmmg_header ALL
    DEPENDS
    ${source_dir}/libparmmg.h )

  COPY_HEADER (
    ${source_dir} libparmmgtypes.h
    ${include_dir} libparmmgtypes.h
    pmmgtypes_header copy_libpmmgtypes )

  COPY_HEADER (
    ${binary_dir} pmmgversion.h
    ${include_dir} pmmgversion.h
    pmmgversion_header copy_pmmgversion )

  COPY_HEADER (
    ${source_dir} libparmmg.h
    ${include_dir} libparmmg.h
    pmmg_header copy_libpmmg
    )

  COPY_HEADER (
    ${binary_dir} libparmmgtypesf.h
    ${include_dir} libparmmgtypesf.h
    pmmgtypes_fortran_header copy_libpmmgtypesf )

  COPY_HEADER (
    ${binary_dir} libparmmgf.h
    ${include_dir} libparmmgf.h
    pmmg_fortran_header copy_libpmmgf
    )

  SET ( tgt_list  copy_libpmmgf copy_libpmmgtypesf
    copy_libpmmg copy_libpmmgtypes copy_pmmgversion)

  IF (NOT WIN32 OR MINGW)
    COPY_HEADER (
      ${binary_dir} git_log_pmmg.h
      ${include_dir} git_log_pmmg.h
      GenerateGitHash copy_pmmggithash )

    LIST ( APPEND tgt_list copy_pmmggithash)
  ENDIF ()

  ADD_CUSTOM_TARGET(copy_pmmg_headers ALL
    DEPENDS ${tgt_list})

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
  add_library( ParMmg::${target_name} ALIAS ${target_name} )

  IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
    INCLUDE_DIRECTORIES ( ${target_name} PRIVATE
      ${COMMON_BINARY_DIR} ${COMMON_SOURCE_DIR} ${CMAKE_BINARY_DIR}/include )
  ELSE ( )
    target_include_directories( ${target_name} PUBLIC
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/>
      $<BUILD_INTERFACE:${COMMON_SOURCE_DIR}>
      $<BUILD_INTERFACE:${COMMON_BINARY_DIR}>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> )
  ENDIF ( )

  SET_TARGET_PROPERTIES ( ${target_name} PROPERTIES
    OUTPUT_NAME ${output_name}
    VERSION ${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}
    SOVERSION ${CMAKE_RELEASE_VERSION_MAJOR} )

  SET_PROPERTY(TARGET ${target_name} PROPERTY C_STANDARD 99)

  TARGET_LINK_LIBRARIES ( ${target_name} ${LIBRARIES} )

  install(TARGETS ${target_name} EXPORT ParMmgTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    Component lib
    )
ENDMACRO ( )

###############################################################################
#####
#####         Add an executable to build and needed include dir, set its
#####         postfix, add link dependencies and the install rule
#####
###############################################################################

MACRO ( ADD_AND_INSTALL_EXECUTABLE
    exec_name lib_files main_file )

  IF ( NOT TARGET lib${exec_name}_a AND NOT TARGET lib${exec_name}_so )
    ADD_EXECUTABLE ( ${exec_name} ${lib_files} ${main_file} )
  ELSE ( )
    ADD_EXECUTABLE ( ${exec_name} ${main_file})

    SET_PROPERTY(TARGET ${exec_name} PROPERTY C_STANDARD 99)

    IF ( NOT TARGET lib${exec_name}_a )
      TARGET_LINK_LIBRARIES(${exec_name} lib${exec_name}_so)
    ELSE ( )
      TARGET_LINK_LIBRARIES(${exec_name} lib${exec_name}_a)
    ENDIF ( )

  ENDIF ( )

  IF ( WIN32 AND NOT MINGW AND SCOTCH_FOUND )
    my_add_link_flags ( ${exec_name} "/SAFESEH:NO")
  ENDIF ( )

 IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
   INCLUDE_DIRECTORIES ( ${exec_name} PUBLIC
     ${COMMON_BINARY_DIR} ${COMMON_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include )
 ELSE ( )
   target_include_directories( ${exec_name} PUBLIC
     $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/>
     $<BUILD_INTERFACE:${COMMON_SOURCE_DIR}>
     $<BUILD_INTERFACE:${COMMON_BINARY_DIR}>
     $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> )
 ENDIF ( )

  TARGET_LINK_LIBRARIES ( ${exec_name} ${LIBRARIES}  )

  INSTALL(TARGETS ${exec_name}
    EXPORT ParMmgTargets RUNTIME DESTINATION bin COMPONENT appli)

  ADD_TARGET_POSTFIX(${exec_name})

ENDMACRO ( )

###############################################################################
#####
#####         Add a target postfix depending on the build type
#####
###############################################################################

MACRO ( ADD_TARGET_POSTFIX target_name )
  IF ( CMAKE_BUILD_TYPE MATCHES "Debug" )
    # in debug mode we name the executable mmgs_debug
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES DEBUG_POSTFIX _debug)
  ELSEIF ( CMAKE_BUILD_TYPE MATCHES "Release" )
    # in Release mode we name the executable mmgs_O3
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES RELEASE_POSTFIX _O3)
  ELSEIF ( CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo" )
    # in RelWithDebInfo mode we name the executable mmgs_O3d
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES RELWITHDEBINFO_POSTFIX _O3d)
  ELSEIF ( CMAKE_BUILD_TYPE MATCHES "MinSizeRel" )
    # in MinSizeRel mode we name the executable mmgs_O3
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES MINSIZEREL_POSTFIX _Os)
  ENDIF ( )
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

  IF ( WIN32 AND ((NOT MINGW) AND SCOTCH_FOUND) )
    MY_ADD_LINK_FLAGS ( ${target_name} "/SAFESEH:NO" )
  ENDIF ( )

  TARGET_LINK_LIBRARIES ( ${target_name}  ${lib_name} )
  INSTALL(TARGETS ${target_name} RUNTIME DESTINATION bin COMPONENT appli )

ENDMACRO ( )
