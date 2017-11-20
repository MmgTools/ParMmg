CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)

PROJECT (parmmg)

###############################################################################
#####
#####         Release version and date
#####
###############################################################################
SET (CMAKE_RELEASE_VERSION_MAJOR "0")
SET (CMAKE_RELEASE_VERSION_MINOR "0")
SET (CMAKE_RELEASE_VERSION_PATCH "0")
SET (CMAKE_RELEASE_DATE "July 2, 2015")

SET (CMAKE_RELEASE_VERSION
  "${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}")

###############################################################################
#####
#####         Main CMake compilation variables
#####
###############################################################################

# To see flags and options of compilation
#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

# Executable path
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
SET(PMMG_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)


# Find "exotic" compilers
IF (${CMAKE_C_COMPILER} MATCHES Clang OR ${CMAKE_C_COMPILER} MATCHES clang)
  # using clang
  SET(CMAKE_COMPILER_IS_CLANG TRUE)
ELSEIF(${CMAKE_C_COMPILER} MATCHES Icc OR ${CMAKE_C_COMPILER} MATCHES icc)
  # using icc
  SET(CMAKE_COMPILER_IS_INTEL TRUE)
ENDIF()

# Hide some options
MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES
  CMAKE_OSX_DEPLOYMENT_TARGET
  CMAKE_OSX_SYSROOT)

IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS " -Wno-char-subscripts ${CMAKE_C_FLAGS}")
  IF(APPLE)
    # Add flags to the compiler to work on old mac
    ADD_DEFINITIONS( -mmacosx-version-min=10.4 -arch x86_64)

    # To avoid pbs with binary files...
    SET(CMAKE_EXE_LINKER_FLAGS "-arch x86_64 ${CMAKE_EXE_LINKER_FLAGS} ")

    # Determine if the processor supports 64bit execution
    EXECUTE_PROCESS(
      COMMAND sysctl hw.cpu64bit_capable
      ERROR_QUIET
      OUTPUT_VARIABLE 64_CMD
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    STRING(REGEX REPLACE "^hw.cpu64bit_capable: (.*)" "\\1" 64_BIT "${64_CMD}")
    #  ELSEIF(UNIX)# UNIX must be after APPLE becaus APPLE is UNIX too
  ENDIF()
ENDIF()

IF(NOT CMAKE_COMPILER_IS_CLANG)
  # Compiler options for profiling... but not possible with clang
  OPTION ( PROFILING "Enable/Disable PROFILING" OFF )
  IF(PROFILING)
    ADD_DEFINITIONS(-pg)
    SET(CMAKE_EXE_LINKER_FLAGS "-pg ${CMAKE_EXE_LINKER_FLAGS}")
  ENDIF(PROFILING)
ENDIF(NOT CMAKE_COMPILER_IS_CLANG)

###############################################################################
#####
#####         Choose executable target to compile
#####
###############################################################################
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  MESSAGE(STATUS "Setting build type to 'Release' as none was specified.")
  seT(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
  SET(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} )
ENDIF()

############################################################################
#####
#####         MPI
#####
############################################################################
# compile parallel version
FIND_PACKAGE( MPI REQUIRED )

SET( CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}" )
SET( CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS}" )
INCLUDE_DIRECTORIES( MPI_INCLUDE_PATH )

SET( CMAKE_C_FLAGS "-DUSE_MPI ${CMAKE_C_FLAGS}" )
MESSAGE( STATUS "Compilation with mpi" )
SET( LIBRARIES ${MPI_LIBRARIES} ${LIBRARIES} )
SET( CMAKE_C_COMPILER  ${MPI_C_COMPILER} )

############################################################################
#####
#####        MMG (for mesh data structure)
#####
############################################################################
SET(MMG_DIR "" CACHE PATH "Installation directory for mmg")

INCLUDE(cmake/modules/FindMmg.cmake)

IF(NOT MMG_DIR OR NOT MMG_FOUND )
  MESSAGE ( FATAL_ERROR "ERROR: The installation directory for mmg is required:"
    "(see https://github.com/MmgTools/mmg and  download the branch develop)."
    "If you have already installed Mmg and want to use it, "
    "please set the CMake variable or environment variable MMG_DIR "
    "to your mmg directory.")
ELSE()
  # Set flags for building test program
  INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIR})
  INCLUDE_DIRECTORIES(${MMG_BUILD_DIR}/src/common)
  INCLUDE_DIRECTORIES(${MMG_DIR}/src/mmg2d)
  INCLUDE_DIRECTORIES(${MMG_DIR}/src/mmg3d)
  INCLUDE_DIRECTORIES(${MMG_DIR}/src/common)

  #SET(CMAKE_REQUIRED_INCLUDES ${MMG_INCLUDE_DIR})

  SET( LIBRARIES ${MMG_LIBRARY} ${LIBRARIES})


ENDIF ( )

############################################################################
#####
#####        Mmg PATTERN?
#####
############################################################################

FILE(READ ${MMG_BUILD_DIR}/CMakeFiles/libmmg_a.dir/flags.make TMPTXT)
STRING(FIND "${TMPTXT}" "PATTERN" matchres)

IF(${matchres} EQUAL -1)
  SET ( PATTERN OFF )
ELSE ( )
  SET ( PATTERN ON )
  SET(CMAKE_C_FLAGS "-DPATTERN ${CMAKE_C_FLAGS}")
ENDIF()

############################################################################
#####
#####        Metis
#####
############################################################################
# add metis library?
SET(METIS_DIR "" CACHE PATH "Installation directory for metis")

FIND_PACKAGE(Metis)

IF( METIS_FOUND)
  OPTION ( USE_METIS "Use METIS TOOL for partioning" ON)
  IF( USE_METIS )
    # Set flags for building test program
    INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIRS})

    SET(CMAKE_C_FLAGS "-DUSE_METIS ${CMAKE_C_FLAGS}")

    MESSAGE(STATUS
      "Compilation with metis: ${METIS_LIBRARIES}")
    SET( LIBRARIES ${METIS_LIBRARIES} ${LIBRARIES})
  ENDIF()
ELSE ( )
  MESSAGE(STATUS ${METIS_INCLUDE_DIR} ${METIS_LIBRARIES} )

  MESSAGE ( ERROR " Metis library not found: "
    "If you have already installed Metis and want to use it, "
    "please set the CMake variable or environment variable METIS_DIR "
    "to your metis directory.")
ENDIF ( )

############################################################################
#####
#####        ParMetis
#####
############################################################################
# add metis library?
SET(PARMETIS_DIR "" CACHE PATH "Installation directory for parmetis")
FIND_PACKAGE(Parmetis)

IF( PARMETIS_FOUND)
  OPTION ( USE_PARMETIS "Use PARMETIS TOOL for partioning" ON)
  IF( USE_PARMETIS )
    # Set flags for building test program
    INCLUDE_DIRECTORIES(${PARMETIS_INCLUDE_DIRS})

    SET(CMAKE_C_FLAGS "-DUSE_PARMETIS ${CMAKE_C_FLAGS}")

    MESSAGE(STATUS
      "Compilation with parmetis: ${PARMETIS_LIBRARIES}")
    SET( LIBRARIES ${PARMETIS_LIBRARIES} ${LIBRARIES})
  ENDIF()
ELSE ( )
  MESSAGE(STATUS ${PARMETIS_INCLUDE_DIR} ${PARMETIS_LIBRARY} )

  MESSAGE ( ERROR " parmetis library not found: "
    "If you have already installed ParMetis and want to use it, "
    "please set the CMake variable or environment variable PARMETIS_DIR "
    "to your metis directory.")
ENDIF ( )

############################################################################
#####
#####         Fortran header: libparmmgf.h
#####
############################################################################

IF (NOT WIN32)

  ADD_EXECUTABLE(genheader ${CMAKE_SOURCE_DIR}/scripts/genheader.c)

  ADD_CUSTOM_COMMAND(OUTPUT ${PMMG_SOURCE_DIR}/libparmmgf.h
    COMMAND genheader ${PMMG_SOURCE_DIR}/libparmmgf.h
    ${PMMG_SOURCE_DIR}/libparmmg.h ${CMAKE_SOURCE_DIR}/scripts/genfort.pl
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS genheader ${PMMG_SOURCE_DIR}/libparmmg.h
    ${CMAKE_SOURCE_DIR}/scripts/genfort.pl
    COMMENT "Generating Fortran header"
    )
ENDIF ()

###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################
# Header files
INCLUDE_DIRECTORIES(${PMMG_SOURCE_DIR})
INCLUDE_DIRECTORIES(${MMG_BUILD_DIR})

# Source files
FILE(
  GLOB
  source_files
  ${PMMG_SOURCE_DIR}/*.c   ${PMMG_SOURCE_DIR}/*.h
  )
LIST(REMOVE_ITEM source_files
  ${PMMG_SOURCE_DIR}/parmmg.c
  ${REMOVE_FILE})
FILE(
  GLOB
  main_file
  ${PMMG_SOURCE_DIR}/parmmg.c
  ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}.c
  ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}f.c
  )
FILE(
  GLOB
  lib_file
  ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}.c
  ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}f.c
  )
IF (NOT WIN32)
  FIND_LIBRARY(M_LIB m)
  SET( LIBRARIES ${LIBRARIES} ${M_LIB})
ENDIF()

############################################################################
#####
#####Set the full RPATH to find libraries independently from LD_LIBRARY_PATH
#####
############################################################################
SET(CMAKE_MACOSX_RPATH TRUE)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system
# directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
  "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

############################################################################
#####
#####         Compile parmmg libraries
#####
############################################################################
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

OPTION ( LIBPARMMG_STATIC "Compile static library" ON)
IF ( LIBPARMMG_STATIC )
  ADD_LIBRARY(${PROJECT_NAME}_a  STATIC
    ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}f.h
    ${source_files} ${lib_file} )
  SET_TARGET_PROPERTIES(${PROJECT_NAME}_a PROPERTIES OUTPUT_NAME
    ${PROJECT_NAME})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME}_a ${LIBRARIES})
  INSTALL(FILES ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}f.h DESTINATION
    ${CMAKE_BINARY_DIR}/include)
  INSTALL(FILES ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}.h DESTINATION
    ${CMAKE_BINARY_DIR}/include)
  INSTALL(TARGETS ${PROJECT_NAME}_a
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
ENDIF()

# Compile shared library
OPTION ( LIBPARMMG_SHARED "Compile dynamic library" OFF)
IF ( LIBPARMMG_SHARED )
  ADD_LIBRARY(${PROJECT_NAME}_so SHARED
    ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}f.h
    ${source_files} ${lib_file})
  SET_TARGET_PROPERTIES(${PROJECT_NAME}_so PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME}_so ${LIBRARIES})
  INSTALL(FILES ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}f.h
    DESTINATION ${CMAKE_BINARY_DIR}/include)
  INSTALL(FILES ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}.h
    DESTINATION ${CMAKE_BINARY_DIR}/include)
  INSTALL(TARGETS ${PROJECT_NAME}_so
    ARCHIVE DESTINATION
    LIBRARY DESTINATION lib)
ENDIF()

############################################################################
#####
#####         Compile program to test library
#####
############################################################################

###############################################################################
#####
#####         Compile PMMG executable
#####
###############################################################################
ADD_EXECUTABLE(${PROJECT_NAME}
  ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}f.h
  ${PMMG_SOURCE_DIR}/lib${PROJECT_NAME}.h
  ${source_files} ${main_file})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${LIBRARIES})
INSTALL(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)

# in debug mode we name the executable parmmg_debug
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX _debug)
# in Release mode we name the executable parmmg_O3
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES RELEASE_POSTFIX _O3)
# in RelWithDebInfo mode we name the executable parmmg_d
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES RELWITHDEBINFO_POSTFIX _O)
# in MinSizeRel mode we name the executable parmmg_Os
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES MINSIZEREL_POSTFIX _Os)

###############################################################################
#####
#####         Create API Documentation
#####
###############################################################################
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)

  # PARMMG Documentation
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile @ONLY)

  ADD_CUSTOM_TARGET(fmg_doc
    COMMAND ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile
    # COMMAND ${CMAKE_COMMAND} -E chdir
    # ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/latex make
    # COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/latex/refman.pdf
    # ${CMAKE_CURRENT_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/
    COMMENT "Generating PARMMG API documentation with Doxygen.
    Open the doc/doxygen/html/index.html file to see it." VERBATIM
    )

  ADD_CUSTOM_TARGET(doc
    DEPENDS fmg_doc
    COMMENT "Generating PARMMG API documentation with Doxygen.
     Open the doc/doxygen/html/index.html file to see it" VERBATIM
    )
ENDIF ( DOXYGEN_FOUND )

OPTION ( BUILD_TESTING "Enable / Disable tests" OFF )
INCLUDE( ${CMAKE_SOURCE_DIR}/cmake/testing/pmmg_tests.cmake )
