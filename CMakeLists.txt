CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)

INCLUDE(CMakeDependentOption)
INCLUDE(cmake/modules/macros.cmake)

PROJECT (parmmg)

###############################################################################
#####
#####         Release version and date
#####
###############################################################################
SET (CMAKE_RELEASE_VERSION_MAJOR "1")
SET (CMAKE_RELEASE_VERSION_MINOR "0")
SET (CMAKE_RELEASE_VERSION_PATCH "0")
SET (CMAKE_RELEASE_DATE "Feb 15, 2019")

SET (CMAKE_RELEASE_VERSION
  "${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}")

###############################################################################
#####
#####         Main CMake compilation variables
#####
###############################################################################

SET(PMMG_BINARY_DIR      ${CMAKE_BINARY_DIR}/src/parmmg)
SET(PMMG_SHRT_INCLUDE    parmmg )
SET(PMMG_INCLUDE         ${CMAKE_BINARY_DIR}/include/${PMMG_SHRT_INCLUDE} )

FILE(MAKE_DIRECTORY ${PMMG_BINARY_DIR})

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
FIND_PACKAGE( MPI COMPONENTS C REQUIRED )

SET( CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}" )
SET( CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS}" )

INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )

SET( CMAKE_C_FLAGS "-DUSE_MPI ${CMAKE_C_FLAGS}" )
MESSAGE( STATUS "Compilation with mpi" )
SET( LIBRARIES ${MPI_C_LIBRARIES} ${LIBRARIES} )

EXECUTE_PROCESS ( COMMAND  ${MPIEXEC} --help mapping
  RESULT_VARIABLE CMD_ERROR
  OUTPUT_VARIABLE CMD_OUTPUT )

STRING(REGEX REPLACE "\"" " " CMD_OUT ${CMD_OUTPUT})

IF ( ${CMD_OUT} MATCHES "oversubscribe" )
  MESSAGE ( STATUS "OK" )
  SET ( MPI_ARGS "-oversubscribe" )
ENDIF()

############################################################################
#####
#####        MMG (for mesh data structure)
#####
############################################################################
SET(MMG_DIR "" CACHE PATH "Installation directory for mmg")
SET(MMG_BUILDDIR "" CACHE PATH "Build directory for mmg")

FIND_PACKAGE(Mmg)

IF(NOT MMG_FOUND )
  MESSAGE ( FATAL_ERROR "ERROR: The installation directory for mmg is required:"
    "(see https://github.com/MmgTools/mmg and  download the branch develop)."
    "If you have already installed Mmg and want to use it, "
    "please set the CMake variable or environment variable MMG_DIR "
    "to your mmg directory and the CMake variable or environment variable"
    " MMG_BUILDDIR "
    "to your mmg build directory.")
ELSE()
  INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIRS})

  MESSAGE(STATUS
    "Compilation with Mmg: ${MMG_LIBRARIES}")
  SET( LIBRARIES ${MMG_LIBRARIES} ${LIBRARIES})

  # Additionnal directories to access the Mmg sources
  INCLUDE_DIRECTORIES(${MMG_BUILDDIR_INTERNAL}/src/common)
  IF ( MMG_DIR )
    INCLUDE_DIRECTORIES(${MMG_DIR}/src/mmg2d)
    INCLUDE_DIRECTORIES(${MMG_DIR}/src/mmg3d)
    INCLUDE_DIRECTORIES(${MMG_DIR}/src/common)
  ELSE ( )
    MESSAGE ( FATAL_ERROR "ERROR: The source directory for mmg is required:"
      "(see https://github.com/MmgTools/mmg and  download the branch develop)."
      "If you have already installed Mmg and want to use it, "
      "please set the CMake variable or environment variable MMG_DIR "
      "to your mmg directory.")
  ENDIF ( )
ENDIF ( )

############################################################################
#####
#####        Mmg PATTERN?
#####
############################################################################

FILE(READ ${MMG_BUILDDIR_INTERNAL}/CMakeFiles/libmmg_a.dir/flags.make TMPTXT)
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
  OPTION ( USE_METIS "Detected METIS Tool" ON)
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

IF ( NOT ${METIS_INTSIZE} EQUAL 4 )
  MESSAGE ( WARNING "Metis integer uses ${METIS_INTSIZE} bytes."
    " Only int32 are supported." )
ENDIF ( )

############################################################################
#####
#####        ParMetis
#####
############################################################################
# add parmetis library?
SET(PARMETIS_DIR "" CACHE PATH "Installation directory for parmetis")
FIND_PACKAGE(Parmetis)

IF( PARMETIS_FOUND)
  OPTION ( USE_PARMETIS "Detected PARMETIS TOOL" ON)
  IF( USE_PARMETIS )
    # Set flags for building test program
    INCLUDE_DIRECTORIES(${PARMETIS_INCLUDE_DIRS})

    SET(CMAKE_C_FLAGS "-DUSE_PARMETIS ${CMAKE_C_FLAGS}")

    MESSAGE(STATUS
      "Compilation with parmetis: ${PARMETIS_LIBRARIES}")
    SET( LIBRARIES ${PARMETIS_LIBRARIES} ${LIBRARIES})
  ENDIF()
ELSE ( )
#  MESSAGE(STATUS ${PARMETIS_INCLUDE_DIR} ${PARMETIS_LIBRARY} )
#
#  MESSAGE ( ERROR " parmetis library not found: "
#    "If you have already installed ParMetis and want to use it, "
#    "please set the CMake variable or environment variable PARMETIS_DIR "
#    "to your parmetis directory.")
ENDIF ( )

############################################################################
#####
#####         Scotch
#####
############################################################################
# Find SCOTCH library?
INCLUDE(cmake/modules/FindScotch.cmake)

SET(SCOTCH_DIR "" CACHE PATH "Installation directory for scotch")
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
FIND_PACKAGE(Scotch)

CMAKE_DEPENDENT_OPTION ( USE_SCOTCH "Use SCOTCH TOOL for renumbering" ON
  "SCOTCH_FOUND" OFF)

IF( USE_SCOTCH )
  # Scotch renumbering doesn't work for now
  #  SET( LIBRARIES ${SCOTCH_LIBRARIES} ${LIBRARIES})
ENDIF()

############################################################################
#####
#####         Fortran header: libparmmgf.h
#####
############################################################################

IF (NOT WIN32)

  ADD_EXECUTABLE(genheader ${CMAKE_SOURCE_DIR}/scripts/genheader.c)
  TARGET_LINK_LIBRARIES ( genheader ${LIBRARIES}  )

  GENERATE_FORTRAN_HEADER ( pmmgtypes
    ${PMMG_SOURCE_DIR} libparmmgtypes.h
    ${PMMG_SHRT_INCLUDE}
    ${PMMG_BINARY_DIR} libparmmgtypesf.h
    )

  GENERATE_FORTRAN_HEADER ( pmmg
    ${PMMG_SOURCE_DIR} libparmmg.h
    ${PMMG_SHRT_INCLUDE}
    ${PMMG_BINARY_DIR} libparmmgf.h
    )

ENDIF ()

###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################
# Header files
INCLUDE_DIRECTORIES(${PMMG_SOURCE_DIR})

# Source files
FILE(
  GLOB
  pmmg_library_files
  ${PMMG_SOURCE_DIR}/*.c
  )
LIST(REMOVE_ITEM pmmg_library_files
  ${PMMG_SOURCE_DIR}/parmmg.c
  ${REMOVE_FILE})
FILE(
  GLOB
  pmmg_main_file
  ${PMMG_SOURCE_DIR}/parmmg.c
  )

IF (NOT WIN32)
  FIND_LIBRARY(M_LIB m)
  SET( LIBRARIES ${LIBRARIES} ${M_LIB})
ENDIF()

############################################################################
#####
##### RPATH for MacOSX
#####
############################################################################

SET(CMAKE_MACOSX_RPATH 1)

############################################################################
#####
#####         Compile parmmg libraries
#####
############################################################################
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

OPTION ( LIBPARMMG_STATIC "Compile static library" ON)
IF ( LIBPARMMG_STATIC )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}_a STATIC
    "${pmmg_library_files}" ${PROJECT_NAME} )
ENDIF()

# Compile shared library
OPTION ( LIBPARMMG_SHARED "Compile dynamic library" OFF)
IF ( LIBPARMMG_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}_so SHARED
    "${pmmg_library_files}" ${PROJECT_NAME} )
ENDIF()

# parmmg header files needed for library
SET( pmmg_headers
  ${PMMG_SOURCE_DIR}/libparmmg.h
  ${PMMG_BINARY_DIR}/libparmmgf.h
  ${PMMG_SOURCE_DIR}/libparmmgtypes.h
  ${PMMG_BINARY_DIR}/libparmmgtypesf.h
  )

# Install header files in /usr/local or equivalent
INSTALL(FILES ${pmmg_headers} DESTINATION include/parmmg COMPONENT headers )

COPY_FORTRAN_HEADER_AND_CREATE_TARGET ( ${PMMG_BINARY_DIR} ${PMMG_INCLUDE} )

# Copy header files in project directory at configuration step
# (generated file don't exists yet or are outdated)
FILE(INSTALL  ${pmmg_headers} DESTINATION ${PMMG_INCLUDE}
  PATTERN "libparmmg*f.h"  EXCLUDE)

###############################################################################
#####
#####         Compile PMMG executable
#####
###############################################################################
ADD_AND_INSTALL_EXECUTABLE ( ${PROJECT_NAME}
  "${pmmg_library_files}" ${pmmg_main_file} )

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
