CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)

INCLUDE(CMakeDependentOption)
INCLUDE(cmake/modules/macros.cmake)

PROJECT (parmmg)

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

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

SET(PMMG_BINARY_DIR      ${PROJECT_BINARY_DIR}/src/parmmg)
SET(PMMG_SHRT_INCLUDE    parmmg )
SET(PMMG_INCLUDE         ${PROJECT_BINARY_DIR}/include/${PMMG_SHRT_INCLUDE} )

FILE(MAKE_DIRECTORY ${PMMG_BINARY_DIR})


# To see flags and options of compilation
#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

# Executable path
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(PMMG_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)


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

# Only allow building of static OR shared libraries
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
FUNCTION (INVERT_BOOL OUTVAR INVAR)
  IF(${INVAR})
    SET(${OUTVAR} OFF PARENT_SCOPE)
  ELSE()
    SET(${OUTVAR} ON PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

OPTION ( LIBPARMMG_STATIC "Compile static library" ON)
OPTION ( LIBPARMMG_SHARED "Compile dynamic library" OFF)
INVERT_BOOL("LIBPARMMG_STATIC" ${LIBPARMMG_SHARED})

############################################################################
#####
#####         MPI
#####
############################################################################
# compile parallel version
ENABLE_LANGUAGE ( Fortran )
FIND_PACKAGE( MPI COMPONENTS C CXX OPTIONAL_COMPONENTS Fortran )

SET( CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}" )
SET( CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS}" )

SET( CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}" )
SET( CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}" )

INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )

# openmpi is bugged (deadlocks) : warn the user
EXECUTE_PROCESS ( COMMAND grep "OMPI_MAJOR_VERSION"  "${MPI_C_INCLUDE_PATH}/mpi.h"
  RESULT_VARIABLE CMD_ERROR
  OUTPUT_VARIABLE CMD_OUTPUT )

IF ( ${CMD_ERROR} MATCHES 0 )
  MESSAGE(WARNING "Possible deadlocks with open-mpi (see https://github.com/open-mpi/ompi/issues/6568 )...")
ENDIF()

SET( CMAKE_C_FLAGS "-DUSE_MPI ${CMAKE_C_FLAGS}" )
SET( CMAKE_CXX_FLAGS "-DUSE_MPI ${CMAKE_CXX_FLAGS}" )
MESSAGE( STATUS "Compilation with mpi" )
SET( LIBRARIES ${MPI_C_LIBRARIES} ${LIBRARIES} )
SET( LIBRARIES ${MPI_CXX_LIBRARIES} ${LIBRARIES} )

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
#####         Scotch
#####
############################################################################
# Find SCOTCH library?
SET(SCOTCH_DIR "" CACHE PATH "Installation directory for scotch")
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
FIND_PACKAGE(Scotch)

IF ( NOT SCOTCH_FOUND )
  MESSAGE ( WARNING "Scotch library not found:"
    "Using scotch reduce the execution time of mmg3d "
    "(see https://gforge.inria.fr/frs/?group_id=248 to download it)."
    "If you have already installed Scotch and want to use it, "
    "please set the CMake variable or environment variable SCOTCH_DIR "
    "to your scotch directory.")
ENDIF ( )

############################################################################
#####
#####         VTK (to parse (p)vtp/(p)vtu files )
#####
############################################################################
FIND_PACKAGE(VTK)

###############################################################################
#####
#####         Add dependent options
#####
###############################################################################

# add Scotch library?
CMAKE_DEPENDENT_OPTION ( USE_SCOTCH "Use SCOTCH TOOL for renumbering" ON
  "SCOTCH_FOUND" OFF)

IF( USE_SCOTCH )
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_SCOTCH")
 MESSAGE(STATUS
    "Compilation with scotch: ${SCOTCH_LIBRARIES}")
  SET( LIBRARIES ${LIBRARIES} ${SCOTCH_LIBRARIES})
ENDIF()

# add the VTK library ?
CMAKE_DEPENDENT_OPTION ( USE_VTK "Use VTK I/O" ON
  "VTK_FOUND" OFF)

IF ( USE_VTK )
  ENABLE_LANGUAGE ( CXX )
  ADD_DEFINITIONS(-DUSE_VTK)
  MESSAGE ( STATUS "Compilation with VTK: add vtp and vtu I/O." )
  INCLUDE ( ${VTK_USE_FILE} )
  SET( LIBRARIES  ${LIBRARIES} "-lstdc++" ${VTK_LIBRARIES} )
ENDIF ( )

############################################################################
#####
#####        MMG (for mesh data structure)
#####
############################################################################
INCLUDE ( ExternalProject )

# forward scotch configuration
SET ( SCOTCH_CFG -DUSE_SCOTCH=${USE_SCOTCH} -DSCOTCH_DIR=${SCOTCH_DIR}
  -DSCOTCH_scotch.h_DIRS=${SCOTCH_scotch.h_DIRS}
  -DSCOTCH_scotch_LIBRARY=${SCOTCH_scotch_LIBRARY}
  -DSCOTCH_scotcherrexit_LIBRARY=${SCOTCH_scotcherrexit_LIBRARY} )

# forward VTK config
SET ( VTK_CFG -DUSE_VTK=${USE_VTK} -DVTK_DIR=${VTK_DIR} )

# forward compiler (for c++ compatibility)
SET ( COMPILER_CFG -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} )

OPTION ( DOWNLOAD_MMG "Download and build automatically Mmg." ON )

IF ( DOWNLOAD_MMG )
  UNSET(MMG_DIR CACHE)
  UNSET(MMG_BUILDDIR CACHE)

  EXTERNALPROJECT_ADD ( Mmg
    GIT_REPOSITORY https://github.com/MmgTools/mmg.git
    GIT_TAG 9bd9ac675c7fef52df46e2cabdd4e9119adad425
    INSTALL_COMMAND  echo "Mmg installation DISABLED"
    CMAKE_ARGS ${MMG_ARGS} -DUSE_ELAS=OFF ${COMPILER_CFG}
    ${SCOTCH_CFG} ${VTK_CFG}
    -DBUILD=MMG -DBUILD_SHARED_LIBS=${LIBPARMMG_SHARED} )

  EXTERNALPROJECT_GET_PROPERTY ( Mmg BINARY_DIR )
  SET ( MMG_BINARY_DIR ${BINARY_DIR} )
  EXTERNALPROJECT_GET_PROPERTY ( Mmg SOURCE_DIR )
  SET ( MMG_SOURCE_DIR ${SOURCE_DIR} )

  EXTERNALPROJECT_ADD_STEP ( Mmg reconfigure
    COMMAND ${CMAKE_COMMAND} -E remove ${MMG_BINARY_DIR}/CMakeCache.txt
    DEPENDEES update
    DEPENDERS configure )

  INCLUDE_DIRECTORIES(${MMG_BINARY_DIR}/include)

  IF( LIBPARMMG_STATIC )
    SET(MMG3D_LIBRARY ${MMG_BINARY_DIR}/lib/libmmg3d${CMAKE_STATIC_LIBRARY_SUFFIX})
  ELSE()
    SET(MMG3D_LIBRARY ${MMG_BINARY_DIR}/lib/libmmg3d${CMAKE_SHARED_LIBRARY_SUFFIX})
  ENDIF()
  MESSAGE(STATUS
    "Compilation with Mmg: ${MMG3D_LIBRARY}")
  SET( LIBRARIES ${MMG3D_LIBRARY} ${LIBRARIES})

  # Additionnal directories to access the Mmg sources
  INCLUDE_DIRECTORIES(${MMG_BINARY_DIR}/src/common)
  INCLUDE_DIRECTORIES(${MMG_SOURCE_DIR}/src/mmg3d)
  INCLUDE_DIRECTORIES(${MMG_SOURCE_DIR}/src/common)

ELSE ( )
  MESSAGE ( STATUS "Manual installation of Mmg: please, specify the MMG_DIR and MMG_BUILDDIR CMake variables" )

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
ENDIF()

############################################################################
#####
#####        Metis / ParMetis
#####
############################################################################
OPTION ( DOWNLOAD_METIS "Download and build automatically Metis/ParMetis." ON )

IF ( DOWNLOAD_METIS )
  UNSET ( METIS_DIR CACHE )

  # Metis
  EXTERNALPROJECT_ADD ( Metis
    URL ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/metis-5.1.0.tar.gz
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ${CMAKE_MAKE_PROGRAM} config prefix=../../../Metis-install
    BUILD_IN_SOURCE 1
    BUILD_COMMAND unset MFLAGS && unset MAKEFLAGS && unset MAKELEVEL && ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND unset MFLAGS && unset MAKEFLAGS && unset MAKELEVEL && ${CMAKE_MAKE_PROGRAM} install
    )

  EXTERNALPROJECT_GET_PROPERTY(Metis SOURCE_DIR)
  SET ( METIS_INSTALL_DIR ${SOURCE_DIR}/../Metis-install)

  INCLUDE_DIRECTORIES(${METIS_INSTALL_DIR}/include)
  SET ( METIS_LIBRARIES ${METIS_INSTALL_DIR}/lib/libmetis.a )

  SET(CMAKE_C_FLAGS "-DUSE_METIS ${CMAKE_C_FLAGS}")

  MESSAGE(STATUS
    "Compilation with metis: ${METIS_LIBRARIES}")
  SET( LIBRARIES ${METIS_LIBRARIES} ${LIBRARIES})

ELSE ( )
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

ENDIF ( )

############################################################################
#####
#####        ParMetis
#####
############################################################################

# add parmetis library?
#SET(PARMETIS_DIR "" CACHE PATH "Installation directory for parmetis")
#FIND_PACKAGE(Parmetis)
#
#IF( PARMETIS_FOUND)
#  OPTION ( USE_PARMETIS "Detected PARMETIS TOOL" ON)
#  IF( USE_PARMETIS )
#    # Set flags for building test program
#    INCLUDE_DIRECTORIES(${PARMETIS_INCLUDE_DIRS})
#
#    SET(CMAKE_C_FLAGS "-DUSE_PARMETIS ${CMAKE_C_FLAGS}")
#
#    MESSAGE(STATUS
#      "Compilation with parmetis: ${PARMETIS_LIBRARIES}")
#    SET( LIBRARIES ${PARMETIS_LIBRARIES} ${LIBRARIES})
#  ENDIF()
#ELSE ( )
#  #  MESSAGE(STATUS ${PARMETIS_INCLUDE_DIR} ${PARMETIS_LIBRARY} )
#  #
#  #  MESSAGE ( ERROR " parmetis library not found: "
#  #    "If you have already installed ParMetis and want to use it, "
#  #    "please set the CMake variable or environment variable PARMETIS_DIR "
#  #    "to your parmetis directory.")
#ENDIF ( )

############################################################################
#####
#####         Fortran header: libparmmgf.h
#####
############################################################################

IF (NOT WIN32)

  ADD_EXECUTABLE(genheader_pmmg ${PROJECT_SOURCE_DIR}/scripts/genheader.c)

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
  ${PMMG_SOURCE_DIR}/*.cpp
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

  IF ( DOWNLOAD_MMG )
    Add_Dependencies(lib${PROJECT_NAME}_a Mmg)
  ENDIF ( )

  IF ( DOWNLOAD_METIS )
    Add_Dependencies(lib${PROJECT_NAME}_a Metis)
  ENDIF ( )

ENDIF()

# Compile shared library
OPTION ( LIBPARMMG_SHARED "Compile dynamic library" OFF)
IF ( LIBPARMMG_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}_so SHARED
    "${pmmg_library_files}" ${PROJECT_NAME} )

  IF ( DOWNLOAD_MMG )
    Add_Dependencies(lib${PROJECT_NAME}_so Mmg)
  ENDIF ( )

  IF ( DOWNLOAD_METIS )
    Add_Dependencies(lib${PROJECT_NAME}_so Metis)
  ENDIF ( )

ENDIF()

# parmmg header files needed for library
SET( pmmg_headers
  ${PMMG_SOURCE_DIR}/libparmmg.h
  ${PMMG_BINARY_DIR}/libparmmgf.h
  ${PMMG_SOURCE_DIR}/libparmmgtypes.h
  ${PMMG_BINARY_DIR}/libparmmgtypesf.h
  )

# Install header files in /usr/local or equivalent
INSTALL(FILES ${pmmg_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/parmmg COMPONENT headers )

COPY_FORTRAN_HEADER_AND_CREATE_TARGET ( ${PMMG_BINARY_DIR} ${PMMG_INCLUDE} )

# Copy header files in project directory at configuration step
# (generated file don't exists yet or are outdated)
FILE(INSTALL  ${pmmg_headers} DESTINATION ${PMMG_INCLUDE}
  PATTERN "libparmmg*f.h"  EXCLUDE)

install(EXPORT ParMmgTargets
  FILE ParMmgTargets.cmake
  NAMESPACE ParMmg::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/parmmg
  )

###############################################################################
#####
#####         Compile PMMG executable
#####
###############################################################################
ADD_AND_INSTALL_EXECUTABLE ( ${PROJECT_NAME}
  "${pmmg_library_files}" ${pmmg_main_file} )

IF ( DOWNLOAD_MMG )
  Add_Dependencies(parmmg Mmg)
ENDIF ( )

IF ( DOWNLOAD_METIS )
  Add_Dependencies(parmmg Metis)
ENDIF ( )

############################################################################
#####
#####         Print git hash to source file
#####
############################################################################

IF (NOT WIN32)

  ADD_CUSTOM_TARGET(GenerateGitHash
    COMMAND ./git_log_pmmg.sh ${PROJECT_SOURCE_DIR} ${PMMG_BINARY_DIR}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/
    COMMENT "Getting git commit hash"
  )
  ADD_DEPENDENCIES(${PROJECT_NAME} GenerateGitHash)
  IF( LIBPARMMG_STATIC )
    ADD_DEPENDENCIES(lib${PROJECT_NAME}_a GenerateGitHash)
  ENDIF ()
  IF( LIBPARMMG_SHARED )
    ADD_DEPENDENCIES(lib${PROJECT_NAME}_so GenerateGitHash)
  ENDIF ()

  INCLUDE_DIRECTORIES(${PMMG_BINARY_DIR})

ENDIF ()

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

  ADD_CUSTOM_TARGET(parmmg_doc
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


  if ( NOT TARGET doc )
    ADD_CUSTOM_TARGET(doc
      DEPENDS parmmg_doc
      COMMENT "Generating PARMMG API documentation with Doxygen.
     Open the doc/doxygen/html/index.html file to see it" VERBATIM
      )
  else()
    add_dependencies(doc parmmg_doc)
  endif()
ENDIF ( DOXYGEN_FOUND )

OPTION ( BUILD_TESTING "Enable / Disable tests" OFF )
CMAKE_DEPENDENT_OPTION (
  ONLY_LIBRARY_TESTS "Enable / Disable tests hosted outside the ParMmg repos
  (non library tests)" OFF  "BUILD_TESTING" OFF )


INCLUDE( ${PROJECT_SOURCE_DIR}/cmake/testing/pmmg_tests.cmake )
