## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
##
##  mmg is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  mmg is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with mmg (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the mmg distribution only if you accept them.
## =============================================================================

############################################################################
#####
#####         Scotch
#####
############################################################################
# Find SCOTCH library?
SET(SCOTCH_DIR "" CACHE PATH "Installation directory for scotch")

# add Scotch library?
SET ( USE_SCOTCH "" CACHE STRING "Use SCOTCH tool for renumbering (ON, OFF or <empty>)" )
SET_PROPERTY(CACHE USE_SCOTCH PROPERTY STRINGS "ON" "OFF" "")

IF ( NOT DEFINED USE_SCOTCH OR USE_SCOTCH STREQUAL "" OR USE_SCOTCH MATCHES " +" )
  # Variable is not provided by user
  FIND_PACKAGE(SCOTCH QUIET)

ELSE ()
  IF ( USE_SCOTCH )
    # User wants to use scotch
    FIND_PACKAGE(SCOTCH)
    IF ( NOT SCOTCH_FOUND )
      MESSAGE ( FATAL_ERROR "Scotch library not found:"
        "Using scotch reduce the execution time of mmg3d "
        "(see https://gforge.inria.fr/frs/?group_id=248 to download it)."
        "If you have already installed Scotch and want to use it, "
        "please set the CMake variable or environment variable SCOTCH_DIR "
        "to your scotch directory.")
    ENDIF ( )
  ENDIF ( )

ENDIF ( )

If ( SCOTCH_FOUND AND NOT USE_SCOTCH MATCHES OFF)
  add_definitions(-DUSE_SCOTCH)

  MESSAGE(STATUS
    "Compilation with scotch: ${SCOTCH_LIBRARIES}")
  SET( LIBRARIES ${SCOTCH_LIBRARIES} ${LIBRARIES})
ENDIF()


############################################################################
#####
#####         LinearElasticity
#####
############################################################################
# add LinearElasticity library?
SET(ELAS_DIR "" CACHE PATH "Installation directory for Elas")

SET ( USE_ELAS "" CACHE STRING "Use the Elas library for lagrangian motion option (ON, OFF or <empty>)" )
SET_PROPERTY(CACHE USE_ELAS PROPERTY STRINGS "ON" "OFF" "")

IF ( NOT DEFINED USE_ELAS OR USE_ELAS STREQUAL "" OR USE_ELAS MATCHES " +" )
  INCLUDE(FindElas)

ELSE()
  IF ( USE_ELAS )
    # User wants to use elas
    INCLUDE(FindElas)
    IF ( NOT ELAS_FOUND )
      MESSAGE ( FATAL_ERROR "Elas is a library to solve the linear elasticity "
        "problem (see https://github.com/ISCDtoolbox/LinearElasticity to"
        " download it). "
        "This library is needed to use the lagrangian motion option. "
        "If you have already installed Elas and want to use it, "
        "please set the CMake variable or environment variable ELAS_DIR "
        "to your Elas directory.")
    ENDIF ( )
  ENDIF ( )

ENDIF ( )

############################################################################
#####
#####         VTK (to parse (p)vtp/(p)vtu files )
#####
############################################################################
# add the VTK library ?
SET ( USE_VTK "" CACHE STRING "Use VTK I/O (ON, OFF or <empty>)" )
SET_PROPERTY(CACHE USE_VTK PROPERTY STRINGS "ON" "OFF" "")

IF ( NOT DEFINED USE_VTK OR USE_VTK STREQUAL "" OR USE_VTK MATCHES " +" OR USE_VTK )
  # Variable is not provided by the user or is setted to on

  # Handle vtk components name change between v8.2 and v9
  # Before v9
  FIND_PACKAGE(VTK QUIET)
  IF ( VTK_FOUND )
    message (STATUS "ParMmg - VTK_VERSION: ${VTK_VERSION}")
    IF (VTK_VERSION VERSION_LESS "9.0.0")
      find_package(VTK  COMPONENTS
        vtkCommonCore
        vtkCommonDataModel
        vtkIOLegacy
        vtkIOParallel
        vtkIOParallelXML
        vtkIOXML
        vtkParallelCore
        vtkParallelMPI
        QUIET)
    ELSE()
      # After v9
      FIND_PACKAGE(VTK COMPONENTS
        CommonCore
        CommonDataModel
        IOLegacy
        IOParallel
        IOParallelXML
        IOXML
        ParallelCore
        ParallelMPI
        QUIET)
    ENDIF()

  ELSEIF ( USE_VTK )
    # USE_VTK is not empty so user explicitely ask for VTK...
    # but it is not found: raise an error
    MESSAGE(FATAL_ERROR "ParMmg - VTK library not found.")
  ENDIF()
ENDIF()

IF ( VTK_FOUND AND NOT USE_VTK MATCHES OFF)

  add_definitions(-DUSE_VTK)

  MESSAGE ( STATUS "ParMmg - Compilation with VTK: add vtk, vtp and vtu I/O." )

  IF( "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}" LESS 8.90 )
    INCLUDE ( ${VTK_USE_FILE} )
  ENDIF()

  SET( LIBRARIES ${VTK_LIBRARIES} ${LIBRARIES} )
ENDIF ( )

############################################################################
#####
#####         HDF5
#####
############################################################################
# Find HDF5 library?
SET ( USE_HDF5 "" CACHE STRING "Use HDF5 I/O (ON, OFF or <empty>)" )
SET_PROPERTY(CACHE USE_HDF5 PROPERTY STRINGS "ON" "OFF" "")

IF ( NOT DEFINED USE_HDF5 OR USE_HDF5 STREQUAL "" OR USE_HDF5 MATCHES " +" )
  # Variable is not provided by user
  SET(HDF5_PREFER_PARALLEL TRUE)

  FIND_PACKAGE(HDF5 QUIET)
  IF ( HDF5_FOUND AND NOT HDF5_IS_PARALLEL )
    UNSET ( HDF5_FOUND )
  ENDIF( )

ELSE ( )
  IF ( USE_HDF5 )
    # User wants to use HDF5 I/O
    SET(HDF5_PREFER_PARALLEL TRUE)

    FIND_PACKAGE(HDF5)

    IF ( NOT HDF5_FOUND )
      MESSAGE ( FATAL_ERROR "HDF5 library not found."
        " Please set the CMake variable USE_HDF5 to OFF to disable HDF5 I/Os.")
    ENDIF ( )

    IF ( NOT HDF5_IS_PARALLEL )
      MESSAGE ( FATAL_ERROR "HDF5 parallel library not found: "
        "Please enable parallel support for HDF5 or set the CMake variable "
        "USE_HDF5 to OFF to disable HDF5 I/Os.")
    ENDIF( )

  ENDIF ( )
ENDIF ( )

IF ( HDF5_FOUND AND NOT USE_HDF5 MATCHES OFF )

  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_HDF5")
  MESSAGE ( STATUS "Compilation with HDF5: add HDF5 I/O." )
  SET( LIBRARIES  ${LIBRARIES} ${HDF5_LIBRARIES} )

ENDIF ( )
