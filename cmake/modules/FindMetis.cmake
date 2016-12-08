## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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

IF ((NOT WIN32) AND (NOT WIN64))
  SET ( METIS_INCLUDE_DIR METIS_INCLUDE_DIR-NOTFOUND )
  SET ( METIS_LIBRARY METIS_LIBRARY-NOTFOUND )
  SET ( PARMETIS_LIBRARY PARMETIS_LIBRARY-NOTFOUND )
ENDIF()

FIND_PATH(METIS_INCLUDE_DIR
  NAMES metis.h
  HINTS ${METIS_INCLUDE_DIR}
  $ENV{METIS_INCLUDE_DIR}
  ${METIS_DIR}/include
  $ENV{METIS_DIR}/include
  PATH_SUFFIXES metis
  DOC "Directory of METIS Header")

# Check for metis
FIND_LIBRARY(METIS_LIBRARY
  NAMES metis metis${METIS_LIB_SUFFIX}
  HINTS ${METIS_LIBRARY}
  $ENV{METIS_LIBRARY}
  ${METIS_DIR}/lib
  $ENV{METIS_DIR}/lib
  ${METIS_DIR}/build/.*/libmetis
  $ENV{METIS_DIR}/build/.*/libmetis
  PATH_SUFFIXES metis
  DOC "The METIS library"
  )

#FIND_LIBRARY(PARMETIS_LIBRARY
#  NAMES parmetis parmetis${METIS_LIB_SUFFIX}
#  HINTS ${PARMETIS_LIBRARY}
#  $ENV{PARMETIS_LIBRARY}
#  ${PARMETIS_DIR}/lib
#  $ENV{PARMETIS_DIR}/lib
#  DOC "The PARMETIS library"
#  )

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(METIS DEFAULT_MSG
  METIS_INCLUDE_DIR METIS_LIBRARY
  #PARMETIS_LIBRARY
  )

IF ((NOT WIN32) AND (NOT WIN64))
  MARK_AS_ADVANCED(METIS_INCLUDE_DIR METIS_LIBRARY PARMETIS_LIBRARY)
ENDIF()
