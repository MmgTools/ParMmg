/* =============================================================================
**  This file is part of the parmmg software package for parallel tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017-
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with parmmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the parmmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file libparmmg_toolsf.c
 * \brief Fortran API functions for PARMMG library.
 * \author Algiane Froehly (Inria)
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \date 09 2020
 * \copyright GNU Lesser General Public License.
 * \note Please, refer to the \ref libparmmg.h file for functions
 * documentation.
 *
 * Define the Fortran API functions for PARMMG library: adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */

#include "parmmg.h"

/**
 * See \ref PMMG_printCommunicator function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_PRINTCOMMUNICATOR,pmmg_printcommunicator,
             (PMMG_pParMesh *parmesh,char* filename, int *strlen,int* retval),
             (parmesh,filename,strlen, retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen+1,char,return);
  strncpy(tmp,filename,*strlen);
  tmp[*strlen] = '\0';

  *retval = PMMG_printCommunicator(*parmesh,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}
