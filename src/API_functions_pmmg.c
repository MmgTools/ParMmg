/**
 * \file API_functions_pmmg.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * C API for PARMMG library.
 *
 */

#include "libparmmg.h"

int PMMG_Init_parMesh( const int starter,... ) {
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = _PMMG_Init_parMesh_var(argptr);

  va_end(argptr);

  return ier;
}
