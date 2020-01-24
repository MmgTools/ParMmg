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
 * \file linkedlist_pmmg.h
 * \brief linkedlist_pmmg.c header file
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#ifndef LINKEDLIST_PMMG_H

#define LINKEDLIST_PMMG_H

#include "parmmg.h"

#define PMMG_LISTSIZE 30

/**
 * \struct PMMG_lnkdCell
 *
 * \brief Cell of a linked list. This cell allow to store 2 values and an index.
 *
 */
typedef struct {
  int     val1; /*!< first value of the cell */
  int     val2; /*!< second value of the cell */
  int     id;   /*!< position of the cell in the array */
  int     nxt;  /*!< position of the next cell in the linked list */
} PMMG_lnkdCell;

/**
 * \struct PMMG_lnkdVal
 *
 * \brief Cell of a linked list. This cell allow to store 1 value.
 *
 */
typedef struct {
  int     val; /*!< value of the cell */
  int     nxt;  /*!< position of the next cell in the linked list */
} PMMG_lnkdVal;

/**
 * \struct PMMG_lnkdList
 *
 * \brief Linked list of Cells of. Each cell can store 2 values.
 *
 */
typedef struct {
  int nitem; /*!< number of used cells in the list (=position to insert the next cell) */
  int nitem_max; /*!< maximal number of item in the list */
  int frst; /*!< position of the first cell of the linked list */
  int id;  /*!< ID of the linked list */
  PMMG_lnkdCell *item; /*!< array of cells */
} PMMG_cellLnkdList;

/**
 * \struct PMMG_lnkdList
 *
 * \brief Linked list of cells. Each cell can store 2 values.
 *
 */
typedef struct {
  int nitem; /*!< number of used cells in the list (=position to insert the next cell) */
  int nitem_max; /*!< maximal number of item in the list */
  int frst; /*!< position of the first cell of the linked list */
  int id;  /*!< ID of the linked list */
  PMMG_lnkdVal *item; /*!< array of Values */
} PMMG_valLnkdList;

int  PMMG_cellLnkdListNew( PMMG_pParMesh parmesh,PMMG_cellLnkdList *list,int,int );
int  PMMG_valLnkdListNew( PMMG_pParMesh parmesh,PMMG_valLnkdList *list,int,int );
void PMMG_reset_cellLnkdList( PMMG_pParMesh parmesh,PMMG_cellLnkdList *list );
int  PMMG_merge_cellLnkdList( PMMG_pParMesh,PMMG_cellLnkdList*,PMMG_cellLnkdList *);
int  PMMG_add_cell2lnkdList( PMMG_pParMesh,PMMG_cellLnkdList*,int,int );
int  PMMG_pop_cell_lnkdList( PMMG_pParMesh,PMMG_cellLnkdList*,int*,int* );
int  PMMG_add_val2lnkdList( PMMG_pParMesh,PMMG_valLnkdList*,int);
int  PMMG_pop_val_lnkdList( PMMG_pParMesh,PMMG_valLnkdList*,int*);
void PMMG_print_cellLnkdList( PMMG_cellLnkdList *list );
int  PMMG_packInArray_cellLnkdList( PMMG_cellLnkdList *list,int *array );
int  PMMG_unpackArray_inCellLnkdList(  PMMG_pParMesh,PMMG_cellLnkdList*,int* );
int  PMMG_compare_cellLnkdList (const void * a, const void * b);
int  PMMG_compare_valLnkdListLen (const void * a, const void * b);
int  PMMG_sort_iarray( PMMG_pParMesh parmesh,int *array1,int *array2,int *oldIdx,
                       int nitem );

#endif
