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
 * \brief Cell of a linked list. This cell allow to store 2 values.
 *
 */
typedef struct {
  int     val1; /*!< first value of the cell */
  int     val2; /*!< second value of the cell */
  int     id;   /*!< position of the cell in the array */
  int     nxt;  /*!< position of the next cell in the linked list */
} PMMG_lnkdCell;

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
  PMMG_lnkdCell *item; /*!< array of cells */
} PMMG_lnkdList;

int  PMMG_lnkdListNew( PMMG_pParMesh parmesh,PMMG_lnkdList *list,int,int );
void PMMG_reset_lnkdList( PMMG_pParMesh parmesh,PMMG_lnkdList *list );
int  PMMG_merge_lnkdList( PMMG_pParMesh,PMMG_lnkdList*,PMMG_lnkdList *);
int  PMMG_add_cell2lnkdList( PMMG_pParMesh,PMMG_lnkdList*,int,int );
int  PMMG_pop_cell_lnkdList( PMMG_pParMesh,PMMG_lnkdList*,int*,int* );
void PMMG_print_lnkdList( PMMG_lnkdList *list );
int  PMMG_packInArray_lnkdList( PMMG_lnkdList *list,int *array );
int  PMMG_unpackArray_inLnkdList(  PMMG_pParMesh,PMMG_lnkdList*,int* );
int  PMMG_compare_lnkdList (const void * a, const void * b);
int  PMMG_sort_iarray( PMMG_pParMesh parmesh,int *array1,int *array2,int *oldIdx,
                       int nitem );

#endif
