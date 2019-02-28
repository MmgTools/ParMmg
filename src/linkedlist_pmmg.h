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

#endif
