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
 * \file linkedlist_pmmg.c
 * \brief functions to manage a linked list
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "linkedlist_pmmg.h"
#include "coorcell_pmmg.h"

/**
 * \param parmesh   pointer toward the parmesh structure.
 * \param list      pointer toward a PMMG_lnkdList structure.
 * \param nitem_max maximal number of cells of the linked list.
 *
 * \return 1 if success, 0 if fail.
 *
 * Initialisation of a linked list of cells.
 *
 */
int PMMG_cellLnkdListNew( PMMG_pParMesh parmesh,PMMG_cellLnkdList *list,int id,int nitem_max )
{

  /* adjust hash table params */
  list->nitem     = 0;
  list->nitem_max = nitem_max;
  list->id        = id;

  PMMG_MALLOC(parmesh,list->item,nitem_max,PMMG_lnkdCell,"linked list array",
              return 0);

  list->frst = PMMG_UNSET;

  return 1;
}

/**
 * \param parmesh   pointer toward the parmesh structure.
 * \param list      pointer toward a PMMG_lnkdList structure.
 * \param nitem_max maximal number of cells of the linked list.
 *
 * \return 1 if success, 0 if fail.
 *
 * Initialisation of a linked list of cells.
 *
 */
int PMMG_valLnkdListNew( PMMG_pParMesh parmesh,PMMG_valLnkdList *list,int id,int nitem_max )
{

  /* adjust hash table params */
  list->nitem     = 0;
  list->nitem_max = nitem_max;
  list->id        = id;

  PMMG_MALLOC(parmesh,list->item,nitem_max,PMMG_lnkdVal,"linked list array",
              return 0);

  list->frst = PMMG_UNSET;

  return 1;
}

/**
 * \param parmesh   pointer toward the parmesh structure.
 * \param list      pointer toward a PMMG_lnkdList structure.
 *
 * Reset a linked list without deallocating it.
 *
 */
void PMMG_reset_cellLnkdList( PMMG_pParMesh parmesh,PMMG_cellLnkdList *list ) {

  /* adjust hash table params */
  list->nitem  = 0;
  list->id     = PMMG_NUL;
  list->frst   = PMMG_UNSET;

  return;
}

/**
 * \param parmesh   pointer toward the parmesh structure.
 * \param list      pointer toward a PMMG_lnkdList structure.
 * \param val1      value to add to the \a val1 field of the new cell.
 * \param val2      value to add to the \a val2 field of the new cell.
 *
 * \return 1 if we insert a new cell, 2 if the cell was already present and 0
 * if fail.
 *
 * Sorted insertion (with unicity) of a new cell with values \a val1 and \a val2
 * in the linked list \a list.
 *
 */
int PMMG_add_cell2lnkdList( PMMG_pParMesh parmesh,PMMG_cellLnkdList *list,
                            int val1,int val2 ) {
  PMMG_lnkdCell *cell;
  int           k,id,prevId,newsize;

  /** Start from the first cell of the list */

  if ( list->frst<0 ) {
    /* Add the first element */
    list->frst    = 0;
    cell          = &list->item[list->nitem++];
    cell->val1    = val1;
    cell->val2    = val2;
    cell->id      = 0;
    cell->nxt     = -1;

    return 1;
  }

  /** Research of the cell position */
  k = 0;
  cell = &list->item[list->frst];
  while ( k<list->nitem &&
          ( (cell->val1<val1)||((cell->val1==val1)&&(cell->val2<val2)) ) ) {
    prevId   = cell->id;
    cell     = &list->item[cell->nxt];
    ++k;
  }

  if ( (k==list->nitem) ||
       ( (k<list->nitem) && ( (cell->val1 != val1) || (cell->val2!=val2) ) ) ) {
    /** The cell has not be founded, add it */
    if ( list->nitem >= list->nitem_max ) {
      newsize = MG_MAX((int)(1.2*list->nitem_max),list->nitem_max+1);
      PMMG_REALLOC(parmesh,list->item,newsize,
                   list->nitem_max,PMMG_lnkdCell,"linked list",return 0);
      list->nitem_max = newsize;
    }
    cell          = &list->item[list->nitem];
    cell->id      = list->nitem++;
    cell->val1    = val1;
    cell->val2    = val2;

    if ( k==0 ) {
      /* Modif of the head of the linked list */
      id         = list->frst;
      list->frst = cell->id;
    }
    else {
      id                      = list->item[prevId].nxt;
      list->item[prevId].nxt  = cell->id;
    }
    cell->nxt  = id;

    return 1;
  }

  /** This cell is alredy listed */
  assert ( cell->val1 == val1 && cell->val2 == val2 );

  return 2;
}

/**
 * \param parmesh   pointer toward the parmesh structure.
 * \param list      pointer toward a PMMG_cellLnkdList structure.
 * \param val1      first Cell value (removed)
 * \param val2      first Cell value (removed)
 *
 * \return 1
 *
 * Remove the first item of a list and fill \a val1 and \a val2 with the values
 * of the removed cell.
 *
 */
int PMMG_pop_cell_lnkdList( PMMG_pParMesh parmesh,PMMG_cellLnkdList *list,
                            int *val1,int *val2 ) {
  PMMG_lnkdCell *cell;

  /** Get first cell */
  cell = &list->item[list->frst];
  *val1 = cell->val1;
  *val2 = cell->val2;

  /** Pop cell from the head of the list */
  list->frst = cell->nxt;
  list->nitem--;

  return 1;
}

/**
 * \param parmesh   pointer toward the parmesh structure.
 * \param list      pointer toward a PMMG_lnkdList structure.
 * \param val1      value to add to the \a val1 field of the new cell.
 * \param val2      value to add to the \a val2 field of the new cell.
 *
 * \return 1 if we insert a new cell and 0 if fail.
 *
 * Non sorted insertion of a new cell with values \a val1 and \a val2
 * in the linked list \a list.
 *
 */
int PMMG_add_val2lnkdList( PMMG_pParMesh parmesh,PMMG_valLnkdList *list,
                           int val ) {
  PMMG_lnkdVal *cell;
  int           k,newsize;

  /** Start from the first cell of the list */

  if ( list->frst<0 ) {
    /* Add the first element */
    list->frst    = 0;
    cell          = &list->item[list->nitem++];
    cell->val     = val;
    cell->nxt     = -1;

    return 1;
  }

  /** Cell insertion at the end of the list */
  if ( list->nitem >= list->nitem_max ) {
    newsize = MG_MAX((int)(1.2*list->nitem_max),list->nitem_max+1);
    PMMG_REALLOC(parmesh,list->item,newsize,
                 list->nitem_max,PMMG_lnkdVal,"linked list",return 0);
    list->nitem_max = newsize;
  }
  k             = list->nitem++;
  cell          = &list->item[k];
  cell->val     = val;
  cell->nxt     = -1;

  list->item[k-1].nxt = k;

  return 1;
}

/**
 * \param parmesh   pointer toward the parmesh structure.
 * \param list      pointer toward a PMMG_valLnkdList structure.
 * \param val       first list value (removed)
 *
 * \return 1
 *
 * Remove the first item of a list and fill \a val with this removed value.
 *
 */
int PMMG_pop_val_lnkdList( PMMG_pParMesh parmesh,PMMG_valLnkdList *list,
                           int *val ) {
  PMMG_lnkdVal *cell;

  /** Get first cell */
  cell = &list->item[list->frst];
  *val = cell->val;

  /** Pop cell from the head of the list */
  list->frst = cell->nxt;
  list->nitem--;

  return 1;
}


/**
 * \param parmesh   pointer toward the parmesh structure.
 * \param list1     pointer toward a PMMG_lnkdList structure.
 * \param list2     pointer toward a PMMG_lnkdList structure.
 *
 * \return 1 if we success to merge the list \a list2 into \a list1, 2 if there
 * is nothing to do (no items in the list \a list2 that are not in \a list1) and
 * 0 if fail.
 *
 * Sorted insertion (with unicity) of a new cell with values \a val1 and \a val2
 * in the linked list \a list.
 *
 */
int PMMG_merge_cellLnkdList( PMMG_pParMesh parmesh,PMMG_cellLnkdList *list1,
                             PMMG_cellLnkdList *list2) {
  PMMG_lnkdCell *cell;
  int           k,ier,minier;

  if ( !list2->nitem ) return 2;

  /** Start from the first cell of the list */
  assert ( list2->frst >=0 );

  /** Insertion of the next cell of list2 into list1 */
  cell = &list2->item[list2->frst];

  minier = 2;
  for ( k=0; k<list2->nitem; ++k ) {
    ier  = PMMG_add_cell2lnkdList( parmesh,list1,cell->val1,cell->val2);
    if ( !ier ) return ier;

    cell = &list2->item[cell->nxt];
    /** if there is at least 1 successfull insertion, we want to return that
     * list 1 has changed */
    minier = MG_MIN(ier,minier);
  }

  return minier;
}

/**
 * \param list      pointer toward a PMMG_lnkdList structure.
 *
 * Print the linked list \a list.
 *
 * \remark for debug purpose
 */
void PMMG_print_cellLnkdList( PMMG_cellLnkdList *list ) {
  PMMG_lnkdCell *cell;
  int           k;

  /** Process the cell list */
  k = 0;
  cell = &list->item[list->frst];
  while ( k<list->nitem ) {
    printf(" %s: value number %d: %d %d\n",__func__,k,cell->val1,cell->val2);
    ++k;
    cell = &list->item[cell->nxt];
  }
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param list    pointer toward a PMMG_lnkdList structure.
 * \param array   pointer toward the array in which to pack the list.
 *
 * \return the number of added items.
 *
 * Copy in an array the size of a linked list followed by the \a val1 and \a
 * val2 fields of each cell.
 *
 */
int PMMG_packInArray_cellLnkdList( PMMG_cellLnkdList *list, int *array ) {
  PMMG_lnkdCell *cell;
  int           k,idx;

  /** Start from the first cell of the list */
  array[0] = list->nitem;

  /** Store each cell into the array */
  idx  = 1;
  k    = 0;
  cell = &list->item[list->frst];

  while ( k<list->nitem ) {
    array[idx++] = cell->val1;
    array[idx++] = cell->val2;
    ++k;
    cell = &list->item[cell->nxt];
  }

  return idx;
}

/**
 * \param list  pointer toward a PMMG_lnkdList structure.
 * \param array pointer toward the array in which to pack the list.
 *
 * \return the number of unpacked items if success, -1 if fail.
 *
 * Copy the values packed by the \a PMMG_unpackInArray_linkedList function in \a
 * array (starting from \a start) into the linked list \a list.
 *
 */
int PMMG_unpackArray_inCellLnkdList(  PMMG_pParMesh parmesh, PMMG_cellLnkdList *list,
                                      int *array ) {
  PMMG_lnkdCell *cell;
  int           k,idx;

  /** Start from the first cell of the list */
  list->nitem = array[0];

  if ( list->nitem >= list->nitem_max ) {
    PMMG_REALLOC(parmesh,list->item,list->nitem,list->nitem_max,
                 PMMG_lnkdCell,"linked list",return -1);
    list->nitem_max = list->nitem;
  }


  /** Research of the cell position */
  list->frst    = 0;
  idx           = 1;
  for ( k=0; k<list->nitem; ++k ) {
    cell          = &list->item[k];
    cell->val1    = array[idx++];
    cell->val2    = array[idx++];
    cell->id      = k;
    cell->nxt     = k+1;
  }

  cell->nxt = -1;

  return idx;
}

/**
 * \param a  pointer toward a PMMG_lnkdList* structure.
 * \param b  pointer toward a PMMG_lnkdList* structure.
 *
 * \return 1 if a is greater than b, -1 if b is greater than 1, 0 if they are
 * equals.
 *
 * Compare 2 linked lists (can be used inside the qsort C fnuction), first on
 * their lengths, second ond their \a val1 field then on their \a val2 one.
 *
 */
int PMMG_compare_cellLnkdList (const void * a, const void * b) {
  PMMG_cellLnkdList *list1,*list2;
  PMMG_lnkdCell     *cell1,*cell2;
  int              k;

  list1 = *(PMMG_cellLnkdList**)a;
  list2 = *(PMMG_cellLnkdList**)b;

  if ( list1->nitem > list2->nitem ) return 1;

  if ( list1->nitem < list2->nitem ) return -1;

  cell1 = &list1->item[list1->frst];
  cell2 = &list2->item[list2->frst];

  for ( k=0; k<list1->nitem; ++k ) {
    if ( cell1->val1 > cell2->val1 ) return 1;
    if ( cell1->val1 < cell2->val1 ) return -1;

    if ( cell1->val2 > cell2->val2 ) return 1;
    if ( cell1->val2 < cell2->val2 ) return -1;

    cell1 = &list1->item[cell1->nxt];
    cell2 = &list2->item[cell2->nxt];
  }

  return 0;
}

/**
 * \param a  pointer toward a PMMG_lnkdList* structure.
 * \param b  pointer toward a PMMG_lnkdList* structure.
 *
 * \return 1 if a is greater than b, -1 if b is greater than 1, 0 if they are
 * equals.
 *
 * Compare lengths of 2 linked lists (can be used inside the qsort C fnuction)
 *
 */
int PMMG_compare_valLnkdListLen (const void * a, const void * b) {
  PMMG_valLnkdList *list1,*list2;

  list1 = *(PMMG_valLnkdList**)a;
  list2 = *(PMMG_valLnkdList**)b;

  if ( list1->nitem > list2->nitem ) return 1;

  if ( list1->nitem < list2->nitem ) return -1;

  return 0;
}

/**
 * \param a  pointer toward a PMMG_lnkdCell structure.
 * \param b  pointer toward a PMMG_lnkdCell structure.
 *
 * \return 1 if val2 of cell \a a is greater than val2 of cell \a b, -1 if
 * smaller 0 if equals.
 *
 * Compare two cells from their second field (\a val2).
 *
 */
int PMMG_compare_cell2 (const void * a, const void * b) {
  PMMG_lnkdCell cell1,cell2;

  cell1 = *(PMMG_lnkdCell*)a;
  cell2 = *(PMMG_lnkdCell*)b;

  if ( cell1.val2 > cell2.val2 ) return 1;
  if ( cell1.val2 < cell2.val2 ) return -1;

  return 0;
}

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param array1  pointer toward an integer array.
 * \param array2  pointer toward another integer array.
 * \param oldIdx  array of the initial item position.
 *
 * \return 1.
 *
 * Sort 2 arrays in the ascending order of the second one.
 * The first array is not mandatory.
 *
 */
int PMMG_sort_iarray( PMMG_pParMesh parmesh,int *array1,int *array2,int *oldIdx,
                      int nitem ) {
  PMMG_lnkdCell *cell;
  int           i;

  /* Initialize and fill lists of int pairs */
  PMMG_CALLOC(parmesh,cell,nitem,PMMG_lnkdCell,"array of cells",return 0);

  if ( array1 ) {
    for( i = 0; i < nitem; i++ ) {
      cell[i].val1 = array1[i];
      cell[i].val2 = array2[i];
      cell[i].id   = i;
    }
  }
  else {
    for( i = 0; i < nitem; i++ ) {
      cell[i].val2 = array2[i];
      cell[i].id   = i;
    }
  }

  /* Sort lists based on values in array2, in ascending order */
  qsort(cell,nitem,sizeof(PMMG_lnkdCell),PMMG_compare_cell2);

  /* Permute arrays and deallocate lists */
  if ( array1 ) {
    for( i = 0; i < nitem; i++ ) {
      array1[i] = cell[i].val1;
      array2[i] = cell[i].val2;
      oldIdx[i] = cell[i].id;
    }
  }
  else {
    for( i = 0; i < nitem; i++ ) {
      array2[i] = cell[i].val2;
      oldIdx[i] = cell[i].id;
    }
  }

  PMMG_DEL_MEM(parmesh,cell,PMMG_lnkdCell,"array of cells");
  return 1;
}
