/*===========================================================================*/
/* Original source from caTools library:                                     */
/* runfunc - running window functions                                        */
/* Copyright (C) 2005 Jarek Tuszynski                                        */
/* Distributed under GNU General Public License version 3                    */
/*===========================================================================*/

/*==================================================*/
/* Index:                                           */
/*  |------------------+------+------+----------|   */
/*  | function         | NaN  | Edge | Underflow|   */
/*  |------------------+------+------+----------|   */
/*  | movmin           | yes  | yes  |   NA     |   */
/*  | movmax           | yes  | yes  |   NA     |   */
/*  |------------------+------+------+----------|   */
/*  NaN - means support for NaN and possibly Inf    */
/*  edge - means calculations are done all the way  */
/*         to the edges                             */
/*  underflow - means at maximum how many numbers   */
/*    are used to store results of addition in case */
/*    of underflow                                  */
/*==================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <float.h>
#include <R.h>
#include <Rinternals.h>

#define notNaN(x)   ((x)==(x))
#define isNaN(x)  (!((x)==(x)))
#define MIN(y,x) ((x)<(y) && (x)==(x) ? (x) : (y))
#define MAX(y,x) ((x)>(y) && (x)==(x) ? (x) : (y))
#define SQR(x) ((x)*(x))


/*==================================================================*/
/* minimum function applied to moving (running) window              */
/* Input :                                                          */
/*   In   - array to run moving window over will remain umchanged   */
/*   Out  - empty space for array to store the results. Out is      */
/*          assumed to have reserved memory for nIn*nProbs elements */
/*   nIn  - size of arrays In and Out                               */
/*   nWin - size of the moving window (odd)                         */
/* Output :                                                         */
/*   Out  - results of runing moving window over array In and       */
/*          colecting window mean                                   */
/*==================================================================*/
void movmin(double *In, double *Out, const int *nIn, const int *nWin)
{ /* full-blown version with NaN's and edge calculation */
  int i, j, k2, n=*nIn, m=*nWin;
  double ptOut, Min, *in, *out, CST = DBL_MAX;
  double NaN = (0.0/0.0);

  k2  = m>>1;               /* right half of window size */
  in  = In;
  out = Out;
  /* --- step 1 - find min of elements 0:(k2-1) */
  Min=CST;                  /* we need to calculate  initial 'Min' */
  for(i=0; i<k2; i++) Min = MIN(Min,in[i]);  /* find minimum over a window of length k2 */
  /* --- step 2 - left edge - start expanding the moving window to the right */
  for(i=k2; i<m-1; i++) {
    Min=MIN(Min,in[i]);     /* cumulative min */
    *(out++) = (Min==CST ? NaN : Min); /* save 'Min' and move window */
  }
  /* --- step 3 - the inner section - window of constant size is moving  */
  ptOut=CST;
  for(i=m-1; i<n; i++) {
    if(ptOut==Min) {        /* if point comining out of the window was window's min than ... */
      Min=CST;              /* we need to recalculate 'Min' */
      for(j=0; j<m; j++)
        Min=MIN(Min,in[j]); /* find minimum over a window of length m */
    } else                  /* if point comining out of the window was NOT window min than min of ... */
      Min=MIN(Min,in[m-1]); /* ... window's first m-1 points is still 'Min', so we have to add a single point */
    ptOut = *(in++);        /* store point comming out of the window for future use and move window */
    *(out++) = (Min==CST ? NaN : Min); /* save 'Min' and move window */
  }
  /* --- step 4 - right edge - right side reached the end and left is shrinking  */
  for(i=0; i<k2; i++) {
    if(ptOut==Min) {        /* if point comining out of the window was window's extreme than ... */
      Min=CST;              /* we need to recalculate 'Min' */
      for(j=0; j<m-i-1; j++)
        Min=MIN(Min,in[j]); /* find minimum over a window of length m */
    }
    ptOut = *(in++);        /* store point comming out of the window for future use and move window */
    *(out++) = (Min==CST ? NaN : Min);  /* and fill the space with window extreme and move window */
  }
}

/*==================================================================*/
/* Maximum function applied to moving (running) window              */
/* Input :                                                          */
/*   In   - array to run moving window over will remain umchanged   */
/*   Out  - empty space for array to store the results. Out is      */
/*          assumed to have reserved memory for nIn*nProbs elements */
/*   nIn  - size of arrays In and Out                               */
/*   nWin - size of the moving window (odd)                         */
/* Output :                                                         */
/*   Out  - results of runing moving window over array In and       */
/*          colecting window mean                                   */
/*==================================================================*/
void movmax(double *In, double *Out, const int *nIn, const int *nWin)
{ /* full-blown version with NaN's and edge calculation */
  int i, j, k2, n=*nIn, m=*nWin;
  double ptOut, Max, *in, *out, CST = -DBL_MAX;
  double NaN = (0.0/0.0);

  k2  = m>>1;               /* right half of window size */
  in  = In;
  out = Out;
  /* step 1 - find max of elements 0:(k2-1) */
  Max= CST;                /* we need to calculate  initial 'Max' */
  for(i=0; i<k2; i++) Max = MAX(Max,in[i]);  /* find maximum over a window of length k2 */
  /* step 2 - left edge - start expanding the moving window to the right */
  for(i=k2; i<m-1; i++) {
    Max=MAX(Max,in[i]);     /* cumulative max */
    *(out++) = (Max==CST ? NaN : Max); /* save 'Max' and move window */
  }
  /* step 3 - the inner section - window of constant size is moving  */
  ptOut=CST;
  for(i=m-1; i<n; i++) {
    if(ptOut==Max) {        /* if point comaxing out of the window was window's max than ... */
      Max=CST;              /* we need to recalculate 'Max' */
      for(j=0; j<m; j++)
        Max=MAX(Max,in[j]); /* find maximum over a window of length m */
    } else                  /* if point comining out of the window was NOT window max than max of ... */
      Max=MAX(Max,in[m-1]); /* ... window's first m-1 points is still 'Max', so we have to add a single point */
    ptOut = *(in++);        /* store point comming out of the window for future use and move window */
    *(out++) = (Max==CST ? NaN : Max); /* save 'Max' and move window */
  }
  /* step 4 - right edge - right side reached the end and left is shrinking  */
  for(i=0; i<k2; i++) {
    if(ptOut==Max) {        /* if point comining out of the window was window's extreme than ... */
      Max=CST;              /* we need to recalculate 'Max' */
      for(j=0; j<m-i-1; j++)
        Max=MAX(Max,in[j]); /* find maximum over a window of length m */
    }
    ptOut = *(in++);        /* store point comming out of the window for future use and move window */
    *(out++) = (Max==CST ? NaN : Max); /* and fill the space with window extreme and move window */
  }
}

