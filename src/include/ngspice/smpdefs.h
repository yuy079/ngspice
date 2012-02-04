#ifndef SMP
#define SMP

typedef  struct MatrixFrame     SMPmatrix;
typedef  struct MatrixElement  *SMPelement;

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000  AlansFixes
**********/

#include <stdio.h>
#include <math.h>
#include "ngspice/complex.h"

#ifdef KLU
#include "ngspice/klu.h"
#endif

#ifdef KLU
void SMPmatrix_CSC (SMPmatrix *, int **, int **, double **, int, double **, double **, double **) ;
void SMPnnz (SMPmatrix *, int *, int *) ;
#endif

int SMPaddElt( SMPmatrix *, int , int , double );
double * SMPmakeElt( SMPmatrix * , int , int );

#ifdef KLU
void SMPcClear( SMPmatrix *, double *, int );
#else
void SMPcClear( SMPmatrix *);
#endif

#ifdef KLU
void SMPclear( SMPmatrix *, double *, int);
#else
void SMPclear( SMPmatrix *);
#endif

#ifdef KLU
int SMPcLUfac( SMPmatrix *, int *, int *, double *, klu_symbolic *, klu_numeric *, klu_common *, double, int );
#else
int SMPcLUfac( SMPmatrix *, double );
#endif

#ifdef KLU
int SMPluFac( SMPmatrix *, int *, int *, double *, klu_symbolic *, klu_numeric *, klu_common *, double **, double, double, int );
#else
int SMPluFac( SMPmatrix *, double , double );
#endif
#ifdef KLU
int SMPcReorder( SMPmatrix *, int *, int *, double *, klu_symbolic **, klu_numeric **, klu_common *, double, double, int *, int );
#else
int SMPcReorder( SMPmatrix * , double , double , int *);
#endif
#ifdef KLU
int SMPreorder( SMPmatrix * , int *, int *, double *, klu_symbolic *, klu_numeric **, klu_common *, double **, double, double, double, int );
#else
int SMPreorder( SMPmatrix * , double , double , double );
#endif
void SMPcaSolve(SMPmatrix *Matrix, double RHS[], double iRHS[],
		double Spare[], double iSpare[]);
#ifdef KLU
void SMPcSolve( SMPmatrix *, klu_symbolic *, klu_numeric *, klu_common *, double [], double [], double [], double [], double [], int );
#else
void SMPcSolve( SMPmatrix *, double [], double [], double [], double []);
#endif
#ifdef KLU
void SMPsolve( SMPmatrix *, klu_symbolic *, klu_numeric *, klu_common *, double [], double [], double [], int);
#else
void SMPsolve( SMPmatrix *, double [], double []);
#endif

int SMPmatSize( SMPmatrix *);
int SMPnewMatrix( SMPmatrix ** );

#ifdef KLU
void SMPdestroy( SMPmatrix *, int **, int **, double **, klu_symbolic **, klu_numeric **, klu_common *, double ***, double ***, double ***, double ***, double **, double **, int );
#else
void SMPdestroy( SMPmatrix *);
#endif

#ifdef KLU
int SMPpreOrder( SMPmatrix *, int *, int *, klu_symbolic **, klu_common *, int );
#else
int SMPpreOrder( SMPmatrix *);
#endif

void SMPprint( SMPmatrix * , FILE *);
void SMPgetError( SMPmatrix *, int *, int *);
int SMPcProdDiag( SMPmatrix *, SPcomplex *, int *);
int SMPcDProd(SMPmatrix *Matrix, SPcomplex *pMantissa, int *pExponent);
SMPelement * SMPfindElt( SMPmatrix *, int , int , int );
int SMPcZeroCol(SMPmatrix *eMatrix, int Col);
int SMPcAddCol(SMPmatrix *eMatrix, int Accum_Col, int Addend_Col);
int SMPzeroRow(SMPmatrix *eMatrix, int Row);

/* Correction for the Spertica's hack */
extern void SMPgmo ( SMPmatrix *, int, double * ) ;
/* End of Correction for the Spertica's hack */

#ifdef PARALLEL_ARCH
void SMPcombine(SMPmatrix *Matrix, double RHS[], double Spare[]);
void SMPcCombine(SMPmatrix *Matrix, double RHS[], double Spare[],
		 double iRHS[], double iSpare[]);
#endif

#endif /*SMP*/
