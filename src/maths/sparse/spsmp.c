/*
 *  Spice3 COMPATIBILITY MODULE
 *
 *  Author:                     Advising professor:
 *     Kenneth S. Kundert           Alberto Sangiovanni-Vincentelli
 *     UC Berkeley
 *
 *  This module contains routines that make Sparse1.3 a direct
 *  replacement for the SMP sparse matrix package in Spice3c1 or Spice3d1.
 *  Sparse1.3 is in general a faster and more robust package than SMP.
 *  These advantages become significant on large circuits.
 *
 *  >>> User accessible functions contained in this file:
 *  SMPaddElt
 *  SMPmakeElt
 *  SMPcClear
 *  SMPclear
 *  SMPcLUfac
 *  SMPluFac
 *  SMPcReorder
 *  SMPreorder
 *  SMPcaSolve
 *  SMPcSolve
 *  SMPsolve
 *  SMPmatSize
 *  SMPnewMatrix
 *  SMPdestroy
 *  SMPpreOrder
 *  SMPprint
 *  SMPgetError
 *  SMPcProdDiag
 *  LoadGmin
 *  SMPfindElt
 *  SMPcombine
 *  SMPcCombine
 */

/*
 *  To replace SMP with Sparse, rename the file spSpice3.h to
 *  spMatrix.h and place Sparse in a subdirectory of SPICE called
 *  `sparse'.  Then on UNIX compile Sparse by executing `make spice'.
 *  If not on UNIX, after compiling Sparse and creating the sparse.a
 *  archive, compile this file (spSMP.c) and spSMP.o to the archive,
 *  then copy sparse.a into the SPICE main directory and rename it
 *  SMP.a.  Finally link SPICE.
 *
 *  To be compatible with SPICE, the following Sparse compiler options
 *  (in spConfig.h) should be set as shown below:
 *
 *      EXPANDABLE                      YES
 *      TRANSLATE                       NO
 *      INITIALIZE                      NO or YES, YES for use with test prog.
 *      DIAGONAL_PIVOTING               YES
 *      MODIFIED_MARKOWITZ              NO
 *      DELETE                          NO
 *      STRIP                           NO
 *      MODIFIED_NODAL                  YES
 *      QUAD_ELEMENT                    NO
 *      TRANSPOSE                       YES
 *      SCALING                         NO
 *      DOCUMENTATION                   YES
 *      MULTIPLICATION                  NO
 *      DETERMINANT                     YES
 *      STABILITY                       NO
 *      CONDITION                       NO
 *      PSEUDOCONDITION                 NO
 *      DEBUG                           YES
 *
 *      spREAL  double
 */

/*
 *  Revision and copyright information.
 *
 *  Copyright (c) 1985,86,87,88,89,90
 *  by Kenneth S. Kundert and the University of California.
 *
 *  Permission to use, copy, modify, and distribute this software and its
 *  documentation for any purpose and without fee is hereby granted, provided
 *  that the above copyright notice appear in all copies and supporting
 *  documentation and that the authors and the University of California
 *  are properly credited.  The authors and the University of California
 *  make no representations as to the suitability of this software for
 *  any purpose.  It is provided `as is', without express or implied warranty.
 */

/*
 *  IMPORTS
 *
 *  >>> Import descriptions:
 *  spMatrix.h
 *     Sparse macros and declarations.
 *  SMPdefs.h
 *     Spice3's matrix macro definitions.
 */

#include "ngspice/config.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "ngspice/spmatrix.h"
#include "spdefs.h"
#include "ngspice/smpdefs.h"

#ifdef KLU
#include "ngspice/klu.h"
#endif

#if defined (_MSC_VER)
extern double scalbn(double, int);
#define logb _logb
extern double logb(double);
#endif

#ifdef KLU
static void LoadGmin_CSC (double **diag, int n, double Gmin) ;
#endif
static void LoadGmin(SMPmatrix *eMatrix, double Gmin);

#ifdef KLU
void SMPmatrix_CSC (SMPmatrix *Matrix, int **Ap, int **Ai, double **Ax, int n, double **Bind_Sparse, double **Bind_KLU, double **Diag) {
    spMatrix_CSC (Matrix, *Ap, *Ai, *Ax, n, Bind_Sparse, Bind_KLU, Diag) ;
    return ;
}
void SMPnnz (SMPmatrix *Matrix, int *CKTkluN, int *CKTklunz) {
    *CKTkluN = spGetSize (Matrix, 1) ;
    *CKTklunz = Matrix->Elements ;
    return ;
}
#endif

/*
 * SMPaddElt()
 */
int
SMPaddElt(SMPmatrix *Matrix, int Row, int Col, double Value)
{
    *spGetElement( Matrix, Row, Col ) = Value;
    return spError( Matrix );
}

/*
 * SMPmakeElt()
 */
double *
SMPmakeElt(SMPmatrix *Matrix, int Row, int Col)
{
    return spGetElement( Matrix, Row, Col );
}

/*
 * SMPcClear()
 */

#ifdef KLU
void
SMPcClear(SMPmatrix *Matrix, double *Ax, int CKTkluMODE)
{
    int i, nz ;
    if (CKTkluMODE) {
	spClear( Matrix ) ;
	if (Ax != NULL) {
	    nz = Matrix->Elements ;
    	    for (i = 0 ; i < 2 * nz ; i++) Ax [i] = 0 ;
	}
    }
    else spClear( Matrix ) ;
}
#else
void
SMPcClear(SMPmatrix *Matrix)
{
    spClear( Matrix );
}
#endif

/*
 * SMPclear()
 */

#ifdef KLU
void
SMPclear(SMPmatrix *Matrix, double *Ax, int CKTkluMODE)
{
    int i, nz ;
    if (CKTkluMODE) {
	spClear( Matrix ) ;
	if (Ax != NULL) {
	    nz = Matrix->Elements ;
    	    for (i = 0 ; i < nz ; i++) Ax [i] = 0 ;
	}
    }
    else spClear( Matrix ) ;
}
#else
void
SMPclear(SMPmatrix *Matrix)
{
    spClear( Matrix );
}
#endif

#define NG_IGNORE(x)  (void)x

/*
 * SMPcLUfac()
 */
/*ARGSUSED*/

#ifdef KLU
int
SMPcLUfac(SMPmatrix *Matrix, int *Ap, int *Ai, double *Ax, klu_symbolic *Symbolic, klu_numeric *Numeric, klu_common *Common, double PivTol, int CKTkluMODE)
{
    int ret ;
    if (CKTkluMODE) {
	NG_IGNORE(PivTol) ;

	spSetComplex( Matrix ) ;
	ret = klu_z_refactor (Ap, Ai, Ax, Symbolic, Numeric, Common) ;
	return (!ret) ;
    }
    else {
	NG_IGNORE(PivTol) ;

	spSetComplex( Matrix ) ;
	return spFactor( Matrix ) ;
    }
}
#else
int
SMPcLUfac(SMPmatrix *Matrix, double PivTol)
{
    NG_IGNORE(PivTol);

    spSetComplex( Matrix );
    return spFactor( Matrix );
}
#endif

/*
 * SMPluFac()
 */
/*ARGSUSED*/

#ifdef KLU
int
SMPluFac(SMPmatrix *Matrix, int *Ap, int *Ai, double *Ax, klu_symbolic *Symbolic, klu_numeric *Numeric, klu_common *Common, double **diag, double PivTol, double Gmin, int CKTkluMODE)
{
    int n, ret ;
    if (CKTkluMODE) {
	NG_IGNORE(PivTol) ;
	spSetReal( Matrix ) ;
	n = spGetSize (Matrix, 1) ;
	LoadGmin_CSC (diag, n, Gmin) ;
	ret = klu_refactor (Ap, Ai, Ax, Symbolic, Numeric, Common) ;
	return (!ret) ;
    }
    else {
	NG_IGNORE(PivTol) ;
	spSetReal( Matrix ) ;
	LoadGmin( Matrix, Gmin ) ;
	return spFactor( Matrix ) ;
    }
}
#else
int
SMPluFac(SMPmatrix *Matrix, double PivTol, double Gmin)
{
    NG_IGNORE(PivTol);
    spSetReal( Matrix );
    LoadGmin( Matrix, Gmin );
    return spFactor( Matrix );
}
#endif

/*
 * SMPcReorder()
 */

#ifdef KLU
int
SMPcReorder(SMPmatrix *Matrix, int *Ap, int *Ai, double *Ax, klu_symbolic **Symbolic, klu_numeric **Numeric, klu_common *Common, double PivTol, double PivRel, int *NumSwaps, int CKTkluMODE)
{
    if (CKTkluMODE) {
	*NumSwaps = 1 ;
	spSetComplex( Matrix ) ;
	klu_z_free_numeric (Numeric, Common) ;
	*Numeric = klu_z_factor (Ap, Ai, Ax, *Symbolic, Common) ;
	if (*Numeric == NULL) return 1 ;
	else return 0 ;
    }
    else {
	*NumSwaps = 1;
	spSetComplex( Matrix );
	return spOrderAndFactor( Matrix, NULL,
				 (spREAL)PivRel, (spREAL)PivTol, YES );
    }
}
#else
int
SMPcReorder(SMPmatrix *Matrix, double PivTol, double PivRel,
	    int *NumSwaps)
{
    *NumSwaps = 1;
    spSetComplex( Matrix );
    return spOrderAndFactor( Matrix, NULL,
                             (spREAL)PivRel, (spREAL)PivTol, YES );
}
#endif

/*
 * SMPreorder()
 */

#ifdef KLU
int
SMPreorder(SMPmatrix *Matrix, int *Ap, int *Ai, double *Ax, klu_symbolic *Symbolic, klu_numeric **Numeric, klu_common *Common, double **diag, double PivTol, double PivRel, double Gmin, int CKTkluMODE)
{
    int n ;
    if (CKTkluMODE) {
	spSetReal( Matrix );
	n = spGetSize (Matrix, 1) ;
	LoadGmin_CSC (diag, n, Gmin) ;
	if (*Numeric != NULL) {
	    klu_free_numeric (Numeric, Common) ;
	    *Numeric = klu_factor (Ap, Ai, Ax, Symbolic, Common) ;
	}
	else {
	    *Numeric = klu_factor (Ap, Ai, Ax, Symbolic, Common) ;
	}
	if (*Numeric == NULL) return 1 ;
	else return 0 ;
    }
    else {
	spSetReal( Matrix );
	LoadGmin( Matrix, Gmin );
	return spOrderAndFactor( Matrix, NULL,
				 (spREAL)PivRel, (spREAL)PivTol, YES );
    }
}
#else
int
SMPreorder(SMPmatrix *Matrix, double PivTol, double PivRel, double Gmin)
{
    spSetReal( Matrix );
    LoadGmin( Matrix, Gmin );
    return spOrderAndFactor( Matrix, NULL,
                             (spREAL)PivRel, (spREAL)PivTol, YES );
}
#endif

/*
 * SMPcaSolve()
 */
void
SMPcaSolve(SMPmatrix *Matrix, double RHS[], double iRHS[],
	   double Spare[], double iSpare[])
{
    printf ("SMPcaSolve\n") ;
    NG_IGNORE(iSpare);
    NG_IGNORE(Spare);

    spSolveTransposed( Matrix, RHS, RHS, iRHS, iRHS );
}

/*
 * SMPcSolve()
 */

#ifdef KLU
void
SMPcSolve(SMPmatrix *Matrix, klu_symbolic *Symbolic, klu_numeric *Numeric, klu_common *Common, double RHS[], double iRHS[], double Intermediate[], double Spare[], double iSpare[], int CKTkluMODE)
{
    int ret, n, i, *pExtOrder ;
    if (CKTkluMODE) {
	n = spGetSize (Matrix, 1) ;

	NG_IGNORE(iSpare);
	NG_IGNORE(Spare);

	pExtOrder = &Matrix->IntToExtRowMap[n];
	for (i = 2 * n - 1 ; i > 0 ; i -= 2) {
	    Intermediate [i] = RHS [*(pExtOrder)] ;
	    Intermediate [i - 1] = iRHS [*(pExtOrder--)] ;
	}

	ret = klu_z_solve (Symbolic, Numeric, n, 1, Intermediate, Common) ;

	pExtOrder = &Matrix->IntToExtColMap[n];
	for (i = 2 * n - 1 ; i > 0 ; i -= 2) {
	    RHS [*(pExtOrder)] = Intermediate [i] ;
	    iRHS [*(pExtOrder--)] = Intermediate [i - 1] ;
	}
    }
    else {
	NG_IGNORE(iSpare);
	NG_IGNORE(Spare);

	spSolve( Matrix, RHS, RHS, iRHS, iRHS );
    }
}
#else
void
SMPcSolve(SMPmatrix *Matrix, double RHS[], double iRHS[],
	  double Spare[], double iSpare[])
{
    NG_IGNORE(iSpare);
    NG_IGNORE(Spare);

    spSolve( Matrix, RHS, RHS, iRHS, iRHS );
}
#endif

/*
 * SMPsolve()
 */

#ifdef KLU
void
SMPsolve(SMPmatrix *Matrix, klu_symbolic *Symbolic, klu_numeric *Numeric, klu_common *Common, double RHS[], double Intermediate[], double Spare[], int CKTkluMODE)
{
    int ret, n, i, *pExtOrder ;
    if (CKTkluMODE) {
	n = spGetSize (Matrix, 1) ;

	NG_IGNORE(Spare);

	pExtOrder = &Matrix->IntToExtRowMap[n];
	for (i = n - 1 ; i >= 0 ; i--) Intermediate [i] = RHS [*(pExtOrder--)] ;

	ret = klu_solve (Symbolic, Numeric, n, 1, Intermediate, Common) ;

	pExtOrder = &Matrix->IntToExtColMap[n];
	for (i = n - 1 ; i >= 0 ; i--) RHS [*(pExtOrder--)] = Intermediate [i] ;
    }
    else {
	NG_IGNORE(Spare);

	spSolve( Matrix, RHS, RHS, NULL, NULL );
    }
}
#else
void
SMPsolve(SMPmatrix *Matrix, double RHS[], double Spare[])
{
    NG_IGNORE(Spare);

    spSolve( Matrix, RHS, RHS, NULL, NULL );
}
#endif

/*
 * SMPmatSize()
 */
int
SMPmatSize(SMPmatrix *Matrix)
{
    return spGetSize( Matrix, 1 );
}

/*
 * SMPnewMatrix()
 */
int
SMPnewMatrix(SMPmatrix **pMatrix)
{
    int Error;
    *pMatrix = spCreate( 0, 1, &Error );
    return Error;
}

/*
 * SMPdestroy()
 */

#ifdef KLU
void
SMPdestroy(SMPmatrix *Matrix, int **Ap, int **Ai, double **Ax, klu_symbolic **Symbolic, klu_numeric **Numeric, klu_common *Common, double ***bind_Sparse_Ptr, double ***bind_KLU_Ptr, double ***bind_KLU_Complex_Ptr, double ***diagPtr, double **Intermediate, double **Intermediate_Complex, int CKTkluMODE)
{
    if (CKTkluMODE) {
	printf("Destroy\n");
	spDestroy( Matrix );
	klu_free_numeric (Numeric, Common) ;
	klu_free_symbolic (Symbolic, Common) ;
	free (*Ap) ;
	free (*Ai) ;
	free (*Ax) ;
	free (*bind_Sparse_Ptr) ;
	free (*bind_KLU_Ptr) ;
	free (*bind_KLU_Complex_Ptr) ;
	free (*diagPtr) ;
	free (*Intermediate) ;
	free (*Intermediate_Complex) ;
    }
    else {
	spDestroy( Matrix );
    }
}
#else
void
SMPdestroy(SMPmatrix *Matrix)
{
    spDestroy( Matrix );
}
#endif

/*
 * SMPpreOrder()
 */

#ifdef KLU
int
SMPpreOrder(SMPmatrix *Matrix, int *Ap, int *Ai, klu_symbolic **Symbolic, klu_common *Common, int CKTkluMODE)
{
    int n ;
    if (CKTkluMODE) {
//	spMNA_Preorder( Matrix );

	n = spGetSize (Matrix, 1) ;
	*Symbolic = klu_analyze (n, Ap, Ai, Common) ;

	return 0 ;
    }
    else {
	spMNA_Preorder( Matrix );
	return spError( Matrix );
    }
}
#else
int
SMPpreOrder(SMPmatrix *Matrix)
{
    spMNA_Preorder( Matrix );
    return spError( Matrix );
}
#endif

/*
 * SMPprint()
 */
/*ARGSUSED*/
void
SMPprint(SMPmatrix *Matrix, FILE *File)
{
    NG_IGNORE(File);

    spPrint( Matrix, 0, 1, 1 );
}

/*
 * SMPgetError()
 */
void
SMPgetError(SMPmatrix *Matrix, int *Col, int *Row)
{
    spWhereSingular( Matrix, Row, Col );
}

/*
 * SMPcProdDiag()
 *    note: obsolete for Spice3d2 and later
 */
int
SMPcProdDiag(SMPmatrix *Matrix, SPcomplex *pMantissa, int *pExponent)
{
    spDeterminant( Matrix, pExponent, &(pMantissa->real),
                                              &(pMantissa->imag) );
    return spError( Matrix );
}

/*
 * SMPcDProd()
 */
int
SMPcDProd(SMPmatrix *Matrix, SPcomplex *pMantissa, int *pExponent)
{
    double	re, im, x, y, z;
    int		p;

    spDeterminant( Matrix, &p, &re, &im);

#ifndef M_LN2
#define M_LN2   0.69314718055994530942
#endif
#ifndef M_LN10
#define M_LN10  2.30258509299404568402
#endif

#ifdef debug_print
    printf("Determinant 10: (%20g,%20g)^%d\n", re, im, p);
#endif

    /* Convert base 10 numbers to base 2 numbers, for comparison */
    y = p * M_LN10 / M_LN2;
    x = (int) y;
    y -= x;

    /* ASSERT
     *	x = integral part of exponent, y = fraction part of exponent
     */

    /* Fold in the fractional part */
#ifdef debug_print
    printf(" ** base10 -> base2 int =  %g, frac = %20g\n", x, y);
#endif
    z = pow(2.0, y);
    re *= z;
    im *= z;
#ifdef debug_print
    printf(" ** multiplier = %20g\n", z);
#endif

    /* Re-normalize (re or im may be > 2.0 or both < 1.0 */
    if (re != 0.0) {
	y = logb(re);
	if (im != 0.0)
	    z = logb(im);
	else
	    z = 0;
    } else if (im != 0.0) {
	z = logb(im);
	y = 0;
    } else {
	/* Singular */
	/*printf("10 -> singular\n");*/
	y = 0;
	z = 0;
    }

#ifdef debug_print
    printf(" ** renormalize changes = %g,%g\n", y, z);
#endif
    if (y < z)
	y = z;

    *pExponent = (int)(x + y);
    x = scalbn(re, (int) -y);
    z = scalbn(im, (int) -y);
#ifdef debug_print
    printf(" ** values are: re %g, im %g, y %g, re' %g, im' %g\n",
	    re, im, y, x, z);
#endif
    pMantissa->real = scalbn(re, (int) -y);
    pMantissa->imag = scalbn(im, (int) -y);

#ifdef debug_print
    printf("Determinant 10->2: (%20g,%20g)^%d\n", pMantissa->real,
	pMantissa->imag, *pExponent);
#endif
    return spError( Matrix );
}



/*
 *  The following routines need internal knowledge of the Sparse data
 *  structures.
 */

/*
 *  LOAD GMIN
 *
 *  This routine adds Gmin to each diagonal element.  Because Gmin is
 *  added to the current diagonal, which may bear little relation to
 *  what the outside world thinks is a diagonal, and because the
 *  elements that are diagonals may change after calling spOrderAndFactor,
 *  use of this routine is not recommended.  It is included here simply
 *  for compatibility with Spice3.
 */


/* Francesco Lannutti */
#ifdef KLU
static void LoadGmin_CSC (double **diag, int n, double Gmin) {

	int i ;
	if (Gmin != 0.0) {
		for (i = 0 ; i < n ; i++) {
			if (diag [i] != NULL) *(diag [i]) += Gmin ;
		}
	}

	return ;
}
#endif
static void
LoadGmin(SMPmatrix *eMatrix, double Gmin)
{
    MatrixPtr Matrix = eMatrix;
    int I;
    ArrayOfElementPtrs Diag;
    ElementPtr diag;

    /* Begin `LoadGmin'. */
    assert( IS_SPARSE( Matrix ) );

    if (Gmin != 0.0) {
	Diag = Matrix->Diag;
	for (I = Matrix->Size; I > 0; I--) {
	    if ((diag = Diag[I]) != NULL)
		diag->Real += Gmin;
	}
    }
    return;
}




/*
 *  FIND ELEMENT
 *
 *  This routine finds an element in the matrix by row and column number.
 *  If the element exists, a pointer to it is returned.  If not, then NULL
 *  is returned unless the CreateIfMissing flag is TRUE, in which case a
 *  pointer to the new element is returned.
 */

SMPelement *
SMPfindElt(SMPmatrix *eMatrix, int Row, int Col, int CreateIfMissing)
{
    MatrixPtr Matrix = eMatrix;
    ElementPtr Element;

    /* Begin `SMPfindElt'. */
    assert( IS_SPARSE( Matrix ) );
    Row = Matrix->ExtToIntRowMap[Row];
    Col = Matrix->ExtToIntColMap[Col];
    Element = Matrix->FirstInCol[Col];
    Element = spcFindElementInCol(Matrix, &Element, Row, Col, CreateIfMissing);
    return (SMPelement *)Element;
}

/* XXX The following should probably be implemented in spUtils */

/*
 * SMPcZeroCol()
 */
int
SMPcZeroCol(SMPmatrix *eMatrix, int Col)
{
    MatrixPtr Matrix = eMatrix;
    ElementPtr	Element;

    Col = Matrix->ExtToIntColMap[Col];

    for (Element = Matrix->FirstInCol[Col];
	Element != NULL;
	Element = Element->NextInCol)
    {
	Element->Real = 0.0;
	Element->Imag = 0.0;
    }

    return spError( Matrix );
}

/*
 * SMPcAddCol()
 */
int
SMPcAddCol(SMPmatrix *eMatrix, int Accum_Col, int Addend_Col)
{
    MatrixPtr Matrix = eMatrix;
    ElementPtr	Accum, Addend, *Prev;

    Accum_Col = Matrix->ExtToIntColMap[Accum_Col];
    Addend_Col = Matrix->ExtToIntColMap[Addend_Col];

    Addend = Matrix->FirstInCol[Addend_Col];
    Prev = &Matrix->FirstInCol[Accum_Col];
    Accum = *Prev;;

    while (Addend != NULL) {
	while (Accum && Accum->Row < Addend->Row) {
	    Prev = &Accum->NextInCol;
	    Accum = *Prev;
	}
	if (!Accum || Accum->Row > Addend->Row) {
	    Accum = spcCreateElement(Matrix, Addend->Row, Accum_Col, Prev, 0);
	}
	Accum->Real += Addend->Real;
	Accum->Imag += Addend->Imag;
	Addend = Addend->NextInCol;
    }

    return spError( Matrix );
}

/*
 * SMPzeroRow()
 */
int
SMPzeroRow(SMPmatrix *eMatrix, int Row)
{
    MatrixPtr Matrix = eMatrix;
    ElementPtr	Element;

    Row = Matrix->ExtToIntColMap[Row];

    if (Matrix->RowsLinked == NO)
	spcLinkRows(Matrix);

    if (Matrix->PreviousMatrixWasComplex || Matrix->Complex) {
	for (Element = Matrix->FirstInRow[Row];
	    Element != NULL;
	    Element = Element->NextInRow)
	{
	    Element->Real = 0.0;
	    Element->Imag = 0.0;
	}
    } else {
	for (Element = Matrix->FirstInRow[Row];
	    Element != NULL;
	    Element = Element->NextInRow)
	{
	    Element->Real = 0.0;
	}
    }

    return spError( Matrix );
}

#ifdef PARALLEL_ARCH
/*
 * SMPcombine()
 */
void
SMPcombine(SMPmatrix *Matrix, double RHS[], double Spare[])
{
    spSetReal( Matrix );
    spCombine( Matrix, RHS, Spare, NULL, NULL );
}

/*
 * SMPcCombine()
 */
void
SMPcCombine(SMPmatrix *Matrix, double RHS[], double Spare[],
	    double iRHS[], double iSpare[])
{
    spSetComplex( Matrix );
    spCombine( Matrix, RHS, Spare, iRHS, iSpare );
}
#endif /* PARALLEL_ARCH */
