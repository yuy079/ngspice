/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified 1999 Emmanuel Rouat
**********/

    /*
     * NIacIter(ckt)
     *
     *  This subroutine performs the actual numerical iteration.
     *  It uses the sparse matrix stored in the NIstruct by NIinit,
     *  along with the matrix loading program, the load data, the
     *  convergence test function, and the convergence parameters
     * - return value is non-zero for convergence failure 
     */

#include "ngspice/ngspice.h"
#include "ngspice/trandefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/sperror.h"
#include "niaciter.h"


int
NIacIter(CKTcircuit *ckt)
{
    int error;
    int ignore;
    double *temp;
    double startTime;

retry:
    ckt->CKTnoncon=0;

    error = CKTacLoad(ckt);
    if(error) return(error);
    
    if(ckt->CKTniState & NIACSHOULDREORDER) {
	startTime = SPfrontEnd->IFseconds();

#ifdef KLU
        error = SMPcReorder(ckt->CKTmatrix, ckt->CKTkluAp, ckt->CKTkluAi, ckt->CKTkluAx,
                &(ckt->CKTkluSymbolic), &(ckt->CKTkluNumeric), ckt->CKTkluCommon,
                ckt->CKTpivotAbsTol, ckt->CKTpivotRelTol, &ignore, ckt->CKTkluMODE);
#else
        error = SMPcReorder(ckt->CKTmatrix,ckt->CKTpivotAbsTol,
                ckt->CKTpivotRelTol,&ignore);
#endif

	ckt->CKTstat->STATreorderTime +=
		SPfrontEnd->IFseconds()- startTime;
        ckt->CKTniState &= ~NIACSHOULDREORDER;
        if(error != 0) {
            /* either singular equations or no memory, in either case,
             * let caller handle problem
             */
            return(error);
        }
    } else {
	startTime = SPfrontEnd->IFseconds();

#ifdef KLU
        error = SMPcLUfac(ckt->CKTmatrix, ckt->CKTkluAp, ckt->CKTkluAi, ckt->CKTkluAx,
                          ckt->CKTkluSymbolic, ckt->CKTkluNumeric, ckt->CKTkluCommon,
                          ckt->CKTpivotAbsTol, ckt->CKTkluMODE);
#else
        error = SMPcLUfac(ckt->CKTmatrix,ckt->CKTpivotAbsTol);
#endif

	ckt->CKTstat->STATdecompTime += 
		SPfrontEnd->IFseconds()-startTime;
        if(error != 0) {
            if(error == E_SINGULAR) {
                /* the problem is that the matrix can't be solved with the
                 * current LU factorization.  Maybe if we reload and
                 * try to reorder again it will help...
                 */
                ckt->CKTniState |= NIACSHOULDREORDER;
                goto retry;
            }
            return(error); /* can't handle E_BADMATRIX, so let caller */
        }
    } 
    startTime = SPfrontEnd->IFseconds();

#ifdef KLU
    SMPcSolve(ckt->CKTmatrix, ckt->CKTkluSymbolic, ckt->CKTkluNumeric, ckt->CKTkluCommon,
              ckt->CKTrhs, ckt->CKTirhs, ckt->CKTkluIntermediate_Complex, ckt->CKTrhsSpare, ckt->CKTirhsSpare, ckt->CKTkluMODE);
#else
    SMPcSolve(ckt->CKTmatrix,ckt->CKTrhs, 
            ckt->CKTirhs, ckt->CKTrhsSpare,
            ckt->CKTirhsSpare);
#endif

    ckt->CKTstat->STATsolveTime += SPfrontEnd->IFseconds() - startTime;

    *ckt->CKTrhs = 0;
    *ckt->CKTrhsSpare = 0;
    *ckt->CKTrhsOld = 0;
    *ckt->CKTirhs = 0;
    *ckt->CKTirhsSpare = 0;
    *ckt->CKTirhsOld = 0;

    temp = ckt->CKTirhsOld;
    ckt->CKTirhsOld = ckt->CKTirhs;
    ckt->CKTirhs = temp;

    temp = ckt->CKTrhsOld;
    ckt->CKTrhsOld = ckt->CKTrhs;
    ckt->CKTrhs = temp;

    return(OK);
}
