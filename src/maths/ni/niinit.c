/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/
/*
 */
    /*
     * NIinit(nistruct,loadpkg)
     *
     *  Initialize the Numerical iteration package to perform Newton-Raphson
     *  iterations on a sparse matrix filled by the specified load package,
     */

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/smpdefs.h"

#ifdef KLU
#include "ngspice/klu.h"
#endif

#include "niinit.h"


int
NIinit(CKTcircuit *ckt)
{
#ifdef SPARSE
/* a concession to Ken Kundert's sparse matrix package - SMP doesn't need this*/
    int Error;
#endif /* SPARSE */

    #ifdef KLU
    ckt->CKTkluCommon = TMALLOC(klu_common, 1);
    ckt->CKTkluSymbolic = NULL ;
    ckt->CKTkluNumeric = NULL ;
    ckt->CKTkluAp = NULL ;
    ckt->CKTkluAi = NULL ;
    ckt->CKTkluAx = NULL ;
    ckt->CKTkluIntermediate = NULL ;
    ckt->CKTkluIntermediate_Complex = NULL ;
    ckt->CKTkluBind_Sparse = NULL ;
    ckt->CKTkluBind_KLU = NULL ;
    ckt->CKTkluBind_KLU_Complex = NULL ;
    ckt->CKTkluDiag = NULL ;
    ckt->CKTkluN = 0 ;
    ckt->CKTklunz = 0 ;
    ckt->CKTkluMODE = CKTkluON ; /* TO BE SUBSTITUTED WITH THE HEURISTICS */

    klu_defaults (ckt->CKTkluCommon) ;
    #endif

    ckt->CKTniState = NIUNINITIALIZED;
    return(SMPnewMatrix( &(ckt->CKTmatrix) ) );
}
