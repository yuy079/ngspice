/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/
/*
 */

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "cntdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"


int
CNTtrunc(GENmodel *inModel, CKTcircuit *ckt, double *timeStep)
{
    CNTmodel *model = (CNTmodel *)inModel;
    CNTinstance *here;

    for( ; model != NULL; model = model->CNTnextModel) {
        for(here=model->CNTinstances;here!=NULL;here = here->CNTnextInstance){
/*	    if (here->CNTowner != ARCHme) continue;*/
        
            CKTterr(here->CNTqgs,ckt,timeStep);
            CKTterr(here->CNTqgd,ckt,timeStep);
            CKTterr(here->CNTqgb,ckt,timeStep);
        }
    }
    return(OK);
}
