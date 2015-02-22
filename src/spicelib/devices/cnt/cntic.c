/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/
/*
 */

#include "ngspice.h"
#include "cktdefs.h"
#include "cntdefs.h"
#include "sperror.h"
#include "suffix.h"


int
CNTgetic(GENmodel *inModel, CKTcircuit *ckt)
{
    CNTmodel *model = (CNTmodel *)inModel;
    CNTinstance *here;
    /*
     * grab initial conditions out of rhs array.   User specified, so use
     * external nodes to get values
     */

    for( ; model ; model = model->CNTnextModel) {
        for(here = model->CNTinstances; here ; here = here->CNTnextInstance) {
	    if (here->CNTowner != ARCHme) continue;
        
            if(!here->CNTicVBSGiven) {
                here->CNTicVBS = 
                        *(ckt->CKTrhs + here->CNTbNode) - 
                        *(ckt->CKTrhs + here->CNTsNode);
            }
            if(!here->CNTicVDSGiven) {
                here->CNTicVDS = 
                        *(ckt->CKTrhs + here->CNTdNode) - 
                        *(ckt->CKTrhs + here->CNTsNode);
            }
            if(!here->CNTicVGSGiven) {
                here->CNTicVGS = 
                        *(ckt->CKTrhs + here->CNTgNode) - 
                        *(ckt->CKTrhs + here->CNTsNode);
            }
        }
    }
    return(OK);
}
