/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "diodefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/cpdefs.h" 

int
DIOaccept(CKTcircuit *ckt, GENmodel *inModel)
        /* make SOA checks after NR has finished */
{
    DIOmodel *model = (DIOmodel *) inModel;
    DIOinstance *here;
    int warn;       /* =1 SOA check should performed */
    double vd;      /* current diode voltage */

    if (!cp_getvar("warn", CP_NUM, &warn))
        warn = 0;

    if (warn > 0) {
        /* loop through all the models */
        for( ; model != NULL; model = model->DIOnextModel ) {
    
            /* loop through all the instances of the model */
            for (here = model->DIOinstances; here != NULL ;
                    here=here->DIOnextInstance) {
    
                vd = *(ckt->CKTrhsOld+here->DIOposPrimeNode)-
                        *(ckt->CKTrhsOld + here->DIOnegNode);
        
                if(!(ckt->CKTmode & (MODETRAN | MODETRANOP))) {
                    /* not transient, so shouldn't be here */
                    return(OK);
                } else {
                    if (vd > model->DIOfv_max)
                        printf("Instance: %s Model: %s Time: %g Vf=%g has exceeded Fv_max=%g\n", 
                        here->DIOname, model->DIOmodName, ckt->CKTtime, vd, model->DIOfv_max);
                    if (-vd > model->DIObv_max)
                        printf("Instance: %s Model: %s Time: %g |Vj|=%g has exceeded Bv_max=%g\n", 
                        here->DIOname, model->DIOmodName, ckt->CKTtime, -vd, model->DIObv_max);
                } // if ... else
    
            } // end instance loop
    
        } // end model loop

    }
    return(OK);
}
