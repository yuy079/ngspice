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
    double vd;         /* actual diode voltage */
    int maxwarns_fv=0, maxwarns_bv=0;
    static int warns_fv=0, warns_bv=0;
   
    if (ckt->CKTsoaCheck > 0) {

        /* loop through all the models */
        for( ; model != NULL; model = model->DIOnextModel ) {
    
            maxwarns_fv = maxwarns_bv = ckt->CKTsoaMaxWarns;

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
                        if (warns_fv < maxwarns_fv) {
                            printf("Instance: %s Model: %s Time: %g Vf=%g has exceeded Fv_max=%g\n", 
                            here->DIOname, model->DIOmodName, ckt->CKTtime, vd, model->DIOfv_max);
                            warns_fv++;
                        }
                    if (-vd > model->DIObv_max)
                        if (warns_bv < maxwarns_bv) {
                            printf("Instance: %s Model: %s Time: %g |Vj|=%g has exceeded Bv_max=%g\n", 
                            here->DIOname, model->DIOmodName, ckt->CKTtime, -vd, model->DIObv_max);
                            warns_bv++;
                        }
                } // if ... else
    
            } // end instance loop
    
        } // end model loop

    }
    return(OK);
}
