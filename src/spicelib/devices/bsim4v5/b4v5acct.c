/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4v5def.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/cpdefs.h" 

extern FILE *slogp;          /* soa log file ('--soa-log file' command line option) */

int
BSIM4v5accept(CKTcircuit *ckt, GENmodel *inModel)
/* make SOA checks after NR has finished */
{
    BSIM4v5model *model = (BSIM4v5model *) inModel;
    BSIM4v5instance *here;
    double vgs, vgd, vgb, vds, vbs, vbd;    /* actual bjt voltages */
    int maxwarns_vgs=0, maxwarns_vgd=0, maxwarns_vgb=0, maxwarns_vds=0, maxwarns_vbs=0, maxwarns_vbd=0;
    static int warns_vgs=0, warns_vgd=0, warns_vgb=0, warns_vds=0, warns_vbs=0, warns_vbd=0;

    if (ckt->CKTsoaCheck > 0) {

        /* loop through all the models */
        for( ; model != NULL; model = model->BSIM4v5nextModel ) {
    
            maxwarns_vgs = maxwarns_vgd = maxwarns_vgb = maxwarns_vds = maxwarns_vbs = maxwarns_vbd = ckt->CKTsoaMaxWarns;

            /* loop through all the instances of the model */
            for (here = model->BSIM4v5instances; here != NULL ;
                    here=here->BSIM4v5nextInstance) {

                vgs = fabs(*(ckt->CKTrhsOld + here->BSIM4v5gNodePrime) 
                           - *(ckt->CKTrhsOld + here->BSIM4v5sNodePrime));
                vgd = fabs(*(ckt->CKTrhsOld + here->BSIM4v5gNodePrime) 
                           - *(ckt->CKTrhsOld + here->BSIM4v5dNodePrime));
                vgb = fabs(*(ckt->CKTrhsOld + here->BSIM4v5gNodePrime) 
                           - *(ckt->CKTrhsOld + here->BSIM4v5bNodePrime));
                vds = fabs(*(ckt->CKTrhsOld + here->BSIM4v5dNodePrime)
                           - *(ckt->CKTrhsOld + here->BSIM4v5sNodePrime));
                vbs = fabs(*(ckt->CKTrhsOld + here->BSIM4v5bNodePrime)
                           - *(ckt->CKTrhsOld + here->BSIM4v5sNodePrime));
                vbd = fabs(*(ckt->CKTrhsOld + here->BSIM4v5bNodePrime)
                           - *(ckt->CKTrhsOld + here->BSIM4v5dNodePrime));

                if(!(ckt->CKTmode & (MODETRAN | MODETRANOP))) {
                    /* not transient, so shouldn't be here */
                    return(OK);
                } else {
                    if (vgs > model->BSIM4v5vgsMax)
                        if (warns_vgs < maxwarns_vgs) {
                            printf("Instance: %s Model: %s Time: %g |Vgs|=%g has exceeded Vgs_max=%g\n", 
                            here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vgs, model->BSIM4v5vgsMax);
                            if (slogp)
                                fprintf(slogp, "Instance: %s Model: %s Time: %g |Vgs|=%g has exceeded Vgs_max=%g\n", 
                                here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vgs, model->BSIM4v5vgsMax);
                            warns_vgs++;
                        }
                    if (vgd > model->BSIM4v5vgdMax)
                        if (warns_vgd < maxwarns_vgd) {
                            printf("Instance: %s Model: %s Time: %g |Vgd|=%g has exceeded Vgd_max=%g\n", 
                            here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vgd, model->BSIM4v5vgdMax);
                            if (slogp)
                                fprintf(slogp, "Instance: %s Model: %s Time: %g |Vgd|=%g has exceeded Vgd_max=%g\n", 
                                here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vgd, model->BSIM4v5vgdMax);
                            warns_vgd++;
                        }
                    if (vgb > model->BSIM4v5vgbMax)
                        if (warns_vgb < maxwarns_vgb) {
                            printf("Instance: %s Model: %s Time: %g |Vgb|=%g has exceeded Vgb_max=%g\n", 
                            here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vgb, model->BSIM4v5vgbMax);
                            if (slogp)
                                fprintf(slogp, "Instance: %s Model: %s Time: %g |Vgb|=%g has exceeded Vgb_max=%g\n", 
                                here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vgb, model->BSIM4v5vgbMax);
                            warns_vgb++;
                        }
                    if (vds > model->BSIM4v5vdsMax)
                        if (warns_vds < maxwarns_vds) {
                            printf("Instance: %s Model: %s Time: %g |Vds|=%g has exceeded Vds_max=%g\n", 
                            here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vds, model->BSIM4v5vdsMax);
                            if (slogp)
                                fprintf(slogp, "Instance: %s Model: %s Time: %g |Vds|=%g has exceeded Vds_max=%g\n", 
                                here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vds, model->BSIM4v5vdsMax);
                            warns_vds++;
                        }
                    if (vbs > model->BSIM4v5vbsMax)
                        if (warns_vbs < maxwarns_vbs) {
                            printf("Instance: %s Model: %s Time: %g |Vbs|=%g has exceeded Vbs_max=%g\n", 
                            here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vbs, model->BSIM4v5vbsMax);
                            if (slogp)
                                fprintf(slogp, "Instance: %s Model: %s Time: %g |Vbs|=%g has exceeded Vbs_max=%g\n", 
                                here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vbs, model->BSIM4v5vbsMax);
                            warns_vbs++;
                        }
                    if (vbd > model->BSIM4v5vbdMax)
                        if (warns_vbd < maxwarns_vbd) {
                            printf("Instance: %s Model: %s Time: %g |Vbd|=%g has exceeded Vbd_max=%g\n", 
                            here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vbd, model->BSIM4v5vbdMax);
                            if (slogp)
                                fprintf(slogp, "Instance: %s Model: %s Time: %g |Vbd|=%g has exceeded Vbd_max=%g\n", 
                                here->BSIM4v5name, model->BSIM4v5modName, ckt->CKTtime, vbd, model->BSIM4v5vbdMax);
                            warns_vbd++;
                        }
                } // if ... else
    
            } // end instance loop
    
        } // end model loop

    }
    return(OK);
}
