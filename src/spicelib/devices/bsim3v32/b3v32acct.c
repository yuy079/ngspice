/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim3v32def.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/cpdefs.h"


extern FILE *slogp;  /* soa log file ('--soa-log file' command line option) */


/* make SOA checks after NR has finished */

int
BSIM3v32accept(CKTcircuit *ckt, GENmodel *inModel)
{
    BSIM3v32model *model = (BSIM3v32model *) inModel;
    BSIM3v32instance *here;
    double vgs, vgd, vgb, vds, vbs, vbd;    /* actual mos voltages */
    int maxwarns_vgs = 0, maxwarns_vgd = 0, maxwarns_vgb = 0, maxwarns_vds = 0, maxwarns_vbs = 0, maxwarns_vbd = 0;
    static int warns_vgs = 0, warns_vgd = 0, warns_vgb = 0, warns_vds = 0, warns_vbs = 0, warns_vbd = 0;

    if (!ckt->CKTsoaCheck)
        return OK;

    if (!(ckt->CKTmode & (MODETRAN | MODETRANOP)))
        return OK;

    for (; model; model = model->BSIM3v32nextModel) {

        maxwarns_vgs = maxwarns_vgd = maxwarns_vgb = maxwarns_vds = maxwarns_vbs = maxwarns_vbd = ckt->CKTsoaMaxWarns;

        for (here = model->BSIM3v32instances; here; here = here->BSIM3v32nextInstance) {

            vgs = fabs(ckt->CKTrhsOld [here->BSIM3v32gNode] -
                       ckt->CKTrhsOld [here->BSIM3v32sNodePrime]);

            vgd = fabs(ckt->CKTrhsOld [here->BSIM3v32gNode] -
                       ckt->CKTrhsOld [here->BSIM3v32dNodePrime]);

            vgb = fabs(ckt->CKTrhsOld [here->BSIM3v32gNode] -
                       ckt->CKTrhsOld [here->BSIM3v32bNode]);

            vds = fabs(ckt->CKTrhsOld [here->BSIM3v32dNodePrime] -
                       ckt->CKTrhsOld [here->BSIM3v32sNodePrime]);

            vbs = fabs(ckt->CKTrhsOld [here->BSIM3v32bNode] -
                       ckt->CKTrhsOld [here->BSIM3v32sNodePrime]);

            vbd = fabs(ckt->CKTrhsOld [here->BSIM3v32bNode] -
                       ckt->CKTrhsOld [here->BSIM3v32dNodePrime]);

            if (vgs > model->BSIM3v32vgsMax)
                if (warns_vgs < maxwarns_vgs) {
                    printf("Instance: %s Model: %s Time: %g |Vgs|=%g has exceeded Vgs_max=%g\n",
                           here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vgs, model->BSIM3v32vgsMax);
                    if (slogp)
                        fprintf(slogp, "Instance: %s Model: %s Time: %g |Vgs|=%g has exceeded Vgs_max=%g\n",
                                here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vgs, model->BSIM3v32vgsMax);
                    warns_vgs++;
                }

            if (vgd > model->BSIM3v32vgdMax)
                if (warns_vgd < maxwarns_vgd) {
                    printf("Instance: %s Model: %s Time: %g |Vgd|=%g has exceeded Vgd_max=%g\n",
                           here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vgd, model->BSIM3v32vgdMax);
                    if (slogp)
                        fprintf(slogp, "Instance: %s Model: %s Time: %g |Vgd|=%g has exceeded Vgd_max=%g\n",
                                here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vgd, model->BSIM3v32vgdMax);
                    warns_vgd++;
                }

            if (vgb > model->BSIM3v32vgbMax)
                if (warns_vgb < maxwarns_vgb) {
                    printf("Instance: %s Model: %s Time: %g |Vgb|=%g has exceeded Vgb_max=%g\n",
                           here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vgb, model->BSIM3v32vgbMax);
                    if (slogp)
                        fprintf(slogp, "Instance: %s Model: %s Time: %g |Vgb|=%g has exceeded Vgb_max=%g\n",
                                here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vgb, model->BSIM3v32vgbMax);
                    warns_vgb++;
                }

            if (vds > model->BSIM3v32vdsMax)
                if (warns_vds < maxwarns_vds) {
                    printf("Instance: %s Model: %s Time: %g |Vds|=%g has exceeded Vds_max=%g\n",
                           here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vds, model->BSIM3v32vdsMax);
                    if (slogp)
                        fprintf(slogp, "Instance: %s Model: %s Time: %g |Vds|=%g has exceeded Vds_max=%g\n",
                                here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vds, model->BSIM3v32vdsMax);
                    warns_vds++;
                }

            if (vbs > model->BSIM3v32vbsMax)
                if (warns_vbs < maxwarns_vbs) {
                    printf("Instance: %s Model: %s Time: %g |Vbs|=%g has exceeded Vbs_max=%g\n",
                           here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vbs, model->BSIM3v32vbsMax);
                    if (slogp)
                        fprintf(slogp, "Instance: %s Model: %s Time: %g |Vbs|=%g has exceeded Vbs_max=%g\n",
                                here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vbs, model->BSIM3v32vbsMax);
                    warns_vbs++;
                }

            if (vbd > model->BSIM3v32vbdMax)
                if (warns_vbd < maxwarns_vbd) {
                    printf("Instance: %s Model: %s Time: %g |Vbd|=%g has exceeded Vbd_max=%g\n",
                           here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vbd, model->BSIM3v32vbdMax);
                    if (slogp)
                        fprintf(slogp, "Instance: %s Model: %s Time: %g |Vbd|=%g has exceeded Vbd_max=%g\n",
                                here->BSIM3v32name, model->BSIM3v32modName, ckt->CKTtime, vbd, model->BSIM3v32vbdMax);
                    warns_vbd++;
                }

        }
    }

    return OK;
}
