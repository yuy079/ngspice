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


/* make SOA checks after NR has finished */

int
BSIM3v32soaCheck(CKTcircuit *ckt, GENmodel *inModel)
{
    BSIM3v32model *model = (BSIM3v32model *) inModel;
    BSIM3v32instance *here;
    double vgs, vgd, vgb, vds, vbs, vbd;    /* actual mos voltages */
    int maxwarns_vgs = 0, maxwarns_vgd = 0, maxwarns_vgb = 0, maxwarns_vds = 0, maxwarns_vbs = 0, maxwarns_vbd = 0;
    static int warns_vgs = 0, warns_vgd = 0, warns_vgb = 0, warns_vds = 0, warns_vbs = 0, warns_vbd = 0;

    if(!(ckt->CKTmode & (MODEDC | MODEDCOP | MODEDCTRANCURVE | MODETRAN | MODETRANOP)))
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
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vgs|=%g has exceeded Vgs_max=%g\n",
                               vgs, model->BSIM3v32vgsMax);
                    warns_vgs++;
                }

            if (vgd > model->BSIM3v32vgdMax)
                if (warns_vgd < maxwarns_vgd) {
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vgd|=%g has exceeded Vgd_max=%g\n",
                               vgd, model->BSIM3v32vgdMax);
                    warns_vgd++;
                }

            if (vgb > model->BSIM3v32vgbMax)
                if (warns_vgb < maxwarns_vgb) {
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vgb|=%g has exceeded Vgb_max=%g\n",
                               vgb, model->BSIM3v32vgbMax);
                    warns_vgb++;
                }

            if (vds > model->BSIM3v32vdsMax)
                if (warns_vds < maxwarns_vds) {
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vds|=%g has exceeded Vds_max=%g\n",
                               vds, model->BSIM3v32vdsMax);
                    warns_vds++;
                }

            if (vbs > model->BSIM3v32vbsMax)
                if (warns_vbs < maxwarns_vbs) {
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vbs|=%g has exceeded Vbs_max=%g\n",
                               vbs, model->BSIM3v32vbsMax);
                    warns_vbs++;
                }

            if (vbd > model->BSIM3v32vbdMax)
                if (warns_vbd < maxwarns_vbd) {
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vbd|=%g has exceeded Vbd_max=%g\n",
                               vbd, model->BSIM3v32vbdMax);
                    warns_vbd++;
                }

        }
    }

    return OK;
}
