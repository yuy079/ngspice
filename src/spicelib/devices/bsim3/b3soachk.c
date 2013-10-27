/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim3def.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/cpdefs.h"


/* make SOA checks after NR has finished */

int
BSIM3soaCheck(CKTcircuit *ckt, GENmodel *inModel)
{
    BSIM3model *model = (BSIM3model *) inModel;
    BSIM3instance *here;
    double vgs, vgd, vgb, vds, vbs, vbd;    /* actual mos voltages */
    int maxwarns_vgs = 0, maxwarns_vgd = 0, maxwarns_vgb = 0, maxwarns_vds = 0, maxwarns_vbs = 0, maxwarns_vbd = 0;
    static int warns_vgs = 0, warns_vgd = 0, warns_vgb = 0, warns_vds = 0, warns_vbs = 0, warns_vbd = 0;

    if(!(ckt->CKTmode & (MODEDC | MODEDCOP | MODEDCTRANCURVE | MODETRAN | MODETRANOP)))
        return OK;

    for (; model; model = model->BSIM3nextModel) {

        maxwarns_vgs = maxwarns_vgd = maxwarns_vgb = maxwarns_vds = maxwarns_vbs = maxwarns_vbd = ckt->CKTsoaMaxWarns;

        for (here = model->BSIM3instances; here; here = here->BSIM3nextInstance) {

            vgs = fabs(ckt->CKTrhsOld [here->BSIM3gNode] -
                       ckt->CKTrhsOld [here->BSIM3sNodePrime]);

            vgd = fabs(ckt->CKTrhsOld [here->BSIM3gNode] -
                       ckt->CKTrhsOld [here->BSIM3dNodePrime]);

            vgb = fabs(ckt->CKTrhsOld [here->BSIM3gNode] -
                       ckt->CKTrhsOld [here->BSIM3bNode]);

            vds = fabs(ckt->CKTrhsOld [here->BSIM3dNodePrime] -
                       ckt->CKTrhsOld [here->BSIM3sNodePrime]);

            vbs = fabs(ckt->CKTrhsOld [here->BSIM3bNode] -
                       ckt->CKTrhsOld [here->BSIM3sNodePrime]);

            vbd = fabs(ckt->CKTrhsOld [here->BSIM3bNode] -
                       ckt->CKTrhsOld [here->BSIM3dNodePrime]);

            if (vgs > model->BSIM3vgsMax)
                if (warns_vgs < maxwarns_vgs) {
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vgs|=%g has exceeded Vgs_max=%g\n",
                               vgs, model->BSIM3vgsMax);
                    warns_vgs++;
                }

            if (vgd > model->BSIM3vgdMax)
                if (warns_vgd < maxwarns_vgd) {
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vgd|=%g has exceeded Vgd_max=%g\n",
                               vgd, model->BSIM3vgdMax);
                    warns_vgd++;
                }

            if (vgb > model->BSIM3vgbMax)
                if (warns_vgb < maxwarns_vgb) {
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vgb|=%g has exceeded Vgb_max=%g\n",
                               vgb, model->BSIM3vgbMax);
                    warns_vgb++;
                }

            if (vds > model->BSIM3vdsMax)
                if (warns_vds < maxwarns_vds) {
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vds|=%g has exceeded Vds_max=%g\n",
                               vds, model->BSIM3vdsMax);
                    warns_vds++;
                }

            if (vbs > model->BSIM3vbsMax)
                if (warns_vbs < maxwarns_vbs) {
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vbs|=%g has exceeded Vbs_max=%g\n",
                               vbs, model->BSIM3vbsMax);
                    warns_vbs++;
                }

            if (vbd > model->BSIM3vbdMax)
                if (warns_vbd < maxwarns_vbd) {
                    soa_printf(ckt, (GENinstance*) here, (GENmodel*) model,
                               "|Vbd|=%g has exceeded Vbd_max=%g\n",
                               vbd, model->BSIM3vbdMax);
                    warns_vbd++;
                }

        }
    }

    return OK;
}
