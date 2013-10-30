/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4def.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/cpdefs.h"


int
BSIM4soaCheck(CKTcircuit *ckt, GENmodel *inModel)
{
    BSIM4model *model = (BSIM4model *) inModel;
    BSIM4instance *here;
    double vgs, vgd, vgb, vds, vbs, vbd;    /* actual mos voltages */
    int maxwarns_vgs = 0, maxwarns_vgd = 0, maxwarns_vgb = 0, maxwarns_vds = 0, maxwarns_vbs = 0, maxwarns_vbd = 0;
    static int warns_vgs = 0, warns_vgd = 0, warns_vgb = 0, warns_vds = 0, warns_vbs = 0, warns_vbd = 0;

    for (; model; model = model->BSIM4nextModel) {

        maxwarns_vgs = maxwarns_vgd = maxwarns_vgb = maxwarns_vds = maxwarns_vbs = maxwarns_vbd = ckt->CKTsoaMaxWarns;

        for (here = model->BSIM4instances; here; here = here->BSIM4nextInstance) {

            vgs = fabs(ckt->CKTrhsOld [here->BSIM4gNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4sNodePrime]);

            vgd = fabs(ckt->CKTrhsOld [here->BSIM4gNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4dNodePrime]);

            vgb = fabs(ckt->CKTrhsOld [here->BSIM4gNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4bNodePrime]);

            vds = fabs(ckt->CKTrhsOld [here->BSIM4dNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4sNodePrime]);

            vbs = fabs(ckt->CKTrhsOld [here->BSIM4bNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4sNodePrime]);

            vbd = fabs(ckt->CKTrhsOld [here->BSIM4bNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4dNodePrime]);

            if (vgs > model->BSIM4vgsMax)
                if (warns_vgs < maxwarns_vgs) {
                    soa_printf(ckt, (GENinstance*) here,
                               "|Vgs|=%g has exceeded Vgs_max=%g\n",
                               vgs, model->BSIM4vgsMax);
                    warns_vgs++;
                }

            if (vgd > model->BSIM4vgdMax)
                if (warns_vgd < maxwarns_vgd) {
                    soa_printf(ckt, (GENinstance*) here,
                               "|Vgd|=%g has exceeded Vgd_max=%g\n",
                               vgd, model->BSIM4vgdMax);
                    warns_vgd++;
                }

            if (vgb > model->BSIM4vgbMax)
                if (warns_vgb < maxwarns_vgb) {
                    soa_printf(ckt, (GENinstance*) here,
                               "|Vgb|=%g has exceeded Vgb_max=%g\n",
                               vgb, model->BSIM4vgbMax);
                    warns_vgb++;
                }

            if (vds > model->BSIM4vdsMax)
                if (warns_vds < maxwarns_vds) {
                    soa_printf(ckt, (GENinstance*) here,
                               "|Vds|=%g has exceeded Vds_max=%g\n",
                               vds, model->BSIM4vdsMax);
                    warns_vds++;
                }

            if (vbs > model->BSIM4vbsMax)
                if (warns_vbs < maxwarns_vbs) {
                    soa_printf(ckt, (GENinstance*) here,
                               "|Vbs|=%g has exceeded Vbs_max=%g\n",
                               vbs, model->BSIM4vbsMax);
                    warns_vbs++;
                }

            if (vbd > model->BSIM4vbdMax)
                if (warns_vbd < maxwarns_vbd) {
                    soa_printf(ckt, (GENinstance*) here,
                               "|Vbd|=%g has exceeded Vbd_max=%g\n",
                               vbd, model->BSIM4vbdMax);
                    warns_vbd++;
                }

        }
    }

    return OK;
}
