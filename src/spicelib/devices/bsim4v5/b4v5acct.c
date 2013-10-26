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

void
soa_printf(GENinstance *, GENmodel *, CKTcircuit *, const char *, ...);

/* make SOA checks after NR has finished */

int
BSIM4v5accept(CKTcircuit *ckt, GENmodel *inModel)
{
    BSIM4v5model *model = (BSIM4v5model *) inModel;
    BSIM4v5instance *here;
    double vgs, vgd, vgb, vds, vbs, vbd;    /* actual mos voltages */
    int maxwarns_vgs = 0, maxwarns_vgd = 0, maxwarns_vgb = 0, maxwarns_vds = 0, maxwarns_vbs = 0, maxwarns_vbd = 0;
    static int warns_vgs = 0, warns_vgd = 0, warns_vgb = 0, warns_vds = 0, warns_vbs = 0, warns_vbd = 0;

    if (!ckt->CKTsoaCheck)
        return OK;

    if(!(ckt->CKTmode & (MODEDC | MODEDCOP | MODEDCTRANCURVE | MODETRAN | MODETRANOP)))
        return OK;

    for (; model; model = model->BSIM4v5nextModel) {

        maxwarns_vgs = maxwarns_vgd = maxwarns_vgb = maxwarns_vds = maxwarns_vbs = maxwarns_vbd = ckt->CKTsoaMaxWarns;

        for (here = model->BSIM4v5instances; here; here = here->BSIM4v5nextInstance) {

            vgs = fabs(ckt->CKTrhsOld [here->BSIM4v5gNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4v5sNodePrime]);

            vgd = fabs(ckt->CKTrhsOld [here->BSIM4v5gNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4v5dNodePrime]);

            vgb = fabs(ckt->CKTrhsOld [here->BSIM4v5gNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4v5bNodePrime]);

            vds = fabs(ckt->CKTrhsOld [here->BSIM4v5dNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4v5sNodePrime]);

            vbs = fabs(ckt->CKTrhsOld [here->BSIM4v5bNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4v5sNodePrime]);

            vbd = fabs(ckt->CKTrhsOld [here->BSIM4v5bNodePrime] -
                       ckt->CKTrhsOld [here->BSIM4v5dNodePrime]);

            if (vgs > model->BSIM4v5vgsMax)
                if (warns_vgs < maxwarns_vgs) {
                    soa_printf((GENinstance*) here, (GENmodel*) model, ckt,
                               "|Vgs|=%g has exceeded Vgs_max=%g\n",
                               vgs, model->BSIM4v5vgsMax);
                    warns_vgs++;
                }

            if (vgd > model->BSIM4v5vgdMax)
                if (warns_vgd < maxwarns_vgd) {
                    soa_printf((GENinstance*) here, (GENmodel*) model, ckt,
                               "|Vgd|=%g has exceeded Vgd_max=%g\n",
                               vgd, model->BSIM4v5vgdMax);
                    warns_vgd++;
                }

            if (vgb > model->BSIM4v5vgbMax)
                if (warns_vgb < maxwarns_vgb) {
                    soa_printf((GENinstance*) here, (GENmodel*) model, ckt,
                               "|Vgb|=%g has exceeded Vgb_max=%g\n",
                               vgb, model->BSIM4v5vgbMax);
                    warns_vgb++;
                }

            if (vds > model->BSIM4v5vdsMax)
                if (warns_vds < maxwarns_vds) {
                    soa_printf((GENinstance*) here, (GENmodel*) model, ckt,
                               "|Vds|=%g has exceeded Vds_max=%g\n",
                               vds, model->BSIM4v5vdsMax);
                    warns_vds++;
                }

            if (vbs > model->BSIM4v5vbsMax)
                if (warns_vbs < maxwarns_vbs) {
                    soa_printf((GENinstance*) here, (GENmodel*) model, ckt,
                               "|Vbs|=%g has exceeded Vbs_max=%g\n",
                               vbs, model->BSIM4v5vbsMax);
                    warns_vbs++;
                }

            if (vbd > model->BSIM4v5vbdMax)
                if (warns_vbd < maxwarns_vbd) {
                    soa_printf((GENinstance*) here, (GENmodel*) model, ckt,
                               "|Vbd|=%g has exceeded Vbd_max=%g\n",
                               vbd, model->BSIM4v5vbdMax);
                    warns_vbd++;
                }

        }
    }

    return OK;
}
