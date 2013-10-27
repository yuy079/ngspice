/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bjtdefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/cpdefs.h"

void
soa_printf(GENinstance *, GENmodel *, CKTcircuit *, const char *, ...);

/* make SOA checks after NR has finished */

int
BJTsoaCheck(CKTcircuit *ckt, GENmodel *inModel)
{
    BJTmodel *model = (BJTmodel *) inModel;
    BJTinstance *here;
    double vbe, vbc, vce;    /* actual bjt voltages */
    int maxwarns_vbe = 0, maxwarns_vbc = 0, maxwarns_vce = 0;
    static int warns_vbe = 0, warns_vbc = 0, warns_vce = 0;

    if(!(ckt->CKTmode & (MODEDC | MODEDCOP | MODEDCTRANCURVE | MODETRAN | MODETRANOP)))
        return OK;

    for (; model; model = model->BJTnextModel) {

        maxwarns_vbe = maxwarns_vbc = maxwarns_vce = ckt->CKTsoaMaxWarns;

        for (here = model->BJTinstances; here; here=here->BJTnextInstance) {

            vbe = fabs(ckt->CKTrhsOld [here->BJTbasePrimeNode] -
                       ckt->CKTrhsOld [here->BJTemitPrimeNode]);
            vbc = fabs(ckt->CKTrhsOld [here->BJTbasePrimeNode] -
                       ckt->CKTrhsOld [here->BJTcolPrimeNode]);
            vce = fabs(vbe - vbc);

            if (vbe > model->BJTvbeMax)
                if (warns_vbe < maxwarns_vbe) {
                    soa_printf((GENinstance*) here, (GENmodel*) model, ckt,
                               "|Vbe|=%g has exceeded Vbe_max=%g\n",
                               vbe, model->BJTvbeMax);
                    warns_vbe++;
                }

            if (vbc > model->BJTvbcMax)
                if (warns_vbc < maxwarns_vbc) {
                    soa_printf((GENinstance*) here, (GENmodel*) model, ckt,
                               "|Vbc|=%g has exceeded Vbc_max=%g\n",
                               vbc, model->BJTvbcMax);
                    warns_vbc++;
                }

            if (vce > model->BJTvceMax)
                if (warns_vce < maxwarns_vce) {
                    soa_printf((GENinstance*) here, (GENmodel*) model, ckt,
                               "|Vce|=%g has exceeded Vce_max=%g\n",
                               vce, model->BJTvceMax);
                    warns_vce++;
                }

        }
    }

    return OK;
}
