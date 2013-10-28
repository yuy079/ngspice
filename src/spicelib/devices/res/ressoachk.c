/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "resdefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/cpdefs.h"


/* make SOA checks after NR has finished */

int
RESsoaCheck(CKTcircuit *ckt, GENmodel *inModel)
{
    RESmodel *model = (RESmodel *) inModel;
    RESinstance *here;
    double vr;  /* current resistor voltage */
    int maxwarns_bv = 0;
    static int warns_bv = 0;

    for (; model; model = model->RESnextModel) {

        maxwarns_bv = ckt->CKTsoaMaxWarns;

        for (here = model->RESinstances; here; here = here->RESnextInstance) {

            vr = fabs(ckt->CKTrhsOld [here->RESposNode] -
                      ckt->CKTrhsOld [here->RESnegNode]);

            if (vr > here->RESbv_max)
                if (warns_bv < maxwarns_bv) {
                    soa_printf(ckt, (GENinstance*) here,
                               "|Vr|=%g has exceeded Bv_max=%g\n",
                               vr, here->RESbv_max);
                    warns_bv++;
                }

        }

    }

    return OK;
}
