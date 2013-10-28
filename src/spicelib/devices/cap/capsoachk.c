/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "capdefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/cpdefs.h"


/* make SOA checks after NR has finished */

int
CAPsoaCheck(CKTcircuit *ckt, GENmodel *inModel)
{
    CAPmodel *model = (CAPmodel *) inModel;
    CAPinstance *here;
    double vc;  /* current capacitor voltage */
    int maxwarns_bv = 0;
    static int warns_bv = 0;

    if (!(ckt->CKTmode & (MODEDC | MODEDCOP | MODEDCTRANCURVE | MODETRAN | MODETRANOP)))
        return OK;

    for (; model; model = model->CAPnextModel) {

        maxwarns_bv = ckt->CKTsoaMaxWarns;

        for (here = model->CAPinstances; here; here = here->CAPnextInstance) {

            vc = ckt->CKTrhsOld [here->CAPposNode] -
                 ckt->CKTrhsOld [here->CAPnegNode];

            if (vc > here->CAPbv_max)
                if (warns_bv < maxwarns_bv) {
                    soa_printf(ckt, (GENinstance*) here,
                               "|Vc|=%g has exceeded Bv_max=%g\n",
                               vc, here->CAPbv_max);
                    warns_bv++;
                }

        }

    }

    return OK;
}
