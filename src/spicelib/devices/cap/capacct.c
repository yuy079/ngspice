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

extern FILE *slogp;  /* soa log file ('--soa-log file' command line option) */

/* make SOA checks after NR has finished */

int
CAPaccept(CKTcircuit *ckt, GENmodel *inModel)
{
    CAPmodel *model = (CAPmodel *) inModel;
    CAPinstance *here;

    int maxwarns_bv = 0;
    static int warns_bv = 0;


    if(!(ckt->CKTmode & (MODETRAN | MODETRANOP)))
        return OK;

    if (!ckt->CKTsoaCheck)
        return OK;

    for(; model; model = model->CAPnextModel)
        for (here = model->CAPinstances; here; here = here->CAPnextInstance) {

            double vc;  /* current capacitor voltage */

            vc = ckt->CKTrhsOld [here->CAPposNode] -
                 ckt->CKTrhsOld [here->CAPnegNode];

            if (vc > model->CAPbv_max)
                if (warns_bv < maxwarns_bv) {
                    printf("Instance: %s Model: %s Time: %g |Vc|=%g has exceeded Bv_max=%g\n",
                           here->CAPname, model->CAPmodName, ckt->CKTtime, vc, model->CAPbv_max);
                    if (slogp)
                        fprintf(slogp, "Instance: %s Model: %s Time: %g |Vc|=%g has exceeded Bv_max=%g\n",
                                here->CAPname, model->CAPmodName, ckt->CKTtime, vc, model->CAPbv_max);
                    warns_bv++;
                }

        }

    return OK;
}
