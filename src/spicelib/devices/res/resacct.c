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

extern FILE *slogp;  /* soa log file ('--soa-log file' command line option) */

/* make SOA checks after NR has finished */

int
RESaccept(CKTcircuit *ckt, GENmodel *inModel)
{
    RESmodel *model = (RESmodel *) inModel;
    RESinstance *here;

    int maxwarns_bv = 0;
    static int warns_bv = 0;


    if(!(ckt->CKTmode & (MODETRAN | MODETRANOP)))
        return OK;

    if (!ckt->CKTsoaCheck)
        return OK;

    for(; model; model = model->RESnextModel)

        maxwarns_bv = ckt->CKTsoaMaxWarns;

        for (here = model->RESinstances; here; here = here->RESnextInstance) {

            double vr;  /* current resistor voltage */

            vr = ckt->CKTrhsOld [here->RESposNode] -
                 ckt->CKTrhsOld [here->RESnegNode];

            if (vr > here->RESbv_max)
                if (warns_bv < maxwarns_bv) {
                    printf("Instance: %s Model: %s Time: %g |Vr|=%g has exceeded Bv_max=%g\n",
                           here->RESname, model->RESmodName, ckt->CKTtime, vr, here->RESbv_max);
                    if (slogp)
                        fprintf(slogp, "Instance: %s Model: %s Time: %g |Vr|=%g has exceeded Bv_max=%g\n",
                                here->RESname, model->RESmodName, ckt->CKTtime, vr, here->RESbv_max);
                    warns_bv++;
                }

        }

    return OK;
}
