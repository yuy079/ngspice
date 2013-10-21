/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "diodefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/cpdefs.h"

extern FILE *slogp;  /* soa log file ('--soa-log file' command line option) */

/* make SOA checks after NR has finished */

int
DIOaccept(CKTcircuit *ckt, GENmodel *inModel)
{
    DIOmodel *model = (DIOmodel *) inModel;
    DIOinstance *here;

    int maxwarns_fv = 0, maxwarns_bv = 0;
    static int warns_fv = 0, warns_bv = 0;


    if(!(ckt->CKTmode & (MODETRAN | MODETRANOP)))
        return OK;

    if (!ckt->CKTsoaCheck)
        return OK;

    for(; model; model = model->DIOnextModel)
        for (here = model->DIOinstances; here; here = here->DIOnextInstance) {

            double vd;  /* current diode voltage */

            vd = ckt->CKTrhsOld [here->DIOposPrimeNode] -
                 ckt->CKTrhsOld [here->DIOnegNode];

            if (vd > model->DIOfv_max)
                if (warns_fv < maxwarns_fv) {
                    printf("Instance: %s Model: %s Time: %g Vf=%g has exceeded Fv_max=%g\n",
                           here->DIOname, model->DIOmodName, ckt->CKTtime, vd, model->DIOfv_max);
                    if (slogp)
                        fprintf(slogp, "Instance: %s Model: %s Time: %g Vf=%g has exceeded Fv_max=%g\n",
                                here->DIOname, model->DIOmodName, ckt->CKTtime, vd, model->DIOfv_max);
                    warns_fv++;
                }

            if (-vd > model->DIObv_max)
                if (warns_bv < maxwarns_bv) {
                    printf("Instance: %s Model: %s Time: %g |Vj|=%g has exceeded Bv_max=%g\n",
                           here->DIOname, model->DIOmodName, ckt->CKTtime, -vd, model->DIObv_max);
                    if (slogp)
                        fprintf(slogp, "Instance: %s Model: %s Time: %g |Vj|=%g has exceeded Bv_max=%g\n",
                                here->DIOname, model->DIOmodName, ckt->CKTtime, -vd, model->DIObv_max);
                    warns_bv++;
                }

        }

    return OK;
}
