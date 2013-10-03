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


/* make SOA checks after NR has finished */

int
DIOaccept(CKTcircuit *ckt, GENmodel *inModel)
{
    DIOmodel *model = (DIOmodel *) inModel;
    DIOinstance *here;

    int warn;   /* whether SOA check should be performed */

    if (!cp_getvar("warn", CP_NUM, &warn))
        warn = 0;

    if (warn <= 0)
        return OK;

    if(!(ckt->CKTmode & (MODETRAN | MODETRANOP)))
        return OK;

    for(; model; model = model->DIOnextModel)
        for (here = model->DIOinstances; here; here = here->DIOnextInstance) {

            double vd;  /* current diode voltage */

            vd = ckt->CKTrhsOld [here->DIOposPrimeNode] -
                ckt->CKTrhsOld [here->DIOnegNode];

            if (vd > model->DIOfv_max)
                printf("Instance: %s Model: %s Time: %g Vf=%g has exceeded Fv_max=%g\n",
                       here->DIOname, model->DIOmodName, ckt->CKTtime, vd, model->DIOfv_max);

            if (-vd > model->DIObv_max)
                printf("Instance: %s Model: %s Time: %g |Vj|=%g has exceeded Bv_max=%g\n",
                       here->DIOname, model->DIOmodName, ckt->CKTtime, -vd, model->DIObv_max);

        }

    return OK;
}
