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

static void
soa_printf(GENinstance *instance, GENmodel *model, CKTcircuit *ckt, const char *fmt, ...)
{
    FILE *fp;

    va_list ap;
    va_start(ap, fmt);

    for (fp = stdout; fp; fp = slogp) {

        if(ckt->CKTmode & MODETRAN)
            fprintf(fp, "Instance: %s Model: %s Time: %g ",
                    instance->GENname, model->GENmodName, ckt->CKTtime);
        else
            fprintf(fp, "Instance: %s Model: %s ",
                    instance->GENname, model->GENmodName);

        vfprintf(fp, fmt, ap);

        if (fp == slogp)
            break;
    }

    va_end(ap);
}


/* make SOA checks after NR has finished */

int
DIOaccept(CKTcircuit *ckt, GENmodel *inModel)
{
    DIOmodel *model = (DIOmodel *) inModel;
    DIOinstance *here;

    int maxwarns_fv = 0, maxwarns_bv = 0;
    static int warns_fv = 0, warns_bv = 0;

    if (!ckt->CKTsoaCheck)
        return OK;

    for(; model; model = model->DIOnextModel) {

        maxwarns_fv = maxwarns_bv = ckt->CKTsoaMaxWarns;

        for (here = model->DIOinstances; here; here = here->DIOnextInstance) {

            double vd;  /* current diode voltage */

            vd = ckt->CKTrhsOld [here->DIOposPrimeNode] -
                 ckt->CKTrhsOld [here->DIOnegNode];

            if (vd > model->DIOfv_max)
                if (warns_fv < maxwarns_fv) {
                    soa_printf((GENinstance*) here, (GENmodel*) model, ckt,
                               "Vf=%g has exceeded Fv_max=%g\n",
                               vd, model->DIOfv_max);
                    warns_fv++;
                }

            if (fabs(vd) > model->DIObv_max)
                if (warns_bv < maxwarns_bv) {
                    soa_printf((GENinstance*) here, (GENmodel*) model, ckt,
                               "|Vj|=%g has exceeded Bv_max=%g\n",
                               vd, model->DIObv_max);
                    warns_bv++;
                }

        }

    }

    return OK;
}
