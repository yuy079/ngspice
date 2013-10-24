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

        if(ckt->CKTmode & (MODEDC | MODEDCOP))
            fprintf(fp, "Instance: %s Model: %s ",
                    instance->GENname, model->GENmodName);
        else
            fprintf(fp, "Instance: %s Model: %s Time: %g ",
                    instance->GENname, model->GENmodName, ckt->CKTtime);

        vfprintf(fp, fmt, ap);

        if (fp == slogp)
            break;
    }

    va_end(ap);
}


#define soa_check(value, limit, warn_counter, message)                  \
    do {                                                                \
        if (fabs(value) > limit  &&  warn_counter < max##warn_counter) { \
            soa_printf((GENinstance*) here, (GENmodel*) model, ckt, message, value, limit); \
            warn_counter++;                                             \
        }                                                               \
    } while (0)


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

            // alternativ
            soa_check(vd, model->DIObv_max, warns_bv, "|Vj|=%g has exceeded Bv_max=%g\n");
        }

    return OK;
}
