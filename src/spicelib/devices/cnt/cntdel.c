/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/
/*
 */

#include "ngspice.h"
#include "cntdefs.h"
#include "sperror.h"
#include "suffix.h"


int
CNTdelete(GENmodel *inModel, IFuid name, GENinstance **inst)
{
    CNTmodel *model = (CNTmodel *)inModel;
    CNTinstance **fast = (CNTinstance **)inst;
    CNTinstance **prev = NULL;
    CNTinstance *here;

    for( ; model ; model = model->CNTnextModel) {
        prev = &(model->CNTinstances);
        for(here = *prev; here ; here = *prev) {
            if(here->CNTname == name || (fast && here==*fast) ) {
                *prev= here->CNTnextInstance;
                FREE(here);
                return(OK);
            }
            prev = &(here->CNTnextInstance);
        }
    }
    return(E_NODEV);
}
