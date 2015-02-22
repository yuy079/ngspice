/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/
/*
 */

#include "ngspice/ngspice.h"
#include "cntdefs.h"
#include "ngspice/suffix.h"


void
CNTdestroy(GENmodel **inModel)
{
    CNTmodel **model = (CNTmodel **)inModel;
    CNTinstance *here;
    CNTinstance *prev = NULL;
    CNTmodel *mod = *model;
    CNTmodel *oldmod = NULL;

    for( ; mod ; mod = mod->CNTnextModel) {
        if(oldmod) FREE(oldmod);
        oldmod = mod;
        prev = (CNTinstance *)NULL;
        for(here = mod->CNTinstances ; here ; here = here->CNTnextInstance) {
            if(prev){
          if(prev->CNTsens) FREE(prev->CNTsens);
          FREE(prev);
        }
            prev = here;
        }
        if(prev) FREE(prev);
    }
    if(oldmod) FREE(oldmod);
    *model = NULL;
}
