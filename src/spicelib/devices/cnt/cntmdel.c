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
CNTmDelete(GENmodel **inModel, IFuid modname, GENmodel *kill)
{
    CNTmodel **model = (CNTmodel **)inModel;
    CNTmodel *modfast = (CNTmodel *)kill;
    CNTinstance *here;
    CNTinstance *prev = NULL;
    CNTmodel **oldmod;
    oldmod = model;
    for( ; *model ; model = &((*model)->CNTnextModel)) {
        if( (*model)->CNTmodName == modname || 
                (modfast && *model == modfast) ) goto delgot;
        oldmod = model;
    }
    return(E_NOMOD);

delgot:
    *oldmod = (*model)->CNTnextModel; /* cut deleted device out of list */
    for(here = (*model)->CNTinstances ; here ; here = here->CNTnextInstance) {
        if(prev) FREE(prev);
        prev = here;
    }
    if(prev) FREE(prev);
    FREE(*model);
    return(OK);

}
