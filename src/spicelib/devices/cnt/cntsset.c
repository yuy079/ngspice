/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles

This function is obsolete (was used by an old sensitivity analysis)
**********/

#include "ngspice.h"
#include "smpdefs.h"
#include "cktdefs.h"
#include "cntdefs.h"
#include "sperror.h"
#include "suffix.h"

int
CNTsSetup(SENstruct *info, GENmodel *inModel)
/* loop through all the devices and 
         * allocate parameter #s to design parameters 
         */
{
    CNTmodel *model = (CNTmodel *)inModel;
    CNTinstance *here;

    /*  loop through all the models */
    for( ; model != NULL; model = model->CNTnextModel ) {

        /* loop through all the instances of the model */
        for (here = model->CNTinstances; here != NULL ;
                here=here->CNTnextInstance) {
	    if (here->CNTowner != ARCHme) continue;

            if(here->CNTsenParmNo){
                if((here->CNTsens_l)&&(here->CNTsens_w)){
                    here->CNTsenParmNo = ++(info->SENparms);
                    ++(info->SENparms);/* CNT has two design parameters */
                }
                else{
                    here->CNTsenParmNo = ++(info->SENparms);
                }
            }
            if((here->CNTsens = (double *)MALLOC(70*sizeof(double))) == NULL) {
                return(E_NOMEM);
            }
            here->CNTsenPertFlag = OFF;

        }
    }
    return(OK);
}


