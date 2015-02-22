/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000 AlansFixes

This function is obsolete (was used by an old sensitivity analysis)
**********/

/* Pretty print the sensitivity info for all 
 * the CNT devices in the circuit.
 */

#include "ngspice.h"
#include "smpdefs.h"
#include "cktdefs.h"
#include "cntdefs.h"
#include "sperror.h"
#include "suffix.h"

void
CNTsPrint(GENmodel *inModel, CKTcircuit *ckt)
/* Pretty print the sensitivity info for all the CNT 
         * devices  in the circuit.
         */
{
    CNTmodel *model = (CNTmodel *)inModel;
    CNTinstance *here;

    printf("LEVEL 1 CNTFETS-----------------\n");
    /*  loop through all the CNT models */
    for( ; model != NULL; model = model->CNTnextModel ) {

        printf("Model name:%s\n",model->CNTmodName);

        /* loop through all the instances of the model */
        for (here = model->CNTinstances; here != NULL ;
                here=here->CNTnextInstance) {
	    if (here->CNTowner != ARCHme) continue;

            printf("    Instance name:%s\n",here->CNTname);
            printf("      Drain, Gate , Source nodes: %s, %s ,%s\n",
            CKTnodName(ckt,here->CNTdNode),CKTnodName(ckt,here->CNTgNode),
            CKTnodName(ckt,here->CNTsNode));
            
            printf("  Multiplier: %g ",here->CNTm);
            printf(here->CNTmGiven ? "(specified)\n" : "(default)\n");
            
            printf("      Length: %g ",here->CNTl);
            printf(here->CNTlGiven ? "(specified)\n" : "(default)\n");
            printf("      Width: %g ",here->CNTw);
            printf(here->CNTwGiven ? "(specified)\n" : "(default)\n");
            if(here->CNTsens_l == 1){
                printf("    CNTsenParmNo:l = %d ",here->CNTsenParmNo);
            }
            else{ 
                printf("    CNTsenParmNo:l = 0 ");
            }
            if(here->CNTsens_w == 1){
                printf("    w = %d \n",here->CNTsenParmNo + here->CNTsens_l);
            }
            else{ 
                printf("    w = 0 \n");
            }


        }
    }
}

