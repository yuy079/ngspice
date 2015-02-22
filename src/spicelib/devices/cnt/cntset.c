/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/

    /* load the CNT device structure with those pointers needed later 
     * for fast matrix loading 
     */

#include "ngspice/ngspice.h"
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "cntdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int
CNTsetup(SMPmatrix *matrix, GENmodel *inModel, CKTcircuit *ckt,
          int *states)
{
    CNTmodel *model = (CNTmodel *)inModel;
    CNTinstance *here;
    int error;
    CKTnode *tmp;

    /*  loop through all the CNT device models */
    for( ; model != NULL; model = model->CNTnextModel ) {

        if(!model->CNTtypeGiven) {
            model->CNTtype = NCNT;
        }
        if(!model->CNTlatDiffGiven) {
            model->CNTlatDiff = 0;
        }
        if(!model->CNTjctSatCurDensityGiven) {
            model->CNTjctSatCurDensity = 0;
        }
        if(!model->CNTjctSatCurGiven) {
            model->CNTjctSatCur = 1e-14;
        }
        if(!model->CNTtransconductanceGiven) {
            model->CNTtransconductance = 2e-5;
        }
        if(!model->CNTgateSourceOverlapCapFactorGiven) {
            model->CNTgateSourceOverlapCapFactor = 0;
        }
        if(!model->CNTgateDrainOverlapCapFactorGiven) {
            model->CNTgateDrainOverlapCapFactor = 0;
        }
        if(!model->CNTgateBulkOverlapCapFactorGiven) {
            model->CNTgateBulkOverlapCapFactor = 0;
        }
        if(!model->CNTvt0Given) {
            model->CNTvt0 = 0;
        }
        if(!model->CNTbulkCapFactorGiven) {
            model->CNTbulkCapFactor = 0;
        }
        if(!model->CNTsideWallCapFactorGiven) {
            model->CNTsideWallCapFactor = 0;
        }
        if(!model->CNTbulkJctPotentialGiven) {
            model->CNTbulkJctPotential = .8;
        }
        if(!model->CNTbulkJctBotGradingCoeffGiven) {
            model->CNTbulkJctBotGradingCoeff = .5;
        }
        if(!model->CNTbulkJctSideGradingCoeffGiven) {
            model->CNTbulkJctSideGradingCoeff = .5;
        }
        if(!model->CNTfwdCapDepCoeffGiven) {
            model->CNTfwdCapDepCoeff = .5;
        }
        if(!model->CNTphiGiven) {
            model->CNTphi = .6;
        }
        if(!model->CNTlambdaGiven) {
            model->CNTlambda = 0;
        }
        if(!model->CNTgammaGiven) {
            model->CNTgamma = 0;
        }
	if(!model->CNTfNcoefGiven) {
	    model->CNTfNcoef = 0;
	}
	if(!model->CNTfNexpGiven) {
	    model->CNTfNexp = 1;
	}

        /* loop through all the instances of the model */
        for (here = model->CNTinstances; here != NULL ;
                here=here->CNTnextInstance) {

	    if (/*here->CNTowner == ARCHme*/1) {
		/* allocate a chunk of the state vector */
		here->CNTstates = *states;
		*states += CNTnumStates;
		if(ckt->CKTsenInfo && (ckt->CKTsenInfo->SENmode & TRANSEN) ){
		    *states += 10 * (ckt->CKTsenInfo->SENparms);
		}
	    }

            if(!here->CNTdrainPerimiterGiven) {
                here->CNTdrainPerimiter = 0;
            }
            if(!here->CNTicVBSGiven) {
                here->CNTicVBS = 0;
            }
            if(!here->CNTicVDSGiven) {
                here->CNTicVDS = 0;
            }
            if(!here->CNTicVGSGiven) {
                here->CNTicVGS = 0;
            }
            if(!here->CNTsourcePerimiterGiven) {
                here->CNTsourcePerimiter = 0;
            }
            if(!here->CNTvdsatGiven) {
                here->CNTvdsat = 0;
            }
            if(!here->CNTvonGiven) {
                here->CNTvon = 0;
            }
	    if(!here->CNTdrainSquaresGiven) {
		here->CNTdrainSquares=1;
	    }
	    if(!here->CNTsourceSquaresGiven) {
		here->CNTsourceSquares=1;
	    }

            if ((model->CNTdrainResistance != 0
		    || (model->CNTsheetResistance != 0
                    && here->CNTdrainSquares != 0) )
		    && here->CNTdNodePrime == 0) {
                error = CKTmkVolt(ckt,&tmp,here->CNTname,"drain");
                if(error) return(error);
                here->CNTdNodePrime = tmp->number;
                
                if (ckt->CKTcopyNodesets) {
		    CKTnode *tmpNode;
		    IFuid tmpName;

                  if (CKTinst2Node(ckt,here,1,&tmpNode,&tmpName)==OK) {
                     if (tmpNode->nsGiven) {
                       tmp->nodeset=tmpNode->nodeset; 
                       tmp->nsGiven=tmpNode->nsGiven; 
                     }
                  }
                }
                
            } else {
                here->CNTdNodePrime = here->CNTdNode;
            }

            if((model->CNTsourceResistance != 0 ||
                    (model->CNTsheetResistance != 0 &&
                     here->CNTsourceSquares != 0) ) &&
                    here->CNTsNodePrime==0) {
                error = CKTmkVolt(ckt,&tmp,here->CNTname,"source");
                if(error) return(error);
                here->CNTsNodePrime = tmp->number;
                
                if (ckt->CKTcopyNodesets) {
		    CKTnode *tmpNode;
		    IFuid tmpName;

                  if (CKTinst2Node(ckt,here,3,&tmpNode,&tmpName)==OK) {
                     if (tmpNode->nsGiven) {
                       tmp->nodeset=tmpNode->nodeset; 
                       tmp->nsGiven=tmpNode->nsGiven; 
                     }
                  }
                }
                
            } else {
                here->CNTsNodePrime = here->CNTsNode;
            }

/* macro to make elements with built in test for out of memory */
#define TSTALLOC(ptr,first,second) \
if((here->ptr = SMPmakeElt(matrix,here->first,here->second))==(double *)NULL){\
    return(E_NOMEM);\
}
            TSTALLOC(CNTDdPtr,CNTdNode,CNTdNode)
            TSTALLOC(CNTGgPtr,CNTgNode,CNTgNode)
            TSTALLOC(CNTSsPtr,CNTsNode,CNTsNode)
            TSTALLOC(CNTBbPtr,CNTbNode,CNTbNode)
            TSTALLOC(CNTDPdpPtr,CNTdNodePrime,CNTdNodePrime)
            TSTALLOC(CNTSPspPtr,CNTsNodePrime,CNTsNodePrime)
            TSTALLOC(CNTDdpPtr,CNTdNode,CNTdNodePrime)
            TSTALLOC(CNTGbPtr,CNTgNode,CNTbNode)
            TSTALLOC(CNTGdpPtr,CNTgNode,CNTdNodePrime)
            TSTALLOC(CNTGspPtr,CNTgNode,CNTsNodePrime)
            TSTALLOC(CNTSspPtr,CNTsNode,CNTsNodePrime)
            TSTALLOC(CNTBdpPtr,CNTbNode,CNTdNodePrime)
            TSTALLOC(CNTBspPtr,CNTbNode,CNTsNodePrime)
            TSTALLOC(CNTDPspPtr,CNTdNodePrime,CNTsNodePrime)
            TSTALLOC(CNTDPdPtr,CNTdNodePrime,CNTdNode)
            TSTALLOC(CNTBgPtr,CNTbNode,CNTgNode)
            TSTALLOC(CNTDPgPtr,CNTdNodePrime,CNTgNode)
            TSTALLOC(CNTSPgPtr,CNTsNodePrime,CNTgNode)
            TSTALLOC(CNTSPsPtr,CNTsNodePrime,CNTsNode)
            TSTALLOC(CNTDPbPtr,CNTdNodePrime,CNTbNode)
            TSTALLOC(CNTSPbPtr,CNTsNodePrime,CNTbNode)
            TSTALLOC(CNTSPdpPtr,CNTsNodePrime,CNTdNodePrime)

        }
    }
    return(OK);
}

int
CNTunsetup(GENmodel *inModel, CKTcircuit *ckt)
{
    CNTmodel *model;
    CNTinstance *here;

    for (model = (CNTmodel *)inModel; model != NULL;
	    model = model->CNTnextModel)
    {
        for (here = model->CNTinstances; here != NULL;
                here=here->CNTnextInstance)
	{
	    if (here->CNTdNodePrime
		    && here->CNTdNodePrime != here->CNTdNode)
	    {
		CKTdltNNum(ckt, here->CNTdNodePrime);
		here->CNTdNodePrime= 0;
	    }
	    if (here->CNTsNodePrime
		    && here->CNTsNodePrime != here->CNTsNode)
	    {
		CKTdltNNum(ckt, here->CNTsNodePrime);
		here->CNTsNodePrime= 0;
	    }
	}
    }
    return OK;
}
