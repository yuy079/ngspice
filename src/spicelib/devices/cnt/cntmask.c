/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1987 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/const.h"
#include "ngspice/ifsim.h"
#include "ngspice/cktdefs.h"
#include "ngspice/devdefs.h"
#include "cntdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

/*ARGSUSED*/
int
CNTmAsk(CKTcircuit *ckt, GENmodel *inst, int which, IFvalue *value)
{
    CNTmodel *model = (CNTmodel *)inst;

    NG_IGNORE(ckt);

    switch(which) {
        case CNT_MOD_TNOM:
            value->rValue = model->CNTtnom-CONSTCtoK;
            return(OK);
        case CNT_MOD_VTO:
            value->rValue = model->CNTvt0;
            return(OK);
        case CNT_MOD_KP:
            value->rValue = model->CNTtransconductance;
            return(OK);
        case CNT_MOD_GAMMA:
            value->rValue = model->CNTgamma;
            return(OK);
        case CNT_MOD_PHI:
            value->rValue = model->CNTphi;
            return(OK);
        case CNT_MOD_LAMBDA:
            value->rValue = model->CNTlambda;
            return(OK);
        case CNT_MOD_RD:
            value->rValue = model->CNTdrainResistance;
            return(OK);
        case CNT_MOD_RS:
            value->rValue = model->CNTsourceResistance;
            return(OK);
        case CNT_MOD_CBD:
            value->rValue = model->CNTcapBD;
            return(OK);
        case CNT_MOD_CBS:
            value->rValue = model->CNTcapBS;
            return(OK);
        case CNT_MOD_IS:
            value->rValue = model->CNTjctSatCur;
            return(OK);
        case CNT_MOD_PB:
            value->rValue = model->CNTbulkJctPotential;
            return(OK);
        case CNT_MOD_CGSO:
            value->rValue = model->CNTgateSourceOverlapCapFactor;
            return(OK);
        case CNT_MOD_CGDO:
            value->rValue = model->CNTgateDrainOverlapCapFactor;
            return(OK);
        case CNT_MOD_CGBO:
            value->rValue = model->CNTgateBulkOverlapCapFactor;
            return(OK);
        case CNT_MOD_CJ:
            value->rValue = model->CNTbulkCapFactor;
            return(OK);
        case CNT_MOD_MJ:
            value->rValue = model->CNTbulkJctBotGradingCoeff;
            return(OK);
        case CNT_MOD_CJSW:
            value->rValue = model->CNTsideWallCapFactor;
            return(OK);
        case CNT_MOD_MJSW:
            value->rValue = model->CNTbulkJctSideGradingCoeff;
            return(OK);
        case CNT_MOD_JS:
            value->rValue = model->CNTjctSatCurDensity;
            return(OK);
        case CNT_MOD_TOX:
            value->rValue = model->CNToxideThickness;
            return(OK);
        case CNT_MOD_LD:
            value->rValue = model->CNTlatDiff;
            return(OK);
        case CNT_MOD_RSH:
            value->rValue = model->CNTsheetResistance;
            return(OK);
        case CNT_MOD_U0:
            value->rValue = model->CNTsurfaceMobility;
            return(OK);
        case CNT_MOD_FC:
            value->rValue = model->CNTfwdCapDepCoeff;
            return(OK);
        case CNT_MOD_NSUB:
            value->rValue = model->CNTsubstrateDoping;
            return(OK);
        case CNT_MOD_TPG:
            value->iValue = model->CNTgateType;
            return(OK);
        case CNT_MOD_NSS:
            value->rValue = model->CNTsurfaceStateDensity;
            return(OK);
        case CNT_MOD_TYPE:
	    if (model->CNTtype > 0)
		value->sValue = "ncnt";
	    else
		value->sValue = "pcnt";
            return(OK);
        default:
            return(E_BADPARM);
    }
    /* NOTREACHED */
}

