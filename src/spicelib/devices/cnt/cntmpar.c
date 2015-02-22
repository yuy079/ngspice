/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/const.h"
#include "ngspice/ifsim.h"
#include "cntdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int
CNTmParam(int param, IFvalue *value, GENmodel *inModel)
{
    CNTmodel *model = (CNTmodel *)inModel;
    switch(param) {
        case CNT_MOD_TNOM:
            model->CNTtnom = value->rValue + CONSTCtoK;
            model->CNTtnomGiven = TRUE;
            break;
        case CNT_MOD_VTO:
            model->CNTvt0 = value->rValue;
            model->CNTvt0Given = TRUE;
            break;
        case CNT_MOD_KP:
            model->CNTtransconductance = value->rValue;
            model->CNTtransconductanceGiven = TRUE;
            break;
        case CNT_MOD_GAMMA:
            model->CNTgamma = value->rValue;
            model->CNTgammaGiven = TRUE;
            break;
        case CNT_MOD_PHI:
            model->CNTphi = value->rValue;
            model->CNTphiGiven = TRUE;
            break;
        case CNT_MOD_LAMBDA:
            model->CNTlambda = value->rValue;
            model->CNTlambdaGiven = TRUE;
            break;
        case CNT_MOD_RD:
            model->CNTdrainResistance = value->rValue;
            model->CNTdrainResistanceGiven = TRUE;
            break;
        case CNT_MOD_RS:
            model->CNTsourceResistance = value->rValue;
            model->CNTsourceResistanceGiven = TRUE;
            break;
        case CNT_MOD_CBD:
            model->CNTcapBD = value->rValue;
            model->CNTcapBDGiven = TRUE;
            break;
        case CNT_MOD_CBS:
            model->CNTcapBS = value->rValue;
            model->CNTcapBSGiven = TRUE;
            break;
        case CNT_MOD_IS:
            model->CNTjctSatCur = value->rValue;
            model->CNTjctSatCurGiven = TRUE;
            break;
        case CNT_MOD_PB:
            model->CNTbulkJctPotential = value->rValue;
            model->CNTbulkJctPotentialGiven = TRUE;
            break;
        case CNT_MOD_CGSO:
            model->CNTgateSourceOverlapCapFactor = value->rValue;
            model->CNTgateSourceOverlapCapFactorGiven = TRUE;
            break;
        case CNT_MOD_CGDO:
            model->CNTgateDrainOverlapCapFactor = value->rValue;
            model->CNTgateDrainOverlapCapFactorGiven = TRUE;
            break;
        case CNT_MOD_CGBO:
            model->CNTgateBulkOverlapCapFactor = value->rValue;
            model->CNTgateBulkOverlapCapFactorGiven = TRUE;
            break;
        case CNT_MOD_CJ:
            model->CNTbulkCapFactor = value->rValue;
            model->CNTbulkCapFactorGiven = TRUE;
            break;
        case CNT_MOD_MJ:
            model->CNTbulkJctBotGradingCoeff = value->rValue;
            model->CNTbulkJctBotGradingCoeffGiven = TRUE;
            break;
        case CNT_MOD_CJSW:
            model->CNTsideWallCapFactor = value->rValue;
            model->CNTsideWallCapFactorGiven = TRUE;
            break;
        case CNT_MOD_MJSW:
            model->CNTbulkJctSideGradingCoeff = value->rValue;
            model->CNTbulkJctSideGradingCoeffGiven = TRUE;
            break;
        case CNT_MOD_JS:
            model->CNTjctSatCurDensity = value->rValue;
            model->CNTjctSatCurDensityGiven = TRUE;
            break;
        case CNT_MOD_TOX:
            model->CNToxideThickness = value->rValue;
            model->CNToxideThicknessGiven = TRUE;
            break;
        case CNT_MOD_LD:
            model->CNTlatDiff = value->rValue;
            model->CNTlatDiffGiven = TRUE;
            break;
        case CNT_MOD_RSH:
            model->CNTsheetResistance = value->rValue;
            model->CNTsheetResistanceGiven = TRUE;
            break;
        case CNT_MOD_U0:
            model->CNTsurfaceMobility = value->rValue;
            model->CNTsurfaceMobilityGiven = TRUE;
            break;
        case CNT_MOD_FC:
            model->CNTfwdCapDepCoeff = value->rValue;
            model->CNTfwdCapDepCoeffGiven = TRUE;
            break;
        case CNT_MOD_NSS:
            model->CNTsurfaceStateDensity = value->rValue;
            model->CNTsurfaceStateDensityGiven = TRUE;
            break;
        case CNT_MOD_NSUB:
            model->CNTsubstrateDoping = value->rValue;
            model->CNTsubstrateDopingGiven = TRUE;
            break;
        case CNT_MOD_TPG:
            model->CNTgateType = value->iValue;
            model->CNTgateTypeGiven = TRUE;
            break;
        case CNT_MOD_NCNT:
            if(value->iValue) {
                model->CNTtype = 1;
                model->CNTtypeGiven = TRUE;
            }
            break;
        case CNT_MOD_PCNT:
            if(value->iValue) {
                model->CNTtype = -1;
                model->CNTtypeGiven = TRUE;
            }
            break;
	case CNT_MOD_KF:
	    model->CNTfNcoef = value->rValue;
	    model->CNTfNcoefGiven = TRUE;
	    break;
	case CNT_MOD_AF:
	    model->CNTfNexp = value->rValue;
	    model->CNTfNexpGiven = TRUE;
	    break;
        default:
            return(E_BADPARM);
    }
    return(OK);
}
