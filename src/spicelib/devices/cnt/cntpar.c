/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/
/*
 */

#include "ngspice/ngspice.h"
#include "ngspice/const.h"
#include "ngspice/ifsim.h"
#include "cntdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"


/* ARGSUSED */
int
CNTparam(int param, IFvalue *value, GENinstance *inst, IFvalue *select)
{
    CNTinstance *here = (CNTinstance *)inst;

    NG_IGNORE(select);

    switch(param) {
        case CNT_TEMP:
            here->CNTtemp = value->rValue+CONSTCtoK;
            here->CNTtempGiven = TRUE;
            break;
        case CNT_DTEMP:
            here->CNTdtemp = value->rValue;
            here->CNTdtempGiven = TRUE;
            break;
        case CNT_M:
            here->CNTm = value->rValue;
            here->CNTmGiven = TRUE;
            break;
        case CNT_W:
            here->CNTw = value->rValue;
            here->CNTwGiven = TRUE;
            break;
        case CNT_L:
            here->CNTl = value->rValue;
            here->CNTlGiven = TRUE;
            break;
        case CNT_DIA:
            here->CNTdia = value->rValue;
            here->CNTdiaGiven = TRUE;
            break;
        case CNT_EF:
            here->CNTef = value->rValue;
            here->CNTefGiven = TRUE;
            break;
        case CNT_AS:
            here->CNTsourceArea = value->rValue;
            here->CNTsourceAreaGiven = TRUE;
            break;
        case CNT_AD:
            here->CNTdrainArea = value->rValue;
            here->CNTdrainAreaGiven = TRUE;
            break;
        case CNT_PS:
            here->CNTsourcePerimiter = value->rValue;
            here->CNTsourcePerimiterGiven = TRUE;
            break;
        case CNT_PD:
            here->CNTdrainPerimiter = value->rValue;
            here->CNTdrainPerimiterGiven = TRUE;
            break;
        case CNT_NRS:
            here->CNTsourceSquares = value->rValue;
            here->CNTsourceSquaresGiven = TRUE;
            break;
        case CNT_NRD:
            here->CNTdrainSquares = value->rValue;
            here->CNTdrainSquaresGiven = TRUE;
            break;
        case CNT_OFF:
            here->CNToff = (value->iValue != 0);
            break;
        case CNT_IC_VBS:
            here->CNTicVBS = value->rValue;
            here->CNTicVBSGiven = TRUE;
            break;
        case CNT_IC_VDS:
            here->CNTicVDS = value->rValue;
            here->CNTicVDSGiven = TRUE;
            break;
        case CNT_IC_VGS:
            here->CNTicVGS = value->rValue;
            here->CNTicVGSGiven = TRUE;
            break;
        case CNT_IC:
            switch(value->v.numValue){
                case 3:
                    here->CNTicVBS = *(value->v.vec.rVec+2);
                    here->CNTicVBSGiven = TRUE;
                case 2:
                    here->CNTicVGS = *(value->v.vec.rVec+1);
                    here->CNTicVGSGiven = TRUE;
                case 1:
                    here->CNTicVDS = *(value->v.vec.rVec);
                    here->CNTicVDSGiven = TRUE;
                    break;
                default:
                    return(E_BADPARM);
            }
            break;
        case CNT_L_SENS:
            if(value->iValue) {
                here->CNTsenParmNo = 1;
                here->CNTsens_l = 1;
            }
            break;
        case CNT_W_SENS:
            if(value->iValue) {
                here->CNTsenParmNo = 1;
                here->CNTsens_w = 1;
            }
            break;
        default:
            return(E_BADPARM);
    }
    return(OK);
}
