/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1987 Thomas L. Quarles
Modified: 2000 AlansFixes
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
CNTask(CKTcircuit *ckt, GENinstance *inst, int which, IFvalue *value,
        IFvalue *select)
{
    CNTinstance *here = (CNTinstance*)inst;
    double vr;
    double vi;
    double sr;
    double si;
    double vm;
    static char *msg = "Current and power not available for ac analysis";
    switch(which) {
        case CNT_TEMP:
            value->rValue = here->CNTtemp - CONSTCtoK;
            return(OK);
        case CNT_DTEMP:
            value->rValue = here->CNTdtemp;
            return(OK);
        case CNT_CGS:
            value->rValue = 2*  *(ckt->CKTstate0 + here->CNTcapgs);
            return(OK);
        case CNT_CGD:
            value->rValue = 2* *(ckt->CKTstate0 + here->CNTcapgd);
            return(OK);   
        case CNT_M:
            value->rValue = here->CNTm;
            return(OK);    
        case CNT_L:
            value->rValue = here->CNTl;
                return(OK);
        case CNT_W:
            value->rValue = here->CNTw;
                return(OK);
        case CNT_AS:
            value->rValue = here->CNTsourceArea;
                return(OK);
        case CNT_AD:
            value->rValue = here->CNTdrainArea;
                return(OK);
        case CNT_PS:
            value->rValue = here->CNTsourcePerimiter;
                return(OK);
        case CNT_PD:
            value->rValue = here->CNTdrainPerimiter;
                return(OK);
        case CNT_NRS:
            value->rValue = here->CNTsourceSquares;
                return(OK);
        case CNT_NRD:
            value->rValue = here->CNTdrainSquares;
                return(OK);
        case CNT_OFF:
            value->rValue = here->CNToff;
                return(OK);
        case CNT_IC_VBS:
            value->rValue = here->CNTicVBS;
                return(OK);
        case CNT_IC_VDS:
            value->rValue = here->CNTicVDS;
                return(OK);
        case CNT_IC_VGS:
            value->rValue = here->CNTicVGS;
                return(OK);
        case CNT_DNODE:
            value->iValue = here->CNTdNode;
            return(OK);
        case CNT_GNODE:
            value->iValue = here->CNTgNode;
            return(OK);
        case CNT_SNODE:
            value->iValue = here->CNTsNode;
            return(OK);
        case CNT_BNODE:
            value->iValue = here->CNTbNode;
            return(OK);
        case CNT_DNODEPRIME:
            value->iValue = here->CNTdNodePrime;
            return(OK);
        case CNT_SNODEPRIME:
            value->iValue = here->CNTsNodePrime;
            return(OK);
        case CNT_SOURCECONDUCT:
            value->rValue = here->CNTsourceConductance;
            return(OK);
        case CNT_SOURCERESIST:
	    if (here->CNTsNodePrime != here->CNTsNode)
		value->rValue = 1.0 / here->CNTsourceConductance;
	    else
		value->rValue = 0.0;
            return(OK);
        case CNT_DRAINCONDUCT:
            value->rValue = here->CNTdrainConductance;
            return(OK);
        case CNT_DRAINRESIST:
	    if (here->CNTdNodePrime != here->CNTdNode)
		value->rValue = 1.0 / here->CNTdrainConductance;
	    else
		value->rValue = 0.0;
            return(OK);
        case CNT_VON:
            value->rValue = here->CNTvon;
            return(OK);
        case CNT_VDSAT:
            value->rValue = here->CNTvdsat;
            return(OK);
        case CNT_SOURCEVCRIT:
            value->rValue = here->CNTsourceVcrit;
            return(OK);
        case CNT_DRAINVCRIT:
            value->rValue = here->CNTdrainVcrit;
            return(OK);
        case CNT_CD:
            value->rValue = here->CNTcd;
            return(OK);
        case CNT_CBS:
            value->rValue = here->CNTcbs;
            return(OK);
        case CNT_CBD:
            value->rValue = here->CNTcbd;
            return(OK);
        case CNT_GMBS:
            value->rValue = here->CNTgmbs;
            return(OK);
        case CNT_GM:
            value->rValue = here->CNTgm;
            return(OK);
        case CNT_GDS:
            value->rValue = here->CNTgds;
            return(OK);
        case CNT_GBD:
            value->rValue = here->CNTgbd;
            return(OK);
        case CNT_GBS:
            value->rValue = here->CNTgbs;
            return(OK);
        case CNT_CAPBD:
            value->rValue = here->CNTcapbd;
            return(OK);
        case CNT_CAPBS:
            value->rValue = here->CNTcapbs;
            return(OK);
        case CNT_CAPZEROBIASBD:
            value->rValue = here->CNTCbd;
            return(OK);
        case CNT_CAPZEROBIASBDSW:
            value->rValue = here->CNTCbdsw;
            return(OK);
        case CNT_CAPZEROBIASBS:
            value->rValue = here->CNTCbs;
            return(OK);
        case CNT_CAPZEROBIASBSSW:
            value->rValue = here->CNTCbssw;
            return(OK);
        case CNT_VBD:
            value->rValue = *(ckt->CKTstate0 + here->CNTvbd);
            return(OK);
        case CNT_VBS:
            value->rValue = *(ckt->CKTstate0 + here->CNTvbs);
            return(OK);
        case CNT_VGS:
            value->rValue = *(ckt->CKTstate0 + here->CNTvgs);
            return(OK);
        case CNT_VDS:
            value->rValue = *(ckt->CKTstate0 + here->CNTvds);
            return(OK);
        case CNT_CAPGS:
            value->rValue = 2* *(ckt->CKTstate0 + here->CNTcapgs);
            /* add overlap capacitance */
            value->rValue += (here->sCNTmodPtr->CNTgateSourceOverlapCapFactor)
                             * here->CNTm
                             * (here->CNTw);
            return(OK);
        case CNT_QGS:
            value->rValue = *(ckt->CKTstate0 + here->CNTqgs);
            return(OK);
        case CNT_CQGS:
            value->rValue = *(ckt->CKTstate0 + here->CNTcqgs);
            return(OK);
        case CNT_CAPGD:
            value->rValue = 2* *(ckt->CKTstate0 + here->CNTcapgd);
            /* add overlap capacitance */
            value->rValue += (here->sCNTmodPtr->CNTgateSourceOverlapCapFactor)
                             * here->CNTm
                             * (here->CNTw);
            return(OK);
        case CNT_QGD:
            value->rValue = *(ckt->CKTstate0 + here->CNTqgd);
            return(OK);
        case CNT_CQGD:
            value->rValue = *(ckt->CKTstate0 + here->CNTcqgd);
            return(OK);
        case CNT_CAPGB:
            value->rValue = 2* *(ckt->CKTstate0 + here->CNTcapgb);
            /* add overlap capacitance */
            value->rValue += (here->sCNTmodPtr->CNTgateBulkOverlapCapFactor)
                             * here->CNTm
                             * (here->CNTl
                                -2*(here->sCNTmodPtr->CNTlatDiff));
            return(OK);
        case CNT_QGB:
            value->rValue = *(ckt->CKTstate0 + here->CNTqgb);
            return(OK);
        case CNT_CQGB:
            value->rValue = *(ckt->CKTstate0 + here->CNTcqgb);
            return(OK);
        case CNT_QBD:
            value->rValue = *(ckt->CKTstate0 + here->CNTqbd);
            return(OK);
        case CNT_CQBD:
            value->rValue = *(ckt->CKTstate0 + here->CNTcqbd);
            return(OK);
        case CNT_QBS:
            value->rValue = *(ckt->CKTstate0 + here->CNTqbs);
            return(OK);
        case CNT_CQBS:
            value->rValue = *(ckt->CKTstate0 + here->CNTcqbs);
            return(OK);
        case CNT_L_SENS_DC:
            if(ckt->CKTsenInfo && here->CNTsens_l){
               value->rValue = *(ckt->CKTsenInfo->SEN_Sap[select->iValue + 1]+
                       here->CNTsenParmNo);
            }
            return(OK);
        case CNT_L_SENS_REAL:
            if(ckt->CKTsenInfo && here->CNTsens_l){
               value->rValue = *(ckt->CKTsenInfo->SEN_RHS[select->iValue + 1]+
                       here->CNTsenParmNo);
            }
            return(OK);
        case CNT_L_SENS_IMAG:
            if(ckt->CKTsenInfo && here->CNTsens_l){
               value->rValue = *(ckt->CKTsenInfo->SEN_iRHS[select->iValue + 1]+
                       here->CNTsenParmNo);
            }
            return(OK);
        case CNT_L_SENS_MAG:
            if(ckt->CKTsenInfo && here->CNTsens_l){
                vr = *(ckt->CKTrhsOld + select->iValue + 1); 
                vi = *(ckt->CKTirhsOld + select->iValue + 1); 
                vm = sqrt(vr*vr + vi*vi);
                if(vm == 0){
                    value->rValue = 0;
                    return(OK);
                }
                sr = *(ckt->CKTsenInfo->SEN_RHS[select->iValue + 1]+
                        here->CNTsenParmNo);
                si = *(ckt->CKTsenInfo->SEN_iRHS[select->iValue + 1]+
                        here->CNTsenParmNo);
                value->rValue = (vr * sr + vi * si)/vm;
            }
            return(OK);
        case CNT_L_SENS_PH:
            if(ckt->CKTsenInfo && here->CNTsens_l){
                vr = *(ckt->CKTrhsOld + select->iValue + 1); 
                vi = *(ckt->CKTirhsOld + select->iValue + 1); 
                vm = vr*vr + vi*vi;
                if(vm == 0){
                    value->rValue = 0;
                    return(OK);
                }
                sr = *(ckt->CKTsenInfo->SEN_RHS[select->iValue + 1]+
                        here->CNTsenParmNo);
                si = *(ckt->CKTsenInfo->SEN_iRHS[select->iValue + 1]+
                        here->CNTsenParmNo);
                value->rValue =  (vr * si - vi * sr)/vm;
            }
            return(OK);
        case CNT_L_SENS_CPLX:
            if(ckt->CKTsenInfo && here->CNTsens_l){
                value->cValue.real= 
                        *(ckt->CKTsenInfo->SEN_RHS[select->iValue + 1]+
                        here->CNTsenParmNo);
                value->cValue.imag= 
                        *(ckt->CKTsenInfo->SEN_iRHS[select->iValue + 1]+
                        here->CNTsenParmNo);
            }
            return(OK);
        case CNT_W_SENS_DC:
            if(ckt->CKTsenInfo && here->CNTsens_w){
                value->rValue = *(ckt->CKTsenInfo->SEN_Sap[select->iValue + 1]+
                        here->CNTsenParmNo + here->CNTsens_l);
            }
            return(OK);
        case CNT_W_SENS_REAL:
            if(ckt->CKTsenInfo && here->CNTsens_w){
                value->rValue = *(ckt->CKTsenInfo->SEN_RHS[select->iValue + 1]+
                        here->CNTsenParmNo + here->CNTsens_l);
            }
             return(OK);
        case CNT_W_SENS_IMAG:
            if(ckt->CKTsenInfo && here->CNTsens_w){
                value->rValue = *(ckt->CKTsenInfo->SEN_iRHS[select->iValue + 1]+
                        here->CNTsenParmNo + here->CNTsens_l);
            }
            return(OK);
        case CNT_W_SENS_MAG:
            if(ckt->CKTsenInfo && here->CNTsens_w){
                vr = *(ckt->CKTrhsOld + select->iValue + 1); 
                vi = *(ckt->CKTirhsOld + select->iValue + 1); 
                vm = sqrt(vr*vr + vi*vi);
                if(vm == 0){
                    value->rValue = 0;
                    return(OK);
                }
                sr = *(ckt->CKTsenInfo->SEN_RHS[select->iValue + 1]+
                        here->CNTsenParmNo + here->CNTsens_l);
                si = *(ckt->CKTsenInfo->SEN_iRHS[select->iValue + 1]+
                        here->CNTsenParmNo + here->CNTsens_l);
                value->rValue = (vr * sr + vi * si)/vm;
            }
            return(OK);
        case CNT_W_SENS_PH:
            if(ckt->CKTsenInfo && here->CNTsens_w){
                vr = *(ckt->CKTrhsOld + select->iValue + 1); 
                vi = *(ckt->CKTirhsOld + select->iValue + 1);     
                vm = vr*vr + vi*vi;
                if(vm == 0){
                    value->rValue = 0;
                    return(OK);
                }
                sr = *(ckt->CKTsenInfo->SEN_RHS[select->iValue + 1]+
                        here->CNTsenParmNo + here->CNTsens_l);
                si = *(ckt->CKTsenInfo->SEN_iRHS[select->iValue + 1]+
                        here->CNTsenParmNo + here->CNTsens_l);
                value->rValue =  (vr * si - vi * sr)/vm;
            }
                    return(OK);
        case CNT_W_SENS_CPLX:
            if(ckt->CKTsenInfo && here->CNTsens_w){
                value->cValue.real= 
                        *(ckt->CKTsenInfo->SEN_RHS[select->iValue + 1]+
                        here->CNTsenParmNo + here->CNTsens_l);
                value->cValue.imag= 
                        *(ckt->CKTsenInfo->SEN_iRHS[select->iValue + 1]+
                        here->CNTsenParmNo + here->CNTsens_l);
            }
            return(OK);
        case CNT_CB :
            if (ckt->CKTcurrentAnalysis & DOING_AC) {
                errMsg = MALLOC(strlen(msg)+1);
                errRtn = "CNTask.c";
                strcpy(errMsg,msg);
                return(E_ASKCURRENT);
            } else {
                value->rValue = here->CNTcbd + here->CNTcbs - *(ckt->CKTstate0
                        + here->CNTcqgb);
            }
            return(OK);
        case CNT_CG :
            if (ckt->CKTcurrentAnalysis & DOING_AC) {
                errMsg = MALLOC(strlen(msg)+1);
                errRtn = "CNTask.c";
                strcpy(errMsg,msg);
                return(E_ASKCURRENT);
            } else if (ckt->CKTcurrentAnalysis & (DOING_DCOP | DOING_TRCV)) {
                value->rValue = 0;
            } else if ((ckt->CKTcurrentAnalysis & DOING_TRAN) && 
                    (ckt->CKTmode & MODETRANOP)) {
                value->rValue = 0;
            } else {
                value->rValue = *(ckt->CKTstate0 + here->CNTcqgb) +
                        *(ckt->CKTstate0 + here->CNTcqgd) + *(ckt->CKTstate0 + 
                        here->CNTcqgs);
            }
            return(OK);
        case CNT_CS :
            if (ckt->CKTcurrentAnalysis & DOING_AC) {
                errMsg = MALLOC(strlen(msg)+1);
                errRtn = "CNTask.c";
                strcpy(errMsg,msg);
                return(E_ASKCURRENT);
            } else {
                value->rValue = -here->CNTcd;
                value->rValue -= here->CNTcbd + here->CNTcbs -
                        *(ckt->CKTstate0 + here->CNTcqgb);
                if ((ckt->CKTcurrentAnalysis & DOING_TRAN) && 
                        !(ckt->CKTmode & MODETRANOP)) {
                    value->rValue -= *(ckt->CKTstate0 + here->CNTcqgb) + 
                            *(ckt->CKTstate0 + here->CNTcqgd) +
                            *(ckt->CKTstate0 + here->CNTcqgs);
                }
            }
            return(OK);
        case CNT_POWER :
            if (ckt->CKTcurrentAnalysis & DOING_AC) {
                errMsg = MALLOC(strlen(msg)+1);
                errRtn = "CNTask.c";
                strcpy(errMsg,msg);
                return(E_ASKPOWER);
            } else {
                double temp;

                value->rValue = here->CNTcd * 
                        *(ckt->CKTrhsOld + here->CNTdNode);
                value->rValue += (here->CNTcbd + here->CNTcbs -
                        *(ckt->CKTstate0 + here->CNTcqgb)) *
                        *(ckt->CKTrhsOld + here->CNTbNode);
                if ((ckt->CKTcurrentAnalysis & DOING_TRAN) && 
                        !(ckt->CKTmode & MODETRANOP)) {
                    value->rValue += (*(ckt->CKTstate0 + here->CNTcqgb) + 
                            *(ckt->CKTstate0 + here->CNTcqgd) +
                            *(ckt->CKTstate0 + here->CNTcqgs)) *
                            *(ckt->CKTrhsOld + here->CNTgNode);
                }
                temp = -here->CNTcd;
                temp -= here->CNTcbd + here->CNTcbs ;
                if ((ckt->CKTcurrentAnalysis & DOING_TRAN) && 
                        !(ckt->CKTmode & MODETRANOP)) {
                    temp -= *(ckt->CKTstate0 + here->CNTcqgb) + 
                            *(ckt->CKTstate0 + here->CNTcqgd) + 
                            *(ckt->CKTstate0 + here->CNTcqgs);
                }
                value->rValue += temp * *(ckt->CKTrhsOld + here->CNTsNode);
            }
            return(OK);
        default:
            return(E_BADPARM);
    }
    /* NOTREACHED */
}

