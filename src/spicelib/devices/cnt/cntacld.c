/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/
/*
 */

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "cntdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"


int
CNTacLoad(GENmodel *inModel, CKTcircuit *ckt)
{
    CNTmodel *model = (CNTmodel*)inModel;
    CNTinstance *here;
    int xnrm;
    int xrev;
    double xgs;
    double xgd;
    double xgb;
    double xbd;
    double xbs;
    double capgs;
    double capgd;
    double capgb;
    double GateBulkOverlapCap;
    double GateDrainOverlapCap;
    double GateSourceOverlapCap;
    double EffectiveLength;

    for( ; model != NULL; model = model->CNTnextModel) {
        for(here = model->CNTinstances; here!= NULL;
                here = here->CNTnextInstance) {
	    /*if (here->CNTowner != ARCHme) continue;*//* Cuidado acabei de mudar!!*/
        
            if (here->CNTmode < 0) {
                xnrm=0;
                xrev=1;
            } else {
                xnrm=1;
                xrev=0;
            }
            /*
             *     meyer's model parameters
             */
            EffectiveLength=here->CNTl - 2*model->CNTlatDiff;
            
            GateSourceOverlapCap = model->CNTgateSourceOverlapCapFactor * 
                    here->CNTm * here->CNTw;
            GateDrainOverlapCap = model->CNTgateDrainOverlapCapFactor * 
                    here->CNTm * here->CNTw;
            GateBulkOverlapCap = model->CNTgateBulkOverlapCapFactor * 
                    here->CNTm * EffectiveLength;
                    
            capgs = ( *(ckt->CKTstate0+here->CNTcapgs)+ 
                      *(ckt->CKTstate0+here->CNTcapgs) +
                      GateSourceOverlapCap );
            capgd = ( *(ckt->CKTstate0+here->CNTcapgd)+ 
                      *(ckt->CKTstate0+here->CNTcapgd) +
                      GateDrainOverlapCap );
            capgb = ( *(ckt->CKTstate0+here->CNTcapgb)+ 
                      *(ckt->CKTstate0+here->CNTcapgb) +
                      GateBulkOverlapCap );
            xgs = capgs * ckt->CKTomega;
            xgd = capgd * ckt->CKTomega;
            xgb = capgb * ckt->CKTomega;
            xbd  = here->CNTcapbd * ckt->CKTomega;
            xbs  = here->CNTcapbs * ckt->CKTomega;
            /*
             *    load matrix
             */

            *(here->CNTGgPtr +1) += xgd+xgs+xgb;
            *(here->CNTBbPtr +1) += xgb+xbd+xbs;
            *(here->CNTDPdpPtr +1) += xgd+xbd;
            *(here->CNTSPspPtr +1) += xgs+xbs;
            *(here->CNTGbPtr +1) -= xgb;
            *(here->CNTGdpPtr +1) -= xgd;
            *(here->CNTGspPtr +1) -= xgs;
            *(here->CNTBgPtr +1) -= xgb;
            *(here->CNTBdpPtr +1) -= xbd;
            *(here->CNTBspPtr +1) -= xbs;
            *(here->CNTDPgPtr +1) -= xgd;
            *(here->CNTDPbPtr +1) -= xbd;
            *(here->CNTSPgPtr +1) -= xgs;
            *(here->CNTSPbPtr +1) -= xbs;
            *(here->CNTDdPtr) += here->CNTdrainConductance;
            *(here->CNTSsPtr) += here->CNTsourceConductance;
            *(here->CNTBbPtr) += here->CNTgbd+here->CNTgbs;
            *(here->CNTDPdpPtr) += here->CNTdrainConductance+
                    here->CNTgds+here->CNTgbd+
                    xrev*(here->CNTgm+here->CNTgmbs);
            *(here->CNTSPspPtr) += here->CNTsourceConductance+
                    here->CNTgds+here->CNTgbs+
                    xnrm*(here->CNTgm+here->CNTgmbs);
            *(here->CNTDdpPtr) -= here->CNTdrainConductance;
            *(here->CNTSspPtr) -= here->CNTsourceConductance;
            *(here->CNTBdpPtr) -= here->CNTgbd;
            *(here->CNTBspPtr) -= here->CNTgbs;
            *(here->CNTDPdPtr) -= here->CNTdrainConductance;
            *(here->CNTDPgPtr) += (xnrm-xrev)*here->CNTgm;
            *(here->CNTDPbPtr) += -here->CNTgbd+(xnrm-xrev)*here->CNTgmbs;
            *(here->CNTDPspPtr) -= here->CNTgds+
                    xnrm*(here->CNTgm+here->CNTgmbs);
            *(here->CNTSPgPtr) -= (xnrm-xrev)*here->CNTgm;
            *(here->CNTSPsPtr) -= here->CNTsourceConductance;
            *(here->CNTSPbPtr) -= here->CNTgbs+(xnrm-xrev)*here->CNTgmbs;
            *(here->CNTSPdpPtr) -= here->CNTgds+
                    xrev*(here->CNTgm+here->CNTgmbs);

        }
    }
    return(OK);
}
