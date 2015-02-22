/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/
/*
 */

#include "ngspice.h"
#include "cktdefs.h"
#include "complex.h"
#include "cntdefs.h"
#include "sperror.h"
#include "suffix.h"


int
CNTpzLoad(GENmodel *inModel, CKTcircuit *ckt, SPcomplex *s)
{
    CNTmodel *model = (CNTmodel *)inModel;
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
	    if (here->CNTowner != ARCHme) continue;
        
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
            
            capgs = ( 2* *(ckt->CKTstate0+here->CNTcapgs)+ 
                      GateSourceOverlapCap );
            capgd = ( 2* *(ckt->CKTstate0+here->CNTcapgd)+ 
                      GateDrainOverlapCap );
            capgb = ( 2* *(ckt->CKTstate0+here->CNTcapgb)+ 
                      GateBulkOverlapCap );
            xgs = capgs;
            xgd = capgd;
            xgb = capgb;
            xbd  = here->CNTcapbd;
            xbs  = here->CNTcapbs;
            /*printf("mos2: xgs=%g, xgd=%g, xgb=%g, xbd=%g, xbs=%g\n",
                    xgs,xgd,xgb,xbd,xbs);*/
            /*
             *    load matrix
             */

            *(here->CNTGgPtr   ) += (xgd+xgs+xgb)*s->real;
            *(here->CNTGgPtr +1) += (xgd+xgs+xgb)*s->imag;
            *(here->CNTBbPtr   ) += (xgb+xbd+xbs)*s->real;
            *(here->CNTBbPtr +1) += (xgb+xbd+xbs)*s->imag;
            *(here->CNTDPdpPtr   ) += (xgd+xbd)*s->real;
            *(here->CNTDPdpPtr +1) += (xgd+xbd)*s->imag;
            *(here->CNTSPspPtr   ) += (xgs+xbs)*s->real;
            *(here->CNTSPspPtr +1) += (xgs+xbs)*s->imag;
            *(here->CNTGbPtr   ) -= xgb*s->real;
            *(here->CNTGbPtr +1) -= xgb*s->imag;
            *(here->CNTGdpPtr   ) -= xgd*s->real;
            *(here->CNTGdpPtr +1) -= xgd*s->imag;
            *(here->CNTGspPtr   ) -= xgs*s->real;
            *(here->CNTGspPtr +1) -= xgs*s->imag;
            *(here->CNTBgPtr   ) -= xgb*s->real;
            *(here->CNTBgPtr +1) -= xgb*s->imag;
            *(here->CNTBdpPtr   ) -= xbd*s->real;
            *(here->CNTBdpPtr +1) -= xbd*s->imag;
            *(here->CNTBspPtr   ) -= xbs*s->real;
            *(here->CNTBspPtr +1) -= xbs*s->imag;
            *(here->CNTDPgPtr   ) -= xgd*s->real;
            *(here->CNTDPgPtr +1) -= xgd*s->imag;
            *(here->CNTDPbPtr   ) -= xbd*s->real;
            *(here->CNTDPbPtr +1) -= xbd*s->imag;
            *(here->CNTSPgPtr   ) -= xgs*s->real;
            *(here->CNTSPgPtr +1) -= xgs*s->imag;
            *(here->CNTSPbPtr   ) -= xbs*s->real;
            *(here->CNTSPbPtr +1) -= xbs*s->imag;
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
