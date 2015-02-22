/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/

#include "ngspice.h"
#include "cktdefs.h"
#include "cntdefs.h"
#include "const.h"
#include "sperror.h"
#include "suffix.h"

int
CNTtemp(GENmodel *inModel, CKTcircuit *ckt)
{
    CNTmodel *model = (CNTmodel *)inModel;
    CNTinstance *here;

    double egfet,egfet1;
    double fact1,fact2;
    double kt,kt1;
    double arg1;
    double ratio,ratio4;
    double phio;
    double pbo;
    double gmanew,gmaold;
    double capfact;
    double pbfact1,pbfact;
    double vt,vtnom;
    double wkfngs;
    double wkfng;
    double fermig;
    double fermis;
    double vfb;
    /* loop through all the resistor models */
    for( ; model != NULL; model = model->CNTnextModel) {
        

        /* perform model defaulting */
        if(!model->CNTtnomGiven) {
            model->CNTtnom = ckt->CKTnomTemp;
        }

        fact1 = model->CNTtnom/REFTEMP;
        vtnom = model->CNTtnom*CONSTKoverQ;
        kt1 = CONSTboltz * model->CNTtnom;
        egfet1 = 1.16-(7.02e-4*model->CNTtnom*model->CNTtnom)/
                (model->CNTtnom+1108);
        arg1 = -egfet1/(kt1+kt1)+1.1150877/(CONSTboltz*(REFTEMP+REFTEMP));
        pbfact1 = -2*vtnom *(1.5*log(fact1)+CHARGE*arg1);

    /* now model parameter preprocessing */

        if(!model->CNToxideThicknessGiven || model->CNToxideThickness == 0) {
            model->CNToxideCapFactor = 0;
        } else {
            model->CNToxideCapFactor = 3.9 * 8.854214871e-12/
                    model->CNToxideThickness;
            if(!model->CNTtransconductanceGiven) {
                if(!model->CNTsurfaceMobilityGiven) {
                    model->CNTsurfaceMobility=600;
                }
                model->CNTtransconductance = model->CNTsurfaceMobility *
                        model->CNToxideCapFactor * 1e-4 /*(m**2/cm**2)*/;
            }
            if(model->CNTsubstrateDopingGiven) {
                if(model->CNTsubstrateDoping*1e6 /*(cm**3/m**3)*/ >1.45e16) {
                    if(!model->CNTphiGiven) {
                        model->CNTphi = 2*vtnom*
                                log(model->CNTsubstrateDoping*
                                1e6/*(cm**3/m**3)*//1.45e16);
                        model->CNTphi = MAX(.1,model->CNTphi);
                    }
                    fermis = model->CNTtype * .5 * model->CNTphi;
                    wkfng = 3.2;
                    if(!model->CNTgateTypeGiven) model->CNTgateType=1;
                    if(model->CNTgateType != 0) {
                        fermig = model->CNTtype *model->CNTgateType*.5*egfet1;
                        wkfng = 3.25 + .5 * egfet1 - fermig;
                    }
                    wkfngs = wkfng - (3.25 + .5 * egfet1 +fermis);
                    if(!model->CNTgammaGiven) {
                        model->CNTgamma = sqrt(2 * 11.70 * 8.854214871e-12 * 
                                CHARGE * model->CNTsubstrateDoping*
                                1e6/*(cm**3/m**3)*/)/
                                model->CNToxideCapFactor;
                    }
                    if(!model->CNTvt0Given) {
                        if(!model->CNTsurfaceStateDensityGiven) 
                                model->CNTsurfaceStateDensity=0;
                        vfb = wkfngs - 
                                model->CNTsurfaceStateDensity * 
                                1e4 /*(cm**2/m**2)*/ * 
                                CHARGE/model->CNToxideCapFactor;
                        model->CNTvt0 = vfb + model->CNTtype * 
                                (model->CNTgamma * sqrt(model->CNTphi)+
                                model->CNTphi);
                    }
                } else {
                    model->CNTsubstrateDoping = 0;
                    (*(SPfrontEnd->IFerror))(ERR_FATAL,
                            "%s: Nsub < Ni",&model->CNTmodName);
                    return(E_BADPARM);
                }
            }
        }

        
        /* loop through all instances of the model */
        for(here = model->CNTinstances; here!= NULL; 
                here = here->CNTnextInstance) {
            double czbd;    /* zero voltage bulk-drain capacitance */
            double czbdsw;  /* zero voltage bulk-drain sidewall capacitance */
            double czbs;    /* zero voltage bulk-source capacitance */
            double czbssw;  /* zero voltage bulk-source sidewall capacitance */
            double arg;     /* 1 - fc */
            double sarg;    /* (1-fc) ^^ (-mj) */
            double sargsw;  /* (1-fc) ^^ (-mjsw) */
	    if (here->CNTowner != ARCHme) continue;

            /* perform the parameter defaulting */
            
            if(!here->CNTdtempGiven) {
                here->CNTdtemp = 0.0;
            }
            if(!here->CNTtempGiven) {
                here->CNTtemp = ckt->CKTtemp + here->CNTdtemp;
            }
            vt = here->CNTtemp * CONSTKoverQ;
            ratio = here->CNTtemp/model->CNTtnom;
            fact2 = here->CNTtemp/REFTEMP;
            kt = here->CNTtemp * CONSTboltz;
            egfet = 1.16-(7.02e-4*here->CNTtemp*here->CNTtemp)/
                    (here->CNTtemp+1108);
            arg = -egfet/(kt+kt)+1.1150877/(CONSTboltz*(REFTEMP+REFTEMP));
            pbfact = -2*vt *(1.5*log(fact2)+CHARGE*arg);

            if(!here->CNTdrainAreaGiven) {
                here->CNTdrainArea = ckt->CKTdefaultMosAD;
            }
            if(!here->CNTmGiven) {
                here->CNTm = ckt->CKTdefaultMosM;
            }
            if(!here->CNTlGiven) {
                here->CNTl = ckt->CKTdefaultMosL;
            }
            if(!here->CNTsourceAreaGiven) {
                here->CNTsourceArea = ckt->CKTdefaultMosAS;
            }
            if(!here->CNTwGiven) {
                here->CNTw = ckt->CKTdefaultMosW;
            }

            if(here->CNTl - 2 * model->CNTlatDiff <=0) {
                (*(SPfrontEnd->IFerror))(ERR_WARNING,
                        "%s: effective channel length less than zero",
                        &(model->CNTmodName));
            }
            ratio4 = ratio * sqrt(ratio);
            here->CNTtTransconductance = model->CNTtransconductance / ratio4;
            here->CNTtSurfMob = model->CNTsurfaceMobility/ratio4;
            phio= (model->CNTphi-pbfact1)/fact1;
            here->CNTtPhi = fact2 * phio + pbfact;
            here->CNTtVbi = 
                    model->CNTvt0 - model->CNTtype * 
                        (model->CNTgamma* sqrt(model->CNTphi))
                    +.5*(egfet1-egfet) 
                    + model->CNTtype*.5* (here->CNTtPhi-model->CNTphi);
            here->CNTtVto = here->CNTtVbi + model->CNTtype * 
                    model->CNTgamma * sqrt(here->CNTtPhi);
            here->CNTtSatCur = model->CNTjctSatCur* 
                    exp(-egfet/vt+egfet1/vtnom);
            here->CNTtSatCurDens = model->CNTjctSatCurDensity *
                    exp(-egfet/vt+egfet1/vtnom);
            pbo = (model->CNTbulkJctPotential - pbfact1)/fact1;
            gmaold = (model->CNTbulkJctPotential-pbo)/pbo;
            capfact = 1/(1+model->CNTbulkJctBotGradingCoeff*
                    (4e-4*(model->CNTtnom-REFTEMP)-gmaold));
            here->CNTtCbd = model->CNTcapBD * capfact;
            here->CNTtCbs = model->CNTcapBS * capfact;
            here->CNTtCj = model->CNTbulkCapFactor * capfact;
            capfact = 1/(1+model->CNTbulkJctSideGradingCoeff*
                    (4e-4*(model->CNTtnom-REFTEMP)-gmaold));
            here->CNTtCjsw = model->CNTsideWallCapFactor * capfact;
            here->CNTtBulkPot = fact2 * pbo+pbfact;
            gmanew = (here->CNTtBulkPot-pbo)/pbo;
            capfact = (1+model->CNTbulkJctBotGradingCoeff*
                    (4e-4*(here->CNTtemp-REFTEMP)-gmanew));
            here->CNTtCbd *= capfact;
            here->CNTtCbs *= capfact;
            here->CNTtCj *= capfact;
            capfact = (1+model->CNTbulkJctSideGradingCoeff*
                    (4e-4*(here->CNTtemp-REFTEMP)-gmanew));
            here->CNTtCjsw *= capfact;
            here->CNTtDepCap = model->CNTfwdCapDepCoeff * here->CNTtBulkPot;
            if( (here->CNTtSatCurDens == 0) ||
                    (here->CNTdrainArea == 0) ||
                    (here->CNTsourceArea == 0) ) {
                here->CNTsourceVcrit = here->CNTdrainVcrit =
                       vt*log(vt/(CONSTroot2*here->CNTm*here->CNTtSatCur));
            } else {
                here->CNTdrainVcrit =
                        vt * log( vt / (CONSTroot2 *
                        here->CNTm *
                        here->CNTtSatCurDens * here->CNTdrainArea));
                here->CNTsourceVcrit =
                        vt * log( vt / (CONSTroot2 *
                        here->CNTm *
                        here->CNTtSatCurDens * here->CNTsourceArea));
            }

            if(model->CNTcapBDGiven) {
                czbd = here->CNTtCbd * here->CNTm;
            } else {
                if(model->CNTbulkCapFactorGiven) {  
                    czbd=here->CNTtCj*here->CNTm*here->CNTdrainArea;
                } else {
                    czbd=0;
                }
            }
            if(model->CNTsideWallCapFactorGiven) {
                czbdsw= here->CNTtCjsw * here->CNTdrainPerimiter *
                     here->CNTm;
            } else {
                czbdsw=0;
            }
            arg = 1-model->CNTfwdCapDepCoeff;
            sarg = exp( (-model->CNTbulkJctBotGradingCoeff) * log(arg) );
            sargsw = exp( (-model->CNTbulkJctSideGradingCoeff) * log(arg) );
            here->CNTCbd = czbd;
            here->CNTCbdsw = czbdsw;
            here->CNTf2d = czbd*(1-model->CNTfwdCapDepCoeff*
                        (1+model->CNTbulkJctBotGradingCoeff))* sarg/arg
                    +  czbdsw*(1-model->CNTfwdCapDepCoeff*
                        (1+model->CNTbulkJctSideGradingCoeff))*
                        sargsw/arg;
            here->CNTf3d = czbd * model->CNTbulkJctBotGradingCoeff * sarg/arg/
                        here->CNTtBulkPot
                    + czbdsw * model->CNTbulkJctSideGradingCoeff * sargsw/arg /
                        here->CNTtBulkPot;
            here->CNTf4d = czbd*here->CNTtBulkPot*(1-arg*sarg)/
                        (1-model->CNTbulkJctBotGradingCoeff)
                    + czbdsw*here->CNTtBulkPot*(1-arg*sargsw)/
                        (1-model->CNTbulkJctSideGradingCoeff)
                    -here->CNTf3d/2*
                        (here->CNTtDepCap*here->CNTtDepCap)
                    -here->CNTtDepCap * here->CNTf2d;
            if(model->CNTcapBSGiven) {
                czbs=here->CNTtCbs * here->CNTm;
            } else {
                if(model->CNTbulkCapFactorGiven) {
                   czbs=here->CNTtCj*here->CNTsourceArea * here->CNTm;
                } else {
                    czbs=0;
                }
            }
            if(model->CNTsideWallCapFactorGiven) {
                czbssw = here->CNTtCjsw * here->CNTsourcePerimiter *
                          here->CNTm;
            } else {
                czbssw=0;
            }
            arg = 1-model->CNTfwdCapDepCoeff;
            sarg = exp( (-model->CNTbulkJctBotGradingCoeff) * log(arg) );
            sargsw = exp( (-model->CNTbulkJctSideGradingCoeff) * log(arg) );
            here->CNTCbs = czbs;
            here->CNTCbssw = czbssw;
            here->CNTf2s = czbs*(1-model->CNTfwdCapDepCoeff*
                        (1+model->CNTbulkJctBotGradingCoeff))* sarg/arg
                    +  czbssw*(1-model->CNTfwdCapDepCoeff*
                        (1+model->CNTbulkJctSideGradingCoeff))*
                        sargsw/arg;
            here->CNTf3s = czbs * model->CNTbulkJctBotGradingCoeff * sarg/arg/
                        here->CNTtBulkPot
                    + czbssw * model->CNTbulkJctSideGradingCoeff * sargsw/arg /
                        here->CNTtBulkPot;
            here->CNTf4s = czbs*here->CNTtBulkPot*(1-arg*sarg)/
                        (1-model->CNTbulkJctBotGradingCoeff)
                    + czbssw*here->CNTtBulkPot*(1-arg*sargsw)/
                        (1-model->CNTbulkJctSideGradingCoeff)
                    -here->CNTf3s/2*
                        (here->CNTtDepCap*here->CNTtDepCap)
                    -here->CNTtDepCap * here->CNTf2s;


            if(model->CNTdrainResistanceGiven) {
                if(model->CNTdrainResistance != 0) {
                   here->CNTdrainConductance = here->CNTm /
                                      model->CNTdrainResistance;
                } else {
                    here->CNTdrainConductance = 0;
                }
            } else if (model->CNTsheetResistanceGiven) {
                if(model->CNTsheetResistance != 0) {
                    here->CNTdrainConductance = 
                       here->CNTm /
                          (model->CNTsheetResistance*here->CNTdrainSquares);
                } else {
                    here->CNTdrainConductance = 0;
                }
            } else {
                here->CNTdrainConductance = 0;
            }
            if(model->CNTsourceResistanceGiven) {
                if(model->CNTsourceResistance != 0) {
                   here->CNTsourceConductance = here->CNTm /
                                         model->CNTsourceResistance;
                } else {
                    here->CNTsourceConductance = 0;
                }
            } else if (model->CNTsheetResistanceGiven) {
                if ((model->CNTsheetResistance != 0) &&
                                   (here->CNTsourceSquares != 0)) {
                    here->CNTsourceConductance = 
                        here->CNTm /
                          (model->CNTsheetResistance*here->CNTsourceSquares);
                } else {
                    here->CNTsourceConductance = 0;
                }
            } else {
                here->CNTsourceConductance = 0;
            }
        }
    }
    return(OK);
}
