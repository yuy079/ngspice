/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles

This function is obsolete (was used by an old sensitivity analysis)
**********/

/* actually load the current sensitivity 
 * information into the  array previously provided 
 */

#include "ngspice/ngspice.h"
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "cntdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int
CNTsLoad(GENmodel *inModel, CKTcircuit *ckt)
{
    CNTmodel *model = (CNTmodel *)inModel;
    CNTinstance *here;
    double   SaveState[44];
    int    save_mode;
    int    i;
    int    iparmno;
    int    error;
    int    flag;
    double A0;
    double DELA;
    double Apert;
    double DELAinv;
    double gspr0;
    double gspr;
    double gdpr0;
    double gdpr;
    double cdpr0;
    double cspr0;
    double cd0;
    double cbd0;
    double cbs0;
    double cd;
    double cbd;
    double cbs;
    double DcdprDp;
    double DcsprDp;
    double DcbDp;
    double DcdDp;
    double DcbsDp;
    double DcbdDp;
    double DcdprmDp;
    double DcsprmDp;
    double qgs0;
    double qgd0;
    double qgb0;
    double qbd0;
    double qbd;
    double qbs0;
    double qbs;
    double DqgsDp;
    double DqgdDp;
    double DqgbDp;
    double DqbdDp;
    double DqbsDp;
    double Osxpgs;
    double Osxpgd;
    double Osxpgb;
    double Osxpbd;
    double Osxpbs;
    double tag0;
    double tag1;
    double arg;
    double sarg;
    double sargsw;
    int offset;
    double EffectiveLength;
    SENstruct *info;

#ifdef SENSDEBUG
    printf("CNTsenload \n");
    printf("CKTtime = %.5e\n",ckt->CKTtime);
    printf("CKTorder = %d\n",ckt->CKTorder);
#endif /* SENSDEBUG */

    info = ckt->CKTsenInfo;
    info->SENstatus = PERTURBATION;

    tag0 = ckt->CKTag[0]; 
    tag1 = ckt->CKTag[1];
    if(ckt->CKTorder == 1){
        tag1 = 0;
    }

    /*  loop through all the models */
    for( ; model != NULL; model = model->CNTnextModel ) {

        /* loop through all the instances of the model */
        for (here = model->CNTinstances; here != NULL ;
                here=here->CNTnextInstance) {
/*	    if (here->CNTowner != ARCHme) continue;*/

#ifdef SENSDEBUG
            printf("senload instance name %s\n",here->CNTname);
            printf("gate = %d ,drain = %d, drainprm = %d\n", 
                    here->CNTgNode,here->CNTdNode,here->CNTdNodePrime);
            printf("source = %d , sourceprm = %d ,body = %d, senparmno = %d\n",
                    here->CNTsNode ,here->CNTsNodePrime,
                    here->CNTbNode,here->CNTsenParmNo);
#endif /* SENSDEBUG */


            /* save the unperturbed values in the state vector */
            for(i=0; i <= 16; i++){
                *(SaveState + i) = *(ckt->CKTstate0 + here->CNTstates + i);
            }

            *(SaveState + 17) = here->CNTsourceConductance;  
            *(SaveState + 18) = here->CNTdrainConductance;  
            *(SaveState + 19) = here->CNTcd;  
            *(SaveState + 20) = here->CNTcbs;  
            *(SaveState + 21) = here->CNTcbd;  
            *(SaveState + 22) = here->CNTgmbs;  
            *(SaveState + 23) = here->CNTgm;  
            *(SaveState + 24) = here->CNTgds;  
            *(SaveState + 25) = here->CNTgbd;  
            *(SaveState + 26) = here->CNTgbs;  
            *(SaveState + 27) = here->CNTcapbd;  
            *(SaveState + 28) = here->CNTcapbs;  
            *(SaveState + 29) = here->CNTCbd;  
            *(SaveState + 30) = here->CNTCbdsw;  
            *(SaveState + 31) = here->CNTCbs;  
            *(SaveState + 32) = here->CNTCbssw;  
            *(SaveState + 33) = here->CNTf2d;  
            *(SaveState + 34) = here->CNTf3d;  
            *(SaveState + 35) = here->CNTf4d;  
            *(SaveState + 36) = here->CNTf2s;  
            *(SaveState + 37) = here->CNTf3s;  
            *(SaveState + 38) = here->CNTf4s;  
            *(SaveState + 39) = here->CNTcgs;  
            *(SaveState + 40) = here->CNTcgd;  
            *(SaveState + 41) = here->CNTcgb;  
            *(SaveState + 42) = here->CNTvdsat;  
            *(SaveState + 43) = here->CNTvon;  
            save_mode  = here->CNTmode;  


            if(here->CNTsenParmNo == 0) goto next1;

#ifdef SENSDEBUG
            printf("without perturbation \n");
            printf("gbd =%.5e\n",here->CNTgbd);
            printf("satCur =%.5e\n",here->CNTtSatCur);
            printf("satCurDens =%.5e\n",here->CNTtSatCurDens);
            printf("vbd =%.5e\n",*(ckt->CKTstate0 + here->CNTvbd));
#endif /* SENSDEBUG */

            cdpr0= here->CNTcd;
            cspr0= -(here->CNTcd + here->CNTcbd + here->CNTcbs);
            if((info->SENmode == TRANSEN) &&
                (ckt->CKTmode & MODEINITTRAN)){
                qgs0 = *(ckt->CKTstate1 + here->CNTqgs);
                qgd0 = *(ckt->CKTstate1 + here->CNTqgd);
                qgb0 = *(ckt->CKTstate1 + here->CNTqgb);
            }
            else{
                qgs0 = *(ckt->CKTstate0 + here->CNTqgs);
                qgd0 = *(ckt->CKTstate0 + here->CNTqgd);
                qgb0 = *(ckt->CKTstate0 + here->CNTqgb);
            }

            here->CNTsenPertFlag = ON;
            error =  CNTload((GENmodel*)model,ckt);
            if(error) return(error);

            cd0 =  here->CNTcd ;
            cbd0 = here->CNTcbd ;
            cbs0 = here->CNTcbs ;
            gspr0= here->CNTsourceConductance ;
            gdpr0= here->CNTdrainConductance ;

            qbs0 = *(ckt->CKTstate0 + here->CNTqbs);
            qbd0 = *(ckt->CKTstate0 + here->CNTqbd);

            for( flag = 0 ; flag <= 1 ; flag++){
                if(here->CNTsens_l == 0)
                    if(flag == 0) goto next2;
                if(here->CNTsens_w == 0)
                    if(flag == 1) goto next2;
                if(flag == 0){
                    A0 = here->CNTl;
                    DELA = info->SENpertfac * A0;
                    DELAinv = 1.0/DELA;
                    Apert = A0 + DELA;
                    here->CNTl = Apert;
                }
                else{
                    A0 = here->CNTw;
                    DELA = info->SENpertfac * A0;
                    DELAinv = 1.0/DELA;
                    Apert = A0 + DELA;
                    here->CNTw = Apert;
                    here->CNTdrainArea *= (1 + info->SENpertfac);
                    here->CNTsourceArea *= (1 + info->SENpertfac);
                    here->CNTCbd *= (1 + info->SENpertfac);
                    here->CNTCbs *= (1 + info->SENpertfac);
                    if(here->CNTdrainPerimiter){
                        here->CNTCbdsw += here->CNTCbdsw *
                            DELA/here->CNTdrainPerimiter;
                    }
                    if(here->CNTsourcePerimiter){
                        here->CNTCbssw += here->CNTCbssw *
                            DELA/here->CNTsourcePerimiter;
                    }
                    if(*(ckt->CKTstate0 + here->CNTvbd) >=
                        here->CNTtDepCap){
                        arg = 1-model->CNTfwdCapDepCoeff;
                        sarg = exp( (-model->CNTbulkJctBotGradingCoeff) * 
                                log(arg) );
                        sargsw = exp( (-model->CNTbulkJctSideGradingCoeff) * 
                                log(arg) );
                        here->CNTf2d = here->CNTCbd*
                                (1-model->CNTfwdCapDepCoeff*
                                (1+model->CNTbulkJctBotGradingCoeff))* sarg/arg
                                +  here->CNTCbdsw*(1-model->CNTfwdCapDepCoeff*
                                (1+model->CNTbulkJctSideGradingCoeff))*
                                sargsw/arg;
                        here->CNTf3d = here->CNTCbd * 
                                model->CNTbulkJctBotGradingCoeff * sarg/arg/
                                here->CNTtBulkPot
                                + here->CNTCbdsw * 
                                model->CNTbulkJctSideGradingCoeff * sargsw/arg/
                                here->CNTtBulkPot;
                        here->CNTf4d = here->CNTCbd*
                                here->CNTtBulkPot*(1-arg*sarg)/
                                (1-model->CNTbulkJctBotGradingCoeff)
                                + here->CNTCbdsw*here->CNTtBulkPot*
                                (1-arg*sargsw)/
                                (1-model->CNTbulkJctSideGradingCoeff)
                                -here->CNTf3d/2*
                                (here->CNTtDepCap*here->CNTtDepCap)
                                -here->CNTtDepCap * here->CNTf2d;
                    }
                    if(*(ckt->CKTstate0 + here->CNTvbs) >=
                        here->CNTtDepCap){
                        arg = 1-model->CNTfwdCapDepCoeff;
                        sarg = exp( (-model->CNTbulkJctBotGradingCoeff) * 
                                log(arg) );
                        sargsw = exp( (-model->CNTbulkJctSideGradingCoeff) * 
                                log(arg) );
                        here->CNTf2s = here->CNTCbs*
                                (1-model->CNTfwdCapDepCoeff*
                                (1+model->CNTbulkJctBotGradingCoeff))* sarg/arg
                                +  here->CNTCbssw*(1-model->CNTfwdCapDepCoeff*
                                (1+model->CNTbulkJctSideGradingCoeff))*
                                sargsw/arg;
                        here->CNTf3s = here->CNTCbs * 
                                model->CNTbulkJctBotGradingCoeff * sarg/arg/
                                here->CNTtBulkPot
                                + here->CNTCbssw * 
                                model->CNTbulkJctSideGradingCoeff * sargsw/arg/
                                here->CNTtBulkPot;
                        here->CNTf4s = here->CNTCbs*
                                here->CNTtBulkPot*(1-arg*sarg)/
                                (1-model->CNTbulkJctBotGradingCoeff)
                                + here->CNTCbssw*here->CNTtBulkPot*
                                (1-arg*sargsw)/
                                (1-model->CNTbulkJctSideGradingCoeff)
                                -here->CNTf3s/2*
                                (here->CNTtDepCap*here->CNTtDepCap)
                                -here->CNTtDepCap * here->CNTf2s;
                    }
                    here->CNTdrainConductance *= Apert/A0;
                    here->CNTsourceConductance *= Apert/A0;
                }


#ifdef SENSDEBUG
                if(flag == 0)
                    printf("perturbation of l\n");
                if(flag == 1)
                    printf("perturbation of w\n");
#endif /* SENSDEBUG */

                error =  CNTload((GENmodel*)model,ckt);
                if(error) return(error);

                if(flag == 0){
                    here->CNTl = A0;
                }
                else{
                    here->CNTw = A0;
                    here->CNTdrainArea /= (1 + info->SENpertfac);
                    here->CNTsourceArea /= (1 + info->SENpertfac);
                    here->CNTdrainConductance *= A0/Apert;
                    here->CNTsourceConductance *= A0/Apert;
                }
                cd =  here->CNTcd ;
                cbd = here->CNTcbd ;
                cbs = here->CNTcbs ;

                gspr= here->CNTsourceConductance ;
                gdpr= here->CNTdrainConductance ;

                DcdDp = (cd - cd0) * DELAinv;
                DcbsDp = (cbs - cbs0) * DELAinv;
                DcbdDp = (cbd - cbd0) * DELAinv;
                DcbDp = ( DcbsDp + DcbdDp );

                DcdprDp = 0;
                DcsprDp = 0;
                if(here->CNTdNode != here->CNTdNodePrime)
                    if(gdpr0) DcdprDp = cdpr0 * (gdpr - gdpr0)/gdpr0 * DELAinv;
                if(here->CNTsNode != here->CNTsNodePrime)
                    if(gspr0) DcsprDp = cspr0 * (gspr - gspr0)/gspr0 * DELAinv;

                DcdprmDp = ( - DcdprDp + DcdDp);
                DcsprmDp = ( - DcbsDp - DcdDp - DcbdDp  - DcsprDp);

                if(flag == 0){
                    EffectiveLength = here->CNTl 
                        - 2*model->CNTlatDiff;
                    if(EffectiveLength == 0){ 
                        DqgsDp = 0;
                        DqgdDp = 0;
                        DqgbDp = 0;
                    }
                    else{
                        DqgsDp = model->CNTtype * qgs0 / EffectiveLength;
                        DqgdDp = model->CNTtype * qgd0 / EffectiveLength;
                        DqgbDp = model->CNTtype * qgb0 / EffectiveLength;
                    }
                }
                else{
                    DqgsDp = model->CNTtype * qgs0 / here->CNTw;
                    DqgdDp = model->CNTtype * qgd0 / here->CNTw;
                    DqgbDp = model->CNTtype * qgb0 / here->CNTw;
                }


                qbd = *(ckt->CKTstate0 + here->CNTqbd);
                qbs = *(ckt->CKTstate0 + here->CNTqbs);

                DqbsDp = model->CNTtype * (qbs - qbs0)*DELAinv;
                DqbdDp = model->CNTtype * (qbd - qbd0)*DELAinv;

                if(flag == 0){
                    *(here->CNTdphigs_dl) = DqgsDp;
                    *(here->CNTdphigd_dl) = DqgdDp;
                    *(here->CNTdphibs_dl) = DqbsDp;
                    *(here->CNTdphibd_dl) = DqbdDp;
                    *(here->CNTdphigb_dl) = DqgbDp;
                }
                else{
                    *(here->CNTdphigs_dw) = DqgsDp;
                    *(here->CNTdphigd_dw) = DqgdDp;
                    *(here->CNTdphibs_dw) = DqbsDp;
                    *(here->CNTdphibd_dw) = DqbdDp;
                    *(here->CNTdphigb_dw) = DqgbDp;
                }


#ifdef SENSDEBUG
                printf("CKTag[0]=%.7e,CKTag[1]=%.7e,flag= %d\n",
                        ckt->CKTag[0],ckt->CKTag[1],flag);
                printf("cd0 = %.7e ,cd = %.7e,\n",cd0,cd); 
                printf("cbs0 = %.7e ,cbs = %.7e,\n",cbs0,cbs); 
                printf("cbd0 = %.7e ,cbd = %.7e,\n",cbd0,cbd); 
                printf("DcdprmDp = %.7e,\n",DcdprmDp); 
                printf("DcsprmDp = %.7e,\n",DcsprmDp); 
                printf("DcdprDp = %.7e,\n",DcdprDp); 
                printf("DcsprDp = %.7e,\n",DcsprDp); 
                printf("qgs0 = %.7e \n",qgs0); 
                printf("qgd0 = %.7e \n",qgd0); 
                printf("qgb0 = %.7e \n",qgb0); 
                printf("qbs0 = %.7e ,qbs = %.7e,\n",qbs0,qbs); 
                printf("qbd0 = %.7e ,qbd = %.7e,\n",qbd0,qbd); 
                printf("DqgsDp = %.7e \n",DqgsDp); 
                printf("DqgdDp = %.7e \n",DqgdDp); 
                printf("DqgbDp = %.7e \n",DqgbDp); 
                printf("DqbsDp = %.7e \n",DqbsDp); 
                printf("DqbdDp = %.7e \n",DqbdDp); 
                printf("EffectiveLength = %.7e \n",EffectiveLength); 
                printf("tdepCap = %.7e \n",here->CNTtDepCap); 
                printf("\n");
#endif /* SENSDEBUG*/
                if((info->SENmode == TRANSEN) &&
                    (ckt->CKTmode & MODEINITTRAN))
                    goto next2;

                /*
                                 *   load RHS matrix
                                 */

                if(flag == 0){
                    *(info->SEN_RHS[here->CNTbNode] + here->CNTsenParmNo) -= 
                            model->CNTtype * DcbDp;
                    *(info->SEN_RHS[here->CNTdNode] + here->CNTsenParmNo) -= 
                            model->CNTtype * DcdprDp;
                    *(info->SEN_RHS[here->CNTdNodePrime] +
                            here->CNTsenParmNo) -= model->CNTtype * DcdprmDp;
                    *(info->SEN_RHS[here->CNTsNode] + here->CNTsenParmNo) -= 
                            model->CNTtype * DcsprDp;
                    *(info->SEN_RHS[here->CNTsNodePrime] + 
                            here->CNTsenParmNo) -= model->CNTtype * DcsprmDp;
                }
                else{  
                    offset = here->CNTsens_l;

                    *(info->SEN_RHS[here->CNTbNode] + here->CNTsenParmNo + 
                            offset) -= model->CNTtype * DcbDp;
                    *(info->SEN_RHS[here->CNTdNode] + here->CNTsenParmNo + 
                            offset) -= model->CNTtype * DcdprDp;
                    *(info->SEN_RHS[here->CNTdNodePrime] + here->CNTsenParmNo
                            + offset) -= model->CNTtype * DcdprmDp;
                    *(info->SEN_RHS[here->CNTsNode] + here->CNTsenParmNo +
                            offset) -= model->CNTtype * DcsprDp;
                    *(info->SEN_RHS[here->CNTsNodePrime] + here->CNTsenParmNo
                            + offset) -= model->CNTtype * DcsprmDp;
                }
#ifdef SENSDEBUG
                printf("after loading\n");
                if(flag == 0){
                    printf("DcbDp=%.7e\n",
                            *(info->SEN_RHS[here->CNTbNode] +
                            here->CNTsenParmNo));
                    printf("DcdprDp=%.7e\n",
                            *(info->SEN_RHS[here->CNTdNode] +
                            here->CNTsenParmNo));
                    printf("DcsprDp=%.7e\n",
                            *(info->SEN_RHS[here->CNTsNode] +
                            here->CNTsenParmNo));
                    printf("DcdprmDp=%.7e\n",
                            *(info->SEN_RHS[here->CNTdNodePrime] +
                            here->CNTsenParmNo));
                    printf("DcsprmDp=%.7e\n",
                            *(info->SEN_RHS[here->CNTsNodePrime] +
                            here->CNTsenParmNo));
                    printf("\n");
                }
                else{
                    printf("DcbDp=%.7e\n",
                            *(info->SEN_RHS[here->CNTbNode] + 
                            here->CNTsenParmNo + here->CNTsens_l));
                    printf("DcdprDp=%.7e\n",
                            *(info->SEN_RHS[here->CNTdNode] + 
                            here->CNTsenParmNo + here->CNTsens_l));
                    printf("DcsprDp=%.7e\n",
                            *(info->SEN_RHS[here->CNTsNode] + 
                            here->CNTsenParmNo + here->CNTsens_l));
                    printf("DcdprmDp=%.7e\n",
                            *(info->SEN_RHS[here->CNTdNodePrime] + 
                            here->CNTsenParmNo + here->CNTsens_l));
                    printf("DcsprmDp=%.7e\n",
                            *(info->SEN_RHS[here->CNTsNodePrime] + 
                            here->CNTsenParmNo + here->CNTsens_l));
                }
#endif /* SENSDEBUG*/
next2:                   
                ;
            }
next1:
            if((info->SENmode == DCSEN) ||
                (ckt->CKTmode&MODETRANOP) ) goto restore;
            if((info->SENmode == TRANSEN) &&
                (ckt->CKTmode & MODEINITTRAN)) goto restore;
            for(iparmno = 1;iparmno<=info->SENparms;iparmno++){
#ifdef SENSDEBUG
                printf("after conductive currents\n");
                printf("iparmno = %d\n",iparmno);
                printf("DcbDp=%.7e\n",
                        *(info->SEN_RHS[here->CNTbNode] + iparmno));
                printf("DcdprDp=%.7e\n",
                        *(info->SEN_RHS[here->CNTdNode] + iparmno));
                printf("DcdprmDp=%.7e\n",
                        *(info->SEN_RHS[here->CNTdNodePrime] + iparmno));
                printf("DcsprDp=%.7e\n",
                        *(info->SEN_RHS[here->CNTsNode] + iparmno));
                printf("DcsprmDp=%.7e\n",
                        *(info->SEN_RHS[here->CNTsNodePrime] + iparmno));
                printf("\n");
#endif /* SENSDEBUG */
                Osxpgs = tag0 * *(ckt->CKTstate1 + here->CNTsensxpgs +
                        10*(iparmno - 1))
                        + tag1 * *(ckt->CKTstate1 + here->CNTsensxpgs +
                        10*(iparmno - 1) + 1);

                Osxpgd = tag0 * *(ckt->CKTstate1 + here->CNTsensxpgd +
                        10*(iparmno - 1))
                        + tag1 * *(ckt->CKTstate1 + here->CNTsensxpgd +
                        10*(iparmno - 1) + 1);

                Osxpbs = tag0 * *(ckt->CKTstate1 + here->CNTsensxpbs +
                        10*(iparmno - 1))
                        + tag1 * *(ckt->CKTstate1 + here->CNTsensxpbs +
                        10*(iparmno - 1) + 1);

                Osxpbd =tag0 * *(ckt->CKTstate1 + here->CNTsensxpbd +
                        10*(iparmno - 1))
                        + tag1 * *(ckt->CKTstate1 + here->CNTsensxpbd +
                        10*(iparmno - 1) + 1);
                Osxpgb = tag0 * *(ckt->CKTstate1 + here->CNTsensxpgb +
                        10*(iparmno - 1))
                        + tag1 * *(ckt->CKTstate1 + here->CNTsensxpgb +
                        10*(iparmno - 1) + 1);

#ifdef SENSDEBUG
                printf("iparmno=%d\n",iparmno);
                printf("sxpgs=%.7e,sdgs=%.7e\n",
                        *(ckt->CKTstate1 + here->CNTsensxpgs + 
                        10*(iparmno - 1)), *(ckt->CKTstate1 + 
                        here->CNTsensxpgs + 10*(iparmno - 1) + 1));
                printf("sxpgd=%.7e,sdgd=%.7e\n",
                        *(ckt->CKTstate1 + here->CNTsensxpgd + 
                        10*(iparmno - 1)), *(ckt->CKTstate1 + 
                        here->CNTsensxpgd + 10*(iparmno - 1) + 1));
                printf("sxpbs=%.7e,sdbs=%.7e\n",
                        *(ckt->CKTstate1 + here->CNTsensxpbs + 
                        10*(iparmno - 1)), *(ckt->CKTstate1 + 
                        here->CNTsensxpbs + 10*(iparmno - 1) + 1));
                printf("sxpbd=%.7e,sdbd=%.7e\n",
                        *(ckt->CKTstate1 + here->CNTsensxpbd + 
                        10*(iparmno - 1)), *(ckt->CKTstate1 + 
                        here->CNTsensxpbd + 10*(iparmno - 1) + 1));
                printf("sxpgb=%.7e,sdgb=%.7e\n",
                        *(ckt->CKTstate1 + here->CNTsensxpgb + 
                        10*(iparmno - 1)), *(ckt->CKTstate1 + 
                        here->CNTsensxpgb + 10*(iparmno - 1) + 1));
                printf("before loading DqDp\n");
                printf("Osxpgs=%.7e,Osxpgd=%.7e\n",Osxpgs,Osxpgd);
                printf("Osxpbs=%.7e,Osxpbd=%.7e,Osxpgb=%.7e\n",
                        Osxpbs,Osxpbd,Osxpgb);
                printf("\n");
#endif /* SENSDEBUG */
                if(here->CNTsens_l && (iparmno == here->CNTsenParmNo)){
                    Osxpgs -= tag0 * *(here->CNTdphigs_dl);
                    Osxpgd -= tag0 * *(here->CNTdphigd_dl);
                    Osxpbs -= tag0 * *(here->CNTdphibs_dl);
                    Osxpbd -= tag0 * *(here->CNTdphibd_dl);
                    Osxpgb -= tag0 * *(here->CNTdphigb_dl);
                }
                if(here->CNTsens_w && 
                        (iparmno == (here->CNTsenParmNo + here->CNTsens_l))){
                    Osxpgs -= tag0 * *(here->CNTdphigs_dw);
                    Osxpgd -= tag0 * *(here->CNTdphigd_dw);
                    Osxpbs -= tag0 * *(here->CNTdphibs_dw);
                    Osxpbd -= tag0 * *(here->CNTdphibd_dw);
                    Osxpgb -= tag0 * *(here->CNTdphigb_dw);
                }
#ifdef SENSDEBUG
                printf("after loading DqDp\n");
                printf("DqgsDp=%.7e",DqgsDp);
                printf("Osxpgs=%.7e,Osxpgd=%.7e\n",Osxpgs,Osxpgd);
                printf("Osxpbs=%.7e,Osxpbd=%.7e,Osxpgb=%.7e\n",
                        Osxpbs,Osxpbd,Osxpgb);
#endif /* SENSDEBUG */

                *(info->SEN_RHS[here->CNTbNode] + iparmno) += 
                        Osxpbs + Osxpbd -Osxpgb;
                *(info->SEN_RHS[here->CNTgNode] + iparmno) += 
                        Osxpgs + Osxpgd + Osxpgb;

                *(info->SEN_RHS[here->CNTdNodePrime] + iparmno) -= 
                        Osxpgd + Osxpbd ;
                *(info->SEN_RHS[here->CNTsNodePrime] + iparmno) -= 
                        Osxpgs + Osxpbs;
#ifdef SENSDEBUG
                printf("after capacitive currents\n");
                printf("DcbDp=%.7e\n",
                        *(info->SEN_RHS[here->CNTbNode] + iparmno));
                printf("DcdprDp=%.7e\n",
                        *(info->SEN_RHS[here->CNTdNode] + iparmno));
                printf("DcdprmDp=%.7e\n",
                        *(info->SEN_RHS[here->CNTdNodePrime] + iparmno));
                printf("DcsprDp=%.7e\n",
                        *(info->SEN_RHS[here->CNTsNode] + iparmno));
                printf("DcsprmDp=%.7e\n",
                        *(info->SEN_RHS[here->CNTsNodePrime] + iparmno));
#endif /* SENSDEBUG */

            }
restore:    /* put the unperturbed values back into the state vector */
            for(i=0; i <= 16; i++)
                *(ckt->CKTstate0 + here->CNTstates + i) = *(SaveState + i);
            here->CNTsourceConductance = *(SaveState + 17) ;   
            here->CNTdrainConductance = *(SaveState + 18) ; 
            here->CNTcd =  *(SaveState + 19) ;  
            here->CNTcbs =  *(SaveState + 20) ;  
            here->CNTcbd =  *(SaveState + 21) ;  
            here->CNTgmbs =  *(SaveState + 22) ;  
            here->CNTgm =  *(SaveState + 23) ;  
            here->CNTgds =  *(SaveState + 24) ;  
            here->CNTgbd =  *(SaveState + 25) ;  
            here->CNTgbs =  *(SaveState + 26) ;  
            here->CNTcapbd =  *(SaveState + 27) ;  
            here->CNTcapbs =  *(SaveState + 28) ;  
            here->CNTCbd =  *(SaveState + 29) ;  
            here->CNTCbdsw =  *(SaveState + 30) ;  
            here->CNTCbs =  *(SaveState + 31) ;  
            here->CNTCbssw =  *(SaveState + 32) ;  
            here->CNTf2d =  *(SaveState + 33) ;  
            here->CNTf3d =  *(SaveState + 34) ;  
            here->CNTf4d =  *(SaveState + 35) ;  
            here->CNTf2s =  *(SaveState + 36) ;  
            here->CNTf3s =  *(SaveState + 37) ;  
            here->CNTf4s =  *(SaveState + 38) ;  
            here->CNTcgs = *(SaveState + 39) ;  
            here->CNTcgd = *(SaveState + 40) ;  
            here->CNTcgb = *(SaveState + 41) ;  
            here->CNTvdsat = *(SaveState + 42) ;  
            here->CNTvon = *(SaveState + 43) ;   
            here->CNTmode = save_mode ;  

            here->CNTsenPertFlag = OFF;

        }
    }
    info->SENstatus = NORMAL;
#ifdef SENSDEBUG
    printf("CNTsenload end\n");
#endif /* SENSDEBUG */
    return(OK);
}

