/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles

This function is obsolete (was used by an old sensitivity analysis)
**********/

#include "ngspice/ngspice.h"
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/const.h"
#include "cntdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

/* actually load the current ac sensitivity 
 * information into the  array previously provided 
 */

int
CNTsAcLoad(GENmodel *inModel, CKTcircuit *ckt)
{
    CNTmodel *model = (CNTmodel*)inModel;
    CNTinstance *here;
    int xnrm;
    int xrev;
    double A0;
    double Apert = 0.0;
    double DELA;
    double DELAinv;
    double gdpr0;
    double gspr0;
    double gds0;
    double gbs0;
    double gbd0;
    double gm0;
    double gmbs0;
    double gdpr;
    double gspr;
    double gds;
    double gbs;
    double gbd;
    double gm;
    double gmbs;
    double xcgs0;
    double xcgd0;
    double xcgb0;
    double xbd0;
    double xbs0;
    double xcgs;
    double xcgd;
    double xcgb;
    double xbd;
    double xbs;
    double vbsOp;
    double vbdOp;
    double vspr;
    double vdpr;
    double vgs;
    double vgd;
    double vgb;
    double vbs;
    double vbd;
    double vds;
    double ivspr;
    double ivdpr;
    double ivgs;
    double ivgd;
    double ivgb;
    double ivbs;
    double ivbd;
    double ivds;
    double cspr;
    double cdpr;
    double cgs;
    double cgd;
    double cgb;
    double cbs;
    double cbd;
    double cds;
    double cs0;
    double csprm0;
    double cd0;
    double cdprm0;
    double cg0;
    double cb0;
    double cs;
    double csprm;
    double cd;
    double cdprm;
    double cg;
    double cb;
    double icspr;
    double icdpr;
    double icgs;
    double icgd;
    double icgb;
    double icbs;
    double icbd;
    double icds;
    double ics0;
    double icsprm0;
    double icd0;
    double icdprm0;
    double icg0;
    double icb0;
    double ics;
    double icsprm;
    double icd;
    double icdprm;
    double icg;
    double icb;
    double DvDp = 0.0;
    int i;
    int flag;
    int error;
    int iparmno;
    double arg;
    double sarg;
    double sargsw;
    double SaveState[44];
    int    save_mode;
    SENstruct *info;

#ifdef SENSDEBUG
    printf("CNTsenacload\n");
    printf("CKTomega = %.5e\n",ckt->CKTomega);
#endif /* SENSDEBUG */
    info = ckt->CKTsenInfo;
    info->SENstatus = PERTURBATION;
    for( ; model != NULL; model = model->CNTnextModel) {
        for(here = model->CNTinstances; here!= NULL;
                here = here->CNTnextInstance) {
	    /*if (here->CNTowner != ARCHme) continue;*/ /*cuidado*/

            /* save the unperturbed values in the state vector */
            for(i=0; i <= 16; i++)
                *(SaveState + i) = *(ckt->CKTstate0 + here->CNTstates + i);

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

            xnrm=1;
            xrev=0;
            if (here->CNTmode < 0) {
                xnrm=0;
                xrev=1;
            }

            vbsOp = model->CNTtype * ( 
            *(ckt->CKTrhsOp+here->CNTbNode) -
                *(ckt->CKTrhsOp+here->CNTsNodePrime));
            vbdOp = model->CNTtype * ( 
            *(ckt->CKTrhsOp+here->CNTbNode) -
                *(ckt->CKTrhsOp+here->CNTdNodePrime));
            vspr = *(ckt->CKTrhsOld + here->CNTsNode) 
                - *(ckt->CKTrhsOld +
                    here->CNTsNodePrime) ;
            ivspr = *(ckt->CKTirhsOld + here->CNTsNode) 
                - *(ckt->CKTirhsOld +
                    here->CNTsNodePrime) ;
            vdpr = *(ckt->CKTrhsOld + here->CNTdNode) 
                - *(ckt->CKTrhsOld +
                    here->CNTdNodePrime) ;
            ivdpr = *(ckt->CKTirhsOld + here->CNTdNode) 
                - *(ckt->CKTirhsOld +
                    here->CNTdNodePrime) ;
            vgb = *(ckt->CKTrhsOld + here->CNTgNode) 
                - *(ckt->CKTrhsOld +
                    here->CNTbNode) ;
            ivgb = *(ckt->CKTirhsOld + here->CNTgNode) 
                - *(ckt->CKTirhsOld +
                    here->CNTbNode) ;
            vbs = *(ckt->CKTrhsOld + here->CNTbNode) 
                - *(ckt->CKTrhsOld +
                    here->CNTsNodePrime) ;
            ivbs = *(ckt->CKTirhsOld + here->CNTbNode) 
                - *(ckt->CKTirhsOld +
                    here->CNTsNodePrime) ;
            vbd = *(ckt->CKTrhsOld + here->CNTbNode) 
                - *(ckt->CKTrhsOld +
                    here->CNTdNodePrime) ;
            ivbd = *(ckt->CKTirhsOld + here->CNTbNode) 
                - *(ckt->CKTirhsOld +
                    here->CNTdNodePrime) ;
            vds = vbs - vbd ;
            ivds = ivbs - ivbd ;
            vgs = vgb + vbs ;
            ivgs = ivgb + ivbs ;
            vgd = vgb + vbd ;
            ivgd = ivgb + ivbd ;

#ifdef SENSDEBUG
            printf("senacload instance name %s\n",here->CNTname);
            printf("gate = %d ,drain = %d, drainprm = %d\n", 
                    here->CNTgNode,here->CNTdNode,here->CNTdNodePrime);
            printf("source = %d , sourceprm = %d ,body = %d, senparmno = %d\n",
                    here->CNTsNode ,here->CNTsNodePrime,here->CNTbNode,
                    here->CNTsenParmNo);
            printf("\n without  perturbation \n");
#endif /* SENSDEBUG */
            /*  without  perturbation  */
            *(ckt->CKTstate0 + here->CNTvbs) = vbsOp;
            *(ckt->CKTstate0 + here->CNTvbd) = vbdOp;

            here->CNTsenPertFlag = ON ;
            if(info->SENacpertflag == 1){
                /* store the  unperturbed values of small signal parameters */
                if((error = CNTload((GENmodel*)model,ckt)))
		  return(error);
                *(here->CNTsenCgs) = here->CNTcgs;
                *(here->CNTsenCgd) = here->CNTcgd;
                *(here->CNTsenCgb) = here->CNTcgb;
                *(here->CNTsenCbd) = here->CNTcapbd;
                *(here->CNTsenCbs) = here->CNTcapbs;
                *(here->CNTsenGds) = here->CNTgds;
                *(here->CNTsenGbs) = here->CNTgbs;
                *(here->CNTsenGbd) = here->CNTgbd;
                *(here->CNTsenGm) = here->CNTgm;
                *(here->CNTsenGmbs) = here->CNTgmbs;

            }
            xcgs0= *(here->CNTsenCgs) * ckt->CKTomega;
            xcgd0= *(here->CNTsenCgd) * ckt->CKTomega;
            xcgb0= *(here->CNTsenCgb) * ckt->CKTomega;
            xbd0= *(here->CNTsenCbd) * ckt->CKTomega;
            xbs0= *(here->CNTsenCbs) * ckt->CKTomega;
            gds0= *(here->CNTsenGds);
            gbs0= *(here->CNTsenGbs);
            gbd0= *(here->CNTsenGbd);
            gm0= *(here->CNTsenGm);
            gmbs0= *(here->CNTsenGmbs);
            gdpr0 = here->CNTdrainConductance;
            gspr0 = here->CNTsourceConductance;


            cspr = gspr0 * vspr ;
            icspr = gspr0 * ivspr ;
            cdpr = gdpr0 * vdpr ;
            icdpr = gdpr0 * ivdpr ;
            cgs = ( - xcgs0 * ivgs );
            icgs =  xcgs0 * vgs ;
            cgd = ( - xcgd0 * ivgd );
            icgd =  xcgd0 * vgd ;
            cgb = ( - xcgb0 * ivgb );
            icgb =  xcgb0 * vgb ;
            cbs = ( gbs0 * vbs  - xbs0 * ivbs );
            icbs = ( xbs0 * vbs  + gbs0 * ivbs );
            cbd = ( gbd0 * vbd  - xbd0 * ivbd );
            icbd = ( xbd0 * vbd  + gbd0 * ivbd );
            cds = ( gds0 * vds  + xnrm * (gm0 * vgs + gmbs0 * vbs) 
                - xrev * (gm0 * vgd + gmbs0 * vbd) );
            icds = ( gds0 * ivds + xnrm * (gm0 * ivgs + gmbs0 * ivbs)
                - xrev * (gm0 * ivgd + gmbs0 * ivbd) );

            cs0 = cspr;
            ics0 = icspr;
            csprm0 = ( -cspr - cgs - cbs - cds ) ;
            icsprm0 = ( -icspr - icgs - icbs - icds ) ;
            cd0 = cdpr;
            icd0 = icdpr;
            cdprm0 = ( -cdpr - cgd - cbd + cds ) ;
            icdprm0 = ( -icdpr - icgd - icbd + icds ) ;
            cg0 = cgs + cgd + cgb ;
            icg0 = icgs + icgd + icgb ;
            cb0 = cbs + cbd - cgb ;
            icb0 = icbs + icbd - icgb ;
#ifdef SENSDEBUG
            printf("gspr0 = %.7e , gdpr0 = %.7e , gds0 = %.7e, gbs0 = %.7e\n",
                    gspr0,gdpr0,gds0,gbs0);
            printf("gbd0 = %.7e , gm0 = %.7e , gmbs0 = %.7e\n",gbd0,gm0,gmbs0);
            printf("xcgs0 = %.7e , xcgd0 = %.7e , xcgb0 = %.7e,",
                    xcgs0,xcgd0,xcgb0);
            printf("xbd0 = %.7e,xbs0 = %.7e\n",xbd0,xbs0);
            printf("vbs = %.7e , vbd = %.7e , vgb = %.7e\n",vbs,vbd,vgb);
            printf("ivbs = %.7e , ivbd = %.7e , ivgb = %.7e\n",ivbs,ivbd,ivgb);
            printf("cbs0 = %.7e , cbd0 = %.7e , cgb0 = %.7e\n",cbs,cbd,cgb);
            printf("cb0 = %.7e , cg0 = %.7e , cs0 = %.7e\n",cb0,cg0,cs0);
            printf("csprm0 = %.7e, cd0 = %.7e, cdprm0 = %.7e\n",
                    csprm0,cd0,cdprm0);
            printf("icb0 = %.7e , icg0 = %.7e , ics0 = %.7e\n",icb0,icg0,ics0);
            printf("icsprm0 = %.7e, icd0 = %.7e, icdprm0 = %.7e\n",
                    icsprm0,icd0,icdprm0);
            printf("\nPerturbation of vbs\n");
#endif /* SENSDEBUG */

            /* Perturbation of vbs */
            flag = 1;
            A0 = vbsOp;
            DELA =  info->SENpertfac * here->CNTtVto ;
            DELAinv = 1.0/DELA;

            if(info->SENacpertflag == 1){
                /* store the  values of small signal parameters 
                 * corresponding to perturbed vbs */
                Apert = A0 + DELA;
                *(ckt->CKTstate0 + here->CNTvbs) = Apert;
                *(ckt->CKTstate0 + here->CNTvbd) = vbdOp;

                if((error = CNTload((GENmodel*)model,ckt)))
		  return(error);

                *(here->CNTsenCgs + 1) = here->CNTcgs;
                *(here->CNTsenCgd + 1) = here->CNTcgd;
                *(here->CNTsenCgb + 1) = here->CNTcgb;
                *(here->CNTsenCbd + 1) = here->CNTcapbd;
                *(here->CNTsenCbs + 1) = here->CNTcapbs;
                *(here->CNTsenGds + 1) = here->CNTgds;
                *(here->CNTsenGbs + 1) = here->CNTgbs;
                *(here->CNTsenGbd + 1) = here->CNTgbd;
                *(here->CNTsenGm + 1) = here->CNTgm;
                *(here->CNTsenGmbs + 1) = here->CNTgmbs;

                *(ckt->CKTstate0 + here->CNTvbs) = A0;


            }

            goto load;


pertvbd:  /* Perturbation of vbd */
#ifdef SENSDEBUG
            printf("\nPerturbation of vbd\n");
#endif /* SENSDEBUG */

            flag = 2;
            A0 = vbdOp;
            DELA =  info->SENpertfac * here->CNTtVto + 1e-8;
            DELAinv = 1.0/DELA;

            if(info->SENacpertflag == 1){
                /* store the  values of small signal parameters 
                 * corresponding to perturbed vbd */
                Apert = A0 + DELA;
                *(ckt->CKTstate0 + here->CNTvbs) = vbsOp;
                *(ckt->CKTstate0 + here->CNTvbd) = Apert;

                if((error = CNTload((GENmodel*)model,ckt)))
		  return(error);

                *(here->CNTsenCgs + 2) = here->CNTcgs;
                *(here->CNTsenCgd + 2) = here->CNTcgd;
                *(here->CNTsenCgb + 2) = here->CNTcgb;
                *(here->CNTsenCbd + 2) = here->CNTcapbd;
                *(here->CNTsenCbs + 2) = here->CNTcapbs;
                *(here->CNTsenGds + 2) = here->CNTgds;
                *(here->CNTsenGbs + 2) = here->CNTgbs;
                *(here->CNTsenGbd + 2) = here->CNTgbd;
                *(here->CNTsenGm + 2) = here->CNTgm;
                *(here->CNTsenGmbs + 2) = here->CNTgmbs;

                *(ckt->CKTstate0 + here->CNTvbd) = A0;

            }

            goto load;


pertvgb:  /* Perturbation of vgb */
#ifdef SENSDEBUG
            printf("\nPerturbation of vgb\n");
#endif /* SENSDEBUG */

            flag = 3;
            A0 = model->CNTtype * (*(ckt->CKTrhsOp + here->CNTgNode) 
                -  *(ckt->CKTrhsOp + here->CNTbNode)); 
            DELA =  info->SENpertfac * A0 + 1e-8;
            DELAinv = model->CNTtype * 1.0/DELA;


            if(info->SENacpertflag == 1){
                /* store the  values of small signal parameters 
                 * corresponding to perturbed vgb */
                *(ckt->CKTstate0 + here->CNTvbs) = vbsOp;
                *(ckt->CKTstate0 + here->CNTvbd) = vbdOp;
                *(ckt->CKTrhsOp + here->CNTbNode) -= DELA; 

                if((error = CNTload((GENmodel*)model,ckt)))
		  return(error);

                *(here->CNTsenCgs + 3) = here->CNTcgs;
                *(here->CNTsenCgd + 3) = here->CNTcgd;
                *(here->CNTsenCgb + 3) = here->CNTcgb;
                *(here->CNTsenCbd + 3) = here->CNTcapbd;
                *(here->CNTsenCbs + 3) = here->CNTcapbs;
                *(here->CNTsenGds + 3) = here->CNTgds;
                *(here->CNTsenGbs + 3) = here->CNTgbs;
                *(here->CNTsenGbd + 3) = here->CNTgbd;
                *(here->CNTsenGm + 3) = here->CNTgm;
                *(here->CNTsenGmbs + 3) = here->CNTgmbs;


                *(ckt->CKTrhsOp + here->CNTbNode) += DELA; 
            }
            goto load;

pertl:    /* Perturbation of length */

            if(here->CNTsens_l == 0){
                goto pertw;
            }
#ifdef SENSDEBUG
            printf("\nPerturbation of length\n");
#endif /* SENSDEBUG */
            flag = 4;
            A0 = here->CNTl;
            DELA =  info->SENpertfac * A0;
            DELAinv = 1.0/DELA;

            if(info->SENacpertflag == 1){
                /* store the  values of small signal parameters 
                 * corresponding to perturbed length */
                Apert = A0 + DELA;
                here->CNTl = Apert;

                *(ckt->CKTstate0 + here->CNTvbs) = vbsOp;
                *(ckt->CKTstate0 + here->CNTvbd) = vbdOp;

                if ((error = CNTload((GENmodel*)model,ckt)))
		  return(error);

                *(here->CNTsenCgs + 4) = here->CNTcgs;
                *(here->CNTsenCgd + 4) = here->CNTcgd;
                *(here->CNTsenCgb + 4) = here->CNTcgb;
                *(here->CNTsenCbd + 4) = here->CNTcapbd;
                *(here->CNTsenCbs + 4) = here->CNTcapbs;
                *(here->CNTsenGds + 4) = here->CNTgds;
                *(here->CNTsenGbs + 4) = here->CNTgbs;
                *(here->CNTsenGbd + 4) = here->CNTgbd;
                *(here->CNTsenGm + 4) = here->CNTgm;
                *(here->CNTsenGmbs + 4) = here->CNTgmbs;

                here->CNTl = A0;

            }

            goto load;

pertw:    /* Perturbation of width */
            if(here->CNTsens_w == 0)
                goto next;
#ifdef SENSDEBUG
            printf("\nPerturbation of width\n");
#endif /* SENSDEBUG */
            flag = 5;
            A0 = here->CNTw;
            DELA = info->SENpertfac * A0;
            DELAinv = 1.0/DELA;
            Apert = A0 + DELA;

            if(info->SENacpertflag == 1){
                /* store the  values of small signal parameters 
                 * corresponding to perturbed width */
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
                if(vbdOp >= here->CNTtDepCap){
                    arg = 1-model->CNTfwdCapDepCoeff;
                    sarg = exp( (-model->CNTbulkJctBotGradingCoeff) * 
                            log(arg) );
                    sargsw = exp( (-model->CNTbulkJctSideGradingCoeff) * 
                            log(arg) );
                    here->CNTf2d = here->CNTCbd*(1-model->CNTfwdCapDepCoeff*
                        (1+model->CNTbulkJctBotGradingCoeff))* sarg/arg
                        +  here->CNTCbdsw*(1-model->CNTfwdCapDepCoeff*
                        (1+model->CNTbulkJctSideGradingCoeff))*
                        sargsw/arg;
                    here->CNTf3d = here->CNTCbd * 
                            model->CNTbulkJctBotGradingCoeff * sarg/arg/
                            here->CNTtBulkPot + here->CNTCbdsw * 
                            model->CNTbulkJctSideGradingCoeff * sargsw/arg /
                            here->CNTtBulkPot;
                    here->CNTf4d = here->CNTCbd*here->CNTtBulkPot*
                            (1-arg*sarg)/ (1-model->CNTbulkJctBotGradingCoeff)
                            + here->CNTCbdsw*here->CNTtBulkPot*(1-arg*sargsw)/
                            (1-model->CNTbulkJctSideGradingCoeff)
                            -here->CNTf3d/2*
                            (here->CNTtDepCap*here->CNTtDepCap)
                            -here->CNTtDepCap * here->CNTf2d;
                }
                if(vbsOp >= here->CNTtDepCap){
                    arg = 1-model->CNTfwdCapDepCoeff;
                    sarg = exp( (-model->CNTbulkJctBotGradingCoeff) *
                            log(arg) );
                    sargsw = exp( (-model->CNTbulkJctSideGradingCoeff) * 
                            log(arg) );
                    here->CNTf2s = here->CNTCbs*(1-model->CNTfwdCapDepCoeff*
                        (1+model->CNTbulkJctBotGradingCoeff))* sarg/arg
                        +  here->CNTCbssw*(1-model->CNTfwdCapDepCoeff*
                        (1+model->CNTbulkJctSideGradingCoeff))*
                        sargsw/arg;
                    here->CNTf3s = here->CNTCbs * 
                            model->CNTbulkJctBotGradingCoeff * sarg/arg/
                            here->CNTtBulkPot + here->CNTCbssw * 
                            model->CNTbulkJctSideGradingCoeff * sargsw/arg /
                            here->CNTtBulkPot;
                    here->CNTf4s = here->CNTCbs*
                            here->CNTtBulkPot*(1-arg*sarg)/
                            (1-model->CNTbulkJctBotGradingCoeff)
                            + here->CNTCbssw*here->CNTtBulkPot*(1-arg*sargsw)/
                            (1-model->CNTbulkJctSideGradingCoeff)
                            -here->CNTf3s/2*
                            (here->CNTtDepCap*here->CNTtDepCap)
                            -here->CNTtDepCap * here->CNTf2s;
                }

                *(ckt->CKTstate0 + here->CNTvbs) = vbsOp;
                *(ckt->CKTstate0 + here->CNTvbd) = vbdOp;

                if ((error = CNTload((GENmodel*)model,ckt)))
		  return(error);

                *(here->CNTsenCgs + 5) = here->CNTcgs;
                *(here->CNTsenCgd + 5) = here->CNTcgd;
                *(here->CNTsenCgb + 5) = here->CNTcgb;
                *(here->CNTsenCbd + 5) = here->CNTcapbd;
                *(here->CNTsenCbs + 5) = here->CNTcapbs;
                *(here->CNTsenGds + 5) = here->CNTgds;
                *(here->CNTsenGbs + 5) = here->CNTgbs;
                *(here->CNTsenGbd + 5) = here->CNTgbd;
                *(here->CNTsenGm + 5) = here->CNTgm;
                *(here->CNTsenGmbs + 5) = here->CNTgmbs;

                here->CNTw = A0;
                here->CNTdrainArea /= (1 + info->SENpertfac);
                here->CNTsourceArea /= (1 + info->SENpertfac);
            }

load:

            gds= *(here->CNTsenGds + flag);
            gbs= *(here->CNTsenGbs + flag);
            gbd= *(here->CNTsenGbd + flag);
            gm= *(here->CNTsenGm + flag);
            gmbs= *(here->CNTsenGmbs + flag);
            if(flag == 5){
                gdpr = here->CNTdrainConductance * Apert/A0;
                gspr = here->CNTsourceConductance * Apert/A0;
            }
            else{
                gdpr = here->CNTdrainConductance;
                gspr = here->CNTsourceConductance;
            }

            xcgs= *(here->CNTsenCgs + flag) * ckt->CKTomega;
            xcgd= *(here->CNTsenCgd + flag) * ckt->CKTomega;
            xcgb= *(here->CNTsenCgb + flag) * ckt->CKTomega;
            xbd= *(here->CNTsenCbd + flag) * ckt->CKTomega;
            xbs= *(here->CNTsenCbs + flag) * ckt->CKTomega;

#ifdef SENSDEBUG
            printf("flag = %d \n",flag);
            printf("gspr = %.7e , gdpr = %.7e , gds = %.7e, gbs = %.7e\n",
                    gspr,gdpr,gds,gbs);
            printf("gbd = %.7e , gm = %.7e , gmbs = %.7e\n",gbd,gm,gmbs);
            printf("xcgs = %.7e , xcgd = %.7e , xcgb = %.7e,", xcgs,xcgd,xcgb);
            printf("xbd = %.7e,xbs = %.7e\n",xbd,xbs);
#endif /* SENSDEBUG */
            cspr = gspr * vspr ;
            icspr = gspr * ivspr ;
            cdpr = gdpr * vdpr ;
            icdpr = gdpr * ivdpr ;
            cgs = ( - xcgs * ivgs );
            icgs =  xcgs * vgs ;
            cgd = ( - xcgd * ivgd );
            icgd =  xcgd * vgd ;
            cgb = ( - xcgb * ivgb );
            icgb =  xcgb * vgb ;
            cbs = ( gbs * vbs  - xbs * ivbs );
            icbs = ( xbs * vbs  + gbs * ivbs );
            cbd = ( gbd * vbd  - xbd * ivbd );
            icbd = ( xbd * vbd  + gbd * ivbd );
            cds = ( gds * vds  + xnrm * (gm * vgs + gmbs * vbs) 
                    - xrev * (gm * vgd + gmbs * vbd) );
            icds = ( gds * ivds + xnrm * (gm * ivgs + gmbs * ivbs)
                    - xrev * (gm * ivgd + gmbs * ivbd) );

            cs = cspr;
            ics = icspr;
            csprm = ( -cspr - cgs - cbs - cds ) ;
            icsprm = ( -icspr - icgs - icbs - icds ) ;
            cd = cdpr;
            icd = icdpr;
            cdprm = ( -cdpr - cgd - cbd + cds ) ;
            icdprm = ( -icdpr - icgd - icbd + icds ) ;
            cg = cgs + cgd + cgb ;
            icg = icgs + icgd + icgb ;
            cb = cbs + cbd - cgb ;
            icb = icbs + icbd - icgb ;

#ifdef SENSDEBUG
            printf("vbs = %.7e , vbd = %.7e , vgb = %.7e\n",vbs,vbd,vgb);
            printf("ivbs = %.7e , ivbd = %.7e , ivgb = %.7e\n",ivbs,ivbd,ivgb);
            printf("cbs = %.7e , cbd = %.7e , cgb = %.7e\n",cbs,cbd,cgb);
            printf("cb = %.7e , cg = %.7e , cs = %.7e\n",cb,cg,cs);
            printf("csprm = %.7e, cd = %.7e, cdprm = %.7e\n",csprm,cd,cdprm);
            printf("icb = %.7e , icg = %.7e , ics = %.7e\n",icb,icg,ics);
            printf("icsprm = %.7e, icd = %.7e, icdprm = %.7e\n",
                    icsprm,icd,icdprm);
#endif /* SENSDEBUG */
            for(iparmno = 1;iparmno<=info->SENparms;iparmno++){
                if((flag == 4) && (iparmno != here->CNTsenParmNo)) continue;
                if( (flag == 5) && (iparmno != (here->CNTsenParmNo +
                        here->CNTsens_l))) continue;

                switch(flag){
                case 1: 
                    DvDp = model->CNTtype * 
                            (info->SEN_Sap[here->CNTbNode][iparmno]
                            -  info->SEN_Sap[here->CNTsNodePrime][iparmno]);
                    break;
                case 2: 
                    DvDp = model->CNTtype * 
                            ( info->SEN_Sap[here->CNTbNode][iparmno]
                            -  info->SEN_Sap[here->CNTdNodePrime][iparmno]);
                    break;
                case 3: 
                    DvDp = model->CNTtype * 
                            ( info->SEN_Sap[here->CNTgNode][iparmno]
                            -  info->SEN_Sap[here->CNTbNode][iparmno]);
                    break;
                case 4: 
                    DvDp = 1;
                    break;
                case 5: 
                    DvDp = 1;
                    break;
                }
                *(info->SEN_RHS[here->CNTbNode] + iparmno) -=  
                        ( cb  - cb0) * DELAinv * DvDp;
                *(info->SEN_iRHS[here->CNTbNode] + iparmno) -=  
                        ( icb  - icb0) * DELAinv * DvDp;

                *(info->SEN_RHS[here->CNTgNode] + iparmno) -=  
                        ( cg  - cg0) * DELAinv * DvDp;
                *(info->SEN_iRHS[here->CNTgNode] + iparmno) -=  
                        ( icg  - icg0) * DELAinv * DvDp;

                if(here->CNTsNode != here->CNTsNodePrime){
                    *(info->SEN_RHS[here->CNTsNode] + iparmno) -=  
                            ( cs  - cs0) * DELAinv * DvDp;
                    *(info->SEN_iRHS[here->CNTsNode] + iparmno) -=  
                            ( ics  - ics0) * DELAinv * DvDp;
                }

                *(info->SEN_RHS[here->CNTsNodePrime] + iparmno) -=  
                        ( csprm  - csprm0) * DELAinv * DvDp;
                *(info->SEN_iRHS[here->CNTsNodePrime] + iparmno) -=  
                        ( icsprm  - icsprm0) * DELAinv * DvDp;

                if(here->CNTdNode != here->CNTdNodePrime){
                    *(info->SEN_RHS[here->CNTdNode] + iparmno) -=  
                            ( cd  - cd0) * DELAinv * DvDp;
                    *(info->SEN_iRHS[here->CNTdNode] + iparmno) -=  
                            ( icd  - icd0) * DELAinv * DvDp;
                }

                *(info->SEN_RHS[here->CNTdNodePrime] + iparmno) -=  
                        ( cdprm  - cdprm0) * DELAinv * DvDp;
                *(info->SEN_iRHS[here->CNTdNodePrime] + iparmno) -=  
                        ( icdprm  - icdprm0) * DELAinv * DvDp;
#ifdef SENSDEBUG
                printf("after loading\n");  
                printf("DvDp = %.5e , DELAinv = %.5e ,flag = %d ,",
                        DvDp,DELAinv,flag);
                printf("iparmno = %d,senparmno = %d\n",
                        iparmno,here->CNTsenParmNo);
                printf("A0 = %.5e , Apert = %.5e ,CONSTvt = %.5e \n",
                        A0,Apert,here->CNTtVto);
                printf("senb = %.7e + j%.7e ",
                        *(info->SEN_RHS[here->CNTbNode] + iparmno),
                        *(info->SEN_iRHS[here->CNTbNode] + iparmno));
                printf("seng = %.7e + j%.7e ",
                        *(info->SEN_RHS[here->CNTgNode] + iparmno),
                        *(info->SEN_iRHS[here->CNTgNode] + iparmno));
                printf("sens = %.7e + j%.7e ",
                        *(info->SEN_RHS[here->CNTsNode] + iparmno),
                        *(info->SEN_iRHS[here->CNTsNode] + iparmno));
                printf("sensprm = %.7e + j%.7e ",
                        *(info->SEN_RHS[here->CNTsNodePrime] + iparmno),
                        *(info->SEN_iRHS[here->CNTsNodePrime] + iparmno));
                printf("send = %.7e + j%.7e ",
                        *(info->SEN_RHS[here->CNTdNode] + iparmno),
                        *(info->SEN_iRHS[here->CNTdNode] + iparmno));
                printf("sendprm = %.7e + j%.7e ",
                        *(info->SEN_RHS[here->CNTdNodePrime] + iparmno),
                        *(info->SEN_iRHS[here->CNTdNodePrime] + iparmno));
#endif /* SENSDEBUG */

            }
            switch(flag){
            case 1: 
                goto pertvbd ;
            case 2: 
                goto pertvgb ; 
            case 3: 
                goto pertl ;
            case 4: 
                goto pertw ;
            case 5: 
                break; 
            }
next:                   
            ;

            /* put the unperturbed values back into the state vector */
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
    printf("CNTsenacload end\n");
#endif /* SENSDEBUG */
    return(OK);
}


