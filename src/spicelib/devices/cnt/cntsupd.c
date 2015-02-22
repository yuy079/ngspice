/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles

This function is obsolete (was used by an old sensitivity analysis)
**********/

#include "ngspice.h"
#include "smpdefs.h"
#include "cktdefs.h"
#include "cntdefs.h"
#include "sperror.h"
#include "suffix.h"

/* update the  charge sensitivities and their derivatives */

int
CNTsUpdate(GENmodel *inModel, CKTcircuit *ckt)
{
    CNTmodel *model = (CNTmodel *)inModel;
    CNTinstance *here;
    int    iparmno;
    double sb;
    double sg;
    double sdprm;
    double ssprm;
    double sxpgs;
    double sxpgd;
    double sxpbs;
    double sxpbd;
    double sxpgb;
    double dummy1;
    double dummy2;
    SENstruct *info;


    if(ckt->CKTtime == 0) return(OK);
    info = ckt->CKTsenInfo;

#ifdef SENSDEBUG
    printf("CNTsenupdate\n");
    printf("CKTtime = %.5e\n",ckt->CKTtime);
#endif /* SENSDEBUG */

    sxpgs = 0;
    sxpgd = 0;
    sxpbs = 0;
    sxpbd = 0;
    sxpgb = 0;
    dummy1 = 0;
    dummy2 = 0;

    /*  loop through all the CNT models */
    for( ; model != NULL; model = model->CNTnextModel ) {

        /* loop through all the instances of the model */
        for (here = model->CNTinstances; here != NULL ;
                here=here->CNTnextInstance) {
	    if (here->CNTowner != ARCHme) continue;

#ifdef SENSDEBUG
            printf("senupdate instance name %s\n",here->CNTname);
            printf("before loading\n");
            printf("CKTag[0] = %.2e,CKTag[1] = %.2e\n",
                    ckt->CKTag[0],ckt->CKTag[1]);
            printf("capgs = %.7e\n",here->CNTcgs);
            printf("capgd = %.7e\n",here->CNTcgd);
            printf("capgb = %.7e\n",here->CNTcgb);
            printf("capbs = %.7e\n",here->CNTcapbs);
            printf("capbd = %.7e\n",here->CNTcapbd);
#endif /* SENSDEBUG */

            for(iparmno = 1;iparmno<=info->SENparms;iparmno++){

                sb = *(info->SEN_Sap[here->CNTbNode] + iparmno);
                sg = *(info->SEN_Sap[here->CNTgNode] + iparmno);
                ssprm = *(info->SEN_Sap[here->CNTsNodePrime] + iparmno);
                sdprm = *(info->SEN_Sap[here->CNTdNodePrime] + iparmno);
#ifdef SENSDEBUG
                printf("iparmno = %d\n",iparmno);
                printf("sb = %.7e,sg = %.7e\n",sb,sg);
                printf("ssprm = %.7e,sdprm = %.7e\n",ssprm,sdprm);
#endif /* SENSDEBUG */

                sxpgs =  (sg - ssprm) * here->CNTcgs ;
                sxpgd =  (sg - sdprm) * here->CNTcgd ;
                sxpgb =  (sg - sb) * here->CNTcgb ;
                sxpbs =  (sb - ssprm) * here->CNTcapbs ;
                sxpbd =  (sb - sdprm) * here->CNTcapbd ;

                if(here->CNTsens_l && (iparmno == here->CNTsenParmNo)){
                    sxpgs += *(here->CNTdphigs_dl);
                    sxpgd += *(here->CNTdphigd_dl);
                    sxpbs += *(here->CNTdphibs_dl);
                    sxpbd += *(here->CNTdphibd_dl);
                    sxpgb += *(here->CNTdphigb_dl);
                }
                if(here->CNTsens_w && 
                        (iparmno == (here->CNTsenParmNo+here->CNTsens_l))){
                    sxpgs += *(here->CNTdphigs_dw);
                    sxpgd += *(here->CNTdphigd_dw);
                    sxpbs += *(here->CNTdphibs_dw);
                    sxpbd += *(here->CNTdphibd_dw);
                    sxpgb += *(here->CNTdphigb_dw);
                }
                if(ckt->CKTmode & MODEINITTRAN) {
                    *(ckt->CKTstate1 + here->CNTsensxpgs + 
                            10 * (iparmno - 1)) = sxpgs;
                    *(ckt->CKTstate1 + here->CNTsensxpgd + 
                            10 * (iparmno - 1)) = sxpgd;
                    *(ckt->CKTstate1 + here->CNTsensxpbs + 
                            10 * (iparmno - 1)) = sxpbs;
                    *(ckt->CKTstate1 + here->CNTsensxpbd + 
                            10 * (iparmno - 1)) = sxpbd;
                    *(ckt->CKTstate1 + here->CNTsensxpgb + 
                            10 * (iparmno - 1)) = sxpgb;
                    *(ckt->CKTstate1 + here->CNTsensxpgs + 
                            10 * (iparmno - 1) + 1) = 0;
                    *(ckt->CKTstate1 + here->CNTsensxpgd + 
                            10 * (iparmno - 1) + 1) = 0;
                    *(ckt->CKTstate1 + here->CNTsensxpbs + 
                            10 * (iparmno - 1) + 1) = 0;
                    *(ckt->CKTstate1 + here->CNTsensxpbd + 
                            10 * (iparmno - 1) + 1) = 0;
                    *(ckt->CKTstate1 + here->CNTsensxpgb + 
                            10 * (iparmno - 1) + 1) = 0;
                    goto next;
                }

                *(ckt->CKTstate0 + here->CNTsensxpgs + 
                        10 * (iparmno - 1)) = sxpgs;
                *(ckt->CKTstate0 + here->CNTsensxpgd + 
                        10 * (iparmno - 1)) = sxpgd;
                *(ckt->CKTstate0 + here->CNTsensxpbs + 
                        10 * (iparmno - 1)) = sxpbs;
                *(ckt->CKTstate0 + here->CNTsensxpbd + 
                        10 * (iparmno - 1)) = sxpbd;
                *(ckt->CKTstate0 + here->CNTsensxpgb + 
                        10 * (iparmno - 1)) = sxpgb;

                NIintegrate(ckt,&dummy1,&dummy2,here->CNTcgs,
                        here->CNTsensxpgs + 10*(iparmno -1));

                NIintegrate(ckt,&dummy1,&dummy2,here->CNTcgd,
                        here->CNTsensxpgd + 10*(iparmno -1));

                NIintegrate(ckt,&dummy1,&dummy2,here->CNTcgb,
                        here->CNTsensxpgb + 10*(iparmno -1));

                NIintegrate(ckt,&dummy1,&dummy2,here->CNTcapbs,
                        here->CNTsensxpbs + 10*(iparmno -1));

                NIintegrate(ckt,&dummy1,&dummy2,here->CNTcapbd,
                        here->CNTsensxpbd + 10*(iparmno -1));
next:   
                ;
#ifdef SENSDEBUG
                printf("after loading\n");
                printf("sxpgs = %.7e,sdotxpgs = %.7e\n",
                        sxpgs,*(ckt->CKTstate0 + here->CNTsensxpgs + 
                        10 * (iparmno - 1) + 1));
                printf("sxpgd = %.7e,sdotxpgd = %.7e\n",
                        sxpgd,*(ckt->CKTstate0 + here->CNTsensxpgd + 
                        10 * (iparmno - 1) + 1));
                printf("sxpgb = %.7e,sdotxpgb = %.7e\n",
                        sxpgb,*(ckt->CKTstate0 + here->CNTsensxpgb + 
                        10 * (iparmno - 1) + 1));
                printf("sxpbs = %.7e,sdotxpbs = %.7e\n",
                        sxpbs,*(ckt->CKTstate0 + here->CNTsensxpbs + 
                        10 * (iparmno - 1) + 1));
                printf("sxpbd = %.7e,sdotxpbd = %.7e\n",
                        sxpbd,*(ckt->CKTstate0 + here->CNTsensxpbd + 
                        10 * (iparmno - 1) + 1));
#endif /* SENSDEBUG */
            }
        }
    }
#ifdef SENSDEBUG
    printf("CNTsenupdate end\n");
#endif /* SENSDEBUG */
    return(OK);

}

