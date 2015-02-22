/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "cntdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int
CNTconvTest(GENmodel *inModel, CKTcircuit *ckt)
{
    CNTmodel *model = (CNTmodel *)inModel;
    CNTinstance *here;
    double delvbs;
    double delvbd;
    double delvgs;
    double delvds;
    double delvgd;
    double cbhat;
    double cdhat;
    double vbs;
    double vbd;
    double vgs;
    double vds;
    double vgd;
    double vgdo;
    double tol;

    for( ; model != NULL; model = model->CNTnextModel) {
        for(here = model->CNTinstances; here!= NULL;
                here = here->CNTnextInstance) {
	    if (here->CNTowner != ARCHme) continue;
        
            vbs = model->CNTtype * ( 
                *(ckt->CKTrhs+here->CNTbNode) -
                *(ckt->CKTrhs+here->CNTsNodePrime));
            vgs = model->CNTtype * ( 
                *(ckt->CKTrhs+here->CNTgNode) -
                *(ckt->CKTrhs+here->CNTsNodePrime));
            vds = model->CNTtype * ( 
                *(ckt->CKTrhs+here->CNTdNodePrime) -
                *(ckt->CKTrhs+here->CNTsNodePrime));
            vbd=vbs-vds;
            vgd=vgs-vds;
            vgdo = *(ckt->CKTstate0 + here->CNTvgs) -
                *(ckt->CKTstate0 + here->CNTvds);
            delvbs = vbs - *(ckt->CKTstate0 + here->CNTvbs);
            delvbd = vbd - *(ckt->CKTstate0 + here->CNTvbd);
            delvgs = vgs - *(ckt->CKTstate0 + here->CNTvgs);
            delvds = vds - *(ckt->CKTstate0 + here->CNTvds);
            delvgd = vgd-vgdo;

            /* these are needed for convergence testing */

            if (here->CNTmode >= 0) {
                cdhat=
                    here->CNTcd-
                    here->CNTgbd * delvbd +
                    here->CNTgmbs * delvbs +
                    here->CNTgm * delvgs + 
                    here->CNTgds * delvds ;
            } else {
                cdhat=
                    here->CNTcd -
                    ( here->CNTgbd -
                    here->CNTgmbs) * delvbd -
                    here->CNTgm * delvgd + 
                    here->CNTgds * delvds ;
            }
            cbhat=
                here->CNTcbs +
                here->CNTcbd +
                here->CNTgbd * delvbd +
                here->CNTgbs * delvbs ;
            /*
             *  check convergence
             */
            tol=ckt->CKTreltol*MAX(fabs(cdhat),fabs(here->CNTcd))+
                    ckt->CKTabstol;
            if (fabs(cdhat-here->CNTcd) >= tol) { 
                ckt->CKTnoncon++;
		ckt->CKTtroubleElt = (GENinstance *) here;
                return(OK); /* no reason to continue, we haven't converged */
            } else {
                tol=ckt->CKTreltol*
                        MAX(fabs(cbhat),fabs(here->CNTcbs+here->CNTcbd))+ 
                        ckt->CKTabstol;
                if (fabs(cbhat-(here->CNTcbs+here->CNTcbd)) > tol) {
                    ckt->CKTnoncon++;
		    ckt->CKTtroubleElt = (GENinstance *) here;
                    return(OK); /* no reason to continue, we haven't converged*/
                }
            }
        }
    }
    return(OK);
}
