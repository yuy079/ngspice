/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/

#include <math.h>
#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/devdefs.h"
#include "cntdefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/const.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

static double CSa[5]={-2.384095e10,7.137675e10,-4.753580e10,0.0,0.0};
static double CSb[5]={-3.576143e10,7.849982e10,-2.852148e10,0.0,0.0};
static double CSc[5]={-1.930062e10,2.640388e10,-5.702512e9,0.0,-1.419906e9};
static double CSd[5]={-3.263224e9,2.830709e9,-3.799295e8,0.0,-2.831046e8};

int
CNTload(GENmodel *inModel, CKTcircuit *ckt)
        /* actually load the current value into the 
         * sparse matrix previously provided 
         */
{
    CNTmodel *model = (CNTmodel *) inModel;
    CNTinstance *here;
    double Beta;
    double DrainSatCur;
    double EffectiveLength;
    double GateBulkOverlapCap;
    double GateDrainOverlapCap;
    double GateSourceOverlapCap;
    double OxideCap;
    double SourceSatCur;
    double arg;
    double cbhat;
    double cdhat;
    double cdrain;
    double cdreq;
    double ceq;
    double ceqbd;
    double ceqbs;
    double ceqgb;
    double ceqgd;
    double ceqgs;
    double delvbd;
    double delvbs;
    double delvds;
    double delvgd;
    double delvgs;
    double evbd;
    double evbs;
    double gcgb;
    double gcgd;
    double gcgs;
    double geq;
    double sarg;
    double sargsw;
    double vbd;
    double vbs;
    double vds;
    double vdsat;
    double vgb1;
    double vgb;
    double vgd1;
    double vgd;
    double vgdo;
    double vgs1;
    double vgs;
    double von;
    double vt;
    double xfact = 0.0;
    int xnrm;
    int xrev;
    double capgs = 0.0;   /* total gate-source capacitance */
    double capgd = 0.0;   /* total gate-drain capacitance */
    double capgb = 0.0;   /* total gate-bulk capacitance */
    int Check;

    double EF;
    double VD;
    double VS;
    double VG;
    double Vdi;
    double Vsi;
    double Vgi;
    double Vdsc;
    double Eft;
    double Efi;
    double N0;
    double c;
    double Cox;
    double Cge;
    double Cse;
    double Cde;
    double Ctot;
    double qC;
    double qCN0;
    double dcnt;
    double xmax=-0.2;
    double xmin=-0.5;
    double inte;
    double Vsc;
    double q=1.6e-19;
    double pi=3.1415926;
    double e0=8.854e-12;
    double t0=1.5e-9;
    double k=3.9;
    double h=6.625e-34;
    double KB=1.38e-23;
    double T=300.0;
    double acc=1.42e-10;
    double Vcc;
    double Eg;
    double D0;
    int QsRange;
    int QdRange;


#ifndef NOBYPASS    
    double tempv;
#endif /*NOBYPASS*/    
    int error;
 #ifdef CAPBYPASS
    int senflag;
#endif /*CAPBYPASS*/           
    int SenCond;



#ifdef CAPBYPASS
    senflag = 0;
    if(ckt->CKTsenInfo && ckt->CKTsenInfo->SENstatus == PERTURBATION &&
        (ckt->CKTsenInfo->SENmode & (ACSEN | TRANSEN))) {
        senflag = 1;
    }
#endif /* CAPBYPASS */ 

    /*  loop through all the CNT device models */
    for( ; model != NULL; model = model->CNTnextModel ) {

        /* loop through all the instances of the model */
        for (here = model->CNTinstances; here != NULL ;
	     here=here->CNTnextInstance) {
	    /*if (here->CNTowner != ARCHme) continue;*//*cuidado*/

            vt = CONSTKoverQ * here->CNTtemp;
            Check=1;
            if(ckt->CKTsenInfo){
#ifdef SENSDEBUG
                printf("CNTload \n");
#endif /* SENSDEBUG */

                if((ckt->CKTsenInfo->SENstatus == PERTURBATION)&&
		   (here->CNTsenPertFlag == OFF))continue;

            }
            SenCond = ckt->CKTsenInfo && here->CNTsenPertFlag;

/*
  
*/

            /* first, we compute a few useful values - these could be
             * pre-computed, but for historical reasons are still done
             * here.  They may be moved at the expense of instance size
             */

            EffectiveLength=here->CNTl - 2*model->CNTlatDiff;
            
            if( (here->CNTtSatCurDens == 0) || 
                    (here->CNTdrainArea == 0) ||
                    (here->CNTsourceArea == 0)) {
                DrainSatCur = here->CNTm * here->CNTtSatCur;
                SourceSatCur = here->CNTm * here->CNTtSatCur;
            } else {
                DrainSatCur = here->CNTtSatCurDens * 
                        here->CNTm * here->CNTdrainArea;
                SourceSatCur = here->CNTtSatCurDens * 
                        here->CNTm * here->CNTsourceArea;
            }
            GateSourceOverlapCap = model->CNTgateSourceOverlapCapFactor * 
                    here->CNTm * here->CNTw;
            GateDrainOverlapCap = model->CNTgateDrainOverlapCapFactor * 
                    here->CNTm * here->CNTw;
            GateBulkOverlapCap = model->CNTgateBulkOverlapCapFactor * 
                    here->CNTm * EffectiveLength;
            Beta = here->CNTtTransconductance * here->CNTm *
                    here->CNTw/EffectiveLength;
            OxideCap = model->CNToxideCapFactor * EffectiveLength * 
                    here->CNTm * here->CNTw;
           
            /* 
             * ok - now to do the start-up operations
             *
             * we must get values for vbs, vds, and vgs from somewhere
             * so we either predict them or recover them from last iteration
             * These are the two most common cases - either a prediction
             * step or the general iteration step and they
             * share some code, so we put them first - others later on
             */

            if(SenCond){
#ifdef SENSDEBUG
                printf("CNTsenPertFlag = ON \n");
#endif /* SENSDEBUG */
                if((ckt->CKTsenInfo->SENmode == TRANSEN) &&
		   (ckt->CKTmode & MODEINITTRAN)) {
                    vgs = *(ckt->CKTstate1 + here->CNTvgs);
                    vds = *(ckt->CKTstate1 + here->CNTvds);
                    vbs = *(ckt->CKTstate1 + here->CNTvbs);
                    vbd = *(ckt->CKTstate1 + here->CNTvbd);
                    vgb = vgs - vbs;
                    vgd = vgs - vds;
                }
                else if (ckt->CKTsenInfo->SENmode == ACSEN){
                    vgb = model->CNTtype * ( 
                        *(ckt->CKTrhsOp+here->CNTgNode) -
                        *(ckt->CKTrhsOp+here->CNTbNode));
                    vbs = *(ckt->CKTstate0 + here->CNTvbs);
                    vbd = *(ckt->CKTstate0 + here->CNTvbd);
                    vgd = vgb + vbd ;
                    vgs = vgb + vbs ;
                    vds = vbs - vbd ;
                }
                else{
                    vgs = *(ckt->CKTstate0 + here->CNTvgs);
                    vds = *(ckt->CKTstate0 + here->CNTvds);
                    vbs = *(ckt->CKTstate0 + here->CNTvbs);
                    vbd = *(ckt->CKTstate0 + here->CNTvbd);
                    vgb = vgs - vbs;
                    vgd = vgs - vds;
                }
#ifdef SENSDEBUG
                printf(" vbs = %.7e ,vbd = %.7e,vgb = %.7e\n",vbs,vbd,vgb);
                printf(" vgs = %.7e ,vds = %.7e,vgd = %.7e\n",vgs,vds,vgd);
#endif /* SENSDEBUG */
                goto next1;
            }


            if((ckt->CKTmode & (MODEINITFLOAT | MODEINITPRED | MODEINITSMSIG
				| MODEINITTRAN)) ||
	       ( (ckt->CKTmode & MODEINITFIX) && (!here->CNToff) )  ) {
#ifndef PREDICTOR
                if(ckt->CKTmode & (MODEINITPRED | MODEINITTRAN) ) {

                    /* predictor step */

                    xfact=ckt->CKTdelta/ckt->CKTdeltaOld[1];
                    *(ckt->CKTstate0 + here->CNTvbs) = 
			*(ckt->CKTstate1 + here->CNTvbs);
                    vbs = (1+xfact)* (*(ckt->CKTstate1 + here->CNTvbs))
			-(xfact * (*(ckt->CKTstate2 + here->CNTvbs)));
                    *(ckt->CKTstate0 + here->CNTvgs) = 
			*(ckt->CKTstate1 + here->CNTvgs);
                    vgs = (1+xfact)* (*(ckt->CKTstate1 + here->CNTvgs))
			-(xfact * (*(ckt->CKTstate2 + here->CNTvgs)));
                    *(ckt->CKTstate0 + here->CNTvds) = 
			*(ckt->CKTstate1 + here->CNTvds);
                    vds = (1+xfact)* (*(ckt->CKTstate1 + here->CNTvds))
			-(xfact * (*(ckt->CKTstate2 + here->CNTvds)));
                    *(ckt->CKTstate0 + here->CNTvbd) = 
			*(ckt->CKTstate0 + here->CNTvbs)-
			*(ckt->CKTstate0 + here->CNTvds);
                } else {
#endif /* PREDICTOR */

                    /* general iteration */

                    vbs = model->CNTtype * ( 
                        *(ckt->CKTrhsOld+here->CNTbNode) -
                        *(ckt->CKTrhsOld+here->CNTsNodePrime));
                    vgs = model->CNTtype * ( 
                        *(ckt->CKTrhsOld+here->CNTgNode) -
                        *(ckt->CKTrhsOld+here->CNTsNodePrime));
                    vds = model->CNTtype * ( 
                        *(ckt->CKTrhsOld+here->CNTdNodePrime) -
                        *(ckt->CKTrhsOld+here->CNTsNodePrime));
#ifndef PREDICTOR
                }
#endif /* PREDICTOR */

                /* now some common crunching for some more useful quantities */

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
  
*/

#ifndef NOBYPASS
                /* now lets see if we can bypass (ugh) */
                tempv = (MAX(fabs(cbhat),
			    fabs(here->CNTcbs + here->CNTcbd)) +
			 ckt->CKTabstol);
                if ((!(ckt->CKTmode &
		       (MODEINITPRED|MODEINITTRAN|MODEINITSMSIG))) &&
		    (ckt->CKTbypass) &&
		    (fabs(cbhat-(here->CNTcbs + 
				 here->CNTcbd)) < ckt->CKTreltol * tempv) &&
		    (fabs(delvbs) < (ckt->CKTreltol *
				     MAX(fabs(vbs),
					 fabs( *(ckt->CKTstate0 +
						 here->CNTvbs))) +
				     ckt->CKTvoltTol)) &&
		    (fabs(delvbd) < (ckt->CKTreltol *
				     MAX(fabs(vbd),
					 fabs(*(ckt->CKTstate0 +
						here->CNTvbd))) +
				     ckt->CKTvoltTol)) &&
		    (fabs(delvgs) < (ckt->CKTreltol *
				     MAX(fabs(vgs),
					 fabs(*(ckt->CKTstate0 +
						here->CNTvgs)))+
				     ckt->CKTvoltTol)) &&
		    (fabs(delvds) < (ckt->CKTreltol *
				     MAX(fabs(vds),
					 fabs(*(ckt->CKTstate0 +
						here->CNTvds))) +
				     ckt->CKTvoltTol)) &&
		    (fabs(cdhat- here->CNTcd) < (ckt->CKTreltol *
						  MAX(fabs(cdhat),
						      fabs(here->CNTcd)) +
						  ckt->CKTabstol))) {
		  /* bypass code */
		  /* nothing interesting has changed since last
		   * iteration on this device, so we just
		   * copy all the values computed last iteration out
		   * and keep going
		   */
		  vbs = *(ckt->CKTstate0 + here->CNTvbs);
		  vbd = *(ckt->CKTstate0 + here->CNTvbd);
		  vgs = *(ckt->CKTstate0 + here->CNTvgs);
		  vds = *(ckt->CKTstate0 + here->CNTvds);
		  vgd = vgs - vds;
		  vgb = vgs - vbs;
		  cdrain = here->CNTmode * (here->CNTcd + here->CNTcbd);
		  if(ckt->CKTmode & (MODETRAN | MODETRANOP)) {
		    capgs = ( *(ckt->CKTstate0+here->CNTcapgs)+ 
			      *(ckt->CKTstate1+here->CNTcapgs) +
			      GateSourceOverlapCap );
		    capgd = ( *(ckt->CKTstate0+here->CNTcapgd)+ 
			      *(ckt->CKTstate1+here->CNTcapgd) +
			      GateDrainOverlapCap );
		    capgb = ( *(ckt->CKTstate0+here->CNTcapgb)+ 
			      *(ckt->CKTstate1+here->CNTcapgb) +
			      GateBulkOverlapCap );
		    
		    if(ckt->CKTsenInfo){
		      here->CNTcgs = capgs;
		      here->CNTcgd = capgd;
		      here->CNTcgb = capgb;
		    }
		  }
		  goto bypass;
		}
#endif /*NOBYPASS*/

/*
  
*/

                /* ok - bypass is out, do it the hard way */

                von = model->CNTtype * here->CNTvon;

#ifndef NODELIMITING
                /* 
                 * limiting
                 *  we want to keep device voltages from changing
                 * so fast that the exponentials churn out overflows
                 * and similar rudeness
                 */

                if(*(ckt->CKTstate0 + here->CNTvds) >=0) {
                    vgs = DEVfetlim(vgs,*(ckt->CKTstate0 + here->CNTvgs)
				    ,von);
                    vds = vgs - vgd;
                    vds = DEVlimvds(vds,*(ckt->CKTstate0 + here->CNTvds));
                    vgd = vgs - vds;
                } else {
                    vgd = DEVfetlim(vgd,vgdo,von);
                    vds = vgs - vgd;
                    if(!(ckt->CKTfixLimit)) {
                        vds = -DEVlimvds(-vds,-(*(ckt->CKTstate0 + 
						  here->CNTvds)));
                    }
                    vgs = vgd + vds;
                }
                if(vds >= 0) {
                    vbs = DEVpnjlim(vbs,*(ckt->CKTstate0 + here->CNTvbs),
				    vt,here->CNTsourceVcrit,&Check);
                    vbd = vbs-vds;
                } else {
                    vbd = DEVpnjlim(vbd,*(ckt->CKTstate0 + here->CNTvbd),
				    vt,here->CNTdrainVcrit,&Check);
                    vbs = vbd + vds;
                }
#endif /*NODELIMITING*/
/*
  
*/

            } else {

                /* ok - not one of the simple cases, so we have to
                 * look at all of the possibilities for why we were
                 * called.  We still just initialize the three voltages
                 */

                if((ckt->CKTmode & MODEINITJCT) && !here->CNToff) {
                    vds= model->CNTtype * here->CNTicVDS;
                    vgs= model->CNTtype * here->CNTicVGS;
                    vbs= model->CNTtype * here->CNTicVBS;
                    if((vds==0) && (vgs==0) && (vbs==0) && 
		       ((ckt->CKTmode & 
			 (MODETRAN|MODEDCOP|MODEDCTRANCURVE)) ||
			(!(ckt->CKTmode & MODEUIC)))) {
                        vbs = -1;
                        vgs = model->CNTtype * here->CNTtVto;
                        vds = 0;
                    }
                } else {
                    vbs=vgs=vds=0;
                } 
            }
/*
  
*/

            /*
             * now all the preliminaries are over - we can start doing the
             * real work
             */
            vbd = vbs - vds;
            vgd = vgs - vds;
            vgb = vgs - vbs;


            /*
             * bulk-source and bulk-drain diodes
             *   here we just evaluate the ideal diode current and the
             *   corresponding derivative (conductance).
             */
next1:     if(vbs <= -3*vt) {
                here->CNTgbs = ckt->CKTgmin;
                here->CNTcbs = here->CNTgbs*vbs-SourceSatCur;
            } else {
                evbs = exp(MIN(MAX_EXP_ARG,vbs/vt));
                here->CNTgbs = SourceSatCur*evbs/vt + ckt->CKTgmin;
                here->CNTcbs = SourceSatCur*(evbs-1) + ckt->CKTgmin*vbs;
            }

			Efi=here->CNTef*q;
            if(Efi>0) {
                here->CNTgbd = ckt->CKTgmin;
                here->CNTcbd = here->CNTgbd*vbd-DrainSatCur;
                evbd = exp(MIN(MAX_EXP_ARG,vbd/vt));
            } else {
                evbd = exp(MIN(MAX_EXP_ARG,vbd/vt));
                here->CNTgbd = DrainSatCur*evbd/vt + ckt->CKTgmin;
                here->CNTcbd = DrainSatCur*(evbd-1) + ckt->CKTgmin*vbd;
            }
            /* now to determine whether the user was able to correctly
             * identify the source and drain of his device
             */
            if(vds >= 0) {
                /* normal mode */
                here->CNTmode = 1;
            } else {
                /* inverse mode */
                here->CNTmode = -1;
            }
/*
  
*/

            {
		/*
		 *     this block of code evaluates the drain current and its 
		 *     derivatives using the shichman-hodges model and the 
		 *     charges associated with the gate, channel and bulk for 
		 *     mosfets
		 *
		 */

		/* the following 4 variables are local to this code block until 
		 * it is obvious that they can be made global 
		 */
		double arg;
		double betap;
		double sarg;
		double vgst;

                if ((here->CNTmode==1?vbs:vbd) <= 0 ) {
                    sarg=sqrt(here->CNTtPhi-(here->CNTmode==1?vbs:vbd));
                } else {
                    sarg=sqrt(here->CNTtPhi);
                    sarg=sarg-(here->CNTmode==1?vbs:vbd)/(sarg+sarg);
                    sarg=MAX(0,sarg);
                }
                von=(here->CNTtVbi*model->CNTtype)+model->CNTgamma*sarg;
                vgst=(here->CNTmode==1?vgs:vgd)-von;
                vdsat=MAX(vgst,0);
                if (sarg <= 0) {
                    arg=0;
                } else {
                    arg=model->CNTgamma/(sarg+sarg);
                }
                if (vgst <= 0) {
                    /*
                     *     cutoff region
                     */
                    cdrain=0;
                    here->CNTgm=0;
                    here->CNTgds=0;
                    here->CNTgmbs=0;
                } else{ 			/* cubic spline */
			dcnt= here->CNTdia;
//printf("%lf\n",dcnt*1000000);
			Vcc=3.0*q;
			D0=8.0/3.0/pi/acc/Vcc;
			Eg=2.0*acc*Vcc/dcnt;
			Cox=2.0*pi*k*e0/log((t0+dcnt/2.0)*2.0/dcnt);
			Cge=Cox;Cse=0.097*Cox;Cde=0.040*Cox;Ctot=Cge+Cse+Cde;
			Efi=here->CNTef*q;
			Vdi= -vbd; Vsi=-vbs; Vgi= vgb;
			if (Efi>0){
				Eft=-Efi;VD=-Vdi;VS=-Vsi;VG=-Vgi;
			} 
			else {
				Eft=Efi;VD=Vdi;VS=Vsi;VG=Vgi;
   			}
			EF=Eft/q;
			N0=Ncnt(D0*q,Eg/q,KB*T/q,EF);
			c=-q*(VG*Cge+VS*Cse+VD*Cde)/Ctot;
			Vdsc=VD-VS;
//printf("%lf %lf %lf\n",VD,VG,Vdsc);
			qC=q*q/Ctot;
			qCN0=qC*N0;
			inte=(xmax-xmin)/(4.0-1.0);
			QsRange=Range1(Vdsc,q,c,qC,qCN0,xmax,xmin,inte,4);
			QdRange=Range2(Vdsc,q,c,qC,qCN0,xmax,xmin,inte,4);
			Vsc=Root(Vdsc,q,c,qC,qCN0,QsRange,QdRange);
printf("%lf %lf %lf %lf\n",VD,VS,VG,Vdsc);
			//if (Efi<=0.0)
				cdrain=4.0*q*KB*T/h*(log(1.0+exp(q*(EF-Vsc)/KB/T))-log(1.0+exp(q*(EF-Vsc-Vdsc)/KB/T)));
			//else 
				//cdrain=-4.0*q*KB*T/h*(log(1.0+exp(q*(EF-Vsc)/KB/T))-log(1.0+exp(q*(EF-Vsc-Vdsc)/KB/T)));
printf("%d %d %lf  %lf\n",QsRange, QdRange, Vsc,cdrain);
			betap=Beta*(1+model->CNTlambda*(vds*here->CNTmode));
                    	here->CNTgm=betap*(vds*here->CNTmode);
                        here->CNTgds=betap*(vgst-(vds*here->CNTmode))+
			    model->CNTlambda*Beta*
			    (vds*here->CNTmode)*
			    (vgst-.5*(vds*here->CNTmode));
                        here->CNTgmbs=here->CNTgm*arg;
                }
                /*
                 *     finished
                 */
            }
/*
  
*/

            /* now deal with n vs p polarity */

            here->CNTvon = model->CNTtype * von;
            here->CNTvdsat = model->CNTtype * vdsat;
            /* line 490 */
            /*
             *  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
             */
            here->CNTcd=here->CNTmode * cdrain - here->CNTcbd;

            if (ckt->CKTmode & (MODETRAN | MODETRANOP | MODEINITSMSIG)) {
                /* 
                 * now we do the hard part of the bulk-drain and bulk-source
                 * diode - we evaluate the non-linear capacitance and
                 * charge
                 *
                 * the basic equations are not hard, but the implementation
                 * is somewhat long in an attempt to avoid log/exponential
                 * evaluations
                 */
                /*
                 *  charge storage elements
                 *
                 *.. bulk-drain and bulk-source depletion capacitances
                 */
#ifdef CAPBYPASS
                if(((ckt->CKTmode & (MODEINITPRED | MODEINITTRAN) ) ||
                        fabs(delvbs) >= ckt->CKTreltol * MAX(fabs(vbs),
                        fabs(*(ckt->CKTstate0+here->CNTvbs)))+
                        ckt->CKTvoltTol)|| senflag)
#endif /*CAPBYPASS*/
		 
                {
                    /* can't bypass the diode capacitance calculations */
#ifdef CAPZEROBYPASS		    
                    if(here->CNTCbs != 0 || here->CNTCbssw != 0 ) {
#endif /*CAPZEROBYPASS*/		    
			if (vbs < here->CNTtDepCap){
			    arg=1-vbs/here->CNTtBulkPot;
			    /*
			     * the following block looks somewhat long and messy,
			     * but since most users use the default grading
			     * coefficients of .5, and sqrt is MUCH faster than an
			     * exp(log()) we use this special case code to buy time.
			     * (as much as 10% of total job time!)
			     */
#ifndef NOSQRT
			    if(model->CNTbulkJctBotGradingCoeff ==
			       model->CNTbulkJctSideGradingCoeff) {
				if(model->CNTbulkJctBotGradingCoeff == .5) {
				    sarg = sargsw = 1/sqrt(arg);
				} else {
				    sarg = sargsw =
                                        exp(-model->CNTbulkJctBotGradingCoeff*
					    log(arg));
				}
			    } else {
				if(model->CNTbulkJctBotGradingCoeff == .5) {
				    sarg = 1/sqrt(arg);
				} else {
#endif /*NOSQRT*/
				    sarg = exp(-model->CNTbulkJctBotGradingCoeff*
					       log(arg));
#ifndef NOSQRT
				}
				if(model->CNTbulkJctSideGradingCoeff == .5) {
				    sargsw = 1/sqrt(arg);
				} else {
#endif /*NOSQRT*/
				    sargsw =exp(-model->CNTbulkJctSideGradingCoeff*
						log(arg));
#ifndef NOSQRT
				}
			    }
#endif /*NOSQRT*/
			    *(ckt->CKTstate0 + here->CNTqbs) =
				here->CNTtBulkPot*(here->CNTCbs*
						    (1-arg*sarg)/(1-model->CNTbulkJctBotGradingCoeff)
						    +here->CNTCbssw*
						    (1-arg*sargsw)/
						    (1-model->CNTbulkJctSideGradingCoeff));
			    here->CNTcapbs=here->CNTCbs*sarg+
                                here->CNTCbssw*sargsw;
			} else {
			    *(ckt->CKTstate0 + here->CNTqbs) = here->CNTf4s +
                                vbs*(here->CNTf2s+vbs*(here->CNTf3s/2));
			    here->CNTcapbs=here->CNTf2s+here->CNTf3s*vbs;
			}
#ifdef CAPZEROBYPASS						
                    } else {
			*(ckt->CKTstate0 + here->CNTqbs) = 0;
                        here->CNTcapbs=0;
                    }
#endif /*CAPZEROBYPASS*/		    
                }
#ifdef CAPBYPASS
                    if(((ckt->CKTmode & (MODEINITPRED | MODEINITTRAN) ) ||
                        fabs(delvbd) >= ckt->CKTreltol * MAX(fabs(vbd),
                        fabs(*(ckt->CKTstate0+here->CNTvbd)))+
                        ckt->CKTvoltTol)|| senflag)
#endif /*CAPBYPASS*/		    	
	
                    /* can't bypass the diode capacitance calculations */
                {
#ifdef CAPZEROBYPASS					
                    if(here->CNTCbd != 0 || here->CNTCbdsw != 0 ) {
#endif /*CAPZEROBYPASS*/		    		    
			if (vbd < here->CNTtDepCap) {
			    arg=1-vbd/here->CNTtBulkPot;
			    /*
			     * the following block looks somewhat long and messy,
			     * but since most users use the default grading
			     * coefficients of .5, and sqrt is MUCH faster than an
			     * exp(log()) we use this special case code to buy time.
			     * (as much as 10% of total job time!)
			     */
#ifndef NOSQRT
			    if(model->CNTbulkJctBotGradingCoeff == .5 &&
			       model->CNTbulkJctSideGradingCoeff == .5) {
				sarg = sargsw = 1/sqrt(arg);
			    } else {
				if(model->CNTbulkJctBotGradingCoeff == .5) {
				    sarg = 1/sqrt(arg);
				} else {
#endif /*NOSQRT*/
				    sarg = exp(-model->CNTbulkJctBotGradingCoeff*
					       log(arg));
#ifndef NOSQRT
				}
				if(model->CNTbulkJctSideGradingCoeff == .5) {
				    sargsw = 1/sqrt(arg);
				} else {
#endif /*NOSQRT*/
				    sargsw =exp(-model->CNTbulkJctSideGradingCoeff*
						log(arg));
#ifndef NOSQRT
				}
			    }
#endif /*NOSQRT*/
			    *(ckt->CKTstate0 + here->CNTqbd) =
				here->CNTtBulkPot*(here->CNTCbd*
						    (1-arg*sarg)
						    /(1-model->CNTbulkJctBotGradingCoeff)
						    +here->CNTCbdsw*
						    (1-arg*sargsw)
						    /(1-model->CNTbulkJctSideGradingCoeff));
			    here->CNTcapbd=here->CNTCbd*sarg+
                                here->CNTCbdsw*sargsw;
			} else {
			    *(ckt->CKTstate0 + here->CNTqbd) = here->CNTf4d +
                                vbd * (here->CNTf2d + vbd * here->CNTf3d/2);
			    here->CNTcapbd=here->CNTf2d + vbd * here->CNTf3d;
			}
#ifdef CAPZEROBYPASS						
		    } else {
			*(ckt->CKTstate0 + here->CNTqbd) = 0;
			here->CNTcapbd = 0;
		    }
#endif /*CAPZEROBYPASS*/		    		    
                }
/*
  
*/

                if(SenCond && (ckt->CKTsenInfo->SENmode==TRANSEN)) goto next2;

                if ( (ckt->CKTmode & MODETRAN) || ( (ckt->CKTmode&MODEINITTRAN)
						    && !(ckt->CKTmode&MODEUIC)) ) {
                    /* (above only excludes tranop, since we're only at this
                     * point if tran or tranop )
                     */

                    /*
                     *    calculate equivalent conductances and currents for
                     *    depletion capacitors
                     */

                    /* integrate the capacitors and save results */

                    error = NIintegrate(ckt,&geq,&ceq,here->CNTcapbd,
					here->CNTqbd);
                    if(error) return(error);
                    here->CNTgbd += geq;
                    here->CNTcbd += *(ckt->CKTstate0 + here->CNTcqbd);
                    here->CNTcd -= *(ckt->CKTstate0 + here->CNTcqbd);
                    error = NIintegrate(ckt,&geq,&ceq,here->CNTcapbs,
					here->CNTqbs);
                    if(error) return(error);
                    here->CNTgbs += geq;
                    here->CNTcbs += *(ckt->CKTstate0 + here->CNTcqbs);
                }
            }
/*
  
*/

            if(SenCond) goto next2;


            /*
             *  check convergence
             */
            if ( (here->CNToff == 0)  || 
		 (!(ckt->CKTmode & (MODEINITFIX|MODEINITSMSIG))) ){
                if (Check == 1) {
                    ckt->CKTnoncon++;
		    ckt->CKTtroubleElt = (GENinstance *) here;
                }
            }
/*
  
*/

            /* save things away for next time */

	next2:      *(ckt->CKTstate0 + here->CNTvbs) = vbs;
            *(ckt->CKTstate0 + here->CNTvbd) = vbd;
            *(ckt->CKTstate0 + here->CNTvgs) = vgs;
            *(ckt->CKTstate0 + here->CNTvds) = vds;

/*
  
*/

            /*
             *     meyer's capacitor model
             */
            if ( ckt->CKTmode & (MODETRAN | MODETRANOP | MODEINITSMSIG) ) {
                /*
                 *     calculate meyer's capacitors
                 */
                /* 
                 * new cmeyer - this just evaluates at the current time,
                 * expects you to remember values from previous time
                 * returns 1/2 of non-constant portion of capacitance
                 * you must add in the other half from previous time
                 * and the constant part
                 */
                if (here->CNTmode > 0){
                    DEVqmeyer (vgs,vgd,vgb,von,vdsat,
			       (ckt->CKTstate0 + here->CNTcapgs),
			       (ckt->CKTstate0 + here->CNTcapgd),
			       (ckt->CKTstate0 + here->CNTcapgb),
			       here->CNTtPhi,OxideCap);
                } else {
                    DEVqmeyer (vgd,vgs,vgb,von,vdsat,
			       (ckt->CKTstate0 + here->CNTcapgd),
			       (ckt->CKTstate0 + here->CNTcapgs),
			       (ckt->CKTstate0 + here->CNTcapgb),
			       here->CNTtPhi,OxideCap);
                }
                vgs1 = *(ckt->CKTstate1 + here->CNTvgs);
                vgd1 = vgs1 - *(ckt->CKTstate1 + here->CNTvds);
                vgb1 = vgs1 - *(ckt->CKTstate1 + here->CNTvbs);
                if(ckt->CKTmode & (MODETRANOP|MODEINITSMSIG)) {
                    capgs =  2 * *(ckt->CKTstate0+here->CNTcapgs)+ 
			GateSourceOverlapCap ;
                    capgd =  2 * *(ckt->CKTstate0+here->CNTcapgd)+ 
			GateDrainOverlapCap ;
                    capgb =  2 * *(ckt->CKTstate0+here->CNTcapgb)+ 
			GateBulkOverlapCap ;
                } else {
                    capgs = ( *(ckt->CKTstate0+here->CNTcapgs)+ 
                              *(ckt->CKTstate1+here->CNTcapgs) +
                              GateSourceOverlapCap );
                    capgd = ( *(ckt->CKTstate0+here->CNTcapgd)+ 
                              *(ckt->CKTstate1+here->CNTcapgd) +
                              GateDrainOverlapCap );
                    capgb = ( *(ckt->CKTstate0+here->CNTcapgb)+ 
                              *(ckt->CKTstate1+here->CNTcapgb) +
                              GateBulkOverlapCap );
                }
                if(ckt->CKTsenInfo){
                    here->CNTcgs = capgs;
                    here->CNTcgd = capgd;
                    here->CNTcgb = capgb;
                }
/*
  
*/

                /*
                 *     store small-signal parameters (for meyer's model)
                 *  all parameters already stored, so done...
                 */
                if(SenCond){
                    if((ckt->CKTsenInfo->SENmode == DCSEN)||
		       (ckt->CKTsenInfo->SENmode == ACSEN)){
                        continue;
                    }
                }

#ifndef PREDICTOR
                if (ckt->CKTmode & (MODEINITPRED | MODEINITTRAN) ) {
                    *(ckt->CKTstate0 + here->CNTqgs) =
                        (1+xfact) * *(ckt->CKTstate1 + here->CNTqgs)
                        - xfact * *(ckt->CKTstate2 + here->CNTqgs);
                    *(ckt->CKTstate0 + here->CNTqgd) =
                        (1+xfact) * *(ckt->CKTstate1 + here->CNTqgd)
                        - xfact * *(ckt->CKTstate2 + here->CNTqgd);
                    *(ckt->CKTstate0 + here->CNTqgb) =
                        (1+xfact) * *(ckt->CKTstate1 + here->CNTqgb)
                        - xfact * *(ckt->CKTstate2 + here->CNTqgb);
                } else {
#endif /*PREDICTOR*/
                    if(ckt->CKTmode & MODETRAN) {
                        *(ckt->CKTstate0 + here->CNTqgs) = (vgs-vgs1)*capgs +
                            *(ckt->CKTstate1 + here->CNTqgs) ;
                        *(ckt->CKTstate0 + here->CNTqgd) = (vgd-vgd1)*capgd +
                            *(ckt->CKTstate1 + here->CNTqgd) ;
                        *(ckt->CKTstate0 + here->CNTqgb) = (vgb-vgb1)*capgb +
                            *(ckt->CKTstate1 + here->CNTqgb) ;
                    } else {
                        /* TRANOP only */
                        *(ckt->CKTstate0 + here->CNTqgs) = vgs*capgs;
                        *(ckt->CKTstate0 + here->CNTqgd) = vgd*capgd;
                        *(ckt->CKTstate0 + here->CNTqgb) = vgb*capgb;
                    }
#ifndef PREDICTOR
                }
#endif /*PREDICTOR*/
            }
	bypass:
            if(SenCond) continue;

            if ( (ckt->CKTmode & (MODEINITTRAN)) || 
		 (! (ckt->CKTmode & (MODETRAN)) )  ) {
                /*
                 *  initialize to zero charge conductances 
                 *  and current
                 */
                gcgs=0;
                ceqgs=0;
                gcgd=0;
                ceqgd=0;
                gcgb=0;
                ceqgb=0;
            } else {
                if(capgs == 0) *(ckt->CKTstate0 + here->CNTcqgs) =0;
                if(capgd == 0) *(ckt->CKTstate0 + here->CNTcqgd) =0;
                if(capgb == 0) *(ckt->CKTstate0 + here->CNTcqgb) =0;
                /*
                 *    calculate equivalent conductances and currents for
                 *    meyer"s capacitors
                 */
                error = NIintegrate(ckt,&gcgs,&ceqgs,capgs,here->CNTqgs);
                if(error) return(error);
                error = NIintegrate(ckt,&gcgd,&ceqgd,capgd,here->CNTqgd);
                if(error) return(error);
                error = NIintegrate(ckt,&gcgb,&ceqgb,capgb,here->CNTqgb);
                if(error) return(error);
                ceqgs=ceqgs-gcgs*vgs+ckt->CKTag[0]* 
		    *(ckt->CKTstate0 + here->CNTqgs);
                ceqgd=ceqgd-gcgd*vgd+ckt->CKTag[0]*
		    *(ckt->CKTstate0 + here->CNTqgd);
                ceqgb=ceqgb-gcgb*vgb+ckt->CKTag[0]*
		    *(ckt->CKTstate0 + here->CNTqgb);
            }
            /*
             *     store charge storage info for meyer's cap in lx table
             */

            /*
             *  load current vector
             */
            ceqbs = model->CNTtype * 
		(here->CNTcbs-(here->CNTgbs)*vbs);
            ceqbd = model->CNTtype * 
		(here->CNTcbd-(here->CNTgbd)*vbd);
            if (here->CNTmode >= 0) {
                xnrm=1;
                xrev=0;
                cdreq=model->CNTtype*(cdrain-here->CNTgds*vds-
				       here->CNTgm*vgs-here->CNTgmbs*vbs);
            } else {
                xnrm=0;
                xrev=1;
                cdreq = -(model->CNTtype)*(cdrain-here->CNTgds*(-vds)-
					    here->CNTgm*vgd-here->CNTgmbs*vbd);
            }
	//if (Efi<=0.0){
            *(ckt->CKTrhs + here->CNTgNode) -= 
                (model->CNTtype * (ceqgs + ceqgb + ceqgd));
            *(ckt->CKTrhs + here->CNTbNode) -=
                (ceqbs + ceqbd - model->CNTtype * ceqgb);
            *(ckt->CKTrhs + here->CNTdNodePrime) +=
		(ceqbd - cdreq + model->CNTtype * ceqgd);
            *(ckt->CKTrhs + here->CNTsNodePrime) += 
		(cdreq + ceqbs + model->CNTtype * ceqgs);
		//}
	//else {
            //*(ckt->CKTrhs + here->CNTgNode) -= 
               // -(model->CNTtype * (ceqgs + ceqgb + ceqgd));
           // *(ckt->CKTrhs + here->CNTbNode) -=
          //      -(ceqbs + ceqbd - model->CNTtype * ceqgb);
           // *(ckt->CKTrhs + here->CNTdNodePrime) +=
	//	-(ceqbd - cdreq + model->CNTtype * ceqgd);
           // *(ckt->CKTrhs + here->CNTsNodePrime) += 
	//	-(cdreq + ceqbs + model->CNTtype * ceqgs);
		//}
printf("%lf %lf %lf\n",cdreq, ceqbs, ceqbd);
            /*
             *  load y matrix
             */

            *(here->CNTDdPtr) += (here->CNTdrainConductance);
            *(here->CNTGgPtr) += ((gcgd+gcgs+gcgb));
            *(here->CNTSsPtr) += (here->CNTsourceConductance);
            *(here->CNTBbPtr) += (here->CNTgbd+here->CNTgbs+gcgb);
            *(here->CNTDPdpPtr) += 
		(here->CNTdrainConductance+here->CNTgds+
		 here->CNTgbd+xrev*(here->CNTgm+here->CNTgmbs)+gcgd);
            *(here->CNTSPspPtr) += 
		(here->CNTsourceConductance+here->CNTgds+
		 here->CNTgbs+xnrm*(here->CNTgm+here->CNTgmbs)+gcgs);
            *(here->CNTDdpPtr) += (-here->CNTdrainConductance);
            *(here->CNTGbPtr) -= gcgb;
            *(here->CNTGdpPtr) -= gcgd;
            *(here->CNTGspPtr) -= gcgs;
            *(here->CNTSspPtr) += (-here->CNTsourceConductance);
            *(here->CNTBgPtr) -= gcgb;
            *(here->CNTBdpPtr) -= here->CNTgbd;
            *(here->CNTBspPtr) -= here->CNTgbs;
            *(here->CNTDPdPtr) += (-here->CNTdrainConductance);
            *(here->CNTDPgPtr) += ((xnrm-xrev)*here->CNTgm-gcgd);
            *(here->CNTDPbPtr) += (-here->CNTgbd+(xnrm-xrev)*here->CNTgmbs);
            *(here->CNTDPspPtr) += (-here->CNTgds-xnrm*
				     (here->CNTgm+here->CNTgmbs));
            *(here->CNTSPgPtr) += (-(xnrm-xrev)*here->CNTgm-gcgs);
            *(here->CNTSPsPtr) += (-here->CNTsourceConductance);
            *(here->CNTSPbPtr) += (-here->CNTgbs-(xnrm-xrev)*here->CNTgmbs);
            *(here->CNTSPdpPtr) += (-here->CNTgds-xrev*
				     (here->CNTgm+here->CNTgmbs));
        }
    }
    return(OK);
}

int
Range1(double Vds1,
    double q1,
    double c1,
    double qC1,
    double qCN01,
    double xmax1,
    double xmin1,
    double inte1,
    int n1)
{
    int y;
    int i;
    int j;

for(i=1;i<=n1;i++)
	{
	if(Root(Vds1,q1,c1,qC1,qCN01,n1+1,n1+1)<=(xmin1-Vds1))
		{
		y=n1+1;
		goto next3;
		}
	else if(Root(Vds1,q1,c1,qC1,qCN01,n1,n1)>xmax1)
		{
		y=n1;
		goto next3;
		}
	else if(Vds1>=(n1-i)*inte1)
		{
		if(i==1)
			{
			if((xmax1-Vds1)<Root(Vds1,q1,c1,qC1,qCN01,n1+1,n1) && Root(Vds1,q1,c1,qC1,qCN01,n1+1,n1)<=xmin1)
				{
				y=n1+1;
				goto next3;
				}
			else
				{
				for(j=1;j<=n1-1;j++)
					{
					if((xmin1-Vds1+(j-1)*inte1)<Root(Vds1,q1,c1,qC1,qCN01,n1+1,j) && Root(Vds1,q1,c1,qC1,qCN01,n1+1,j)<=(xmin1-Vds1+j*inte1))
						{
						y=n1+1;
						goto next3;
						}
					else if((xmin1+(j-1)*inte1)<Root(Vds1,q1,c1,qC1,qCN01,j,n1) && Root(Vds1,q1,c1,qC1,qCN01,j,n1)<=(xmin1+j*inte1))
						{
						y=j;
						goto next3;
						}
					}
				}
			}
		else if(i>1)
			{
			if((xmin1-Vds1+(n1-i)*inte1)<Root(Vds1,q1,c1,qC1,qCN01,n1+1,n1-i+1) && Root(Vds1,q1,c1,qC1,qCN01,n1+1,n1-i+1)<=xmin1)
				{
				y=n1+1;
				goto next3;
				}
			else
				{
				for(j=1;j<=i-1;j++)
					{
					if((xmin1+(j-1)*inte1)<Root(Vds1,q1,c1,qC1,qCN01,j,n1-i+j) && Root(Vds1,q1,c1,qC1,qCN01,j,n1-i+j)<=(xmin1-Vds1+(n1-i+j)*inte1))
						{
						y=j;
						goto next3;
						}
					else if((xmin1-Vds1+(n1-i+j)*inte1)<Root(Vds1,q1,c1,qC1,qCN01,j,n1-i+j+1) && Root(Vds1,q1,c1,qC1,qCN01,j,n1-i+j+1)<=(xmin1+j*inte1))
						{
						y=j;
						goto next3;
						}
					}
				if(i<n1)
					{
					for(j=1;j<n1-i;j++)
						{
						if((xmin1-Vds1+(j-1)*inte1)<Root(Vds1,q1,c1,qC1,qCN01,n1+1,j) && Root(Vds1,q1,c1,qC1,qCN01,n1+1,j)<=(xmin1-Vds1+j*inte1))
							{
							y=n1+1;
							goto next3;
							}
						else if((xmin1+(i+j-2)*inte1)<Root(Vds1,q1,c1,qC1,qCN01,i+j-1,n1) && Root(Vds1,q1,c1,qC1,qCN01,i+j-1,n1)<=(xmin1+(i+j-1)*inte1));
							{
							y=i+j-1;
							goto next3;
							}
						}
					}
				}
			}
		}
	else{y=4;}
	}

next3:	return y;

}

int
Range2(double Vds2,
    double q2,
    double c2,
    double qC2,
    double qCN02,
    double xmax2,
    double xmin2,
    double inte2,
    int n2)
{
    int y;
    int i;
    int j;

for(i=1;i<=n2;i++)
	{
	if(Root(Vds2,q2,c2,qC2,qCN02,n2+1,n2+1)<=(xmin2-Vds2))
		{
		y=n2+1;
		goto next4;
		}
	else if(Root(Vds2,q2,c2,qC2,qCN02,n2,n2)>xmax2)
		{
		y=n2;
		goto next4;
		}
	else if(Vds2>=(n2-i)*inte2)
		{
		if(i==1)
			{
			if((xmax2-Vds2)<Root(Vds2,q2,c2,qC2,qCN02,n2+1,n2) && Root(Vds2,q2,c2,qC2,qCN02,n2+1,n2)<=xmin2)
				{
				y=n2;
				goto next4;
				}
			else
				{
				for(j=1;j<=n2-1;j++)
					{
					if((xmin2-Vds2+(j-1)*inte2)<Root(Vds2,q2,c2,qC2,qCN02,n2+1,j) && Root(Vds2,q2,c2,qC2,qCN02,n2+1,j)<=(xmin2-Vds2+j*inte2))
						{
						y=j;
						goto next4;
						}
					else if((xmin2+(j-1)*inte2)<Root(Vds2,q2,c2,qC2,qCN02,j,n2) && Root(Vds2,q2,c2,qC2,qCN02,j,n2)<=(xmin2+j*inte2))
						{
						y=n2;
						goto next4;
						}
					}
				}
			}
		else if(i>1)
			{
			if((xmin2-Vds2+(n2-i)*inte2)<Root(Vds2,q2,c2,qC2,qCN02,n2+1,n2-i+1) && Root(Vds2,q2,c2,qC2,qCN02,n2+1,n2-i+1)<=xmin2)
				{
				y=n2-i+1;
				goto next4;
				}
			else
				{
				for(j=1;j<=i-1;j++)
					{
					if((xmin2+(j-1)*inte2)<Root(Vds2,q2,c2,qC2,qCN02,j,n2-i+j) && Root(Vds2,q2,c2,qC2,qCN02,j,n2-i+j)<=(xmin2-Vds2+(n2-i+j)*inte2))
						{
						y=n2-i+j;
						goto next4;
						}
					else if((xmin2-Vds2+(n2-i+j)*inte2)<Root(Vds2,q2,c2,qC2,qCN02,j,n2-i+j+1) && Root(Vds2,q2,c2,qC2,qCN02,j,n2-i+j+1)<=(xmin2+j*inte2))
						{
						y=n2-i+j+1;
						goto next4;
						}
					}
				if(i<n2)
					{
					for(j=1;j<n2-i;j++)
						{
						if((xmin2-Vds2+(j-1)*inte2)<Root(Vds2,q2,c2,qC2,qCN02,n2+1,j) && Root(Vds2,q2,c2,qC2,qCN02,n2+1,j)<=(xmin2-Vds2+j*inte2))
							{
							y=j;
							goto next4;
							}
						else if((xmin2+(i+j-2)*inte2)<Root(Vds2,q2,c2,qC2,qCN02,i+j-1,n2) && Root(Vds2,q2,c2,qC2,qCN02,i+j-1,n2)<=(xmin2+(i+j-1)*inte2));
							{
							y=n2;
							goto next4;
							}
						}
					}
				}
			}
		}
	else{y=4;}
	}

next4:	return y;

}

double
Root(double Vds,
    double q,
    double c,
    double qC,
    double qCN0,
    int coeS,
    int coeD)
{
    
    double aa;
    double bb;
    double cc;
    double dd;
    double p;
    double qq;
    double r;
    double u;
    double v;
    double g;
    double hh;
    double tmp;
    double fai;
    double pi=3.1415926;
    double realy1;
    double realy2;
    double realy3;
    double imagy1;
    double imagy2;
    double imagy3;
    double rt;

aa=qC*(Coeff(coeS,1)+Coeff(coeD,1));
bb=qC*(Coeff(coeD,1)*3.0*Vds+Coeff(coeS,2)+Coeff(coeD,2));
cc=qC*(Coeff(coeD,1)*3.0*Vds*Vds+Coeff(coeD,2)*2.0*Vds+Coeff(coeS,3)+Coeff(coeD,3))-q;
dd=qC*(Coeff(coeD,1)*Vds*Vds*Vds+Coeff(coeD,2)*Vds*Vds+Coeff(coeD,3)*Vds+Coeff(coeS,4)+Coeff(coeD,4))+c-qCN0;
//printf("abcd %.20lf %.20lf %.20lf %.20lf\n",aa,bb,cc,dd);
if(aa==0.0)
	{
	if(bb==0.0)
		{
		rt=-dd/cc;
		}
	else
		{
		rt=(-cc+sqrt(cc*cc-4.0*bb*dd))/2.0*bb;
		}
	}
else
	{
	p=(3.0*aa*cc-bb*bb)/3.0/aa/aa;
	qq=(2.0*bb*bb*bb-9.0*aa*bb*cc+27.0*aa*aa*dd)/(27*aa*aa*aa);
	r=bb/3.0/aa;
	hh=(qq/2.0)*(qq/2.0)+(p/3.0)*(p/3.0)*(p/3.0);
//printf("pqrh %.20lf %.20lf %.20lf %.20lf\n",p,qq,r,hh);

	if(hh>=0.0)
		{
		g=sqrt(hh);
		
		if(-qq/2.0+g<0.0)
			{
			tmp=fabs(-qq/2.0+g);
			u=-pow(tmp,0.3333);
			}
		else
			{
			tmp=-qq/2.0+g;
			u=pow(tmp,0.3333);
			}

		if(-qq/2.0-g<0.0)
			v=-pow((fabs(-qq/2.0-g)),0.3333);
		else
			v=pow((-qq/2.0-g),0.3333);
//printf("uvtmp %.20lf %.20lf %.20lf\n",u,v,tmp);

		if(hh==0.0)
			{
			realy1=u+v-r;
			realy2=-(u+v)/2.0-r;
			realy3=-(u+v)/2.0-r;
			imagy1=0.0;
			imagy2=0.0;
			imagy3=0.0;
			}
		else
			{
			realy1=u+v-r;
			realy2=-(u+v)/2.0;
			realy3=-(u+v)/2.0;
			imagy1=0.0;
			imagy2=sqrt(3.0)*(u-v)/2.0;
			imagy3=-sqrt(3.0)*(u-v)/2.0;
			}
		}
	else
		{
		fai=acos((-qq/2.0)/(sqrt(pow(fabs(p),3.0)/27.0)));
		realy1=2.0*sqrt(fabs(p)/3.0)*cos(fai/3.0)-r;

/* FRA - Original Code */
//		realy1=-2.0*sqrt(fabs(p)/3.0)*cos((fai+pi)/3.0)-r;
//		realy1=-2.0*sqrt(fabs(p)/3.0)*cos((fai+pi)/3.0)-r;
/***********************/
/* FRA - New Code */
		realy2=-2.0*sqrt(fabs(p)/3.0)*cos((fai+pi)/3.0)-r;
		realy3=-2.0*sqrt(fabs(p)/3.0)*cos((fai+pi)/3.0)-r;
/***********************/

		imagy1=0.0;
		imagy2=0.0;
		imagy3=0.0;
		}
		//printf("r123 %lf %lf %lf %lf\n",fai,realy1,realy2,realy3);
	if(-0.5<=realy3 && realy3<=-0.2)
		rt=realy3;
	else if(-0.5<=realy2 && realy2<=-0.2)
		rt=realy2;
	else if(-0.5<=realy1 && realy1<=-0.2)
		rt=realy1;
	else
		rt=0.0;
	}
//printf("rt:%lf\n",rt);
return rt;

}



double
Coeff(int cor,int coc)
{
double co;

if(coc==1)
	co=CSa[cor-1];
else if(coc==2)
	co=CSb[cor-1];
else if(coc==3)
	co=CSc[cor-1];
else if(coc==4)
	co=CSd[cor-1];
else
	co=0.0;

return co;

}

double
Ncnt(double ND0,double NEG,double NKT,double NEF)
{
double n0;
double En;
double zn;
double z0;
double z1=0.0;
double f0;
double fn=0.0;
int i;

if(NEF<2.0*NKT)
	En=10.0*NKT+NEG/2.0;
else
	En=(8.0*NKT+NEF)+NEG/2.0;
zn=sqrt(En*En-(NEG/2.0)*(NEG/2.0));
z0=zn/1000.0;
for(i=0;i<=1000;i++)
	{
	f0=1.0/(1.0+exp((sqrt(z1*z1+(NEG/2.0)*(NEG/2.0))-NEG/2.0-NEF)/NKT));
	fn=fn+f0;
	z1=z1+z0;
	}

n0=ND0*fn*z0;

return n0;

}
