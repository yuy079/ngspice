/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/

#ifndef CNT
#define CNT

#include "ngspice/ifsim.h"
#include "ngspice/cktdefs.h"
#include "ngspice/gendefs.h"
#include "ngspice/complex.h"
#include "ngspice/noisedef.h"

/* declarations for level 1 CNTFETs */

/* information needed for each instance */

typedef struct sCNTinstance {
    struct sCNTmodel *sCNTmodPtr; /* backpointer to model */
    struct sCNTinstance *CNTnextInstance;  /* pointer to next instance of
                                              *current model*/
    IFuid CNTname; /* pointer to character string naming this instance */
    int CNTowner;  /* number of owner process */
    int CNTstates;     /* index into state table for this device */

    int CNTdNode;  /* number of the gate node of the cntfet */
    int CNTgNode;  /* number of the gate node of the cntfet */
    int CNTsNode;  /* number of the source node of the cntfet */
    int CNTbNode;  /* number of the bulk node of the cntfet */
    int CNTdNodePrime; /* number of the internal drain node of the cntfet */
    int CNTsNodePrime; /* number of the internal source node of the cntfet */
    
    double CNTm;   /* parallel device multiplier */

    double CNTl;   /* the length of the channel region */
    double CNTw;   /* the width of the channel region */
    double CNTdia;   /* the diameter of the cnt channel */
    double CNTef;   /* the fermi-level of the cnt */
    double CNTdrainArea;   /* the area of the drain diffusion */
    double CNTsourceArea;  /* the area of the source diffusion */
    double CNTdrainSquares;    /* the length of the drain in squares */
    double CNTsourceSquares;   /* the length of the source in squares */
    double CNTdrainPerimiter;
    double CNTsourcePerimiter;
    double CNTsourceConductance;   /*conductance of source(or 0):set in setup*/
    double CNTdrainConductance;    /*conductance of drain(or 0):set in setup*/
    double CNTtemp;    /* operating temperature of this instance */
    double CNTdtemp;   /* operating temperature of the instance relative to circuit temperature*/

    double CNTtTransconductance;   /* temperature corrected transconductance*/
    double CNTtSurfMob;            /* temperature corrected surface mobility */
    double CNTtPhi;                /* temperature corrected Phi */
    double CNTtVto;                /* temperature corrected Vto */
    double CNTtSatCur;             /* temperature corrected saturation Cur. */
    double CNTtSatCurDens; /* temperature corrected saturation Cur. density*/
    double CNTtCbd;                /* temperature corrected B-D Capacitance */
    double CNTtCbs;                /* temperature corrected B-S Capacitance */
    double CNTtCj;         /* temperature corrected Bulk bottom Capacitance */
    double CNTtCjsw;       /* temperature corrected Bulk side Capacitance */
    double CNTtBulkPot;    /* temperature corrected Bulk potential */
    double CNTtDepCap;     /* temperature adjusted transition point in */
                            /* the cureve matching Fc * Vj */
    double CNTtVbi;        /* temperature adjusted Vbi */

    double CNTicVBS;   /* initial condition B-S voltage */
    double CNTicVDS;   /* initial condition D-S voltage */
    double CNTicVGS;   /* initial condition G-S voltage */
    double CNTvon;
    double CNTvdsat;
    double CNTsourceVcrit; /* Vcrit for pos. vds */
    double CNTdrainVcrit;  /* Vcrit for pos. vds */
    double CNTcd;
    double CNTcbs;
    double CNTcbd;
    double CNTgmbs;
    double CNTgm;
    double CNTgds;
    double CNTgbd;
    double CNTgbs;
    double CNTcapbd;
    double CNTcapbs;
    double CNTCbd;
    double CNTCbdsw;
    double CNTCbs;
    double CNTCbssw;
    double CNTf2d;
    double CNTf3d;
    double CNTf4d;
    double CNTf2s;
    double CNTf3s;
    double CNTf4s;

/*
 * naming convention:
 * x = vgs
 * y = vbs
 * z = vds
 * cdr = cdrain
 */

#define	CNTNDCOEFFS	30

#ifndef NODISTO
	double CNTdCoeffs[CNTNDCOEFFS];
#else /* NODISTO */
	double *CNTdCoeffs;
#endif /* NODISTO */

#ifndef CONFIG

#define	capbs2		CNTdCoeffs[0]
#define	capbs3		CNTdCoeffs[1]
#define	capbd2		CNTdCoeffs[2]
#define	capbd3		CNTdCoeffs[3]
#define	gbs2		CNTdCoeffs[4]
#define	gbs3		CNTdCoeffs[5]
#define	gbd2		CNTdCoeffs[6]
#define	gbd3		CNTdCoeffs[7]
#define	capgb2		CNTdCoeffs[8]
#define	capgb3		CNTdCoeffs[9]
#define	cdr_x2		CNTdCoeffs[10]
#define	cdr_y2		CNTdCoeffs[11]
#define	cdr_z2		CNTdCoeffs[12]
#define	cdr_xy		CNTdCoeffs[13]
#define	cdr_yz		CNTdCoeffs[14]
#define	cdr_xz		CNTdCoeffs[15]
#define	cdr_x3		CNTdCoeffs[16]
#define	cdr_y3		CNTdCoeffs[17]
#define	cdr_z3		CNTdCoeffs[18]
#define	cdr_x2z		CNTdCoeffs[19]
#define	cdr_x2y		CNTdCoeffs[20]
#define	cdr_y2z		CNTdCoeffs[21]
#define	cdr_xy2		CNTdCoeffs[22]
#define	cdr_xz2		CNTdCoeffs[23]
#define	cdr_yz2		CNTdCoeffs[24]
#define	cdr_xyz		CNTdCoeffs[25]
#define	capgs2		CNTdCoeffs[26]
#define	capgs3		CNTdCoeffs[27]
#define	capgd2		CNTdCoeffs[28]
#define	capgd3		CNTdCoeffs[29]

#endif

#define CNTRDNOIZ	0
#define CNTRSNOIZ   1
#define CNTIDNOIZ       2
#define CNTFLNOIZ 3
#define CNTTOTNOIZ    4

#define CNTNSRCS     5     /* the number of CNTFET noise sources*/

#ifndef NONOISE
    double CNTnVar[NSTATVARS][CNTNSRCS];
#else /* NONOISE */
	double **CNTnVar;
#endif /* NONOISE */

    int CNTmode;       /* device mode : 1 = normal, -1 = inverse */


/* FRA - Original Code */
//    unsigned CNToff:1;  /* non-zero to indicate device is off for dc analysis*/
/***********************/
/* FRA - New Code */
    int CNToff;  /* non-zero to indicate device is off for dc analysis*/
/******************/

    unsigned CNTtempGiven :1;  /* instance temperature specified */
    unsigned CNTdtempGiven :1;  /* instance delta temperature specified */
    unsigned CNTmGiven :1;
    unsigned CNTlGiven :1;
    unsigned CNTwGiven :1;
    unsigned CNTdiaGiven :1;  /*cnt diameter*/
    unsigned CNTefGiven :1;  /*cnt fermi-level*/
    unsigned CNTdrainAreaGiven :1;
    unsigned CNTsourceAreaGiven    :1;
    unsigned CNTdrainSquaresGiven  :1;
    unsigned CNTsourceSquaresGiven :1;
    unsigned CNTdrainPerimiterGiven    :1;
    unsigned CNTsourcePerimiterGiven   :1;
    unsigned CNTdNodePrimeSet  :1;
    unsigned CNTsNodePrimeSet  :1;
    unsigned CNTicVBSGiven :1;
    unsigned CNTicVDSGiven :1;
    unsigned CNTicVGSGiven :1;
    unsigned CNTvonGiven   :1;
    unsigned CNTvdsatGiven :1;
    unsigned CNTmodeGiven  :1;


    double *CNTDdPtr;      /* pointer to sparse matrix element at
                                     * (Drain node,drain node) */
    double *CNTGgPtr;      /* pointer to sparse matrix element at
                                     * (gate node,gate node) */
    double *CNTSsPtr;      /* pointer to sparse matrix element at
                                     * (source node,source node) */
    double *CNTBbPtr;      /* pointer to sparse matrix element at
                                     * (bulk node,bulk node) */
    double *CNTDPdpPtr;    /* pointer to sparse matrix element at
                                     * (drain prime node,drain prime node) */
    double *CNTSPspPtr;    /* pointer to sparse matrix element at
                                     * (source prime node,source prime node) */
    double *CNTDdpPtr;     /* pointer to sparse matrix element at
                                     * (drain node,drain prime node) */
    double *CNTGbPtr;      /* pointer to sparse matrix element at
                                     * (gate node,bulk node) */
    double *CNTGdpPtr;     /* pointer to sparse matrix element at
                                     * (gate node,drain prime node) */
    double *CNTGspPtr;     /* pointer to sparse matrix element at
                                     * (gate node,source prime node) */
    double *CNTSspPtr;     /* pointer to sparse matrix element at
                                     * (source node,source prime node) */
    double *CNTBdpPtr;     /* pointer to sparse matrix element at
                                     * (bulk node,drain prime node) */
    double *CNTBspPtr;     /* pointer to sparse matrix element at
                                     * (bulk node,source prime node) */
    double *CNTDPspPtr;    /* pointer to sparse matrix element at
                                     * (drain prime node,source prime node) */
    double *CNTDPdPtr;     /* pointer to sparse matrix element at
                                     * (drain prime node,drain node) */
    double *CNTBgPtr;      /* pointer to sparse matrix element at
                                     * (bulk node,gate node) */
    double *CNTDPgPtr;     /* pointer to sparse matrix element at
                                     * (drain prime node,gate node) */

    double *CNTSPgPtr;     /* pointer to sparse matrix element at
                                     * (source prime node,gate node) */
    double *CNTSPsPtr;     /* pointer to sparse matrix element at
                                     * (source prime node,source node) */
    double *CNTDPbPtr;     /* pointer to sparse matrix element at
                                     * (drain prime node,bulk node) */
    double *CNTSPbPtr;     /* pointer to sparse matrix element at
                                     * (source prime node,bulk node) */
    double *CNTSPdpPtr;    /* pointer to sparse matrix element at
                                     * (source prime node,drain prime node) */

    int  CNTsenParmNo;   /* parameter # for sensitivity use;
            set equal to 0  if  neither length
            nor width of the cntfet is a design
            parameter */
    unsigned CNTsens_l :1;   /* field which indicates whether  
            length of the cntfet is a design
            parameter or not */
    unsigned CNTsens_w :1;   /* field which indicates whether  
            width of the cntfet is a design
            parameter or not */
    unsigned CNTsenPertFlag :1; /* indictes whether the the 
            parameter of the particular instance is 
            to be perturbed */
    double CNTcgs;
    double CNTcgd;
    double CNTcgb;
    double *CNTsens;

#define CNTsenCgs CNTsens /* contains pertured values of cgs */
#define CNTsenCgd CNTsens + 6 /* contains perturbed values of cgd*/
#define CNTsenCgb CNTsens + 12 /* contains perturbed values of cgb*/
#define CNTsenCbd CNTsens + 18 /* contains perturbed values of cbd*/
#define CNTsenCbs CNTsens + 24 /* contains perturbed values of cbs*/
#define CNTsenGds CNTsens + 30 /* contains perturbed values of gds*/
#define CNTsenGbs CNTsens + 36 /* contains perturbed values of gbs*/
#define CNTsenGbd CNTsens + 42 /* contains perturbed values of gbd*/
#define CNTsenGm CNTsens + 48 /* contains perturbed values of gm*/
#define CNTsenGmbs CNTsens + 54 /* contains perturbed values of gmbs*/
#define CNTdphigs_dl CNTsens + 60
#define CNTdphigd_dl CNTsens + 61
#define CNTdphigb_dl CNTsens + 62
#define CNTdphibs_dl CNTsens + 63
#define CNTdphibd_dl CNTsens + 64
#define CNTdphigs_dw CNTsens + 65
#define CNTdphigd_dw CNTsens + 66
#define CNTdphigb_dw CNTsens + 67
#define CNTdphibs_dw CNTsens + 68
#define CNTdphibd_dw CNTsens + 69

} CNTinstance ;

#define CNTvbd CNTstates+ 0   /* bulk-drain voltage */
#define CNTvbs CNTstates+ 1   /* bulk-source voltage */
#define CNTvgs CNTstates+ 2   /* gate-source voltage */
#define CNTvds CNTstates+ 3   /* drain-source voltage */

#define CNTcapgs CNTstates+4  /* gate-source capacitor value */
#define CNTqgs CNTstates+ 5   /* gate-source capacitor charge */
#define CNTcqgs CNTstates+ 6  /* gate-source capacitor current */

#define CNTcapgd CNTstates+ 7 /* gate-drain capacitor value */
#define CNTqgd CNTstates+ 8   /* gate-drain capacitor charge */
#define CNTcqgd CNTstates+ 9  /* gate-drain capacitor current */

#define CNTcapgb CNTstates+10 /* gate-bulk capacitor value */
#define CNTqgb CNTstates+ 11  /* gate-bulk capacitor charge */
#define CNTcqgb CNTstates+ 12 /* gate-bulk capacitor current */

#define CNTqbd CNTstates+ 13  /* bulk-drain capacitor charge */
#define CNTcqbd CNTstates+ 14 /* bulk-drain capacitor current */

#define CNTqbs CNTstates+ 15  /* bulk-source capacitor charge */
#define CNTcqbs CNTstates+ 16 /* bulk-source capacitor current */

#define CNTnumStates 17

#define CNTsensxpgs CNTstates+17 /* charge sensitivities and 
          their derivatives.  +18 for the derivatives:
          pointer to the beginning of the array */
#define CNTsensxpgd  CNTstates+19
#define CNTsensxpgb  CNTstates+21
#define CNTsensxpbs  CNTstates+23
#define CNTsensxpbd  CNTstates+25

#define CNTnumSenStates 10


/* per model data */

    /* NOTE:  parameters marked 'input - use xxxx' are paramters for
     * which a temperature correction is applied in CNTtemp, thus
     * the CNTxxxx value in the per-instance structure should be used
     * instead in all calculations 
     */


typedef struct sCNTmodel {       /* model structure for a resistor */
    int CNTmodType;    /* type index to this device type */
    struct sCNTmodel *CNTnextModel;    /* pointer to next possible model 
                                          *in linked list */
    CNTinstance * CNTinstances; /* pointer to list of instances 
                                   * that have this model */
    IFuid CNTmodName;       /* pointer to character string naming this model */
    int CNTtype;       /* device type : 1 = ncnt,  -1 = pcnt */
    double CNTtnom;        /* temperature at which parameters measured */
    double CNTlatDiff;
    double CNTjctSatCurDensity;    /* input - use tSatCurDens */
    double CNTjctSatCur;   /* input - use tSatCur */
    double CNTdrainResistance;
    double CNTsourceResistance;
    double CNTsheetResistance;
    double CNTtransconductance;    /* input - use tTransconductance */
    double CNTgateSourceOverlapCapFactor;
    double CNTgateDrainOverlapCapFactor;
    double CNTgateBulkOverlapCapFactor;
    double CNToxideCapFactor;
    double CNTvt0; /* input - use tVto */
    double CNTcapBD;   /* input - use tCbd */
    double CNTcapBS;   /* input - use tCbs */
    double CNTbulkCapFactor;   /* input - use tCj */
    double CNTsideWallCapFactor;   /* input - use tCjsw */
    double CNTbulkJctPotential;    /* input - use tBulkPot */
    double CNTbulkJctBotGradingCoeff;
    double CNTbulkJctSideGradingCoeff;
    double CNTfwdCapDepCoeff;
    double CNTphi; /* input - use tPhi */
    double CNTgamma;
    double CNTlambda;
    double CNTsubstrateDoping;
    int CNTgateType;
    double CNTsurfaceStateDensity;
    double CNToxideThickness;
    double CNTsurfaceMobility; /* input - use tSurfMob */
    double CNTfNcoef;
    double CNTfNexp;

    unsigned CNTtypeGiven  :1;
    unsigned CNTlatDiffGiven   :1;
    unsigned CNTjctSatCurDensityGiven  :1;
    unsigned CNTjctSatCurGiven :1;
    unsigned CNTdrainResistanceGiven   :1;
    unsigned CNTsourceResistanceGiven  :1;
    unsigned CNTsheetResistanceGiven   :1;
    unsigned CNTtransconductanceGiven  :1;
    unsigned CNTgateSourceOverlapCapFactorGiven    :1;
    unsigned CNTgateDrainOverlapCapFactorGiven :1;
    unsigned CNTgateBulkOverlapCapFactorGiven  :1;
    unsigned CNTvt0Given   :1;
    unsigned CNTcapBDGiven :1;
    unsigned CNTcapBSGiven :1;
    unsigned CNTbulkCapFactorGiven :1;
    unsigned CNTsideWallCapFactorGiven   :1;
    unsigned CNTbulkJctPotentialGiven  :1;
    unsigned CNTbulkJctBotGradingCoeffGiven    :1;
    unsigned CNTbulkJctSideGradingCoeffGiven   :1;
    unsigned CNTfwdCapDepCoeffGiven    :1;
    unsigned CNTphiGiven   :1;
    unsigned CNTgammaGiven :1;
    unsigned CNTlambdaGiven    :1;
    unsigned CNTsubstrateDopingGiven   :1;
    unsigned CNTgateTypeGiven  :1;
    unsigned CNTsurfaceStateDensityGiven   :1;
    unsigned CNToxideThicknessGiven    :1;
    unsigned CNTsurfaceMobilityGiven   :1;
    unsigned CNTtnomGiven  :1;
    unsigned CNTfNcoefGiven  :1;
    unsigned CNTfNexpGiven   :1;

} CNTmodel;

#ifndef NCNT
#define NCNT 1
#define PCNT -1
#endif /*NCNT*/

/* device parameters */
#define CNT_W 1
#define CNT_L 2
#define CNT_AS 3
#define CNT_AD 4
#define CNT_PS 5
#define CNT_PD 6
#define CNT_NRS 7
#define CNT_NRD 8
#define CNT_OFF 9
#define CNT_IC 10
#define CNT_IC_VBS 11
#define CNT_IC_VDS 12
#define CNT_IC_VGS 13
#define CNT_W_SENS 14
#define CNT_L_SENS 15
#define CNT_CB 16
#define CNT_CG 17
#define CNT_CS 18
#define CNT_POWER 19
#define CNT_TEMP 20
#define CNT_M 21
#define CNT_DTEMP 22
#define CNT_DIA                23  /*cnt diameter*/
#define CNT_EF                24  /*cnt fermi-level*/
/* model paramerers */
#define CNT_MOD_VTO 101
#define CNT_MOD_KP 102
#define CNT_MOD_GAMMA 103
#define CNT_MOD_PHI 104
#define CNT_MOD_LAMBDA 105
#define CNT_MOD_RD 106
#define CNT_MOD_RS 107
#define CNT_MOD_CBD 108
#define CNT_MOD_CBS 109
#define CNT_MOD_IS 110
#define CNT_MOD_PB 111
#define CNT_MOD_CGSO 112
#define CNT_MOD_CGDO 113
#define CNT_MOD_CGBO 114
#define CNT_MOD_CJ 115
#define CNT_MOD_MJ 116
#define CNT_MOD_CJSW 117
#define CNT_MOD_MJSW 118
#define CNT_MOD_JS 119
#define CNT_MOD_TOX 120
#define CNT_MOD_LD 121
#define CNT_MOD_RSH 122
#define CNT_MOD_U0 123
#define CNT_MOD_FC 124
#define CNT_MOD_NSUB 125
#define CNT_MOD_TPG 126
#define CNT_MOD_NSS 127
#define CNT_MOD_NCNT 128
#define CNT_MOD_PCNT 129
#define CNT_MOD_TNOM 130
#define CNT_MOD_KF 131
#define CNT_MOD_AF 132
#define CNT_MOD_TYPE 133

/* device questions */
#define CNT_CGS                201
#define CNT_CGD                202
#define CNT_DNODE              203
#define CNT_GNODE              204
#define CNT_SNODE              205
#define CNT_BNODE              206
#define CNT_DNODEPRIME         207
#define CNT_SNODEPRIME         208
#define CNT_SOURCECONDUCT      209
#define CNT_DRAINCONDUCT       210
#define CNT_VON                211
#define CNT_VDSAT              212
#define CNT_SOURCEVCRIT        213
#define CNT_DRAINVCRIT         214
#define CNT_CD                 215
#define CNT_CBS                216
#define CNT_CBD                217
#define CNT_GMBS               218
#define CNT_GM                 219
#define CNT_GDS                220
#define CNT_GBD                221
#define CNT_GBS                222
#define CNT_CAPBD              223
#define CNT_CAPBS              224
#define CNT_CAPZEROBIASBD      225
#define CNT_CAPZEROBIASBDSW    226
#define CNT_CAPZEROBIASBS      227
#define CNT_CAPZEROBIASBSSW    228
#define CNT_VBD                229
#define CNT_VBS                230
#define CNT_VGS                231
#define CNT_VDS                232
#define CNT_CAPGS              233
#define CNT_QGS                234
#define CNT_CQGS               235
#define CNT_CAPGD              236
#define CNT_QGD                237
#define CNT_CQGD               238
#define CNT_CAPGB              239
#define CNT_QGB                240
#define CNT_CQGB               241
#define CNT_QBD                242
#define CNT_CQBD               243
#define CNT_QBS                244
#define CNT_CQBS               245
#define CNT_L_SENS_REAL               246
#define CNT_L_SENS_IMAG               247
#define CNT_L_SENS_MAG                248 
#define CNT_L_SENS_PH                 249 
#define CNT_L_SENS_CPLX               250
#define CNT_W_SENS_REAL               251
#define CNT_W_SENS_IMAG               252
#define CNT_W_SENS_MAG                253 
#define CNT_W_SENS_PH                 254 
#define CNT_W_SENS_CPLX               255
#define CNT_L_SENS_DC                 256
#define CNT_W_SENS_DC                 257
#define CNT_SOURCERESIST      258
#define CNT_DRAINRESIST       259

/* model questions */

#include "cntext.h"

#endif /*CNT*/

