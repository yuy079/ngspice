/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1987 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/

#include "ngspice/ngspice.h"
#include "ngspice/devdefs.h"
#include "ngspice/ifsim.h"
#include "cntdefs.h"
#include "ngspice/suffix.h"

IFparm CNTpTable[] = { /* parameters */ 
 IOPU("m",            CNT_M,          IF_REAL   , "Multiplier"),
 IOPU("l",            CNT_L,          IF_REAL   , "Length"),
 IOPU("w",            CNT_W,          IF_REAL   , "Width"),
 IOPU("dia",       CNT_DIA,     IF_REAL   , "Diameter"), /*cnt diameter*/
 IOPU("ef",       CNT_EF,     IF_REAL   , "Fermi-level"), /*cnt fermi-level*/
 IOPU("ad",           CNT_AD,         IF_REAL   , "Drain area"),
 IOPU("as",           CNT_AS,         IF_REAL   , "Source area"),
 IOPU("pd",           CNT_PD,         IF_REAL   , "Drain perimeter"),
 IOPU("ps",           CNT_PS,         IF_REAL   , "Source perimeter"),
 IOPU("nrd",          CNT_NRD,        IF_REAL   , "Drain squares"),
 IOPU("nrs",          CNT_NRS,        IF_REAL   , "Source squares"),
 IP("off",           CNT_OFF,        IF_FLAG   , "Device initially off"),
 IOPU("icvds",        CNT_IC_VDS,     IF_REAL   , "Initial D-S voltage"),
 IOPU("icvgs",        CNT_IC_VGS,     IF_REAL   , "Initial G-S voltage"),
 IOPU("icvbs",        CNT_IC_VBS,     IF_REAL   , "Initial B-S voltage"),
 IOPU("temp",         CNT_TEMP,       IF_REAL,    "Instance temperature"),
 IOPU("dtemp",         CNT_DTEMP,       IF_REAL,    "Instance temperature difference"),
 IP( "ic",           CNT_IC,  IF_REALVEC, "Vector of D-S, G-S, B-S voltages"),
 IP( "sens_l", CNT_L_SENS, IF_FLAG, "flag to request sensitivity WRT length"),
 IP( "sens_w", CNT_W_SENS, IF_FLAG, "flag to request sensitivity WRT width"),

 OP( "id",           CNT_CD,         IF_REAL,    "Drain current"),
 OP( "is",           CNT_CS,         IF_REAL,    "Source current"),
 OP( "ig",           CNT_CG,         IF_REAL,    "Gate current "),
 OP( "ib",           CNT_CB,         IF_REAL,    "Bulk current "),
 OPU( "ibd",      CNT_CBD,    IF_REAL,    "B-D junction current"),
 OPU( "ibs",      CNT_CBS,    IF_REAL,    "B-S junction current"),
 OP( "vgs",          CNT_VGS,        IF_REAL,    "Gate-Source voltage"),
 OP( "vds",          CNT_VDS,        IF_REAL,    "Drain-Source voltage"),
 OP( "vbs",          CNT_VBS,        IF_REAL,    "Bulk-Source voltage"),
 OPU( "vbd",          CNT_VBD,        IF_REAL,    "Bulk-Drain voltage"),
 /*
 OP( "cgs",          CNT_CGS,        IF_REAL   , "Gate-Source capacitance"),
 OP( "cgd",          CNT_CGD,        IF_REAL   , "Gate-Drain capacitance"),
 */

 OPU( "dnode",      CNT_DNODE,      IF_INTEGER, "Number of the drain node "),
 OPU( "gnode",      CNT_GNODE,      IF_INTEGER, "Number of the gate node "),
 OPU( "snode",      CNT_SNODE,      IF_INTEGER, "Number of the source node "),
 OPU( "bnode",      CNT_BNODE,      IF_INTEGER, "Number of the node "),
 OPU( "dnodeprime", CNT_DNODEPRIME, IF_INTEGER, "Number of int. drain node"),
 OPU( "snodeprime", CNT_SNODEPRIME, IF_INTEGER, "Number of int. source node "),

 OP( "von",          CNT_VON,        IF_REAL,    " "),
 OP( "vdsat",        CNT_VDSAT,      IF_REAL,    "Saturation drain voltage"),
 OPU( "sourcevcrit",  CNT_SOURCEVCRIT,IF_REAL,    "Critical source voltage"),
 OPU( "drainvcrit",   CNT_DRAINVCRIT, IF_REAL,    "Critical drain voltage"),
 OP( "rs", CNT_SOURCERESIST, IF_REAL, "Source resistance"),
 OPU("sourceconductance", CNT_SOURCECONDUCT, IF_REAL, "Conductance of source"),
 OP( "rd",  CNT_DRAINRESIST,  IF_REAL, "Drain conductance"),
 OPU("drainconductance",  CNT_DRAINCONDUCT,  IF_REAL, "Conductance of drain"),

 OP( "gm",           CNT_GM,         IF_REAL,    "Transconductance"),
 OP( "gds",          CNT_GDS,        IF_REAL,    "Drain-Source conductance"),
 OP( "gmb",     CNT_GMBS,   IF_REAL,    "Bulk-Source transconductance"),
 OPR( "gmbs",     CNT_GMBS,   IF_REAL,    ""),
 OPU( "gbd",          CNT_GBD,        IF_REAL,    "Bulk-Drain conductance"),
 OPU( "gbs",          CNT_GBS,        IF_REAL,    "Bulk-Source conductance"),

 OP( "cbd",        CNT_CAPBD,      IF_REAL,    "Bulk-Drain capacitance"),
 OP( "cbs",        CNT_CAPBS,      IF_REAL,    "Bulk-Source capacitance"),
 OP( "cgs",        CNT_CAPGS,      IF_REAL,    "Gate-Source capacitance"),
 OP( "cgd",        CNT_CAPGD,      IF_REAL,    "Gate-Drain capacitance"),
 OP( "cgb",        CNT_CAPGB,      IF_REAL,    "Gate-Bulk capacitance"),

 OPU( "cqgs",CNT_CQGS,IF_REAL,"Capacitance due to gate-source charge storage"),
 OPU( "cqgd",CNT_CQGD,IF_REAL,"Capacitance due to gate-drain charge storage"),
 OPU( "cqgb",CNT_CQGB,IF_REAL,"Capacitance due to gate-bulk charge storage"),
 OPU( "cqbd",CNT_CQBD,IF_REAL,"Capacitance due to bulk-drain charge storage"),
 OPU( "cqbs",CNT_CQBS,IF_REAL,"Capacitance due to bulk-source charge storage"),

 OP( "cbd0", CNT_CAPZEROBIASBD, IF_REAL, "Zero-Bias B-D junction capacitance"),
 OP( "cbdsw0",        CNT_CAPZEROBIASBDSW, IF_REAL,    " "),
 OP( "cbs0", CNT_CAPZEROBIASBS, IF_REAL, "Zero-Bias B-S junction capacitance"),
 OP( "cbssw0",        CNT_CAPZEROBIASBSSW, IF_REAL,    " "),

 OPU( "qgs",         CNT_QGS,        IF_REAL,    "Gate-Source charge storage"),
 OPU( "qgd",         CNT_QGD,        IF_REAL,    "Gate-Drain charge storage"),
 OPU( "qgb",         CNT_QGB,        IF_REAL,    "Gate-Bulk charge storage"),
 OPU( "qbd",         CNT_QBD,        IF_REAL,    "Bulk-Drain charge storage"),
 OPU( "qbs",         CNT_QBS,        IF_REAL,    "Bulk-Source charge storage"),
 OPU( "p",            CNT_POWER,      IF_REAL,    "Instaneous power"),
 OPU( "sens_l_dc",    CNT_L_SENS_DC,  IF_REAL,    "dc sensitivity wrt length"),
 OPU( "sens_l_real", CNT_L_SENS_REAL,IF_REAL,
        "real part of ac sensitivity wrt length"),
 OPU( "sens_l_imag",  CNT_L_SENS_IMAG,IF_REAL,    
        "imag part of ac sensitivity wrt length"),
 OPU( "sens_l_mag",   CNT_L_SENS_MAG, IF_REAL,    
        "sensitivity wrt l of ac magnitude"),
 OPU( "sens_l_ph",    CNT_L_SENS_PH,  IF_REAL,    
        "sensitivity wrt l of ac phase"),
 OPU( "sens_l_cplx",  CNT_L_SENS_CPLX,IF_COMPLEX, "ac sensitivity wrt length"),
 OPU( "sens_w_dc",    CNT_W_SENS_DC,  IF_REAL,    "dc sensitivity wrt width"),
 OPU( "sens_w_real",  CNT_W_SENS_REAL,IF_REAL,    
        "real part of ac sensitivity wrt width"),
 OPU( "sens_w_imag",  CNT_W_SENS_IMAG,IF_REAL,    
        "imag part of ac sensitivity wrt width"),
 OPU( "sens_w_mag",   CNT_W_SENS_MAG, IF_REAL,    
        "sensitivity wrt w of ac magnitude"),
 OPU( "sens_w_ph",    CNT_W_SENS_PH,  IF_REAL,    
        "sensitivity wrt w of ac phase"),
 OPU( "sens_w_cplx",  CNT_W_SENS_CPLX,IF_COMPLEX, "ac sensitivity wrt width")
};

IFparm CNTmPTable[] = { /* model parameters */
 OP("type",   CNT_MOD_TYPE,  IF_STRING, "N-channel or P-channel CNT"),
 IOP("vto",   CNT_MOD_VTO,   IF_REAL   ,"Threshold voltage"),
 IOPR("vt0",  CNT_MOD_VTO,   IF_REAL   ,"Threshold voltage"),
 IOP("kp",    CNT_MOD_KP,    IF_REAL   ,"Transconductance parameter"),
 IOP("gamma", CNT_MOD_GAMMA, IF_REAL   ,"Bulk threshold parameter"),
 IOP("phi",   CNT_MOD_PHI,   IF_REAL   ,"Surface potential"),
 IOP("lambda",CNT_MOD_LAMBDA,IF_REAL   ,"Channel length modulation"),
 IOP("rd",    CNT_MOD_RD,    IF_REAL   ,"Drain ohmic resistance"),
 IOP("rs",    CNT_MOD_RS,    IF_REAL   ,"Source ohmic resistance"),
 IOPA("cbd",  CNT_MOD_CBD,   IF_REAL   ,"B-D junction capacitance"),
 IOPA("cbs",  CNT_MOD_CBS,   IF_REAL   ,"B-S junction capacitance"),
 IOP("is",    CNT_MOD_IS,    IF_REAL   ,"Bulk junction sat. current"),
 IOP("pb",    CNT_MOD_PB,    IF_REAL   ,"Bulk junction potential"),
 IOPA("cgso", CNT_MOD_CGSO,  IF_REAL   ,"Gate-source overlap cap."),
 IOPA("cgdo", CNT_MOD_CGDO,  IF_REAL   ,"Gate-drain overlap cap."),
 IOPA("cgbo", CNT_MOD_CGBO,  IF_REAL   ,"Gate-bulk overlap cap."),
 IOP("rsh",   CNT_MOD_RSH,   IF_REAL   ,"Sheet resistance"),
 IOPA("cj",   CNT_MOD_CJ,    IF_REAL   ,"Bottom junction cap per area"),
 IOP("mj",    CNT_MOD_MJ,    IF_REAL   ,"Bottom grading coefficient"),
 IOPA("cjsw", CNT_MOD_CJSW,  IF_REAL   ,"Side junction cap per area"),
 IOP("mjsw",  CNT_MOD_MJSW,  IF_REAL   ,"Side grading coefficient"),
 IOP("js",    CNT_MOD_JS,    IF_REAL   ,"Bulk jct. sat. current density"),
 IOP("tox",   CNT_MOD_TOX,   IF_REAL   ,"Oxide thickness"),
 IOP("ld",    CNT_MOD_LD,    IF_REAL   ,"Lateral diffusion"),
 IOP("u0",    CNT_MOD_U0,    IF_REAL   ,"Surface mobility"),
 IOPR("uo",   CNT_MOD_U0,    IF_REAL   ,"Surface mobility"),
 IOP("fc",    CNT_MOD_FC,    IF_REAL   ,"Forward bias jct. fit parm."),
 IP("ncnt",   CNT_MOD_NCNT,  IF_FLAG   ,"N type CNT Transistor model"),
 IP("pcnt",   CNT_MOD_PCNT,  IF_FLAG   ,"P type CNT Transistor model"),
 IOP("nsub",  CNT_MOD_NSUB,  IF_REAL   ,"Substrate doping"),
 IOP("tpg",   CNT_MOD_TPG,   IF_INTEGER,"Gate type"),
 IOP("nss",   CNT_MOD_NSS,   IF_REAL   ,"Surface state density"),
 IOP("tnom",  CNT_MOD_TNOM,  IF_REAL   ,"Parameter measurement temperature"),
 IOP("kf",     CNT_MOD_KF,    IF_REAL   ,"Flicker noise coefficient"),
 IOP("af",     CNT_MOD_AF,    IF_REAL   ,"Flicker noise exponent")
};

char *CNTnames[] = {
    "Drain",
    "Gate",
    "Source",
    "Bulk"
};

int	CNTnSize = NUMELEMS(CNTnames);
int	CNTpTSize = NUMELEMS(CNTpTable);
int	CNTmPTSize = NUMELEMS(CNTmPTable);
int	CNTiSize = sizeof(CNTinstance);
int	CNTmSize = sizeof(CNTmodel);
