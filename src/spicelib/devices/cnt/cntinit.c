#include "config.h"

#include "devdefs.h"

#include "cntitf.h"
#include "cntext.h"
#include "cntinit.h"


SPICEdev CNTinfo = {
    {
        "Cnt",
        "CNTfet model with Meyer capacitance model",

        &CNTnSize,
        &CNTnSize,
        CNTnames,

        &CNTpTSize,
        CNTpTable,

        &CNTmPTSize,
        CNTmPTable,

#ifdef XSPICE
/*----  Fixed by SDB 5.2.2003 to enable XSPICE/tclspice integration  -----*/
        NULL,  /* This is a SPICE device, it has no MIF info data */

        0,     /* This is a SPICE device, it has no MIF info data */
        NULL,  /* This is a SPICE device, it has no MIF info data */

        0,     /* This is a SPICE device, it has no MIF info data */
        NULL,  /* This is a SPICE device, it has no MIF info data */

        0,     /* This is a SPICE device, it has no MIF info data */
        NULL,  /* This is a SPICE device, it has no MIF info data */
/*---------------------------  End of SDB fix   -------------------------*/
#endif

	DEV_DEFAULT
    },

 /* DEVparam      */ CNTparam,
 /* DEVmodParam   */ CNTmParam,
 /* DEVload       */ CNTload,
 /* DEVsetup      */ CNTsetup,
 /* DEVunsetup    */ CNTunsetup,
 /* DEVpzSetup    */ CNTsetup,
 /* DEVtemperature*/ CNTtemp,
 /* DEVtrunc      */ CNTtrunc,
 /* DEVfindBranch */ NULL,
 /* DEVacLoad     */ CNTacLoad,
 /* DEVaccept     */ NULL,
 /* DEVdestroy    */ CNTdestroy,
 /* DEVmodDelete  */ CNTmDelete,
 /* DEVdelete     */ CNTdelete,
 /* DEVsetic      */ CNTgetic,
 /* DEVask        */ CNTask,
 /* DEVmodAsk     */ CNTmAsk,
 /* DEVpzLoad     */ CNTpzLoad,
 /* DEVconvTest   */ CNTconvTest,
 /* DEVsenSetup   */ CNTsSetup,
 /* DEVsenLoad    */ CNTsLoad,
 /* DEVsenUpdate  */ CNTsUpdate,
 /* DEVsenAcLoad  */ CNTsAcLoad,
 /* DEVsenPrint   */ CNTsPrint,
 /* DEVsenTrunc   */ NULL,
 /* DEVdisto      */ CNTdisto,
 /* DEVnoise      */ CNTnoise,
#ifdef CIDER
 /* DEVdump       */ NULL,
 /* DEVacct       */ NULL,
#endif    
 /* DEVinstSize   */ &CNTiSize,
 /* DEVmodSize    */ &CNTmSize
};


SPICEdev *
get_cnt_info(void)
{
    return &CNTinfo;
}
