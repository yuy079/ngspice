/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/

extern int CNTacLoad(GENmodel*,CKTcircuit*);
extern int CNTask(CKTcircuit*,GENinstance*,int,IFvalue*,IFvalue*);
extern int CNTmAsk(CKTcircuit*,GENmodel*,int,IFvalue*);
extern int CNTconvTest(GENmodel*,CKTcircuit*);
extern int CNTdelete(GENmodel*,IFuid,GENinstance**);
extern void CNTdestroy(GENmodel**);
extern int CNTgetic(GENmodel*,CKTcircuit*);
extern int CNTload(GENmodel*,CKTcircuit*);
extern int CNTmDelete(GENmodel**,IFuid,GENmodel*);
extern int CNTmParam(int,IFvalue*,GENmodel*);
extern int CNTparam(int,IFvalue*,GENinstance*,IFvalue*);
extern int CNTpzLoad(GENmodel*,CKTcircuit*,SPcomplex*);
extern int CNTsAcLoad(GENmodel*,CKTcircuit*);
extern int CNTsLoad(GENmodel*,CKTcircuit*);
extern void CNTsPrint(GENmodel*,CKTcircuit*);
extern int CNTsSetup(SENstruct*,GENmodel*);
extern int CNTsUpdate(GENmodel*,CKTcircuit*);
extern int CNTsetup(SMPmatrix*,GENmodel*,CKTcircuit*,int*);
extern int CNTunsetup(GENmodel*,CKTcircuit*);
extern int CNTtemp(GENmodel*,CKTcircuit*);
extern int CNTtrunc(GENmodel*,CKTcircuit*,double*);
extern int CNTdisto(int,GENmodel*,CKTcircuit*);
extern int CNTnoise(int,int,GENmodel*,CKTcircuit*,Ndata*,double*);
extern int CNTdSetup(GENmodel*,CKTcircuit*);

extern double Coeff();
extern double Root();
extern int Range1();
extern int Range2();
extern double Ncnt();
