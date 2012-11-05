#ifndef _RHS_H_
#define _RHS_H_

/* constants for the single cell detailed circadian model in Kim 2012 */ 
#define	trPo 25.9201
#define	trPt 44.854
#define	trRo 23.0747
#define	trRt 39.9409
#define	trB 46.1038
#define	trRev 102.923
#define	trNp 0.329749
#define	tlp 1.81031
#define	tlr 5.03882
#define	tlb 0.530436
#define	tlrev 8.90744
#define	tlc 4.64589
#define	tlnp 1.25099
#define	agp 1.3962
#define	dg 2.93521
#define	ac 0.0456572
#define	dc 0.108072
#define	ar 0.0235285
#define	dr 0.605268
#define	cbin 0.0454894
#define	uncbin 7.27215
#define	bbin 6.92686
#define	unbbin 0.130196
#define	cbbin 6.59924
#define	uncbbin 0.304176
#define	ag 0.162392
#define	bin 6.97166
#define	unbin 0.255032
#define	binrev 0.0120525
#define	unbinrev 10.9741
#define	binr 6.15445
#define	unbinr 2.91009
#define	binc 0.280863
#define	unbinc 0.00886752
#define	binrevb 0.00626588
#define	unbinrevb 5.30559
#define	tmc 0.16426
#define	tmcrev 9.2631
#define	nl 0.643086
#define	ne 0.0269078
#define	nlrev 9.63702
#define	nerev 0.0152514
#define	lne 0.594609
#define	nlbc 5.26501
#define	hoo 0.527453
#define	hto 2.45584
#define	phos 0.291429
#define	lono 0.205813
#define	lont 0.396392
#define	lta 0.607387
#define	ltb 0.013
#define	trgto 0.644602
#define	ugto 0.0625777
#define	Nf 3.35063
#define	up 3.537
#define	uro 0.17491
#define	urt 0.481895
#define	umNp 0.369493
#define	umPo 0.766962
#define	umPt 0.58892
#define	umRo 0.403425
#define	umRt 0.455544
#define	ub 0.0188002
#define	uc 0.0251651
#define	ubc 0.348829
#define	upu 0.0700322
#define	urev 1.64876
#define	uprev 0.517303
#define	umB 0.795402
#define	umRev 1.51019

/* define the functions used in the differential equations*/ 
float rhsGR(float GR, float G, float x01011, float x02011);
float rhsG(float GR, float G, float x00011);
float rhsGrR(float G, float GrR, float Gr, float x01011, float x02011);
float rhsGr(float G, float GrR, float Gr, float x00011);
float rhsGcR(float G, float GcR, float Gc, float x01011, float x02011);
float rhsGc(float G, float GcR, float Gc, float x00011);
float rhsGBR(float G, float GBR, float GB, float B, float revn, float revng, float revngp, float revnp);
float rhsGB(float G, float GBR, float GB, float B, float revn, float revng, float revngp, float revnp);
float rhsGBRb(float G, float GBR, float GB, float GBRb, float GBb, float B, float revn, float revng, float revngp, float revnp);
float rhsGBb(float G, float GBR, float GB, float GBRb, float GBb, float B, float revn, float revng, float revngp, float revnp);
float rhsMnPo(float G, float MnPo);
float rhsMcPo(float MnPo, float McPo);
float rhsMnPt(float G, float MnPt);
float rhsMcPt(float MnPt, float McPt);
float rhsMnRt(float G, float Gc, float MnRt);
float rhsMcRt(float MnRt, float McRt);
float rhsMnRev(float G, float Gr, float MnRev, float x00011);
float rhsMcRev(float MnRev, float McRev);
float rhsMnRo(float G, float GB, float MnRo, float B);
float rhsMcRo(float MnRo, float McRo);
float rhsMnB(float G, float GB, float GBb, float MnB, float B);
float rhsMcB(float MnB, float McB, float B);
float rhsMnNp(float G, float GB, float MnNp, float B);
float rhsMcNp(float MnNp, float McNp);
float rhsB(float McB, float B, float Cl, float BC);
float rhsCl(float McNp, float B, float Cl, float BC);
float rhsBC(float B, float Cl, float BC);
float rhscyrev(float McRev, float cyrev, float revn, float cyrevg, float x00200);
float rhsrevn(float cyrev, float revn, float revng, float x00210);
float rhscyrevg(float cyrev, float revn, float cyrevg, float revng, float gto, float x00200);
float rhsrevng(float cyrev, float revn, float cyrevg, float revng, float gto, float x00210);
float rhscyrevgp(float cyrev, float revn, float cyrevg, float revng, float cyrevgp, float revngp, float gto);
float rhsrevngp(float cyrev, float revn, float cyrevg, float revng, float cyrevgp, float revngp, float gto);
float rhscyrevp(float cyrev, float revn, float cyrevg, float cyrevgp, float cyrevp, float revnp);
float rhsrevnp(float cyrev, float revn, float revng, float revngp, float cyrevp, float revnp);
float rhsgto(float G, float GB, float B, float gto);
float rhsx00001(float B, float BC, float x00001);
float rhsx00011(float x00001, float x00011, float x01010, float x01011, float x02010, float x02011, float x20010, float x20011, float x20110, float x20111, float x21010, float x21011, float x21110, float x21111, float x22010, float x22011, float x22110, float x22111, float x40010, float x40011, float x40110, float x40111, float x40210, float x40211, float x40310, float x40311, float x41010, float x41011, float x41110, float x41111, float x41210, float x41211, float x41310, float x41311, float x42010, float x42011, float x42110, float x42111, float x42210, float x42211, float x42310, float x42311, float x50010, float x50011, float x50110, float x50111, float x50210, float x50211, float x50310, float x50311, float x51010, float x51011, float x51110, float x51111, float x51210, float x51211, float x51310, float x51311, float x52010, float x52011, float x52110, float x52111, float x52210, float x52211, float x52310, float x52311, float x60010, float x60011, float x60110, float x60111, float x60210, float x60211, float x60310, float x60311, float x61010, float x61011, float x61110, float x61111, float x61210, float x61211, float x61310, float x61311, float x62010, float x62011, float x62110, float x62111, float x62210, float x62211, float x62310, float x62311);
float rhsx00100(float x00100, float x00110, float x10000, float x10100, float x20000, float x20100, float x21000, float x21100, float x22000, float x22100, float x30000, float x30100, float x30200, float x30300, float x40000, float x40100, float x40200, float x40300, float x41000, float x41100, float x41200, float x41300, float x42000, float x42100, float x42200, float x42300, float x50000, float x50100, float x50200, float x50300, float x51000, float x51100, float x51200, float x51300, float x52000, float x52100, float x52200, float x52300, float x60000, float x60100, float x60200, float x60300, float x61000, float x61100, float x61200, float x61300, float x62000, float x62100, float x62200, float x62300);
float rhsx00110(float x00110, float x20010, float x20011, float x20110, float x20111, float x21010, float x21011, float x21110, float x21111, float x22010, float x22011, float x22110, float x22111, float x40010, float x40011, float x40110, float x40111, float x40210, float x40211, float x40310, float x40311, float x41010, float x41011, float x41110, float x41111, float x41210, float x41211, float x41310, float x41311, float x42010, float x42011, float x42110, float x42111, float x42210, float x42211, float x42310, float x42311, float x50010, float x50011, float x50110, float x50111, float x50210, float x50211, float x50310, float x50311, float x51010, float x51011, float x51110, float x51111, float x51210, float x51211, float x51310, float x51311, float x52010, float x52011, float x52110, float x52111, float x52210, float x52211, float x52310, float x52311, float x60010, float x60011, float x60110, float x60111, float x60210, float x60211, float x60310, float x60311, float x61010, float x61011, float x61110, float x61111, float x61210, float x61211, float x61310, float x61311, float x62010, float x62011, float x62110, float x62111, float x62210, float x62211, float x62310, float x62311);
float rhsx00200(float cyrev, float cyrevg, float cyrevgp, float x00200, float x00210, float x30000, float x30100, float x30200, float x30300, float x40000, float x40100, float x40200, float x40300, float x41000, float x41100, float x41200, float x41300, float x42000, float x42100, float x42200, float x42300, float x50000, float x50100, float x50200, float x50300, float x51000, float x51100, float x51200, float x51300, float x52000, float x52100, float x52200, float x52300, float x60000, float x60100, float x60200, float x60300, float x61000, float x61100, float x61200, float x61300, float x62000, float x62100, float x62200, float x62300);
float rhsx00210(float revn, float revng, float revngp, float x00210, float x40010, float x40011, float x40110, float x40111, float x40210, float x40211, float x40310, float x40311, float x41010, float x41011, float x41110, float x41111, float x41210, float x41211, float x41310, float x41311, float x42010, float x42011, float x42110, float x42111, float x42210, float x42211, float x42310, float x42311, float x50010, float x50011, float x50110, float x50111, float x50210, float x50211, float x50310, float x50311, float x51010, float x51011, float x51110, float x51111, float x51210, float x51211, float x51310, float x51311, float x52010, float x52011, float x52110, float x52111, float x52210, float x52211, float x52310, float x52311, float x60010, float x60011, float x60110, float x60111, float x60210, float x60211, float x60310, float x60311, float x61010, float x61011, float x61110, float x61111, float x61210, float x61211, float x61310, float x61311, float x62010, float x62011, float x62110, float x62111, float x62210, float x62211, float x62310, float x62311);
float rhsx01000(float McRo, float x01000, float x20000, float x20100, float x21000, float x21100, float x40000, float x40100, float x40200, float x40300, float x41000, float x41100, float x41200, float x41300, float x50000, float x50100, float x50200, float x50300, float x51000, float x51100, float x51200, float x51300, float x60000, float x60100, float x60200, float x60300, float x61000, float x61100, float x61200, float x61300);
float rhsx01010(float x00011, float x01010, float x01011, float x20010, float x20011, float x20110, float x20111, float x21010, float x21011, float x21110, float x21111, float x40010, float x40011, float x40110, float x40111, float x40210, float x40211, float x40310, float x40311, float x41010, float x41011, float x41110, float x41111, float x41210, float x41211, float x41310, float x41311, float x50010, float x50011, float x50110, float x50111, float x50210, float x50211, float x50310, float x50311, float x51010, float x51011, float x51110, float x51111, float x51210, float x51211, float x51310, float x51311, float x60010, float x60011, float x60110, float x60111, float x60210, float x60211, float x60310, float x60311, float x61010, float x61011, float x61110, float x61111, float x61210, float x61211, float x61310, float x61311);
float rhsx01011(float x00011, float x01010, float x01011, float x20010, float x20110, float x21011, float x21111, float x40010, float x40110, float x40210, float x40310, float x41011, float x41111, float x41211, float x41311, float x50010, float x50110, float x50210, float x50310, float x51011, float x51111, float x51211, float x51311, float x60010, float x60110, float x60210, float x60310, float x61011, float x61111, float x61211, float x61311);
float rhsx02000(float McRt, float x02000, float x20000, float x20100, float x22000, float x22100, float x40000, float x40100, float x40200, float x40300, float x42000, float x42100, float x42200, float x42300, float x50000, float x50100, float x50200, float x50300, float x52000, float x52100, float x52200, float x52300, float x60000, float x60100, float x60200, float x60300, float x62000, float x62100, float x62200, float x62300);
float rhsx02010(float x00011, float x02010, float x02011, float x20010, float x20011, float x20110, float x20111, float x22010, float x22011, float x22110, float x22111, float x40010, float x40011, float x40110, float x40111, float x40210, float x40211, float x40310, float x40311, float x42010, float x42011, float x42110, float x42111, float x42210, float x42211, float x42310, float x42311, float x50010, float x50011, float x50110, float x50111, float x50210, float x50211, float x50310, float x50311, float x52010, float x52011, float x52110, float x52111, float x52210, float x52211, float x52310, float x52311, float x60010, float x60011, float x60110, float x60111, float x60210, float x60211, float x60310, float x60311, float x62010, float x62011, float x62110, float x62111, float x62210, float x62211, float x62310, float x62311);
float rhsx02011(float x00011, float x02010, float x02011, float x20010, float x20110, float x22011, float x22111, float x40010, float x40110, float x40210, float x40310, float x42011, float x42111, float x42211, float x42311, float x50010, float x50110, float x50210, float x50310, float x52011, float x52111, float x52211, float x52311, float x60010, float x60110, float x60210, float x60310, float x62011, float x62111, float x62211, float x62311);
float rhsx10000(float McPo, float x00100, float x10000, float x10100);
float rhsx10100(float x00100, float x10000, float x10100);
float rhsx20000(float x00100, float x01000, float x02000, float x20000, float x20010, float x20100, float x21000, float x22000);
float rhsx20010(float x00011, float x00110, float x01010, float x01011, float x02010, float x02011, float x20000, float x20010, float x20011, float x20110, float x21010, float x21011, float x22010, float x22011);
float rhsx20011(float x00011, float x00110, float x01010, float x02010, float x20010, float x20011, float x20111, float x21011, float x22011);
float rhsx20100(float x00100, float x01000, float x02000, float x10100, float x20000, float x20100, float x20110, float x21100, float x22100);
float rhsx20110(float x00011, float x00110, float x01010, float x01011, float x02010, float x02011, float x20010, float x20100, float x20110, float x20111, float x21110, float x21111, float x22110, float x22111);
float rhsx20111(float x00011, float x00110, float x01010, float x02010, float x20011, float x20110, float x20111, float x21111, float x22111);
float rhsx21000(float x00100, float x01000, float x20000, float x21000, float x21010, float x21100);
float rhsx21010(float x00011, float x00110, float x01010, float x20010, float x21000, float x21010, float x21011, float x21110);
float rhsx21011(float x00011, float x00110, float x01010, float x01011, float x20010, float x20011, float x21010, float x21011, float x21111);
float rhsx21100(float x00100, float x01000, float x20100, float x21000, float x21100, float x21110);
float rhsx21110(float x00011, float x00110, float x01010, float x20110, float x21010, float x21100, float x21110, float x21111);
float rhsx21111(float x00011, float x00110, float x01010, float x01011, float x20110, float x20111, float x21011, float x21110, float x21111);
float rhsx22000(float x00100, float x02000, float x20000, float x22000, float x22010, float x22100);
float rhsx22010(float x00011, float x00110, float x02010, float x20010, float x22000, float x22010, float x22011, float x22110);
float rhsx22011(float x00011, float x00110, float x02010, float x02011, float x20010, float x20011, float x22010, float x22011, float x22111);
float rhsx22100(float x00100, float x02000, float x20100, float x22000, float x22100, float x22110);
float rhsx22110(float x00011, float x00110, float x02010, float x20110, float x22010, float x22100, float x22110, float x22111);
float rhsx22111(float x00011, float x00110, float x02010, float x02011, float x20110, float x20111, float x22011, float x22110, float x22111);
float rhsx30000(float McPt, float x00100, float x00200, float x30000, float x30100, float x30200);
float rhsx30100(float x00100, float x00200, float x30000, float x30100, float x30300);
float rhsx30200(float gto, float x00100, float x00200, float x30000, float x30200, float x30300);
float rhsx30300(float gto, float x00100, float x00200, float x30100, float x30200, float x30300);
float rhsx40000(float x00100, float x00200, float x01000, float x02000, float x40000, float x40010, float x40100, float x40200, float x41000, float x42000);
float rhsx40010(float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x40000, float x40010, float x40011, float x40110, float x40210, float x41010, float x41011, float x42010, float x42011);
float rhsx40011(float x00011, float x00110, float x00210, float x01010, float x02010, float x40010, float x40011, float x40111, float x40211, float x41011, float x42011);
float rhsx40100(float x00100, float x00200, float x01000, float x02000, float x30100, float x40000, float x40100, float x40110, float x40300, float x41100, float x42100);
float rhsx40110(float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x40010, float x40100, float x40110, float x40111, float x40310, float x41110, float x41111, float x42110, float x42111);
float rhsx40111(float x00011, float x00110, float x00210, float x01010, float x02010, float x40011, float x40110, float x40111, float x40311, float x41111, float x42111);
float rhsx40200(float gto, float x00100, float x00200, float x01000, float x02000, float x40000, float x40200, float x40210, float x40300, float x41200, float x42200);
float rhsx40210(float gto, float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x40010, float x40200, float x40210, float x40211, float x40310, float x41210, float x41211, float x42210, float x42211);
float rhsx40211(float gto, float x00011, float x00110, float x00210, float x01010, float x02010, float x40011, float x40210, float x40211, float x40311, float x41211, float x42211);
float rhsx40300(float gto, float x00100, float x00200, float x01000, float x02000, float x30300, float x40100, float x40200, float x40300, float x40310, float x41300, float x42300);
float rhsx40310(float gto, float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x40110, float x40210, float x40300, float x40310, float x40311, float x41310, float x41311, float x42310, float x42311);
float rhsx40311(float gto, float x00011, float x00110, float x00210, float x01010, float x02010, float x40111, float x40211, float x40310, float x40311, float x41311, float x42311);
float rhsx41000(float x00100, float x00200, float x01000, float x40000, float x41000, float x41010, float x41100, float x41200);
float rhsx41010(float x00011, float x00110, float x00210, float x01010, float x40010, float x41000, float x41010, float x41011, float x41110, float x41210);
float rhsx41011(float x00011, float x00110, float x00210, float x01010, float x01011, float x40010, float x40011, float x41010, float x41011, float x41111, float x41211);
float rhsx41100(float x00100, float x00200, float x01000, float x40100, float x41000, float x41100, float x41110, float x41300);
float rhsx41110(float x00011, float x00110, float x00210, float x01010, float x40110, float x41010, float x41100, float x41110, float x41111, float x41310);
float rhsx41111(float x00011, float x00110, float x00210, float x01010, float x01011, float x40110, float x40111, float x41011, float x41110, float x41111, float x41311);
float rhsx41200(float gto, float x00100, float x00200, float x01000, float x40200, float x41000, float x41200, float x41210, float x41300);
float rhsx41210(float gto, float x00011, float x00110, float x00210, float x01010, float x40210, float x41010, float x41200, float x41210, float x41211, float x41310);
float rhsx41211(float gto, float x00011, float x00110, float x00210, float x01010, float x01011, float x40210, float x40211, float x41011, float x41210, float x41211, float x41311);
float rhsx41300(float gto, float x00100, float x00200, float x01000, float x40300, float x41100, float x41200, float x41300, float x41310);
float rhsx41310(float gto, float x00011, float x00110, float x00210, float x01010, float x40310, float x41110, float x41210, float x41300, float x41310, float x41311);
float rhsx41311(float gto, float x00011, float x00110, float x00210, float x01010, float x01011, float x40310, float x40311, float x41111, float x41211, float x41310, float x41311);
float rhsx42000(float x00100, float x00200, float x02000, float x40000, float x42000, float x42010, float x42100, float x42200);
float rhsx42010(float x00011, float x00110, float x00210, float x02010, float x40010, float x42000, float x42010, float x42011, float x42110, float x42210);
float rhsx42011(float x00011, float x00110, float x00210, float x02010, float x02011, float x40010, float x40011, float x42010, float x42011, float x42111, float x42211);
float rhsx42100(float x00100, float x00200, float x02000, float x40100, float x42000, float x42100, float x42110, float x42300);
float rhsx42110(float x00011, float x00110, float x00210, float x02010, float x40110, float x42010, float x42100, float x42110, float x42111, float x42310);
float rhsx42111(float x00011, float x00110, float x00210, float x02010, float x02011, float x40110, float x40111, float x42011, float x42110, float x42111, float x42311);
float rhsx42200(float gto, float x00100, float x00200, float x02000, float x40200, float x42000, float x42200, float x42210, float x42300);
float rhsx42210(float gto, float x00011, float x00110, float x00210, float x02010, float x40210, float x42010, float x42200, float x42210, float x42211, float x42310);
float rhsx42211(float gto, float x00011, float x00110, float x00210, float x02010, float x02011, float x40210, float x40211, float x42011, float x42210, float x42211, float x42311);
float rhsx42300(float gto, float x00100, float x00200, float x02000, float x40300, float x42100, float x42200, float x42300, float x42310);
float rhsx42310(float gto, float x00011, float x00110, float x00210, float x02010, float x40310, float x42110, float x42210, float x42300, float x42310, float x42311);
float rhsx42311(float gto, float x00011, float x00110, float x00210, float x02010, float x02011, float x40310, float x40311, float x42111, float x42211, float x42310, float x42311);
float rhsx50000(float x00100, float x00200, float x01000, float x02000, float x50000, float x50010, float x50100, float x50200, float x51000, float x52000);
float rhsx50010(float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x50000, float x50010, float x50011, float x50110, float x50210, float x51010, float x51011, float x52010, float x52011);
float rhsx50011(float x00011, float x00110, float x00210, float x01010, float x02010, float x50010, float x50011, float x50111, float x50211, float x51011, float x52011);
float rhsx50100(float x00100, float x00200, float x01000, float x02000, float x50000, float x50100, float x50110, float x50300, float x51100, float x52100);
float rhsx50110(float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x50010, float x50100, float x50110, float x50111, float x50310, float x51110, float x51111, float x52110, float x52111);
float rhsx50111(float x00011, float x00110, float x00210, float x01010, float x02010, float x50011, float x50110, float x50111, float x50311, float x51111, float x52111);
float rhsx50200(float gto, float x00100, float x00200, float x01000, float x02000, float x30200, float x50000, float x50200, float x50210, float x50300, float x51200, float x52200);
float rhsx50210(float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x50010, float x50200, float x50210, float x50211, float x50310, float x51210, float x51211, float x52210, float x52211);
float rhsx50211(float x00011, float x00110, float x00210, float x01010, float x02010, float x50011, float x50210, float x50211, float x50311, float x51211, float x52211);
float rhsx50300(float gto, float x00100, float x00200, float x01000, float x02000, float x30300, float x50100, float x50200, float x50300, float x50310, float x51300, float x52300);
float rhsx50310(float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x50110, float x50210, float x50300, float x50310, float x50311, float x51310, float x51311, float x52310, float x52311);
float rhsx50311(float x00011, float x00110, float x00210, float x01010, float x02010, float x50111, float x50211, float x50310, float x50311, float x51311, float x52311);
float rhsx51000(float x00100, float x00200, float x01000, float x50000, float x51000, float x51010, float x51100, float x51200);
float rhsx51010(float x00011, float x00110, float x00210, float x01010, float x50010, float x51000, float x51010, float x51011, float x51110, float x51210);
float rhsx51011(float x00011, float x00110, float x00210, float x01010, float x01011, float x50010, float x50011, float x51010, float x51011, float x51111, float x51211);
float rhsx51100(float x00100, float x00200, float x01000, float x50100, float x51000, float x51100, float x51110, float x51300);
float rhsx51110(float x00011, float x00110, float x00210, float x01010, float x50110, float x51010, float x51100, float x51110, float x51111, float x51310);
float rhsx51111(float x00011, float x00110, float x00210, float x01010, float x01011, float x50110, float x50111, float x51011, float x51110, float x51111, float x51311);
float rhsx51200(float x00100, float x00200, float x01000, float x50200, float x51000, float x51200, float x51210, float x51300);
float rhsx51210(float x00011, float x00110, float x00210, float x01010, float x50210, float x51010, float x51200, float x51210, float x51211, float x51310);
float rhsx51211(float x00011, float x00110, float x00210, float x01010, float x01011, float x50210, float x50211, float x51011, float x51210, float x51211, float x51311);
float rhsx51300(float x00100, float x00200, float x01000, float x50300, float x51100, float x51200, float x51300, float x51310);
float rhsx51310(float x00011, float x00110, float x00210, float x01010, float x50310, float x51110, float x51210, float x51300, float x51310, float x51311);
float rhsx51311(float x00011, float x00110, float x00210, float x01010, float x01011, float x50310, float x50311, float x51111, float x51211, float x51310, float x51311);
float rhsx52000(float x00100, float x00200, float x02000, float x50000, float x52000, float x52010, float x52100, float x52200);
float rhsx52010(float x00011, float x00110, float x00210, float x02010, float x50010, float x52000, float x52010, float x52011, float x52110, float x52210);
float rhsx52011(float x00011, float x00110, float x00210, float x02010, float x02011, float x50010, float x50011, float x52010, float x52011, float x52111, float x52211);
float rhsx52100(float x00100, float x00200, float x02000, float x50100, float x52000, float x52100, float x52110, float x52300);
float rhsx52110(float x00011, float x00110, float x00210, float x02010, float x50110, float x52010, float x52100, float x52110, float x52111, float x52310);
float rhsx52111(float x00011, float x00110, float x00210, float x02010, float x02011, float x50110, float x50111, float x52011, float x52110, float x52111, float x52311);
float rhsx52200(float x00100, float x00200, float x02000, float x50200, float x52000, float x52200, float x52210, float x52300);
float rhsx52210(float x00011, float x00110, float x00210, float x02010, float x50210, float x52010, float x52200, float x52210, float x52211, float x52310);
float rhsx52211(float x00011, float x00110, float x00210, float x02010, float x02011, float x50210, float x50211, float x52011, float x52210, float x52211, float x52311);
float rhsx52300(float x00100, float x00200, float x02000, float x50300, float x52100, float x52200, float x52300, float x52310);
float rhsx52310(float x00011, float x00110, float x00210, float x02010, float x50310, float x52110, float x52210, float x52300, float x52310, float x52311);
float rhsx52311(float x00011, float x00110, float x00210, float x02010, float x02011, float x50310, float x50311, float x52111, float x52211, float x52310, float x52311);
float rhsx60000(float x00100, float x00200, float x01000, float x02000, float x60000, float x60010, float x60100, float x60200, float x61000, float x62000);
float rhsx60010(float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x60000, float x60010, float x60011, float x60110, float x60210, float x61010, float x61011, float x62010, float x62011);
float rhsx60011(float x00011, float x00110, float x00210, float x01010, float x02010, float x60010, float x60011, float x60111, float x60211, float x61011, float x62011);
float rhsx60100(float x00100, float x00200, float x01000, float x02000, float x50100, float x60000, float x60100, float x60110, float x60300, float x61100, float x62100);
float rhsx60110(float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x50110, float x60010, float x60100, float x60110, float x60111, float x60310, float x61110, float x61111, float x62110, float x62111);
float rhsx60111(float x00011, float x00110, float x00210, float x01010, float x02010, float x50111, float x60011, float x60110, float x60111, float x60311, float x61111, float x62111);
float rhsx60200(float gto, float x00100, float x00200, float x01000, float x02000, float x40200, float x60000, float x60200, float x60210, float x60300, float x61200, float x62200);
float rhsx60210(float gto, float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x40210, float x60010, float x60200, float x60210, float x60211, float x60310, float x61210, float x61211, float x62210, float x62211);
float rhsx60211(float gto, float x00011, float x00110, float x00210, float x01010, float x02010, float x40211, float x60011, float x60210, float x60211, float x60311, float x61211, float x62211);
float rhsx60300(float gto, float x00100, float x00200, float x01000, float x02000, float x40300, float x50300, float x60100, float x60200, float x60300, float x60310, float x61300, float x62300);
float rhsx60310(float gto, float x00011, float x00110, float x00210, float x01010, float x01011, float x02010, float x02011, float x40310, float x50310, float x60110, float x60210, float x60300, float x60310, float x60311, float x61310, float x61311, float x62310, float x62311);
float rhsx60311(float gto, float x00011, float x00110, float x00210, float x01010, float x02010, float x40311, float x50311, float x60111, float x60211, float x60310, float x60311, float x61311, float x62311);
float rhsx61000(float x00100, float x00200, float x01000, float x60000, float x61000, float x61010, float x61100, float x61200);
float rhsx61010(float x00011, float x00110, float x00210, float x01010, float x60010, float x61000, float x61010, float x61011, float x61110, float x61210);
float rhsx61011(float x00011, float x00110, float x00210, float x01010, float x01011, float x60010, float x60011, float x61010, float x61011, float x61111, float x61211);
float rhsx61100(float x00100, float x00200, float x01000, float x60100, float x61000, float x61100, float x61110, float x61300);
float rhsx61110(float x00011, float x00110, float x00210, float x01010, float x60110, float x61010, float x61100, float x61110, float x61111, float x61310);
float rhsx61111(float x00011, float x00110, float x00210, float x01010, float x01011, float x60110, float x60111, float x61011, float x61110, float x61111, float x61311);
float rhsx61200(float gto, float x00100, float x00200, float x01000, float x41200, float x60200, float x61000, float x61200, float x61210, float x61300);
float rhsx61210(float gto, float x00011, float x00110, float x00210, float x01010, float x41210, float x60210, float x61010, float x61200, float x61210, float x61211, float x61310);
float rhsx61211(float gto, float x00011, float x00110, float x00210, float x01010, float x01011, float x41211, float x60210, float x60211, float x61011, float x61210, float x61211, float x61311);
float rhsx61300(float gto, float x00100, float x00200, float x01000, float x41300, float x60300, float x61100, float x61200, float x61300, float x61310);
float rhsx61310(float gto, float x00011, float x00110, float x00210, float x01010, float x41310, float x60310, float x61110, float x61210, float x61300, float x61310, float x61311);
float rhsx61311(float gto, float x00011, float x00110, float x00210, float x01010, float x01011, float x41311, float x60310, float x60311, float x61111, float x61211, float x61310, float x61311);
float rhsx62000(float x00100, float x00200, float x02000, float x60000, float x62000, float x62010, float x62100, float x62200);
float rhsx62010(float x00011, float x00110, float x00210, float x02010, float x60010, float x62000, float x62010, float x62011, float x62110, float x62210);
float rhsx62011(float x00011, float x00110, float x00210, float x02010, float x02011, float x60010, float x60011, float x62010, float x62011, float x62111, float x62211);
float rhsx62100(float x00100, float x00200, float x02000, float x60100, float x62000, float x62100, float x62110, float x62300);
float rhsx62110(float x00011, float x00110, float x00210, float x02010, float x60110, float x62010, float x62100, float x62110, float x62111, float x62310);
float rhsx62111(float x00011, float x00110, float x00210, float x02010, float x02011, float x60110, float x60111, float x62011, float x62110, float x62111, float x62311);
float rhsx62200(float gto, float x00100, float x00200, float x02000, float x42200, float x60200, float x62000, float x62200, float x62210, float x62300);
float rhsx62210(float gto, float x00011, float x00110, float x00210, float x02010, float x42210, float x60210, float x62010, float x62200, float x62210, float x62211, float x62310);
float rhsx62211(float gto, float x00011, float x00110, float x00210, float x02010, float x02011, float x42211, float x60210, float x60211, float x62011, float x62210, float x62211, float x62311);
float rhsx62300(float gto, float x00100, float x00200, float x02000, float x42300, float x60300, float x62100, float x62200, float x62300, float x62310);
float rhsx62310(float gto, float x00011, float x00110, float x00210, float x02010, float x42310, float x60310, float x62110, float x62210, float x62300, float x62310, float x62311);
float rhsx62311(float gto, float x00011, float x00110, float x00210, float x02010, float x02011, float x42311, float x60310, float x60311, float x62111, float x62211, float x62310, float x62311);

#endif
