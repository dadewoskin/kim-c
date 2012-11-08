#ifndef _RHS_H_
#define _RHS_H_

#define step_size 0.0001 /* step size (in hours) */

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
double rhsGR(double GR, double G, double x01011, double x02011);
double rhsG(double GR, double G, double x00011);
double rhsGrR(double G, double GrR, double Gr, double x01011, double x02011);
double rhsGr(double G, double GrR, double Gr, double x00011);
double rhsGcR(double G, double GcR, double Gc, double x01011, double x02011);
double rhsGc(double G, double GcR, double Gc, double x00011);
double rhsGBR(double G, double GBR, double GB, double B, double revn, double revng, double revngp, double revnp);
double rhsGB(double G, double GBR, double GB, double B, double revn, double revng, double revngp, double revnp);
double rhsGBRb(double G, double GBR, double GB, double GBRb, double GBb, double B, double revn, double revng, double revngp, double revnp);
double rhsGBb(double G, double GBR, double GB, double GBRb, double GBb, double B, double revn, double revng, double revngp, double revnp);
double rhsMnPo(double G, double MnPo);
double rhsMcPo(double MnPo, double McPo);
double rhsMnPt(double G, double MnPt);
double rhsMcPt(double MnPt, double McPt);
double rhsMnRt(double G, double Gc, double MnRt);
double rhsMcRt(double MnRt, double McRt);
double rhsMnRev(double G, double Gr, double MnRev, double x00011);
double rhsMcRev(double MnRev, double McRev);
double rhsMnRo(double G, double GB, double MnRo, double B);
double rhsMcRo(double MnRo, double McRo);
double rhsMnB(double G, double GB, double GBb, double MnB, double B);
double rhsMcB(double MnB, double McB, double B);
double rhsMnNp(double G, double GB, double MnNp, double B);
double rhsMcNp(double MnNp, double McNp);
double rhsB(double McB, double B, double Cl, double BC);
double rhsCl(double McNp, double B, double Cl, double BC);
double rhsBC(double B, double Cl, double BC);
double rhscyrev(double McRev, double cyrev, double revn, double cyrevg, double x00200);
double rhsrevn(double cyrev, double revn, double revng, double x00210);
double rhscyrevg(double cyrev, double revn, double cyrevg, double revng, double gto, double x00200);
double rhsrevng(double cyrev, double revn, double cyrevg, double revng, double gto, double x00210);
double rhscyrevgp(double cyrev, double revn, double cyrevg, double revng, double cyrevgp, double revngp, double gto);
double rhsrevngp(double cyrev, double revn, double cyrevg, double revng, double cyrevgp, double revngp, double gto);
double rhscyrevp(double cyrev, double revn, double cyrevg, double cyrevgp, double cyrevp, double revnp);
double rhsrevnp(double cyrev, double revn, double revng, double revngp, double cyrevp, double revnp);
double rhsgto(double G, double GB, double B, double gto);
double rhsx00001(double B, double BC, double x00001);
double rhsx00011(double x00001, double x00011, double x01010, double x01011, double x02010, double x02011, double x20010, double x20011, double x20110, double x20111, double x21010, double x21011, double x21110, double x21111, double x22010, double x22011, double x22110, double x22111, double x40010, double x40011, double x40110, double x40111, double x40210, double x40211, double x40310, double x40311, double x41010, double x41011, double x41110, double x41111, double x41210, double x41211, double x41310, double x41311, double x42010, double x42011, double x42110, double x42111, double x42210, double x42211, double x42310, double x42311, double x50010, double x50011, double x50110, double x50111, double x50210, double x50211, double x50310, double x50311, double x51010, double x51011, double x51110, double x51111, double x51210, double x51211, double x51310, double x51311, double x52010, double x52011, double x52110, double x52111, double x52210, double x52211, double x52310, double x52311, double x60010, double x60011, double x60110, double x60111, double x60210, double x60211, double x60310, double x60311, double x61010, double x61011, double x61110, double x61111, double x61210, double x61211, double x61310, double x61311, double x62010, double x62011, double x62110, double x62111, double x62210, double x62211, double x62310, double x62311);
double rhsx00100(double x00100, double x00110, double x10000, double x10100, double x20000, double x20100, double x21000, double x21100, double x22000, double x22100, double x30000, double x30100, double x30200, double x30300, double x40000, double x40100, double x40200, double x40300, double x41000, double x41100, double x41200, double x41300, double x42000, double x42100, double x42200, double x42300, double x50000, double x50100, double x50200, double x50300, double x51000, double x51100, double x51200, double x51300, double x52000, double x52100, double x52200, double x52300, double x60000, double x60100, double x60200, double x60300, double x61000, double x61100, double x61200, double x61300, double x62000, double x62100, double x62200, double x62300);
double rhsx00110(double x00110, double x20010, double x20011, double x20110, double x20111, double x21010, double x21011, double x21110, double x21111, double x22010, double x22011, double x22110, double x22111, double x40010, double x40011, double x40110, double x40111, double x40210, double x40211, double x40310, double x40311, double x41010, double x41011, double x41110, double x41111, double x41210, double x41211, double x41310, double x41311, double x42010, double x42011, double x42110, double x42111, double x42210, double x42211, double x42310, double x42311, double x50010, double x50011, double x50110, double x50111, double x50210, double x50211, double x50310, double x50311, double x51010, double x51011, double x51110, double x51111, double x51210, double x51211, double x51310, double x51311, double x52010, double x52011, double x52110, double x52111, double x52210, double x52211, double x52310, double x52311, double x60010, double x60011, double x60110, double x60111, double x60210, double x60211, double x60310, double x60311, double x61010, double x61011, double x61110, double x61111, double x61210, double x61211, double x61310, double x61311, double x62010, double x62011, double x62110, double x62111, double x62210, double x62211, double x62310, double x62311);
double rhsx00200(double cyrev, double cyrevg, double cyrevgp, double x00200, double x00210, double x30000, double x30100, double x30200, double x30300, double x40000, double x40100, double x40200, double x40300, double x41000, double x41100, double x41200, double x41300, double x42000, double x42100, double x42200, double x42300, double x50000, double x50100, double x50200, double x50300, double x51000, double x51100, double x51200, double x51300, double x52000, double x52100, double x52200, double x52300, double x60000, double x60100, double x60200, double x60300, double x61000, double x61100, double x61200, double x61300, double x62000, double x62100, double x62200, double x62300);
double rhsx00210(double revn, double revng, double revngp, double x00210, double x40010, double x40011, double x40110, double x40111, double x40210, double x40211, double x40310, double x40311, double x41010, double x41011, double x41110, double x41111, double x41210, double x41211, double x41310, double x41311, double x42010, double x42011, double x42110, double x42111, double x42210, double x42211, double x42310, double x42311, double x50010, double x50011, double x50110, double x50111, double x50210, double x50211, double x50310, double x50311, double x51010, double x51011, double x51110, double x51111, double x51210, double x51211, double x51310, double x51311, double x52010, double x52011, double x52110, double x52111, double x52210, double x52211, double x52310, double x52311, double x60010, double x60011, double x60110, double x60111, double x60210, double x60211, double x60310, double x60311, double x61010, double x61011, double x61110, double x61111, double x61210, double x61211, double x61310, double x61311, double x62010, double x62011, double x62110, double x62111, double x62210, double x62211, double x62310, double x62311);
double rhsx01000(double McRo, double x01000, double x20000, double x20100, double x21000, double x21100, double x40000, double x40100, double x40200, double x40300, double x41000, double x41100, double x41200, double x41300, double x50000, double x50100, double x50200, double x50300, double x51000, double x51100, double x51200, double x51300, double x60000, double x60100, double x60200, double x60300, double x61000, double x61100, double x61200, double x61300);
double rhsx01010(double x00011, double x01010, double x01011, double x20010, double x20011, double x20110, double x20111, double x21010, double x21011, double x21110, double x21111, double x40010, double x40011, double x40110, double x40111, double x40210, double x40211, double x40310, double x40311, double x41010, double x41011, double x41110, double x41111, double x41210, double x41211, double x41310, double x41311, double x50010, double x50011, double x50110, double x50111, double x50210, double x50211, double x50310, double x50311, double x51010, double x51011, double x51110, double x51111, double x51210, double x51211, double x51310, double x51311, double x60010, double x60011, double x60110, double x60111, double x60210, double x60211, double x60310, double x60311, double x61010, double x61011, double x61110, double x61111, double x61210, double x61211, double x61310, double x61311);
double rhsx01011(double x00011, double x01010, double x01011, double x20010, double x20110, double x21011, double x21111, double x40010, double x40110, double x40210, double x40310, double x41011, double x41111, double x41211, double x41311, double x50010, double x50110, double x50210, double x50310, double x51011, double x51111, double x51211, double x51311, double x60010, double x60110, double x60210, double x60310, double x61011, double x61111, double x61211, double x61311);
double rhsx02000(double McRt, double x02000, double x20000, double x20100, double x22000, double x22100, double x40000, double x40100, double x40200, double x40300, double x42000, double x42100, double x42200, double x42300, double x50000, double x50100, double x50200, double x50300, double x52000, double x52100, double x52200, double x52300, double x60000, double x60100, double x60200, double x60300, double x62000, double x62100, double x62200, double x62300);
double rhsx02010(double x00011, double x02010, double x02011, double x20010, double x20011, double x20110, double x20111, double x22010, double x22011, double x22110, double x22111, double x40010, double x40011, double x40110, double x40111, double x40210, double x40211, double x40310, double x40311, double x42010, double x42011, double x42110, double x42111, double x42210, double x42211, double x42310, double x42311, double x50010, double x50011, double x50110, double x50111, double x50210, double x50211, double x50310, double x50311, double x52010, double x52011, double x52110, double x52111, double x52210, double x52211, double x52310, double x52311, double x60010, double x60011, double x60110, double x60111, double x60210, double x60211, double x60310, double x60311, double x62010, double x62011, double x62110, double x62111, double x62210, double x62211, double x62310, double x62311);
double rhsx02011(double x00011, double x02010, double x02011, double x20010, double x20110, double x22011, double x22111, double x40010, double x40110, double x40210, double x40310, double x42011, double x42111, double x42211, double x42311, double x50010, double x50110, double x50210, double x50310, double x52011, double x52111, double x52211, double x52311, double x60010, double x60110, double x60210, double x60310, double x62011, double x62111, double x62211, double x62311);
double rhsx10000(double McPo, double x00100, double x10000, double x10100);
double rhsx10100(double x00100, double x10000, double x10100);
double rhsx20000(double x00100, double x01000, double x02000, double x20000, double x20010, double x20100, double x21000, double x22000);
double rhsx20010(double x00011, double x00110, double x01010, double x01011, double x02010, double x02011, double x20000, double x20010, double x20011, double x20110, double x21010, double x21011, double x22010, double x22011);
double rhsx20011(double x00011, double x00110, double x01010, double x02010, double x20010, double x20011, double x20111, double x21011, double x22011);
double rhsx20100(double x00100, double x01000, double x02000, double x10100, double x20000, double x20100, double x20110, double x21100, double x22100);
double rhsx20110(double x00011, double x00110, double x01010, double x01011, double x02010, double x02011, double x20010, double x20100, double x20110, double x20111, double x21110, double x21111, double x22110, double x22111);
double rhsx20111(double x00011, double x00110, double x01010, double x02010, double x20011, double x20110, double x20111, double x21111, double x22111);
double rhsx21000(double x00100, double x01000, double x20000, double x21000, double x21010, double x21100);
double rhsx21010(double x00011, double x00110, double x01010, double x20010, double x21000, double x21010, double x21011, double x21110);
double rhsx21011(double x00011, double x00110, double x01010, double x01011, double x20010, double x20011, double x21010, double x21011, double x21111);
double rhsx21100(double x00100, double x01000, double x20100, double x21000, double x21100, double x21110);
double rhsx21110(double x00011, double x00110, double x01010, double x20110, double x21010, double x21100, double x21110, double x21111);
double rhsx21111(double x00011, double x00110, double x01010, double x01011, double x20110, double x20111, double x21011, double x21110, double x21111);
double rhsx22000(double x00100, double x02000, double x20000, double x22000, double x22010, double x22100);
double rhsx22010(double x00011, double x00110, double x02010, double x20010, double x22000, double x22010, double x22011, double x22110);
double rhsx22011(double x00011, double x00110, double x02010, double x02011, double x20010, double x20011, double x22010, double x22011, double x22111);
double rhsx22100(double x00100, double x02000, double x20100, double x22000, double x22100, double x22110);
double rhsx22110(double x00011, double x00110, double x02010, double x20110, double x22010, double x22100, double x22110, double x22111);
double rhsx22111(double x00011, double x00110, double x02010, double x02011, double x20110, double x20111, double x22011, double x22110, double x22111);
double rhsx30000(double McPt, double x00100, double x00200, double x30000, double x30100, double x30200);
double rhsx30100(double x00100, double x00200, double x30000, double x30100, double x30300);
double rhsx30200(double gto, double x00100, double x00200, double x30000, double x30200, double x30300);
double rhsx30300(double gto, double x00100, double x00200, double x30100, double x30200, double x30300);
double rhsx40000(double x00100, double x00200, double x01000, double x02000, double x40000, double x40010, double x40100, double x40200, double x41000, double x42000);
double rhsx40010(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40000, double x40010, double x40011, double x40110, double x40210, double x41010, double x41011, double x42010, double x42011);
double rhsx40011(double x00011, double x00110, double x00210, double x01010, double x02010, double x40010, double x40011, double x40111, double x40211, double x41011, double x42011);
double rhsx40100(double x00100, double x00200, double x01000, double x02000, double x30100, double x40000, double x40100, double x40110, double x40300, double x41100, double x42100);
double rhsx40110(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40010, double x40100, double x40110, double x40111, double x40310, double x41110, double x41111, double x42110, double x42111);
double rhsx40111(double x00011, double x00110, double x00210, double x01010, double x02010, double x40011, double x40110, double x40111, double x40311, double x41111, double x42111);
double rhsx40200(double gto, double x00100, double x00200, double x01000, double x02000, double x40000, double x40200, double x40210, double x40300, double x41200, double x42200);
double rhsx40210(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40010, double x40200, double x40210, double x40211, double x40310, double x41210, double x41211, double x42210, double x42211);
double rhsx40211(double gto, double x00011, double x00110, double x00210, double x01010, double x02010, double x40011, double x40210, double x40211, double x40311, double x41211, double x42211);
double rhsx40300(double gto, double x00100, double x00200, double x01000, double x02000, double x30300, double x40100, double x40200, double x40300, double x40310, double x41300, double x42300);
double rhsx40310(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40110, double x40210, double x40300, double x40310, double x40311, double x41310, double x41311, double x42310, double x42311);
double rhsx40311(double gto, double x00011, double x00110, double x00210, double x01010, double x02010, double x40111, double x40211, double x40310, double x40311, double x41311, double x42311);
double rhsx41000(double x00100, double x00200, double x01000, double x40000, double x41000, double x41010, double x41100, double x41200);
double rhsx41010(double x00011, double x00110, double x00210, double x01010, double x40010, double x41000, double x41010, double x41011, double x41110, double x41210);
double rhsx41011(double x00011, double x00110, double x00210, double x01010, double x01011, double x40010, double x40011, double x41010, double x41011, double x41111, double x41211);
double rhsx41100(double x00100, double x00200, double x01000, double x40100, double x41000, double x41100, double x41110, double x41300);
double rhsx41110(double x00011, double x00110, double x00210, double x01010, double x40110, double x41010, double x41100, double x41110, double x41111, double x41310);
double rhsx41111(double x00011, double x00110, double x00210, double x01010, double x01011, double x40110, double x40111, double x41011, double x41110, double x41111, double x41311);
double rhsx41200(double gto, double x00100, double x00200, double x01000, double x40200, double x41000, double x41200, double x41210, double x41300);
double rhsx41210(double gto, double x00011, double x00110, double x00210, double x01010, double x40210, double x41010, double x41200, double x41210, double x41211, double x41310);
double rhsx41211(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x40210, double x40211, double x41011, double x41210, double x41211, double x41311);
double rhsx41300(double gto, double x00100, double x00200, double x01000, double x40300, double x41100, double x41200, double x41300, double x41310);
double rhsx41310(double gto, double x00011, double x00110, double x00210, double x01010, double x40310, double x41110, double x41210, double x41300, double x41310, double x41311);
double rhsx41311(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x40310, double x40311, double x41111, double x41211, double x41310, double x41311);
double rhsx42000(double x00100, double x00200, double x02000, double x40000, double x42000, double x42010, double x42100, double x42200);
double rhsx42010(double x00011, double x00110, double x00210, double x02010, double x40010, double x42000, double x42010, double x42011, double x42110, double x42210);
double rhsx42011(double x00011, double x00110, double x00210, double x02010, double x02011, double x40010, double x40011, double x42010, double x42011, double x42111, double x42211);
double rhsx42100(double x00100, double x00200, double x02000, double x40100, double x42000, double x42100, double x42110, double x42300);
double rhsx42110(double x00011, double x00110, double x00210, double x02010, double x40110, double x42010, double x42100, double x42110, double x42111, double x42310);
double rhsx42111(double x00011, double x00110, double x00210, double x02010, double x02011, double x40110, double x40111, double x42011, double x42110, double x42111, double x42311);
double rhsx42200(double gto, double x00100, double x00200, double x02000, double x40200, double x42000, double x42200, double x42210, double x42300);
double rhsx42210(double gto, double x00011, double x00110, double x00210, double x02010, double x40210, double x42010, double x42200, double x42210, double x42211, double x42310);
double rhsx42211(double gto, double x00011, double x00110, double x00210, double x02010, double x02011, double x40210, double x40211, double x42011, double x42210, double x42211, double x42311);
double rhsx42300(double gto, double x00100, double x00200, double x02000, double x40300, double x42100, double x42200, double x42300, double x42310);
double rhsx42310(double gto, double x00011, double x00110, double x00210, double x02010, double x40310, double x42110, double x42210, double x42300, double x42310, double x42311);
double rhsx42311(double gto, double x00011, double x00110, double x00210, double x02010, double x02011, double x40310, double x40311, double x42111, double x42211, double x42310, double x42311);
double rhsx50000(double x00100, double x00200, double x01000, double x02000, double x50000, double x50010, double x50100, double x50200, double x51000, double x52000);
double rhsx50010(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x50000, double x50010, double x50011, double x50110, double x50210, double x51010, double x51011, double x52010, double x52011);
double rhsx50011(double x00011, double x00110, double x00210, double x01010, double x02010, double x50010, double x50011, double x50111, double x50211, double x51011, double x52011);
double rhsx50100(double x00100, double x00200, double x01000, double x02000, double x50000, double x50100, double x50110, double x50300, double x51100, double x52100);
double rhsx50110(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x50010, double x50100, double x50110, double x50111, double x50310, double x51110, double x51111, double x52110, double x52111);
double rhsx50111(double x00011, double x00110, double x00210, double x01010, double x02010, double x50011, double x50110, double x50111, double x50311, double x51111, double x52111);
double rhsx50200(double gto, double x00100, double x00200, double x01000, double x02000, double x30200, double x50000, double x50200, double x50210, double x50300, double x51200, double x52200);
double rhsx50210(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x50010, double x50200, double x50210, double x50211, double x50310, double x51210, double x51211, double x52210, double x52211);
double rhsx50211(double x00011, double x00110, double x00210, double x01010, double x02010, double x50011, double x50210, double x50211, double x50311, double x51211, double x52211);
double rhsx50300(double gto, double x00100, double x00200, double x01000, double x02000, double x30300, double x50100, double x50200, double x50300, double x50310, double x51300, double x52300);
double rhsx50310(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x50110, double x50210, double x50300, double x50310, double x50311, double x51310, double x51311, double x52310, double x52311);
double rhsx50311(double x00011, double x00110, double x00210, double x01010, double x02010, double x50111, double x50211, double x50310, double x50311, double x51311, double x52311);
double rhsx51000(double x00100, double x00200, double x01000, double x50000, double x51000, double x51010, double x51100, double x51200);
double rhsx51010(double x00011, double x00110, double x00210, double x01010, double x50010, double x51000, double x51010, double x51011, double x51110, double x51210);
double rhsx51011(double x00011, double x00110, double x00210, double x01010, double x01011, double x50010, double x50011, double x51010, double x51011, double x51111, double x51211);
double rhsx51100(double x00100, double x00200, double x01000, double x50100, double x51000, double x51100, double x51110, double x51300);
double rhsx51110(double x00011, double x00110, double x00210, double x01010, double x50110, double x51010, double x51100, double x51110, double x51111, double x51310);
double rhsx51111(double x00011, double x00110, double x00210, double x01010, double x01011, double x50110, double x50111, double x51011, double x51110, double x51111, double x51311);
double rhsx51200(double x00100, double x00200, double x01000, double x50200, double x51000, double x51200, double x51210, double x51300);
double rhsx51210(double x00011, double x00110, double x00210, double x01010, double x50210, double x51010, double x51200, double x51210, double x51211, double x51310);
double rhsx51211(double x00011, double x00110, double x00210, double x01010, double x01011, double x50210, double x50211, double x51011, double x51210, double x51211, double x51311);
double rhsx51300(double x00100, double x00200, double x01000, double x50300, double x51100, double x51200, double x51300, double x51310);
double rhsx51310(double x00011, double x00110, double x00210, double x01010, double x50310, double x51110, double x51210, double x51300, double x51310, double x51311);
double rhsx51311(double x00011, double x00110, double x00210, double x01010, double x01011, double x50310, double x50311, double x51111, double x51211, double x51310, double x51311);
double rhsx52000(double x00100, double x00200, double x02000, double x50000, double x52000, double x52010, double x52100, double x52200);
double rhsx52010(double x00011, double x00110, double x00210, double x02010, double x50010, double x52000, double x52010, double x52011, double x52110, double x52210);
double rhsx52011(double x00011, double x00110, double x00210, double x02010, double x02011, double x50010, double x50011, double x52010, double x52011, double x52111, double x52211);
double rhsx52100(double x00100, double x00200, double x02000, double x50100, double x52000, double x52100, double x52110, double x52300);
double rhsx52110(double x00011, double x00110, double x00210, double x02010, double x50110, double x52010, double x52100, double x52110, double x52111, double x52310);
double rhsx52111(double x00011, double x00110, double x00210, double x02010, double x02011, double x50110, double x50111, double x52011, double x52110, double x52111, double x52311);
double rhsx52200(double x00100, double x00200, double x02000, double x50200, double x52000, double x52200, double x52210, double x52300);
double rhsx52210(double x00011, double x00110, double x00210, double x02010, double x50210, double x52010, double x52200, double x52210, double x52211, double x52310);
double rhsx52211(double x00011, double x00110, double x00210, double x02010, double x02011, double x50210, double x50211, double x52011, double x52210, double x52211, double x52311);
double rhsx52300(double x00100, double x00200, double x02000, double x50300, double x52100, double x52200, double x52300, double x52310);
double rhsx52310(double x00011, double x00110, double x00210, double x02010, double x50310, double x52110, double x52210, double x52300, double x52310, double x52311);
double rhsx52311(double x00011, double x00110, double x00210, double x02010, double x02011, double x50310, double x50311, double x52111, double x52211, double x52310, double x52311);
double rhsx60000(double x00100, double x00200, double x01000, double x02000, double x60000, double x60010, double x60100, double x60200, double x61000, double x62000);
double rhsx60010(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x60000, double x60010, double x60011, double x60110, double x60210, double x61010, double x61011, double x62010, double x62011);
double rhsx60011(double x00011, double x00110, double x00210, double x01010, double x02010, double x60010, double x60011, double x60111, double x60211, double x61011, double x62011);
double rhsx60100(double x00100, double x00200, double x01000, double x02000, double x50100, double x60000, double x60100, double x60110, double x60300, double x61100, double x62100);
double rhsx60110(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x50110, double x60010, double x60100, double x60110, double x60111, double x60310, double x61110, double x61111, double x62110, double x62111);
double rhsx60111(double x00011, double x00110, double x00210, double x01010, double x02010, double x50111, double x60011, double x60110, double x60111, double x60311, double x61111, double x62111);
double rhsx60200(double gto, double x00100, double x00200, double x01000, double x02000, double x40200, double x60000, double x60200, double x60210, double x60300, double x61200, double x62200);
double rhsx60210(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40210, double x60010, double x60200, double x60210, double x60211, double x60310, double x61210, double x61211, double x62210, double x62211);
double rhsx60211(double gto, double x00011, double x00110, double x00210, double x01010, double x02010, double x40211, double x60011, double x60210, double x60211, double x60311, double x61211, double x62211);
double rhsx60300(double gto, double x00100, double x00200, double x01000, double x02000, double x40300, double x50300, double x60100, double x60200, double x60300, double x60310, double x61300, double x62300);
double rhsx60310(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40310, double x50310, double x60110, double x60210, double x60300, double x60310, double x60311, double x61310, double x61311, double x62310, double x62311);
double rhsx60311(double gto, double x00011, double x00110, double x00210, double x01010, double x02010, double x40311, double x50311, double x60111, double x60211, double x60310, double x60311, double x61311, double x62311);
double rhsx61000(double x00100, double x00200, double x01000, double x60000, double x61000, double x61010, double x61100, double x61200);
double rhsx61010(double x00011, double x00110, double x00210, double x01010, double x60010, double x61000, double x61010, double x61011, double x61110, double x61210);
double rhsx61011(double x00011, double x00110, double x00210, double x01010, double x01011, double x60010, double x60011, double x61010, double x61011, double x61111, double x61211);
double rhsx61100(double x00100, double x00200, double x01000, double x60100, double x61000, double x61100, double x61110, double x61300);
double rhsx61110(double x00011, double x00110, double x00210, double x01010, double x60110, double x61010, double x61100, double x61110, double x61111, double x61310);
double rhsx61111(double x00011, double x00110, double x00210, double x01010, double x01011, double x60110, double x60111, double x61011, double x61110, double x61111, double x61311);
double rhsx61200(double gto, double x00100, double x00200, double x01000, double x41200, double x60200, double x61000, double x61200, double x61210, double x61300);
double rhsx61210(double gto, double x00011, double x00110, double x00210, double x01010, double x41210, double x60210, double x61010, double x61200, double x61210, double x61211, double x61310);
double rhsx61211(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x41211, double x60210, double x60211, double x61011, double x61210, double x61211, double x61311);
double rhsx61300(double gto, double x00100, double x00200, double x01000, double x41300, double x60300, double x61100, double x61200, double x61300, double x61310);
double rhsx61310(double gto, double x00011, double x00110, double x00210, double x01010, double x41310, double x60310, double x61110, double x61210, double x61300, double x61310, double x61311);
double rhsx61311(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x41311, double x60310, double x60311, double x61111, double x61211, double x61310, double x61311);
double rhsx62000(double x00100, double x00200, double x02000, double x60000, double x62000, double x62010, double x62100, double x62200);
double rhsx62010(double x00011, double x00110, double x00210, double x02010, double x60010, double x62000, double x62010, double x62011, double x62110, double x62210);
double rhsx62011(double x00011, double x00110, double x00210, double x02010, double x02011, double x60010, double x60011, double x62010, double x62011, double x62111, double x62211);
double rhsx62100(double x00100, double x00200, double x02000, double x60100, double x62000, double x62100, double x62110, double x62300);
double rhsx62110(double x00011, double x00110, double x00210, double x02010, double x60110, double x62010, double x62100, double x62110, double x62111, double x62310);
double rhsx62111(double x00011, double x00110, double x00210, double x02010, double x02011, double x60110, double x60111, double x62011, double x62110, double x62111, double x62311);
double rhsx62200(double gto, double x00100, double x00200, double x02000, double x42200, double x60200, double x62000, double x62200, double x62210, double x62300);
double rhsx62210(double gto, double x00011, double x00110, double x00210, double x02010, double x42210, double x60210, double x62010, double x62200, double x62210, double x62211, double x62310);
double rhsx62211(double gto, double x00011, double x00110, double x00210, double x02010, double x02011, double x42211, double x60210, double x60211, double x62011, double x62210, double x62211, double x62311);
double rhsx62300(double gto, double x00100, double x00200, double x02000, double x42300, double x60300, double x62100, double x62200, double x62300, double x62310);
double rhsx62310(double gto, double x00011, double x00110, double x00210, double x02010, double x42310, double x60310, double x62110, double x62210, double x62300, double x62310, double x62311);
double rhsx62311(double gto, double x00011, double x00110, double x00210, double x02010, double x02011, double x42311, double x60310, double x60311, double x62111, double x62211, double x62310, double x62311);

#endif
