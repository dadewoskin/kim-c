/* simulate one cell using the equations from Kim 2012 */
#include <stdio.h> 
#include <math.h> 
#include <stdlib.h>
#include "rhs.h"
//#include <string.h> 

double rhsGR(double GR, double G, double x01011, double x02011)
{
	return -(unbin*GR)+bin*(1-G-GR)*(x01011+x02011);
};
double rhsG(double GR, double G, double x00011)
{
	return -(unbin*G)+bin*(1-G-GR)*x00011;
};
double rhsGrR(double G, double GrR, double Gr, double x01011, double x02011)
{
	return -(unbinr*GrR)+binr*(1-Gr-GrR)*(x01011+x02011);
};
double rhsGr(double G, double GrR, double Gr, double x00011)
{
	return -(unbinr*Gr)+binr*(1-Gr-GrR)*x00011;
};
double rhsGcR(double G, double GcR, double Gc, double x01011, double x02011)
{
	return -(unbinc*GcR)+binc*(1-Gc-GcR)*(x01011+x02011);
};
double rhsGc(double G, double GcR, double Gc, double x00011)
{
	return -(unbinc*Gc)+binc*(1-Gc-GcR)*x00011;
};
double rhsGBR(double G, double GBR, double GB, double B, double revn, double revng, double revngp, double revnp)
{
	return -(unbinrev*GBR)+binrev*GB*(revn+revng+revngp+revnp);
};
double rhsGB(double G, double GBR, double GB, double B, double revn, double revng, double revngp, double revnp)
{
	return unbinrev*GBR-binrev*GB*(revn+revng+revngp+revnp);
};
double rhsGBRb(double G, double GBR, double GB, double GBRb, double GBb, double B, double revn, double revng, double revngp, double revnp)
{
	return -(unbinrevb*GBRb)+binrevb*GBb*(revn+revng+revngp+revnp);
};
double rhsGBb(double G, double GBR, double GB, double GBRb, double GBb, double B, double revn, double revng, double revngp, double revnp)
{
	return unbinrevb*GBRb-binrevb*GBb*(revn+revng+revngp+revnp);
};
double rhsMnPo(double G, double MnPo)
{
	return trPo*G-tmc*MnPo-umPo*MnPo;
};
double rhsMcPo(double MnPo, double McPo)
{
	return -(umPo*McPo)+tmc*MnPo;
};
double rhsMnPt(double G, double MnPt)
{
	return trPt*G-tmc*MnPt-umPt*MnPt;
};
double rhsMcPt(double MnPt, double McPt)
{
	return -(umPt*McPt)+tmc*MnPt;
};
double rhsMnRt(double G, double Gc, double MnRt)
{
	return trRt*Gc-tmc*MnRt-umRt*MnRt;
};
double rhsMcRt(double MnRt, double McRt)
{
	return -(umRt*McRt)+tmc*MnRt;
};
double rhsMnRev(double G, double Gr, double MnRev, double x00011)
{
	return -(tmcrev*MnRev)-umRev*MnRev+trRev*Gr*x00011;
};
double rhsMcRev(double MnRev, double McRev)
{
	return -(umRev*McRev)+tmcrev*MnRev;
};
double rhsMnRo(double G, double GB, double MnRo, double B)
{
	return trRo*G*GB-tmc*MnRo-umRo*MnRo;
};
double rhsMcRo(double MnRo, double McRo)
{
	return -(umRo*McRo)+tmc*MnRo;
};
double rhsMnB(double G, double GB, double GBb, double MnB, double B)
{
	return trB*GBb-tmc*MnB-umB*MnB;
};
double rhsMcB(double MnB, double McB, double B)
{
	return -(umB*McB)+tmc*MnB;
};
double rhsMnNp(double G, double GB, double MnNp, double B)
{
	return trNp*GB-tmc*MnNp-umNp*MnNp;
};
double rhsMcNp(double MnNp, double McNp)
{
	return -(umNp*McNp)+tmc*MnNp;
};
double rhsB(double McB, double B, double Cl, double BC)
{
	return -(ub*B)+uncbin*BC-cbin*B*Cl+tlb*McB;
};
double rhsCl(double McNp, double B, double Cl, double BC)
{
	return tlc+uncbin*BC-uc*Cl-cbin*B*Cl+tlnp*McNp;
};
double rhsBC(double B, double Cl, double BC)
{
	return -(phos*BC)-ubc*BC-uncbin*BC+cbin*B*Cl;
};
double rhscyrev(double McRev, double cyrev, double revn, double cyrevg, double x00200)
{
	return -((nlrev+urev)*cyrev)+dg*cyrevg+tlrev*McRev+nerev*revn-ag*cyrev*x00200;
};
double rhsrevn(double cyrev, double revn, double revng, double x00210)
{
	return nlrev*cyrev+(-nerev-urev)*revn+dg*revng-ag*Nf*revn*x00210;
};
double rhscyrevg(double cyrev, double revn, double cyrevg, double revng, double gto, double x00200)
{
	return -(cyrevg*(dg+nlrev+urev+gto))+nerev*revng+ag*cyrev*x00200;
};
double rhsrevng(double cyrev, double revn, double cyrevg, double revng, double gto, double x00210)
{
	return nlrev*cyrevg-(dg+nerev+urev+gto)*revng+ag*Nf*revn*x00210;
};
double rhscyrevgp(double cyrev, double revn, double cyrevg, double revng, double cyrevgp, double revngp, double gto)
{
	return -((dg+nlrev+uprev)*cyrevgp)+cyrevg*gto+nerev*revngp;
};
double rhsrevngp(double cyrev, double revn, double cyrevg, double revng, double cyrevgp, double revngp, double gto)
{
	return nlrev*cyrevgp+gto*revng-(dg+nerev+uprev)*revngp;
};
double rhscyrevp(double cyrev, double revn, double cyrevg, double cyrevgp, double cyrevp, double revnp)
{
	return dg*cyrevgp-(nlrev+uprev)*cyrevp+nerev*revnp;
};
double rhsrevnp(double cyrev, double revn, double revng, double revngp, double cyrevp, double revnp)
{
	return nlrev*cyrevp+dg*revngp-(nerev+uprev)*revnp;
};
double rhsgto(double G, double GB, double B, double gto)
{
	return trgto*G*GB-ugto*gto;
};
double rhsx00001(double B, double BC, double x00001)
{
	return phos*BC-nlbc*x00001-ubc*x00001;
};
double rhsx00011(double x00001, double x00011, double x01010, double x01011, double x02010, double x02011, double x20010, double x20011, double x20110, double x20111, double x21010, double x21011, double x21110, double x21111, double x22010, double x22011, double x22110, double x22111, double x40010, double x40011, double x40110, double x40111, double x40210, double x40211, double x40310, double x40311, double x41010, double x41011, double x41110, double x41111, double x41210, double x41211, double x41310, double x41311, double x42010, double x42011, double x42110, double x42111, double x42210, double x42211, double x42310, double x42311, double x50010, double x50011, double x50110, double x50111, double x50210, double x50211, double x50310, double x50311, double x51010, double x51011, double x51110, double x51111, double x51210, double x51211, double x51310, double x51311, double x52010, double x52011, double x52110, double x52111, double x52210, double x52211, double x52310, double x52311, double x60010, double x60011, double x60110, double x60111, double x60210, double x60211, double x60310, double x60311, double x61010, double x61011, double x61110, double x61111, double x61210, double x61211, double x61310, double x61311, double x62010, double x62011, double x62110, double x62111, double x62210, double x62211, double x62310, double x62311)
{
	return nlbc*x00001-ubc*x00011+uro*x01011-cbbin*Nf*x00011*(x01010+x02010)+urt*x02011+uncbbin*(x01011+x02011)+upu*(x50011+x50111+x50211+x50311)+up*(x20011+x20111+x40011+x40111+x40211+x40311+x60011+x60111+x60211+x60311)-bbin*Nf*x00011*(x20010+x20110+x21010+x21110+x22010+x22110+x40010+x40110+x40210+x40310+x41010+x41110+x41210+x41310+x42010+x42110+x42210+x42310+x50010+x50110+x50210+x50310+x51010+x51110+x51210+x51310+x52010+x52110+x52210+x52310+x60010+x60110+x60210+x60310+x61010+x61110+x61210+x61310+x62010+x62110+x62210+x62310)+unbbin*(x20011+x20111+x21011+x21111+x22011+x22111+x40011+x40111+x40211+x40311+x41011+x41111+x41211+x41311+x42011+x42111+x42211+x42311+x50011+x50111+x50211+x50311+x51011+x51111+x51211+x51311+x52011+x52111+x52211+x52311+x60011+x60111+x60211+x60311+x61011+x61111+x61211+x61311+x62011+x62111+x62211+x62311);
};
double rhsx00100(double x00100, double x00110, double x10000, double x10100, double x20000, double x20100, double x21000, double x21100, double x22000, double x22100, double x30000, double x30100, double x30200, double x30300, double x40000, double x40100, double x40200, double x40300, double x41000, double x41100, double x41200, double x41300, double x42000, double x42100, double x42200, double x42300, double x50000, double x50100, double x50200, double x50300, double x51000, double x51100, double x51200, double x51300, double x52000, double x52100, double x52200, double x52300, double x60000, double x60100, double x60200, double x60300, double x61000, double x61100, double x61200, double x61300, double x62000, double x62100, double x62200, double x62300)
{
	return lne*x00110+upu*(x10100+x30100+x30300+x50100+x50300)+up*(x20100+x40100+x40300+x60100+x60300)-ac*x00100*(x10000+x20000+x21000+x22000+x30000+x40000+x41000+x42000+x50000+x51000+x52000+x60000+x61000+x62000)+dc*(x10100+x20100+x21100+x22100+x30100+x40100+x41100+x42100+x50100+x51100+x52100+x60100+x61100+x62100)-ac*x00100*(x30200+x40200+x41200+x42200+x50200+x51200+x52200+x60200+x61200+x62200)+dc*(x30300+x40300+x41300+x42300+x50300+x51300+x52300+x60300+x61300+x62300);
};
double rhsx00110(double x00110, double x20010, double x20011, double x20110, double x20111, double x21010, double x21011, double x21110, double x21111, double x22010, double x22011, double x22110, double x22111, double x40010, double x40011, double x40110, double x40111, double x40210, double x40211, double x40310, double x40311, double x41010, double x41011, double x41110, double x41111, double x41210, double x41211, double x41310, double x41311, double x42010, double x42011, double x42110, double x42111, double x42210, double x42211, double x42310, double x42311, double x50010, double x50011, double x50110, double x50111, double x50210, double x50211, double x50310, double x50311, double x51010, double x51011, double x51110, double x51111, double x51210, double x51211, double x51310, double x51311, double x52010, double x52011, double x52110, double x52111, double x52210, double x52211, double x52310, double x52311, double x60010, double x60011, double x60110, double x60111, double x60210, double x60211, double x60310, double x60311, double x61010, double x61011, double x61110, double x61111, double x61210, double x61211, double x61310, double x61311, double x62010, double x62011, double x62110, double x62111, double x62210, double x62211, double x62310, double x62311)
{
	return -(lne*x00110)+upu*(x50110+x50111+x50310+x50311)+up*(x20110+x20111+x40110+x40111+x40310+x40311+x60110+x60111+x60310+x60311)-ac*Nf*x00110*(x20010+x21010+x22010+x40010+x41010+x42010+x50010+x51010+x52010+x60010+x61010+x62010)-ac*Nf*x00110*(x20011+x21011+x22011+x40011+x41011+x42011+x50011+x51011+x52011+x60011+x61011+x62011)+dc*(x20110+x21110+x22110+x40110+x41110+x42110+x50110+x51110+x52110+x60110+x61110+x62110)+dc*(x20111+x21111+x22111+x40111+x41111+x42111+x50111+x51111+x52111+x60111+x61111+x62111)-ac*Nf*x00110*(x40210+x41210+x42210+x50210+x51210+x52210+x60210+x61210+x62210)-ac*Nf*x00110*(x40211+x41211+x42211+x50211+x51211+x52211+x60211+x61211+x62211)+dc*(x40310+x41310+x42310+x50310+x51310+x52310+x60310+x61310+x62310)+dc*(x40311+x41311+x42311+x50311+x51311+x52311+x60311+x61311+x62311);
};
double rhsx00200(double cyrev, double cyrevg, double cyrevgp, double x00200, double x00210, double x30000, double x30100, double x30200, double x30300, double x40000, double x40100, double x40200, double x40300, double x41000, double x41100, double x41200, double x41300, double x42000, double x42100, double x42200, double x42300, double x50000, double x50100, double x50200, double x50300, double x51000, double x51100, double x51200, double x51300, double x52000, double x52100, double x52200, double x52300, double x60000, double x60100, double x60200, double x60300, double x61000, double x61100, double x61200, double x61300, double x62000, double x62100, double x62200, double x62300)
{
	return dg*cyrevg+urev*cyrevg+dg*cyrevgp+uprev*cyrevgp-ag*cyrev*x00200+lne*x00210+upu*(x30200+x30300+x50200+x50300)+up*(x40200+x40300+x60200+x60300)-agp*x00200*(x30000+x30100+x40000+x40100+x41000+x41100+x42000+x42100+x50000+x50100+x51000+x51100+x52000+x52100+x60000+x60100+x61000+x61100+x62000+x62100)+dg*(x30200+x30300+x40200+x40300+x41200+x41300+x42200+x42300+x50200+x50300+x51200+x51300+x52200+x52300+x60200+x60300+x61200+x61300+x62200+x62300);
};
double rhsx00210(double revn, double revng, double revngp, double x00210, double x40010, double x40011, double x40110, double x40111, double x40210, double x40211, double x40310, double x40311, double x41010, double x41011, double x41110, double x41111, double x41210, double x41211, double x41310, double x41311, double x42010, double x42011, double x42110, double x42111, double x42210, double x42211, double x42310, double x42311, double x50010, double x50011, double x50110, double x50111, double x50210, double x50211, double x50310, double x50311, double x51010, double x51011, double x51110, double x51111, double x51210, double x51211, double x51310, double x51311, double x52010, double x52011, double x52110, double x52111, double x52210, double x52211, double x52310, double x52311, double x60010, double x60011, double x60110, double x60111, double x60210, double x60211, double x60310, double x60311, double x61010, double x61011, double x61110, double x61111, double x61210, double x61211, double x61310, double x61311, double x62010, double x62011, double x62110, double x62111, double x62210, double x62211, double x62310, double x62311)
{
	return dg*revng+urev*revng+dg*revngp+uprev*revngp-lne*x00210-ag*Nf*revn*x00210+upu*(x50210+x50211+x50310+x50311)+up*(x40210+x40211+x40310+x40311+x60210+x60211+x60310+x60311)-agp*Nf*x00210*(x40010+x40011+x40110+x40111+x41010+x41011+x41110+x41111+x42010+x42011+x42110+x42111+x50010+x50011+x50110+x50111+x51010+x51011+x51110+x51111+x52010+x52011+x52110+x52111+x60010+x60011+x60110+x60111+x61010+x61011+x61110+x61111+x62010+x62011+x62110+x62111)+dg*(x40210+x40211+x40310+x40311+x41210+x41211+x41310+x41311+x42210+x42211+x42310+x42311+x50210+x50211+x50310+x50311+x51210+x51211+x51310+x51311+x52210+x52211+x52310+x52311+x60210+x60211+x60310+x60311+x61210+x61211+x61310+x61311+x62210+x62211+x62310+x62311);
};
double rhsx01000(double McRo, double x01000, double x20000, double x20100, double x21000, double x21100, double x40000, double x40100, double x40200, double x40300, double x41000, double x41100, double x41200, double x41300, double x50000, double x50100, double x50200, double x50300, double x51000, double x51100, double x51200, double x51300, double x60000, double x60100, double x60200, double x60300, double x61000, double x61100, double x61200, double x61300)
{
	return tlr*McRo-uro*x01000-ar*x01000*(x20000+x20100+x40000+x40100+x40200+x40300+x50000+x50100+x50200+x50300+x60000+x60100+x60200+x60300)+dr*(x21000+x21100+x41000+x41100+x41200+x41300+x51000+x51100+x51200+x51300+x61000+x61100+x61200+x61300);
};
double rhsx01010(double x00011, double x01010, double x01011, double x20010, double x20011, double x20110, double x20111, double x21010, double x21011, double x21110, double x21111, double x40010, double x40011, double x40110, double x40111, double x40210, double x40211, double x40310, double x40311, double x41010, double x41011, double x41110, double x41111, double x41210, double x41211, double x41310, double x41311, double x50010, double x50011, double x50110, double x50111, double x50210, double x50211, double x50310, double x50311, double x51010, double x51011, double x51110, double x51111, double x51210, double x51211, double x51310, double x51311, double x60010, double x60011, double x60110, double x60111, double x60210, double x60211, double x60310, double x60311, double x61010, double x61011, double x61110, double x61111, double x61210, double x61211, double x61310, double x61311)
{
	return -(uro*x01010)-cbbin*Nf*x00011*x01010+uncbbin*x01011-ar*Nf*x01010*(x20010+x20110+x40010+x40110+x40210+x40310+x50010+x50110+x50210+x50310+x60010+x60110+x60210+x60310)-ar*Nf*x01010*(x20011+x20111+x40011+x40111+x40211+x40311+x50011+x50111+x50211+x50311+x60011+x60111+x60211+x60311)+dr*(x21010+x21110+x41010+x41110+x41210+x41310+x51010+x51110+x51210+x51310+x61010+x61110+x61210+x61310)+dr*(x21011+x21111+x41011+x41111+x41211+x41311+x51011+x51111+x51211+x51311+x61011+x61111+x61211+x61311);
};
double rhsx01011(double x00011, double x01010, double x01011, double x20010, double x20110, double x21011, double x21111, double x40010, double x40110, double x40210, double x40310, double x41011, double x41111, double x41211, double x41311, double x50010, double x50110, double x50210, double x50310, double x51011, double x51111, double x51211, double x51311, double x60010, double x60110, double x60210, double x60310, double x61011, double x61111, double x61211, double x61311)
{
	return cbbin*Nf*x00011*x01010-uncbbin*x01011-uro*x01011-ar*Nf*x01011*(x20010+x20110+x40010+x40110+x40210+x40310+x50010+x50110+x50210+x50310+x60010+x60110+x60210+x60310)+dr*(x21011+x21111+x41011+x41111+x41211+x41311+x51011+x51111+x51211+x51311+x61011+x61111+x61211+x61311);
};
double rhsx02000(double McRt, double x02000, double x20000, double x20100, double x22000, double x22100, double x40000, double x40100, double x40200, double x40300, double x42000, double x42100, double x42200, double x42300, double x50000, double x50100, double x50200, double x50300, double x52000, double x52100, double x52200, double x52300, double x60000, double x60100, double x60200, double x60300, double x62000, double x62100, double x62200, double x62300)
{
	return tlr*McRt-urt*x02000-ar*x02000*(x20000+x20100+x40000+x40100+x40200+x40300+x50000+x50100+x50200+x50300+x60000+x60100+x60200+x60300)+dr*(x22000+x22100+x42000+x42100+x42200+x42300+x52000+x52100+x52200+x52300+x62000+x62100+x62200+x62300);
};
double rhsx02010(double x00011, double x02010, double x02011, double x20010, double x20011, double x20110, double x20111, double x22010, double x22011, double x22110, double x22111, double x40010, double x40011, double x40110, double x40111, double x40210, double x40211, double x40310, double x40311, double x42010, double x42011, double x42110, double x42111, double x42210, double x42211, double x42310, double x42311, double x50010, double x50011, double x50110, double x50111, double x50210, double x50211, double x50310, double x50311, double x52010, double x52011, double x52110, double x52111, double x52210, double x52211, double x52310, double x52311, double x60010, double x60011, double x60110, double x60111, double x60210, double x60211, double x60310, double x60311, double x62010, double x62011, double x62110, double x62111, double x62210, double x62211, double x62310, double x62311)
{
	return -(urt*x02010)-cbbin*Nf*x00011*x02010+uncbbin*x02011-ar*Nf*x02010*(x20010+x20110+x40010+x40110+x40210+x40310+x50010+x50110+x50210+x50310+x60010+x60110+x60210+x60310)-ar*Nf*x02010*(x20011+x20111+x40011+x40111+x40211+x40311+x50011+x50111+x50211+x50311+x60011+x60111+x60211+x60311)+dr*(x22010+x22110+x42010+x42110+x42210+x42310+x52010+x52110+x52210+x52310+x62010+x62110+x62210+x62310)+dr*(x22011+x22111+x42011+x42111+x42211+x42311+x52011+x52111+x52211+x52311+x62011+x62111+x62211+x62311);
};
double rhsx02011(double x00011, double x02010, double x02011, double x20010, double x20110, double x22011, double x22111, double x40010, double x40110, double x40210, double x40310, double x42011, double x42111, double x42211, double x42311, double x50010, double x50110, double x50210, double x50310, double x52011, double x52111, double x52211, double x52311, double x60010, double x60110, double x60210, double x60310, double x62011, double x62111, double x62211, double x62311)
{
	return cbbin*Nf*x00011*x02010-uncbbin*x02011-urt*x02011-ar*Nf*x02011*(x20010+x20110+x40010+x40110+x40210+x40310+x50010+x50110+x50210+x50310+x60010+x60110+x60210+x60310)+dr*(x22011+x22111+x42011+x42111+x42211+x42311+x52011+x52111+x52211+x52311+x62011+x62111+x62211+x62311);
};
double rhsx10000(double McPo, double x00100, double x10000, double x10100)
{
	return tlp*McPo-upu*x10000-ac*x00100*x10000+dc*x10100;
};
double rhsx10100(double x00100, double x10000, double x10100)
{
	return ac*x00100*x10000-dc*x10100-hoo*x10100-upu*x10100;
};
double rhsx20000(double x00100, double x01000, double x02000, double x20000, double x20010, double x20100, double x21000, double x22000)
{
	return -(nl*x20000)-up*x20000-ac*x00100*x20000-ar*(x01000+x02000)*x20000+ne*x20010+dc*x20100+dr*(x21000+x22000);
};
double rhsx20010(double x00011, double x00110, double x01010, double x01011, double x02010, double x02011, double x20000, double x20010, double x20011, double x20110, double x21010, double x21011, double x22010, double x22011)
{
	return nl*x20000-ne*x20010-up*x20010-bbin*Nf*x00011*x20010-ac*Nf*x00110*x20010-ar*Nf*(x01010+x02010)*x20010-ar*Nf*(x01011+x02011)*x20010+ubc*x20011+unbbin*x20011+dc*x20110+dr*(x21010+x22010)+dr*(x21011+x22011);
};
double rhsx20011(double x00011, double x00110, double x01010, double x02010, double x20010, double x20011, double x20111, double x21011, double x22011)
{
	return bbin*Nf*x00011*x20010-ubc*x20011-unbbin*x20011-up*x20011-ac*Nf*x00110*x20011-ar*Nf*(x01010+x02010)*x20011+dc*x20111+dr*(x21011+x22011);
};
double rhsx20100(double x00100, double x01000, double x02000, double x10100, double x20000, double x20100, double x20110, double x21100, double x22100)
{
	return hoo*x10100+ac*x00100*x20000-dc*x20100-nl*x20100-up*x20100-ar*(x01000+x02000)*x20100+ne*x20110+dr*(x21100+x22100);
};
double rhsx20110(double x00011, double x00110, double x01010, double x01011, double x02010, double x02011, double x20010, double x20100, double x20110, double x20111, double x21110, double x21111, double x22110, double x22111)
{
	return ac*Nf*x00110*x20010+nl*x20100-dc*x20110-ne*x20110-up*x20110-bbin*Nf*x00011*x20110-ar*Nf*(x01010+x02010)*x20110-ar*Nf*(x01011+x02011)*x20110+ubc*x20111+unbbin*x20111+dr*(x21110+x22110)+dr*(x21111+x22111);
};
double rhsx20111(double x00011, double x00110, double x01010, double x02010, double x20011, double x20110, double x20111, double x21111, double x22111)
{
	return ac*Nf*x00110*x20011+bbin*Nf*x00011*x20110-dc*x20111-ubc*x20111-unbbin*x20111-up*x20111-ar*Nf*(x01010+x02010)*x20111+dr*(x21111+x22111);
};
double rhsx21000(double x00100, double x01000, double x20000, double x21000, double x21010, double x21100)
{
	return ar*x01000*x20000-dr*x21000-nl*x21000-ac*x00100*x21000+ne*x21010+dc*x21100;
};
double rhsx21010(double x00011, double x00110, double x01010, double x20010, double x21000, double x21010, double x21011, double x21110)
{
	return ar*Nf*x01010*x20010+nl*x21000-dr*x21010-ne*x21010-bbin*Nf*x00011*x21010-ac*Nf*x00110*x21010+unbbin*x21011+dc*x21110;
};
double rhsx21011(double x00011, double x00110, double x01010, double x01011, double x20010, double x20011, double x21010, double x21011, double x21111)
{
	return ar*Nf*x01011*x20010+ar*Nf*x01010*x20011+bbin*Nf*x00011*x21010-2*dr*x21011-unbbin*x21011-ac*Nf*x00110*x21011+dc*x21111;
};
double rhsx21100(double x00100, double x01000, double x20100, double x21000, double x21100, double x21110)
{
	return ar*x01000*x20100+ac*x00100*x21000-dc*x21100-dr*x21100-nl*x21100+ne*x21110;
};
double rhsx21110(double x00011, double x00110, double x01010, double x20110, double x21010, double x21100, double x21110, double x21111)
{
	return ar*Nf*x01010*x20110+ac*Nf*x00110*x21010+nl*x21100-dc*x21110-dr*x21110-ne*x21110-bbin*Nf*x00011*x21110+unbbin*x21111;
};
double rhsx21111(double x00011, double x00110, double x01010, double x01011, double x20110, double x20111, double x21011, double x21110, double x21111)
{
	return ar*Nf*x01011*x20110+ar*Nf*x01010*x20111+ac*Nf*x00110*x21011+bbin*Nf*x00011*x21110-dc*x21111-2*dr*x21111-unbbin*x21111;
};
double rhsx22000(double x00100, double x02000, double x20000, double x22000, double x22010, double x22100)
{
	return ar*x02000*x20000-dr*x22000-nl*x22000-ac*x00100*x22000+ne*x22010+dc*x22100;
};
double rhsx22010(double x00011, double x00110, double x02010, double x20010, double x22000, double x22010, double x22011, double x22110)
{
	return ar*Nf*x02010*x20010+nl*x22000-dr*x22010-ne*x22010-bbin*Nf*x00011*x22010-ac*Nf*x00110*x22010+unbbin*x22011+dc*x22110;
};
double rhsx22011(double x00011, double x00110, double x02010, double x02011, double x20010, double x20011, double x22010, double x22011, double x22111)
{
	return ar*Nf*x02011*x20010+ar*Nf*x02010*x20011+bbin*Nf*x00011*x22010-2*dr*x22011-unbbin*x22011-ac*Nf*x00110*x22011+dc*x22111;
};
double rhsx22100(double x00100, double x02000, double x20100, double x22000, double x22100, double x22110)
{
	return ar*x02000*x20100+ac*x00100*x22000-dc*x22100-dr*x22100-nl*x22100+ne*x22110;
};
double rhsx22110(double x00011, double x00110, double x02010, double x20110, double x22010, double x22100, double x22110, double x22111)
{
	return ar*Nf*x02010*x20110+ac*Nf*x00110*x22010+nl*x22100-dc*x22110-dr*x22110-ne*x22110-bbin*Nf*x00011*x22110+unbbin*x22111;
};
double rhsx22111(double x00011, double x00110, double x02010, double x02011, double x20110, double x20111, double x22011, double x22110, double x22111)
{
	return ar*Nf*x02011*x20110+ar*Nf*x02010*x20111+ac*Nf*x00110*x22011+bbin*Nf*x00011*x22110-dc*x22111-2*dr*x22111-unbbin*x22111;
};
double rhsx30000(double McPt, double x00100, double x00200, double x30000, double x30100, double x30200)
{
	return tlp*McPt-upu*x30000-ac*x00100*x30000-agp*x00200*x30000+dc*x30100+dg*x30200;
};
double rhsx30100(double x00100, double x00200, double x30000, double x30100, double x30300)
{
	return ac*x00100*x30000-dc*x30100-hto*x30100-upu*x30100-agp*x00200*x30100+dg*x30300;
};
double rhsx30200(double gto, double x00100, double x00200, double x30000, double x30200, double x30300)
{
	return agp*x00200*x30000-dg*x30200-upu*x30200-gto*x30200-ac*x00100*x30200+dc*x30300;
};
double rhsx30300(double gto, double x00100, double x00200, double x30100, double x30200, double x30300)
{
	return agp*x00200*x30100+ac*x00100*x30200-dc*x30300-dg*x30300-hto*x30300-upu*x30300-gto*x30300;
};
double rhsx40000(double x00100, double x00200, double x01000, double x02000, double x40000, double x40010, double x40100, double x40200, double x41000, double x42000)
{
	return -(nl*x40000)-up*x40000-ac*x00100*x40000-agp*x00200*x40000-ar*(x01000+x02000)*x40000+ne*x40010+dc*x40100+dg*x40200+dr*(x41000+x42000);
};
double rhsx40010(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40000, double x40010, double x40011, double x40110, double x40210, double x41010, double x41011, double x42010, double x42011)
{
	return nl*x40000-ne*x40010-up*x40010-bbin*Nf*x00011*x40010-ac*Nf*x00110*x40010-agp*Nf*x00210*x40010-ar*Nf*(x01010+x02010)*x40010-ar*Nf*(x01011+x02011)*x40010+ubc*x40011+unbbin*x40011+dc*x40110+dg*x40210+dr*(x41010+x42010)+dr*(x41011+x42011);
};
double rhsx40011(double x00011, double x00110, double x00210, double x01010, double x02010, double x40010, double x40011, double x40111, double x40211, double x41011, double x42011)
{
	return bbin*Nf*x00011*x40010-ubc*x40011-unbbin*x40011-up*x40011-ac*Nf*x00110*x40011-agp*Nf*x00210*x40011-ar*Nf*(x01010+x02010)*x40011+dc*x40111+dg*x40211+dr*(x41011+x42011);
};
double rhsx40100(double x00100, double x00200, double x01000, double x02000, double x30100, double x40000, double x40100, double x40110, double x40300, double x41100, double x42100)
{
	return hto*x30100+ac*x00100*x40000-dc*x40100-nl*x40100-up*x40100-agp*x00200*x40100-ar*(x01000+x02000)*x40100+ne*x40110+dg*x40300+dr*(x41100+x42100);
};
double rhsx40110(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40010, double x40100, double x40110, double x40111, double x40310, double x41110, double x41111, double x42110, double x42111)
{
	return ac*Nf*x00110*x40010+nl*x40100-dc*x40110-ne*x40110-up*x40110-bbin*Nf*x00011*x40110-agp*Nf*x00210*x40110-ar*Nf*(x01010+x02010)*x40110-ar*Nf*(x01011+x02011)*x40110+ubc*x40111+unbbin*x40111+dg*x40310+dr*(x41110+x42110)+dr*(x41111+x42111);
};
double rhsx40111(double x00011, double x00110, double x00210, double x01010, double x02010, double x40011, double x40110, double x40111, double x40311, double x41111, double x42111)
{
	return ac*Nf*x00110*x40011+bbin*Nf*x00011*x40110-dc*x40111-ubc*x40111-unbbin*x40111-up*x40111-agp*Nf*x00210*x40111-ar*Nf*(x01010+x02010)*x40111+dg*x40311+dr*(x41111+x42111);
};
double rhsx40200(double gto, double x00100, double x00200, double x01000, double x02000, double x40000, double x40200, double x40210, double x40300, double x41200, double x42200)
{
	return agp*x00200*x40000-dg*x40200-nl*x40200-up*x40200-gto*x40200-ac*x00100*x40200-ar*(x01000+x02000)*x40200+ne*x40210+dc*x40300+dr*(x41200+x42200);
};
double rhsx40210(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40010, double x40200, double x40210, double x40211, double x40310, double x41210, double x41211, double x42210, double x42211)
{
	return agp*Nf*x00210*x40010+nl*x40200-dg*x40210-ne*x40210-up*x40210-gto*x40210-bbin*Nf*x00011*x40210-ac*Nf*x00110*x40210-ar*Nf*(x01010+x02010)*x40210-ar*Nf*(x01011+x02011)*x40210+ubc*x40211+unbbin*x40211+dc*x40310+dr*(x41210+x42210)+dr*(x41211+x42211);
};
double rhsx40211(double gto, double x00011, double x00110, double x00210, double x01010, double x02010, double x40011, double x40210, double x40211, double x40311, double x41211, double x42211)
{
	return agp*Nf*x00210*x40011+bbin*Nf*x00011*x40210-dg*x40211-ubc*x40211-unbbin*x40211-up*x40211-gto*x40211-ac*Nf*x00110*x40211-ar*Nf*(x01010+x02010)*x40211+dc*x40311+dr*(x41211+x42211);
};
double rhsx40300(double gto, double x00100, double x00200, double x01000, double x02000, double x30300, double x40100, double x40200, double x40300, double x40310, double x41300, double x42300)
{
	return hto*x30300+agp*x00200*x40100+ac*x00100*x40200-dc*x40300-dg*x40300-nl*x40300-up*x40300-gto*x40300-ar*(x01000+x02000)*x40300+ne*x40310+dr*(x41300+x42300);
};
double rhsx40310(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40110, double x40210, double x40300, double x40310, double x40311, double x41310, double x41311, double x42310, double x42311)
{
	return agp*Nf*x00210*x40110+ac*Nf*x00110*x40210+nl*x40300-dc*x40310-dg*x40310-ne*x40310-up*x40310-gto*x40310-bbin*Nf*x00011*x40310-ar*Nf*(x01010+x02010)*x40310-ar*Nf*(x01011+x02011)*x40310+ubc*x40311+unbbin*x40311+dr*(x41310+x42310)+dr*(x41311+x42311);
};
double rhsx40311(double gto, double x00011, double x00110, double x00210, double x01010, double x02010, double x40111, double x40211, double x40310, double x40311, double x41311, double x42311)
{
	return agp*Nf*x00210*x40111+ac*Nf*x00110*x40211+bbin*Nf*x00011*x40310-dc*x40311-dg*x40311-ubc*x40311-unbbin*x40311-up*x40311-gto*x40311-ar*Nf*(x01010+x02010)*x40311+dr*(x41311+x42311);
};
double rhsx41000(double x00100, double x00200, double x01000, double x40000, double x41000, double x41010, double x41100, double x41200)
{
	return ar*x01000*x40000-dr*x41000-nl*x41000-ac*x00100*x41000-agp*x00200*x41000+ne*x41010+dc*x41100+dg*x41200;
};
double rhsx41010(double x00011, double x00110, double x00210, double x01010, double x40010, double x41000, double x41010, double x41011, double x41110, double x41210)
{
	return ar*Nf*x01010*x40010+nl*x41000-dr*x41010-ne*x41010-bbin*Nf*x00011*x41010-ac*Nf*x00110*x41010-agp*Nf*x00210*x41010+unbbin*x41011+dc*x41110+dg*x41210;
};
double rhsx41011(double x00011, double x00110, double x00210, double x01010, double x01011, double x40010, double x40011, double x41010, double x41011, double x41111, double x41211)
{
	return ar*Nf*x01011*x40010+ar*Nf*x01010*x40011+bbin*Nf*x00011*x41010-2*dr*x41011-unbbin*x41011-ac*Nf*x00110*x41011-agp*Nf*x00210*x41011+dc*x41111+dg*x41211;
};
double rhsx41100(double x00100, double x00200, double x01000, double x40100, double x41000, double x41100, double x41110, double x41300)
{
	return ar*x01000*x40100+ac*x00100*x41000-dc*x41100-dr*x41100-nl*x41100-agp*x00200*x41100+ne*x41110+dg*x41300;
};
double rhsx41110(double x00011, double x00110, double x00210, double x01010, double x40110, double x41010, double x41100, double x41110, double x41111, double x41310)
{
	return ar*Nf*x01010*x40110+ac*Nf*x00110*x41010+nl*x41100-dc*x41110-dr*x41110-ne*x41110-bbin*Nf*x00011*x41110-agp*Nf*x00210*x41110+unbbin*x41111+dg*x41310;
};
double rhsx41111(double x00011, double x00110, double x00210, double x01010, double x01011, double x40110, double x40111, double x41011, double x41110, double x41111, double x41311)
{
	return ar*Nf*x01011*x40110+ar*Nf*x01010*x40111+ac*Nf*x00110*x41011+bbin*Nf*x00011*x41110-dc*x41111-2*dr*x41111-unbbin*x41111-agp*Nf*x00210*x41111+dg*x41311;
};
double rhsx41200(double gto, double x00100, double x00200, double x01000, double x40200, double x41000, double x41200, double x41210, double x41300)
{
	return ar*x01000*x40200+agp*x00200*x41000-dg*x41200-dr*x41200-nl*x41200-gto*x41200-ac*x00100*x41200+ne*x41210+dc*x41300;
};
double rhsx41210(double gto, double x00011, double x00110, double x00210, double x01010, double x40210, double x41010, double x41200, double x41210, double x41211, double x41310)
{
	return ar*Nf*x01010*x40210+agp*Nf*x00210*x41010+nl*x41200-dg*x41210-dr*x41210-ne*x41210-gto*x41210-bbin*Nf*x00011*x41210-ac*Nf*x00110*x41210+unbbin*x41211+dc*x41310;
};
double rhsx41211(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x40210, double x40211, double x41011, double x41210, double x41211, double x41311)
{
	return ar*Nf*x01011*x40210+ar*Nf*x01010*x40211+agp*Nf*x00210*x41011+bbin*Nf*x00011*x41210-dg*x41211-2*dr*x41211-unbbin*x41211-gto*x41211-ac*Nf*x00110*x41211+dc*x41311;
};
double rhsx41300(double gto, double x00100, double x00200, double x01000, double x40300, double x41100, double x41200, double x41300, double x41310)
{
	return ar*x01000*x40300+agp*x00200*x41100+ac*x00100*x41200-dc*x41300-dg*x41300-dr*x41300-nl*x41300-gto*x41300+ne*x41310;
};
double rhsx41310(double gto, double x00011, double x00110, double x00210, double x01010, double x40310, double x41110, double x41210, double x41300, double x41310, double x41311)
{
	return ar*Nf*x01010*x40310+agp*Nf*x00210*x41110+ac*Nf*x00110*x41210+nl*x41300-dc*x41310-dg*x41310-dr*x41310-ne*x41310-gto*x41310-bbin*Nf*x00011*x41310+unbbin*x41311;
};
double rhsx41311(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x40310, double x40311, double x41111, double x41211, double x41310, double x41311)
{
	return ar*Nf*x01011*x40310+ar*Nf*x01010*x40311+agp*Nf*x00210*x41111+ac*Nf*x00110*x41211+bbin*Nf*x00011*x41310-dc*x41311-dg*x41311-2*dr*x41311-unbbin*x41311-gto*x41311;
};
double rhsx42000(double x00100, double x00200, double x02000, double x40000, double x42000, double x42010, double x42100, double x42200)
{
	return ar*x02000*x40000-dr*x42000-nl*x42000-ac*x00100*x42000-agp*x00200*x42000+ne*x42010+dc*x42100+dg*x42200;
};
double rhsx42010(double x00011, double x00110, double x00210, double x02010, double x40010, double x42000, double x42010, double x42011, double x42110, double x42210)
{
	return ar*Nf*x02010*x40010+nl*x42000-dr*x42010-ne*x42010-bbin*Nf*x00011*x42010-ac*Nf*x00110*x42010-agp*Nf*x00210*x42010+unbbin*x42011+dc*x42110+dg*x42210;
};
double rhsx42011(double x00011, double x00110, double x00210, double x02010, double x02011, double x40010, double x40011, double x42010, double x42011, double x42111, double x42211)
{
	return ar*Nf*x02011*x40010+ar*Nf*x02010*x40011+bbin*Nf*x00011*x42010-2*dr*x42011-unbbin*x42011-ac*Nf*x00110*x42011-agp*Nf*x00210*x42011+dc*x42111+dg*x42211;
};
double rhsx42100(double x00100, double x00200, double x02000, double x40100, double x42000, double x42100, double x42110, double x42300)
{
	return ar*x02000*x40100+ac*x00100*x42000-dc*x42100-dr*x42100-nl*x42100-agp*x00200*x42100+ne*x42110+dg*x42300;
};
double rhsx42110(double x00011, double x00110, double x00210, double x02010, double x40110, double x42010, double x42100, double x42110, double x42111, double x42310)
{
	return ar*Nf*x02010*x40110+ac*Nf*x00110*x42010+nl*x42100-dc*x42110-dr*x42110-ne*x42110-bbin*Nf*x00011*x42110-agp*Nf*x00210*x42110+unbbin*x42111+dg*x42310;
};
double rhsx42111(double x00011, double x00110, double x00210, double x02010, double x02011, double x40110, double x40111, double x42011, double x42110, double x42111, double x42311)
{
	return ar*Nf*x02011*x40110+ar*Nf*x02010*x40111+ac*Nf*x00110*x42011+bbin*Nf*x00011*x42110-dc*x42111-2*dr*x42111-unbbin*x42111-agp*Nf*x00210*x42111+dg*x42311;
};
double rhsx42200(double gto, double x00100, double x00200, double x02000, double x40200, double x42000, double x42200, double x42210, double x42300)
{
	return ar*x02000*x40200+agp*x00200*x42000-dg*x42200-dr*x42200-nl*x42200-gto*x42200-ac*x00100*x42200+ne*x42210+dc*x42300;
};
double rhsx42210(double gto, double x00011, double x00110, double x00210, double x02010, double x40210, double x42010, double x42200, double x42210, double x42211, double x42310)
{
	return ar*Nf*x02010*x40210+agp*Nf*x00210*x42010+nl*x42200-dg*x42210-dr*x42210-ne*x42210-gto*x42210-bbin*Nf*x00011*x42210-ac*Nf*x00110*x42210+unbbin*x42211+dc*x42310;
};
double rhsx42211(double gto, double x00011, double x00110, double x00210, double x02010, double x02011, double x40210, double x40211, double x42011, double x42210, double x42211, double x42311)
{
	return ar*Nf*x02011*x40210+ar*Nf*x02010*x40211+agp*Nf*x00210*x42011+bbin*Nf*x00011*x42210-dg*x42211-2*dr*x42211-unbbin*x42211-gto*x42211-ac*Nf*x00110*x42211+dc*x42311;
};
double rhsx42300(double gto, double x00100, double x00200, double x02000, double x40300, double x42100, double x42200, double x42300, double x42310)
{
	return ar*x02000*x40300+agp*x00200*x42100+ac*x00100*x42200-dc*x42300-dg*x42300-dr*x42300-nl*x42300-gto*x42300+ne*x42310;
};
double rhsx42310(double gto, double x00011, double x00110, double x00210, double x02010, double x40310, double x42110, double x42210, double x42300, double x42310, double x42311)
{
	return ar*Nf*x02010*x40310+agp*Nf*x00210*x42110+ac*Nf*x00110*x42210+nl*x42300-dc*x42310-dg*x42310-dr*x42310-ne*x42310-gto*x42310-bbin*Nf*x00011*x42310+unbbin*x42311;
};
double rhsx42311(double gto, double x00011, double x00110, double x00210, double x02010, double x02011, double x40310, double x40311, double x42111, double x42211, double x42310, double x42311)
{
	return ar*Nf*x02011*x40310+ar*Nf*x02010*x40311+agp*Nf*x00210*x42111+ac*Nf*x00110*x42211+bbin*Nf*x00011*x42310-dc*x42311-dg*x42311-2*dr*x42311-unbbin*x42311-gto*x42311;
};
double rhsx50000(double x00100, double x00200, double x01000, double x02000, double x50000, double x50010, double x50100, double x50200, double x51000, double x52000)
{
	return -(nl*x50000)-upu*x50000-ac*x00100*x50000-agp*x00200*x50000-ar*(x01000+x02000)*x50000+ne*x50010+dc*x50100+dg*x50200+dr*(x51000+x52000);
};
double rhsx50010(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x50000, double x50010, double x50011, double x50110, double x50210, double x51010, double x51011, double x52010, double x52011)
{
	return nl*x50000-ne*x50010-upu*x50010-bbin*Nf*x00011*x50010-ac*Nf*x00110*x50010-agp*Nf*x00210*x50010-ar*Nf*(x01010+x02010)*x50010-ar*Nf*(x01011+x02011)*x50010+ubc*x50011+unbbin*x50011+dc*x50110+dg*x50210+dr*(x51010+x52010)+dr*(x51011+x52011);
};
double rhsx50011(double x00011, double x00110, double x00210, double x01010, double x02010, double x50010, double x50011, double x50111, double x50211, double x51011, double x52011)
{
	return bbin*Nf*x00011*x50010-ubc*x50011-unbbin*x50011-upu*x50011-ac*Nf*x00110*x50011-agp*Nf*x00210*x50011-ar*Nf*(x01010+x02010)*x50011+dc*x50111+dg*x50211+dr*(x51011+x52011);
};
double rhsx50100(double x00100, double x00200, double x01000, double x02000, double x50000, double x50100, double x50110, double x50300, double x51100, double x52100)
{
	return ac*x00100*x50000-dc*x50100-hto*x50100-nl*x50100-upu*x50100-agp*x00200*x50100-ar*(x01000+x02000)*x50100+ne*x50110+dg*x50300+dr*(x51100+x52100);
};
double rhsx50110(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x50010, double x50100, double x50110, double x50111, double x50310, double x51110, double x51111, double x52110, double x52111)
{
	return ac*Nf*x00110*x50010+nl*x50100-dc*x50110-hto*x50110-ne*x50110-upu*x50110-bbin*Nf*x00011*x50110-agp*Nf*x00210*x50110-ar*Nf*(x01010+x02010)*x50110-ar*Nf*(x01011+x02011)*x50110+ubc*x50111+unbbin*x50111+dg*x50310+dr*(x51110+x52110)+dr*(x51111+x52111);
};
double rhsx50111(double x00011, double x00110, double x00210, double x01010, double x02010, double x50011, double x50110, double x50111, double x50311, double x51111, double x52111)
{
	return ac*Nf*x00110*x50011+bbin*Nf*x00011*x50110-dc*x50111-hto*x50111-ubc*x50111-unbbin*x50111-upu*x50111-agp*Nf*x00210*x50111-ar*Nf*(x01010+x02010)*x50111+dg*x50311+dr*(x51111+x52111);
};
double rhsx50200(double gto, double x00100, double x00200, double x01000, double x02000, double x30200, double x50000, double x50200, double x50210, double x50300, double x51200, double x52200)
{
	return gto*x30200+agp*x00200*x50000-dg*x50200-nl*x50200-upu*x50200-ac*x00100*x50200-ar*(x01000+x02000)*x50200+ne*x50210+dc*x50300+dr*(x51200+x52200);
};
double rhsx50210(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x50010, double x50200, double x50210, double x50211, double x50310, double x51210, double x51211, double x52210, double x52211)
{
	return agp*Nf*x00210*x50010+nl*x50200-dg*x50210-ne*x50210-upu*x50210-bbin*Nf*x00011*x50210-ac*Nf*x00110*x50210-ar*Nf*(x01010+x02010)*x50210-ar*Nf*(x01011+x02011)*x50210+ubc*x50211+unbbin*x50211+dc*x50310+dr*(x51210+x52210)+dr*(x51211+x52211);
};
double rhsx50211(double x00011, double x00110, double x00210, double x01010, double x02010, double x50011, double x50210, double x50211, double x50311, double x51211, double x52211)
{
	return agp*Nf*x00210*x50011+bbin*Nf*x00011*x50210-dg*x50211-ubc*x50211-unbbin*x50211-upu*x50211-ac*Nf*x00110*x50211-ar*Nf*(x01010+x02010)*x50211+dc*x50311+dr*(x51211+x52211);
};
double rhsx50300(double gto, double x00100, double x00200, double x01000, double x02000, double x30300, double x50100, double x50200, double x50300, double x50310, double x51300, double x52300)
{
	return gto*x30300+agp*x00200*x50100+ac*x00100*x50200-dc*x50300-dg*x50300-hto*x50300-nl*x50300-upu*x50300-ar*(x01000+x02000)*x50300+ne*x50310+dr*(x51300+x52300);
};
double rhsx50310(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x50110, double x50210, double x50300, double x50310, double x50311, double x51310, double x51311, double x52310, double x52311)
{
	return agp*Nf*x00210*x50110+ac*Nf*x00110*x50210+nl*x50300-dc*x50310-dg*x50310-hto*x50310-ne*x50310-upu*x50310-bbin*Nf*x00011*x50310-ar*Nf*(x01010+x02010)*x50310-ar*Nf*(x01011+x02011)*x50310+ubc*x50311+unbbin*x50311+dr*(x51310+x52310)+dr*(x51311+x52311);
};
double rhsx50311(double x00011, double x00110, double x00210, double x01010, double x02010, double x50111, double x50211, double x50310, double x50311, double x51311, double x52311)
{
	return agp*Nf*x00210*x50111+ac*Nf*x00110*x50211+bbin*Nf*x00011*x50310-dc*x50311-dg*x50311-hto*x50311-ubc*x50311-unbbin*x50311-upu*x50311-ar*Nf*(x01010+x02010)*x50311+dr*(x51311+x52311);
};
double rhsx51000(double x00100, double x00200, double x01000, double x50000, double x51000, double x51010, double x51100, double x51200)
{
	return ar*x01000*x50000-dr*x51000-nl*x51000-ac*x00100*x51000-agp*x00200*x51000+ne*x51010+dc*x51100+dg*x51200;
};
double rhsx51010(double x00011, double x00110, double x00210, double x01010, double x50010, double x51000, double x51010, double x51011, double x51110, double x51210)
{
	return ar*Nf*x01010*x50010+nl*x51000-dr*x51010-ne*x51010-bbin*Nf*x00011*x51010-ac*Nf*x00110*x51010-agp*Nf*x00210*x51010+unbbin*x51011+dc*x51110+dg*x51210;
};
double rhsx51011(double x00011, double x00110, double x00210, double x01010, double x01011, double x50010, double x50011, double x51010, double x51011, double x51111, double x51211)
{
	return ar*Nf*x01011*x50010+ar*Nf*x01010*x50011+bbin*Nf*x00011*x51010-2*dr*x51011-unbbin*x51011-ac*Nf*x00110*x51011-agp*Nf*x00210*x51011+dc*x51111+dg*x51211;
};
double rhsx51100(double x00100, double x00200, double x01000, double x50100, double x51000, double x51100, double x51110, double x51300)
{
	return ar*x01000*x50100+ac*x00100*x51000-dc*x51100-dr*x51100-nl*x51100-agp*x00200*x51100+ne*x51110+dg*x51300;
};
double rhsx51110(double x00011, double x00110, double x00210, double x01010, double x50110, double x51010, double x51100, double x51110, double x51111, double x51310)
{
	return ar*Nf*x01010*x50110+ac*Nf*x00110*x51010+nl*x51100-dc*x51110-dr*x51110-ne*x51110-bbin*Nf*x00011*x51110-agp*Nf*x00210*x51110+unbbin*x51111+dg*x51310;
};
double rhsx51111(double x00011, double x00110, double x00210, double x01010, double x01011, double x50110, double x50111, double x51011, double x51110, double x51111, double x51311)
{
	return ar*Nf*x01011*x50110+ar*Nf*x01010*x50111+ac*Nf*x00110*x51011+bbin*Nf*x00011*x51110-dc*x51111-2*dr*x51111-unbbin*x51111-agp*Nf*x00210*x51111+dg*x51311;
};
double rhsx51200(double x00100, double x00200, double x01000, double x50200, double x51000, double x51200, double x51210, double x51300)
{
	return ar*x01000*x50200+agp*x00200*x51000-dg*x51200-dr*x51200-nl*x51200-ac*x00100*x51200+ne*x51210+dc*x51300;
};
double rhsx51210(double x00011, double x00110, double x00210, double x01010, double x50210, double x51010, double x51200, double x51210, double x51211, double x51310)
{
	return ar*Nf*x01010*x50210+agp*Nf*x00210*x51010+nl*x51200-dg*x51210-dr*x51210-ne*x51210-bbin*Nf*x00011*x51210-ac*Nf*x00110*x51210+unbbin*x51211+dc*x51310;
};
double rhsx51211(double x00011, double x00110, double x00210, double x01010, double x01011, double x50210, double x50211, double x51011, double x51210, double x51211, double x51311)
{
	return ar*Nf*x01011*x50210+ar*Nf*x01010*x50211+agp*Nf*x00210*x51011+bbin*Nf*x00011*x51210-dg*x51211-2*dr*x51211-unbbin*x51211-ac*Nf*x00110*x51211+dc*x51311;
};
double rhsx51300(double x00100, double x00200, double x01000, double x50300, double x51100, double x51200, double x51300, double x51310)
{
	return ar*x01000*x50300+agp*x00200*x51100+ac*x00100*x51200-dc*x51300-dg*x51300-dr*x51300-nl*x51300+ne*x51310;
};
double rhsx51310(double x00011, double x00110, double x00210, double x01010, double x50310, double x51110, double x51210, double x51300, double x51310, double x51311)
{
	return ar*Nf*x01010*x50310+agp*Nf*x00210*x51110+ac*Nf*x00110*x51210+nl*x51300-dc*x51310-dg*x51310-dr*x51310-ne*x51310-bbin*Nf*x00011*x51310+unbbin*x51311;
};
double rhsx51311(double x00011, double x00110, double x00210, double x01010, double x01011, double x50310, double x50311, double x51111, double x51211, double x51310, double x51311)
{
	return ar*Nf*x01011*x50310+ar*Nf*x01010*x50311+agp*Nf*x00210*x51111+ac*Nf*x00110*x51211+bbin*Nf*x00011*x51310-dc*x51311-dg*x51311-2*dr*x51311-unbbin*x51311;
};
double rhsx52000(double x00100, double x00200, double x02000, double x50000, double x52000, double x52010, double x52100, double x52200)
{
	return ar*x02000*x50000-dr*x52000-nl*x52000-ac*x00100*x52000-agp*x00200*x52000+ne*x52010+dc*x52100+dg*x52200;
};
double rhsx52010(double x00011, double x00110, double x00210, double x02010, double x50010, double x52000, double x52010, double x52011, double x52110, double x52210)
{
	return ar*Nf*x02010*x50010+nl*x52000-dr*x52010-ne*x52010-bbin*Nf*x00011*x52010-ac*Nf*x00110*x52010-agp*Nf*x00210*x52010+unbbin*x52011+dc*x52110+dg*x52210;
};
double rhsx52011(double x00011, double x00110, double x00210, double x02010, double x02011, double x50010, double x50011, double x52010, double x52011, double x52111, double x52211)
{
	return ar*Nf*x02011*x50010+ar*Nf*x02010*x50011+bbin*Nf*x00011*x52010-2*dr*x52011-unbbin*x52011-ac*Nf*x00110*x52011-agp*Nf*x00210*x52011+dc*x52111+dg*x52211;
};
double rhsx52100(double x00100, double x00200, double x02000, double x50100, double x52000, double x52100, double x52110, double x52300)
{
	return ar*x02000*x50100+ac*x00100*x52000-dc*x52100-dr*x52100-nl*x52100-agp*x00200*x52100+ne*x52110+dg*x52300;
};
double rhsx52110(double x00011, double x00110, double x00210, double x02010, double x50110, double x52010, double x52100, double x52110, double x52111, double x52310)
{
	return ar*Nf*x02010*x50110+ac*Nf*x00110*x52010+nl*x52100-dc*x52110-dr*x52110-ne*x52110-bbin*Nf*x00011*x52110-agp*Nf*x00210*x52110+unbbin*x52111+dg*x52310;
};
double rhsx52111(double x00011, double x00110, double x00210, double x02010, double x02011, double x50110, double x50111, double x52011, double x52110, double x52111, double x52311)
{
	return ar*Nf*x02011*x50110+ar*Nf*x02010*x50111+ac*Nf*x00110*x52011+bbin*Nf*x00011*x52110-dc*x52111-2*dr*x52111-unbbin*x52111-agp*Nf*x00210*x52111+dg*x52311;
};
double rhsx52200(double x00100, double x00200, double x02000, double x50200, double x52000, double x52200, double x52210, double x52300)
{
	return ar*x02000*x50200+agp*x00200*x52000-dg*x52200-dr*x52200-nl*x52200-ac*x00100*x52200+ne*x52210+dc*x52300;
};
double rhsx52210(double x00011, double x00110, double x00210, double x02010, double x50210, double x52010, double x52200, double x52210, double x52211, double x52310)
{
	return ar*Nf*x02010*x50210+agp*Nf*x00210*x52010+nl*x52200-dg*x52210-dr*x52210-ne*x52210-bbin*Nf*x00011*x52210-ac*Nf*x00110*x52210+unbbin*x52211+dc*x52310;
};
double rhsx52211(double x00011, double x00110, double x00210, double x02010, double x02011, double x50210, double x50211, double x52011, double x52210, double x52211, double x52311)
{
	return ar*Nf*x02011*x50210+ar*Nf*x02010*x50211+agp*Nf*x00210*x52011+bbin*Nf*x00011*x52210-dg*x52211-2*dr*x52211-unbbin*x52211-ac*Nf*x00110*x52211+dc*x52311;
};
double rhsx52300(double x00100, double x00200, double x02000, double x50300, double x52100, double x52200, double x52300, double x52310)
{
	return ar*x02000*x50300+agp*x00200*x52100+ac*x00100*x52200-dc*x52300-dg*x52300-dr*x52300-nl*x52300+ne*x52310;
};
double rhsx52310(double x00011, double x00110, double x00210, double x02010, double x50310, double x52110, double x52210, double x52300, double x52310, double x52311)
{
	return ar*Nf*x02010*x50310+agp*Nf*x00210*x52110+ac*Nf*x00110*x52210+nl*x52300-dc*x52310-dg*x52310-dr*x52310-ne*x52310-bbin*Nf*x00011*x52310+unbbin*x52311;
};
double rhsx52311(double x00011, double x00110, double x00210, double x02010, double x02011, double x50310, double x50311, double x52111, double x52211, double x52310, double x52311)
{
	return ar*Nf*x02011*x50310+ar*Nf*x02010*x50311+agp*Nf*x00210*x52111+ac*Nf*x00110*x52211+bbin*Nf*x00011*x52310-dc*x52311-dg*x52311-2*dr*x52311-unbbin*x52311;
};
double rhsx60000(double x00100, double x00200, double x01000, double x02000, double x60000, double x60010, double x60100, double x60200, double x61000, double x62000)
{
	return -(nl*x60000)-up*x60000-ac*x00100*x60000-agp*x00200*x60000-ar*(x01000+x02000)*x60000+ne*x60010+dc*x60100+dg*x60200+dr*(x61000+x62000);
};
double rhsx60010(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x60000, double x60010, double x60011, double x60110, double x60210, double x61010, double x61011, double x62010, double x62011)
{
	return nl*x60000-ne*x60010-up*x60010-bbin*Nf*x00011*x60010-ac*Nf*x00110*x60010-agp*Nf*x00210*x60010-ar*Nf*(x01010+x02010)*x60010-ar*Nf*(x01011+x02011)*x60010+ubc*x60011+unbbin*x60011+dc*x60110+dg*x60210+dr*(x61010+x62010)+dr*(x61011+x62011);
};
double rhsx60011(double x00011, double x00110, double x00210, double x01010, double x02010, double x60010, double x60011, double x60111, double x60211, double x61011, double x62011)
{
	return bbin*Nf*x00011*x60010-ubc*x60011-unbbin*x60011-up*x60011-ac*Nf*x00110*x60011-agp*Nf*x00210*x60011-ar*Nf*(x01010+x02010)*x60011+dc*x60111+dg*x60211+dr*(x61011+x62011);
};
double rhsx60100(double x00100, double x00200, double x01000, double x02000, double x50100, double x60000, double x60100, double x60110, double x60300, double x61100, double x62100)
{
	return hto*x50100+ac*x00100*x60000-dc*x60100-nl*x60100-up*x60100-agp*x00200*x60100-ar*(x01000+x02000)*x60100+ne*x60110+dg*x60300+dr*(x61100+x62100);
};
double rhsx60110(double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x50110, double x60010, double x60100, double x60110, double x60111, double x60310, double x61110, double x61111, double x62110, double x62111)
{
	return hto*x50110+ac*Nf*x00110*x60010+nl*x60100-dc*x60110-ne*x60110-up*x60110-bbin*Nf*x00011*x60110-agp*Nf*x00210*x60110-ar*Nf*(x01010+x02010)*x60110-ar*Nf*(x01011+x02011)*x60110+ubc*x60111+unbbin*x60111+dg*x60310+dr*(x61110+x62110)+dr*(x61111+x62111);
};
double rhsx60111(double x00011, double x00110, double x00210, double x01010, double x02010, double x50111, double x60011, double x60110, double x60111, double x60311, double x61111, double x62111)
{
	return hto*x50111+ac*Nf*x00110*x60011+bbin*Nf*x00011*x60110-dc*x60111-ubc*x60111-unbbin*x60111-up*x60111-agp*Nf*x00210*x60111-ar*Nf*(x01010+x02010)*x60111+dg*x60311+dr*(x61111+x62111);
};
double rhsx60200(double gto, double x00100, double x00200, double x01000, double x02000, double x40200, double x60000, double x60200, double x60210, double x60300, double x61200, double x62200)
{
	return gto*x40200+agp*x00200*x60000-dg*x60200-nl*x60200-up*x60200-ac*x00100*x60200-ar*(x01000+x02000)*x60200+ne*x60210+dc*x60300+dr*(x61200+x62200);
};
double rhsx60210(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40210, double x60010, double x60200, double x60210, double x60211, double x60310, double x61210, double x61211, double x62210, double x62211)
{
	return gto*x40210+agp*Nf*x00210*x60010+nl*x60200-dg*x60210-ne*x60210-up*x60210-bbin*Nf*x00011*x60210-ac*Nf*x00110*x60210-ar*Nf*(x01010+x02010)*x60210-ar*Nf*(x01011+x02011)*x60210+ubc*x60211+unbbin*x60211+dc*x60310+dr*(x61210+x62210)+dr*(x61211+x62211);
};
double rhsx60211(double gto, double x00011, double x00110, double x00210, double x01010, double x02010, double x40211, double x60011, double x60210, double x60211, double x60311, double x61211, double x62211)
{
	return gto*x40211+agp*Nf*x00210*x60011+bbin*Nf*x00011*x60210-dg*x60211-ubc*x60211-unbbin*x60211-up*x60211-ac*Nf*x00110*x60211-ar*Nf*(x01010+x02010)*x60211+dc*x60311+dr*(x61211+x62211);
};
double rhsx60300(double gto, double x00100, double x00200, double x01000, double x02000, double x40300, double x50300, double x60100, double x60200, double x60300, double x60310, double x61300, double x62300)
{
	return gto*x40300+hto*x50300+agp*x00200*x60100+ac*x00100*x60200-dc*x60300-dg*x60300-nl*x60300-up*x60300-ar*(x01000+x02000)*x60300+ne*x60310+dr*(x61300+x62300);
};
double rhsx60310(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x02010, double x02011, double x40310, double x50310, double x60110, double x60210, double x60300, double x60310, double x60311, double x61310, double x61311, double x62310, double x62311)
{
	return gto*x40310+hto*x50310+agp*Nf*x00210*x60110+ac*Nf*x00110*x60210+nl*x60300-dc*x60310-dg*x60310-ne*x60310-up*x60310-bbin*Nf*x00011*x60310-ar*Nf*(x01010+x02010)*x60310-ar*Nf*(x01011+x02011)*x60310+ubc*x60311+unbbin*x60311+dr*(x61310+x62310)+dr*(x61311+x62311);
};
double rhsx60311(double gto, double x00011, double x00110, double x00210, double x01010, double x02010, double x40311, double x50311, double x60111, double x60211, double x60310, double x60311, double x61311, double x62311)
{
	return gto*x40311+hto*x50311+agp*Nf*x00210*x60111+ac*Nf*x00110*x60211+bbin*Nf*x00011*x60310-dc*x60311-dg*x60311-ubc*x60311-unbbin*x60311-up*x60311-ar*Nf*(x01010+x02010)*x60311+dr*(x61311+x62311);
};
double rhsx61000(double x00100, double x00200, double x01000, double x60000, double x61000, double x61010, double x61100, double x61200)
{
	return ar*x01000*x60000-dr*x61000-nl*x61000-ac*x00100*x61000-agp*x00200*x61000+ne*x61010+dc*x61100+dg*x61200;
};
double rhsx61010(double x00011, double x00110, double x00210, double x01010, double x60010, double x61000, double x61010, double x61011, double x61110, double x61210)
{
	return ar*Nf*x01010*x60010+nl*x61000-dr*x61010-ne*x61010-bbin*Nf*x00011*x61010-ac*Nf*x00110*x61010-agp*Nf*x00210*x61010+unbbin*x61011+dc*x61110+dg*x61210;
};
double rhsx61011(double x00011, double x00110, double x00210, double x01010, double x01011, double x60010, double x60011, double x61010, double x61011, double x61111, double x61211)
{
	return ar*Nf*x01011*x60010+ar*Nf*x01010*x60011+bbin*Nf*x00011*x61010-2*dr*x61011-unbbin*x61011-ac*Nf*x00110*x61011-agp*Nf*x00210*x61011+dc*x61111+dg*x61211;
};
double rhsx61100(double x00100, double x00200, double x01000, double x60100, double x61000, double x61100, double x61110, double x61300)
{
	return ar*x01000*x60100+ac*x00100*x61000-dc*x61100-dr*x61100-nl*x61100-agp*x00200*x61100+ne*x61110+dg*x61300;
};
double rhsx61110(double x00011, double x00110, double x00210, double x01010, double x60110, double x61010, double x61100, double x61110, double x61111, double x61310)
{
	return ar*Nf*x01010*x60110+ac*Nf*x00110*x61010+nl*x61100-dc*x61110-dr*x61110-ne*x61110-bbin*Nf*x00011*x61110-agp*Nf*x00210*x61110+unbbin*x61111+dg*x61310;
};
double rhsx61111(double x00011, double x00110, double x00210, double x01010, double x01011, double x60110, double x60111, double x61011, double x61110, double x61111, double x61311)
{
	return ar*Nf*x01011*x60110+ar*Nf*x01010*x60111+ac*Nf*x00110*x61011+bbin*Nf*x00011*x61110-dc*x61111-2*dr*x61111-unbbin*x61111-agp*Nf*x00210*x61111+dg*x61311;
};
double rhsx61200(double gto, double x00100, double x00200, double x01000, double x41200, double x60200, double x61000, double x61200, double x61210, double x61300)
{
	return gto*x41200+ar*x01000*x60200+agp*x00200*x61000-dg*x61200-dr*x61200-nl*x61200-ac*x00100*x61200+ne*x61210+dc*x61300;
};
double rhsx61210(double gto, double x00011, double x00110, double x00210, double x01010, double x41210, double x60210, double x61010, double x61200, double x61210, double x61211, double x61310)
{
	return gto*x41210+ar*Nf*x01010*x60210+agp*Nf*x00210*x61010+nl*x61200-dg*x61210-dr*x61210-ne*x61210-bbin*Nf*x00011*x61210-ac*Nf*x00110*x61210+unbbin*x61211+dc*x61310;
};
double rhsx61211(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x41211, double x60210, double x60211, double x61011, double x61210, double x61211, double x61311)
{
	return gto*x41211+ar*Nf*x01011*x60210+ar*Nf*x01010*x60211+agp*Nf*x00210*x61011+bbin*Nf*x00011*x61210-dg*x61211-2*dr*x61211-unbbin*x61211-ac*Nf*x00110*x61211+dc*x61311;
};
double rhsx61300(double gto, double x00100, double x00200, double x01000, double x41300, double x60300, double x61100, double x61200, double x61300, double x61310)
{
	return gto*x41300+ar*x01000*x60300+agp*x00200*x61100+ac*x00100*x61200-dc*x61300-dg*x61300-dr*x61300-nl*x61300+ne*x61310;
};
double rhsx61310(double gto, double x00011, double x00110, double x00210, double x01010, double x41310, double x60310, double x61110, double x61210, double x61300, double x61310, double x61311)
{
	return gto*x41310+ar*Nf*x01010*x60310+agp*Nf*x00210*x61110+ac*Nf*x00110*x61210+nl*x61300-dc*x61310-dg*x61310-dr*x61310-ne*x61310-bbin*Nf*x00011*x61310+unbbin*x61311;
};
double rhsx61311(double gto, double x00011, double x00110, double x00210, double x01010, double x01011, double x41311, double x60310, double x60311, double x61111, double x61211, double x61310, double x61311)
{
	return gto*x41311+ar*Nf*x01011*x60310+ar*Nf*x01010*x60311+agp*Nf*x00210*x61111+ac*Nf*x00110*x61211+bbin*Nf*x00011*x61310-dc*x61311-dg*x61311-2*dr*x61311-unbbin*x61311;
};
double rhsx62000(double x00100, double x00200, double x02000, double x60000, double x62000, double x62010, double x62100, double x62200)
{
	return ar*x02000*x60000-dr*x62000-nl*x62000-ac*x00100*x62000-agp*x00200*x62000+ne*x62010+dc*x62100+dg*x62200;
};
double rhsx62010(double x00011, double x00110, double x00210, double x02010, double x60010, double x62000, double x62010, double x62011, double x62110, double x62210)
{
	return ar*Nf*x02010*x60010+nl*x62000-dr*x62010-ne*x62010-bbin*Nf*x00011*x62010-ac*Nf*x00110*x62010-agp*Nf*x00210*x62010+unbbin*x62011+dc*x62110+dg*x62210;
};
double rhsx62011(double x00011, double x00110, double x00210, double x02010, double x02011, double x60010, double x60011, double x62010, double x62011, double x62111, double x62211)
{
	return ar*Nf*x02011*x60010+ar*Nf*x02010*x60011+bbin*Nf*x00011*x62010-2*dr*x62011-unbbin*x62011-ac*Nf*x00110*x62011-agp*Nf*x00210*x62011+dc*x62111+dg*x62211;
};
double rhsx62100(double x00100, double x00200, double x02000, double x60100, double x62000, double x62100, double x62110, double x62300)
{
	return ar*x02000*x60100+ac*x00100*x62000-dc*x62100-dr*x62100-nl*x62100-agp*x00200*x62100+ne*x62110+dg*x62300;
};
double rhsx62110(double x00011, double x00110, double x00210, double x02010, double x60110, double x62010, double x62100, double x62110, double x62111, double x62310)
{
	return ar*Nf*x02010*x60110+ac*Nf*x00110*x62010+nl*x62100-dc*x62110-dr*x62110-ne*x62110-bbin*Nf*x00011*x62110-agp*Nf*x00210*x62110+unbbin*x62111+dg*x62310;
};
double rhsx62111(double x00011, double x00110, double x00210, double x02010, double x02011, double x60110, double x60111, double x62011, double x62110, double x62111, double x62311)
{
	return ar*Nf*x02011*x60110+ar*Nf*x02010*x60111+ac*Nf*x00110*x62011+bbin*Nf*x00011*x62110-dc*x62111-2*dr*x62111-unbbin*x62111-agp*Nf*x00210*x62111+dg*x62311;
};
double rhsx62200(double gto, double x00100, double x00200, double x02000, double x42200, double x60200, double x62000, double x62200, double x62210, double x62300)
{
	return gto*x42200+ar*x02000*x60200+agp*x00200*x62000-dg*x62200-dr*x62200-nl*x62200-ac*x00100*x62200+ne*x62210+dc*x62300;
};
double rhsx62210(double gto, double x00011, double x00110, double x00210, double x02010, double x42210, double x60210, double x62010, double x62200, double x62210, double x62211, double x62310)
{
	return gto*x42210+ar*Nf*x02010*x60210+agp*Nf*x00210*x62010+nl*x62200-dg*x62210-dr*x62210-ne*x62210-bbin*Nf*x00011*x62210-ac*Nf*x00110*x62210+unbbin*x62211+dc*x62310;
};
double rhsx62211(double gto, double x00011, double x00110, double x00210, double x02010, double x02011, double x42211, double x60210, double x60211, double x62011, double x62210, double x62211, double x62311)
{
	return gto*x42211+ar*Nf*x02011*x60210+ar*Nf*x02010*x60211+agp*Nf*x00210*x62011+bbin*Nf*x00011*x62210-dg*x62211-2*dr*x62211-unbbin*x62211-ac*Nf*x00110*x62211+dc*x62311;
};
double rhsx62300(double gto, double x00100, double x00200, double x02000, double x42300, double x60300, double x62100, double x62200, double x62300, double x62310)
{
	return gto*x42300+ar*x02000*x60300+agp*x00200*x62100+ac*x00100*x62200-dc*x62300-dg*x62300-dr*x62300-nl*x62300+ne*x62310;
};
double rhsx62310(double gto, double x00011, double x00110, double x00210, double x02010, double x42310, double x60310, double x62110, double x62210, double x62300, double x62310, double x62311)
{
     return gto*x42310+ar*Nf*x02010*x60310+agp*Nf*x00210*x62110+ac*Nf*x00110*x62210+nl*x62300-dc*x62310-dg*x62310-dr*x62310-ne*x62310-bbin*Nf*x00011*x62310+unbbin*x62311;
};

double rhsx62311(double gto, double x00011, double x00110, double x00210, double x02010, double x02011, double x42311, double x60310, double x60311, double x62111, double x62211, double x62310, double x62311)
{
     return gto*x42311+ar*Nf*x02011*x60310+ar*Nf*x02010*x60311+agp*Nf*x00210*x62111+ac*Nf*x00110*x62211+bbin*Nf*x00011*x62310-dc*x62311-dg*x62311-2*dr*x62311-unbbin*x62311;
};
