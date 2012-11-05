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
void rhs(float t, float statevector[], float output[]); 

#endif
