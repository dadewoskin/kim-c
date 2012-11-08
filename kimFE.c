// PROGRAM: kimFE
// MODEL: Kim & Forger, MSB 2012 single cell circadian clock model
// METHOD: solved using forward Euler method - OPENMP version

// INPUT: 1 argument = number of threads to use
// OUTPUT: text file with time series for all variables
// FILES: rhs contains user supplied routines that compute right hand derivatives (Eulers method)

#include <stdio.h> 
#include <math.h> 
#include <stdlib.h>
#include "rhs.h"
#include "omp.h"
#include <time.h> 
#include <string.h> 

int main(int argc, char *argv[])
{
	/* Read in number of threads to use */
	int n_th;
	if (argc != 2){
		printf("usage: %s num_threads\n",argv[0]);
		return -1;
	}
	else {
		n_th = atoi(argv[1]);
	} 

	/* Set initial conditions */
	double GR = 0.787955;
	double G = 0.210048;
	double GrR = 0.558565;
	double Gr = 0.416217;
	double GcR = 0.804404;
	double Gc = 0.193886;
	double GBR = 0.39359;
	double GB = 0.599965;
	double GBRb = 0.404247;
	double GBb = 0.589308;
	double MnPo = 4.18947;
	double McPo = 0.579108;
	double MnPt = 8.34963;
	double McPt = 1.41168;
	double MnRt = 12.3735;
	double McRt = 4.47054;
	double MnRev = 31.869;
	double McRev = 161.144;
	double MnRo = 3.76085;
	double McRo = 1.09547;
	double MnB = 35.1842;
	double McB = 8.56043;
	double MnNp = 0.520367;
	double McNp = 0.25596;
	double B = 23.7741;
	double Cl = 48.0977;
	double BC = 6.5815;
	double cyrev = 123.95;
	double revn = 467.844;
	double cyrevg = 0.193289;
	double revng = 58.2391;
	double cyrevgp = 0.0397137;
	double revngp = 19.0662;
	double cyrevp = 0.11012;
	double revnp = 67.6175;
	double gto = 1.13837;
	double x00001 = 0.342061;
	double x00011 = 8.17977;
	double x00100 = 50.5794;
	double x00110 = 3.08486;
	double x00200 = 0.0985802;
	double x00210 = 1.30846;
	double x01000 = 41.1658;
	double x01010 = 0.0156084;
	double x01011 = 8.72089;
	double x02000 = 46.5721;
	double x02010 = 0.0031953;
	double x02011 = 1.36799;
	double x10000 = 0.422898;
	double x10100 = 1.09254;
	double x20000 = 0.00181722;
	double x20010 = 5.19698e-05;
	double x20011 = 0.00512223;
	double x20100 = 0.106575;
	double x20110 = 0.000833073;
	double x20111 = 0.0566636;
	double x21000 = 0.00329174;
	double x21010 = 1.4996e-05;
	double x21011 = 0.00543134;
	double x21100 = 0.0876909;
	double x21110 = 0.000333744;
	double x21111 = 0.0548927;
	double x22000 = 0.00348546;
	double x22010 = 1.54312e-05;
	double x22011 = 0.00510457;
	double x22100 = 0.0938758;
	double x22110 = 0.000353611;
	double x22111 = 0.0535918;
	double x30000 = 0.914365;
	double x30100 = 0.664039;
	double x30200 = 0.0226162;
	double x30300 = 0.0241531;
	double x40000 = 0.00356643;
	double x40010 = 1.94469e-05;
	double x40011 = 0.000959363;
	double x40100 = 0.258442;
	double x40110 = 0.000985386;
	double x40111 = 0.0249991;
	double x40200 = 0.000157103;
	double x40210 = 1.1242e-05;
	double x40211 = 0.00132337;
	double x40300 = 0.0105914;
	double x40310 = 0.00023062;
	double x40311 = 0.025311;
	double x41000 = 0.00390694;
	double x41010 = 1.34882e-05;
	double x41011 = 0.000809309;
	double x41100 = 0.124712;
	double x41110 = 0.000416813;
	double x41111 = 0.0148919;
	double x41200 = 0.000188833;
	double x41210 = 1.79979e-06;
	double x41211 = 0.00109293;
	double x41300 = 0.00584769;
	double x41310 = 4.35837e-05;
	double x41311 = 0.0166107;
	double x42000 = 0.00428454;
	double x42010 = 1.47716e-05;
	double x42011 = 0.00086132;
	double x42100 = 0.137696;
	double x42110 = 0.000459811;
	double x42111 = 0.0160349;
	double x42200 = 0.0002068;
	double x42210 = 1.94018e-06;
	double x42211 = 0.0011516;
	double x42300 = 0.00644626;
	double x42310 = 4.76349e-05;
	double x42311 = 0.0177021;
	double x50000 = 0.00569806;
	double x50010 = 0.000190832;
	double x50011 = 0.0566834;
	double x50100 = 0.0287569;
	double x50110 = 0.000413832;
	double x50111 = 0.0494651;
	double x50200 = 0.00414505;
	double x50210 = 0.000386223;
	double x50211 = 0.124464;
	double x50300 = 0.00740796;
	double x50310 = 0.000654378;
	double x50311 = 0.104089;
	double x51000 = 0.00666974;
	double x51010 = 2.54221e-05;
	double x51011 = 0.00484524;
	double x51100 = 0.0828541;
	double x51110 = 0.000294564;
	double x51111 = 0.0323573;
	double x51200 = 0.0012142;
	double x51210 = 1.1508e-05;
	double x51211 = 0.00976858;
	double x51300 = 0.0087474;
	double x51310 = 7.98494e-05;
	double x51311 = 0.0607878;
	double x52000 = 0.00706385;
	double x52010 = 2.64292e-05;
	double x52011 = 0.00442251;
	double x52100 = 0.0852703;
	double x52110 = 0.000300906;
	double x52111 = 0.0301143;
	double x52200 = 0.00132641;
	double x52210 = 1.12497e-05;
	double x52211 = 0.00878991;
	double x52300 = 0.0091505;
	double x52310 = 7.78605e-05;
	double x52311 = 0.055586;
	double x60000 = 0.000792825;
	double x60010 = 2.57804e-05;
	double x60011 = 0.00334582;
	double x60100 = 0.0274125;
	double x60110 = 0.000385498;
	double x60111 = 0.0589937;
	double x60200 = 0.000125512;
	double x60210 = 4.99769e-05;
	double x60211 = 0.00735186;
	double x60300 = 0.00550233;
	double x60310 = 0.000631987;
	double x60311 = 0.12768;
	double x61000 = 0.00252101;
	double x61010 = 1.00779e-05;
	double x61011 = 0.00248103;
	double x61100 = 0.0581233;
	double x61110 = 0.000207429;
	double x61111 = 0.0239137;
	double x61200 = 0.000305307;
	double x61210 = 4.90217e-06;
	double x61211 = 0.00521949;
	double x61300 = 0.00753943;
	double x61310 = 6.39585e-05;
	double x61311 = 0.0471435;
	double x62000 = 0.00259772;
	double x62010 = 1.01709e-05;
	double x62011 = 0.00224807;
	double x62100 = 0.0602479;
	double x62110 = 0.000213449;
	double x62111 = 0.0225783;
	double x62200 = 0.000320067;
	double x62210 = 4.60372e-06;
	double x62211 = 0.004701;
	double x62300 = 0.00797183;
	double x62310 = 6.35261e-05;
	double x62311 = 0.0440794;
	///////////////////////////
	double GR_dot, G_dot, GrR_dot, Gr_dot, GcR_dot, Gc_dot, GBR_dot, GB_dot, GBRb_dot, GBb_dot, MnPo_dot, McPo_dot, MnPt_dot, McPt_dot, MnRt_dot, McRt_dot, MnRev_dot, McRev_dot, MnRo_dot, McRo_dot, MnB_dot, McB_dot, MnNp_dot, McNp_dot, B_dot, Cl_dot, BC_dot, cyrev_dot, revn_dot, cyrevg_dot, revng_dot, cyrevgp_dot, revngp_dot, cyrevp_dot, revnp_dot, gto_dot, x00001_dot, x00011_dot, x00100_dot, x00110_dot, x00200_dot, x00210_dot, x01000_dot, x01010_dot, x01011_dot, x02000_dot, x02010_dot, x02011_dot, x10000_dot, x10100_dot, x20000_dot, x20010_dot, x20011_dot, x20100_dot, x20110_dot, x20111_dot, x21000_dot, x21010_dot, x21011_dot, x21100_dot, x21110_dot, x21111_dot, x22000_dot, x22010_dot, x22011_dot, x22100_dot, x22110_dot, x22111_dot, x30000_dot, x30100_dot, x30200_dot, x30300_dot, x40000_dot, x40010_dot, x40011_dot, x40100_dot, x40110_dot, x40111_dot, x40200_dot, x40210_dot, x40211_dot, x40300_dot, x40310_dot, x40311_dot, x41000_dot, x41010_dot, x41011_dot, x41100_dot, x41110_dot, x41111_dot, x41200_dot, x41210_dot, x41211_dot, x41300_dot, x41310_dot, x41311_dot, x42000_dot, x42010_dot, x42011_dot, x42100_dot, x42110_dot, x42111_dot, x42200_dot, x42210_dot, x42211_dot, x42300_dot, x42310_dot, x42311_dot, x50000_dot, x50010_dot, x50011_dot, x50100_dot, x50110_dot, x50111_dot, x50200_dot, x50210_dot, x50211_dot, x50300_dot, x50310_dot, x50311_dot, x51000_dot, x51010_dot, x51011_dot, x51100_dot, x51110_dot, x51111_dot, x51200_dot, x51210_dot, x51211_dot, x51300_dot, x51310_dot, x51311_dot, x52000_dot, x52010_dot, x52011_dot, x52100_dot, x52110_dot, x52111_dot, x52200_dot, x52210_dot, x52211_dot, x52300_dot, x52310_dot, x52311_dot, x60000_dot, x60010_dot, x60011_dot, x60100_dot, x60110_dot, x60111_dot, x60200_dot, x60210_dot, x60211_dot, x60300_dot, x60310_dot, x60311_dot, x61000_dot, x61010_dot, x61011_dot, x61100_dot, x61110_dot, x61111_dot, x61200_dot, x61210_dot, x61211_dot, x61300_dot, x61310_dot, x61311_dot, x62000_dot, x62010_dot, x62011_dot, x62100_dot, x62110_dot, x62111_dot, x62200_dot, x62210_dot, x62211_dot, x62300_dot, x62310_dot, x62311_dot;

	time_t start; 
	start=time(NULL); 
	
	double t1=0; /* initial time */
	double t2=500; /* end time */
	double t=t1;
	int nstep=(t2-t1)/step_size; /* number of identical steps to take */

	int i,j,t_step;
	double record_size=0.1; /* record data to file every record_size hours (if step_size <= record_size)*/
	int record; /* record data to file every record time steps */
	if (step_size<=record_size) record = round(record_size/step_size);
	else record = 1;

//////////////////////// For Data Output /////////////////////////////
	char filename1[50]; 
	FILE *outfi1; 

	sprintf(filename1,"./Results/kim%d.txt", n_th); /* records output every record timesteps */
	remove(filename1); 
	if ((outfi1 = fopen(filename1, "w")) == NULL) { 
	  printf("can not open output file"); 
	};/* opens an output file*/ 

	/* Output initial conditions */
	//Header
	fprintf(outfi1, "number of threads: %d\n", n_th);
	fprintf(outfi1, "t\t GR\t G\t GrR\t Gr\t GcR\t Gc\t GBR\t GB\t GBRb\t GBb\t MnPo\t McPo\t MnPt\t McPt\t MnRt\t McRt\t MnRev\t McRev\t MnRo\t McRo\t MnB\t McB\t MnNp\t McNp\t B\t Cl\t BC\t cyrev\t revn\t cyrevg\t revng\t cyrevgp\t revngp\t cyrevp\t revnp\t gto\t x00001\t x00011\t x00100\t x00110\t x00200\t x00210\t x01000\t x01010\t x01011\t x02000\t x02010\t x02011\t x10000\t x10100\t x20000\t x20010\t x20011\t x20100\t x20110\t x20111\t x21000\t x21010\t x21011\t x21100\t x21110\t x21111\t x22000\t x22010\t x22011\t x22100\t x22110\t x22111\t x30000\t x30100\t x30200\t x30300\t x40000\t x40010\t x40011\t x40100\t x40110\t x40111\t x40200\t x40210\t x40211\t x40300\t x40310\t x40311\t x41000\t x41010\t x41011\t x41100\t x41110\t x41111\t x41200\t x41210\t x41211\t x41300\t x41310\t x41311\t x42000\t x42010\t x42011\t x42100\t x42110\t x42111\t x42200\t x42210\t x42211\t x42300\t x42310\t x42311\t x50000\t x50010\t x50011\t x50100\t x50110\t x50111\t x50200\t x50210\t x50211\t x50300\t x50310\t x50311\t x51000\t x51010\t x51011\t x51100\t x51110\t x51111\t x51200\t x51210\t x51211\t x51300\t x51310\t x51311\t x52000\t x52010\t x52011\t x52100\t x52110\t x52111\t x52200\t x52210\t x52211\t x52300\t x52310\t x52311\t x60000\t x60010\t x60011\t x60100\t x60110\t x60111\t x60200\t x60210\t x60211\t x60300\t x60310\t x60311\t x61000\t x61010\t x61011\t x61100\t x61110\t x61111\t x61200\t x61210\t x61211\t x61300\t x61310\t x61311\t x62000\t x62010\t x62011\t x62100\t x62110\t x62111\t x62200\t x62210\t x62211\t x62300\t x62310\t x62311\n");
	//Initial data
	fprintf(outfi1, "%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f\n", t, GR, G, GrR, Gr, GcR, Gc, GBR, GB, GBRb, GBb, MnPo, McPo, MnPt, McPt, MnRt, McRt, MnRev, McRev, MnRo, McRo, MnB, McB, MnNp, McNp, B, Cl, BC, cyrev, revn, cyrevg, revng, cyrevgp, revngp, cyrevp, revnp, gto, x00001, x00011, x00100, x00110, x00200, x00210, x01000, x01010, x01011, x02000, x02010, x02011, x10000, x10100, x20000, x20010, x20011, x20100, x20110, x20111, x21000, x21010, x21011, x21100, x21110, x21111, x22000, x22010, x22011, x22100, x22110, x22111, x30000, x30100, x30200, x30300, x40000, x40010, x40011, x40100, x40110, x40111, x40200, x40210, x40211, x40300, x40310, x40311, x41000, x41010, x41011, x41100, x41110, x41111, x41200, x41210, x41211, x41300, x41310, x41311, x42000, x42010, x42011, x42100, x42110, x42111, x42200, x42210, x42211, x42300, x42310, x42311, x50000, x50010, x50011, x50100, x50110, x50111, x50200, x50210, x50211, x50300, x50310, x50311, x51000, x51010, x51011, x51100, x51110, x51111, x51200, x51210, x51211, x51300, x51310, x51311, x52000, x52010, x52011, x52100, x52110, x52111, x52200, x52210, x52211, x52300, x52310, x52311, x60000, x60010, x60011, x60100, x60110, x60111, x60200, x60210, x60211, x60300, x60310, x60311, x61000, x61010, x61011, x61100, x61110, x61111, x61200, x61210, x61211, x61300, x61310, x61311, x62000, x62010, x62011, x62100, x62110, x62111, x62200, x62210, x62211, x62300, x62310, x62311);
	
/////////////////////// Run computation //////////////////////////////
//	omp_set_num_threads(n_th);
	for (t_step=1; t_step<=nstep; t_step++) {
		t=t_step*step_size;

	/* Calculate derivatives */
		#pragma omp parallel
		{
			#pragma omp sections
			{
				#pragma omp section
				{
					GR_dot=rhsGR(GR, G, x01011, x02011);
				}
				#pragma omp section
				{
					G_dot=rhsG(GR, G, x00011);
				}
				#pragma omp section
				{
					GrR_dot=rhsGrR(G, GrR, Gr, x01011, x02011);
				}
				#pragma omp section
				{
					Gr_dot=rhsGr(G, GrR, Gr, x00011);
				}
				#pragma omp section
				{
					GcR_dot=rhsGcR(G, GcR, Gc, x01011, x02011);
				}
				#pragma omp section
				{
					Gc_dot=rhsGc(G, GcR, Gc, x00011);
				}
				#pragma omp section
				{
					GBR_dot=rhsGBR(G, GBR, GB, B, revn, revng, revngp, revnp);
				}
				#pragma omp section
				{
					GB_dot=rhsGB(G, GBR, GB, B, revn, revng, revngp, revnp);
				}
				#pragma omp section
				{
					GBRb_dot=rhsGBRb(G, GBR, GB, GBRb, GBb, B, revn, revng, revngp, revnp);
				}
				#pragma omp section
				{
					GBb_dot=rhsGBb(G, GBR, GB, GBRb, GBb, B, revn, revng, revngp, revnp);
				}
				#pragma omp section
				{
					MnPo_dot=rhsMnPo(G, MnPo);
				}
				#pragma omp section
				{
					McPo_dot=rhsMcPo(MnPo, McPo);
				}
				#pragma omp section
				{
					MnPt_dot=rhsMnPt(G, MnPt);
				}
				#pragma omp section
				{
					McPt_dot=rhsMcPt(MnPt, McPt);
				}
				#pragma omp section
				{
					MnRt_dot=rhsMnRt(G, Gc, MnRt);
				}
				#pragma omp section
				{
					McRt_dot=rhsMcRt(MnRt, McRt);
				}
	#pragma omp section
				{
					MnRev_dot=rhsMnRev(G, Gr, MnRev, x00011);
				}
	#pragma omp section
				{
					McRev_dot=rhsMcRev(MnRev, McRev);
				}
	#pragma omp section
				{
					MnRo_dot=rhsMnRo(G, GB, MnRo, B);
				}
	#pragma omp section
				{
					McRo_dot=rhsMcRo(MnRo, McRo);
				}
	#pragma omp section
				{
					MnB_dot=rhsMnB(G, GB, GBb, MnB, B);
				}
	#pragma omp section
				{
					McB_dot=rhsMcB(MnB, McB, B);
				}
	#pragma omp section
				{
					MnNp_dot=rhsMnNp(G, GB, MnNp, B);
				}
	#pragma omp section
				{
					McNp_dot=rhsMcNp(MnNp, McNp);
				}
	#pragma omp section
				{
					B_dot=rhsB(McB, B, Cl, BC);
				}
	#pragma omp section
				{
					Cl_dot=rhsCl(McNp, B, Cl, BC);
				}
	#pragma omp section
				{
					BC_dot=rhsBC(B, Cl, BC);
				}
	#pragma omp section
				{
					cyrev_dot=rhscyrev(McRev, cyrev, revn, cyrevg, x00200);
				}
	#pragma omp section
				{
					revn_dot=rhsrevn(cyrev, revn, revng, x00210);
				}
	#pragma omp section
				{
					cyrevg_dot=rhscyrevg(cyrev, revn, cyrevg, revng, gto, x00200);
				}
	#pragma omp section
				{
					revng_dot=rhsrevng(cyrev, revn, cyrevg, revng, gto, x00210);
				}
	#pragma omp section
				{
					cyrevgp_dot=rhscyrevgp(cyrev, revn, cyrevg, revng, cyrevgp, revngp, gto);
				}
	#pragma omp section
				{
					revngp_dot=rhsrevngp(cyrev, revn, cyrevg, revng, cyrevgp, revngp, gto);
				}
	#pragma omp section
				{
					cyrevp_dot=rhscyrevp(cyrev, revn, cyrevg, cyrevgp, cyrevp, revnp);
				}
	#pragma omp section
				{
					revnp_dot=rhsrevnp(cyrev, revn, revng, revngp, cyrevp, revnp);
				}
	#pragma omp section
				{
					gto_dot=rhsgto(G, GB, B, gto);
				}
	#pragma omp section
				{
					x00001_dot=rhsx00001(B, BC, x00001);
				}
	#pragma omp section
				{
					x00011_dot=rhsx00011(x00001, x00011, x01010, x01011, x02010, x02011, x20010, x20011, x20110, x20111, x21010, x21011, x21110, x21111, x22010, x22011, x22110, x22111, x40010, x40011, x40110, x40111, x40210, x40211, x40310, x40311, x41010, x41011, x41110, x41111, x41210, x41211, x41310, x41311, x42010, x42011, x42110, x42111, x42210, x42211, x42310, x42311, x50010, x50011, x50110, x50111, x50210, x50211, x50310, x50311, x51010, x51011, x51110, x51111, x51210, x51211, x51310, x51311, x52010, x52011, x52110, x52111, x52210, x52211, x52310, x52311, x60010, x60011, x60110, x60111, x60210, x60211, x60310, x60311, x61010, x61011, x61110, x61111, x61210, x61211, x61310, x61311, x62010, x62011, x62110, x62111, x62210, x62211, x62310, x62311);
				}
	#pragma omp section
				{
					x00100_dot=rhsx00100(x00100, x00110, x10000, x10100, x20000, x20100, x21000, x21100, x22000, x22100, x30000, x30100, x30200, x30300, x40000, x40100, x40200, x40300, x41000, x41100, x41200, x41300, x42000, x42100, x42200, x42300, x50000, x50100, x50200, x50300, x51000, x51100, x51200, x51300, x52000, x52100, x52200, x52300, x60000, x60100, x60200, x60300, x61000, x61100, x61200, x61300, x62000, x62100, x62200, x62300);
				}
	#pragma omp section
				{
					x00110_dot=rhsx00110(x00110, x20010, x20011, x20110, x20111, x21010, x21011, x21110, x21111, x22010, x22011, x22110, x22111, x40010, x40011, x40110, x40111, x40210, x40211, x40310, x40311, x41010, x41011, x41110, x41111, x41210, x41211, x41310, x41311, x42010, x42011, x42110, x42111, x42210, x42211, x42310, x42311, x50010, x50011, x50110, x50111, x50210, x50211, x50310, x50311, x51010, x51011, x51110, x51111, x51210, x51211, x51310, x51311, x52010, x52011, x52110, x52111, x52210, x52211, x52310, x52311, x60010, x60011, x60110, x60111, x60210, x60211, x60310, x60311, x61010, x61011, x61110, x61111, x61210, x61211, x61310, x61311, x62010, x62011, x62110, x62111, x62210, x62211, x62310, x62311);
				}
	#pragma omp section
				{
					x00200_dot=rhsx00200(cyrev, cyrevg, cyrevgp, x00200, x00210, x30000, x30100, x30200, x30300, x40000, x40100, x40200, x40300, x41000, x41100, x41200, x41300, x42000, x42100, x42200, x42300, x50000, x50100, x50200, x50300, x51000, x51100, x51200, x51300, x52000, x52100, x52200, x52300, x60000, x60100, x60200, x60300, x61000, x61100, x61200, x61300, x62000, x62100, x62200, x62300);
				}
	#pragma omp section
				{
					x00210_dot=rhsx00210(revn, revng, revngp, x00210, x40010, x40011, x40110, x40111, x40210, x40211, x40310, x40311, x41010, x41011, x41110, x41111, x41210, x41211, x41310, x41311, x42010, x42011, x42110, x42111, x42210, x42211, x42310, x42311, x50010, x50011, x50110, x50111, x50210, x50211, x50310, x50311, x51010, x51011, x51110, x51111, x51210, x51211, x51310, x51311, x52010, x52011, x52110, x52111, x52210, x52211, x52310, x52311, x60010, x60011, x60110, x60111, x60210, x60211, x60310, x60311, x61010, x61011, x61110, x61111, x61210, x61211, x61310, x61311, x62010, x62011, x62110, x62111, x62210, x62211, x62310, x62311);
				}
	#pragma omp section
				{
					x01000_dot=rhsx01000(McRo, x01000, x20000, x20100, x21000, x21100, x40000, x40100, x40200, x40300, x41000, x41100, x41200, x41300, x50000, x50100, x50200, x50300, x51000, x51100, x51200, x51300, x60000, x60100, x60200, x60300, x61000, x61100, x61200, x61300);
				}
	#pragma omp section
				{
					x01010_dot=rhsx01010(x00011, x01010, x01011, x20010, x20011, x20110, x20111, x21010, x21011, x21110, x21111, x40010, x40011, x40110, x40111, x40210, x40211, x40310, x40311, x41010, x41011, x41110, x41111, x41210, x41211, x41310, x41311, x50010, x50011, x50110, x50111, x50210, x50211, x50310, x50311, x51010, x51011, x51110, x51111, x51210, x51211, x51310, x51311, x60010, x60011, x60110, x60111, x60210, x60211, x60310, x60311, x61010, x61011, x61110, x61111, x61210, x61211, x61310, x61311);
				}
	#pragma omp section
				{
					x01011_dot=rhsx01011(x00011, x01010, x01011, x20010, x20110, x21011, x21111, x40010, x40110, x40210, x40310, x41011, x41111, x41211, x41311, x50010, x50110, x50210, x50310, x51011, x51111, x51211, x51311, x60010, x60110, x60210, x60310, x61011, x61111, x61211, x61311);
				}
	#pragma omp section
				{
					x02000_dot=rhsx02000(McRt, x02000, x20000, x20100, x22000, x22100, x40000, x40100, x40200, x40300, x42000, x42100, x42200, x42300, x50000, x50100, x50200, x50300, x52000, x52100, x52200, x52300, x60000, x60100, x60200, x60300, x62000, x62100, x62200, x62300);
				}
	#pragma omp section
				{
					x02010_dot=rhsx02010(x00011, x02010, x02011, x20010, x20011, x20110, x20111, x22010, x22011, x22110, x22111, x40010, x40011, x40110, x40111, x40210, x40211, x40310, x40311, x42010, x42011, x42110, x42111, x42210, x42211, x42310, x42311, x50010, x50011, x50110, x50111, x50210, x50211, x50310, x50311, x52010, x52011, x52110, x52111, x52210, x52211, x52310, x52311, x60010, x60011, x60110, x60111, x60210, x60211, x60310, x60311, x62010, x62011, x62110, x62111, x62210, x62211, x62310, x62311);
				}
	#pragma omp section
				{
					x02011_dot=rhsx02011(x00011, x02010, x02011, x20010, x20110, x22011, x22111, x40010, x40110, x40210, x40310, x42011, x42111, x42211, x42311, x50010, x50110, x50210, x50310, x52011, x52111, x52211, x52311, x60010, x60110, x60210, x60310, x62011, x62111, x62211, x62311);
				}
	#pragma omp section
				{
					x10000_dot=rhsx10000(McPo, x00100, x10000, x10100);
				}
	#pragma omp section
				{
					x10100_dot=rhsx10100(x00100, x10000, x10100);
				}
	#pragma omp section
				{
					x20000_dot=rhsx20000(x00100, x01000, x02000, x20000, x20010, x20100, x21000, x22000);
				}
	#pragma omp section
				{
					x20010_dot=rhsx20010(x00011, x00110, x01010, x01011, x02010, x02011, x20000, x20010, x20011, x20110, x21010, x21011, x22010, x22011);
				}
	#pragma omp section
				{
					x20011_dot=rhsx20011(x00011, x00110, x01010, x02010, x20010, x20011, x20111, x21011, x22011);
				}
	#pragma omp section
				{
					x20100_dot=rhsx20100(x00100, x01000, x02000, x10100, x20000, x20100, x20110, x21100, x22100);
				}
	#pragma omp section
				{
					x20110_dot=rhsx20110(x00011, x00110, x01010, x01011, x02010, x02011, x20010, x20100, x20110, x20111, x21110, x21111, x22110, x22111);
				}
	#pragma omp section
				{
					x20111_dot=rhsx20111(x00011, x00110, x01010, x02010, x20011, x20110, x20111, x21111, x22111);
				}
	#pragma omp section
				{
					x21000_dot=rhsx21000(x00100, x01000, x20000, x21000, x21010, x21100);
				}
	#pragma omp section
				{
					x21010_dot=rhsx21010(x00011, x00110, x01010, x20010, x21000, x21010, x21011, x21110);
				}
	#pragma omp section
				{
					x21011_dot=rhsx21011(x00011, x00110, x01010, x01011, x20010, x20011, x21010, x21011, x21111);
				}
	#pragma omp section
				{
					x21100_dot=rhsx21100(x00100, x01000, x20100, x21000, x21100, x21110);
				}
	#pragma omp section
				{
					x21110_dot=rhsx21110(x00011, x00110, x01010, x20110, x21010, x21100, x21110, x21111);
				}
	#pragma omp section
				{
					x21111_dot=rhsx21111(x00011, x00110, x01010, x01011, x20110, x20111, x21011, x21110, x21111);
				}
	#pragma omp section
				{
					x22000_dot=rhsx22000(x00100, x02000, x20000, x22000, x22010, x22100);
				}
	#pragma omp section
				{
					x22010_dot=rhsx22010(x00011, x00110, x02010, x20010, x22000, x22010, x22011, x22110);
				}
	#pragma omp section
				{
					x22011_dot=rhsx22011(x00011, x00110, x02010, x02011, x20010, x20011, x22010, x22011, x22111);
				}
	#pragma omp section
				{
					x22100_dot=rhsx22100(x00100, x02000, x20100, x22000, x22100, x22110);
				}
	#pragma omp section
				{
					x22110_dot=rhsx22110(x00011, x00110, x02010, x20110, x22010, x22100, x22110, x22111);
				}
	#pragma omp section
				{
					x22111_dot=rhsx22111(x00011, x00110, x02010, x02011, x20110, x20111, x22011, x22110, x22111);
				}
	#pragma omp section
				{
					x30000_dot=rhsx30000(McPt, x00100, x00200, x30000, x30100, x30200);
				}
	#pragma omp section
				{
					x30100_dot=rhsx30100(x00100, x00200, x30000, x30100, x30300);
				}
	#pragma omp section
				{
					x30200_dot=rhsx30200(gto, x00100, x00200, x30000, x30200, x30300);
				}
	#pragma omp section
				{
					x30300_dot=rhsx30300(gto, x00100, x00200, x30100, x30200, x30300);
				}
	#pragma omp section
				{
					x40000_dot=rhsx40000(x00100, x00200, x01000, x02000, x40000, x40010, x40100, x40200, x41000, x42000);
				}
	#pragma omp section
				{
					x40010_dot=rhsx40010(x00011, x00110, x00210, x01010, x01011, x02010, x02011, x40000, x40010, x40011, x40110, x40210, x41010, x41011, x42010, x42011);
				}
	#pragma omp section
				{
					x40011_dot=rhsx40011(x00011, x00110, x00210, x01010, x02010, x40010, x40011, x40111, x40211, x41011, x42011);
				}
	#pragma omp section
				{
					x40100_dot=rhsx40100(x00100, x00200, x01000, x02000, x30100, x40000, x40100, x40110, x40300, x41100, x42100);
				}
	#pragma omp section
				{
					x40110_dot=rhsx40110(x00011, x00110, x00210, x01010, x01011, x02010, x02011, x40010, x40100, x40110, x40111, x40310, x41110, x41111, x42110, x42111);
				}
	#pragma omp section
				{
					x40111_dot=rhsx40111(x00011, x00110, x00210, x01010, x02010, x40011, x40110, x40111, x40311, x41111, x42111);
				}
	#pragma omp section
				{
					x40200_dot=rhsx40200(gto, x00100, x00200, x01000, x02000, x40000, x40200, x40210, x40300, x41200, x42200);
				}
	#pragma omp section
				{
					x40210_dot=rhsx40210(gto, x00011, x00110, x00210, x01010, x01011, x02010, x02011, x40010, x40200, x40210, x40211, x40310, x41210, x41211, x42210, x42211);
				}
	#pragma omp section
				{
					x40211_dot=rhsx40211(gto, x00011, x00110, x00210, x01010, x02010, x40011, x40210, x40211, x40311, x41211, x42211);
				}
	#pragma omp section
				{
					x40300_dot=rhsx40300(gto, x00100, x00200, x01000, x02000, x30300, x40100, x40200, x40300, x40310, x41300, x42300);
				}
	#pragma omp section
				{
					x40310_dot=rhsx40310(gto, x00011, x00110, x00210, x01010, x01011, x02010, x02011, x40110, x40210, x40300, x40310, x40311, x41310, x41311, x42310, x42311);
				}
	#pragma omp section
				{
					x40311_dot=rhsx40311(gto, x00011, x00110, x00210, x01010, x02010, x40111, x40211, x40310, x40311, x41311, x42311);
				}
	#pragma omp section
				{
					x41000_dot=rhsx41000(x00100, x00200, x01000, x40000, x41000, x41010, x41100, x41200);
				}
	#pragma omp section
				{
					x41010_dot=rhsx41010(x00011, x00110, x00210, x01010, x40010, x41000, x41010, x41011, x41110, x41210);
				}
	#pragma omp section
				{
					x41011_dot=rhsx41011(x00011, x00110, x00210, x01010, x01011, x40010, x40011, x41010, x41011, x41111, x41211);
				}
	#pragma omp section
				{
					x41100_dot=rhsx41100(x00100, x00200, x01000, x40100, x41000, x41100, x41110, x41300);
				}
	#pragma omp section
				{
					x41110_dot=rhsx41110(x00011, x00110, x00210, x01010, x40110, x41010, x41100, x41110, x41111, x41310);
				}
	#pragma omp section
				{
					x41111_dot=rhsx41111(x00011, x00110, x00210, x01010, x01011, x40110, x40111, x41011, x41110, x41111, x41311);
				}
	#pragma omp section
				{
					x41200_dot=rhsx41200(gto, x00100, x00200, x01000, x40200, x41000, x41200, x41210, x41300);
				}
	#pragma omp section
				{
					x41210_dot=rhsx41210(gto, x00011, x00110, x00210, x01010, x40210, x41010, x41200, x41210, x41211, x41310);
				}
	#pragma omp section
				{
					x41211_dot=rhsx41211(gto, x00011, x00110, x00210, x01010, x01011, x40210, x40211, x41011, x41210, x41211, x41311);
				}
	#pragma omp section
				{
					x41300_dot=rhsx41300(gto, x00100, x00200, x01000, x40300, x41100, x41200, x41300, x41310);
				}
	#pragma omp section
				{
					x41310_dot=rhsx41310(gto, x00011, x00110, x00210, x01010, x40310, x41110, x41210, x41300, x41310, x41311);
				}
	#pragma omp section
				{
					x41311_dot=rhsx41311(gto, x00011, x00110, x00210, x01010, x01011, x40310, x40311, x41111, x41211, x41310, x41311);
				}
	#pragma omp section
				{
					x42000_dot=rhsx42000(x00100, x00200, x02000, x40000, x42000, x42010, x42100, x42200);
				}
	#pragma omp section
				{
					x42010_dot=rhsx42010(x00011, x00110, x00210, x02010, x40010, x42000, x42010, x42011, x42110, x42210);
				}
	#pragma omp section
				{
					x42011_dot=rhsx42011(x00011, x00110, x00210, x02010, x02011, x40010, x40011, x42010, x42011, x42111, x42211);
				}
	#pragma omp section
				{
					x42100_dot=rhsx42100(x00100, x00200, x02000, x40100, x42000, x42100, x42110, x42300);
				}
	#pragma omp section
				{
					x42110_dot=rhsx42110(x00011, x00110, x00210, x02010, x40110, x42010, x42100, x42110, x42111, x42310);
				}
	#pragma omp section
				{
					x42111_dot=rhsx42111(x00011, x00110, x00210, x02010, x02011, x40110, x40111, x42011, x42110, x42111, x42311);
				}
	#pragma omp section
				{
					x42200_dot=rhsx42200(gto, x00100, x00200, x02000, x40200, x42000, x42200, x42210, x42300);
				}
	#pragma omp section
				{
					x42210_dot=rhsx42210(gto, x00011, x00110, x00210, x02010, x40210, x42010, x42200, x42210, x42211, x42310);
				}
	#pragma omp section
				{
					x42211_dot=rhsx42211(gto, x00011, x00110, x00210, x02010, x02011, x40210, x40211, x42011, x42210, x42211, x42311);
				}
	#pragma omp section
				{
					x42300_dot=rhsx42300(gto, x00100, x00200, x02000, x40300, x42100, x42200, x42300, x42310);
				}
	#pragma omp section
				{
					x42310_dot=rhsx42310(gto, x00011, x00110, x00210, x02010, x40310, x42110, x42210, x42300, x42310, x42311);
				}
	#pragma omp section
				{
					x42311_dot=rhsx42311(gto, x00011, x00110, x00210, x02010, x02011, x40310, x40311, x42111, x42211, x42310, x42311);
				}
	#pragma omp section
				{
					x50000_dot=rhsx50000(x00100, x00200, x01000, x02000, x50000, x50010, x50100, x50200, x51000, x52000);
				}
	#pragma omp section
				{
					x50010_dot=rhsx50010(x00011, x00110, x00210, x01010, x01011, x02010, x02011, x50000, x50010, x50011, x50110, x50210, x51010, x51011, x52010, x52011);
				}
	#pragma omp section
				{
					x50011_dot=rhsx50011(x00011, x00110, x00210, x01010, x02010, x50010, x50011, x50111, x50211, x51011, x52011);
				}
	#pragma omp section
				{
					x50100_dot=rhsx50100(x00100, x00200, x01000, x02000, x50000, x50100, x50110, x50300, x51100, x52100);
				}
	#pragma omp section
				{
					x50110_dot=rhsx50110(x00011, x00110, x00210, x01010, x01011, x02010, x02011, x50010, x50100, x50110, x50111, x50310, x51110, x51111, x52110, x52111);
				}
	#pragma omp section
				{
					x50111_dot=rhsx50111(x00011, x00110, x00210, x01010, x02010, x50011, x50110, x50111, x50311, x51111, x52111);
				}
	#pragma omp section
				{
					x50200_dot=rhsx50200(gto, x00100, x00200, x01000, x02000, x30200, x50000, x50200, x50210, x50300, x51200, x52200);
				}
	#pragma omp section
				{
					x50210_dot=rhsx50210(x00011, x00110, x00210, x01010, x01011, x02010, x02011, x50010, x50200, x50210, x50211, x50310, x51210, x51211, x52210, x52211);
				}
	#pragma omp section
				{
					x50211_dot=rhsx50211(x00011, x00110, x00210, x01010, x02010, x50011, x50210, x50211, x50311, x51211, x52211);
				}
	#pragma omp section
				{
					x50300_dot=rhsx50300(gto, x00100, x00200, x01000, x02000, x30300, x50100, x50200, x50300, x50310, x51300, x52300);
				}
	#pragma omp section
				{
					x50310_dot=rhsx50310(x00011, x00110, x00210, x01010, x01011, x02010, x02011, x50110, x50210, x50300, x50310, x50311, x51310, x51311, x52310, x52311);
				}
	#pragma omp section
				{
					x50311_dot=rhsx50311(x00011, x00110, x00210, x01010, x02010, x50111, x50211, x50310, x50311, x51311, x52311);
				}
	#pragma omp section
				{
					x51000_dot=rhsx51000(x00100, x00200, x01000, x50000, x51000, x51010, x51100, x51200);
				}
	#pragma omp section
				{
					x51010_dot=rhsx51010(x00011, x00110, x00210, x01010, x50010, x51000, x51010, x51011, x51110, x51210);
				}
	#pragma omp section
				{
					x51011_dot=rhsx51011(x00011, x00110, x00210, x01010, x01011, x50010, x50011, x51010, x51011, x51111, x51211);
				}
	#pragma omp section
				{
					x51100_dot=rhsx51100(x00100, x00200, x01000, x50100, x51000, x51100, x51110, x51300);
				}
	#pragma omp section
				{
					x51110_dot=rhsx51110(x00011, x00110, x00210, x01010, x50110, x51010, x51100, x51110, x51111, x51310);
				}
	#pragma omp section
				{
					x51111_dot=rhsx51111(x00011, x00110, x00210, x01010, x01011, x50110, x50111, x51011, x51110, x51111, x51311);
				}
	#pragma omp section
				{
					x51200_dot=rhsx51200(x00100, x00200, x01000, x50200, x51000, x51200, x51210, x51300);
				}
	#pragma omp section
				{
					x51210_dot=rhsx51210(x00011, x00110, x00210, x01010, x50210, x51010, x51200, x51210, x51211, x51310);
				}
	#pragma omp section
				{
					x51211_dot=rhsx51211(x00011, x00110, x00210, x01010, x01011, x50210, x50211, x51011, x51210, x51211, x51311);
				}
	#pragma omp section
				{
					x51300_dot=rhsx51300(x00100, x00200, x01000, x50300, x51100, x51200, x51300, x51310);
				}
	#pragma omp section
				{
					x51310_dot=rhsx51310(x00011, x00110, x00210, x01010, x50310, x51110, x51210, x51300, x51310, x51311);
				}
	#pragma omp section
				{
					x51311_dot=rhsx51311(x00011, x00110, x00210, x01010, x01011, x50310, x50311, x51111, x51211, x51310, x51311);
				}
	#pragma omp section
				{
					x52000_dot=rhsx52000(x00100, x00200, x02000, x50000, x52000, x52010, x52100, x52200);
				}
	#pragma omp section
				{
					x52010_dot=rhsx52010(x00011, x00110, x00210, x02010, x50010, x52000, x52010, x52011, x52110, x52210);
				}
	#pragma omp section
				{
					x52011_dot=rhsx52011(x00011, x00110, x00210, x02010, x02011, x50010, x50011, x52010, x52011, x52111, x52211);
				}
	#pragma omp section
				{
					x52100_dot=rhsx52100(x00100, x00200, x02000, x50100, x52000, x52100, x52110, x52300);
				}
	#pragma omp section
				{
					x52110_dot=rhsx52110(x00011, x00110, x00210, x02010, x50110, x52010, x52100, x52110, x52111, x52310);
				}
	#pragma omp section
				{
					x52111_dot=rhsx52111(x00011, x00110, x00210, x02010, x02011, x50110, x50111, x52011, x52110, x52111, x52311);
				}
	#pragma omp section
				{
					x52200_dot=rhsx52200(x00100, x00200, x02000, x50200, x52000, x52200, x52210, x52300);
				}
	#pragma omp section
				{
					x52210_dot=rhsx52210(x00011, x00110, x00210, x02010, x50210, x52010, x52200, x52210, x52211, x52310);
				}
	#pragma omp section
				{
					x52211_dot=rhsx52211(x00011, x00110, x00210, x02010, x02011, x50210, x50211, x52011, x52210, x52211, x52311);
				}
	#pragma omp section
				{
					x52300_dot=rhsx52300(x00100, x00200, x02000, x50300, x52100, x52200, x52300, x52310);
				}
	#pragma omp section
				{
					x52310_dot=rhsx52310(x00011, x00110, x00210, x02010, x50310, x52110, x52210, x52300, x52310, x52311);
				}
	#pragma omp section
				{
					x52311_dot=rhsx52311(x00011, x00110, x00210, x02010, x02011, x50310, x50311, x52111, x52211, x52310, x52311);
				}
	#pragma omp section
				{
					x60000_dot=rhsx60000(x00100, x00200, x01000, x02000, x60000, x60010, x60100, x60200, x61000, x62000);
				}
	#pragma omp section
				{
					x60010_dot=rhsx60010(x00011, x00110, x00210, x01010, x01011, x02010, x02011, x60000, x60010, x60011, x60110, x60210, x61010, x61011, x62010, x62011);
				}
	#pragma omp section
				{
					x60011_dot=rhsx60011(x00011, x00110, x00210, x01010, x02010, x60010, x60011, x60111, x60211, x61011, x62011);
				}
	#pragma omp section
				{
					x60100_dot=rhsx60100(x00100, x00200, x01000, x02000, x50100, x60000, x60100, x60110, x60300, x61100, x62100);
				}
	#pragma omp section
				{
					x60110_dot=rhsx60110(x00011, x00110, x00210, x01010, x01011, x02010, x02011, x50110, x60010, x60100, x60110, x60111, x60310, x61110, x61111, x62110, x62111);
				}
	#pragma omp section
				{
					x60111_dot=rhsx60111(x00011, x00110, x00210, x01010, x02010, x50111, x60011, x60110, x60111, x60311, x61111, x62111);
				}
	#pragma omp section
				{
					x60200_dot=rhsx60200(gto, x00100, x00200, x01000, x02000, x40200, x60000, x60200, x60210, x60300, x61200, x62200);
				}
	#pragma omp section
				{
					x60210_dot=rhsx60210(gto, x00011, x00110, x00210, x01010, x01011, x02010, x02011, x40210, x60010, x60200, x60210, x60211, x60310, x61210, x61211, x62210, x62211);
				}
	#pragma omp section
				{
					x60211_dot=rhsx60211(gto, x00011, x00110, x00210, x01010, x02010, x40211, x60011, x60210, x60211, x60311, x61211, x62211);
				}
	#pragma omp section
				{
					x60300_dot=rhsx60300(gto, x00100, x00200, x01000, x02000, x40300, x50300, x60100, x60200, x60300, x60310, x61300, x62300);
				}
	#pragma omp section
				{
					x60310_dot=rhsx60310(gto, x00011, x00110, x00210, x01010, x01011, x02010, x02011, x40310, x50310, x60110, x60210, x60300, x60310, x60311, x61310, x61311, x62310, x62311);
				}
	#pragma omp section
				{
					x60311_dot=rhsx60311(gto, x00011, x00110, x00210, x01010, x02010, x40311, x50311, x60111, x60211, x60310, x60311, x61311, x62311);
				}
	#pragma omp section
				{
					x61000_dot=rhsx61000(x00100, x00200, x01000, x60000, x61000, x61010, x61100, x61200);
				}
	#pragma omp section
				{
					x61010_dot=rhsx61010(x00011, x00110, x00210, x01010, x60010, x61000, x61010, x61011, x61110, x61210);
				}
	#pragma omp section
				{
					x61011_dot=rhsx61011(x00011, x00110, x00210, x01010, x01011, x60010, x60011, x61010, x61011, x61111, x61211);
				}
	#pragma omp section
				{
					x61100_dot=rhsx61100(x00100, x00200, x01000, x60100, x61000, x61100, x61110, x61300);
				}
	#pragma omp section
				{
					x61110_dot=rhsx61110(x00011, x00110, x00210, x01010, x60110, x61010, x61100, x61110, x61111, x61310);
				}
	#pragma omp section
				{
					x61111_dot=rhsx61111(x00011, x00110, x00210, x01010, x01011, x60110, x60111, x61011, x61110, x61111, x61311);
				}
	#pragma omp section
				{
					x61200_dot=rhsx61200(gto, x00100, x00200, x01000, x41200, x60200, x61000, x61200, x61210, x61300);
				}
	#pragma omp section
				{
					x61210_dot=rhsx61210(gto, x00011, x00110, x00210, x01010, x41210, x60210, x61010, x61200, x61210, x61211, x61310);
				}
	#pragma omp section
				{
					x61211_dot=rhsx61211(gto, x00011, x00110, x00210, x01010, x01011, x41211, x60210, x60211, x61011, x61210, x61211, x61311);
				}
	#pragma omp section
				{
					x61300_dot=rhsx61300(gto, x00100, x00200, x01000, x41300, x60300, x61100, x61200, x61300, x61310);
				}
	#pragma omp section
				{
					x61310_dot=rhsx61310(gto, x00011, x00110, x00210, x01010, x41310, x60310, x61110, x61210, x61300, x61310, x61311);
				}
	#pragma omp section
				{
					x61311_dot=rhsx61311(gto, x00011, x00110, x00210, x01010, x01011, x41311, x60310, x60311, x61111, x61211, x61310, x61311);
				}
	#pragma omp section
				{
					x62000_dot=rhsx62000(x00100, x00200, x02000, x60000, x62000, x62010, x62100, x62200);
				}
	#pragma omp section
				{
					x62010_dot=rhsx62010(x00011, x00110, x00210, x02010, x60010, x62000, x62010, x62011, x62110, x62210);
				}
	#pragma omp section
				{
					x62011_dot=rhsx62011(x00011, x00110, x00210, x02010, x02011, x60010, x60011, x62010, x62011, x62111, x62211);
				}
	#pragma omp section
				{
					x62100_dot=rhsx62100(x00100, x00200, x02000, x60100, x62000, x62100, x62110, x62300);
				}
	#pragma omp section
				{
					x62110_dot=rhsx62110(x00011, x00110, x00210, x02010, x60110, x62010, x62100, x62110, x62111, x62310);
				}
	#pragma omp section
				{
					x62111_dot=rhsx62111(x00011, x00110, x00210, x02010, x02011, x60110, x60111, x62011, x62110, x62111, x62311);
				}
	#pragma omp section
				{
					x62200_dot=rhsx62200(gto, x00100, x00200, x02000, x42200, x60200, x62000, x62200, x62210, x62300);
				}
	#pragma omp section
				{
					x62210_dot=rhsx62210(gto, x00011, x00110, x00210, x02010, x42210, x60210, x62010, x62200, x62210, x62211, x62310);
				}
	#pragma omp section
				{
					x62211_dot=rhsx62211(gto, x00011, x00110, x00210, x02010, x02011, x42211, x60210, x60211, x62011, x62210, x62211, x62311);
				}
	#pragma omp section
				{
					x62300_dot=rhsx62300(gto, x00100, x00200, x02000, x42300, x60300, x62100, x62200, x62300, x62310);
				}
	#pragma omp section
				{
					x62310_dot=rhsx62310(gto, x00011, x00110, x00210, x02010, x42310, x60310, x62110, x62210, x62300, x62310, x62311);
				}
	#pragma omp section
				{
					x62311_dot=rhsx62311(gto, x00011, x00110, x00210, x02010, x02011, x42311, x60310, x60311, x62111, x62211, x62310, x62311);
				}
			}
			#pragma omp sections
			{
				#pragma omp section
				{
					GR += step_size*GR_dot; /* Explicit Forward Euler method */
					G += step_size*G_dot; /* Explicit Forward Euler method */
					GrR += step_size*GrR_dot; /* Explicit Forward Euler method */
				}
				#pragma omp section
				{
					Gr += step_size*Gr_dot; /* Explicit Forward Euler method */
					GcR += step_size*GcR_dot; /* Explicit Forward Euler method */
					Gc += step_size*Gc_dot; /* Explicit Forward Euler method */
				}
				#pragma omp section
				{
					GBR += step_size*GBR_dot; /* Explicit Forward Euler method */
					GB += step_size*GB_dot; /* Explicit Forward Euler method */
					GBRb += step_size*GBRb_dot; /* Explicit Forward Euler method */
				}
				#pragma omp section
				{
					GBb += step_size*GBb_dot; /* Explicit Forward Euler method */
					MnPo += step_size*MnPo_dot; /* Explicit Forward Euler method */
					McPo += step_size*McPo_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					MnPt += step_size*MnPt_dot; /* Explicit Forward Euler method */
					McPt += step_size*McPt_dot; /* Explicit Forward Euler method */
					MnRt += step_size*MnRt_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					McRt += step_size*McRt_dot; /* Explicit Forward Euler method */
					MnRev += step_size*MnRev_dot; /* Explicit Forward Euler method */
					McRev += step_size*McRev_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					MnRo += step_size*MnRo_dot; /* Explicit Forward Euler method */
					McRo += step_size*McRo_dot; /* Explicit Forward Euler method */
					MnB += step_size*MnB_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					McB += step_size*McB_dot; /* Explicit Forward Euler method */
					MnNp += step_size*MnNp_dot; /* Explicit Forward Euler method */
					McNp += step_size*McNp_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					B += step_size*B_dot; /* Explicit Forward Euler method */
					Cl += step_size*Cl_dot; /* Explicit Forward Euler method */
					BC += step_size*BC_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					cyrev += step_size*cyrev_dot; /* Explicit Forward Euler method */
					revn += step_size*revn_dot; /* Explicit Forward Euler method */
					cyrevg += step_size*cyrevg_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					revng += step_size*revng_dot; /* Explicit Forward Euler method */
					cyrevgp += step_size*cyrevgp_dot; /* Explicit Forward Euler method */
					revngp += step_size*revngp_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					cyrevp += step_size*cyrevp_dot; /* Explicit Forward Euler method */
					revnp += step_size*revnp_dot; /* Explicit Forward Euler method */
					gto += step_size*gto_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x00001 += step_size*x00001_dot; /* Explicit Forward Euler method */
					x00011 += step_size*x00011_dot; /* Explicit Forward Euler method */
					x00100 += step_size*x00100_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x00110 += step_size*x00110_dot; /* Explicit Forward Euler method */
					x00200 += step_size*x00200_dot; /* Explicit Forward Euler method */
					x00210 += step_size*x00210_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x01000 += step_size*x01000_dot; /* Explicit Forward Euler method */
					x01010 += step_size*x01010_dot; /* Explicit Forward Euler method */
					x01011 += step_size*x01011_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x02000 += step_size*x02000_dot; /* Explicit Forward Euler method */
					x02010 += step_size*x02010_dot; /* Explicit Forward Euler method */
					x02011 += step_size*x02011_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x10000 += step_size*x10000_dot; /* Explicit Forward Euler method */
					x10100 += step_size*x10100_dot; /* Explicit Forward Euler method */
					x20000 += step_size*x20000_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x20010 += step_size*x20010_dot; /* Explicit Forward Euler method */
					x20011 += step_size*x20011_dot; /* Explicit Forward Euler method */
					x20100 += step_size*x20100_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x20110 += step_size*x20110_dot; /* Explicit Forward Euler method */
					x20111 += step_size*x20111_dot; /* Explicit Forward Euler method */
					x21000 += step_size*x21000_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x21010 += step_size*x21010_dot; /* Explicit Forward Euler method */
					x21011 += step_size*x21011_dot; /* Explicit Forward Euler method */
					x21100 += step_size*x21100_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x21110 += step_size*x21110_dot; /* Explicit Forward Euler method */
					x21111 += step_size*x21111_dot; /* Explicit Forward Euler method */
					x22000 += step_size*x22000_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x22010 += step_size*x22010_dot; /* Explicit Forward Euler method */
					x22011 += step_size*x22011_dot; /* Explicit Forward Euler method */
					x22100 += step_size*x22100_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x22110 += step_size*x22110_dot; /* Explicit Forward Euler method */
					x22111 += step_size*x22111_dot; /* Explicit Forward Euler method */
					x30000 += step_size*x30000_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x30100 += step_size*x30100_dot; /* Explicit Forward Euler method */
					x30200 += step_size*x30200_dot; /* Explicit Forward Euler method */
					x30300 += step_size*x30300_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x40000 += step_size*x40000_dot; /* Explicit Forward Euler method */
					x40010 += step_size*x40010_dot; /* Explicit Forward Euler method */
					x40011 += step_size*x40011_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x40100 += step_size*x40100_dot; /* Explicit Forward Euler method */
					x40110 += step_size*x40110_dot; /* Explicit Forward Euler method */
					x40111 += step_size*x40111_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x40200 += step_size*x40200_dot; /* Explicit Forward Euler method */
					x40210 += step_size*x40210_dot; /* Explicit Forward Euler method */
					x40211 += step_size*x40211_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x40300 += step_size*x40300_dot; /* Explicit Forward Euler method */
					x40310 += step_size*x40310_dot; /* Explicit Forward Euler method */
					x40311 += step_size*x40311_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x41000 += step_size*x41000_dot; /* Explicit Forward Euler method */
					x41010 += step_size*x41010_dot; /* Explicit Forward Euler method */
					x41011 += step_size*x41011_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x41100 += step_size*x41100_dot; /* Explicit Forward Euler method */
					x41110 += step_size*x41110_dot; /* Explicit Forward Euler method */
					x41111 += step_size*x41111_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x41200 += step_size*x41200_dot; /* Explicit Forward Euler method */
					x41210 += step_size*x41210_dot; /* Explicit Forward Euler method */
					x41211 += step_size*x41211_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x41300 += step_size*x41300_dot; /* Explicit Forward Euler method */
					x41310 += step_size*x41310_dot; /* Explicit Forward Euler method */
					x41311 += step_size*x41311_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x42000 += step_size*x42000_dot; /* Explicit Forward Euler method */
					x42010 += step_size*x42010_dot; /* Explicit Forward Euler method */
					x42011 += step_size*x42011_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x42100 += step_size*x42100_dot; /* Explicit Forward Euler method */
					x42110 += step_size*x42110_dot; /* Explicit Forward Euler method */
					x42111 += step_size*x42111_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x42200 += step_size*x42200_dot; /* Explicit Forward Euler method */
					x42210 += step_size*x42210_dot; /* Explicit Forward Euler method */
					x42211 += step_size*x42211_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x42300 += step_size*x42300_dot; /* Explicit Forward Euler method */
					x42310 += step_size*x42310_dot; /* Explicit Forward Euler method */
					x42311 += step_size*x42311_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x50000 += step_size*x50000_dot; /* Explicit Forward Euler method */
					x50010 += step_size*x50010_dot; /* Explicit Forward Euler method */
					x50011 += step_size*x50011_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x50100 += step_size*x50100_dot; /* Explicit Forward Euler method */
					x50110 += step_size*x50110_dot; /* Explicit Forward Euler method */
					x50111 += step_size*x50111_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x50200 += step_size*x50200_dot; /* Explicit Forward Euler method */
					x50210 += step_size*x50210_dot; /* Explicit Forward Euler method */
					x50211 += step_size*x50211_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x50300 += step_size*x50300_dot; /* Explicit Forward Euler method */
					x50310 += step_size*x50310_dot; /* Explicit Forward Euler method */
					x50311 += step_size*x50311_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x51000 += step_size*x51000_dot; /* Explicit Forward Euler method */
					x51010 += step_size*x51010_dot; /* Explicit Forward Euler method */
					x51011 += step_size*x51011_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x51100 += step_size*x51100_dot; /* Explicit Forward Euler method */
					x51110 += step_size*x51110_dot; /* Explicit Forward Euler method */
					x51111 += step_size*x51111_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x51200 += step_size*x51200_dot; /* Explicit Forward Euler method */
					x51210 += step_size*x51210_dot; /* Explicit Forward Euler method */
					x51211 += step_size*x51211_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x51300 += step_size*x51300_dot; /* Explicit Forward Euler method */
					x51310 += step_size*x51310_dot; /* Explicit Forward Euler method */
					x51311 += step_size*x51311_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x52000 += step_size*x52000_dot; /* Explicit Forward Euler method */
					x52010 += step_size*x52010_dot; /* Explicit Forward Euler method */
					x52011 += step_size*x52011_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x52100 += step_size*x52100_dot; /* Explicit Forward Euler method */
					x52110 += step_size*x52110_dot; /* Explicit Forward Euler method */
					x52111 += step_size*x52111_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x52200 += step_size*x52200_dot; /* Explicit Forward Euler method */
					x52210 += step_size*x52210_dot; /* Explicit Forward Euler method */
					x52211 += step_size*x52211_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x52300 += step_size*x52300_dot; /* Explicit Forward Euler method */
					x52310 += step_size*x52310_dot; /* Explicit Forward Euler method */
					x52311 += step_size*x52311_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x60000 += step_size*x60000_dot; /* Explicit Forward Euler method */
					x60010 += step_size*x60010_dot; /* Explicit Forward Euler method */
					x60011 += step_size*x60011_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x60100 += step_size*x60100_dot; /* Explicit Forward Euler method */
					x60110 += step_size*x60110_dot; /* Explicit Forward Euler method */
					x60111 += step_size*x60111_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x60200 += step_size*x60200_dot; /* Explicit Forward Euler method */
					x60210 += step_size*x60210_dot; /* Explicit Forward Euler method */
					x60211 += step_size*x60211_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x60300 += step_size*x60300_dot; /* Explicit Forward Euler method */
					x60310 += step_size*x60310_dot; /* Explicit Forward Euler method */
					x60311 += step_size*x60311_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x61000 += step_size*x61000_dot; /* Explicit Forward Euler method */
					x61010 += step_size*x61010_dot; /* Explicit Forward Euler method */
					x61011 += step_size*x61011_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x61100 += step_size*x61100_dot; /* Explicit Forward Euler method */
					x61110 += step_size*x61110_dot; /* Explicit Forward Euler method */
					x61111 += step_size*x61111_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x61200 += step_size*x61200_dot; /* Explicit Forward Euler method */
					x61210 += step_size*x61210_dot; /* Explicit Forward Euler method */
					x61211 += step_size*x61211_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x61300 += step_size*x61300_dot; /* Explicit Forward Euler method */
					x61310 += step_size*x61310_dot; /* Explicit Forward Euler method */
					x61311 += step_size*x61311_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x62000 += step_size*x62000_dot; /* Explicit Forward Euler method */
					x62010 += step_size*x62010_dot; /* Explicit Forward Euler method */
					x62011 += step_size*x62011_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x62100 += step_size*x62100_dot; /* Explicit Forward Euler method */
					x62110 += step_size*x62110_dot; /* Explicit Forward Euler method */
					x62111 += step_size*x62111_dot; /* Explicit Forward Euler method */
				}
#pragma omp section
				{
					x62200 += step_size*x62200_dot; /* Explicit Forward Euler method */
					x62210 += step_size*x62210_dot; /* Explicit Forward Euler method */
					x62211 += step_size*x62211_dot; /* Explicit Forward Euler method */
				}
				#pragma omp section
				{
					x62300 += step_size*x62300_dot; /* Explicit Forward Euler method */
					x62310 += step_size*x62310_dot; /* Explicit Forward Euler method */
					x62311 += step_size*x62311_dot; /* Explicit Forward Euler method */
				}
			}
		}

		if (t_step%record == 0) {
			fprintf(outfi1, "%f", t);
			fprintf(outfi1, "\t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f\n", GR, G, GrR, Gr, GcR, Gc, GBR, GB, GBRb, GBb, MnPo, McPo, MnPt, McPt, MnRt, McRt, MnRev, McRev, MnRo, McRo, MnB, McB, MnNp, McNp, B, Cl, BC, cyrev, revn, cyrevg, revng, cyrevgp, revngp, cyrevp, revnp, gto, x00001, x00011, x00100, x00110, x00200, x00210, x01000, x01010, x01011, x02000, x02010, x02011, x10000, x10100, x20000, x20010, x20011, x20100, x20110, x20111, x21000, x21010, x21011, x21100, x21110, x21111, x22000, x22010, x22011, x22100, x22110, x22111, x30000, x30100, x30200, x30300, x40000, x40010, x40011, x40100, x40110, x40111, x40200, x40210, x40211, x40300, x40310, x40311, x41000, x41010, x41011, x41100, x41110, x41111, x41200, x41210, x41211, x41300, x41310, x41311, x42000, x42010, x42011, x42100, x42110, x42111, x42200, x42210, x42211, x42300, x42310, x42311, x50000, x50010, x50011, x50100, x50110, x50111, x50200, x50210, x50211, x50300, x50310, x50311, x51000, x51010, x51011, x51100, x51110, x51111, x51200, x51210, x51211, x51300, x51310, x51311, x52000, x52010, x52011, x52100, x52110, x52111, x52200, x52210, x52211, x52300, x52310, x52311, x60000, x60010, x60011, x60100, x60110, x60111, x60200, x60210, x60211, x60300, x60310, x60311, x61000, x61010, x61011, x61100, x61110, x61111, x61200, x61210, x61211, x61300, x61310, x61311, x62000, x62010, x62011, x62100, x62110, x62111, x62200, x62210, x62211, x62300, x62310, x62311);
		}
	}
	
//////////////////////////////////////////////////////////////////////

    time_t end; 
    end=time(NULL); 
    printf("Runtime: %f\n",difftime(end,start)/60); 

    fprintf(outfi1,"\nRuntime: %f\n",difftime(end,start)/60); 
    fclose(outfi1); 
	return 1;
} 
