/* simulate one cell using the equations from Diekman & Forger */
#include <stdio.h> 
#include <math.h> 
#include <stdlib.h>
//#include "nrutil.h"
#include "rhs.h"
#include <time.h> 
#include <string.h> 

int main(int argc, char *argv[])
{ 
	int nvar = 181; /* number of variables (+1 since statevector[0] is unused)

	/* Set initial conditions */
	float statevector[]={0, 0.787955, 0.210048, 0.558565, 0.416217, 0.804404, 0.193886, 0.39359, 0.599965, 0.404247, 0.589308, 
	4.18947, 0.579108, 8.34963, 1.41168, 12.3735, 4.47054, 31.869, 161.144, 3.76085, 1.09547, 
	35.1842, 8.56043, 0.520367, 0.25596, 23.7741, 48.0977, 6.5815, 123.95, 467.844, 0.193289, 
	58.2391, 0.0397137, 19.0662, 0.11012, 67.6175, 1.13837, 0.342061, 8.17977, 50.5794, 3.08486, 
	0.0985802, 1.30846, 41.1658, 0.0156084, 8.72089, 46.5721, 0.0031953, 1.36799, 0.422898, 1.09254, 
	0.00181722, .0000519698, 0.00512223, 0.106575, 0.000833073, 0.0566636, 0.00329174, .000014996, 0.00543134, 0.0876909, 
	0.000333744, 0.0548927, 0.00348546, .0000154312, 0.00510457, 0.0938758, 0.000353611, 0.0535918, 0.914365, 0.664039, 
	0.0226162, 0.0241531, 0.00356643, .0000194469, 0.000959363, 0.258442, 0.000985386, 0.0249991, 0.000157103, .000011242, 
	0.00132337, 0.0105914, 0.00023062, 0.025311, 0.00390694, .0000134882, 0.000809309, 0.124712, 0.000416813, 0.0148919, 
	0.000188833, .00000179979, 0.00109293, 0.00584769, .0000435837, 0.0166107, 0.00428454, .0000147716, 0.00086132, 0.137696, 
	0.000459811, 0.0160349, 0.0002068, .00000194018, 0.0011516, 0.00644626, .0000476349, 0.0177021, 0.00569806, 0.000190832, 
	0.0566834, 0.0287569, 0.000413832, 0.0494651, 0.00414505, 0.000386223, 0.124464, 0.00740796, 0.000654378, 0.104089, 
	0.00666974, .0000254221, 0.00484524, 0.0828541, 0.000294564, 0.0323573, 0.0012142, .000011508, 0.00976858, 0.0087474, 
	.0000798494, 0.0607878, 0.00706385, .0000264292, 0.00442251, 0.0852703, 0.000300906, 0.0301143, 0.00132641, .0000112497, 
	0.00878991, 0.0091505, .0000778605, 0.055586, 0.000792825, .0000257804, 0.00334582, 0.0274125, 0.000385498, 0.0589937, 
	0.000125512, .0000499769, 0.00735186, 0.00550233, 0.000631987, 0.12768, 0.00252101, .0000100779, 0.00248103, 0.0581233, 
	0.000207429, 0.0239137, 0.000305307, .00000490217, 0.00521949, 0.00753943, .0000639585, 0.0471435, 0.00259772, .0000101709, 
	0.00224807, 0.0602479, 0.000213449, 0.0225783, 0.000320067, .00000460372, 0.004701, 0.00797183, .0000635261, 0.0440794};
	float output[nvar];
	time_t start; 
	start=time(NULL); 
	
	float t1=0; /* initial time */
	float t2=500; /* end time */
	float t=t1;
	float step_size=0.0001; /* step size (in hours) */
	int nstep=(t2-t1)/step_size; /* number of identical steps to take */

	int i,j,t_step;
	float record_size=0.1; /* record data to file every record_size hours (if step_size <= record_size)*/
	int record; /* record data to file every record time steps */
	if (step_size<=record_size) record = round(record_size/step_size);
	else record = 1;

	float *results[nvar][nstep/record];
		
	/* statevector[0...nvar] is the vector of all state variables (1 through nvar) at the beginning of a timestep */
	/* output[0...nvar] is the vector of all state variables (1 through nvar) at the end of a timestep */
	/* rhs is the user supplied routine that computes right hand derivatives using statevector, and returns output */
	

//////////////////////// For Data Output /////////////////////////////
	char filename1[50]; 
	FILE *outfi1; 

	sprintf(filename1,"./Results/Kim.txt"); /* records output every record timesteps */
	remove(filename1); 
	if ((outfi1 = fopen(filename1, "w")) == NULL) { 
	  printf("can not open output file"); 
	};/* opens an output file*/ 

	/* Output initial conditions */
	fprintf(outfi1, "%f", t);
	for (i=1; i<=nvar; i++) {
	  fprintf(outfi1, "\t%f", statevector[i]);
	}
	fprintf(outfi1, "\n");
	
/////////////////////// Run computation //////////////////////////////
	for (t_step=1; t_step<=nstep; t_step++) {
	  t=t_step*step_size;

	  /* Calculate derivatives */
	  rhs(t,statevector, output);
	  
	  for (i=1; i<=nvar; i++) {
	    statevector[i] += step_size*output[i]; /* Explicit Forward Euler method */
	  }
	  if (t_step%record == 0) {
	    fprintf(outfi1, "%f", t);
	    for (i=1; i<=nvar; i++) {
	    fprintf(outfi1, "\t%f", statevector[i]);
	    }
	    fprintf(outfi1, "\n");
	  }
	}
    fclose(outfi1); 
//////////////////////////////////////////////////////////////////////

    time_t end; 
    end=time(NULL); 
    printf("Runtime: %f\n",difftime(end,start)/60); 

	return 1;
} 