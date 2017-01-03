/*	compile: gcc 20150113_DU1.c -o 20150113_DU1
	usage: ./20150113_DU1 "LEVEL" "X(image columns)" "Y(image rows)" "MAXVALUE(Int16)" "MINVALUE(Int16)" "STEP(Int16)"
	ex) for gap-filled MODIS and landsat, LEVEL:0~30(MNDWI -0.50~1.00), use for loop in command line.
	confirm corrent directory has spatial-base&temporal-base binary files and lists of their filenames (spatial_binlist.txt & temporal_binlist.txt).
	confirm output binary file list (lookup_output_binlist.txt) and enough HDD capacity.
	image data format: Bytes=Int16, null=-20000. range=-10000~10000.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "define.h"

int main(int argc, char *argv[]){
	if(argc != 2){
		fprintf(stderr,"usage: ./DBUX MODE(0,1,2)\n");
		exit(1);
	}
	printf("+-----DBUX ver 1.0----------------------------------------+\n");
	printf("| __/ Compile date : %s                      __/ | \n",__DATE__);
	printf("| __/         time : %s                         __/ | \n",__TIME__);
	printf("| __/ Program name : %s                            __/ | \n",__FILE__);
	printf("| __/ Made by : Mizuochi, H., University of Tsukuba. __/ | \n");
	printf("+---------------------------------------------------------+\n");

	/*variable definition and memory allocation*/
  int count=0,MODE=0;
	char date_listname[]="../input/datelist.txt";
	char s_listname[]="../input/spatial_binlist.txt";
	char t_listname[]="../input/temporal_binlist.txt";
	short *t_input[SAMPLESIZE],*s_input[SAMPLESIZE];
	sscanf(argv[1],"%d",&MODE);

	for(count=0;count<SAMPLESIZE;count++){
		if((t_input[count]=(short*)malloc(COL*ROW*sizeof(short)))==NULL||(s_input[count]=(short*)malloc(COL*ROW*sizeof(short)))==NULL){
			fprintf(stderr,"main: can't allocate memory\n");
			exit(1);
		}
	}

	/*initial check for parameters*/
	if(TPRANGE<TNRANGE){
		fprintf(stderr,"main: TPRANGE must be greater than TNRANGE\n");
		exit(1);
	}
	if(SPRANGE<SNRANGE){
		fprintf(stderr,"main: SPRANGE must be greater than SNRANGE\n");
		exit(1);
	}
	if(((TNRANGE<=NVALUE)&&(NVALUE<=TPRANGE))||((SNRANGE<=NVALUE)&&(NVALUE<=SPRANGE))){
		fprintf(stderr,"main: NVALUE must NOT be between NRANGE and PRANGE\n");
		exit(1);
	}

	/*load input datasets*/
	printf("input data loading...\n");
	if(ReadST(s_listname, t_listname, t_input, s_input)!=0){
		fprintf(stderr,"main: ReadST error!\n");
		exit(1);
	}
	printf("input data loaded! execute DBUX...\n");
	if(MODE==0){
		if(ExecDBUX(date_listname, t_input, s_input)!=0){
			fprintf(stderr,"main: ExecDBUX error!\n");
			exit(1);
		}
	}
	for(count=0;count<SAMPLESIZE;count++){
		free(t_input[count]);
		free(s_input[count]);
	}
	return 0;
}
