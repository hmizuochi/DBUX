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
	printf("+-----DBUX ver 1.0----------------------------------------+\n");
	printf("| __/ Compile date : %s                      __/ | \n",__DATE__);
	printf("| __/         time : %s                         __/ | \n",__TIME__);
	printf("| __/ Program name : %s                            __/ | \n",__FILE__);
	printf("| __/ Made by : Mizuochi, H., University of Tsukuba. __/ | \n");
	printf("+---------------------------------------------------------+\n");

	/*variable definition and memory allocation*/
  int count=0;
	char pred_listname[]="../input/predlist.txt";
	char s_listname[]="../input/spatial_pairlist.txt";
	char t_listname[]="../input/temporal_pairlist.txt";
	short *t_pair[PAIRSIZE],*s_pair[PAIRSIZE];

	for(count=0;count<PAIRSIZE;count++){
		if((t_pair[count]=(short*)malloc(COL*ROW*sizeof(short)))==NULL||(s_pair[count]=(short*)malloc(COL*ROW*sizeof(short)))==NULL){
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
	if(ReadST(s_listname, t_listname, t_pair, s_pair)!=0){
		fprintf(stderr,"main: ReadST error!\n");
		exit(1);
	}
	printf("input data loaded! execute DBUX...\n");
	if(ExecDBUX(pred_listname, t_pair, s_pair)!=0){
		fprintf(stderr,"main: ExecDBUX error!\n");
		exit(1);
	}

	for(count=0;count<PAIRSIZE;count++){
		free(t_pair[count]);
		free(s_pair[count]);
	}
	return 0;
}
