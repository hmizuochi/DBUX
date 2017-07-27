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
        int count=0,i=0;
	char pred_listname[]="../input/predlist.txt";
	char s_listname[]="../input/spatial_pairlist.txt";
	char t_listname[]="../input/temporal_pairlist.txt";
	//short *t_pair[PAIRSIZE],*s_pair[PAIRSIZE];
	short **t_pair,**s_pair;
	t_pair=malloc(sizeof(short *)*PAIRSIZE);
	s_pair=malloc(sizeof(short *)*PAIRSIZE);
	for(count=0;count<PAIRSIZE;count++){
		if((t_pair[count]=malloc(COL*ROW*sizeof(short)))==NULL||(s_pair[count]=malloc(COL*ROW*sizeof(short)))==NULL){
			fprintf(stderr,"main: can't allocate memory\n");
			exit(1);
		}
	}
	//*initialization
	for(count=0;count<PAIRSIZE;count++){
		for(i=0;i<ROW*COL;i++){
			t_pair[count][i]=0;
			s_pair[count][i]=0;
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
	/*execute DBUX*/
	printf("input data loaded! execute DBUX...\n");
	if(ExecDBUX(pred_listname, t_pair, s_pair)!=0){
		fprintf(stderr,"main: ExecDBUX error!\n");
		exit(1);
	}
	/*memory release*/
	for(count=0;count<PAIRSIZE;count++){
		free(t_pair[count]);
		free(s_pair[count]);
	}
	return 0;
}
