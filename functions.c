#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "define.h"

int ReadST(char *s_listname, char *t_listname, short **t_input, short **s_input){
  FILE *fp0,*fp1,*fp2;
  char buf[MAXTEXT],s_filename[MAXTEXT],t_filename[MAXTEXT];
  int count=0;
  if((fp0=fopen(s_listname,"r"))==NULL){
    fprintf(stderr,"ReadST: %s was not opened successfully\n",s_listname);
    exit(1);
  }
  if((fp1=fopen(t_listname,"r"))==NULL){
    fprintf(stderr,"ReadST: %s was not opened successfully\n",t_listname);
    exit(1);
  }
  while(fgets(buf,sizeof(buf),fp0) != NULL){
    strtok(buf,"\n\0");
    snprintf(s_filename,MAXTEXT,"../input/%s",buf);
  	fgets(buf,sizeof(buf),fp1);
  	strtok(buf,"\n\0");
    snprintf(t_filename,MAXTEXT,"../input/%s",buf);
    printf("input s_filename is %s\n",s_filename);
  	printf("input t_filename is %s\n",t_filename);
  	if((fp2=fopen(s_filename,"rb"))==NULL){
  		fprintf(stderr,"ReadST: %s was not opened successfully\n",s_filename);
      exit(1);
  	}
  	if((fread(s_input[count],sizeof(short),COL*ROW,fp2)) != (unsigned int)(COL*ROW)) {
  		fprintf(stderr, "ReadST: can't read s file\n");
      exit(1);
  	}
  	fclose(fp2);
  	if((fp2=fopen(t_filename,"rb"))==NULL){
  		fprintf(stderr,"ReadST: %s was not opened successfully\n",t_filename);
      exit(1);
  	}
  	if((fread(t_input[count],sizeof(short),COL*ROW,fp2)) != (unsigned int)(COL*ROW)) {
  		fprintf(stderr, "ReadST: can't read t file\n");
      exit(1);
  	}
  	fclose(fp2);
  	count++;
  }
  if(count!=PAIRSIZE){
  	fprintf(stderr,"ReadST Warning: PAIRSIZE and the number of read files were not consistent\n");
  }
  return 0;
}

int ExecDBUX(char *date_listname, short **t_input, short **s_input){
  int i=0,LEVEL=0;
  short *tmin_input, *lookup[MAXLEVEL+1], *lookup_ave[MAXLEVEL+1], *snum[MAXLEVEL+1], *snum_ave[MAXLEVEL+1];
  FILE *fp;
  for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
		if((lookup[LEVEL]=(short*)malloc(COL*ROW*sizeof(short)))==NULL||(lookup_ave[LEVEL]=(short*)malloc(COL*ROW*sizeof(short)))==NULL||
    (snum[LEVEL]=(short*)malloc(COL*ROW*sizeof(short)))==NULL||(snum_ave[LEVEL]=(short*)malloc(COL*ROW*sizeof(short)))==NULL){
			fprintf(stderr,"ExecDBUX: can't allocate memory\n");
      exit(1);
		}
    for(i=0;i<COL*ROW;i++){
      lookup[LEVEL][i]=0;
      lookup_ave[LEVEL][i]=0;
			snum[LEVEL][i]=0;
      snum_ave[LEVEL][i]=0;
    }
	}
  if((tmin_input=(short*)malloc(COL*ROW*sizeof(short)))==NULL){
    fprintf(stderr,"ExecDBUX: can't allocate memory\n");
    exit(1);
  }
  if((fp=fopen("../input/tmin.raw","rb"))==NULL){
		fprintf(stderr,"ExecDBUX: %s was not opened successfully\n","tmin.raw");
    exit(1);
	}
	if((fread(tmin_input,sizeof(short),COL*ROW,fp)) != (unsigned int)(COL*ROW)){
		fprintf(stderr,"ExecDBUX: can't read tmin\n");
    exit(1);
	}
	fclose(fp);
  printf("LUT generating...\n");
  if(GenLUT(tmin_input, t_input, s_input, lookup, snum)!=0){
    fprintf(stderr,"ExecDBUX: GenLUT error!\n");
    exit(1);
  } //input:tmin_input,t_input,s_input    output:lookup,snum    generate lookup maps.
  printf("LUT generated! applying moving average with size %d...\n",MWSIZE);
  if(AveLUT(lookup, snum, lookup_ave, snum_ave)!=0){
    fprintf(stderr,"ExecDBUX: AveLUT error!\n");
    exit(1);
  } //input:lookup,snum   output:lookup_ave,snum_ave    generate averaged lookup maps.
  printf("do prediction...\n");
  if(PredDBUX(date_listname, tmin_input, lookup_ave, snum_ave)!=0){
    fprintf(stderr,"ExecDBUX: PredDBUX error!\n");
    exit(1);
  } //input:date_listname,tmin_input,lookup_ave,snum_ave    output:None,    generate predected maps.
  for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
    free(lookup[LEVEL]);
    free(lookup_ave[LEVEL]);
    free(snum[LEVEL]);
    free(snum_ave[LEVEL]);
  }
  free(tmin_input);
  return 0;
}

int GenLUT(short *tmin_input, short **t_input, short **s_input, short **lookup_output, short **snum){
  short *s_image;
  int *lookup[MAXLEVEL+1];
  int LEVEL=0,count=0,i=0;
  char lookup_output_filename[256];
  FILE *fp;
  for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
		if((lookup[LEVEL]=(int*)malloc(COL*ROW*sizeof(int)))==NULL){
			fprintf(stderr,"GenLUT: can't allocate memory\n");
      exit(1);
		}
	}
  if((s_image=(short*)malloc(COL*ROW*sizeof(short)))==NULL){
    fprintf(stderr,"GenLUT: can't allocate memory\n");
    exit(1);
  }
  for(i=0;i<COL*ROW;i++){
    for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
    	lookup[LEVEL][i]=0;
	  }
  }
  for(count=0;count<PAIRSIZE;count++){
    for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
      for(i=0;i<COL*ROW;i++){
        s_image[i]=NVALUE;
        if(LEVEL==0){
          if((TNRANGE<=t_input[count][i])&&(t_input[count][i]<=tmin_input[i]+(LEVEL+1)*STEP)){
            s_image[i]=s_input[count][i];
          }else{
            s_image[i]=NVALUE;
          }
        }else if(LEVEL==MAXLEVEL){
          if((tmin_input[i]+LEVEL*STEP<t_input[count][i])&&(t_input[count][i]<=TPRANGE)){
            s_image[i]=s_input[count][i];
          }else{
            s_image[i]=NVALUE;
          }
        }else{
          if(tmin_input[i]+LEVEL*STEP<t_input[count][i]&&t_input[count][i]<=tmin_input[i]+(LEVEL+1)*STEP){
            s_image[i]=s_input[count][i];
          }else{
            s_image[i]=NVALUE;
          }
        }
        if(SNRANGE<=s_image[i]&&s_image[i]<=SPRANGE){
          snum[LEVEL][i]+=1;
        }else{
          s_image[i]=0;
        }
        lookup[LEVEL][i]=lookup[LEVEL][i]+s_image[i];
      }
    }
  }
  for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
    for(i=0;i<COL*ROW;i++){
      if(snum[LEVEL][i]==0){
        lookup[LEVEL][i]=NVALUE;
      }else{
        lookup[LEVEL][i]=lookup[LEVEL][i]/snum[LEVEL][i];
      }
    }
  }
  /*output lookup maps.*/
  for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
    snprintf(lookup_output_filename,MAXTEXT,"../output/LOOKUP%d.raw",LEVEL);
    for(i=0;i<COL*ROW;i++){
      lookup_output[LEVEL][i]=(short)lookup[LEVEL][i];
    }
    if((fp=fopen(lookup_output_filename,"wb"))==NULL){
      fprintf(stderr,"GenLUT: can't open output file\n");
      exit(1);
    }
    if((fwrite(lookup_output[LEVEL],sizeof(short),COL*ROW,fp)) != (unsigned int)(COL*ROW)){
      fprintf(stderr,"GenLUT: can't write output file\n");
      exit(1);
    }
    fclose(fp);
    snprintf(lookup_output_filename,MAXTEXT,"../output/LOOKUP%d_rel.bin",LEVEL);
    if((fp=fopen(lookup_output_filename,"wb"))==NULL){
      fprintf(stderr,"GenLUT: can't open output file\n");
      exit(1);
    }
    if((fwrite(snum[LEVEL],sizeof(short),COL*ROW,fp)) != (unsigned int)(COL*ROW)){
      fprintf(stderr,"GenLUT: can't write output file\n");
      exit(1);
    }
    fclose(fp);
  }
  for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
    free(lookup[LEVEL]);
  }
  free(s_image);
  return 0;
}

int AveLUT(short **lookup_input, short **snum_input, short **lookup_output, short **snum_output){
  char output_lookup_filename[256],output_snum_filename[256];
  int *lookup[MAXLEVEL+1];
  int LEVEL=0,i=0,m=0,count=0;
  FILE *fp;
  for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
		if((lookup[LEVEL]=(int*)malloc(COL*ROW*sizeof(int)))==NULL){
			fprintf(stderr,"GenLUT: can't allocate memory\n");
      exit(1);
		}
	}
	for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
    for(i=0;i<ROW*COL;i++){
		  lookup_output[LEVEL][i]=0;
      lookup[LEVEL][i]=0;
      snum_output[LEVEL][i]=0;
		}
	}

	for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
	   for(i=0;i<ROW*COL;i++){
	      if((int)MWSIZE/2<=LEVEL&&LEVEL<=MAXLEVEL-(int)MWSIZE/2){
          count=0;
          for(m=0;m<MWSIZE;m++){
            if(lookup_input[LEVEL+(int)(m-MWSIZE/2)][i]<=SNRANGE||SPRANGE<=lookup_input[LEVEL+(int)(m-MWSIZE/2)][i]){
              count+=1;
            }else{
				      lookup[LEVEL][i]+=lookup_input[LEVEL+(int)(m-MWSIZE/2)][i];
              snum_output[LEVEL][i]+=snum_input[LEVEL+(int)(m-MWSIZE/2)][i];
				    }
          }
          if(count==MWSIZE){
            lookup[LEVEL][i]=NVALUE;
					  snum_output[LEVEL][i]=0;
          }else{
            lookup[LEVEL][i]=lookup[LEVEL][i]/(MWSIZE-count);
          }
        }else{
				  lookup[LEVEL][i]=lookup_input[LEVEL][i];
				  snum_output[LEVEL][i]=snum_input[LEVEL][i];
		    }
		}
  }
  for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
	  snprintf(output_lookup_filename,MAXTEXT,"../output/LOOKUP%d_ave.raw",LEVEL);
    for(i=0;i<COL*ROW;i++){
      lookup_output[LEVEL][i]=(short)lookup[LEVEL][i];
    }
	  if((fp=fopen(output_lookup_filename,"wb"))==NULL){
		  fprintf(stderr,"AveLUT: can't open output file\n");
      exit(1);
	  }
	  if((fwrite(lookup_output[LEVEL],sizeof(short),COL*ROW,fp)) != (unsigned int)(COL*ROW)){
		  fprintf(stderr,"AveLUT: can't write output file\n");
      exit(1);
	  }
	  fclose(fp);
	  snprintf(output_snum_filename,MAXTEXT,"../output/LOOKUP%d_ave_rel.raw",LEVEL);
    if((fp=fopen(output_snum_filename,"wb"))==NULL){
      fprintf(stderr,"AveLUT: can't open output file\n");
      exit(1);
    }
    if((fwrite(snum_output[LEVEL],sizeof(short),COL*ROW,fp)) != (unsigned int)(COL*ROW)){
      fprintf(stderr,"AveLUT: can't write output file\n");
      exit(1);
    }
    fclose(fp);
  }
  return 0;
}

int PredDBUX(char *date_listname, short *tmin_input, short **lookup, short **snum){
	FILE *fp,*fp2;
	int i=0,count=0,LEVEL=0;
	char date[MAXTEXT],s_filename[MAXTEXT],t_filename[MAXTEXT],output_filename[MAXTEXT];
	short *t_input, *s_input, *s_output;
	int *image,*snum_output;

  if(((t_input=(short*)malloc(COL*ROW*sizeof(short)))==NULL)||
  ((s_input=(short*)malloc(COL*ROW*sizeof(short)))==NULL)||((s_output=(short*)malloc(COL*ROW*sizeof(short)))==NULL)||
  ((image=(int*)malloc(COL*ROW*sizeof(int)))==NULL)||((snum_output=(int*)malloc(COL*ROW*sizeof(int)))==NULL)){
    fprintf(stderr,"PredDBUX: can't allocate memory\n");
    exit(1);
  }
  //read prediction date.
	if((fp=fopen(date_listname,"r"))==NULL){
    fprintf(stderr,"PredDBUX: %s was not opened successfully\n",date_listname);
    exit(1);
  }
  //loop for each prediction date.
  for(count=0;count<PREDSIZE;count++){
	  if(fgets(date,sizeof(date),fp)==NULL){
      fprintf(stderr,"PredDBUX Warning: datelist EOF\n");
    }
		strtok(date,"\n\0");
		printf("calculation date is %s\n",date);
    snprintf(s_filename,MAXTEXT,"../input/spatial_%s.raw",date);
    snprintf(t_filename,MAXTEXT,"../input/temporal_%s.raw",date);
    if((fp2=fopen(s_filename,"rb"))==NULL){
      for(i=0;i<COL*ROW;i++)
        s_input[i]=NVALUE;
    }else{
      if((fread(s_input,sizeof(short),COL*ROW,fp2)) != (unsigned int)(COL*ROW)){
        fprintf(stderr,"PredDBUX: can't read s file\n");
        exit(1);
      }
      fclose(fp2);
    }
    if((fp2=fopen(t_filename,"rb"))==NULL){
      fprintf(stderr,"PredDBUX Warning: can't open %s, continue with nullmap\n",t_filename);
      for(i=0;i<COL*ROW;i++)
        t_input[i]=NVALUE;
    }else{
      if((fread(t_input,sizeof(short),COL*ROW,fp2)) != (unsigned int)(COL*ROW)){
        fprintf(stderr,"PredDBUX: can't read t file\n");
        exit(1);
      }
      fclose(fp2);
    }
		/* image initialize */
		for(i=0;i<COL*ROW;i++){
			image[i]=NVALUE;
      snum_output[i]=NVALUE;
      s_output[i]=NVALUE;
		}
		/* make gap-filled map. */
		for(i=0;i<COL*ROW;i++){
			if((t_input[i]<TNRANGE)||(t_input[i]>TPRANGE)){
				image[i]=s_input[i]; //no temporal-base map.
				snum_output[i]=-88;
			}else{
				if((s_input[i]<SNRANGE)||(s_input[i]>SPRANGE)){
					//exist temporal-base map, no spatial-base map -> gap-filling
					for(LEVEL=0;LEVEL<=MAXLEVEL;LEVEL++){
        		if(LEVEL==0){
						  if((TNRANGE<=t_input[i])&&(t_input[i]<=tmin_input[i]+(LEVEL+1)*STEP)){
								image[i]=lookup[LEVEL][i];
								snum_output[i]=snum[LEVEL][i];
							}
						}else if(LEVEL==MAXLEVEL){
							if((tmin_input[i]+STEP*LEVEL<t_input[i])&&(t_input[i]<=TPRANGE)){
								image[i]=lookup[LEVEL][i];
								snum_output[i]=snum[LEVEL][i];
							}
						}else{
							if((tmin_input[i]+STEP*LEVEL<t_input[i])&&(t_input[i]<=tmin_input[i]+STEP*(LEVEL+1))){
								image[i]=lookup[LEVEL][i];
								snum_output[i]=snum[LEVEL][i];
							}
						}
					}
				}else{
					//exist temporal-base map, exist spatial-base map -> spatial-base map
					image[i]=s_input[i];
					snum_output[i]=-99;
				}
			}
		}

		//output gap-filled map.
		snprintf(output_filename,MAXTEXT,"../output/spatial_%s_comp.raw",date);
		for(i=0;i<COL*ROW;i++){
      s_output[i]=(short)image[i];
    }
    if((fp2=fopen(output_filename,"wb"))==NULL){
      fprintf(stderr,"PredDBUX: can't open output file\n");
      exit(1);
    }
    if((fwrite(s_output,sizeof(short),COL*ROW,fp2)) != (unsigned int)(COL*ROW)){
      fprintf(stderr,"PredDBUX: can't write output file\n");
      exit(1);
    }else{
			printf("%s generated!\n",output_filename);
		}
    fclose(fp2);
    snprintf(output_filename,MAXTEXT,"../output/spatial_%s_rel.raw",date);
    for(i=0;i<COL*ROW;i++){
      s_output[i]=(short)snum_output[i];
    }
    if((fp2=fopen(output_filename,"wb"))==NULL){
      fprintf(stderr,"PredDBUX: can't open output file\n");
      exit(1);
    }
    if((fwrite(s_output,sizeof(short),COL*ROW,fp2)) != (unsigned int)(COL*ROW)){
      fprintf(stderr,"PredDBUX: can't write output file\n");
      exit(1);
    }else{
      printf("%s generated!\n",output_filename);
    }
    fclose(fp2);
	}
  fclose(fp);
  free(snum_output);
  free(image);
  free(t_input);
  free(s_input);
  free(s_output);
  return 0;
}
