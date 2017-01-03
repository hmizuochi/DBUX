#define NVALUE -20000
#define TPRANGE 10000
#define TNRANGE -10000
#define SPRANGE 10000
#define SNRANGE -10000
#define SAMPLESIZE 5113
#define COL 180
#define ROW 180
#define MAXTEXT 256
#define MWSIZE 3 //option for LUT mode ('0')
#define MAXLEVEL 40 //option for LUT mode ('0')
#define STEP 420 //option for LUT mode ('0')
int ReadST(char *s_listname, char *t_listname, short **t_input, short **s_input);
int ExecDBUX(char *date_listname, short **t_input, short **s_input);
int GenLUT(short *tmin_input, short **t_input, short **s_input, short **lookup_output, short **snum);
int AveLUT(short **lookup_input, short **snum_input, short **lookup_output, short **snum_output);
int PredDBUX(char *date_listname, short *tmin_input, short **t_input, short **s_input, short **lookup, short **snum);
