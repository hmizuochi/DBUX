#define NVALUE -20000 //null value. must not be between TNRANGE and TPRANGE, nor SNRANGE and SPRANGE.
#define TPRANGE 10000 //potentially maximum value of temporally frequent maps ("T maps").
#define TNRANGE -10000 //potentially minimum value of temporally frequent maps.
#define SPRANGE 10000 //potentially maximum value of spatially fine maps ("S maps").
#define SNRANGE -10000 //potentially minumum value of spatially fine maps.
#define PAIRSIZE 36 //the number of match-up pairs.
#define PREDSIZE 1 //the number of dates in which DBUX will make prediction.
#define COL 180 //columns of a map. must be common between T maps and S maps.
#define ROW 180 //rows of a map must be common between T maps and S maps.
#define MAXTEXT 256 //maximum text buffer
#define MWSIZE 3 //moving window size for LUT. if not apply moving average for LUT, set 1.
#define STEP 420 //slicing step of LUT.
#define LIMITLEVEL 40 //maximum level of LUT. DBUX calulates it automatically when '0' was set.
int ReadST(char *s_listname, char *t_listname, short **t_input, short **s_input);
int ExecDBUX(char *date_listname, short **t_input, short **s_input);
int GenLUT(short *tmin_input, short *tmax_input, short **t_input, short **s_input, short **lookup_output, short **snum, int MAXLEVEL);
int AveLUT(short **lookup_input, short **snum_input, short **lookup_output, short **snum_output, int MAXLEVEL);
int PredDBUX(char *date_listname, short *tmin_input, short **lookup, short **snum, int MAXLEVEL);
int VisualizeLUT(short *tmin_input, short **lookup, int MAXLEVEL);
int StatsCalc(short **t_input, short *tmin_input, short *tmax_input);
