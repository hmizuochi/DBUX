Description
=======================
This is a stable program of Database Unmixing (Mizuochi et al., 2014), which is LUT-based data fusion algorithm, written in C and python languages. They have been tested in Linux system (for C, using gcc compiler).

### Revision History:
First release on 2017/03/01: original version (Mizuochi et al., 2014)
revised on 2018/02/01: LUT gap-filling implemented  
revised on 2018/05/01: add python version, which includes spatial-smoothing and uncertainty estimation


### License:
This program is provided free of charge, without restriction of use. For the full license information, see "LICENSE.txt". Publications, models and data products that make use of this program must include proper acknowledgement, including citing the journal article as in the following references.

### References:
1. Mizuochi, H., Hiyama, T., Ohta, T., Nasahara, K. N. (2014): Evaluation of the surface water distribution in north-central Namibia based on MODIS and AMSR series. Remote Sensing. 6(8), pp. 7660-7682.
2. Mizuochi, H., Hiyama, T., Ohta, T., Fujioka, Y., Kambatuku, J. R., Iijima, M., Nasahara, K. N. (2017): Development and evaluation of a lookup-table-based approach to data fusion for seasonal wetlands monitoring: An integrated use of AMSR series, MODIS, and Landsat. Remote Sensing of Environment. 199c, pp. 370-388.


How to Use
=====================
### Extract:
	$ tar zxvf DBUX.tar.gz  
	Sub-directories: src, input, output, sample_data  

### Preparation:
Before running program, the following preparation is required.  
##### A) Put the following input maps and text files under the input directory.

        - temporally frequent maps ("T maps" hereafter) and spatially fine maps ("S maps" hereafter).
                 data format is 2 bytes Integer. filename must be:
                 "spatial_YYYYDOY.raw", "temporal_YYYYDOY.raw"; YYYYDOY is indicated by the following text files.

        - namelist of S maps which are used for LUT generation - "spatial_pairlist.txt"
                        e.g.
                        spatial_2002001.raw
                        spatial_2002002.raw
                        spatial_2002005.raw
                        spatial_2002007.raw
                        ...

        - namelist of T maps which are used for LUT generation - "temporal_pairlist.txt".
                 dates and order must be the same as spatial_pairlist.txt.
                        e.g.
                        temporal_2002001.raw
                        temporal_2002002.raw
                        temporal_2002005.raw
                        temporal_2002007.raw
                        ...

        - list of prediction dates - "predlist.txt"
                        e.g.
                        2002001
                        2002002
                        2002003
                        ...


##### B) Give parameters.
*C program*
Edit the following parameters in "src/define.h", and compile the program:

        $cd src
        $make clean
        $make #require gcc
please do this every time when you revise "src/define.h".

*python program*
Give the parameters as arguments when you run the program.

	NVALUE //null value. must NOT be between TNRANGE and TPRANGE, nor SNRANGE and SPRANGE.  
	TPRANGE //potentially maximum value of temporally frequent maps ("T maps").  
	TNRANGE //potentially minimum value of temporally frequent maps.  
	SPRANGE //potentially maximum value of spatially fine maps ("S maps").
	SNRANGE //potentially minumum value of spatially fine maps.
	PAIRSIZE //the number of match-up pairs.
	PREDSIZE //the number of dates in which DBUX will make prediction.
	COL //columns of a map. it must be common between T maps and S maps.
	ROW //rows of a map. it must be common between T maps and S maps.
	LUT_MWSIZE //moving window size for LUT (it must be an odd number). if not apply moving average for LUT, set 1.
	STEP //slicing step of LUT. please determine it by uncertainty calculation (See Mizuochi et al., 2017), or by trial and error.
	MAP_MWSIZE //moving window size for spatial smoothing (it must be an odd number). it can be used only with python program.


### Execute:
*C program*

        $./DBUX.exe

*python program*

	$python 20180424_DBUX.py NVALUE TPRANGE TNRANGE SPRANGE SNRANGE PAIRSIZE PREDSIZE COL ROW LUT_MWSIZE STEP MAP_MWSIZE

lookup maps, reliability maps (See Mizuochi et al., 2017), predicted maps will be generated under output directory. if you use python program, uncertainty maps, which indicate uncertainty of the prediction for each pixel, are also created.  
predicted maps: "spatial_YYYYDOY_comp.raw", reliability maps: "spatial_YYYYDOY_rel.raw", uncertainty maps: "spatial_YYYYDOY_uncertainty.raw"


Sample data
===========

sample_data directory includes 9 Smaps and 366 Tmaps (i.e. 9 match-up pairs are available, and 366 prediction can be executable) over seasonal wetlands in north-central Namibia. (upper-left 17:27S,15:21E; lower-right 17:30S,15:24E; resolution 1 arcsec)
Smap is modified normalized water index (MNDWI) derived from Landsat TM and ETM+, and Tmap is that derived from MODIS, with scalefactor = 10000.
For accurate prediction, more match-up pairs are desirable.
