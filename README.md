Description (ver. 1.0)
======================
This is a stable program of Database Unmixing (Mizuochi et al., 2014), which is LUT-based data fusion algorithm, written in C language. It has been tested in Linux system using gcc compiler.

### ver 1.0: original version (first release)

### License:
This program is provided free of charge, without restriction of use. For the full license information, see "LICENSE.txt". Publications, models and data products that make use of this program must include proper acknowledgement, including citing the journal article as in the following references.

### References:
1. Mizuochi, H., Hiyama, T., Ohta, T., Nasahara, K. N. (2014): Evaluation of the surface water distribution in north-central Namibia based on MODIS and AMSR series. Remote Sensing. 6(8), pp. 7660-7682.
2. Mizuochi, H., Hiyama, T., Ohta, T., Fujioka, Y., Kambatuku, J. R., Iijima, M., Nasahara, K. N. (2017): Development and evaluation of a lookup-table-based approach to data fusion for seasonal wetlands monitoring: An integrated use of AMSR series, MODIS, and Landsat. Remote Sensing of Environment. 199c, pp. 370-388.


Usage of this program
=====================
### Extract:
	$ tar zxvf DBUX.tar.gz  
	Sub-directories: src, input, output, sample_data  

### Preparation:
Before running program, the following preparation is required.  
A) Edit parameters in "src/define.h".  
B) Compile program.

	$cd src  
	$make #require gcc  
when you revise define.h, please use

	$make clean
	$make

C) Put the following input maps and text files under the input directory.  

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

### Execute:

        $./DBUX.exe

lookup maps, reliability maps (See Mizuochi et al., 2017), predicted maps will be generated under output directory.  
predicted maps: "spatial_YYYYDOY_comp.raw", reliability maps: "spatial_YYYYDOY_rel.raw"


Sample data
===========

sample_data directory includes 9 Smaps and 366 Tmaps (i.e. 9 match-up pairs are available, and 366 prediction can be executable) over seasonal wetlands in north-central Namibia.
Smap is modified normalized water index (MNDWI) derived from Landsat TM and ETM+, and Tmap is that derived from MODIS, with scalefactor = 10000.
For accurate prediction, more match-up pairs are desirable.
