This is a stable program of Database Unmixing, which is LUT-based data fusion algorithm.
Developed by Mizuochi, H.

References:

Mizuochi, H., Hiyama, T., Ohta, T., Nasahara, K. N. (2014): Evaluation of the surface water distribution in north-central Namibia based on MODIS and AMSR series. Remote Sens. 6(8), pp. 7660-7682.

Mizuochi, H., Hiyama, T., Ohta, T., Fujioka, Y., Kambatuku, J. R., Iijima, M., Nasahara, K. N. (2017): Development and evaluation of lookup-table-based approach of data fusion for seasonal wetlands monitoring: Integrated use of AMSR series, MODIS and Landsat. Remote Sens. Environ. Submitted.

Extract: $ tar zxvf DBUX.tar.gz

Sub-directories: src, input, output, sample_data

Before running program, the following preparation is required.

A) Edit parameters in "src/define.h".

B) Compile program.

	$cd src

	$make #require gcc

C) put the following input maps and text files under input directory.

	1) temporally frequent maps ("T maps") and spatially fine maps ("S maps").
		 data format is 2 bytes Integer. filename format is:
		 "spatial_YYYYDOY.raw", "temporal_YYYYDOY.raw"

	2) map which indicates minimum value of T maps for each pixel - "tmin.raw"

	3) namelist of S maps which are used for LUT generation - "spatial_pairlist.txt"
			e.g.
			spatial_2002001.raw
			spatial_2002002.raw
			spatial_2002005.raw
			spatial_2002007.raw
			...

	4) namelist of T maps which are used for LUT generation - "temporal_pairlist.txt".
		 dates and order must be the same as spatial_pairlist.txt.
			e.g.
			temporal_2002001.raw
			temporal_2002002.raw
			temporal_2002005.raw
			temporal_2002007.raw
			...

	5) list of prediction dates in YYYYDOY format - "predlist.txt"
			e.g.
			2002001
			2002002
			2002003
			2002004
			...

Execute:
        $./DBUX

lookup maps, reliability maps, predicted maps will be written under output directory.

predicted maps: "spatial_YYYYDOY_comp.raw", reliability maps: "spatial_YYYYDOY_rel.raw"

sample_data includes 9 Smaps and 366 Tmaps (i.e. 9 match-up pairs are available).
Smap is modified normalized water index (MNDWI) derived from Landsat TM and ETM+, and Tmap is that derived from MODIS, with scalefactor = 10000.
For accurate prediction, more match-up pairs are desirable.
