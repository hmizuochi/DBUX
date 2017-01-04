This is a stable program of Database Unmixing, which is LUT-based data fusion algorithm.
Developed by Mizuochi, H. (See Mizuochi et al., 2014; Mizuochi et al., 2017)

Sub-directories: src, input, output

Compile:

	$cd src

	$make #require gcc

Before running program, the following preparation is required.

A) Edit parameters in "define.h".

B) put the following input maps and text files under input directory.

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

predicted maps: "spatial_YYYYDOY_comp.raw"
reliability maps: "spatial_YYYYDOY_rel.raw"
