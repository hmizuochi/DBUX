This is a beta program of Database Unmixing with its addtional modes.
Developed by Mizuochi, H. (See Mizuochi et al., 2014; Mizuochi et al., 2017)

Sub-directories: src, input, output
Compile:
	$cd src
	$make
the src directory includes original programs and executable files. Edit parameters of "define.h".

input directory must include
	1) fine and coarse maps. each pixel format must be 2 bytes Integer.
		ex)
		spatial_2002001.bin,spatial_2002002.bin...
		temporal_2002001.bin,temporal_2002002.bin...

	2) map which indicates minimum coarse value for each pixel - "tmin.bin"
	3) namelist of fine maps - "spatial_binlist.txt"
	4) namelist of coarse maps - "temporal_binlist.txt"
	5) datelist in YYYYDOY format - "datelist.txt"

Execute:
        $./DBUX 0 #prediction via lookup table (LUT).

lookup maps, reliability maps, predicted maps will be written under output directory.

