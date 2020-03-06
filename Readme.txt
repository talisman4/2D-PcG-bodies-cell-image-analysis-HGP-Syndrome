The algorithm has been implemented and tested with MATLAB versions 8.4 (2014a) and 9.1 (2016b) on Linux.
1) For the segmentation of cell nuclei, we use an algorithm that implements one of the most efficient minimization methods for a non convex segmentation problem.
To this purpose we use a C code to be turned in MEX file so that it can be called directly from MATLAB:
Download the files ac_mex.c and utils.h from https://github.com/talisman4/2D-PcG-bodies-cell-image-analysis-HGP-Syndrome to your current working directory.
Alternatively you can download the entire ZIP archive so you'll have already all the necessary files in the working directory.
Start MATLAB and run
mex ac_mex.c
2) 8-bit 2D fluorescence cell image acquired counterstaining nuclei, lamin and PcG bodies.
Example images, used to test the code, have pixel size of 0.24 micron in x and y.
You can download the example images from the repository: https://github.com/talisman4/2D-PcG-bodies-cell-image-analysis-HGP-Syndrome
3) Create a Matlab script file to set all the necessary input parameter to run the main segmentation program.
You can find an example script: launch_script.m
4) Run
launch_script
5) It will create the directory specified in population(1).name. In the example test this directory is 2019-03-08_C001Ezh2
In the subdirectory CSV it will create a csv file with the number of nuclei for each image;
In the subdirectories Series* it will create the trinary segmented (background, nuclei, PcGs) images for each input image;
In the subdirectory SeparatedNuclei it will create a segmented image for each nucleus;
In the subdirectory STCmat it will create a Matlab file with relevant information to complete the analysis.
6) Create a Matlab script file to set all the necessary input parameters to run f11 and f16 analysis functions. 
You can find an example script: launch_script_analysis.m
7) Run
launch_script_analysis
8) In the directory specified in outpath
In the subdirectory CSV it will create csv files containing:
the number of the PcG bodies, and the area of any PcG body and all the necessary information to compute the proximity of any PcG body to nuclear periphery;
It will create the graphs related to the above measures.
