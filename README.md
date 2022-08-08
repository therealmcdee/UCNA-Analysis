# UCNA-Analysis
analysis code for UCNA+ prototype electron detector developed at ETSU 

---------------------------------------------------------------------------------------------------------------------------------
INPUT file number as argument "python SCRIPT INT"

base_analysis_v* == python based analysis of DT5550W csv readout (pandas)

multi_read.py & multi_totoal_HG.py == "" "" capable of reading multiple input files (only used for data sets 15-19 when 24 channels were active)

---------------------------------------------------------------------------------------------------------------------------------
Used for position reconstruction in analysis

sipm_map.txt == SiPM #, x-coord, y-coord from G4 GEARS geometry file

sipm_pos.txt == the actual G4 GEARS file of sipm locations

-----------------------------------------------------------------------------------------------------------------------------------------

stepper.ino == Arduino control software (only four positions; very early try)

stepper.py == python script to command Arduino

------------------------------------------------------------------------------------------------------------------------------------

parse_binary.cpp == C++ code to parse DT5550W binary output files

hist_disp.cpp == C++ ROOT analysis code (work in progress)
