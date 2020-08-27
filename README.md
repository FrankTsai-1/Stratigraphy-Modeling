Louisiana State University, Baton Rouge, LA 70803

February 23, 2020

This example code distribution accompanies the manuscript "Constructing large-scale complex aquifer systems with big well log data: Louisiana model", Hamid Vahdat-Aboueshagh and Frank T.-C. Tsai, submitted to Computers and Geosciences.  The code demonstrates the ideas expressed in the manuscript on constructing stratigraphy of sand-clay facies for a tile in Louisiana given a data file of well logs (WellLogs.csv) and location and land surface elevation of two-dimensional (2D) grid cells (GridTopo.csv). A stratigraphy model for the entire Louisiana domain can be constructed by running the example code for all tiles serially or parallelly in computers. 

The code includes the main file in Fortran 90 file (Lithology.F90) developed by Louisiana State University and the natural neighbor interpolation method in Fortran 77 (nn_int.f) developed by University of Bath. 

The code "Lithology.F90" is copyright of the Louisiana State University and the code "nn_int.f" is copyright of University of Bath as indicated in the source files, and is released under version 3 of the GNU General Public License.  Please see LICENSE.txt for the full text of the license.

File Descriptions
-------------------
Source code files:
	Lithology.F90  : main program that read input files, run the natural neighbor interpolation method, and produce output files.
	nn_int.f       : natural neighbor interpolation method

Input data files
	GridTopo.csv   : input file, center of cell locations (x,y) in meters and its land surface (feet). 
	WellLogs.csv   : input file, well log data file.
	Geo_model.inp  : parameter values for the program

Output files
	geomodel.prn   : A documment of model information 
	Geo-model.dat  : Modeled data file
	Geo-model_no_thin.dat  : Modeled data file with no thin facies less than ThickMin in Geo_model.inp

Compiling the Code
-----------------------

The code is simple Fortran 90 and Fortran 77, and should be able to be compiled on most modern Fortran compilers.  The primary development system was Microsoft Visual Studio 2008 and Intel Fortran Compiler (version 11.1.035) on Windows 10.

Building the code is simply a matter of compiling and linking. See the instruction file (Instruction.pdf). 

Place all input files and source code files in the same folder before compilation. Output files are produced in the same folder.

Test Data and Results
--------------------------

Examples of two datasets, corresponding to two tiles in Louisiana, are provided for testing the code, along with the expected outputs.  Test data is provided in ./TestData, and the corresponding results in ./TestResults. One example is for a tile in north Louisiana. The other example is for a tile in southeast Louisiana

../TestData/Input_example
				--GridTopo.csv
				--WellLogs.csv
				--Geo_model.inp

../TestData/Input_example2
				--GridTopo.csv
				--WellLogs.csv
				--Geo_model.inp
				
../TestData/Output_example
				--geomodel.prn
				--Geo-model.dat
				--Geo-model_no_thin.dat
				
../TestData/Output_example2
				--geomodel.prn
				--Geo-model.dat
				--Geo-model_no_thin.dat

--------------------------
The VTU files are created for the modeled data, Geo-model_no_thin.dat for the two examples. The VTU files can be visualized by the open-source program, ParaView.

../vtuFiles/
			--example.vtu
			--example2.vtu
