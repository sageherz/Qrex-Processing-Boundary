# Qrex-Processing-Boundary
Repository containing codes and directions on how to process Qrex boundary tests (ground, wall, corner, partial). Processing Qreex data from near-boundary flights is NOT a fully autonomous procedure, you will need to understand the provided code and run sections at a time to fully process the data. All processing codes are included in the repository, but are also included on the OneDrive (AFCAD Research Group --> VLRCOE --> Qrex Post Processing --> Indoor Boundary Tests). Do not make any large changes to the codes on the OneDrive!!! If you need to make changes either create a copy file or process the data in a copied directory and not in the folder on the OneDrive.

# Processing Steps
Below is the outlined order that the processing codes should be executed.

## Step #1: Preliminary Processing 
Use the code ```processDataQREX_boundary.mat```. All near-boundary flight tests (ground, wall, partial) will use this code as the first preliminary processing code. The first section of the code will allow the user to select the boundary test type that will dictate more specifically how the data is processed in this and the following codes. 

This code will preliminarily process the recorded data, including extracting data from the .bag or .db3 file from the pi and save as a .mat file, apply scale factors, and align data. Run this code section by section and enter information where indicated. Include the following functions in your working library: 
  * ```air_density```: Computes local air density from input air pressure, temperature, and percent humidity (RECORD THESE VALUES EVERY TIME YOU CONDUCT A FLIGHT TEST).
  * ```saveData_indoor``` and ```saveData_ROS2```: Extracts data from pi and converts to structure called ```rawData```, and saves all raw data as a .mat file.
  * ```Ulog_Parser_Function```: Extracts data from Pixhawk log file and saves as a structure called ```logData```.
  * ```unbiasThrust_constant``` and ```unbiasThrust_linear```: Unbiases thrust data by either one constant value or linear fit of two values at beginning and end of data. This is based on whether you select 1 or 2 unbiasing points when prompted in MATLAB (recommended to choose two points at beginning and end of flight test to do a linear unbias).

Run this code for every flight test data file you have. If multiple flights were completed in one data file, re-run the second half of the code (after line 258) for each flight test in the file. This second half will save data for each flight test condition. 

## Step #2: Normalizing Data
Use the code ```normalizeDataQREX_boundary.mat```.  All near-boundary flight tests (ground, wall, partial) will use this code as the first preliminary processing code. The test type selection made in Step #1 will indicate which type of boundary test and how to specifically process certain data. This code will normalize thrust, RPM, CT, etc. data for boundary tests. Run this code section by section and enter information where indicated. 

For all tests except partial boundary tests with a constant altitude (z/R), ground/wall distance windows will have to be input by the user. In the second section of the code, select windows (input into variables ```Dist_norm_low``` and ```Dist_norm_high```) based on clusters of hover data points. These values are input as a function of rotor radius (R). 

For the partial boundary tests with constant altitude (z/R), averaging windows are selected earlier in Step #1 (lines 363-370 in ```processDataQREX_boundary.mat```) based on clusters of wall distance data. Remember, for these tests, the laser range finder is pointed away from the box at the wall, so the wall distance values will not correlate to the x/R values away from the box. It is important to note down the x/R data points tested in each constant z/R partial boundary test flight in order to correlate the raw wall distance data to distance from the box. 

Within this code there is options to normalize data by the "far boundary" point or hover out of boundary effect "HOBE." If you want to normalize by HOBE, you will need to define start and end HOBE points (more than one point will average HOBE sections of flight). 

The results of running this code will be plots showing Thrust, RPM, CT, Power, Total Power, and Kiel probe Pressure (if qPressure was set to 1 in Step #1). 

## Step #3: Processed Data Plots and Repeatability Comparisons
Use the code ```plottingDataQREX_boundary.mat```. Use this code when you want to combine data from multiple boundary tests or multiple trials. All near-boundary flight tests (ground, wall, partial) will use this code as for collecting all test data and preliminary plotting. In the first section of this code, the user will define what type of boundary tests they are wanting to process: Ground (w/ and w/o pressure probes), Wall/Corner (w/ and w/o pressure probes), and Partial. 

The code will automatically find and process all the tests from a selected boundary test types based on the file information provided in the EXCEL file ("GOOD_BOUNDARY_TESTS" located in the VLRCOE --> Qrex Post Processing --> Indoor Boundary Tests). This spreadsheet has a different sheet for each type of boundary test. To add additional data files to be included in the overall averaging, fill out the spreadsheet with filepath information based on the given structure. Basically, this file will direct the MATLAB code ```plottingDataQREX_boundary.mat``` to which flight test files to include in overall averaging and plotting. 

The code will create matrices containing all data from all trials to be processed in Step #4. The matrices are 4-dimensional and have the following format: ```Thrust_all(i,j,k,n)```. The index positions are described as follows:
1. ```i```: the first index denotes the z/R distance
2. ```j```: the second index denotes the x/R distance (if applicable)
3. ```k```: the third index represents either the arm number, x-y-z axis, or pressure probe port number
4. ```n```: the fourth index represents the trial number (if repeated test conditions).

These matrices containing all data from all tests is saved to a sub-directory of each boundary directory titled "ALL_PROCESSED_TESTS"). The path to save each of these all data matrices is included in lines 164-193 in ```plottingDataQREX_boundary.mat``` for each boundary test type. 

The rest of the code includes plotting to allow for comparisons between different trials of each boudnary test. The section headers indicate which type of boundary test to use for comparison plotting.

## Step #4: Total Averaging and Test-Specific Plotting
The last step is boundary-test specific and uses a different MATLAB code for each boundary type:
* Use ```plottingDataQREX_ground.mat``` for ground effect tests.
* Use ```plottingDataQREX_wall.mat``` for wall and corner effect tests.
* Use ```plottingDataQREX_partial.mat``` for partial boundary tests.

Each of these codes will compute the total averages for thrust, RPM, power, etc. from the 4D matrices created in Step #3. Each of these will also create boundary test specific plots showing the overall averaged rotor and pressure (when applicable) trends. 




