# MCTA
developed by Carine P. Beatrici and Fabricio A. B. Silva

Detailed tutorial is available in the document MCTA.pdf. 

How to run the program:

In a linux or MacOS terminal:

1) After the download, compile the code using "make"; 

2) To execute the program, it is necessary to provide the configuration file in.dat. See an example below:

#file information                                        example

a) name of the fcs file                               Sample_38.fcs
b) number of considered data columns in the fcs file  11
c) number of fluorescent channels                     7
d) columns to copy into memory                        1 5 6 7 9 10 11 12 13 14 15
first the fsc, ssc and additional data (optional)
after the fluorescence channels
e) wavelength values for the channel colors           425 475 525 575 625 675 725
f) Background values                                   20.4  1000.0  12.09  19.33  1000.0  24.01  12.53
g) Compensation values                 100.00    1.63     6.10    0.00    0.70    0.00    0.00
                                        22.36  100.00     6.50    0.00    0.00    0.00    0.00
                                        13.82   29.27   100.00    6.84   37.41    0.00    0.00
                                         2.45   15.04    28.43  100.00   12.94    0.00    0.00
                                         0.00    2.44    10.84   15.73  100.00    0.00    8.74
                                         0.00    0.41     4.90    1.68    0.00  100.00   11.34
                                         0.00    0.00     0.00    0.00    0.00    6.10  100.00
h) Fluorescent channels used in the color resultant calculation       0 0 0 0 0 0 0

Helpful observations:
The in.dat must have the above informations, in order;

For itens b, c and d, to discover the number of columns of the fcs file, just execute the program with the 
in.dat file incomplete and choose from the options displayed by the program output. The counting starts at 1.

3) the wavelength values (item e) can be related or not to real filters values. 
The MCTA application automatically distributes colors evenly along the hue circle, 
ordered according to the  wavelength values informed.

4)The background labeling (item f) corresponds to the maximum fluorescence value of negative events. 
Any positive event for a given channel must have superior fluorescence intensity, ]
when compared with the background reference value.
This value must be in agreement with the data.

5) The values of the conventional compensation table (item g) must be between [0,100] 
and the value 100 is required when the same channel is considered for both dimensions.

6) To calculate the color tendency for all channels, set this value as zeros 
(one zero for each color channel - see item h) or
complete with the sequence according to the number of channels (e.g 1 2 3 4 5 6). 
By doing this, all the channels are taken into account to define the resultant color.

To calculate the resultant value only for positive events on a subset of channels, type the corresponding 
channel numbers followed by zeros, e.g. 2 4 0 0 0 0. 
In this case, only positives events for both the second and fourth channels will be considered 
and will receive a color, all others are going to stay black.

To calculate the resultant value for negative events in a set of channels, just put the number of channels
as negative numbers, e.g 2 -4 0 0 0 0. 
Using this configuration, the events are going to receive a color if and only if they are positive 
for the second channel and negative for the fourth channel (single positive events).

The cytometry files must be present in the same folder as the executable file.

This software was developed to run in Linux and MacOS based systems.

Dependencies:

1 - make; 
2 - gawk;
3 - gfortran compiler;
4 - gnuplot;
5 - R compiler;
6 - the Bioconductor package flowcore


Files in this package:

FCS2CSV.R -------------------------- R script that converts the FCS file into the text format,
                                     used internally by the MCTA application.
Makefile --------------------------- Compiles the code into an executable file.
Sample_37.fcs ---------------------- Flow cytometry experimental data used
                                     as example (control group of mice).
Sample_39.fcs ---------------------- Flow cytometry experimental data used
                                     as example (infected group of mice).
in.dat ----------------------------- This configuration file uses the Sample_39.fcs as data source. 
                                     The MCTA application will generate the color 
                                     dot-plot according to this configuration file.    
in_Sample_37.dat ------------------- Another example of configuration file, this time using 
                                     the Sample_37.fcs as input file.
MCTA.f90 - Fortran program that generates the colored dot-plot.
MCTA.pdf - Tutorial of the package.
script-gnu-color ------------------- Gnuplot script, this script uses the information in the 
                                     temporary files to generate the colored dot-plot.                                   
script-gnu-bw ---------------------- Gnuplot script, this script generates the trust map 
                                     of the colored dot-plot.
script_hue_filter------------------- Interactive dot plot script. 
                                     This script filters resultant colors based on hue values (0.0 - 360). 
                                     It requires the fluo.dat file generated by 
                                     MCTA.
script_trust_filter----------------- Interactive dot plot script. 
                                     This script filters resultant colors based on trust map values (0.0 - 1). 
                                     It requires the fluo.dat file generated by 
                                     MCTA.

