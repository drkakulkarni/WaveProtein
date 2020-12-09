# WaveProtein
**WaveProtein** is a Wavelet transformation based tool to analyze propagation of conformational changes across the polypeptide chain in Small GTPase structures. 
This tool can be extended for the analysis of other protein structures.

The package contains two codes; “PDB2RCO.c” and “WC_new_cdc42.py”

The codes can be downloaded from https://github.com/drkakulkarni/WaveProtein

### Prerequisites

PDB2RCO is a c code to convert 3D structures of protein in pdb file format to a 1D RCO signal. First you need to compile the program using following command
```
gcc –o PDB2RCO PRB2RCO.c –lm
```

Small_GTPases_WC.py requires ```PyCWT```, ```NumPy```, ```SciPy```, ```tqdm``` and ```matplotlib```.

You can install NumPy, SciPy, tqdm and matplotlib using pip or conda as
```
pip install numpy scipy tqdm matplotlib
```
Before you install PyCWT make a director “.pycwt” in your home directory as
```
mkdir ~/.pycwt
 ```
PyCWT can be download from https://github.com/regeirk/pycwt

After downloading changes the 568 and 631 lines of ```wavelet.py``` file present in the ```pycwt``` directory as 

Line 568:
```
dat = np.loadtxt(' ”your home directory”/.pycwt/kl.txt')
```
Line 631:
```
np.savetxt(' ”your home directory”/.pycwt/kl.txt', sig95)
```
Install PyCWT by running the following command in the ```pycwt``` directory
```
python setup.py install
```
Before you run WC_new_cdc42.py, you must edit the lines from 41 to 49 of the code as stated in the code. Change lines 145 & 231 as

Line 145:
```
os.rename(sim_file , '”your home directory”/.pycwt/kl.txt')
```
Line 231:
```
os.rename(' ”your home directory”/.pycwt/kl.txt', sim_file )
```

## To run:
First you need to generate the RCO file using PDB2RCO program.

```
./PDB2RCO my.pdb chain_ID starting_residue end_residue
```
Example:
```
./PDB2RCO 1an0.pdb A 5 179
```

Successful running of the program will generate two files; 1an0.txt and 1an0.csv

After generating CSV files for the two structures to be analyzed you have to run the second code for wavelet coherence analysis.
Next to run the Small_GTPases_WC.py script:
```
./WC_new_cdc42.py first_pdb_file_name second_pdb_file_name
```
Example:
```
./WC_new_cdc42.py 1an0 2qrz
```
Successful running of the script will generate two plots (png files), one showing the RCOs of both the structures and wavelet coherence plot; the other shows the resultant phase vectors of regions of interest (in this case SWI, SWII and P-loop). 

## Disclaimer:
This software is released under a BSD-style open source license. Please read the license file for further information. This routine is provided as is without any express or implied warranties whatsoever.


