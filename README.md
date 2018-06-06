# MassComp

## About MassComp
MassComp is a loseless compressor for mass spectrometry data. It compresses the mass-to-charge ratio and intensity pairs in mzXML files efficiently by calculating the hexadecimal difference of consecutive m/z values, and by searching for parts of the intensity values that match previous ones. The remaining parts of the mzXML (e.g., metadata associated to the experiments) is compressed with the general compression algorithm gzip. 


## Getting Started
Download the full project.

## Run MassComp
### Linux system
To compile:

`g++ -o masscomp masscomp.cpp tinyxml2.cpp`

To compress:

`./masscomp -c fileOri.mzXML fileMasscomp`

To decompress

`./masscomp -d fileMasscomp fileDecomp.mzXML`

To compare

`./masscomp -cmp fileOri.mzXML fileDecomp.mzXML`

### Windows system
Current implementation of the code can be run by visual studio on windows system.

Here's an example of this. Folder 'MSV000080896' is downloaded from MassIVE with id MSV000080896 and contains two mzXML files.

Run the executable file MassComp in the project.
With the hint "please input the path of files to be compressing:", input the folder path "\MSV000080896\peak\Data_mzXML" to start compressing.

With the hint "please input the path of files to be decompressing:", input the folder path "\output\MSV000080896\peak\Data_mzXML" to start decompressing.

Note: Running the application in Windows uses gzip and requires installation of Cygwin.

## Datasets
Datasets of mass spectrometry data can be downloaded from MassIVE https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp

## Authors
MassComp was created by Ruochen Yang, Xi Chen, and Idoia Ochoa at University of Illinois at Urbana-Champaign.

## Contact
If you have any problem, please email Ruochen Yang (rcyang624@126.com), Xi Chen (xichen30@illinois.edu) or Idoia Ochoa (idoia@illinois.edu).
