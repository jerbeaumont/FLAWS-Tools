# FLAWS-Tools

Open source tools for image processing and image reconstruction methods dedicated to the fluid and white matter suppression sequence (FLAWS) (1). 

The main website for FLAWS-Tools is here: https://github.com/jerbeaumont/FLAWS-Tools

This file describes the content of the FLAWS-Tools source code. The process used to compile the code is also provided. The main website is likely to contain more up-to-date information.

## Contents

- Description
- Citation
- References
- Acknowledgements
- License
- Requirements
- Third party components
- Build instruction

## Description

The FLAWS-Tools software implements the image reconstruction methods described in (1-3) for FLAWS multi-contrast imaging.

The permanent citable DOI link to the original source code used in the production of this paper is here: LINK TO THE DOI ONCE OBTAINED. 

This software has been developed by a team of researchers from [CSIRO](http://www.csiro.au/)'s [The Australian E-Health Research Centre](http://aehrc.com/) and from the [LTSI](www.ltsi.univ-rennes1.fr). See AUTHORS.txt for more details.

## Citation

If you use this program for scientific research, we would appreciate if you could cite the following papers:

- for any use of the *FLAWS_Min* image, please cite (1)
- for any use of the high contrast images (*FLAWS_HC, FLAWS_HCO, FLAWS_HC-Den, FLAWS_HCO-Den*), please cite (2-3)
- for any use of the signed inversion contrasts (*FLAWS_INV1_Signed, FLAWS_INV2_Signed*), please cite (4)
- for any use of the *FLAWS_Mip* image, please cite (4)
- for any use of the FLAWS T1 mapping (*FLAWS_T1Map*) and any use of the FLAWS B1+ corrected contrasts (*FLAWS_HC_B1PlusCorr, FLAWS_HCO_B1PlusCorr*), please cite (4)

## References

(1) Tanner et al. Fluid and white matter suppression with the MP2RAGE sequence. *J Magn Reson Imaging* 2012; 35:1063-1070.
(2) Beaumont et al. High contrast T1-weighted MRI with fluid and white matter suppression using MP2RAGE. In *IEEE Int Symp Biomed Imaging*, Venice, Italy; 2019.
(3) Beaumont et al. Multi T1-weighted contrast MRI with fluid and white matter suppression at 1.5T. *Magn Reson Imaging* 2019; 63:217-225.
(4) Paper submitted to MRM. To update before the release of the code.

## Aknowledgements

The developers aknowledge the "Region Bretagne" which partially funded this project.

## License

Copyright (c) 2020 CSIRO / LTSI-INSERM Université Rennes 1. All rights reserved.

For complete copyright, license and disclaimer of warranty information see LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the above copyright notice for more information.

## Requirements

- Git
- CMake (>= 3.1.0)
- A linux environment (the compilation on OSX and Windows environment has not been tested and is not guaranteed by the developers)

## Third party components

FLAWS-Tools are built using CMake. In addition, it depends on two other libraries:

- The Insight Toolkit, with powerful pipelining, multithreading and
  medical image IO/manipulation routines, and large helpful user community.
    https://itk.org/

- The Templatized C++ Command Line Parser Library
    http://tclap.sourceforge.net/

These libraries were selected both because they are useful and because they have BSD/MIT style licenses.

All of the libraries will be automatically downloaded and installed during the compilation of the FLAWS-Tools software. An option is available for the user that would prefer to use their own version of ITK to build the code (See the build instructions).

## Build instructions

FLAWS-Tools C++ programs are available as a superproject. To build these programs, you must follow the next instructions:
- Clone the FLAWS-Tools repository from GitHub (git clone https://github.com/jerbeaumont/FLAWS-Tools.git).
- Create a folder named "build" in the folder FLAWS-Tools.
- Run CMake in the new build folder. If you wish to change the default compilation (which automatically downloads and compiles all dependencies and tools), use ccmake.
- Build using your environment (building works with make or ninja).

An example of command lines to build FLAWS-Tools is provided below (to run in the folder in which you want to install FLAWS-Tools):

```bash
git clone https://github.com/jerbeaumont/FLAWS-Tools.git
cd FLAWS-Tools
mkdir build
cd build
cmake ../src -DCMAKE_BUILD_TYPE=Release
make
```
