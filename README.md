# Local Isosurface Modeler

[![Build status][s1]][av] [![License: MIT][s3]][li]

[s1]: https://ci.appveyor.com/api/projects/status/uoft6p4a8gb0f0um?svg=true
[s3]: https://img.shields.io/badge/License-MIT-orange.svg

[av]: https://ci.appveyor.com/project/korfriend/localisosurfacemodeler
[li]: https://opensource.org/licenses/MIT

<br/>
This project is for a high-fidelity surface extraction technique for volumetric (CT) scans, which is based on the paper entitled "Confidence-controlled Local Isosurfacing" (under review. minor revisions in TVCG).
Detailed "get started" will be prepared with a sample code project. 
You can download the code by using Git and cloning the repository, or downloading it as zip. This will give you the full C++ source code that you must build for yourself. 

### NOTICE:
If you need binary files of built version of this source code, please contact me (korfriend@gmail.com).

### Platforms:
- Windows PC Desktop (x64)

### Requirements:

- Windows 10
- Visual Studio 2017

### Dependencies:

- Eign 3.3.4 (included)
- GL math (included)
- VisMotive (included)
- modified Adaptive multigrid solver (ver. 10.04, included) http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.04/
- Nanoflann (ver. 1.3.0, included) https://github.com/jlblancoc/nanoflann

### Build Environments
Current build environment assumes the following structure of the developement folders. As external dependencies, our VisMotive-based projects use the core APIs and libraries (https://github.com/korfriend/VisMotive-CoreAPIs/) for most of the volumetric and polygonal processing tasks. To be clear your folder structure should be something quite similar to:

    yourdevfolder/
     |
     ├──bin (built files are available in https://github.com/korfriend/VisMotive-BuiltBinary)
     │   ├──X64_Debug
     │   └──X64_Release
     └──External Projects
         ├──LocalIsosurfaceModeler (this module project https://github.com/korfriend/LocalIsosurfaceModeler/)
         │   ├──ct_modeler
         │   └──PoissonRecon (ref : http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.04/)
         ├──other module projects
         ├──...
         └──...

build folders (e.g., X64_Debug and X64_Release) should include the following dll files
- CommonUnits.dll
- vismtv_morphfilters.dll

### What Next?!
- Sample code for "get started" will be available via VisMotive framework.
