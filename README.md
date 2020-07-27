# Local Isosurface Modeler

<br/>
This project is for a high-fidelity surface extraction technique for volumetric (CT) scans, which is based on the paper entitled "Confidence-controlled Local Isosurfacing" (under review. minor revisions in TVCG).
Detailed "getting started" will be prepared with an example code project. 
You can download the code by using Git and cloning the repository, or downloading it as zip. This will give you the full C++ source code that you must build for yourself. 

### NOTICE:
We are uploading and updating the codes and corresponding libraries with certain build information.
Full source code with build environment will be available at this Aug or Sep.
If you need binary files of built version of this source code with interface header for the purpose of testing this code, please contact me (korfriend@gmail.com).

### Platforms:
- Windows PC Desktop (x64)

### Requirements:

- Windows 10
- Visual Studio 2017

### Dependencies:

- Eign 3.3.4 (included)
- GL math (included)
- Adaptive multigrid solver (ver. 10.04, included) http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.04/
- Nanoflann (ver. 1.3.0, included) https://github.com/jlblancoc/nanoflann
- VisMotive Engine (using data structure and specific modules of morphological filters and marching cubes)
- OpenCV 4.0 or higher
