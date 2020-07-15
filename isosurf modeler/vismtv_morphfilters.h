#include "../../EngineCores/CommonUnits/VimCommon.h"

__vmstatic bool MoFiInitializeDx11();
__vmstatic bool MoFiDeinitializeDx11();

__vmstatic void MoFiSetSampleTF(int iWindowingMin, int iWindowingMax, int iMaxValue);

__vmstatic bool MorphGaussianBlur2D(ushort* pusSliceIn, ushort* pusSliceOut, vmint2 i2SizePixels, int iKernelSizeHalf, float fSigma);
__vmstatic bool MorphMedianFilter2D(ushort* pusSliceIn, ushort* pusSliceOut, vmint2 i2SizePixels, int iKernelSizeHalf);
__vmstatic bool MorphMeanFilter2D(ushort* pusSliceIn, ushort* pusSliceOut, vmint2 i2SizePixels, int iKernelSizeHalf);
__vmstatic bool MorphAdaptiveGaussianFilter2D(ushort* pusSliceIn, ushort* pusSliceOut, vmint2 i2SizePixels, int iKernelSizeHalf, float fSigma, float fThreshold);
__vmstatic bool MorphErosionFilter2D(ushort* pusSliceIn, ushort* pusSliceOut, vmint2 i2SizePixels, int iKernelSizeHalf);
__vmstatic bool MorphDilationFilter2D(ushort* pusSliceIn, ushort* pusSliceOut, vmint2 i2SizePixels, int iKernelSizeHalf);
__vmstatic bool MorphLaplacianGaussianFilter2D(ushort* pusSliceIn, ushort* pusSliceOut, vmint2 i2SizePixels, int iKernelSizeHalf, float fSigma);
// To Do
__vmstatic bool MorphAnisotropicDiffusionFilter2D(ushort* pusSliceIn, ushort* pusSliceOut, vmint2 i2SizePixels, float fKappa, float fLambda, int iOptionForEq);

// 3D Filters
__vmstatic bool MorphGaussianBlur3D(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, int iKernelSizeHalf, float fSigma, LocalProgress* _progress);
__vmstatic bool MorphMedianFilter3D(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, int iKernelSizeHalf, LocalProgress* _progress);
__vmstatic bool MorphMeanFilter3D(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, int iKernelSizeHalf, LocalProgress* _progress);
__vmstatic bool MorphAdaptiveGaussianFilter3D(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, int iKernelSizeHalf, float fSigma, float fThreshold, LocalProgress* _progress);
__vmstatic bool MorphErosionFilter3D(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, int iKernelSizeHalf, LocalProgress* _progress);
__vmstatic bool MorphDilationFilter3D(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, int iKernelSizeHalf, LocalProgress* _progress);
__vmstatic bool MorphLaplacianGaussianFilter3D(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, int iKernelSizeHalf, float fSigma, LocalProgress* _progress);
// // To Do
// __vmstatic bool MorphAnisotropicDiffusionFilter3D(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, float fKappa, float fLambda, int iOptionForEq, LocalProgress* _progress);
