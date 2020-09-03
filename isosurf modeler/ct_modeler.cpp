#include "ct_modeler.h"
#include "helpers.hpp"
//#include "vismtv_morphfilters.h"
#include "../adaptive solver/AdaptiveSolvers.h"
#pragma warning (disable:4756)

map <string, VmHMODULE> dll_import;
template <typename T>
T LoadDLL(string dll_name, string function_name)
{
	VmHMODULE hMouleLib = NULL;
	auto it = dll_import.find(dll_name);
	if (it == dll_import.end())
	{
		VMLOADLIBRARY(hMouleLib, dll_name.c_str());
		if (hMouleLib != NULL)
			dll_import[dll_name] = hMouleLib;
	}
	else
		hMouleLib = it->second;
	if (hMouleLib == NULL)
	{
		cout << "DLL import ERROR : " << dll_name << endl;
		return NULL;
	}

	T lpdll_function = NULL;
	lpdll_function = (T)VMGETPROCADDRESS(hMouleLib, function_name.c_str());
	if (lpdll_function == NULL)
		cout << "FUNC import ERROR : " << function_name << endl;
	return lpdll_function;
}

void ReleaseDlls()
{
	for (auto it = dll_import.begin(); it != dll_import.end(); it++)
		VMFREELIBRARY(it->second);
	dll_import.clear();
}

int load_foreground(const std::string& file_name, const vmobjects::VmVObjectVolume& main_volume, vmobjects::VmVObjectVolume& fore_volume,
	const int iso_value, const int downsacled_factor, const bool flip_x, const bool flip_y, const bool flip_z)
{
	VmVObjectVolume* pmain_volume = (VmVObjectVolume*)&main_volume;
	const VolumeData* main_vol_data = (const VolumeData*)(pmain_volume->GetVolumeData());
	//VolumeData* main_vol_data = main_volume.GetVolumeData();
	VolumeData* fore_vol_data = fore_volume.GetVolumeData();

	if (main_vol_data == NULL || fore_vol_data == NULL || fore_vol_data->vol_size != main_vol_data->vol_size)
		return -1;

	vmint3 vol_size = main_vol_data->vol_size;
	vmint3 bnd_size = main_vol_data->bnd_size;
	// control??

#define LOAD_VALUE 3

	if (iso_value < 0)
	{
		// when the external segmention result is processed with downsampling scheme, downsample_scale > 1
		const int downsample_scale = downsacled_factor;
		vmint3 load_mask_size = main_vol_data->vol_size / downsample_scale;
		vmint3 manual_shift = vmint3();

		const int crop_boundary = 2;
		vmfloat3 sign_scale = vmfloat3(flip_x ? -1.f : 1.f, flip_y ? -1.f : 1.f, flip_z ? -1.f : 1.f);
		vmmat44f mat_mvs2ts, mat_s, mat_t, mat_t1, mat_s1, mat_ts2mvs;
		fMatrixTranslation(&mat_t, &vmfloat3(0.5f));
		fMatrixScaling(&mat_s, &vmfloat3(1.f / (float)(load_mask_size.x - manual_shift.x), 1.f / (float)(load_mask_size.y - manual_shift.y), 1.f / (float)(load_mask_size.z - manual_shift.z)));
		fMatrixTranslation(&mat_t1, &vmfloat3(-0.5f));
		fMatrixScaling(&mat_s1, &sign_scale);
		mat_mvs2ts = mat_t * mat_s * mat_t1 * mat_s1;
		fMatrixInverse(&mat_ts2mvs, &mat_mvs2ts);

		fMatrixScaling(&mat_s, &vmfloat3(1.f / (float)(vol_size.x - 0), 1.f / (float)(vol_size.y - 0), 1.f / (float)(vol_size.z - 0)));
		vmmat44f mat_vs2ts = mat_t * mat_s * mat_t1;
		vmmat44f mat_vs2mvs = mat_vs2ts * mat_ts2mvs;

		FILE* pFile;
		fopen_s(&pFile, file_name.c_str(), "rb");
		if (pFile == NULL)
		{
			___debugout("File Open Error : " + file_name, -1);
			return -2;
		}

		typedef float LOAD_MASK; /// float
		LOAD_MASK** temp_slices;
		vmhelpers::AllocateVoidPointer2D((void***)&temp_slices, load_mask_size.z, load_mask_size.x*load_mask_size.y * sizeof(LOAD_MASK)); /// vxrDataTypeFLOAT
		for (int z = 0; z < load_mask_size.z; z++)
			fread(temp_slices[z], sizeof(LOAD_MASK), load_mask_size.x*load_mask_size.y, pFile);
		fclose(pFile);

		auto read_pixel = [&](vmint3& p) -> int
		{
			if (p.x < 0 || p.x >= load_mask_size.x || p.y < 0 || p.y >= load_mask_size.y || p.z < 0 || p.z >= load_mask_size.z) return 0;
			return temp_slices[p.z][p.x + p.y * load_mask_size.x] < 0 ? LOAD_VALUE : 0;
		};

		auto trilinearInterpolation = [](float v_0, float v_1, float v_2, float v_3, float v_4, float v_5, float v_6, float v_7,
			const vmfloat3& ratio)
		{
			float v01 = v_0 * (1.f - ratio.x) + v_1 * ratio.x;
			float v23 = v_2 * (1.f - ratio.x) + v_3 * ratio.x;
			float v0123 = v01 * (1.f - ratio.y) + v23 * ratio.y;
			float v45 = v_4 * (1.f - ratio.x) + v_5 * ratio.x;
			float v67 = v_6 * (1.f - ratio.x) + v_7 * ratio.x;
			float v4567 = v45 * (1.f - ratio.y) + v67 * ratio.y;
			return v0123 * (1.f - ratio.z) + v4567 * ratio.z;
		};

		int w = vol_size.x + bnd_size.x * 2;
		for (int z = crop_boundary; z < vol_size.z - crop_boundary; z++)
		{
			for (int y = crop_boundary; y < vol_size.y - crop_boundary; y++)
			{
				for (int x = crop_boundary; x < vol_size.x - crop_boundary; x++)
				{
					vmfloat3 pos_vs = vmfloat3((float)x, (float)y, (float)z), pos_mvs;
					fTransformPoint(&pos_mvs, &pos_vs, &mat_vs2mvs);

					vmint3 pos_idx = vmint3((int)pos_mvs.x, (int)pos_mvs.y, (int)pos_mvs.z);
					int v00 = read_pixel(pos_idx + vmint3(0, 0, 0));
					int v10 = read_pixel(pos_idx + vmint3(1, 0, 0));
					int v20 = read_pixel(pos_idx + vmint3(0, 1, 0));
					int v30 = read_pixel(pos_idx + vmint3(1, 1, 0));
					int v01 = read_pixel(pos_idx + vmint3(0, 0, 1));
					int v11 = read_pixel(pos_idx + vmint3(1, 0, 1));
					int v21 = read_pixel(pos_idx + vmint3(0, 1, 1));
					int v31 = read_pixel(pos_idx + vmint3(1, 1, 1));

					int xy_addr = x + bnd_size.x + (y + bnd_size.y) * (vol_size.x + bnd_size.x * 2);

					vmfloat3 pos_idx_f(_f3_(pos_idx.x, pos_idx.y, pos_idx.z));
					float v = trilinearInterpolation((float)v00, (float)v10, (float)v20, (float)v30, (float)v01, (float)v11, (float)v21, (float)v31, pos_mvs - pos_idx_f);
					((byte**)fore_vol_data->vol_slices)[z + bnd_size.z][xy_addr] = v >= 1.5f ? (byte)LOAD_VALUE : (byte)0;

					//int v = max(max(max(v00, v10), max(v20, v30)), max(max(v01, v11), max(v21, v31)));
					//((byte**)vol_archive_out->vol_slices)[z + vol_archive_in->bnd_size.z][xy_addr] = (byte)v;
				}
			}
		}

		// test //
		// search 101 in ui code
		/*{
			__VolSampleInfo<ushort> mask_filtered_info = Get_volsample_info<ushort>(pmain_volume, 1.f);

			int w = vol_size.x + bnd_size.x * 2;
			int h = vol_size.y + bnd_size.x * 2;
			int d = vol_size.z + bnd_size.x * 2;
			vmhelpers::AllocateVoidPointer2D((void***)&mask_filtered_info.vol_slices, d, w*h * sizeof(ushort));
			for (int z = 0; z < vol_size.z; z++)
			{
				for (int y = 0; y < vol_size.y; y++)
				{
					for (int x = 0; x < vol_size.x; x++)
					{
						int xy_addr = x + bnd_size.x + (y + bnd_size.y) * (vol_size.x + bnd_size.x * 2);
						byte __v = ((byte**)fore_vol_data->vol_slices)[z + bnd_size.z][xy_addr];
						__WriteVoxel(__v > 0 ? (ushort)60000 : (ushort)0, vmint3(x, y, z), mask_filtered_info.width_slice, EXB, (ushort**)mask_filtered_info.vol_slices);
					}
				}
			}
			LocalProgress _progress;
			MorphGaussianBlur3D(mask_filtered_info.vol_slices, mask_filtered_info.vol_slices, vmint3(w, h, d), vmint2(), 2, 1.4f, &_progress);
			for (int z = 0; z < vol_size.z; z++)
			{
				for (int y = 0; y < vol_size.y; y++)
				{
					for (int x = 0; x < vol_size.x; x++)
					{
						ushort __v = __ReadVoxel(vmint3(x, y, z), mask_filtered_info.width_slice, EXB, mask_filtered_info.vol_slices);
						int xy_addr = x + bnd_size.x + (y + bnd_size.y) * (vol_size.x + bnd_size.x * 2);
						((byte**)fore_vol_data->vol_slices)[z + bnd_size.z][xy_addr] = __v >= 30000 ? (byte)LOAD_VALUE : (byte)0;
					}
				}
			}
		}/**/

		VMSAFE_DELETE2DARRAY(temp_slices, load_mask_size.z);
	}
	else
	{
		// isovalue thresolding
		const ushort** ppdata_in = (const ushort**)main_vol_data->vol_slices;
		byte** ppdata_out = (byte**)fore_vol_data->vol_slices;

		int w = main_vol_data->vol_size.x + 2 * 2;
		int h = main_vol_data->vol_size.y + 2 * 2;
		int d = main_vol_data->vol_size.z + 2 * 2;
		for (int z = 0; z < d; z++)
			for (int xy = 0; xy < w*h; xy++)
				ppdata_out[z][xy] = ppdata_in[z][xy] >= iso_value ? LOAD_VALUE : 0;

	}

	return 0;
}

int processing_stage1(vmobjects::VmVObjectPrimitive& candidate_pts, vmobjects::VmVObjectVolume& filtered_volume, 
	const vmobjects::VmVObjectVolume& main_volume, const vmobjects::VmVObjectVolume& fore_volume, const int t_in, const int t_out,
	const float min_sample_v, const float max_sample_v/*additional feature since 20200215*/,
	const int gaussian_kernel_width/*must be odd, recommend 1 (no filter) or 3 or 5*/,
	const float simplify_grid_length/* = 1.f /*voxel unit*/,
	const float geometric_complexity_kernel_ratio/* = 0.01 /*w.r.t. bounding box of target geometry*/)
{
	if (t_in <= 0 && t_in <= 0) return -1;

	cout << ">>> # of CPU Cores : " << omp_get_num_procs() << endl;

	VmVObjectVolume* pmain_volume = (VmVObjectVolume*)&main_volume;
	VmVObjectVolume* pfore_volume = (VmVObjectVolume*)&fore_volume;
	VolumeData* main_vol_data = pmain_volume->GetVolumeData();
	VolumeData* fore_vol_data = pfore_volume->GetVolumeData();

	float sample_dist_scale = 1.f;
	__VolSampleInfo<ushort> vol_original_info = Get_volsample_info<ushort>(pmain_volume, sample_dist_scale);
	__VolSampleInfo<ushort> mask_info_inside = Get_volsample_info<ushort>(pmain_volume, sample_dist_scale);
	__VolSampleInfo<ushort> mask_info_outside = Get_volsample_info<ushort>(pmain_volume, sample_dist_scale);
	__VolSampleInfo<ushort> mask_info_mid = Get_volsample_info<ushort>(pmain_volume, sample_dist_scale);
	__VolSampleInfo<char> band_info = Get_volsample_info<char>(pmain_volume, sample_dist_scale);

	__VolSampleInfo<ushort> vol_filtered_info = Get_volsample_info<ushort>(&filtered_volume, sample_dist_scale);

	int w = mask_info_mid.vol_size.x + 2 * 2;
	int h = mask_info_mid.vol_size.y + 2 * 2;
	int d = mask_info_mid.vol_size.z + 2 * 2;

	DWORD time_stage1_start = timeGetTime();
	DWORD time_stage1_band = 0, time_stage1_rt = 0, time_stage1_sp = 0, time_stage1_gc = 0;

	{
		typedef bool(*LPDLL_CALL_MoFiInitializeDx11)();
		typedef bool(*LPDLL_CALL_MoFiDeinitializeDx11)();
		typedef bool(*LPDLL_CALL_MorphGaussianBlur3D)(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, int iKernelSizeHalf, float fSigma, LocalProgress* _progress);
		typedef bool(*LPDLL_CALL_MorphErosionFilter3D)(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, int iKernelSizeHalf, LocalProgress* _progress);
		typedef bool(*LPDLL_CALL_MorphDilationFilter3D)(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, int iKernelSizeHalf, LocalProgress* _progress);
		LPDLL_CALL_MoFiInitializeDx11 lpdll_init = LoadDLL<LPDLL_CALL_MoFiInitializeDx11>("vismtv_morphfilters", "MoFiInitializeDx11");
		LPDLL_CALL_MoFiDeinitializeDx11 lpdll_deinit = LoadDLL<LPDLL_CALL_MoFiDeinitializeDx11>("vismtv_morphfilters", "MoFiDeinitializeDx11");
		LPDLL_CALL_MorphGaussianBlur3D lpdll_gaussianblur = LoadDLL<LPDLL_CALL_MorphGaussianBlur3D>("vismtv_morphfilters", "MorphGaussianBlur3D");
		LPDLL_CALL_MorphErosionFilter3D lpdll_erosion = LoadDLL<LPDLL_CALL_MorphErosionFilter3D>("vismtv_morphfilters", "MorphErosionFilter3D");
		LPDLL_CALL_MorphDilationFilter3D lpdll_dilation = LoadDLL<LPDLL_CALL_MorphDilationFilter3D>("vismtv_morphfilters", "MorphDilationFilter3D");
		if (lpdll_init == NULL || lpdll_deinit == NULL || lpdll_gaussianblur == NULL || lpdll_erosion == NULL || lpdll_dilation == NULL)
		{
			cout << "fail to load morphological filter module" << endl;
			return -1;
		}

		LocalProgress _progress;
		lpdll_init();
		if (gaussian_kernel_width > 1)
		{
			//{
			//
			//	__VolSampleInfo<ushort> vol_resample_info = Get_volsample_info<ushort>(pmain_volume, sample_dist_scale);
			//	vmhelpers::AllocateVoidPointer2D((void***)&vol_resample_info.vol_slices, d, w*h * sizeof(ushort));
			//	for (int z = 0; z < d; z++)
			//		for (int xy = 0; xy < w*h; xy++)
			//			vol_resample_info.vol_slices[z][xy] = vol_original_info.vol_slices[z][xy] > 17000 ? vol_original_info.vol_slices[z][xy] : 0;
			//	lpdll_gaussianblur(vol_resample_info.vol_slices, vol_filtered_info.vol_slices, vmint3(w, h, d), vmint2(), gaussian_kernel_width / 2, 1.4f, &_progress);
			//	VMSAFE_DELETE2DARRAY(vol_resample_info.vol_slices, d);
			//}
			//vmhelpers::AllocateVoidPointer2D((void***)&vol_filtered_info.vol_slices, d, w*h * sizeof(ushort));
			lpdll_gaussianblur(vol_original_info.vol_slices, vol_filtered_info.vol_slices, vmint3(w, h, d), vmint2(), gaussian_kernel_width / 2, 1.4f, &_progress);
		}
		else
		{
			for (int z = 0; z < d; z++)
				memcpy(vol_filtered_info.vol_slices[z], vol_original_info.vol_slices[z], sizeof(ushort) * w * h);
		}

		mask_info_inside.vol_slices = NULL;
		mask_info_outside.vol_slices = NULL;
		mask_info_mid.vol_slices = NULL;
		band_info.vol_slices = NULL;
		time_stage1_band = timeGetTime();
		if (t_in > 0) vmhelpers::AllocateVoidPointer2D((void***)&mask_info_inside.vol_slices, d, w*h * sizeof(ushort));
		if (t_out > 0) vmhelpers::AllocateVoidPointer2D((void***)&mask_info_outside.vol_slices, d, w*h * sizeof(ushort));
		vmhelpers::AllocateVoidPointer2D((void***)&mask_info_mid.vol_slices, d, w*h * sizeof(ushort));
		vmhelpers::AllocateVoidPointer2D((void***)&band_info.vol_slices, d, w*h * sizeof(char));
		for (int z = 0; z < d; z++)
			for (int xy = 0; xy < w*h; xy++)
				mask_info_mid.vol_slices[z][xy] = ((byte**)fore_vol_data->vol_slices)[z][xy] > 0 ? 60000 : 0;

		DWORD time_morph = timeGetTime();
		lpdll_gaussianblur(mask_info_mid.vol_slices, mask_info_mid.vol_slices, vmint3(w, h, d), vmint2(), 2, 1.4f, &_progress);
		
		if (t_in > 0) lpdll_erosion(mask_info_mid.vol_slices, mask_info_inside.vol_slices, vmint3(w, h, d), vmint2(), t_in, &_progress);
		if (t_out > 0) lpdll_dilation(mask_info_mid.vol_slices, mask_info_outside.vol_slices, vmint3(w, h, d), vmint2(), t_out, &_progress);

		//cv::Mat img1(h, w, CV_8UC1);
		//cv::Mat img2(h, w, CV_8UC1);
		////cv::Mat img3(h, w, CV_8UC1);
		//for(int i = 0; i < w * h; i++)
		//{
		//	((byte*)img1.data)[i] = ((ushort**)mask_info_mid.vol_slices)[32][i] > 0 ? 255 : 0; 
		//	((byte*)img2.data)[i] = ((byte**)fore_vol_data->vol_slices)[32][i] > 0 ? 255 : 0;
		//	//if ((ushort**)mask_info_outside.vol_slices[32][i] > 0) ((byte*)img3.data)[i] = 255;
		//}
		//cv::imwrite("d:__test0.png", img1);
		//cv::imwrite("d:__test1.png", img2);
		////cv::imwrite("d:__test2.png", img3);

		if (t_in <= 0 || t_out <= 0)
		{
			if (t_in <= 0) mask_info_inside = mask_info_mid;
			if (t_out <= 0) mask_info_outside = mask_info_mid;
			mask_info_mid.vol_slices = NULL;
		}

		if (mask_info_mid.vol_slices != NULL) lpdll_gaussianblur(mask_info_mid.vol_slices, mask_info_mid.vol_slices, vmint3(w, h, d), vmint2(), 1, 1.4f, &_progress);
		lpdll_gaussianblur(mask_info_outside.vol_slices, mask_info_outside.vol_slices, vmint3(w, h, d), vmint2(), 1, 1.4f, &_progress);
		lpdll_gaussianblur(mask_info_inside.vol_slices, mask_info_inside.vol_slices, vmint3(w, h, d), vmint2(), 1, 1.4f, &_progress);
		lpdll_deinit(); // 일단 생략... MPR 때문
		make_band(band_info, (int)30000, (int)30000, mask_info_outside, mask_info_inside);
		cout << ">>> Time for band computation : " << timeGetTime() - time_stage1_band << " ms, (morph : " << timeGetTime() - time_morph << " ms)" << endl;

		ReleaseDlls();
	}

	// test // 
	//for (int z = 0; z < band_info.vol_size.z; z++)
	//	for (int y = 0; y < band_info.vol_size.y; y++)
	//		for (int x = 0; x < band_info.vol_size.x; x++)
	//		{
	//			if (__ReadVoxel(vmint3(x, y, z), band_info.width_slice, EXB, band_info.vol_slices) == 0)
	//				__WriteVoxel((byte)1, vmint3(x, y, z), band_info.width_slice, EXB, (byte**)fore_vol_data->vol_slices);
	//			else 
	//				__WriteVoxel((byte)0, vmint3(x, y, z), band_info.width_slice, EXB, (byte**)fore_vol_data->vol_slices);
	//		}

	std::vector<vmfloat3> lmax_pts;
	{
		time_stage1_rt = timeGetTime();

		vector<vmfloat3> pos_iso_pts_inside, pos_iso_pts_outside, pos_iso_pts_mid;
		// iso surfaces
		set<tuple<float, float, float>> voxel_iso_pts_inside, voxel_iso_pts_outside, voxel_iso_pts_mid;
		get_isosurface_vtx(voxel_iso_pts_inside, 30000.f, 1, mask_info_inside);
		get_isosurface_vtx(voxel_iso_pts_outside, 30000.f, 1, mask_info_outside);
		if (mask_info_mid.vol_slices != NULL) get_isosurface_vtx(voxel_iso_pts_mid, 30000.f, 1, mask_info_mid);

		for (auto itr = voxel_iso_pts_inside.begin(); itr != voxel_iso_pts_inside.end(); itr++)
		{
			vmfloat3 posVS(std::get<0>(*itr), std::get<1>(*itr), std::get<2>(*itr)), posWS;
			fTransformPoint(&posWS, &posVS, &vol_filtered_info.mat_vs2ws);
			pos_iso_pts_inside.push_back(*(vmfloat3*)&posWS);
		}
		voxel_iso_pts_inside.clear();

		for (auto itr = voxel_iso_pts_outside.begin(); itr != voxel_iso_pts_outside.end(); itr++)
		{
			vmfloat3 posVS(std::get<0>(*itr), std::get<1>(*itr), std::get<2>(*itr)), posWS;
			fTransformPoint(&posWS, &posVS, &vol_filtered_info.mat_vs2ws);
			pos_iso_pts_outside.push_back(*(vmfloat3*)&posWS);
		}
		voxel_iso_pts_outside.clear();
		if (mask_info_mid.vol_slices != NULL)
		{
			for (auto itr = voxel_iso_pts_mid.begin(); itr != voxel_iso_pts_mid.end(); itr++)
			{
				vmfloat3 posVS(std::get<0>(*itr), std::get<1>(*itr), std::get<2>(*itr)), posWS;
				fTransformPoint(&posWS, &posVS, &vol_filtered_info.mat_vs2ws);
				pos_iso_pts_mid.push_back(*(vmfloat3*)&posWS);
			}
			voxel_iso_pts_mid.clear();
		}

		vector<vmfloat3> nrl_iso_pts_inside(pos_iso_pts_inside.size());
		vector<vmfloat3> nrl_iso_pts_outside(pos_iso_pts_outside.size());
		vector<vmfloat3> nrl_iso_pts_mid(pos_iso_pts_mid.size());

		//auto __grad_normal_2 = [](const vmfloat3 pos, const __VolSampleInfo<ushort>& vol_info__, const __VolSampleInfo<ushort>& mask_info__)
		//{
		//	vmfloat3 pos_vx;
		//	fTransformPoint(&pos_vx, &pos, &vol_info__.mat_ws2vs);
		//	vmfloat3 vol_nrl = -__Safe_Gradient_by_Samples(pos_vx, vol_info__.vol_size,
		//		vol_info__.vec_grad_dirs[0], vol_info__.vec_grad_dirs[1], vol_info__.vec_grad_dirs[2],
		//		vol_info__.width_slice, EXB, vol_info__.vol_slices);
		//	vmfloat3 mask_nrl = -__Safe_Gradient_by_Samples(pos_vx, mask_info__.vol_size,
		//		mask_info__.vec_grad_dirs[0], mask_info__.vec_grad_dirs[1], mask_info__.vec_grad_dirs[2],
		//		mask_info__.width_slice, EXB, mask_info__.vol_slices);
		//
		//	if (fDotVector(&vol_nrl, &mask_nrl) < 0)
		//		vol_nrl *= -1.f;
		//
		//	float len = fLengthVector(&vol_nrl);
		//	if (len < FLT_EPSILON) return vmfloat3();
		//	return vol_nrl / len;
		//};

		auto __grad_normal = [](const vmfloat3 pos, const __VolSampleInfo<ushort>& vol_info__, const __VolSampleInfo<ushort>& mask_info__)
		{
			vmfloat3 pos_vx;
			fTransformPoint(&pos_vx, &pos, &mask_info__.mat_ws2vs);
			vmfloat3 mask_nrl = -__Safe_Gradient_by_Samples(pos_vx, mask_info__.vol_size,
				mask_info__.vec_grad_dirs[0] * 2.f, mask_info__.vec_grad_dirs[1] * 2.f, mask_info__.vec_grad_dirs[2] * 2.f,
				mask_info__.width_slice, EXB, mask_info__.vol_slices);
			float len = fLengthVector(&mask_nrl);
			if (len < FLT_EPSILON) return vmfloat3();
			return mask_nrl / len;
		};

#pragma omp parallel for num_threads( omp_get_num_procs() ) // 
		for (int i = 0; i < (int)pos_iso_pts_inside.size(); i++)
			nrl_iso_pts_inside[i] = __grad_normal(pos_iso_pts_inside[i], vol_filtered_info, mask_info_inside);
#pragma omp parallel for num_threads( omp_get_num_procs() ) // 
		for (int i = 0; i < (int)pos_iso_pts_outside.size(); i++)
			nrl_iso_pts_outside[i] = __grad_normal(pos_iso_pts_outside[i], vol_filtered_info, mask_info_outside);
#pragma omp parallel for num_threads( omp_get_num_procs() ) // 
		for (int i = 0; i < (int)pos_iso_pts_mid.size(); i++)
			nrl_iso_pts_mid[i] = __grad_normal(pos_iso_pts_mid[i], vol_filtered_info, mask_info_mid);

		cout << "          # of iso surfaces pts (in, out, mid) : " << pos_iso_pts_inside.size() << ", " << pos_iso_pts_outside.size() << ", " << pos_iso_pts_mid.size() << endl;
		cout << "          initial points preparation time : " << timeGetTime() - time_stage1_rt << " ms" << endl;

		//auto __check_nrl_length = [](vector<vmfloat3>& nrl)
		//{
		//	int count = 0;
		//	for (int i = 0; i < (int)nrl.size(); i++)
		//		if (fLengthVectorSq(&nrl[i]) < FLT_EPSILON)
		//			count++;
		//	cout << "******* nrl with zero-length : " << count << ", " << nrl.size() << endl;
		//};
		//__check_nrl_length(pos_iso_pts_outside);
		//__check_nrl_length(pos_iso_pts_inside);
		//__check_nrl_length(pos_iso_pts_mid);

		float _min_v = min_sample_v; //18600.f;// 
		float _max_v = max_sample_v; //27200.f;// 
		std::vector<vmfloat3> lmax_in_pts, lmax_midin_pts, lmax_midout_pts;
		float sample_dist = vol_filtered_info.min_sample_dist * 0.5f;
		raytraversal_gradmax(lmax_pts, pos_iso_pts_outside, nrl_iso_pts_outside, 100, -sample_dist, _min_v, _max_v, band_info, vol_filtered_info);
		raytraversal_gradmax(lmax_in_pts, pos_iso_pts_inside, nrl_iso_pts_inside, 100, sample_dist, _min_v, _max_v, band_info, vol_filtered_info);
		lmax_pts.insert(lmax_pts.end(), lmax_in_pts.begin(), lmax_in_pts.end());

		if (mask_info_mid.vol_slices != NULL)
		{
			raytraversal_gradmax(lmax_midout_pts, pos_iso_pts_mid, nrl_iso_pts_mid, 100, sample_dist, _min_v, _max_v, band_info, vol_filtered_info);
			raytraversal_gradmax(lmax_midin_pts, pos_iso_pts_mid, nrl_iso_pts_mid, 100, -sample_dist, _min_v, _max_v, band_info, vol_filtered_info);
			if (lmax_midout_pts.size() != lmax_midin_pts.size())
			{
				std::cout << "raycast_local_maximum ERROR " << lmax_midout_pts.size() << ", " << lmax_midin_pts.size() << std::endl;
				return -1;
			}
			for (int i = 0; i < (int)lmax_midout_pts.size(); i++)
			{
				vmfloat3 pos_in_vs, pos_out_vs;
				fTransformPoint(&pos_in_vs, &lmax_midin_pts[i], &vol_filtered_info.mat_ws2vs);
				fTransformPoint(&pos_out_vs, &lmax_midout_pts[i], &vol_filtered_info.mat_ws2vs);
				if (fLengthVectorSq(&SGRAD(pos_out_vs, vol_filtered_info)) > fLengthVectorSq(&SGRAD(pos_in_vs, vol_filtered_info)))
					lmax_pts.push_back(lmax_midout_pts[i]);
				else
					lmax_pts.push_back(lmax_midin_pts[i]);
			}
		}/**/

		cout << ">>> Time for ray traversal : " << timeGetTime() - time_stage1_rt << " ms" << endl;
	}

	// simplification //
	vector<vmfloat3> pos_candidate_pts;
	vmfloat3 aabb_diff;
	{
		time_stage1_sp = timeGetTime();
		simplify_points_ugrid(pos_candidate_pts, aabb_diff, lmax_pts, vol_filtered_info.min_sample_dist * simplify_grid_length);
		cout << "          # of pts (before simplification, after simplification) : " << lmax_pts.size() << ", " << pos_candidate_pts.size() << endl;
		cout << ">>> Time for simplification : " << timeGetTime() - time_stage1_sp << " ms" << endl;
		lmax_pts.clear();
	}
	//calculate gc
	{
		time_stage1_gc = timeGetTime();
		float diag_length = fLengthVector(&aabb_diff);
		const float kernel_radius = max(geometric_complexity_kernel_ratio * diag_length, vol_filtered_info.min_sample_dist * 3.f);
		cout << "          kernel radius for geometric complexity : " << kernel_radius << " mm"<< endl;

		vector<vmfloat3> nrl_candidate_pts;
		vector<vmfloat2> gm_gc_candidate_pts;
		compute_geometry_info(nrl_candidate_pts, gm_gc_candidate_pts, pos_candidate_pts, kernel_radius, vol_filtered_info);

		cout << ">>> Time for computing GC : " << timeGetTime() - time_stage1_gc << " ms" << endl;
		
		register_pobj(candidate_pts, pos_candidate_pts, nrl_candidate_pts, NULL);

		vector<float> metric_gm(pos_candidate_pts.size());
		vector<float> metric_gc(pos_candidate_pts.size());
		float min_gm(FLT_MAX), max_gm(0);
		for (uint i = 0; i < (uint)pos_candidate_pts.size(); i++)
		{
			metric_gm[i] = gm_gc_candidate_pts[i].x;
			min_gm = min(metric_gm[i], min_gm);
			max_gm = max(metric_gm[i], max_gm);
			metric_gc[i] = gm_gc_candidate_pts[i].y;
		}
		VmLObject* lobj = candidate_pts.GetBufferObject();
		lobj->ReplaceOrAddBufferPtr("_vlist_float_gm", (float*)&metric_gm[0], (int)metric_gm.size(), sizeof(float));
		lobj->ReplaceOrAddBufferPtr("_vlist_float_gc", (float*)&metric_gc[0], (int)metric_gc.size(), sizeof(float));
		candidate_pts.RegisterCustomParameter("_double_minpitch", (double)vol_original_info.min_sample_dist);

		vector<int> jet_color_map(COLORMAP_SIZE);
		fill_jet_colormap(&jet_color_map[0], COLORMAP_SIZE);
		lobj->ReplaceOrAddBufferPtr("_vlist_int_jetmap", (int*)&jet_color_map[0], (int)jet_color_map.size(), sizeof(int));
		vector<int> cool_color_map(COLORMAP_SIZE);
		fill_cool_colormap(&cool_color_map[0], COLORMAP_SIZE);
		lobj->ReplaceOrAddBufferPtr("_vlist_int_coolmap", (int*)&cool_color_map[0], (int)cool_color_map.size(), sizeof(int));

		auto ___otsu_threshold = [](const vector<float>& metric_gm, const float min_v, const float max_v, const int bin_num)->float
		{
			int _bin_num = (int)(max_v - min_v) > bin_num ? bin_num : (int)(max_v - min_v);
			ullong* hgm_raw_trajectories = new ullong[_bin_num];
			ZeroMemory(hgm_raw_trajectories, sizeof(ullong) * _bin_num);

			float _max_v = min_v == max_v ? max_v + 1 : max_v;

			// Histogram of Gradient Magnitude (HGM)
			for (size_t i = 0; i < metric_gm.size(); i++)
			{
				float normalized = (metric_gm[i] - min_v) / (_max_v - min_v);
				normalized = max(min(normalized, 1.f), 0.f);
				int idx = (int)(normalized * (float)(_bin_num - 1));
				hgm_raw_trajectories[idx]++;
			}

			// threshold_strong and threshold_weak
			// try using otsu0
			int otsu_idx = OtsuThresholdValue(hgm_raw_trajectories, _bin_num, 0, _bin_num - 1);
			delete[] hgm_raw_trajectories;
			float thd = (float)otsu_idx / (float)(_bin_num - 1) * (_max_v - min_v) + min_v;
			//printf("** otsu thresholds : idx(%d), thd(%f)\n", otsu_idx, thd);

			if (otsu_idx == 0)
			{
				cout << "** insufficient samples of histogram value range ==> mean value" << endl;
				return (_max_v + min_v) / 2.f;
			}

			return thd;
		};

		float otsu_gm = ___otsu_threshold(metric_gm, min_gm, max_gm, 10000);
		otsu_gm = ___otsu_threshold(metric_gm, min_gm, otsu_gm, 10000);
		candidate_pts.RegisterCustomParameter("_double_mingm", (double)min_gm);
		candidate_pts.RegisterCustomParameter("_double_maxgm", (double)max_gm);
		candidate_pts.RegisterCustomParameter("_double_otsugm", (double)otsu_gm);
	}

	//if (gaussian_kernel_width > 1)
	//	VMSAFE_DELETE2DARRAY(vol_filtered_info.vol_slices, d);
	VMSAFE_DELETE2DARRAY(mask_info_inside.vol_slices, d);
	VMSAFE_DELETE2DARRAY(mask_info_outside.vol_slices, d);
	VMSAFE_DELETE2DARRAY(mask_info_mid.vol_slices, d);
	VMSAFE_DELETE2DARRAY(band_info.vol_slices, d);
	return 0;
}

#define __STRONG_ONLY_GM
int update_strongweak_points(vmobjects::VmVObjectPrimitive& candidate_pts, const float g_h, const float mu_u)
{
	cout << "     * g_h  : " << g_h << endl;
	cout << "     * mu_u : " << mu_u << endl;

	double min_pitch = 0;
	candidate_pts.GetCustomParameter("_double_minpitch", data_type::dtype<double>(), &min_pitch);
	int* jet_color_map_ptr;
	size_t _t_size;
	int num_ele_jet_color_map = 0;
	VmLObject* lobj = candidate_pts.GetBufferObject();
	lobj->LoadBufferPtr("_vlist_int_jetmap", (void**)&jet_color_map_ptr, _t_size, &num_ele_jet_color_map);

	float* metric_gm_ptr;
	float* metric_gc_ptr;
	int num_ele_metric = 0;
	lobj->LoadBufferPtr("_vlist_float_gm", (void**)&metric_gm_ptr, _t_size, &num_ele_metric);
	lobj->LoadBufferPtr("_vlist_float_gc", (void**)&metric_gc_ptr, _t_size, &num_ele_metric);

	PrimitiveData* prim_data = candidate_pts.GetPrimitiveData();
	vmfloat3* pos_vtx = prim_data->GetVerticeDefinition("POSITION");
	vmfloat3* clr_vtx = prim_data->GetVerticeDefinition("TEXCOORD0");

	float r_sq = (float)(min_pitch * 2.0); // epsilon
	cout << "     * r   : " << r_sq << endl;
	r_sq *= r_sq;
	float g_l = (float)g_h * 0.05f;

	// test : Fig.6
	int* cool_color_map_ptr;
	lobj->LoadBufferPtr("_vlist_int_coolmap", (void**)&cool_color_map_ptr, _t_size, &num_ele_jet_color_map);
	// test : Fig.6 (a)
	/*{
		float _g_h = 10000; // Piston model : 10000, Knee : 90
		float _g_l = _g_h * 0.05;
		for (uint i = 0; i < prim_data->num_vtx; i++)
		{
			float gm = metric_gm_ptr[i];
			gm = min(gm, _g_h);
			gm = max(gm, _g_l);
			int clr_idx = (int)((gm - _g_l) / (_g_h - _g_l) * (float)(COLORMAP_SIZE - 1));
			int color = cool_color_map_ptr[clr_idx];
			convert_int_to_float3(color, clr_vtx[i]);
		}
		return 0;
	}/**/

	// test : Fig.6 (b)
	/*{
		for (uint i = 0; i < prim_data->num_vtx; i++)
		{
			float gc = metric_gc_ptr[i];
			int clr_idx = (int)(gc * (float)(COLORMAP_SIZE - 1));
			int color = jet_color_map_ptr[clr_idx];
			convert_int_to_float3(color, clr_vtx[i]);
			//clr_vtx[i] = vmfloat3(1, 1, 0); // Fig.13 candidate points (displaying yellow points)
		}
		return 0;
	}/**/

	for (uint i = 0; i < prim_data->num_vtx; i++)
	{
		float gm = metric_gm_ptr[i];
		if (gm >= g_h)
			clr_vtx[i] = vmfloat3(1);
		else if (gm <= g_l)
			clr_vtx[i] = vmfloat3(0);
		else
		{
			gm = min(g_h, max(g_l, gm));
			int clr_idx = (int)((gm - g_l) / (g_h - g_l) * (float)(COLORMAP_SIZE - 1));
			//int color = jet_color_map_ptr[clr_idx];
			int color = cool_color_map_ptr[clr_idx]; // Fig.6 (c) (d)
			convert_int_to_float3(color, clr_vtx[i]);
		}
	
	}
	// test : Fig.6 (c)
	//return 0;

	for (uint i = 0; i < prim_data->num_vtx; i++)
	{
		if (metric_gc_ptr[i] >= (float)mu_u)
		{
			vmfloat3 pos_src = pos_vtx[i];
			std::vector<std::pair<size_t, float>> ret_matches;
			const int nMatches = (int)candidate_pts.KDTSearchRadius(pos_src, r_sq, false, ret_matches);
			for (int j = 0; j < nMatches; j++)
			{
				int idx_neighbor = (int)ret_matches[j].first;
#ifdef __STRONG_ONLY_GM
				if(metric_gm_ptr[idx_neighbor] < g_h)
#endif
					clr_vtx[idx_neighbor] = vmfloat3();
			}
		}
	}

	return 0;
}

int processing_stage2(vmobjects::VmVObjectPrimitive& reliable_pts, vmobjects::VmVObjectPrimitive& holefilling_pts,
	const vmobjects::VmVObjectVolume& sample_volume, const vmobjects::VmVObjectPrimitive& candidate_pts, const float g_h, const float mu_u,
	const float epsilon /*= 2.f, voxel unit, connecting radius*/,
	const int connectivity_criterion /*= 5, number of connecting points for the connectivity test*/,
	const float epsilon_b /*= 3.f, voxel unit, for detecting the hole boundary*/,
	const float angle_criterion /*= 90.f, degree, for detecting the hole boundary*/)
{
	//std::vector<int>& cluster_cntpts ==> holefilling_pts
	float sample_dist_scale = 1.f;
	VmVObjectVolume* pmain_volume = (VmVObjectVolume*)&sample_volume;
	__VolSampleInfo<ushort> vol_info = Get_volsample_info<ushort>(pmain_volume, sample_dist_scale);
	float real_epsilon = vol_info.min_sample_dist * epsilon;
	float real_epsilon_b = vol_info.min_sample_dist * epsilon_b;

	VmLObject* lobj = ((VmVObjectPrimitive*)&candidate_pts)->GetBufferObject();
	size_t _t_size;
	float* metric_gm_ptr;
	float* metric_gc_ptr;
	int num_ele_metric = 0;
	lobj->LoadBufferPtr("_vlist_float_gm", (void**)&metric_gm_ptr, _t_size, &num_ele_metric);
	lobj->LoadBufferPtr("_vlist_float_gc", (void**)&metric_gc_ptr, _t_size, &num_ele_metric);

	VmVObjectPrimitive* pcandidate_pts = (VmVObjectPrimitive*)&candidate_pts;
	PrimitiveData* prim_data = pcandidate_pts->GetPrimitiveData();
	vmfloat3* pos_cand_vtx = prim_data->GetVerticeDefinition("POSITION");
	vmfloat3* nrl_cand_vtx = prim_data->GetVerticeDefinition("NORMAL");

	DWORD time_stage2_start = timeGetTime();
	DWORD time_stage2_strongweak = 0, time_stage2_confidence = 0, time_stage2_holefilling = 0;

	float e_u_sq = (float)(vol_info.min_sample_dist * 2.0);
	e_u_sq *= e_u_sq;
	vector<char> state_strongweak(prim_data->num_vtx, 0);
	for (uint i = 0; i < prim_data->num_vtx; i++)
	{
		if (metric_gc_ptr[i] >= (float)mu_u
#ifdef __STRONG_ONLY_GM
			&& metric_gm_ptr[i] < g_h
#endif
			) // optional ... this results in better..
		{
			vmfloat3 pos_src = pos_cand_vtx[i];
			std::vector<std::pair<size_t, float>> ret_matches;
			const int nMatches = (int)pcandidate_pts->KDTSearchRadius(pos_src, e_u_sq, false, ret_matches);
			for (int j = 0; j < nMatches; j++)
			{
				int idx_neighbor = (int)ret_matches[j].first;
				state_strongweak[idx_neighbor] = -1;
			}
		}
	}

	float g_l = (float)g_h * 0.05f;
	vector<vmfloat3> pos_strong_pts, nrl_strong_pts, pos_weak_pts, nrl_weak_pts;
	vector<vmfloat3> pos_remain_pts, nrl_remain_pts;
	{
		time_stage2_strongweak = timeGetTime();
		// identifying strong and weak
		for (uint i = 0; i < prim_data->num_vtx; i++)
		{
			if (state_strongweak[i] == -1)
			{
				pos_remain_pts.push_back(pos_cand_vtx[i]);
				nrl_remain_pts.push_back(nrl_cand_vtx[i]);
			}
			else
			{
				float gm = metric_gm_ptr[i];
				if (gm >= g_h)
				{
					pos_strong_pts.push_back(pos_cand_vtx[i]);
					nrl_strong_pts.push_back(nrl_cand_vtx[i]);
				}
				else if (gm > g_l)
				{
					pos_weak_pts.push_back(pos_cand_vtx[i]);
					nrl_weak_pts.push_back(nrl_cand_vtx[i]);
				}
			}
		}
		cout << "          # of pts (strong, weak) : " << pos_strong_pts.size() << ", " << pos_weak_pts.size() << endl;
		cout << ">>> Time for strong/weak pts identification: " << timeGetTime() - time_stage2_strongweak << " ms" << endl;
	}
	state_strongweak.clear(); // recover RAM

	//vector<vmfloat3> bnd_pts; // test

	auto compute_num_boundary_pts = [](std::map<int, std::vector<int>>& map_clusters, std::map<int, int>& cluster_connectivity,
		const vector<vmfloat3>& pos_ref_pts, const vector<vmfloat3>& nrl_ref_pts, const vector<vmfloat3>& pos_dst_pts, const vector<vmfloat3>& nrl_dst_pts,
		float epsilon, float epsilon_b, vector<vmfloat3>* pos_ref_bnd_ptr)
	{
		float e_sq = epsilon * epsilon;
		PointCloud<float, vmfloat3> pc_ref(pos_ref_pts);
		kd_tree_t kdt_ref(3, pc_ref, nanoflann::KDTreeSingleIndexAdaptorParams(10));
		kdt_ref.buildIndex();
		vector<int> idx_bnd_ref_pts;
		detecting_boundary(idx_bnd_ref_pts, pos_ref_pts, nrl_ref_pts, kdt_ref, epsilon_b, true);
		vector<vmfloat3> pos_ref_bnd(idx_bnd_ref_pts.size());
		for (int i = 0; i < (int)idx_bnd_ref_pts.size(); i++)
			pos_ref_bnd[i] = pos_ref_pts[idx_bnd_ref_pts[i]];
		if (pos_ref_bnd_ptr) // test
			*pos_ref_bnd_ptr = pos_ref_bnd;

		// identifying confidence and rest
		PointCloud<float, vmfloat3> pc_weak(pos_dst_pts);
		kd_tree_t kdt_weak(3, pc_weak, nanoflann::KDTreeSingleIndexAdaptorParams(10));
		kdt_weak.buildIndex();

		clustering(map_clusters, pos_dst_pts, kdt_weak, nrl_dst_pts, epsilon);
		cout << "          # of clusters : " << map_clusters.size() << endl;

		PointCloud<float, vmfloat3> pc_ref_bnd(pos_ref_bnd);
		kd_tree_t kdt_ref_bnd(3, pc_ref_bnd, nanoflann::KDTreeSingleIndexAdaptorParams(10));
		kdt_ref_bnd.buildIndex();
		for (auto it = map_clusters.begin(); it != map_clusters.end(); it++)
		{
			const std::vector<int>& idx_cluster = it->second;

			int connecting_cnt = 0;
			// NumConnectingBoundaryPoints(idx_cluster, Ps's bnd, epsilon)
			for (int i = 0; i < (int)idx_cluster.size(); i++)
			{
				const vmfloat3& pos_src = pos_dst_pts[idx_cluster[i]];
				size_t nb_idx;
				float nb_dist_sq;
				kdt_ref_bnd.knnSearch((float*)&pos_src, 1, &nb_idx, &nb_dist_sq);
				if (nb_dist_sq < e_sq) connecting_cnt++;
			}

			cluster_connectivity[it->first] = connecting_cnt;
		}
	};

	vector<vmfloat3> pos_reliable_pts = pos_strong_pts, nrl_reliable_pts = nrl_strong_pts;
	{
		time_stage2_confidence = timeGetTime();

		map<int, std::vector<int>> map_weak_clusters;
		map<int, int> cluster_connectivity;
		compute_num_boundary_pts(map_weak_clusters, cluster_connectivity, pos_strong_pts, nrl_strong_pts, pos_weak_pts, nrl_weak_pts, real_epsilon, real_epsilon_b, NULL);
		// test
		//holefilling_pts.GetBufferObject()->ReplaceOrAddBufferPtr("_buffer_pos_spheres", &bnd_pts[0], (int)bnd_pts.size(), sizeof(vmfloat3));

		for (auto it = map_weak_clusters.begin(); it != map_weak_clusters.end(); it++)
		{
			int connecting_cnt = cluster_connectivity[it->first];
			const std::vector<int>& idx_cluster = it->second;
			if (connecting_cnt >= connectivity_criterion)
			{
				for (int i = 0; i < (int)idx_cluster.size(); i++)
				{
					pos_reliable_pts.push_back(pos_weak_pts[idx_cluster[i]]);
					nrl_reliable_pts.push_back(nrl_weak_pts[idx_cluster[i]]);
				}
			}
			else
			{
				for (int i = 0; i < (int)idx_cluster.size(); i++)
				{
					pos_remain_pts.push_back(pos_weak_pts[idx_cluster[i]]);
					nrl_remain_pts.push_back(nrl_weak_pts[idx_cluster[i]]);
				}
			}
		}
		cout << "          # of pts (reliable, remain) : " << pos_reliable_pts.size() << ", " << pos_remain_pts.size() << endl;
		cout << ">>> Time for reliable pts identification: " << timeGetTime() - time_stage2_confidence << " ms" << endl;
	}
	pos_strong_pts.clear(); // recover RAM
	nrl_strong_pts.clear(); // recover RAM
	pos_weak_pts.clear(); // recover RAM
	nrl_weak_pts.clear(); // recover RAM

	vector<vmfloat3> pos_holefill_pts, nrl_holefill_pts;
	vector<int> size_clusters;
	vector<int> connectivity_clusters;
	int min_bnd_cnt(100000000), max_bnd_cnt(0);
	{
		time_stage2_holefilling = timeGetTime();

		map<int, std::vector<int>> map_remain_clusters;
		map<int, int> cluster_connectivity;
		compute_num_boundary_pts(map_remain_clusters, cluster_connectivity, pos_reliable_pts, nrl_reliable_pts, pos_remain_pts, nrl_remain_pts, real_epsilon, real_epsilon_b, NULL);

		for (auto it = map_remain_clusters.begin(); it != map_remain_clusters.end(); it++)
		{
			int connecting_cnt = cluster_connectivity[it->first];
			const std::vector<int>& idx_cluster = it->second;
			if (connecting_cnt >= connectivity_criterion)
			{
				size_clusters.push_back((int)idx_cluster.size());
				for (int i = 0; i < (int)idx_cluster.size(); i++)
				{
					pos_holefill_pts.push_back(pos_remain_pts[idx_cluster[i]]);
					nrl_holefill_pts.push_back(nrl_remain_pts[idx_cluster[i]]);
				}
			}
		}
		map_remain_clusters.clear(); // recover RAM
		cluster_connectivity.clear();

		vector<vmfloat3> pos_rh_pts = pos_reliable_pts, nrl_rh_pts = nrl_reliable_pts;
		pos_rh_pts.insert(pos_rh_pts.end(), pos_holefill_pts.begin(), pos_holefill_pts.end());
		nrl_rh_pts.insert(nrl_rh_pts.end(), nrl_holefill_pts.begin(), nrl_holefill_pts.end());

		PointCloud<float, vmfloat3> pc_rh(pos_rh_pts);
		kd_tree_t kdt_rh(3, pc_rh, nanoflann::KDTreeSingleIndexAdaptorParams(10));
		kdt_rh.buildIndex();
		vector<int> idx_bnd_rh_pts;
		detecting_boundary(idx_bnd_rh_pts, pos_rh_pts, nrl_rh_pts, kdt_rh, real_epsilon_b * 1.5f, true);
		
		{ // becomes official ui
			vector<vmfloat3> pos_rh_hole_pts(idx_bnd_rh_pts.size());
			for (int i = 0; i < (int)idx_bnd_rh_pts.size(); i++)
				pos_rh_hole_pts[i] = pos_rh_pts[idx_bnd_rh_pts[i]];
			holefilling_pts.GetBufferObject()->ReplaceOrAddBufferPtr("_buffer_pos_spheres", &pos_rh_hole_pts[0], (int)idx_bnd_rh_pts.size(), sizeof(vmfloat3));
		}		
		// test
		//{
		//	PointCloud<float, vmfloat3> pc_r(pos_reliable_pts);
		//	kd_tree_t kdt_r(3, pc_r, nanoflann::KDTreeSingleIndexAdaptorParams(10));
		//	kdt_r.buildIndex();
		//	vector<int> idx_bnd_r_pts;
		//	detecting_boundary(idx_bnd_r_pts, pos_reliable_pts, nrl_reliable_pts, kdt_r, real_epsilon_b, true);
		//	vector<vmfloat3> pos_r_hole_pts(idx_bnd_r_pts.size());
		//	for (int i = 0; i < (int)idx_bnd_r_pts.size(); i++)
		//		pos_r_hole_pts[i] = pos_reliable_pts[idx_bnd_r_pts[i]];
		//	holefilling_pts.GetBufferObject()->ReplaceOrAddBufferPtr("_buffer_pos_spheres", &pos_r_hole_pts[0], (int)idx_bnd_r_pts.size(), sizeof(vmfloat3));
		//}

		vector<bool> bndstate_holefill_pts(pos_holefill_pts.size(), false);
		for (int i = 0; i < (int)idx_bnd_rh_pts.size(); i++)
		{
			int idx = idx_bnd_rh_pts[i] - (int)pos_reliable_pts.size();
			if (idx >= 0)
				bndstate_holefill_pts[idx] = true;
		}
		int prev_idx_offset = 0;
		connectivity_clusters.assign((int)size_clusters.size(), 0);
		for (int i = 0; i < (int)size_clusters.size(); i++)
		{
			int num_cluster_pts = size_clusters[i];
			int bnd_count = 0;
			for (int j = 0; j < num_cluster_pts; j++)
			{
				int idx = prev_idx_offset + j;
				if (bndstate_holefill_pts[idx]) bnd_count++;
			}
			min_bnd_cnt = min(min_bnd_cnt, bnd_count);
			max_bnd_cnt = max(max_bnd_cnt, bnd_count);
			connectivity_clusters[i] = bnd_count;
			prev_idx_offset += num_cluster_pts;
		}

		cout << "          # of bnd pts of a cluster (min, max) : " << min_bnd_cnt << ", " << max_bnd_cnt << endl;
		cout << "          # of pts (remain, holefilling) : " << pos_remain_pts.size() << ", " << pos_holefill_pts.size() << endl;
		cout << ">>> Time for holefilling pts identification: " << timeGetTime() - time_stage2_holefilling << " ms" << endl;
	}
	pos_remain_pts.clear(); // recover RAM
	nrl_remain_pts.clear(); // recover RAM

	{
		if (pos_reliable_pts.size() == 0)
		{
			cout << "at least a reliable point is necessary!" << endl;
			return -1;
		}

		{
			register_pobj(reliable_pts, pos_reliable_pts, nrl_reliable_pts, NULL);
			reliable_pts.UpdateKDTree();
			vector<float> value_r_pts(pos_reliable_pts.size());
#pragma omp parallel for num_threads( omp_get_num_procs() ) // 
			for (int i = 0; i < (int)pos_reliable_pts.size(); i++)
			{
				vmfloat3 pos_sample;
				fTransformPoint(&pos_sample, &pos_reliable_pts[i], &vol_info.mat_ws2vs);
				float sample_v = __Safe_TrilinearSample<ushort>(pos_sample, vol_info.vol_size, vol_info.width_slice, EXB, vol_info.vol_slices, 0);
				value_r_pts[i] = sample_v;
			}
			reliable_pts.GetBufferObject()->ReplaceOrAddBufferPtr("_vlist_float_values", (float*)&value_r_pts[0], (int)value_r_pts.size(), sizeof(float));;
		}
		{
			if (pos_holefill_pts.size() == 0)
			{
				return -1;
			}
			register_pobj(holefilling_pts, pos_holefill_pts, nrl_holefill_pts, NULL);
			holefilling_pts.UpdateKDTree();

			vector<int> jet_color_map(COLORMAP_SIZE);
			fill_jet_colormap(&jet_color_map[0], COLORMAP_SIZE);
			VmLObject* lobj = holefilling_pts.GetBufferObject();

			lobj->ReplaceOrAddBufferPtr("_vlist_int_jetmap", (int*)&jet_color_map[0], (int)jet_color_map.size(), sizeof(int));
			lobj->ReplaceOrAddBufferPtr("_vlist_int_numclusters", (int*)&size_clusters[0], (int)size_clusters.size(), sizeof(int));
			lobj->ReplaceOrAddBufferPtr("_vlist_int_connectivityclusters", (int*)&connectivity_clusters[0], (int)connectivity_clusters.size(), sizeof(int));

			// confidence //
			vector<float> confidence_hf_pts(pos_holefill_pts.size());
			vector<float> value_hf_pts(pos_holefill_pts.size());
#pragma omp parallel for num_threads( omp_get_num_procs() ) // 
			for (int i = 0; i < (int)pos_holefill_pts.size(); i++)
			{
				vmfloat3 pos_sample;
				fTransformPoint(&pos_sample, &pos_holefill_pts[i], &vol_info.mat_ws2vs);
				float sample_v = __Safe_TrilinearSample<ushort>(pos_sample, vol_info.vol_size, vol_info.width_slice, EXB, vol_info.vol_slices, 0);
				vmfloat3 grad = __Safe_Gradient_by_Samples(pos_sample, vol_info.vol_size,
					vol_info.vec_grad_dirs[0], vol_info.vec_grad_dirs[1], vol_info.vec_grad_dirs[2],
					vol_info.width_slice, EXB, vol_info.vol_slices);
				float gm = fLengthVector(&grad);
				confidence_hf_pts[i] = min(max((gm - g_l) / (g_h - g_l), 0.f), 1.f);
				value_hf_pts[i] = sample_v;
			}
			lobj->ReplaceOrAddBufferPtr("_vlist_float_confidences", (float*)&confidence_hf_pts[0], (int)confidence_hf_pts.size(), sizeof(float));
			lobj->ReplaceOrAddBufferPtr("_vlist_float_values", (float*)&value_hf_pts[0], (int)value_hf_pts.size(), sizeof(float));

			holefilling_pts.RegisterCustomParameter("_double_minpitch", (double)vol_info.min_sample_dist);
			holefilling_pts.RegisterCustomParameter("_int_mineta", min_bnd_cnt);
			holefilling_pts.RegisterCustomParameter("_int_maxeta", max_bnd_cnt);

			holefilling_pts.RegisterCustomParameter("_double2_gm_lowhigh", vmdouble2(g_l, g_h));
		}

	}

	return 0;
}

int update_holefill_points(vmobjects::VmVObjectPrimitive& holefill_pts, const int eta, const float m, const bool show_nc)
{
	cout << "     * eta     : " << eta << endl;
	cout << "     * m       : " << m << endl;
	cout << "     * show_nc : " << show_nc << endl;

	PrimitiveData* prim_data = holefill_pts.GetPrimitiveData();
	if (prim_data == NULL)
	{
		cout << "     no hole-filling points!!" << eta << endl;
		return -1;
	}

	double min_pitch = 0;
	holefill_pts.GetCustomParameter("_double_minpitch", data_type::dtype<double>(), &min_pitch);
	vmdouble2 gm_lh;
	holefill_pts.GetCustomParameter("_double2_gm_lowhigh", data_type::dtype<vmdouble2>(), &gm_lh);

	VmLObject* lobj = holefill_pts.GetBufferObject();
	size_t _t_size;
	int* jet_color_map_ptr;
	int num_ele_jet_color_map = 0;
	lobj->LoadBufferPtr("_vlist_int_jetmap", (void**)&jet_color_map_ptr, _t_size, &num_ele_jet_color_map);

	int* size_hf_clusters;
	int num_ele_numclusters = 0;
	int* connectivity_hf_clusters;
	int num_ele_clusters = 0;
	float* confidence_hf_pts;
	int num_ele_confidences = 0;
	lobj->LoadBufferPtr("_vlist_int_numclusters", (void**)&size_hf_clusters, _t_size, &num_ele_numclusters);
	lobj->LoadBufferPtr("_vlist_int_connectivityclusters", (void**)&connectivity_hf_clusters, _t_size, &num_ele_clusters);
	lobj->LoadBufferPtr("_vlist_float_confidences", (void**)&confidence_hf_pts, _t_size, &num_ele_confidences);

	vmfloat3* clr_vtx = prim_data->GetVerticeDefinition("TEXCOORD0");

	int prev_offset = 0;
	for (int i = 0; i < num_ele_clusters; i++)
	{
		int cnnt = connectivity_hf_clusters[i];
		float _m = cnnt > eta ? m : 1.f;

		int count_vtx = size_hf_clusters[i];
		for (int j = 0; j < count_vtx; j++)
		{
			if (show_nc)
			{
				if (_m > 1.f)
					convert_int_to_float3(0xFF6347, clr_vtx[j + prev_offset]);
				else
					convert_int_to_float3(0x87CEEB, clr_vtx[j + prev_offset]);
			}
			else
			{
				float conf = pow(confidence_hf_pts[j + prev_offset], _m);
				int clr_idx = (int)(conf * (COLORMAP_SIZE - 1) + 0.5f);
				int color = jet_color_map_ptr[clr_idx];
				convert_int_to_float3(color, clr_vtx[j + prev_offset]);
			}
		}
		prev_offset += count_vtx;
	}

	return 0;
}

int localisosurface_points(vmobjects::VmVObjectPrimitive& sm_reliable_pts, vmobjects::VmVObjectPrimitive& sm_holefill_pts,
	const vmobjects::VmVObjectVolume& original_volume, const vmobjects::VmVObjectVolume& filtered_volume,
	const vmobjects::VmVObjectPrimitive& reliable_pts, const vmobjects::VmVObjectPrimitive& holefill_pts, const int e_s, const int eta)
{
	cout << "     * epilon_s  : " << e_s << endl;

	float sample_dist_scale = 1.f;
	VmVObjectVolume* poriginal_volume = (VmVObjectVolume*)&original_volume;
	VmVObjectVolume* pfiltered_volume = (VmVObjectVolume*)&filtered_volume;
	__VolSampleInfo<ushort> vol_original_info = Get_volsample_info<ushort>(poriginal_volume, sample_dist_scale);
	__VolSampleInfo<ushort> vol_filtered_info = Get_volsample_info<ushort>(pfiltered_volume, sample_dist_scale);
	VmVObjectPrimitive* preliable_pts = (VmVObjectPrimitive*)&reliable_pts;
	VmVObjectPrimitive* pholefill_pts = (VmVObjectPrimitive*)&holefill_pts;
	PrimitiveData* prim_reliable_data = preliable_pts->GetPrimitiveData();
	PrimitiveData* prim_holefill_data = pholefill_pts->GetPrimitiveData();

	auto __copy_pobj = [](vmobjects::VmVObjectPrimitive& pobj_src, vmobjects::VmVObjectPrimitive& pobj_dst)
	{
		PrimitiveData* prim_src_data = pobj_src.GetPrimitiveData();
		PrimitiveData prim_dst_data = *prim_src_data;
		prim_dst_data.ClearVertexDefinitionContainer();
		vmfloat3* pos_dst_pts = new vmfloat3[prim_src_data->num_vtx];
		vmfloat3* nrl_dst_pts = new vmfloat3[prim_src_data->num_vtx];
		vmfloat3* clr_dst_pts = new vmfloat3[prim_src_data->num_vtx];
		memcpy(pos_dst_pts, prim_src_data->GetVerticeDefinition("POSITION"), sizeof(vmfloat3) * prim_src_data->num_vtx);
		memcpy(nrl_dst_pts, prim_src_data->GetVerticeDefinition("NORMAL"), sizeof(vmfloat3) * prim_src_data->num_vtx);
		memcpy(clr_dst_pts, prim_src_data->GetVerticeDefinition("TEXCOORD0"), sizeof(vmfloat3) * prim_src_data->num_vtx);
		prim_dst_data.ReplaceOrAddVerticeDefinition("POSITION", pos_dst_pts);
		prim_dst_data.ReplaceOrAddVerticeDefinition("NORMAL", nrl_dst_pts);
		prim_dst_data.ReplaceOrAddVerticeDefinition("TEXCOORD0", clr_dst_pts);
		prim_dst_data.ComputeOrthoBoundingBoxWithCurrentValues();
		pobj_dst.RegisterPrimitiveData(prim_dst_data);
	};
	if(prim_reliable_data)
		__copy_pobj(*preliable_pts, sm_reliable_pts);
	if(prim_holefill_data)
		__copy_pobj(*pholefill_pts, sm_holefill_pts);

	if (e_s <= 0 || (prim_reliable_data == NULL && prim_holefill_data == NULL))
		return 0;

	auto compute_normal = [&vol_filtered_info](const vmfloat3& pos_pt, const vmfloat3& nrl_ref_pt)
	{
		__VolSampleInfo<ushort>& vol_info = vol_filtered_info;
		vmfloat3 pos_sample;
		fTransformPoint(&pos_sample, &pos_pt, &vol_info.mat_ws2vs);
		vmfloat3 grad = __Safe_Gradient_by_Samples(pos_sample, vol_info.vol_size,
			vol_info.vec_grad_dirs[0], vol_info.vec_grad_dirs[1], vol_info.vec_grad_dirs[2],
			vol_info.width_slice, EXB, vol_info.vol_slices);
		float gm = fLengthVector(&grad);
		vmfloat3 nrl = grad / gm;
		if (fDotVector(&nrl_ref_pt, &nrl) < 0) nrl *= -1.f;
		return nrl;
	};

	auto relocate_points = [&e_s, &eta, &preliable_pts, &pholefill_pts, &compute_normal](vmobjects::VmVObjectPrimitive& pobj_src, vmobjects::VmVObjectPrimitive& pobj_dst,
		__VolSampleInfo<ushort>& vol_info)
	{
		PrimitiveData* prim_src_data = pobj_src.GetPrimitiveData();

		VmLObject* lobj = pobj_src.GetBufferObject();
		size_t _t_size;
		int* size_hf_clusters;
		int num_ele_numclusters = 0;
		int* connectivity_hf_clusters;
		int num_ele_clusters = 0;
		float* confidence_hf_pts;
		int num_ele_confidences = 0;
		lobj->LoadBufferPtr("_vlist_int_numclusters", (void**)&size_hf_clusters, _t_size, &num_ele_numclusters);
		lobj->LoadBufferPtr("_vlist_int_connectivityclusters", (void**)&connectivity_hf_clusters, _t_size, &num_ele_clusters);
		lobj->LoadBufferPtr("_vlist_float_confidences", (void**)&confidence_hf_pts, _t_size, &num_ele_confidences);

		vector<bool> valid_src_pts(prim_src_data->num_vtx, true);
		if (size_hf_clusters != NULL)
		{
			int prev_offset = 0;
			for (int i = 0; i < num_ele_clusters; i++)
			{
				int cnnt = connectivity_hf_clusters[i];
				int count_vtx = size_hf_clusters[i];
				if (cnnt > eta)
				{
					for (int j = 0; j < count_vtx; j++)
						valid_src_pts[j + prev_offset] = false;
				}
				prev_offset += count_vtx;
			}
		}

		float* value_r_pts;
		float* value_hf_pts;
		int num_ele_pts = 0;
		preliable_pts->GetBufferObject()->LoadBufferPtr("_vlist_float_values", (void**)&value_r_pts, _t_size, &num_ele_pts);
		pholefill_pts->GetBufferObject()->LoadBufferPtr("_vlist_float_values", (void**)&value_hf_pts, _t_size, &num_ele_pts);

		vector<float> local_isovalue_pts(prim_src_data->num_prims, 0);
		vmfloat3* pos_src_pts = prim_src_data->GetVerticeDefinition("POSITION");
		vmfloat3* nrl_src_pts = prim_src_data->GetVerticeDefinition("NORMAL");

		float eps_s = e_s * vol_info.min_sample_dist;
		float r_sq = eps_s * eps_s;
#pragma omp parallel for num_threads( omp_get_num_procs())
		for (int i = 0; i < (int)prim_src_data->num_prims; i++)
		{
			if (!valid_src_pts[i])
			{
				local_isovalue_pts[i] = 7777777;
				continue;
			}
			const vmfloat3& pos_src = pos_src_pts[i];
			std::vector<std::pair<size_t, float>> r_ret_matches, hf_ret_matches;
			preliable_pts->KDTSearchRadius(pos_src, r_sq, false, r_ret_matches);
			pholefill_pts->KDTSearchRadius(pos_src, r_sq, false, hf_ret_matches);

			double sum_sample_v = 0, cnt = 0;
			auto count_sample = [&sum_sample_v, &cnt](float* buf_value_pts, std::vector<std::pair<size_t, float>>& ret_matches)
			{
				for (int j = 0; j < (int)ret_matches.size(); j++)
				{
					int nb_idx = (int)ret_matches[j].first;
					sum_sample_v += buf_value_pts[nb_idx];
					cnt++;
				}
			};
			count_sample(value_r_pts, r_ret_matches);
			count_sample(value_hf_pts, hf_ret_matches);
			assert(cnt > 0);
			local_isovalue_pts[i] = (float)(sum_sample_v / cnt); 
		}

		PrimitiveData* prim_dst_data = pobj_dst.GetPrimitiveData();
		vmfloat3* pos_dst_pts = prim_dst_data->GetVerticeDefinition("POSITION");

		relocate_to_target_densities(pos_dst_pts, pos_src_pts, nrl_src_pts, &local_isovalue_pts[0], prim_src_data->num_prims, vol_info.min_sample_dist * 0.5f, vol_info);

		// normal update //
		vmfloat3* nrl_dst_pts = prim_dst_data->GetVerticeDefinition("NORMAL");
#pragma omp parallel for num_threads( omp_get_num_procs() ) // 
		for (int i = 0; i < (int)prim_src_data->num_prims; i++)
		{
			nrl_dst_pts[i] = compute_normal(pos_dst_pts[i], nrl_src_pts[i]);
		}
	};

	DWORD time_stage_localis = timeGetTime();
	if (prim_reliable_data)
		relocate_points(*preliable_pts, sm_reliable_pts, vol_original_info);
	if (prim_holefill_data)
		relocate_points(*pholefill_pts, sm_holefill_pts, vol_original_info);
	cout << ">>> Time for local-isosurfacing : " << timeGetTime() - time_stage_localis << " ms" << endl;

	return 0;
}

auto ___compute_tr_solver = [](vmmat44f& mat_os2ws, const int res,
	const vector<__float3>& pos_c_pts, const vector<__float3>& pos_a_pts)
{
	vmfloat3 min_pos(FLT_MAX, FLT_MAX, FLT_MAX), max_pos(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	auto minmax_pos = [&](const vector<__float3>& p)
	{
		for (int i = 0; i < p.size(); i++)
		{
			const __float3& pos = p[i];
			min_pos.x = min(min_pos.x, pos.x);
			min_pos.y = min(min_pos.y, pos.y);
			min_pos.z = min(min_pos.z, pos.z);
			max_pos.x = max(max_pos.x, pos.x);
			max_pos.y = max(max_pos.y, pos.y);
			max_pos.z = max(max_pos.z, pos.z);
		}
	};
	minmax_pos(pos_c_pts);
	minmax_pos(pos_a_pts);

	vmfloat3 pos_cen = (min_pos + max_pos) / 2.f;
	vmfloat3 diff = max_pos - min_pos;
	float scale = max(max(diff.x, diff.y), diff.z) * 1.1f;
	pos_cen -= vmfloat3(scale, scale, scale) / 2.f;
	// [0,1] --> [WS]
	vmmat44f mat_s, mat_t, matd_os2ws;
	fMatrixScaling(&mat_s, &vmfloat3(1.f / scale, 1.f / scale, 1.f / scale));
	fMatrixTranslation(&mat_t, &(-pos_cen));
	matd_os2ws = mat_t * mat_s; // row major
	fMatrixInverse(&matd_os2ws, &matd_os2ws);
	// [0, res-1] = [-0.5, res - 0.5] --> [0, res] --> [0,1]
	float inv_scale = (float)(res);
	fMatrixScaling(&mat_s, &vmfloat3(1.f / inv_scale, 1.f / inv_scale, 1.f / inv_scale));
	fMatrixTranslation(&mat_t, &vmfloat3(0.5f));
	mat_os2ws = mat_t * mat_s * matd_os2ws;
};
enum DataType {
	__DataTypeUNDEFINED = 0,/*!< No Flip */
	__DataTypeBYTE,/*!< 1byte, unsigned char, defined as byte */
	__DataTypeUNSIGNEDSHORT,/*!< 2bytes, unsigned short, defined as ushort*/
	__DataTypeUNSIGNEDINT,/*!< 4bytes, unsigned int, defined as uint*/
	__DataTypeFLOAT,/*!< 4bytes, float, defined as float*/
	__DataTypeDOUBLE,/*!< 8bytes, float, defined as double*/
}; 

auto ___mesh_extract = [](__ProcBuffers<__float3>& trisBuf, const vmmat44f& mat_os2ws, const float iso_value,
	const void** slices, const void* chunk, DataType dt, const __int3 size, const bool inverse_normal)
{
	// with trisBuf.level_set, trisBuf.res
	typedef bool(*LPDLL_MC_MESH)(__ProcBuffers<__float3>& proc_out, const float iso_value, const void** vol_slices, const void* vol_chunk, const __int3& vol_size, const __int3& vol_size_ex, const DataType dt);

	VmHMODULE hMouleLib_MC = NULL;
	LPDLL_MC_MESH lpdll_functionMC = NULL;
	VMLOADLIBRARY(hMouleLib_MC, "vismtv_meshextraction");
	if (hMouleLib_MC == NULL)
	{
		cout << "ERROR import mesh_extract_launcher!!" << endl;
		return true;
	}

	lpdll_functionMC = (LPDLL_MC_MESH)VMGETPROCADDRESS(hMouleLib_MC, "extract_mc_mesh");
	if (lpdll_functionMC == NULL)
	{
		cout << "ERROR import dll function!!" << endl;
		VMFREELIBRARY(hMouleLib_MC);
		return true;
	}

	if (slices == NULL && chunk == NULL && trisBuf.level_set != NULL)
		lpdll_functionMC(trisBuf, trisBuf.iso_value, NULL, trisBuf.level_set, __int3(trisBuf.res, trisBuf.res, trisBuf.res), __int3(), DataType::__DataTypeFLOAT);
	else if (slices != NULL || chunk != NULL)
		lpdll_functionMC(trisBuf, iso_value, slices, chunk, size, __int3(), dt);
	else
		___debugout("___mesh_extract error! ", -1);

#define _cst_vx *(vmfloat3*)&
#define _icst_vx *(__float3*)&
#pragma omp parallel for num_threads( omp_get_num_procs() )
	for (int i = 0; i < (int)trisBuf.num_pts; i++)
	{
		fTransformPoint((vmfloat3*)&trisBuf.pos_pt[i], (vmfloat3*)&trisBuf.pos_pt[i], &mat_os2ws);
		fTransformVector((vmfloat3*)&trisBuf.nrl_pt[i], (vmfloat3*)&trisBuf.nrl_pt[i], &mat_os2ws);
		if (inverse_normal) trisBuf.nrl_pt[i] = _icst_vx(-(_cst_vx trisBuf.nrl_pt[i]));
	}
	VMFREELIBRARY(hMouleLib_MC);

	return true;
};

int processing_final(vmobjects::VmVObjectPrimitive& surface_mesh,
	const vmobjects::VmVObjectPrimitive& sm_reliable_pts, const vmobjects::VmVObjectPrimitive& sm_holefill_pts, const vmobjects::VmVObjectPrimitive& holefill_pts, 
	const int eta, const float m, const int otlev)
{
	cout << "     * eta  : " << eta << endl;
	cout << "     * m    : " << m << endl;

	// poisson //
	const bool use_our_mc = false;
	vector<string> str_cmd_sps_params;
	{
		str_cmd_sps_params.push_back("--depth");
		str_cmd_sps_params.push_back(std::to_string(otlev)); // cube - 1 
		str_cmd_sps_params.push_back("--fullDepth");
		str_cmd_sps_params.push_back(std::to_string(0));
		str_cmd_sps_params.push_back("--cgDepth");
		str_cmd_sps_params.push_back(std::to_string(5));
		str_cmd_sps_params.push_back("--baseDepth");
		str_cmd_sps_params.push_back(std::to_string(0));
		str_cmd_sps_params.push_back("--pointWeight");
		str_cmd_sps_params.push_back(std::to_string(2.0)); // alpha 
		str_cmd_sps_params.push_back("--scale");
		str_cmd_sps_params.push_back(std::to_string(1.1));
		str_cmd_sps_params.push_back("--samplesPerNode");
		str_cmd_sps_params.push_back(std::to_string(2.0));
		str_cmd_sps_params.push_back("--iters");
		str_cmd_sps_params.push_back(std::to_string(8)); // gaussian relaxation
		//if (param.DensityFlag) str_cmd_sps_params.push_back("--density");
		str_cmd_sps_params.push_back("--bType");
		str_cmd_sps_params.push_back(std::to_string(3)); // # 1 : FREE, 2 : DIRICHLET, 3 : NEUMANN
		str_cmd_sps_params.push_back("--confidence");
		str_cmd_sps_params.push_back(std::to_string(1.0));
		str_cmd_sps_params.push_back("--confidenceBias");
		str_cmd_sps_params.push_back(std::to_string(0));
		//if (lvs_normal) str_cmd_sps_params.push_back("--normals");
		//bool lvs_map = _fncontainer.GetParamValue("_bool_LevelsetMap", false);
		//if (lvs_map)
		//{
		if (use_our_mc)
		{
			str_cmd_sps_params.push_back("--levelSetMap");
			str_cmd_sps_params.push_back("--skipContouring");
		}
		//}
	}
	std::vector<char*> cmd_sps_params;
	std::for_each(str_cmd_sps_params.begin(), str_cmd_sps_params.end(), [&cmd_sps_params](string &str) {
		cmd_sps_params.push_back((char*)str.c_str()); }
	);
	std::vector<void*> additional_params; // deprecated

	vector<__float3> pos_r_pts, nrl_r_pts;
	vector<__float3> pos_hf_pts, nrl_hf_pts;
	{
		VmVObjectPrimitive* psm_reliable_pts = (VmVObjectPrimitive*)&sm_reliable_pts;
		VmVObjectPrimitive* psm_holefill_pts = (VmVObjectPrimitive*)&sm_holefill_pts;
		PrimitiveData* prim_sm_reliable_data = psm_reliable_pts->GetPrimitiveData();
		if (prim_sm_reliable_data == NULL) 
			return -1;
		pos_r_pts.assign(prim_sm_reliable_data->num_vtx, __float3());
		nrl_r_pts.assign(prim_sm_reliable_data->num_vtx, __float3());
		PrimitiveData* prim_sm_holefill_data = psm_holefill_pts->GetPrimitiveData();
		memcpy(&pos_r_pts[0], prim_sm_reliable_data->GetVerticeDefinition("POSITION"), prim_sm_reliable_data->num_vtx * sizeof(vmfloat3));
		memcpy(&nrl_r_pts[0], prim_sm_reliable_data->GetVerticeDefinition("NORMAL"), prim_sm_reliable_data->num_vtx * sizeof(vmfloat3));

		if (prim_sm_holefill_data)
		{
			pos_hf_pts.assign(prim_sm_holefill_data->num_vtx, __float3());
			nrl_hf_pts.assign(prim_sm_holefill_data->num_vtx, __float3());
			memcpy(&pos_hf_pts[0], prim_sm_holefill_data->GetVerticeDefinition("POSITION"), prim_sm_holefill_data->num_vtx * sizeof(vmfloat3));
			memcpy(&nrl_hf_pts[0], prim_sm_holefill_data->GetVerticeDefinition("NORMAL"), prim_sm_holefill_data->num_vtx * sizeof(vmfloat3));
			VmLObject* lobj = ((VmVObjectPrimitive*)&holefill_pts)->GetBufferObject();
			size_t _t_size;

			int* size_hf_clusters;
			int num_ele_numclusters = 0;
			int* connectivity_hf_clusters;
			int num_ele_clusters = 0;
			float* confidence_hf_pts;
			int num_ele_confidences = 0;
			lobj->LoadBufferPtr("_vlist_int_numclusters", (void**)&size_hf_clusters, _t_size, &num_ele_numclusters);
			lobj->LoadBufferPtr("_vlist_int_connectivityclusters", (void**)&connectivity_hf_clusters, _t_size, &num_ele_clusters);
			lobj->LoadBufferPtr("_vlist_float_confidences", (void**)&confidence_hf_pts, _t_size, &num_ele_confidences);

			int prev_offset = 0;
			for (int i = 0; i < num_ele_clusters; i++)
			{
				int cnnt = connectivity_hf_clusters[i];
				int count_vtx = size_hf_clusters[i];
				if (cnnt > eta)
				{
					for (int j = 0; j < count_vtx; j++)
					{
						float conf = confidence_hf_pts[j + prev_offset];
						vmfloat3 n = *(vmfloat3*)&nrl_hf_pts[j + prev_offset] * pow(conf, m);
						nrl_hf_pts[j + prev_offset] = *(__float3*)&n;
					}
				}
				else
				{
					for (int j = 0; j < count_vtx; j++)
					{
						float conf = confidence_hf_pts[j + prev_offset];
						vmfloat3 n = *(vmfloat3*)&nrl_hf_pts[j + prev_offset] * conf;
						nrl_hf_pts[j + prev_offset] = *(__float3*)&n;
					}
				}
				prev_offset += count_vtx;
			}
		}

		// test //
		psm_reliable_pts->GetBufferObject()->ReplaceOrAddBufferPtr("_vlist_FLOAT3_PosField", &pos_hf_pts[0], (int)pos_hf_pts.size(), sizeof(vmfloat3));
		psm_reliable_pts->GetBufferObject()->ReplaceOrAddBufferPtr("_vlist_FLOAT3_VecField", &nrl_hf_pts[0], (int)nrl_hf_pts.size(), sizeof(vmfloat3));
	}

	// test //
	//auto ___test = [](vector<__float3>& n_pts)
	//{
	//	float __minn(10000000), __maxn(0);
	//	int n_max = 0, n_min = 0;
	//	for (int i = 0; i < (int)n_pts.size(); i++)
	//	{
	//		float l = fLengthVector((vmfloat3*)&n_pts[i]);
	//		__minn = min(__minn, l);
	//		__maxn = max(__maxn, l);
	//		if (l == 0) n_min++;
	//		if (l > 1) n_max++;
	//	}
	//
	//	cout << "*** __minn : " << __minn << endl;
	//	cout << "*** __maxn : " << __maxn << endl;
	//	cout << "***  n_min : " << n_min << endl;
	//	cout << "***  n_max : " << n_max << endl;
	//};
	//___test(nrl_r_pts);
	//___test(nrl_hf_pts);

	__ProcBuffers<__float3> trisbuf;
	ScreenedPoissonSurface3D(&trisbuf, pos_r_pts, nrl_r_pts, pos_hf_pts, nrl_hf_pts, cmd_sps_params, additional_params);
	if(use_our_mc)
	{
		vmmat44f mat_os2ws;
		___compute_tr_solver(mat_os2ws, trisbuf.res, pos_r_pts, pos_hf_pts);

		bool use_smooth_filter = true;
		if (use_smooth_filter)
		{
			
			ushort** vol_ushort = new ushort*[trisbuf.res];
			for (int i = 0; i < trisbuf.res; i++) vol_ushort[i] = new ushort[trisbuf.res * trisbuf.res];
			float min_v = FLT_MAX, max_v = -FLT_MAX;
			for (int i = 0; i < trisbuf.res; i++)
				for(int j = 0; j < trisbuf.res*trisbuf.res; j++)
				{
					float v = trisbuf.level_set[i][j];
					min_v = min(min_v, v);
					max_v = max(max_v, v);
				}
			for (int i = 0; i < trisbuf.res; i++)
				for (int j = 0; j < trisbuf.res*trisbuf.res; j++)
				{
					float v = trisbuf.level_set[i][j]; 
					vol_ushort[i][j] = (ushort)((v - min_v) / (max_v - min_v) * 65535.f);
				}
			
			typedef bool(*LPDLL_CALL_MoFiInitializeDx11)();
			typedef bool(*LPDLL_CALL_MoFiDeinitializeDx11)();
			typedef bool(*LPDLL_CALL_MorphGaussianBlur3D)(ushort** ppusVolumeIn, ushort** ppusVolumeOut, vmint3 i3SizeVolume, vmint2 i2OffsetZ, int iKernelSizeHalf, float fSigma, LocalProgress* _progress);
			LPDLL_CALL_MoFiInitializeDx11 lpdll_init = LoadDLL<LPDLL_CALL_MoFiInitializeDx11>("vismtv_morphfilters", "MoFiInitializeDx11");
			LPDLL_CALL_MoFiDeinitializeDx11 lpdll_deinit = LoadDLL<LPDLL_CALL_MoFiDeinitializeDx11>("vismtv_morphfilters", "MoFiDeinitializeDx11");
			LPDLL_CALL_MorphGaussianBlur3D lpdll_gaussianblur = LoadDLL<LPDLL_CALL_MorphGaussianBlur3D>("vismtv_morphfilters", "MorphGaussianBlur3D");
			if (lpdll_init == NULL || lpdll_deinit == NULL || lpdll_gaussianblur == NULL)
			{
				cout << "fail to load morphological filter module" << endl;
				return -1;
			}

			lpdll_init();
			LocalProgress _progress;
			lpdll_gaussianblur(vol_ushort, vol_ushort, vmint3(trisbuf.res, trisbuf.res, trisbuf.res), vmint2(), 3, 1.4f, &_progress);
			___mesh_extract(trisbuf, mat_os2ws, ((trisbuf.iso_value - min_v) / (max_v - min_v) * 65535.f), (const void**)vol_ushort, NULL,
				DataType::__DataTypeUNSIGNEDSHORT, __int3(trisbuf.res, trisbuf.res, trisbuf.res), false);
			VMSAFE_DELETE2DARRAY(vol_ushort, trisbuf.res);
			lpdll_deinit();
			ReleaseDlls();
		}
		else
		{
			___mesh_extract(trisbuf, mat_os2ws, trisbuf.iso_value, (const void**)trisbuf.level_set, NULL,
				DataType::__DataTypeFLOAT, __int3(trisbuf.res, trisbuf.res, trisbuf.res), false);
		}
	}

	auto register_trisBuf_to_pobj = [](__ProcBuffers<__float3>& trisbuf, vmobjects::VmVObjectPrimitive& pobj)
	{
		PrimitiveData mesh_data;
		mesh_data.is_ccw = true;
		mesh_data.is_stripe = false;
		mesh_data.check_redundancy = true;
		mesh_data.ptype = PrimitiveTypeTRIANGLE;
		mesh_data.idx_stride = 3;
		mesh_data.num_vtx = trisbuf.num_pts;
		mesh_data.num_prims = trisbuf.num_primitives;
		mesh_data.num_vidx = trisbuf.num_primitives * 3;
		mesh_data.vidx_buffer = trisbuf.index_buffer;

		vmfloat3* pos_pts = (vmfloat3*)trisbuf.pos_pt;
		mesh_data.ReplaceOrAddVerticeDefinition("POSITION", pos_pts);
		// compute normal
		{
			vmfloat3* nrl_pts = new vmfloat3[mesh_data.num_vtx];
			ZeroMemory(nrl_pts, sizeof(vmfloat3) * mesh_data.num_vtx);
			for (uint i = 0; i < mesh_data.num_prims; i++)
			{
				int iIndex0, iIndex1, iIndex2;
				iIndex0 = mesh_data.vidx_buffer[3 * i + 0];
				iIndex1 = mesh_data.vidx_buffer[3 * i + 1];
				iIndex2 = mesh_data.vidx_buffer[3 * i + 2];

				const vmfloat3& f3Pos0 = pos_pts[iIndex0];
				const vmfloat3& f3Pos1 = pos_pts[iIndex1];
				const vmfloat3& f3Pos2 = pos_pts[iIndex2];

				vmfloat3 f3Vec01 = f3Pos1 - f3Pos0;
				vmfloat3 f3Vec02 = f3Pos2 - f3Pos0;

				vmfloat3 f3VecTriangleNormal;
				fCrossDotVector(&f3VecTriangleNormal, &f3Vec01, &f3Vec02);
				fNormalizeVector(&f3VecTriangleNormal, &f3VecTriangleNormal);

				nrl_pts[iIndex0] += f3VecTriangleNormal;
				nrl_pts[iIndex1] += f3VecTriangleNormal;
				nrl_pts[iIndex2] += f3VecTriangleNormal;
			}
			for (uint i = 0; i < mesh_data.num_vtx; i++)
			{
				fNormalizeVector(&nrl_pts[i], &nrl_pts[i]);
			}
			mesh_data.ReplaceOrAddVerticeDefinition("NORMAL", nrl_pts);
		}

		mesh_data.ComputeOrthoBoundingBoxWithCurrentValues();
		pobj.RegisterPrimitiveData(mesh_data);
		pobj.RegisterCustomParameter("_bool_ApplyShadingFactors", true);

		trisbuf.pos_pt = NULL;
		trisbuf.nrl_pt = NULL;
		trisbuf.index_buffer = NULL;
	};

	register_trisBuf_to_pobj(trisbuf, surface_mesh);

	return 0;
}

bool function_launcher(const std::vector<vmobjects::VmObject*>& io_objs, const std::map<std::string, std::any>& parameters)
{
	using namespace std;
	std::map<std::string, std::any>& _parameters = (std::map<std::string, std::any>&)parameters;
	const std::string& function_name = any_cast<string>(_parameters["function_name"]);
	if (function_name == "load_foreground")
	{
		if (io_objs.size() != 2) return false;
		VmVObjectVolume& main_volume = *(VmVObjectVolume*)io_objs[0];
		VmVObjectVolume& fore_volume = *(VmVObjectVolume*)io_objs[1];

		load_foreground(any_cast<string>(_parameters["file_name"]), main_volume, fore_volume,
			any_cast<int>(_parameters["iso_value"]), any_cast<int>(_parameters["downsacled_factor"]),
			any_cast<bool>(_parameters["flip_x"]), any_cast<bool>(_parameters["flip_y"]), any_cast<bool>(_parameters["flip_z"]));
	}
	else if (function_name == "processing_stage1")
	{
		if (io_objs.size() != 4) return false;
		VmVObjectPrimitive& cand_pts = *(VmVObjectPrimitive*)io_objs[0];
		VmVObjectVolume& filtered_volume = *(VmVObjectVolume*)io_objs[1];
		VmVObjectVolume& main_volume = *(VmVObjectVolume*)io_objs[2];
		VmVObjectVolume& fore_volume = *(VmVObjectVolume*)io_objs[3];

		processing_stage1(cand_pts, filtered_volume, main_volume, fore_volume,
			any_cast<int>(_parameters["t_in"]), any_cast<int>(_parameters["t_out"]),
			any_cast<float>(_parameters["min_sample_v"]), any_cast<float>(_parameters["max_sample_v"]),
			any_cast<int>(_parameters["gaussian_kernel_width"]), any_cast<float>(_parameters["simplify_grid_length"]), any_cast<float>(_parameters["geometric_complexity_kernel_ratio"]));
	}
	else if (function_name == "update_strongweak_points")
	{
		if (io_objs.size() != 1) return false;
		VmVObjectPrimitive& cand_pts = *(VmVObjectPrimitive*)io_objs[0];

		update_strongweak_points(cand_pts, any_cast<float>(_parameters["g_h"]), any_cast<float>(_parameters["mu_u"]));
	}
	else if (function_name == "processing_stage2")
	{
		if (io_objs.size() != 4) return false;
		VmVObjectPrimitive& reliable_pts = *(VmVObjectPrimitive*)io_objs[0];
		VmVObjectPrimitive& hf_pts = *(VmVObjectPrimitive*)io_objs[1];
		VmVObjectVolume& sample_volume = *(VmVObjectVolume*)io_objs[2];
		VmVObjectPrimitive& cand_pts = *(VmVObjectPrimitive*)io_objs[3];

		processing_stage2(reliable_pts, hf_pts, sample_volume, cand_pts, any_cast<float>(_parameters["g_h"]), any_cast<float>(_parameters["mu_u"]),
			any_cast<float>(_parameters["epsilon"]), any_cast<int>(_parameters["connectivity_criterion"]), any_cast<float>(_parameters["epsilon_b"]), any_cast<float>(_parameters["angle_criterion"]));
	}
	else if (function_name == "update_holefill_points")
	{
		if (io_objs.size() != 1) return false;
		VmVObjectPrimitive& hf_pts = *(VmVObjectPrimitive*)io_objs[0];

		update_holefill_points(hf_pts, any_cast<int>(_parameters["eta"]), any_cast<float>(_parameters["m"]), any_cast<bool>(_parameters["show_nc"]));
	}
	else if (function_name == "localisosurface_points")
	{
		if (io_objs.size() != 6) return false;
		VmVObjectPrimitive& ls_reliable_pts = *(VmVObjectPrimitive*)io_objs[0];
		VmVObjectPrimitive& ls_hf_pts = *(VmVObjectPrimitive*)io_objs[1];
		VmVObjectVolume& ct_vol = *(VmVObjectVolume*)io_objs[2];
		VmVObjectVolume& filtered_vol = *(VmVObjectVolume*)io_objs[3];
		VmVObjectPrimitive& reliable_pts = *(VmVObjectPrimitive*)io_objs[4];
		VmVObjectPrimitive& hf_pts = *(VmVObjectPrimitive*)io_objs[5];

		localisosurface_points(ls_reliable_pts, ls_hf_pts, ct_vol, filtered_vol, reliable_pts, hf_pts, any_cast<int>(_parameters["epsilon_s"]), any_cast<int>(_parameters["eta"]));
	}
	else if (function_name == "processing_final")
	{
		if (io_objs.size() != 4) return false;
		VmVObjectPrimitive& final_mesh = *(VmVObjectPrimitive*)io_objs[0];
		VmVObjectPrimitive& ls_reliable_pts = *(VmVObjectPrimitive*)io_objs[1];
		VmVObjectPrimitive& ls_hf_pts = *(VmVObjectPrimitive*)io_objs[2];
		VmVObjectPrimitive& hf_pts = *(VmVObjectPrimitive*)io_objs[3];

		processing_final(final_mesh, ls_reliable_pts, ls_hf_pts, hf_pts, any_cast<int>(_parameters["eta"]), any_cast<float>(_parameters["m"]), any_cast<int>(_parameters["otlev"]));
	}
	return false;
}