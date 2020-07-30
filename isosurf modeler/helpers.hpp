#pragma once
#include "VimCommon.h"
#pragma warning (disable:4756)

#include <iostream>
#include <queue>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#ifndef _OPENMP
int omp_get_num_procs(void) { return 1; }
int omp_get_thread_num(void) { return 0; }
#endif // _OPENMP

#define _f3_(x, y, z) (float)(x), (float)(y), (float)(z)
#define _i3_(x, y, z) (int)(x), (int)(y), (int)(z)
#define _fp_ float*

#define __v3minmax(__INOUT, __IN, __OP) {\
	__INOUT.x = __OP(__INOUT.x, __IN.x);\
	__INOUT.y = __OP(__INOUT.y, __IN.y);\
	__INOUT.z = __OP(__INOUT.z, __IN.z);\
};

#define EXB 2

using namespace vmobjects;
using namespace vmmath;
using namespace std;

template <typename T>
inline T __ReadVoxel(const vmint3& idx, const int width_slice, const int exb, T** vol_slices)
{
	return vol_slices[idx.z + exb][idx.x + exb + (idx.y + exb) * width_slice];
};

template <typename T>
inline void __WriteVoxel(T v, const vmint3& idx, const int width_slice, const int exb, T** vol_slices)
{
	vol_slices[idx.z + exb][idx.x + exb + (idx.y + exb) * width_slice] = v;
};

inline vmint3 __MultInt3(const vmint3* pi3_0, const vmint3* pi3_1)
{
	return vmint3(pi3_0->x * pi3_1->x, pi3_0->y * pi3_1->y, pi3_0->z * pi3_1->z);
}

inline bool __SafeCheck(const vmint3& idx, const vmint3& vol_size)
{
	vmint3 max_dir = idx - (vol_size - vmint3(1, 1, 1));
	vmint3 mult_dot = __MultInt3(&max_dir, &idx);
	return !(mult_dot.x > 0 || mult_dot.y > 0 || mult_dot.z > 0);
}

template <typename T>
inline T __Safe__ReadVoxel(const vmint3& idx, const vmint3& vol_size, const int width_slice, const int exb, T** vol_slices, const T bnd_v = 0)
{
	return __SafeCheck(idx, vol_size) ? __ReadVoxel<T>(idx, width_slice, exb, vol_slices) : (T)bnd_v;
};

inline float __TrilinearInterpolation(float v_0, float v_1, float v_2, float v_3, float v_4, float v_5, float v_6, float v_7,
	const vmfloat3& ratio)
{
	float v01 = v_0 * (1.f - ratio.x) + v_1 * ratio.x;
	float v23 = v_2 * (1.f - ratio.x) + v_3 * ratio.x;
	float v0123 = v01 * (1.f - ratio.y) + v23 * ratio.y;
	float v45 = v_4 * (1.f - ratio.x) + v_5 * ratio.x;
	float v67 = v_6 * (1.f - ratio.x) + v_7 * ratio.x;
	float v4567 = v45 * (1.f - ratio.y) + v67 * ratio.y;
	return v0123 * (1.f - ratio.z) + v4567 * ratio.z;
}

template <typename T>
inline float __Safe_TrilinearSample(const vmfloat3& pos_sample, const vmint3& vol_size, const int width_slice, const int exb, T** vol_slices, const T bnd_v = 0)
{
	vmfloat3 __pos_sample = pos_sample + vmfloat3(1000.f, 1000.f, 1000.f); // SAFE //
	vmint3 idx_sample = vmint3((int)__pos_sample.x, (int)__pos_sample.y, (int)__pos_sample.z);
	vmfloat3 ratio = vmfloat3((__pos_sample.x) - (float)(idx_sample.x), (__pos_sample.y) - (float)(idx_sample.y), (__pos_sample.z) - (float)(idx_sample.z));
	idx_sample -= vmint3(1000, 1000, 1000);

	float v0, v1, v2, v3, v4, v5, v6, v7;
	v0 = (float)__Safe__ReadVoxel<T>(idx_sample + vmint3(0, 0, 0), vol_size, width_slice, exb, vol_slices, bnd_v);
	v1 = (float)__Safe__ReadVoxel<T>(idx_sample + vmint3(1, 0, 0), vol_size, width_slice, exb, vol_slices, bnd_v);
	v2 = (float)__Safe__ReadVoxel<T>(idx_sample + vmint3(0, 1, 0), vol_size, width_slice, exb, vol_slices, bnd_v);
	v3 = (float)__Safe__ReadVoxel<T>(idx_sample + vmint3(1, 1, 0), vol_size, width_slice, exb, vol_slices, bnd_v);
	v4 = (float)__Safe__ReadVoxel<T>(idx_sample + vmint3(0, 0, 1), vol_size, width_slice, exb, vol_slices, bnd_v);
	v5 = (float)__Safe__ReadVoxel<T>(idx_sample + vmint3(1, 0, 1), vol_size, width_slice, exb, vol_slices, bnd_v);
	v6 = (float)__Safe__ReadVoxel<T>(idx_sample + vmint3(0, 1, 1), vol_size, width_slice, exb, vol_slices, bnd_v);
	v7 = (float)__Safe__ReadVoxel<T>(idx_sample + vmint3(1, 1, 1), vol_size, width_slice, exb, vol_slices, bnd_v);
	return __TrilinearInterpolation(v0, v1, v2, v3, v4, v5, v6, v7, ratio);
};

template <typename T>
inline vmfloat3 __Safe_Gradient_by_Samples(const vmfloat3& pos_sample, const vmint3& vol_size, const vmfloat3& dir_x, const vmfloat3& dir_y, const vmfloat3& dir_z, const int width_slice, const int exb, T** vol_slices)
{
	float v_XL = __Safe_TrilinearSample<T>(pos_sample - dir_x, vol_size, width_slice, exb, vol_slices);
	float v_XR = __Safe_TrilinearSample<T>(pos_sample + dir_x, vol_size, width_slice, exb, vol_slices);
	float v_YL = __Safe_TrilinearSample<T>(pos_sample - dir_y, vol_size, width_slice, exb, vol_slices);
	float v_YR = __Safe_TrilinearSample<T>(pos_sample + dir_y, vol_size, width_slice, exb, vol_slices);
	float v_ZL = __Safe_TrilinearSample<T>(pos_sample - dir_z, vol_size, width_slice, exb, vol_slices);
	float v_ZR = __Safe_TrilinearSample<T>(pos_sample + dir_z, vol_size, width_slice, exb, vol_slices);
	float g_x = v_XR - v_XL;
	float g_y = v_YR - v_YL;
	float g_z = v_ZR - v_ZL;
	return vmfloat3(g_x, g_y, g_z) * 0.5f;
}

#define SGRAD(p, vinfo) __Safe_Gradient_by_Samples(p, vinfo.vol_size, \
	vinfo.vec_grad_dirs[0], vinfo.vec_grad_dirs[1], vinfo.vec_grad_dirs[2], \
		vinfo.width_slice, EXB, vinfo.vol_slices)

template <typename T>
void ___debugout(std::string str, T v)
{
	std::cout << str.c_str() << std::to_string(v).c_str() << std::endl;
}

template <typename T>
struct __VolSampleInfo
{
	//ushort** vol_slices_filtered;
	T** vol_slices;
	vmint3 vol_size;
	int width_slice;
	vmmat44f mat_ws2vs;
	vmmat44f mat_vs2ws;
	vmfloat3 vec_grad_dirs[3];
	float min_sample_dist;

	__VolSampleInfo() : vol_slices(NULL), //vol_slices_origin(NULL),
		vol_size(vmint3()), width_slice(0), mat_ws2vs(NULL), mat_vs2ws(NULL) {};

	__VolSampleInfo(T** _vol_slices, vmint3 _vol_size,
		int _width_slice, vmmat44f _mat_ws2vs, vmmat44f _mat_vs2ws,
		vmfloat3 _vec_grad_dirs[3], float _min_sample_dist) : vol_slices(_vol_slices),
		vol_size(_vol_size), width_slice(_width_slice), mat_ws2vs(_mat_ws2vs), mat_vs2ws(_mat_vs2ws), min_sample_dist(_min_sample_dist)
	{
		//vol_slices_filtered = _vol_slices;
		vec_grad_dirs[0] = _vec_grad_dirs[0];
		vec_grad_dirs[1] = _vec_grad_dirs[1];
		vec_grad_dirs[2] = _vec_grad_dirs[2];
	}
};

template <typename T>
__VolSampleInfo<T> Get_volsample_info(VmVObjectVolume* pCVolume, const double sample_dist_scale)
{
	if (pCVolume == NULL) return __VolSampleInfo<T>();

	vmfloat3 vecVoxelGradDirs[3];// = { vmfloat3(1, 0, 0), vmfloat3(0, 1, 0), vmfloat3(0, 0, 1) };
	VolumeData* volArchive = pCVolume->GetVolumeData();
	float min_dist_sample = (float)min(min(volArchive->vox_pitch.x, volArchive->vox_pitch.y), volArchive->vox_pitch.z);
	fTransformVector(&vecVoxelGradDirs[0], &vmfloat3(min_dist_sample, 0, 0), &pCVolume->GetMatrixWS2OSf());
	fTransformVector(&vecVoxelGradDirs[1], &vmfloat3(0, min_dist_sample, 0), &pCVolume->GetMatrixWS2OSf());
	fTransformVector(&vecVoxelGradDirs[2], &vmfloat3(0, 0, min_dist_sample), &pCVolume->GetMatrixWS2OSf());

	vmfloat3 vecVoxelGradDirs_unit[3] = {
		vecVoxelGradDirs[0] * (float)sample_dist_scale,
		vecVoxelGradDirs[1] * (float)sample_dist_scale,
		vecVoxelGradDirs[2] * (float)sample_dist_scale };

	VolumeData *volArchiveSample = pCVolume->GetVolumeData();
	vmint3 volSize = volArchiveSample->vol_size;
	vmint3 volSizeEx = volArchiveSample->bnd_size;
	int widthSamplePitch = volSize.x + volSizeEx.x * 2;

	return __VolSampleInfo<T>((T**)volArchiveSample->vol_slices,
		volSize, widthSamplePitch, pCVolume->GetMatrixWS2OSf(), pCVolume->GetMatrixOS2WSf(), vecVoxelGradDirs_unit, min_dist_sample);
};

#include "../nanoflann.hpp"
template <typename T, typename TT>
struct PointCloud
{
	//public:
	const TT* pts;
	const size_t num_pts;
	PointCloud(const TT* _pts, const size_t _num_pts) : pts(_pts), num_pts(_num_pts) { }
	PointCloud(const vector<TT>& vtr_pts) : pts(&vtr_pts[0]), num_pts(vtr_pts.size()) { }

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return num_pts; }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline T kdtree_distance(const T *p1, const size_t idx_p2, size_t) const
	{
		const T d0 = p1[0] - pts[idx_p2].x;
		const T d1 = p1[1] - pts[idx_p2].y;
		const T d2 = p1[2] - pts[idx_p2].z;
		return d0 * d0 + d1 * d1 + d2 * d2;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim == 0) return pts[idx].x;
		else if (dim == 1) return pts[idx].y;
		else return pts[idx].z;
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<float, PointCloud<float, vmfloat3> >,
	PointCloud<float, vmfloat3>,
	3 // dim 
> kd_tree_t;

void make_band(__VolSampleInfo<char>& band_info, const int v_low, const int v_high, const __VolSampleInfo<ushort>& vol_info_outside, const __VolSampleInfo<ushort>& vol_info_inside)
{
	const int clip_bnd = 2;
#pragma omp parallel for num_threads( omp_get_num_procs() )
	for (int z = 0; z < band_info.vol_size.z; z++)
		for (int y = 0; y < band_info.vol_size.y; y++)
			for (int x = 0; x < band_info.vol_size.x; x++)
			{
				if (x < clip_bnd || y < clip_bnd || z < clip_bnd
					|| x >= band_info.vol_size.x - clip_bnd
					|| y >= band_info.vol_size.y - clip_bnd
					|| z >= band_info.vol_size.z - clip_bnd)
					continue;
				int v_out = __ReadVoxel(vmint3(x, y, z), vol_info_outside.width_slice, EXB, vol_info_outside.vol_slices);
				int v_in = __ReadVoxel(vmint3(x, y, z), vol_info_inside.width_slice, EXB, vol_info_inside.vol_slices);
				if (v_in > v_high) __WriteVoxel((char)-1, vmint3(x, y, z), band_info.width_slice, EXB, band_info.vol_slices);
				else if (v_out < v_low) __WriteVoxel((char)1, vmint3(x, y, z), band_info.width_slice, EXB, band_info.vol_slices);
			}
}

template <typename T>
void get_isosurface_vtx(set<tuple<float, float, float>>& voxel_pts_set,
	const float isovalue, const uint sampleoffset, const __VolSampleInfo<T>& vol_info)
{
	auto __read_v = [&](vmint3& idx)
	{
		idx.x = min(idx.x, vol_info.vol_size.x - 1);
		idx.y = min(idx.y, vol_info.vol_size.y - 1);
		idx.z = min(idx.z, vol_info.vol_size.z - 1);
		return (float)__ReadVoxel(idx, vol_info.width_slice, EXB, vol_info.vol_slices);
	};

	const float offset = (float)sampleoffset;
	for (int z = 0; z <= vol_info.vol_size.z; z += sampleoffset)
	{
		for (int y = 0; y <= vol_info.vol_size.y; y += sampleoffset)
		{
			for (int x = 0; x <= vol_info.vol_size.x; x += sampleoffset)
			{
				// sample at 8 nodes of a cell
				float samplevalues[4] = {
					__read_v(vmint3(x + 0, y + 0, z + 0)),
					__read_v(vmint3(x + sampleoffset, y + 0, z + 0)),
					__read_v(vmint3(x + 0, y + sampleoffset, z + 0)),
					__read_v(vmint3(x + 0, y + 0, z + sampleoffset)),
				};

				bool isOriSmall = samplevalues[0] < isovalue;
				bool isEdgeXSmall = samplevalues[1] < isovalue;
				bool isEdgeYSmall = samplevalues[2] < isovalue;
				bool isEdgeZSmall = samplevalues[3] < isovalue;

				if (isOriSmall != isEdgeXSmall)
				{
					float ratio = (float)(isovalue - samplevalues[0]) / (float)(samplevalues[1] - samplevalues[0]);
					vmfloat3 pos = vmfloat3((float)x + ratio * offset, (float)y, (float)z);
					voxel_pts_set.insert(tuple<float, float, float>(pos.x, pos.y, pos.z));
				}
				if (isOriSmall != isEdgeYSmall)
				{
					float ratio = (float)(isovalue - samplevalues[0]) / (float)(samplevalues[2] - samplevalues[0]);
					vmfloat3 pos = vmfloat3((float)x, (float)y + ratio * offset, (float)z);
					voxel_pts_set.insert(tuple<float, float, float>(pos.x, pos.y, pos.z));
				}
				if (isOriSmall != isEdgeZSmall)
				{
					float ratio = (float)(isovalue - samplevalues[0]) / (float)(samplevalues[3] - samplevalues[0]);
					vmfloat3 pos = vmfloat3((float)x, (float)y, (float)z + ratio * offset);
					voxel_pts_set.insert(tuple<float, float, float>(pos.x, pos.y, pos.z));
				}
			} // for x
		} // for y
	} //for z
}

void simplify_points_ugrid(vector<vmfloat3>& pos_simplified_pts, vmfloat3& aabb_diff, const vector<vmfloat3>& pos_pts, const float grid_length)
{
	uint num_pts = (uint)pos_pts.size();
	vmfloat3 aabb_min(FLT_MAX), aabb_max(-FLT_MAX);
	for (uint i = 0; i < num_pts; i++)
	{
		const vmfloat3& pos_pt = pos_pts[i];
		__v3minmax(aabb_min, pos_pt, min);
		__v3minmax(aabb_max, pos_pt, max);
	}
	aabb_diff = aabb_max - aabb_min;
	vmint3 aabb_size = vmint3(aabb_diff / grid_length) + vmint3(1);
	uint* index_map = new uint[aabb_size.x * aabb_size.y * aabb_size.z];
	memset(index_map, 0, sizeof(uint) * aabb_size.x * aabb_size.y * aabb_size.z);

	uint count = 0;
	for (uint i = 0; i < num_pts; i++)
	{
		const vmfloat3& pos_pt = pos_pts[i];
		vmfloat3 pos_cell = (pos_pt - aabb_min) / grid_length;
		vmint3 idx_cell = pos_cell;
		uint addr = (uint)idx_cell.x + (uint)(idx_cell.y * aabb_size.x) + (uint)idx_cell.z * (uint)(aabb_size.x * aabb_size.y);
		uint prev_idx = index_map[addr];
		if (prev_idx == 0)
		{
			index_map[addr] = i + 1;
			count++;
		}
		else
		{
			const vmfloat3& pos_prev_pt = pos_pts[prev_idx];
			vmfloat3 pos_prev_cell = (pos_prev_pt - aabb_min) / grid_length;
			vmint3 idx_prev_cell = pos_prev_cell;

			vmfloat3 cell_pos = pos_cell - (vmfloat3)idx_cell - vmfloat3(0.5f);
			vmfloat3 prev_cell_pos = pos_prev_cell - (vmfloat3)idx_prev_cell - vmfloat3(0.5f);
			if (fLengthVectorSq(&cell_pos) < fLengthVectorSq(&prev_cell_pos))
				index_map[addr] = i + 1;
		}
	}

	pos_simplified_pts.assign(count, vmfloat3());
	uint count_idx = 0;
	for (uint i = 0; i < (uint)aabb_size.x * (uint)aabb_size.y * (uint)aabb_size.z; i++)
	{
		uint idx = index_map[i];
		if (idx > 0)
			pos_simplified_pts[count_idx++] = pos_pts[idx - 1];
	}

	VMSAFE_DELETEARRAY(index_map);
}

int OtsuThresholdValue(ullong* pullHistogramValues, int histo_size, int begin_idx, int end_idx)
{
	double total_elements = 0;
	double sum1 = 0;

	//for (int i = begin_idx; i <= end_idx; i++)
	for (int i = 0; i < histo_size; i++)
	{
		total_elements += (double)pullHistogramValues[i];
		sum1 += (double)i * (double)pullHistogramValues[i];
	}

	double sumB = 0;
	double wB = 0;
	double _maximum = 0;
	begin_idx = min(begin_idx, histo_size - 1);
	end_idx = min(end_idx, histo_size - 1);

	int otsuThreshold = 0;
	for (int i = begin_idx; i <= end_idx; i++)
	{
		wB += (double)pullHistogramValues[i];
		double wF = total_elements - wB;
		if (wB == 0 || wF == 0)
			continue;

		sumB += (double)i * (double)pullHistogramValues[i];
		double mF = (sum1 - sumB) / wF;
		double btn = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
		if (btn >= _maximum)
		{
			otsuThreshold = i;
			_maximum = btn;
		}
	}

	return otsuThreshold;
}

void raytraversal_gradmax(std::vector<vmfloat3>& lmax_pts,
	const std::vector<vmfloat3>& pos_pts, const std::vector<vmfloat3>& dir_pts,
	const int num_maxrc_steps,
	const float dir_scale, /*const bool use_quadinterpolation, const bool use_consistent_dir,*/
	const float min_v, const float max_v,
	const __VolSampleInfo<char>& band_info,
	const __VolSampleInfo<ushort>& vol_info)
{
	auto __gradient = [&](const vmfloat3& pos_vs)
	{
		return SGRAD(pos_vs, vol_info);
		//return -__Safe_TrilinearSample_T(pos_vs, mask_Info.vol_size, mask_Info.width_slice, 0, mask_Info.grad_slices, vmfloat3());
	};

	int num_pts = (int)pos_pts.size();
	lmax_pts.assign(num_pts, vmfloat3(0, 0, 0));
	//std::vector<__float3> tmp_lmax_pts(num_pts, __float3(FLT_MAX, FLT_MAX, FLT_MAX));
	vmfloat3 vol_size_f(_f3_(vol_info.vol_size.x, vol_info.vol_size.y, vol_info.vol_size.z));
	const float safe_bnd = 2.f;

#pragma omp parallel for num_threads( omp_get_num_procs() )
	for (int i = 0; i < num_pts; i++)
	{
		const vmfloat3& dir = dir_pts[i];
		if (fLengthVectorSq(&dir) < FLT_EPSILON)
		{
			lmax_pts[i] = pos_pts[i];
			continue;
		}

		float max_gmsq = -1;
		const vmfloat3& pos_start = pos_pts[i];
		vmfloat3 pos_lmax = pos_start;

		for (int j = 0; j < num_maxrc_steps; j++)
		{
			vmfloat3 pos_cur = pos_start + dir * dir_scale * (float)j, pos_cur_vs;
			fTransformPoint(&pos_cur_vs, &pos_cur, &vol_info.mat_ws2vs);

			if (pos_cur_vs.x < safe_bnd || pos_cur_vs.y < safe_bnd || pos_cur_vs.z < safe_bnd
				|| pos_cur_vs.x >= vol_size_f.x - safe_bnd
				|| pos_cur_vs.y >= vol_size_f.y - safe_bnd
				|| pos_cur_vs.z >= vol_size_f.z - safe_bnd)
				continue;

			if (j > 1 && band_info.vol_slices) // j > 1 is for safe zone
			{
				float mask_value = __Safe_TrilinearSample<char>(pos_cur_vs, band_info.vol_size, band_info.width_slice, EXB, band_info.vol_slices, -1);
				if (mask_value == 1.f || mask_value == -1.f)
					break;
			}

			vmfloat3 grad = SGRAD(pos_cur_vs, vol_info);
			float sample_v = __Safe_TrilinearSample<ushort>(pos_cur_vs, vol_info.vol_size, vol_info.width_slice, EXB, vol_info.vol_slices, -1);
			if (fDotVector(&dir, &grad) >= 0 || sample_v  < min_v || sample_v > max_v)
				continue;

			float gmsq = fLengthVectorSq(&grad);
			if (max_gmsq < gmsq)
			{
				max_gmsq = gmsq;
				pos_lmax = pos_cur;
			}
		}

//#define QUADRATIC
#ifdef QUADRATIC
		auto quadratic_maxsampler = [](const vmfloat3& pos_a, const vmfloat3& pos_c, const float gm_a, const float gm_b, const float gm_c)
		{
			vmfloat3 pos_max;
			float _div_ = (gm_a + gm_c) - 2.f * gm_b;
			if (_div_ == 0)
			{
				pos_max = pos_a;
			}
			else
			{
				float t = (3.f * gm_a + gm_c - 4.f * gm_b) / (4.f * _div_);
				pos_max = pos_a + min(max(t, 0.f), 1.f) * (pos_c - pos_a);
			}
			return pos_max;
		};
		if (max_gmsq > 0)
		{
			vmfloat3 pos_a = pos_lmax - dir * dir_scale, pos_a_vs;
			vmfloat3 pos_c = pos_lmax + dir * dir_scale, pos_c_vs;
			fTransformPoint(&pos_a_vs, &pos_a, &vol_info.mat_ws2vs);
			fTransformPoint(&pos_c_vs, &pos_c, &vol_info.mat_ws2vs);
			float edge_a = fLengthVector(&SGRAD(pos_a_vs, vol_info));
			float edge_c = fLengthVector(&SGRAD(pos_c_vs, vol_info));
			lmax_pts[i] = quadratic_maxsampler(pos_a, pos_c, edge_a, sqrt(max_gmsq), edge_c);
		}
		else
			lmax_pts[i] = pos_lmax;
#else
		lmax_pts[i] = pos_lmax;
#endif
	}
}

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <Eigen/Householder>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#define _DIM_ 3
auto __construct_covariance_matrix = [](const Eigen::VectorXd& cov) -> Eigen::MatrixXd
{
	Eigen::MatrixXd m(_DIM_, _DIM_);

	for (std::size_t i = 0; i < _DIM_; ++i)
	{
		for (std::size_t j = i; j < _DIM_; ++j)
		{
			m(i, j) = static_cast<float>(cov[(_DIM_ * i) + j - ((i * (i + 1)) / 2)]);

			if (i != j)
				m(j, i) = m(i, j);
		}
	}

	return m;
};
auto __diagonalize_selfadjoint_matrix = [](Eigen::MatrixXd& m, Eigen::MatrixXd& eigenvectors, Eigen::VectorXd& eigenvalues) -> bool
{
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;

	//eigensolver.computeDirect(m);
	eigensolver.compute(m);

	if (eigensolver.info() != Eigen::Success)
		return false;

	eigenvalues = eigensolver.eigenvalues();
	eigenvectors = eigensolver.eigenvectors();

	return true;
};
auto __diagonalize_selfadjoint_covariance_matrix = [&]
(const Eigen::VectorXd& cov, float* eigenvalues, vmfloat3* eigenvectors)
{
	Eigen::MatrixXd m = __construct_covariance_matrix(cov);

	// Diagonalizing the matrix
	Eigen::VectorXd eigenvalues_;
	Eigen::MatrixXd eigenvectors_;
	bool res = __diagonalize_selfadjoint_matrix(m, eigenvectors_, eigenvalues_);

	if (res)
	{
		for (std::size_t i = 0; i < _DIM_; ++i)
		{
			eigenvalues[i] = static_cast<float>(eigenvalues_[i]);

			for (std::size_t j = 0; j < _DIM_; ++j)
				((float*)eigenvectors)[_DIM_*i + j] = static_cast<float>(eigenvectors_(j, i));
		}
	}

	return res;
};

float estimate_shape_variation(const vmfloat3& pos_src, const vmfloat3& nrl_src, const vmfloat3* pos_pts, const vmfloat3* nrl_pts, const kd_tree_t& kdt_index, const float kernel_radius)
{
	float sv = 1.f;
	const float r_sq = kernel_radius * kernel_radius;

	static nanoflann::SearchParams params;
	params.sorted = false;

	std::vector<std::pair<size_t, float>> ret_matches;
	const int nMatches = (int)kdt_index.radiusSearch((float*)&pos_src, r_sq, ret_matches, params);

	std::vector<int> valid_idx_pts;
	for (int k = 0; k < nMatches; k++)
	{
		int idx_neighbor = (int)ret_matches[k].first;
		vmfloat3 nrl_nb = nrl_pts[idx_neighbor];
		if (nrl_src.x * nrl_nb.x + nrl_src.y * nrl_nb.y + nrl_src.z * nrl_nb.z >= 0)
			valid_idx_pts.push_back(idx_neighbor);
	}

	const int num_valid_nb_pts = (int)valid_idx_pts.size();
	if (num_valid_nb_pts < 3)
		return 1; // treated as an outlier

	vmdouble3 pos_centroid = vmdouble3();
	double ft_num_ptns = (double)nMatches;
	for (int k = 0; k < nMatches; k++)
	{
		int idx_nb = (int)ret_matches[k].first;
		vmdouble3 pos_nb = pos_pts[idx_nb];
		pos_centroid.x += ((double)pos_nb.x / ft_num_ptns);
		pos_centroid.y += ((double)pos_nb.y / ft_num_ptns);
		pos_centroid.z += ((double)pos_nb.z / ft_num_ptns);
	}

	vmfloat3 pos_centroid_f = vmfloat3((float)pos_centroid.x, (float)pos_centroid.y, (float)pos_centroid.z);
	//Eigen::VectorXd evecs(6); evecs << 0, 0, 0, 0, 0, 0;
	Eigen::VectorXd evecs = Eigen::VectorXd::Zero(6);
	for (int k = 0; k < nMatches; k++)
	{
		int idx_nb = (int)ret_matches[k].first;
		vmfloat3 pos_nb = pos_pts[idx_nb];
		vmfloat3 diff = vmfloat3(pos_centroid_f.x - pos_nb.x, pos_centroid_f.y - pos_nb.y, pos_centroid_f.z - pos_nb.z);
		evecs(0) += diff.x * diff.x;
		evecs(1) += diff.x * diff.y;
		evecs(2) += diff.x * diff.z;
		evecs(3) += diff.y * diff.y;
		evecs(4) += diff.y * diff.z;
		evecs(5) += diff.z * diff.z;
	}

	float eigenvalues[3];
	vmfloat3 eigenvectors[3];
	if (__diagonalize_selfadjoint_covariance_matrix(evecs, eigenvalues, eigenvectors))
	{
		// eigenvalues[0] is smallest
		float sum_egv = eigenvalues[0] + eigenvalues[1] + eigenvalues[2];
		sv = eigenvalues[0] / sum_egv * 3.f;
	}

	return sv;
}

float estimate_orient_variation(const vmfloat3& pos_src, const __VolSampleInfo<ushort>& vol_info, const float kernel_radius)
{
	auto __sample_v = [&](const vmfloat3& pos_vs)
	{
		return __Safe_TrilinearSample(pos_vs, vol_info.vol_size, vol_info.width_slice, EXB, vol_info.vol_slices);
	};

	float ov = 1.f;
	float voxel_kernel = kernel_radius / vol_info.min_sample_dist;
	vmfloat3 dirs[3] = { vol_info.vec_grad_dirs[0] * 0.5f * voxel_kernel, vol_info.vec_grad_dirs[1] * 0.5f * voxel_kernel, vol_info.vec_grad_dirs[2] * 0.5f * voxel_kernel };

	vmfloat3 pos_sample_vs;
	fTransformPoint(&pos_sample_vs, &pos_src, &vol_info.mat_ws2vs);

	float v = __sample_v(pos_sample_vs);
	float v_XXR = __sample_v(pos_sample_vs + 2.f * dirs[0]);
	float v_XXL = __sample_v(pos_sample_vs - 2.f * dirs[0]);
	float v_YYR = __sample_v(pos_sample_vs + 2.f * dirs[1]);
	float v_YYL = __sample_v(pos_sample_vs - 2.f * dirs[1]);
	float v_ZZR = __sample_v(pos_sample_vs + 2.f * dirs[2]);
	float v_ZZL = __sample_v(pos_sample_vs - 2.f * dirs[2]);
	float v_XR = __sample_v(pos_sample_vs + dirs[0]);
	float v_XL = __sample_v(pos_sample_vs - dirs[0]);
	float v_YR = __sample_v(pos_sample_vs + dirs[1]);
	float v_YL = __sample_v(pos_sample_vs - dirs[1]);
	float v_ZR = __sample_v(pos_sample_vs + dirs[2]);
	float v_ZL = __sample_v(pos_sample_vs - dirs[2]);
	float v_XRYR = __sample_v(pos_sample_vs + dirs[0] + dirs[1]);
	float v_XRYL = __sample_v(pos_sample_vs + dirs[0] - dirs[1]);
	float v_XLYR = __sample_v(pos_sample_vs - dirs[0] + dirs[1]);
	float v_XLYL = __sample_v(pos_sample_vs - dirs[0] - dirs[1]);
	float v_YRZR = __sample_v(pos_sample_vs + dirs[1] + dirs[2]);
	float v_YRZL = __sample_v(pos_sample_vs + dirs[1] - dirs[2]);
	float v_YLZR = __sample_v(pos_sample_vs - dirs[1] + dirs[2]);
	float v_YLZL = __sample_v(pos_sample_vs - dirs[1] - dirs[2]);
	float v_XRZR = __sample_v(pos_sample_vs + dirs[0] + dirs[2]);
	float v_XRZL = __sample_v(pos_sample_vs + dirs[0] - dirs[2]);
	float v_XLZR = __sample_v(pos_sample_vs - dirs[0] + dirs[2]);
	float v_XLZL = __sample_v(pos_sample_vs - dirs[0] - dirs[2]);

	vmmat44f H;

	vmfloat3 g = vmfloat3(v_XR - v_XL, v_YR - v_YL, v_ZR - v_ZL);
	H[0][0] = v_XXR - 2.f * v + v_XXL; // f_xx
	H[0][1] = (v_XRYR - v_XLYR - v_XRYL + v_XLYL);	//f_xy
	H[0][2] = (v_XRZR - v_XLZR - v_XRZL + v_XLZL);	//f_xz
	H[1][0] = H[0][1];
	H[1][1] = v_YYR - 2.f * v + v_YYL;								//f_yy
	H[1][2] = (v_YRZR - v_YLZR - v_YRZL + v_YLZL);	//f_yz
	H[2][0] = H[0][2];
	H[2][1] = H[1][2];
	H[2][2] = v_ZZR - 2.f * v + v_ZZL;								//f_zz
	H[3][3] = 1;
	//
	float gm = fLengthVector(&g);
	if (gm > FLT_EPSILON)
	{
		vmfloat3 n = -g / gm;
		vmmat44f nnT, P;
		nnT[0][0] = n.x*n.x;
		nnT[0][1] = n.x*n.y;
		nnT[0][2] = n.x*n.z;
		nnT[1][0] = n.y*n.x;
		nnT[1][1] = n.y*n.y;
		nnT[1][2] = n.y*n.z;
		nnT[2][0] = n.z*n.x;
		nnT[2][1] = n.z*n.y;
		nnT[2][2] = n.z*n.z;
		nnT[3][3] = 0;
		P[0][0] = 1.0f;
		P[1][1] = 1.0f;
		P[2][2] = 1.0f;
		P[3][3] = 1.0f;
		P = P - nnT;
		/////////////////////////
		vmmat44f FF = (-P * H) / (float)gm * nnT;
		float fFlowCurv = sqrt(FF[0][0] * FF[0][0] + FF[0][1] * FF[0][1] + FF[0][2] * FF[0][2]
			+ FF[1][0] * FF[1][0] + FF[1][1] * FF[1][1] + FF[1][2] * FF[1][2]
			+ FF[2][0] * FF[2][0] + FF[2][1] * FF[2][1] + FF[2][2] * FF[2][2]);
		ov = min(fFlowCurv, 1.0f);
	}
	return ov;
}

void compute_geometry_info(vector<vmfloat3>& nrl_pts, vector<vmfloat2>& gm_gc_pts, vector<vmfloat3>& pos_pts, const float kernel_radius, const __VolSampleInfo<ushort>& vol_info)
{
	int num_pts = (int)pos_pts.size();
	nrl_pts.assign(num_pts, vmfloat3());
	gm_gc_pts.assign(num_pts, vmfloat2());

#pragma omp parallel for num_threads( omp_get_num_procs() )
	for (int i = 0; i < num_pts; i++)
	{
		vmfloat3 pos_src = pos_pts[i], pos_vs;
		fTransformPoint(&pos_vs, &pos_src, &vol_info.mat_ws2vs);
		nrl_pts[i] = -SGRAD(pos_vs, vol_info);
		float leng = fLengthVector(&nrl_pts[i]);
		nrl_pts[i] = leng > FLT_EPSILON ? nrl_pts[i] / leng : vmfloat3();
		gm_gc_pts[i].x = leng;
	}

	PointCloud<float, vmfloat3> pc_kdt(pos_pts);
	kd_tree_t kdt(3, pc_kdt, nanoflann::KDTreeSingleIndexAdaptorParams(10));
	kdt.buildIndex();

#pragma omp parallel for num_threads( omp_get_num_procs() )
	for (int i = 0; i < num_pts; i++)
	{
		vmfloat3 pos_src = pos_pts[i];
		float sv = min(estimate_shape_variation(pos_src, nrl_pts[i], &pos_pts[0], &nrl_pts[0], kdt, kernel_radius), 1.f);
		float ov = min(estimate_orient_variation(pos_src, vol_info, kernel_radius), 1.f);
		gm_gc_pts[i].y = max(sv, ov);
	}
}

void clustering(std::map<int, std::vector<int>>& map_clusters, const std::vector<vmfloat3>& pos_pts, const kd_tree_t& kdt, const std::vector<vmfloat3>& nrl_pts, const float e_c)
{
	int num_pts = (int)pos_pts.size();
	if (num_pts == 0) return;

	nanoflann::SearchParams params;
	params.sorted = false;
	const float r_sq = e_c * e_c;

	auto is_side_angle = [](const vmfloat3& p, const vmfloat3& q, const float eps)
	{
		float angle = std::acos(max(min(fDotVector(&p, &q), 1.f), -1.f)); // 0 to PI
		return (angle > VM_PI / 2.f - eps) && (angle < VM_PI / 2.f + eps);
	};
	auto subs_norm = [](const vmfloat3& p, const vmfloat3& q)
	{
		vmfloat3 v = vmfloat3(p.x - q.x, p.y - q.y, p.z - q.z);
		fNormalizeVector(&v, &v);
		return v;
	};

	std::vector<int> cluster_map_pts;
	cluster_map_pts.assign(num_pts, 0);

	int count_id_cluster = 1;
	for (int i = 0; i < num_pts; i++)
	{
		if (cluster_map_pts[i] != 0) continue;

		std::vector<int> index_cluster_pts;
		std::queue<int> que_for_process;

		que_for_process.push(i);
		index_cluster_pts.push_back(i);
		cluster_map_pts[i] = count_id_cluster;

		int error_prop_count = 0;
		while (!que_for_process.empty())
		{
			if (error_prop_count++ > 50000000)
			{
				cout << "ERROR, num conf. clusters : " << map_clusters.size() << endl;
				return;
			}

			int source_idx = que_for_process.front();
			que_for_process.pop();

			const vmfloat3& pos_src = pos_pts[source_idx];
			const vmfloat3& nrl_src = nrl_pts[source_idx];

			std::vector<std::pair<size_t, float> >   ret_matches;
			const int nMatches = (int)kdt.radiusSearch((float*)&pos_src, r_sq, ret_matches, params);

			for (int j = 0; j < nMatches; j++)
			{
				int idx_nb = (int)ret_matches[j].first;
				vmfloat3 pos_nb = pos_pts[idx_nb];
				if (cluster_map_pts[idx_nb] == 0
					&& is_side_angle(nrl_src, subs_norm(pos_nb, pos_src), VM_fPI / 6.f)
					)
				{
					que_for_process.push(idx_nb);
					index_cluster_pts.push_back(idx_nb);
					cluster_map_pts[idx_nb] = count_id_cluster;
				}
			}
		}

		if (index_cluster_pts.size() > 0)
			map_clusters[count_id_cluster++] = index_cluster_pts;
	}
}

void estimating_hole(bool& is_boundary, bool& is_src_vector_set, vmfloat3& vector_o,
	double(&btw_angles_max)[4], double(&btw_angles_min)[4],
	const std::vector<std::pair<size_t, float> > ret_matches,
	const std::vector<vmfloat3>& pos_pts, const vmfloat3& pos_src, const vmfloat3& nrl_src,
	const float thres_susceptibility_sq)
{
	const int nMatches = (int)ret_matches.size();
	for (int j = 0; j < nMatches; j++)
	{
		int idx_nb = (int)ret_matches[j].first;

		const vmfloat3& pos_nb = pos_pts[idx_nb];

		// plane is defined by pos_src and nrl_src
		vmfloat3 _v = pos_nb - pos_src;
		vmfloat3 _tv = nrl_src * fDotVector(&_v, &nrl_src);
		vmfloat3 pos_nb_on_plane = pos_nb - _tv;
		vmfloat3 vec_on_plane = pos_nb_on_plane - pos_src;
		double length_sq = fLengthVectorSq(&vec_on_plane);

		if (length_sq > thres_susceptibility_sq)
		{
			vec_on_plane /= (float)std::sqrt(length_sq);

			if (!is_src_vector_set)
			{
				is_src_vector_set = true;
				vector_o = vec_on_plane;
			}
			else
			{
				double btw_angle = std::acos(max(min(fDotVector(&vector_o, &vec_on_plane), 1.f), -1.f)); // 0 to PI
				vmfloat3 cross_v;
				fCrossDotVector(&cross_v, &vec_on_plane, &vector_o);
				if (fDotVector(&nrl_src, &cross_v) < 0)
					btw_angle = 2 * VM_fPI - btw_angle;
				//ordered_angles.insert((float)btw_angle);

				if (btw_angle < VM_fPI * 0.5)
				{
					btw_angles_max[0] = max(btw_angle, btw_angles_max[0]);
					//btw_angles_min[0] = min(btw_angle, btw_angles_min[0]);//
				}
				else if (btw_angle < VM_fPI)
				{
					btw_angles_max[1] = max(btw_angle, btw_angles_max[1]);
					btw_angles_min[1] = min(btw_angle, btw_angles_min[1]);
				}
				else if (btw_angle < VM_fPI * 1.5)
				{
					btw_angles_max[2] = max(btw_angle, btw_angles_max[2]);
					btw_angles_min[2] = min(btw_angle, btw_angles_min[2]);
				}
				else
				{
					//btw_angles_max[3] = max(btw_angle, btw_angles_max[3]); //
					btw_angles_min[3] = min(btw_angle, btw_angles_min[3]);
				}

				if (btw_angles_min[1] - btw_angles_max[0] < VM_fPI / 2.
					&& btw_angles_min[2] - btw_angles_max[1] < VM_fPI / 2.
					&& btw_angles_min[3] - btw_angles_max[2] < VM_fPI / 2.)
				{
					is_boundary = false;
					break; // for (int j = 0; j < nMatches; j++)
				}
			}
		}
	}
}

void detecting_boundary(std::vector<int>& bnd_id_pts, const std::vector<vmfloat3>& pos_pts, const std::vector<vmfloat3>& nrl_pts,
	const kd_tree_t& kdt_index, const float e_b, const bool consistent_dir_neighbors/*projecting except opposite directional points of the target poin*/)
{
	float thres_susceptibility_sq = e_b * 0.01f;
	thres_susceptibility_sq *= thres_susceptibility_sq;
	int num_pts = (int)pos_pts.size();

	nanoflann::SearchParams params;
	params.sorted = false;
	const float r_sq = e_b * e_b;
	// params.eps

	bool* hole_flags = new bool[num_pts];
	ZeroMemory(hole_flags, sizeof(bool) * num_pts);

#pragma omp parallel for num_threads( omp_get_num_procs() )
	for (int i = 0; i < num_pts; i++)
	{
		const vmfloat3& pos_src = pos_pts[i];
		const vmfloat3& nrl_src = nrl_pts[i];

		std::vector<std::pair<size_t, float> >   ret_matches;
		const int nMatches = (int)kdt_index.radiusSearch((float*)&pos_src, r_sq, ret_matches, params);

		//std::set<float> ordered_angles;
		bool is_src_vector_set = false;
		vmfloat3 vector_o;

		double btw_angles_max[4] = { -100000., -100000., -100000., -100000. };
		double btw_angles_min[4] = { 100000., 100000., 100000., 100000. };
		bool is_boundary = true;
		// criterion 2 angle criterion
		std::vector<std::pair<size_t, float> >   dir_matches;
		{
			for (int j = 0; j < nMatches; j++)
			{
				vmfloat3 nrl_nb = nrl_pts[ret_matches[j].first];
				if (fDotVector((vmfloat3*)&nrl_src, (vmfloat3*)&nrl_nb) >= 0 || !consistent_dir_neighbors)
					dir_matches.push_back(ret_matches[j]);
			}
		}

		estimating_hole(is_boundary, is_src_vector_set, vector_o,
			btw_angles_max, btw_angles_min, dir_matches, pos_pts, pos_src, nrl_src, thres_susceptibility_sq);

		// store boundary result
		if (is_boundary) hole_flags[i] = true;
	}

	for (int i = 0; i < num_pts; i++)
		if (hole_flags[i]) bnd_id_pts.push_back(i);

	delete[] hole_flags;
}


#define RANGE_SURF_STEP_MAX_NUM 10 
#define SURFACE_REFINEMENT_NUM 5 
void relocate_to_target_densities(vmfloat3* pos_relocated_pts,
	const vmfloat3* pos_pts, const vmfloat3* nrl_pts, const float* densities_pts, const int num_pts,
	const float sample_dist, const __VolSampleInfo<ushort>& vol_info)
{
	auto __sample_v = [&](const vmfloat3& pos_vs)
	{
		return __Safe_TrilinearSample(pos_vs, vol_info.vol_size, vol_info.width_slice, EXB, vol_info.vol_slices);
	};

	auto __traverse_serach_dstv = [&__sample_v, &vol_info](vmfloat3& pos_dst, const vmfloat3& pos_start, const vmfloat3& vec_sample, const float dst_v, const bool from_lower)
	{
		for (int j = 1; j < RANGE_SURF_STEP_MAX_NUM; j++)
		{
			vmfloat3 pos_sample = pos_start + vec_sample * (float)j, pos_sample_vs;
			fTransformPoint(&pos_sample_vs, &pos_sample, &vol_info.mat_ws2vs);
			float sample_v = __sample_v(pos_sample_vs);

			if (
				(sample_v > dst_v && from_lower)
				|| (sample_v < dst_v && !from_lower)
				)
			{
				// requires boundary-fitting
				vmfloat3 pos_S = pos_sample - vec_sample;
				vmfloat3 pos_E = pos_sample;
				vmfloat3 pos_O = pos_S;
				for (uint k = 0; k < SURFACE_REFINEMENT_NUM; k++)
				{
					vmfloat3 pos_bis = (pos_S + pos_E) * 0.5f, pos_bis_vs;
					fTransformPoint(&pos_bis_vs, &pos_bis, &vol_info.mat_ws2vs);
					sample_v = __sample_v(pos_bis_vs);

					if (
						(sample_v > dst_v && from_lower)
						|| (sample_v < dst_v && !from_lower))
						pos_E = pos_bis;
					else
						pos_S = pos_bis;
				}

				pos_dst = pos_S;
				return true;
			}
		} // for (uint j = 1; j < RANGE_SURF_STEP_MAX_NUM; j++)
		return false;
	};

	//pos_relocated_pts.clear();
	//pos_relocated_pts.assign(num_pts, vmfloat3());

#pragma omp parallel for num_threads( omp_get_num_procs() )
	for (int i = 0; i < num_pts; i++)
	{
		const vmfloat3& pos_sample = pos_pts[i];
		const vmfloat3& nrl_sample = nrl_pts[i];
		vmfloat3 pos_sample_vs;
		float dst_v = densities_pts[i];
		if (dst_v == 7777777)
		{
			pos_relocated_pts[i] = pos_sample;
			continue; // magic number exception
		}
		fTransformPoint(&pos_sample_vs, &pos_sample, &vol_info.mat_ws2vs);
		float sample_v = __sample_v(pos_sample_vs);

		vmfloat3 pos_dst_0, pos_dst_1;
		bool ret_0 = __traverse_serach_dstv(pos_dst_0, pos_sample, nrl_sample * sample_dist, dst_v, sample_v < dst_v);
		bool ret_1 = __traverse_serach_dstv(pos_dst_1, pos_sample, -nrl_sample * sample_dist, dst_v, sample_v < dst_v);
		if (ret_0 && ret_1)
		{
			pos_relocated_pts[i] = fLengthVectorSq(&(pos_sample - pos_dst_0)) < fLengthVectorSq(&(pos_sample - pos_dst_1)) ? pos_dst_0 : pos_dst_1;
		}
		else if (ret_0)
		{
			pos_relocated_pts[i] = pos_dst_0;
		}
		else if (ret_1)
		{
			pos_relocated_pts[i] = pos_dst_1;
		}
		else
		{
			pos_relocated_pts[i] = pos_sample;
		}
	} // for (uint i = 0; i < numVoxels; i++)
}

//#include <opencv2/highgui.hpp>
// ui work
void register_pobj(VmVObjectPrimitive& pobj, const vector<vmfloat3>& pos_pts, const vector<vmfloat3>& nrl_pts, const vector<vmfloat3>* ptr_clr_pts)
{
	PrimitiveData vtx_data;
	vtx_data.is_ccw = true; // NA //
	vtx_data.is_stripe = false; // NA//
	vtx_data.check_redundancy = true; // NA //
	vtx_data.ptype = PrimitiveTypePOINT;
	vtx_data.idx_stride = 1;
	vtx_data.num_vtx = (uint)pos_pts.size();
	vtx_data.num_prims = vtx_data.num_vtx;

	vmfloat3* buf_vtx_pos = new vmfloat3[vtx_data.num_vtx];
	vmfloat3* buf_nrl_pos = new vmfloat3[vtx_data.num_vtx];
	vmfloat3* buf_clr_pos = new vmfloat3[vtx_data.num_vtx];
	memcpy(buf_vtx_pos, &pos_pts[0], sizeof(vmfloat3) * vtx_data.num_vtx);
	memcpy(buf_nrl_pos, &nrl_pts[0], sizeof(vmfloat3) * vtx_data.num_vtx);
	if (ptr_clr_pts)
	{
		memcpy(buf_clr_pos, &ptr_clr_pts->at(0), sizeof(vmfloat3) * vtx_data.num_vtx);
	}
	else
	{
		for (uint i = 0; i < vtx_data.num_vtx; i++)
			buf_clr_pos[i] = vmfloat3(1);
	}

	vtx_data.ReplaceOrAddVerticeDefinition("POSITION", buf_vtx_pos);
	vtx_data.ReplaceOrAddVerticeDefinition("NORMAL", buf_nrl_pos);
	vtx_data.ReplaceOrAddVerticeDefinition("TEXCOORD0", buf_clr_pos); // special case //
	vtx_data.ComputeOrthoBoundingBoxWithCurrentValues();

	pobj.RegisterPrimitiveData(vtx_data);
	pobj.UpdateKDTree();
	pobj.RegisterCustomParameter("_bool_ApplyShadingFactors", true);
	pobj.RegisterCustomParameter("_bool_CtModelerStep", true);
}

#define COLORMAP_SIZE 1024
void fill_jet_colormap(int* colorarray, int ary_size)
{
	int prev_idx = 0;

	float gap = (float)ary_size / 4.f;
	vmfloat3 redf = vmfloat3(1.f, 0, 0);
	vmfloat3 greenf = vmfloat3(0, 1.f, 0);
	vmfloat3 bluef = vmfloat3(0, 0, 1.f);

	int index_1 = (int)(gap / 2.f + 0.5f);
	int index_2 = (int)(gap / 2.f + gap + 0.5f);
	int index_3 = (int)(gap / 2.f + 2 * gap + 0.5f);
	int index_4 = (int)(gap / 2.f + 3 * gap + 0.5f);
	int index_5 = ary_size;

	auto float3_2_int = [&](vmfloat3& c)
	{
		byte r = (byte)__min(c.x * 255.f, 255.f);
		byte g = (byte)__min(c.y * 255.f, 255.f);
		byte b = (byte)__min(c.z * 255.f, 255.f);
		return (r << 16) | (g << 8) | (b);
	};

	for (int i = 0; i <= index_1; i++)
	{
		vmfloat3 clr = bluef * (0.5f + 0.5f * (float)i / (float)index_1);
		colorarray[i] = float3_2_int(clr);
	}
	for (int i = index_1; i <= index_2; i++)
	{
		float ratio = (float)(i - index_1) / (float)(index_2 - index_1);
		vmfloat3 clr = bluef + greenf * ratio;
		colorarray[i] = float3_2_int(clr);
	}
	for (int i = index_2; i <= index_3; i++)
	{
		float ratio = (float)(i - index_2) / (float)(index_3 - index_2);
		vmfloat3 clr = greenf + bluef * (1.f - ratio) + redf * ratio;
		colorarray[i] = float3_2_int(clr);
	}
	for (int i = index_3; i <= index_4; i++)
	{
		float ratio = (float)(i - index_3) / (float)(index_4 - index_3);
		vmfloat3 clr = redf + greenf * (1.f - ratio);
		colorarray[i] = float3_2_int(clr);
	}
	for (int i = index_4; i < index_5; i++)
	{
		vmfloat3 clr = redf * (1.f - 0.5f * (float)(i - index_4) / (float)(index_5 - index_4));
		colorarray[i] = float3_2_int(clr);
	}
}

void fill_cool_colormap(int* colorarray, int ary_size)
{
	int prev_idx = 0;

	float gap = (float)ary_size / 4.f;
	vmfloat3 pinkf = vmfloat3(1.f, 0, 1.f);
	vmfloat3 skyf = vmfloat3(0, 1.f, 1.f);

	auto float3_2_int = [&](vmfloat3& c)
	{
		byte r = (byte)__min(c.x * 255.f, 255.f);
		byte g = (byte)__min(c.y * 255.f, 255.f);
		byte b = (byte)__min(c.z * 255.f, 255.f);
		return (r << 16) | (g << 8) | (b);
	};

	for (int i = 0; i < ary_size; i++)
	{
		float ratio = (float)(i) / (float)(ary_size);
		vmfloat3 clr = skyf * (1 - ratio) + pinkf * ratio;
		colorarray[i] = float3_2_int(clr);
	}
}

inline void convert_int_to_float3(const int _color, vmfloat3& rgb)
{
	rgb.r = ((_color >> 16) & 0xFF) / 255.f;
	rgb.g = ((_color >> 8) & 0xFF) / 255.f;
	rgb.b = ((_color >> 0) & 0xFF) / 255.f;
}


//auto is_inside_box = [](vmfloat3& p, vmfloat3& bp, float bs)
//{
//	vmfloat3 pos_min = bp - vmfloat3(bs);
//	vmfloat3 pos_max = bp + vmfloat3(bs);
//	return p.x >= pos_min.x && p.x <= pos_max.x
//		&& p.y >= pos_min.y && p.y <= pos_max.y
//		&& p.z >= pos_min.z && p.z <= pos_max.z;
//};
//auto cout_f3 = [](string prefix, vmfloat3& p)
//{
//	cout << prefix << p.x << ", " << p.y << ", " << p.z << endl;
//};