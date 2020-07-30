#pragma once

#ifndef __dojo_export
#define __dojo_export extern "C" __declspec(dllexport)
#endif

#include "VimCommon.h"

// OpenGL defaults to column major 'order' matrices based on col-major computation, which means access is of the format mat[col][row].
// https://en.wikipedia.org/wiki/Row-_and_column-major_order

// "return 0" means success

// input
__dojo_export int load_foreground(const std::string& file_name, const vmobjects::VmVObjectVolume& main_volume, vmobjects::VmVObjectVolume& fore_volume,
	const bool flip_x = false, const bool flip_y = false, const bool flip_z = false);

__dojo_export int processing_stage1(vmobjects::VmVObjectPrimitive& candidate_pts, vmobjects::VmVObjectVolume& filtered_volume,
	const vmobjects::VmVObjectVolume& main_volume, const vmobjects::VmVObjectVolume& fore_volume, const int t_in /*voxel unit*/, const int t_out /*voxel unit*/,
	const float min_sample_v = 0, const float max_sample_v = 65535.f/*additional feature since 20200215*/,
	const int gaussian_kernel_width = 3 /*must be odd, recommend 1 (no filter) or 3 or 5*/, 
	const float simplify_grid_length = 1.f /*voxel unit*/,
	const float geometric_complexity_kernel_ratio = 0.01 /*w.r.t. bounding box of target geometry*/);

__dojo_export int processing_stage2(vmobjects::VmVObjectPrimitive& reliable_pts, vmobjects::VmVObjectPrimitive& holefilling_pts,
	const vmobjects::VmVObjectVolume& sample_volume, const vmobjects::VmVObjectPrimitive& candidate_pts, const float g_h, const float mu_u,
	const float epsilon = 2.f /*voxel unit, connecting radius*/,
	const int connectivity_criterion = 5, /*number of connecting points for the connectivity test*/
	const float epsilon_b = 3.f /*voxel unit, for detecting the hole boundary*/,
	const float angle_criterion = 90.f /*degree, for detecting the hole boundary*/);// ,
//	std::vector<vmfloat3>& vtx_strong_bnd, std::vector<vmfloat3>& vtx_confidence_bnd);

// user interaction at stage 2
__dojo_export int update_strongweak_points(vmobjects::VmVObjectPrimitive& candidate_pts, const float g_h, const float mu_u);

// user interaction at stage 3
__dojo_export int update_holefill_points(vmobjects::VmVObjectPrimitive& holefill_pts, const int eta, const float m, const bool show_nc);

// user interaction at stage 4
__dojo_export int localisosurface_points(vmobjects::VmVObjectPrimitive& sm_reliable_pts, vmobjects::VmVObjectPrimitive& sm_holefill_pts,
	const vmobjects::VmVObjectVolume& original_volume, const vmobjects::VmVObjectVolume& filtered_volume,
	const vmobjects::VmVObjectPrimitive& reliable_pts, const vmobjects::VmVObjectPrimitive& holefill_pts, const int e_s, const int eta);

// output
__dojo_export int processing_final(vmobjects::VmVObjectPrimitive& surface_mesh,
	const vmobjects::VmVObjectPrimitive& sm_reliable_pts, const vmobjects::VmVObjectPrimitive& sm_holefill_pts, const vmobjects::VmVObjectPrimitive& holefill_pts, 
	const int eta, const float m, const int otlev);