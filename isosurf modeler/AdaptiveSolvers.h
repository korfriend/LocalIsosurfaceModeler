#pragma once

#include <vector>
#include <string>
#include <iostream>

typedef struct __float3__
{
	float x, y, z;
	__float3__() { x = y = z = 0; };
	__float3__(float _x, float _y, float _z) { x = _x; y = _y; z = _z; };
	__float3__(int _color) {
		x = (float)((_color >> 16) & 0xFF) / 255.f;
		y = (float)((_color >> 8) & 0xFF) / 255.f;
		z = (float)((_color >> 0) & 0xFF) / 255.f;
	};
}  __float3;

typedef struct __float2__
{
	float x, y;
	__float2__() { x = y = 0; };
	__float2__(float _x, float _y) { x = _x; y = _y; };
}  __float2;

template <typename T>
struct __ProcBuffers
{
	T* pos_pt;
	T* nrl_pt;
	unsigned int* index_buffer;
	unsigned int num_pts;
	unsigned int num_primitives;

	__ProcBuffers() {
		pos_pt = 0;
		nrl_pt = 0;
		index_buffer = 0;
		num_pts = num_primitives = 0;
		res = 0;
		level_set = NULL;
		iso_value = 0;
	};
	~__ProcBuffers() { Delete(); };

	float iso_value;
	float** level_set;
	int res;

	void Delete()
	{
		if (pos_pt != 0) delete[] pos_pt;
		if (nrl_pt != 0) delete[] nrl_pt;
		if (index_buffer != 0) delete[] index_buffer;
		if (level_set != 0) delete[] level_set;
		pos_pt = NULL;
		nrl_pt = NULL;
		index_buffer = NULL;
		level_set = NULL;
		iso_value = 0;
	}

	void GetPointsList(std::vector<T>& _pos_pts)
	{
		for (unsigned int i = 0; i < num_pts; i++)
			_pos_pts.push_back(pos_pt[i]);
	}
};

extern "C" __declspec(dllexport) bool ScreenedPoissonSurface2D(__ProcBuffers<__float2>* pOut,
	const std::vector<__float2>& pos_pts, const std::vector<__float2>& nrl_pts,
	const std::vector<__float2>& pos_aux_pts, const std::vector<__float2>& nrl_aux_pts,
	std::vector<char*>& cmd_sps_params, std::vector<void*>& additional_params);

extern "C" __declspec(dllexport) bool ScreenedPoissonSurface3D(__ProcBuffers< __float3>* pOut,
	const std::vector<__float3>& pos_pts, const std::vector<__float3>& nrl_pts,
	const std::vector<__float3>& pos_aux_pts, const std::vector<__float3>& nrl_aux_pts,
	std::vector<char*>& cmd_sps_params, std::vector<void*>& additional_params);
