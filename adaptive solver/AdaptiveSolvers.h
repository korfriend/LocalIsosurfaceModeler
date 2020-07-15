#pragma once

#include "InteropHeader.h"

extern "C" __declspec(dllexport) bool ScreenedPoissonSurface2D(__ProcBuffers<__float2>* pOut,
	const std::vector<__float2>& pos_pts, const std::vector<__float2>& nrl_pts,
	const std::vector<__float2>& pos_aux_pts, const std::vector<__float2>& nrl_aux_pts,
	std::vector<char*>& cmd_sps_params, std::vector<void*>& additional_params);

extern "C" __declspec(dllexport) bool ScreenedPoissonSurface3D(__ProcBuffers< __float3>* pOut,
	const std::vector<__float3>& pos_pts, const std::vector<__float3>& nrl_pts,
	const std::vector<__float3>& pos_aux_pts, const std::vector<__float3>& nrl_aux_pts,
	std::vector<char*>& cmd_sps_params, std::vector<void*>& additional_params);
