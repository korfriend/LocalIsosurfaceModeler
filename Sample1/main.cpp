#include <string>
#include <iostream>
#include <any>

#include "VisMtvApi.h"

//#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtx/transform.hpp>
//#include <glm/gtc/constants.hpp>
//#include <glm/glm.hpp>
//#include <glm/gtc/type_ptr.hpp>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

int main()
{
#if (_MSVC_LANG < 201703L)
	# this program uses C++ standard ver 17
#endif

	vzm::InitEngineLib();

	int ct_vol_id = 0;
	vzm::LoadModelFile("D:\\Data\\Experiments 2019\\synthesis\\testcube_64x64x64_1.0x1.0x1.0_ushort.den", ct_vol_id);
	int foreground_vol_id = 0;
	vzm::GenerateEmptyVolume(foreground_vol_id, ct_vol_id, "BYTE", 0, 255, 0);

	vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { ct_vol_id , foreground_vol_id }, 
		{ {"function_name", std::string("load_foreground")}, {"file_name", std::string("")}, {"iso_value", 50000}, {"downsacled_factor", -1},
		{"flip_x", false}, {"flip_y", false}, {"flip_z", false} });

	int candidate_pts_id = 0, filtered_vol_id = 0;
	vzm::GenerateEmptyVolume(filtered_vol_id, ct_vol_id, "USHORT", 0, 65535, 0);
	vzm::GenerateEmptyPrimitive(candidate_pts_id);

	int t_in = 1; // control parameter (vox)
	int t_out = 7; // control parameter (vox)
	vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { candidate_pts_id , filtered_vol_id, ct_vol_id, foreground_vol_id },
		{ {"function_name", std::string("processing_stage1")}, {"t_in", t_in}, {"t_out", t_out}, {"min_sample_v", 0.f}, {"max_sample_v", 65535.f},
		{"gaussian_kernel_width", 3}, {"simplify_grid_length", 1.f}, {"geometric_complexity_kernel_ratio", 0.01f} });

	// w/ interactive point set rendering
	float g_h = 3000.f; // control parameter (ct value)
	float mu_u = 0.3f; // control parameter (0.0 ~ 1.0)
	{
		vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { candidate_pts_id },
			{ {"function_name", std::string("update_strongweak_points")}, {"g_h", g_h}, {"mu_u", mu_u} });
	}

	int reliable_pts_id = 0, holefilling_pts_id = 0;
	vzm::GenerateEmptyPrimitive(reliable_pts_id);
	vzm::GenerateEmptyPrimitive(holefilling_pts_id);
	vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { reliable_pts_id , holefilling_pts_id, filtered_vol_id, candidate_pts_id },
		{ {"function_name", std::string("processing_stage2")}, {"g_h", g_h}, {"mu_u", mu_u}, {"epsilon", 2.f}, {"connectivity_criterion", 5},
		{"epsilon_b", 3.f}, {"angle_criterion", 90.f /*fixed*/} });
	
	// w/ interactive point set rendering
	bool show_nc = true; // 
	int eta = 5; // control parameter (to avoid noise clusters)
	float m = 4.f; // control parameter (to suppress noise effects)
	{
		vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { holefilling_pts_id },
			{ {"function_name", std::string("update_holefill_points")}, {"eta", eta}, {"m", m}, {"show_nc", show_nc} });
	}

	int ls_reliable_pts_id = 0, ls_holefilling_pts_id = 0;
	vzm::GenerateEmptyPrimitive(ls_reliable_pts_id);
	vzm::GenerateEmptyPrimitive(ls_holefilling_pts_id);
	// w/ semi-interactive point set rendering
	int epsilon_s = 10; // control parameter (to get local isovalues)
	{
		vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { ls_reliable_pts_id, ls_holefilling_pts_id, ct_vol_id, filtered_vol_id, reliable_pts_id, holefilling_pts_id },
			{ {"function_name", std::string("localisosurface_points")}, {"epsilon_s", epsilon_s}, {"eta", eta} });
	}

	int final_mesh_id = 0;
	vzm::GenerateEmptyPrimitive(final_mesh_id);
	vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { final_mesh_id, ls_reliable_pts_id, ls_holefilling_pts_id, holefilling_pts_id },
		{ {"function_name", std::string("processing_final")}, {"eta", eta}, {"m", m}, {"otlev", -1 /*-1 means automatic decision for octree level, unless set 5 to 11 (int)*/} });

	vzm::DeinitEngineLib();
	return 0;
}