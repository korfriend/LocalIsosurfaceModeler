#include <string>
#include <iostream>
#include <any>

#include "VisMtvApi.h"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

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
	float g_h = 5000.f; // control parameter (ct value)
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

	std::map<std::string, int> prim_ids = { {"candidate_pts_id", candidate_pts_id} , {"reliable_pts_id", reliable_pts_id} , {"holefilling_pts_id", holefilling_pts_id} , 
	{"ls_reliable_pts_id", ls_reliable_pts_id} , {"ls_holefilling_pts_id", ls_holefilling_pts_id} , {"final_mesh_id", final_mesh_id} };

	// test rendering example
#define __cv3__ *(glm::fvec3*)
#define __cv4__ *(glm::fvec4*)
#define __cm4__ *(glm::fmat4x4*)
	vzm::CameraParameters cam_params;
	__cv3__ cam_params.pos = glm::fvec3(0, 0, 150);
	__cv3__ cam_params.up = glm::fvec3(0, 1.f, 0);
	__cv3__ cam_params.view = glm::fvec3(0, 0, -1.f);
	cam_params.fov_y = 3.141592654f / 4.f;
	cam_params.aspect_ratio = 640.f / 480.f;
	cam_params.projection_mode = 2;
	cam_params.w = 640;
	cam_params.h = 480;

	vzm::SceneEnvParameters scn_env_params;
	scn_env_params.is_on_camera = true;
	scn_env_params.is_pointlight = true;
	__cv3__ scn_env_params.pos_light = __cv3__ cam_params.pos;
	__cv3__ scn_env_params.dir_light = __cv3__ cam_params.view;

	vzm::ObjStates obj_state;
	obj_state.emission = 0.4f;
	obj_state.diffusion = 0.6f;
	obj_state.specular = 0.2f;
	obj_state.sp_pow = 30.f;
	__cv4__ obj_state.color = glm::fvec4(0.8f, 1.f, 1.f, 0.7f);
	__cm4__ obj_state.os2ws = glm::fmat4x4(); // identity
	
	int scene_id = 0;
	auto show_window = [](const std::string& title, const int scene_id, const int cam_id)
	{
		static uint count = 0;
		static double times_sum = 0;
		vzm::RenderScene(scene_id, cam_id);
		std::cout << "rendering sec : " << times_sum / (double)(++count) << " s" << std::endl;
		unsigned char* ptr_rgba;
		float* ptr_zdepth;
		int w, h;
		if (vzm::GetRenderBufferPtrs(scene_id, &ptr_rgba, &ptr_zdepth, &w, &h, cam_id))
		{
			cv::Mat cvmat(h, w, CV_8UC4, ptr_rgba);
			//show the image
			cv::imshow(title, cvmat);
		}
	};
	for (auto& it : prim_ids)
	{
		vzm::SetSceneEnvParameters(scene_id, scn_env_params);
		vzm::SetCameraParameters(scene_id, cam_params, 0);
		vzm::ReplaceOrAddSceneObject(scene_id, it.second, obj_state);

		cv::namedWindow(it.first, 1);
		show_window(it.first, scene_id, 0);
		scene_id++;
	}

	cv::waitKey(0);

	vzm::DeinitEngineLib();
	return 0;
}