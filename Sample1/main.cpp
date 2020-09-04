#include <string>
#include <iostream>
#include <chrono>
#include <thread>

#include "VisMtvApi.h"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#define __cv3__ *(glm::fvec3*)
#define __cv4__ *(glm::fvec4*)
#define __cm4__ *(glm::fmat4x4*)

static int scene_id_cnt = 0;
static std::map<int, std::string> scene_name;
auto show_window = [](const std::string& title, const int scene_id, const int cam_id)
{
	vzm::RenderScene(scene_id, cam_id);
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
void CallBackFunc_Mouse(int event, int x, int y, int flags, void* userdata)
{
	using namespace cv;

	static int x_old = x;
	static int y_old = y;

	//if (x_old == x && y_old == y) return;
	//x_old = x;
	//y_old = y;

	vzm::CameraParameters cam_params;
	vzm::GetCameraParameters(0, cam_params, 0);

	static helpers::arcball aball_vr;
	if (event == EVENT_LBUTTONDOWN || event == EVENT_RBUTTONDOWN)
	{
		aball_vr.intializer((float*)&glm::fvec3(), 500);

		helpers::cam_pose arc_cam_pose;
		glm::fvec3 pos = __cv3__ arc_cam_pose.pos = __cv3__ cam_params.pos;
		__cv3__ arc_cam_pose.up = __cv3__ cam_params.up;
		__cv3__ arc_cam_pose.view = __cv3__ cam_params.view;
		aball_vr.start((int*)&glm::ivec2(x, y), (float*)&glm::fvec2(cam_params.w / 2, cam_params.h / 2), arc_cam_pose);
	}
	else if (event == EVENT_MBUTTONDOWN)
	{
	}
	else if (event == EVENT_MOUSEHWHEEL)
	{
	}
	else if (event == EVENT_MOUSEMOVE)
	{
		if ((flags & EVENT_FLAG_LBUTTON) || (flags & EVENT_FLAG_RBUTTON))
		{
			helpers::cam_pose arc_cam_pose;
			if (flags & EVENT_FLAG_LBUTTON)
				aball_vr.pan_move((int*)&glm::ivec2(x, y), arc_cam_pose);
			else if (flags & EVENT_FLAG_RBUTTON)
				aball_vr.move((int*)&glm::ivec2(x, y), arc_cam_pose);

			__cv3__ cam_params.pos = __cv3__ arc_cam_pose.pos;
			__cv3__ cam_params.up = __cv3__ arc_cam_pose.up;
			__cv3__ cam_params.view = __cv3__ arc_cam_pose.view;

			for (auto& it : scene_name)
			{
				vzm::SetCameraParameters(it.first, cam_params);
			}
		}
	}
}

int main()
{
#if (_MSVC_LANG < 201703L)
	# this program uses C++ standard ver 17
#endif

	vzm::InitEngineLib();

	int ct_vol_id = 0;

	// notice : all control parameters and initial parameters can be manipulated while exploring the dataset with interactive 3D rendering.
	// (refer to the video files)

	//vzm::LoadModelFile("..\\data\\data_1.x3d", ct_vol_id);
	//const int thres_value = 41500;
	//const int t_in_init = 6;
	//const int t_out_init = 2;
	//const int g_h_init = 4500;

	vzm::LoadModelFile("..\\data\\data_2.x3d", ct_vol_id);
	const int thres_value = 18800;
	const int t_in_init = 8;
	const int t_out_init = 0;
	const int g_h_init = 3600;
	
	int foreground_vol_id = 0;
	vzm::GenerateEmptyVolume(foreground_vol_id, ct_vol_id, "BYTE", 0, 255, 0);

	vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { ct_vol_id , foreground_vol_id }, 
		{ {"function_name", std::string("load_foreground")}, {"file_name", std::string("")}, {"iso_value", thres_value}, {"downsacled_factor", -1},
		{"flip_x", false}, {"flip_y", false}, {"flip_z", false} });

	int candidate_pts_id = 0, filtered_vol_id = 0;
	vzm::GenerateEmptyVolume(filtered_vol_id, ct_vol_id, "USHORT", 0, 65535, 0);
	vzm::GenerateEmptyPrimitive(candidate_pts_id);

	int t_in = t_in_init; // control parameter (vox)
	int t_out = t_out_init; // control parameter (vox)
	const float min_ct_v = 0.f, max_ct_v = 65535.f; // using full range of CT values... these can be used as control parameters as well.
	vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { candidate_pts_id , filtered_vol_id, ct_vol_id, foreground_vol_id },
		{ {"function_name", std::string("processing_stage1")}, {"t_in", t_in}, {"t_out", t_out}, {"min_sample_v", min_ct_v}, {"max_sample_v", max_ct_v},
		{"gaussian_kernel_width", 5}, {"simplify_grid_length", 1.f}, {"geometric_complexity_kernel_ratio", 0.01f} });

	// w/ interactive point set rendering
<<<<<<< HEAD
	float g_h = g_h_init; // control parameter (ct value)
=======
	float g_h = 5000.f; // control parameter (ct value)
>>>>>>> d5179c6e148847296050ad00246f81df2c9d5e49
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
	int eta = 30; // control parameter (to avoid noise clusters)
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
		{ {"function_name", std::string("processing_final")}, {"eta", eta}, {"m", m}, {"otlev", -1 /*-1 means automatic decision for octree level, unless, set 5 to 11 (int)*/} });

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
	obj_state.point_thickness = 0.f; // control for visible point size
	__cv4__ obj_state.color = glm::fvec4(0.8f, 1.f, 1.f, 1.f); 
	obj_state.color[3] = 0.9f; // control for transparency
	__cm4__ obj_state.os2ws = glm::fmat4x4(); // identity

	scene_name[0] = "candidate points";
	scene_name[1] = "confidence points";
	scene_name[2] = "final mesh";
	std::map<std::string, std::list<int>> prim_ids = { 
		{scene_name[0], {candidate_pts_id}} ,
		{scene_name[1], {reliable_pts_id, holefilling_pts_id}} ,
		{scene_name[2], {final_mesh_id}}
	};

	for (int i = 0; i < (int)scene_name.size(); i++)
	{
		vzm::SetSceneEnvParameters(i, scn_env_params);
		vzm::SetCameraParameters(i, cam_params, 0);

		auto it = prim_ids.find(scene_name[i]);
		for (int& obj_id : it->second)
			vzm::ReplaceOrAddSceneObject(i, obj_id, obj_state);

		cv::namedWindow(scene_name[i], cv::WINDOW_AUTOSIZE);
		cv::setMouseCallback(scene_name[i], CallBackFunc_Mouse, NULL);// &scenes[i]);
		show_window(scene_name[i], i, 0);
	}

	int key = -1;
	while (key != 'q')
	{
		for (auto& it : scene_name)
		{
			show_window(it.second, it.first, 0);
		}
		key = cv::waitKey(1);
	}

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