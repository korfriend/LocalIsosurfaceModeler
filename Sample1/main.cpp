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

	const int base_scene_id = 10;
	vzm::CameraParameters cam_params;
	vzm::GetCameraParameters(base_scene_id, cam_params, 0);

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
	else if (event == EVENT_MOUSEWHEEL && (flags & EVENT_FLAG_CTRLKEY))
	{
		glm::fvec3 pos_min, pos_max;
		vzm::GetSceneBoundingBox({}, base_scene_id, (float*)&pos_min, (float*)&pos_max);
		float scene_scale = glm::length(pos_max - pos_min);

		vzm::CameraParameters cam_params;
		vzm::GetCameraParameters(base_scene_id, cam_params, 0);

		if (getMouseWheelDelta(flags) > 0)
			__cv3__ cam_params.pos += scene_scale * 0.01f * (__cv3__ cam_params.view);
		else
			__cv3__ cam_params.pos -= scene_scale * 0.01f * (__cv3__ cam_params.view);

		for (auto& it : scene_name)
		{
			vzm::SetCameraParameters(base_scene_id + it.first, cam_params);
		}
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
				vzm::SetCameraParameters(base_scene_id + it.first, cam_params);
			}
		}
	}
}

void CallBackFunc_SingleStage_Mouse(int event, int x, int y, int flags, void* userdata)
{
	using namespace cv;

	static helpers::arcball aball_vr;
	if (event == EVENT_LBUTTONDOWN || event == EVENT_RBUTTONDOWN)
	{
		std::pair<std::string, int> scene_pair = *(std::pair<std::string, int>*)userdata;
		vzm::CameraParameters cam_params;
		vzm::GetCameraParameters(scene_pair.second, cam_params, 0);

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
	else if (event == EVENT_MOUSEWHEEL && (flags & EVENT_FLAG_CTRLKEY))
	{
		std::pair<std::string, int> scene_pair = *(std::pair<std::string, int>*)userdata;
		glm::fvec3 pos_min, pos_max;
		vzm::GetSceneBoundingBox({}, scene_pair.second, (float*)&pos_min, (float*)&pos_max);
		float scene_scale = glm::length(pos_max - pos_min);

		vzm::CameraParameters cam_params;
		vzm::GetCameraParameters(scene_pair.second, cam_params, 0);

		if (getMouseWheelDelta(flags) > 0)
			__cv3__ cam_params.pos += scene_scale * 0.01f * (__cv3__ cam_params.view);
		else
			__cv3__ cam_params.pos -= scene_scale * 0.01f * (__cv3__ cam_params.view);
		vzm::SetCameraParameters(scene_pair.second, cam_params, 0);
		vzm::RenderScene(scene_pair.second, 0);
	}
	else if (event == EVENT_MOUSEMOVE)
	{
		if ((flags & EVENT_FLAG_LBUTTON) || (flags & EVENT_FLAG_RBUTTON))
		{
			std::pair<std::string, int> scene_pair = *(std::pair<std::string, int>*)userdata;
			vzm::CameraParameters cam_params;
			vzm::GetCameraParameters(scene_pair.second, cam_params, 0);

			helpers::cam_pose arc_cam_pose;
			if (flags & EVENT_FLAG_LBUTTON)
				aball_vr.pan_move((int*)&glm::ivec2(x, y), arc_cam_pose);
			else if (flags & EVENT_FLAG_RBUTTON)
				aball_vr.move((int*)&glm::ivec2(x, y), arc_cam_pose);

			__cv3__ cam_params.pos = __cv3__ arc_cam_pose.pos;
			__cv3__ cam_params.up = __cv3__ arc_cam_pose.up;
			__cv3__ cam_params.view = __cv3__ arc_cam_pose.view;

			vzm::SetCameraParameters(scene_pair.second, cam_params, 0);
			vzm::RenderScene(scene_pair.second, 0);
		}
	}
}

int main()
{
#if (_MSVC_LANG < 201703L)
	# this program uses C++ standard ver 17
#endif

	vzm::InitEngineLib();

	// test rendering example
#define __cv3__ *(glm::fvec3*)
#define __cv4__ *(glm::fvec4*)
#define __cm4__ *(glm::fmat4x4*)
	vzm::CameraParameters cam_params;
	__cv3__ cam_params.pos = glm::fvec3(0, 150, 0);
	__cv3__ cam_params.up = glm::fvec3(0, 0, -1.f);
	__cv3__ cam_params.view = glm::fvec3(0, -1, 0);
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

	auto update_alpha = [](int v, void* params)
	{
		int scene_id = ((int*)params)[0];
		int num_objs = ((int*)params)[1];
		for (int i = 0; i < num_objs; i++)
		{
			int obj_id = ((int*)params)[2 + i];
			vzm::ObjStates objstate;
			vzm::GetSceneObjectState(scene_id, obj_id, objstate);
			objstate.color[3] = v / 100.f;
			vzm::ReplaceOrAddSceneObject(scene_id, obj_id, objstate);
		}
		vzm::RenderScene(scene_id, 0);
	};

	int ct_vol_id = 0;

	// notice : all control parameters and initial parameters can be manipulated while exploring the dataset with interactive 3D rendering.
	// (refer to the video files)

	//vzm::LoadModelFile("..\\data\\data_1.x3d", ct_vol_id);
	//const int thres_value = 41500;
	//const int t_in_init = 6;
	//const int t_out_init = 2;
	//const int g_h_init = 4500;

	vzm::LoadModelFile("..\\data\\data_2.x3d", ct_vol_id);
	glm::fvec3 vol_pitch;
	vzm::GetVolumeInfo(ct_vol_id, NULL, NULL, (float*)&vol_pitch, NULL);
	obj_state.surfel_size = std::max(std::max(vol_pitch.x, vol_pitch.y), vol_pitch.z) * 2.f;

	const int thres_value = 18800;
	const int t_in_init = 8;
	const int t_out_init = 0;
	const int g_h_init = 3600;
	
	int foreground_vol_id = 0;
	vzm::GenerateEmptyVolume(foreground_vol_id, ct_vol_id, "BYTE", 0, 255, 0);;

	std::map<std::string, int> scenes = { {"stage 1", 0}, {"stage 2", 1}, {"stage 3", 2}, {"stage 4", 3} };

	vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { ct_vol_id , foreground_vol_id }, 
		{ {"function_name", std::string("load_foreground")}, {"file_name", std::string("")}, {"iso_value", thres_value}, {"downsacled_factor", -1},
		{"flip_x", false}, {"flip_y", false}, {"flip_z", false} });

	int candidate_pts_id = 0, filtered_vol_id = 0;
	vzm::GenerateEmptyVolume(filtered_vol_id, ct_vol_id, "USHORT", 0, 65535, 0);
	vzm::GenerateEmptyPrimitive(candidate_pts_id);

	int t_in = t_in_init; // control parameter (vox)
	int t_out = t_out_init; // control parameter (vox)
	const float min_ct_v = 0.f, max_ct_v = 65535.f; // using full range of CT values... these can be used as control parameters as well.
	{
		// interactive MPR for selecting t_in and t_out
	}
	vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { candidate_pts_id , filtered_vol_id, ct_vol_id, foreground_vol_id },
		{ {"function_name", std::string("processing_stage1")}, {"t_in", t_in}, {"t_out", t_out}, {"min_sample_v", min_ct_v}, {"max_sample_v", max_ct_v},
		{"gaussian_kernel_width", 5}, {"simplify_grid_length", 1.f}, {"geometric_complexity_kernel_ratio", 0.01f} });

	// w/ interactive point set rendering
	float g_h = g_h_init, mu_u = 0.3f;
	{
		vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { candidate_pts_id },
			{ {"function_name", std::string("update_strongweak_points")}, {"g_h", g_h}, {"mu_u", mu_u} });

		auto it = scenes.find("stage 2");
		vzm::SetSceneEnvParameters(it->second, scn_env_params);
		vzm::SetCameraParameters(it->second, cam_params, 0);
		vzm::ReplaceOrAddSceneObject(it->second, candidate_pts_id, obj_state);

		cv::namedWindow(it->first, cv::WINDOW_AUTOSIZE);
		static std::pair<std::string, int> scene_pair = { it->first, it->second };
		cv::setMouseCallback(it->first, CallBackFunc_SingleStage_Mouse, &scene_pair);

		struct param_info
		{
			int candidate_pts_id;
			float g_h;
			float mu_u;
			std::pair<std::string, int> scene_pair;
		};
		float g_h = g_h_init; // control parameter (ct value)
		float mu_u = 0.3f; // control parameter (0.0 ~ 1.0)

		param_info pram_info_v;
		pram_info_v.g_h = g_h;
		pram_info_v.mu_u = mu_u;
		pram_info_v.candidate_pts_id = candidate_pts_id;
		pram_info_v.scene_pair = scene_pair;

#define CALL_MODULE_CANDIDATE {vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { pinfo->candidate_pts_id },\
		{ {"function_name", std::string("update_strongweak_points")}, { "g_h", pinfo->g_h }, { "mu_u", pinfo->mu_u } }); \
		vzm::RenderScene(pinfo->scene_pair.second, 0); \
		}

		auto update_g_h = [](int v, void* params)
		{
			param_info* pinfo = (param_info*)params;
			pinfo->g_h = (float)v;
			CALL_MODULE_CANDIDATE;
		};
		auto update_mu_u = [](int v, void* params)
		{
			param_info* pinfo = (param_info*)params;
			pinfo->mu_u = v / 100.f;
			CALL_MODULE_CANDIDATE;
		};
		int control_g_h = g_h;
		cv::createTrackbar("g_h", it->first, &control_g_h, 15000, update_g_h, &pram_info_v);
		int control_mu = 30;
		cv::createTrackbar("mu_h x100", it->first, &control_mu, 100, update_mu_u, &pram_info_v);
		int control_alpha = 100;
		int opacity_param_info[3] = { it->second, 1, candidate_pts_id };
		cv::createTrackbar("opacity x100", it->first, &control_alpha, 100, update_alpha, opacity_param_info);

		show_window(it->first, it->second, 0);

		int key = -1;
		while (key != 'n')
		{
			unsigned char* ptr_rgba;
			float* ptr_zdepth;
			int w, h;
			size_t render_count;
			if (vzm::GetRenderBufferPtrs(it->second, &ptr_rgba, &ptr_zdepth, &w, &h, 0, &render_count))
			{
				cv::Mat cvmat(h, w, CV_8UC4, ptr_rgba);
				static size_t prev_count = 0;
				if (prev_count != render_count)
				{
					cv::putText(cvmat, "press 'q' ==> exit program", cv::Point(10, 20), cv::FONT_HERSHEY_DUPLEX, 0.7, CV_RGB(255, 100, 125), 1, cv::LineTypes::LINE_AA);
					cv::putText(cvmat, "press 'n' ==> go to next stage", cv::Point(10, 50), cv::FONT_HERSHEY_DUPLEX, 0.7, CV_RGB(255, 100, 125), 1, cv::LineTypes::LINE_AA);
					prev_count = render_count;
				}
				cv::imshow(it->first, cvmat);
			}
			if (key == 'q')
			{
				vzm::DeinitEngineLib();
				return 0;
			}

			key = cv::waitKey(1);
		}
		cv::destroyWindow(it->first);
		g_h = pram_info_v.g_h;
		mu_u = pram_info_v.mu_u;
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

		auto prev_it = scenes.find("stage 2");
		auto it = scenes.find("stage 3");
		vzm::SetSceneEnvParameters(it->second, scn_env_params);

		vzm::GetCameraParameters(prev_it->second, cam_params, 0);
		vzm::SetCameraParameters(it->second, cam_params, 0);

		vzm::ReplaceOrAddSceneObject(it->second, reliable_pts_id, obj_state);
		vzm::ReplaceOrAddSceneObject(it->second, holefilling_pts_id, obj_state);

		cv::namedWindow(it->first, cv::WINDOW_AUTOSIZE);
		std::pair<std::string, int> scene_pair = { it->first, it->second };
		cv::setMouseCallback(it->first, CallBackFunc_SingleStage_Mouse, &scene_pair);

		struct param_info
		{
			int holefilling_pts_id;
			bool show_nc;
			float m;
			int eta;
			std::pair<std::string, int> scene_pair;
		};

		param_info pram_info_v;
		pram_info_v.holefilling_pts_id = holefilling_pts_id;
		pram_info_v.show_nc = show_nc;
		pram_info_v.m = m;
		pram_info_v.eta = eta;
		pram_info_v.scene_pair = scene_pair;

		auto update_eta = [](int v, void* params)
		{
			param_info* pinfo = (param_info*)params;
			pinfo->eta = v;
			vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { pinfo->holefilling_pts_id },
				{ {"function_name", std::string("update_holefill_points")}, {"eta", pinfo->eta}, {"m", pinfo->m}, {"show_nc", pinfo->show_nc} });
			vzm::RenderScene(pinfo->scene_pair.second, 0);
		};

		auto update_m = [](int v, void* params)
		{
			param_info* pinfo = (param_info*)params;
			pinfo->m = (float)v;
			vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { pinfo->holefilling_pts_id },
				{ {"function_name", std::string("update_holefill_points")}, {"eta", pinfo->eta}, {"m", pinfo->m}, {"show_nc", pinfo->show_nc} });
			vzm::RenderScene(pinfo->scene_pair.second, 0);
		};

		int control_eta = eta;
		cv::createTrackbar("eta", it->first, &control_eta, 500, update_eta, &pram_info_v);
		int control_m = m;
		cv::createTrackbar("m", it->first, &control_m, 10, update_m, &pram_info_v);
		int control_alpha = 100;
		int opacity_param_info[4] = { it->second, 2, reliable_pts_id, holefilling_pts_id };
		cv::createTrackbar("opacity x100", it->first, &control_alpha, 100, update_alpha, opacity_param_info);

		show_window(it->first, it->second, 0);

		int key = -1;
		while (key != 'n')
		{
			unsigned char* ptr_rgba;
			float* ptr_zdepth;
			int w, h;
			size_t render_count;
			if (vzm::GetRenderBufferPtrs(it->second, &ptr_rgba, &ptr_zdepth, &w, &h, 0, &render_count))
			{
				cv::Mat cvmat(h, w, CV_8UC4, ptr_rgba);

				static size_t prev_count = 0;
				if (prev_count != render_count)
				{
					cv::putText(cvmat, "press 'q' ==> exit program", cv::Point(10, 20), cv::FONT_HERSHEY_DUPLEX, 0.7, CV_RGB(255, 100, 125), 1, cv::LineTypes::LINE_AA);
					cv::putText(cvmat, "press 'n' ==> go to next stage", cv::Point(10, 50), cv::FONT_HERSHEY_DUPLEX, 0.7, CV_RGB(255, 100, 125), 1, cv::LineTypes::LINE_AA);
					cv::putText(cvmat, "press 'c' ==> " + std::string(pram_info_v.show_nc ? "show confidence" : "show noise-shaped clusters"), cv::Point(10, 80), cv::FONT_HERSHEY_DUPLEX, 0.7, CV_RGB(255, 100, 125), 1, cv::LineTypes::LINE_AA);
					prev_count = render_count;
				}

				cv::imshow(it->first, cvmat);
			}

			if (key == 'c')
			{
				pram_info_v.show_nc = !pram_info_v.show_nc;
				vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { pram_info_v.holefilling_pts_id },
					{ {"function_name", std::string("update_holefill_points")}, {"eta", pram_info_v.eta}, {"m", pram_info_v.m}, {"show_nc", pram_info_v.show_nc} });
				vzm::RenderScene(pram_info_v.scene_pair.second, 0);
			}

			if (key == 'q')
			{
				vzm::DeinitEngineLib();
				return 0;
			}

			key = cv::waitKey(1);
		}
		cv::destroyWindow(it->first);
		eta = pram_info_v.eta;
		m = pram_info_v.m;
		show_nc = pram_info_v.show_nc;
	}
	
	int ls_reliable_pts_id = 0, ls_holefilling_pts_id = 0;
	vzm::GenerateEmptyPrimitive(ls_reliable_pts_id);
	vzm::GenerateEmptyPrimitive(ls_holefilling_pts_id);
	// w/ semi-interactive point set rendering
	int epsilon_s = 10; // control parameter (to get local isovalues)
	int ot_lev = 0; /*0 or less means automatic decision for octree level, unless, set 5 to 11 (int)*/
	{
		vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { ls_reliable_pts_id, ls_holefilling_pts_id, ct_vol_id, filtered_vol_id, reliable_pts_id, holefilling_pts_id },
			{ {"function_name", std::string("localisosurface_points")}, {"epsilon_s", epsilon_s}, {"eta", eta} });

		auto prev_it = scenes.find("stage 3");
		auto it = scenes.find("stage 4");
		vzm::SetSceneEnvParameters(it->second, scn_env_params);

		vzm::GetCameraParameters(prev_it->second, cam_params, 0);
		vzm::SetCameraParameters(it->second, cam_params, 0);

		vzm::ReplaceOrAddSceneObject(it->second, ls_reliable_pts_id, obj_state);
		vzm::ReplaceOrAddSceneObject(it->second, ls_holefilling_pts_id, obj_state);

		cv::namedWindow(it->first, cv::WINDOW_AUTOSIZE);
		std::pair<std::string, int> scene_pair = { it->first, it->second };
		cv::setMouseCallback(it->first, CallBackFunc_SingleStage_Mouse, &scene_pair);

		struct param_info
		{
			int ls_reliable_pts_id;
			int ls_holefilling_pts_id;
			int ct_vol_id;
			int filtered_vol_id;
			int reliable_pts_id;
			int holefilling_pts_id;
			int epsilon_s;
			int eta;
			std::pair<std::string, int> scene_pair;
		};

		param_info pram_info_v;
		pram_info_v.ls_reliable_pts_id = ls_reliable_pts_id;
		pram_info_v.ls_holefilling_pts_id = ls_holefilling_pts_id;
		pram_info_v.ct_vol_id = ct_vol_id;
		pram_info_v.filtered_vol_id = filtered_vol_id;
		pram_info_v.reliable_pts_id = reliable_pts_id;
		pram_info_v.holefilling_pts_id = holefilling_pts_id;
		pram_info_v.epsilon_s = epsilon_s;
		pram_info_v.eta = eta;
		pram_info_v.scene_pair = scene_pair;

		auto update_epsilon_s = [](int v, void* params)
		{
			param_info* pinfo = (param_info*)params;
			pinfo->epsilon_s = v;
			vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { pinfo->ls_reliable_pts_id, pinfo->ls_holefilling_pts_id, \
				pinfo->ct_vol_id, pinfo->filtered_vol_id, pinfo->reliable_pts_id, pinfo->holefilling_pts_id },
				{ {"function_name", std::string("localisosurface_points")}, {"epsilon_s", pinfo->epsilon_s}, {"eta", pinfo->eta} });
			vzm::RenderScene(pinfo->scene_pair.second, 0);
		};
		int control_epsilon_s = epsilon_s;
		cv::createTrackbar("eta", it->first, &control_epsilon_s, 20, update_epsilon_s, &pram_info_v);
		cv::createTrackbar("lev", it->first, &ot_lev, 12, NULL, &pram_info_v);
		int control_alpha = 100;
		int opacity_param_info[4] = { it->second, 2, ls_reliable_pts_id, ls_holefilling_pts_id };
		cv::createTrackbar("opacity x100", it->first, &control_alpha, 100, update_alpha, opacity_param_info);

		show_window(it->first, it->second, 0);

		int key = -1;
		while (key != 'n')
		{
			unsigned char* ptr_rgba;
			float* ptr_zdepth;
			int w, h;
			size_t render_count;
			if (vzm::GetRenderBufferPtrs(it->second, &ptr_rgba, &ptr_zdepth, &w, &h, 0, &render_count))
			{
				cv::Mat cvmat(h, w, CV_8UC4, ptr_rgba);
				static size_t prev_count = 0;
				if (prev_count != render_count)
				{
					cv::putText(cvmat, "press 'q' ==> exit program", cv::Point(10, 20), cv::FONT_HERSHEY_DUPLEX, 0.7, CV_RGB(255, 100, 125), 1, cv::LineTypes::LINE_AA);
					cv::putText(cvmat, "press 'n' ==> go to next stage", cv::Point(10, 50), cv::FONT_HERSHEY_DUPLEX, 0.7, CV_RGB(255, 100, 125), 1, cv::LineTypes::LINE_AA);
					cv::putText(cvmat, "if octree level is 0, automatically set an optimal fitting level", cv::Point(10, 80), cv::FONT_HERSHEY_DUPLEX, 0.5, CV_RGB(115, 100, 200), 1, cv::LineTypes::LINE_AA);
					prev_count = render_count;
				}

				cv::imshow(it->first, cvmat);
			}

			key = cv::waitKey(1);
		}
		cv::destroyWindow(it->first);
		epsilon_s = pram_info_v.epsilon_s;
	}

	int final_mesh_id = 0;
	vzm::GenerateEmptyPrimitive(final_mesh_id);
	vzm::ExecuteModule2("ct_modeler.dll", "function_launcher", { final_mesh_id, ls_reliable_pts_id, ls_holefilling_pts_id, holefilling_pts_id },
		{ {"function_name", std::string("processing_final")}, {"eta", eta}, {"m", m}, {"otlev", ot_lev } });


	scene_name[0] = "candidate points";
	scene_name[1] = "confidence points";
	scene_name[2] = "final mesh";
	std::map<std::string, std::list<int>> prim_ids = { 
		{scene_name[0], {candidate_pts_id}} ,
		{scene_name[1], {reliable_pts_id, holefilling_pts_id}} ,
		{scene_name[2], {final_mesh_id}}
	};

	const int base_scn_id = 10;
	for (int i = 0; i < (int)scene_name.size(); i++)
	{
		vzm::SetSceneEnvParameters(base_scn_id + i, scn_env_params);
		vzm::SetCameraParameters(base_scn_id + i, cam_params, 0);

		auto it = prim_ids.find(scene_name[i]);
		for (int& obj_id : it->second)
			vzm::ReplaceOrAddSceneObject(base_scn_id + i, obj_id, obj_state);

		cv::namedWindow(scene_name[i], cv::WINDOW_AUTOSIZE);
		cv::setMouseCallback(scene_name[i], CallBackFunc_Mouse, NULL);// &scenes[i]);
		show_window(scene_name[i], base_scn_id + i, 0);
	}

	int key = -1;
	while (key != 'q')
	{
		for (auto& it : scene_name)
		{
			show_window(it.second, base_scn_id + it.first, 0);
		}
		key = cv::waitKey(1);
	}

	vzm::DeinitEngineLib();
	return 0;
}