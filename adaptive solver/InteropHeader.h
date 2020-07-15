#pragma once

#define IH_VERSION 0.1.190220

#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <windows.h>
#include <time.h>
#include <map>
#include <set>

#define _cst_vx *(vxfloat3*)&
#define _icst_vx *(__float3*)&
#define _f3_(x, y, z) (float)(x), (float)(y), (float)(z)
#define _i3_(x, y, z) (int)(x), (int)(y), (int)(z)
#define _fp_ (float*)

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

typedef struct __int3__
{
	int x, y, z;
	__int3__() { x = y = z = 0; };
	__int3__(float _x, float _y, float _z) { x = (int)_x; y = (int)_y; z = (int)_z; };
	__int3__(int _x, int _y, int _z) { x = _x; y = _y; z = _z; };
}  __int3;

typedef struct __int2__
{
	int x, y;
	__int2__() { x = y = 0; };
	__int2__(float _x, float _y) { x = (int)_x; y = (int)_y; };
	__int2__(int _x, int _y) { x = _x; y = _y; };
}  __int2;

//typedef std::unordered_map<const __int3, data, key_hash, key_equal> map_t;

template <typename T> void convert_type_stdvector(std::vector<__float3>& _pts, const std::vector<T>& pts)
{
	_pts.assign(pts.size(), __float3());
	for (int i = 0; i < (int)pts.size(); i++)
		_pts[i] = *(__float3*)&pts[i];
}

template <typename TIN, typename TOUT> void convert_type_stdvector(std::vector<TOUT>& _pts, const std::vector<TIN>& pts)
{
	_pts.assign(pts.size(), TOUT());
	for (int i = 0; i < (int)pts.size(); i++)
		_pts[i] = *(TOUT*)&pts[i];
}
//typedef 

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

enum DataType {
	__DataTypeUNDEFINED = 0,/*!< No Flip */
	__DataTypeBYTE,/*!< 1byte, unsigned char, defined as byte */
	__DataTypeUNSIGNEDSHORT,/*!< 2bytes, unsigned short, defined as ushort*/
	__DataTypeUNSIGNEDINT,/*!< 4bytes, unsigned int, defined as uint*/
	__DataTypeFLOAT,/*!< 4bytes, float, defined as float*/
	__DataTypeDOUBLE,/*!< 8bytes, float, defined as double*/
};

#include "../nanoflann.hpp"

template <typename T, typename TT>
struct PointCloud
{
	//public:
	const std::vector<TT>& pts;
	PointCloud(const std::vector<TT>& _pts) : pts(_pts) { }

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline T kdtree_distance(const T *p1, const size_t idx_p2, size_t) const
	{
		int dim = (int)sizeof(TT) / (int)sizeof(T);
		T d_sq = 0;
		for (int i = 0; i < dim; i++)
		{
			const T d = p1[i] - ((const T*)&(pts[idx_p2]))[i];
			d_sq += d * d;
		}
		return d_sq;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		return ((const T*)&(pts[idx]))[dim];
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const { return false; }
};

static double __kd_tree_time = 0;
template <typename Real, typename T, unsigned int DIM>
struct point_info
{
	typedef nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<Real, PointCloud<Real, T> >,
		PointCloud<Real, T>,
		DIM // dim 
	> kd_tree;

	PointCloud<Real, T>* pc_kdt;
	kd_tree* kdt;

	std::vector<T> pos_pts;
	std::vector<T> nrl_pts;
	std::vector<T> lfn_pts; // local frame normal
	std::vector<Real> iv_pts; // image value
	std::vector<Real> cf_pts; // confidence
	std::vector<Real> fc_pts; // tmp : flow-curvature-based uncertainty
	std::vector<Real> pv_pts; // tmp : pca-based points variance
	std::vector<Real> sv_pts; // shape variation
	std::vector<Real> ri_pts; // regional info. for confidence weight

	point_info() {
		pc_kdt = NULL;
		kdt = NULL;
	}
	~point_info() {
		clear();
	}

	void build_kdt(const std::string& ment)
	{
		//DWORD _time = timeGetTime();
		build_kdt();
		//std::cout << ment.c_str() << timeGetTime() - _time << " ms" << std::endl;
	}
	void build_kdt()
	{
		if (pos_pts.size() == 0) return;
		if (kdt) {
			delete pc_kdt, delete kdt;
		}

		//DWORD _time = timeGetTime();
		pc_kdt = new PointCloud<Real, T>(pos_pts);
		kdt = new kd_tree(DIM, *pc_kdt, nanoflann::KDTreeSingleIndexAdaptorParams(10));
		kdt->buildIndex();
		//__kd_tree_time += timeGetTime() - _time;
	};
	void clear()
	{
		pos_pts.clear();
		nrl_pts.clear();
		lfn_pts.clear();
		iv_pts.clear(); // ct value
		//gm_pts.clear(); // grad mag
		cf_pts.clear(); // confidence
		fc_pts.clear(); // flow curvature
		pv_pts.clear();
		sv_pts.clear(); // local shape variation
		ri_pts.clear();
		if (pc_kdt) delete pc_kdt;
		if (kdt) delete kdt;
		pc_kdt = NULL;
		kdt = NULL;
	};
	void add_pt(const point_info& src_pts, int idx)
	{
		if (idx < src_pts.pos_pts.size()) pos_pts.push_back(src_pts.pos_pts[idx]);
		if (idx < src_pts.nrl_pts.size()) nrl_pts.push_back(src_pts.nrl_pts[idx]);
		if (idx < src_pts.lfn_pts.size()) lfn_pts.push_back(src_pts.lfn_pts[idx]);
		if (idx < src_pts.iv_pts.size())  iv_pts.push_back(src_pts.iv_pts[idx]);
		//if (idx < src_pts.gm_pts.size())  gm_pts.push_back(src_pts.gm_pts[idx]);
		if (idx < src_pts.cf_pts.size())  cf_pts.push_back(src_pts.cf_pts[idx]);
		if (idx < src_pts.fc_pts.size())  fc_pts.push_back(src_pts.fc_pts[idx]);
		if (idx < src_pts.pv_pts.size())  pv_pts.push_back(src_pts.pv_pts[idx]);
		if (idx < src_pts.sv_pts.size())  sv_pts.push_back(src_pts.sv_pts[idx]);
		if (idx < src_pts.ri_pts.size())  ri_pts.push_back(src_pts.ri_pts[idx]);
	}

	void insert_back(const point_info& src_pts)
	{
		if (src_pts.pos_pts.size() > 0) pos_pts.insert(pos_pts.end(), src_pts.pos_pts.begin(), src_pts.pos_pts.end());
		if (src_pts.nrl_pts.size() > 0) nrl_pts.insert(nrl_pts.end(), src_pts.nrl_pts.begin(), src_pts.nrl_pts.end());
		if (src_pts.lfn_pts.size() > 0) lfn_pts.insert(lfn_pts.end(), src_pts.lfn_pts.begin(), src_pts.lfn_pts.end());
		if (src_pts.iv_pts.size() > 0)  iv_pts.insert(iv_pts.end(), src_pts.iv_pts.begin(), src_pts.iv_pts.end());
		//if (src_pts.gm_pts.size() > 0)  gm_pts.insert(gm_pts.end(), src_pts.gm_pts.begin(), src_pts.gm_pts.end());
		if (src_pts.cf_pts.size() > 0)  cf_pts.insert(cf_pts.end(), src_pts.cf_pts.begin(), src_pts.cf_pts.end());
		if (src_pts.fc_pts.size() > 0)  fc_pts.insert(fc_pts.end(), src_pts.fc_pts.begin(), src_pts.fc_pts.end());
		if (src_pts.pv_pts.size() > 0)  pv_pts.insert(pv_pts.end(), src_pts.pv_pts.begin(), src_pts.pv_pts.end());
		if (src_pts.sv_pts.size() > 0)  sv_pts.insert(sv_pts.end(), src_pts.sv_pts.begin(), src_pts.sv_pts.end());
		if (src_pts.ri_pts.size() > 0)  ri_pts.insert(ri_pts.end(), src_pts.ri_pts.begin(), src_pts.ri_pts.end());
	}

	void copy(const point_info& src_pts)
	{
		pos_pts = src_pts.pos_pts;
		nrl_pts = src_pts.nrl_pts;
		lfn_pts = src_pts.lfn_pts;
		iv_pts = src_pts.iv_pts;
		//gm_pts = src_pts.gm_pts;
		cf_pts = src_pts.cf_pts;
		fc_pts = src_pts.fc_pts;
		pv_pts = src_pts.pv_pts;
		sv_pts = src_pts.sv_pts;
		ri_pts = src_pts.ri_pts;
	}

	size_t size() const { return pos_pts.size(); }

	point_info& operator = (const point_info& src_pts) {
		clear();
		this->pos_pts = src_pts.pos_pts;
		this->nrl_pts = src_pts.nrl_pts;
		this->lfn_pts = src_pts.lfn_pts;
		this->iv_pts = src_pts.iv_pts;
		//this->gm_pts = src_pts.gm_pts;
		this->cf_pts = src_pts.cf_pts;
		this->fc_pts = src_pts.fc_pts;
		this->pv_pts = src_pts.pv_pts;
		this->sv_pts = src_pts.sv_pts;
		this->ri_pts = src_pts.ri_pts;
		std::cout << "copy call in point_info" << std::endl;
		return *this;
	};
};

struct debug_info
{
	// lamda functions
	std::function< void(const std::vector<__float3>&, const __float3&, std::wstring) > ___show_test_pts_with_global_color;
	std::function< void(const std::vector<__float3>&, const std::vector<__float3>&, std::wstring) > ___show_test_pts_with_color_map;
	std::function< std::vector<__float3>(const std::map<int, std::vector<int>>&, const std::set<int>&, const int, const bool) > ___get_colormap_clusters_pts;
};

#pragma region // RGB TABLE
#define __maroon					0x800000
#define __darkred					0x8B0000
#define __brown					0xA52A2A
#define __firebrick				0xB22222
#define __crimson					0xDC143C
#define __red						0xFF0000
#define __tomato					0xFF6347
#define __coral					0xFF7F50
#define __indianred				0xCD5C5C
#define __lightcoral				0xF08080
#define __darksalmon				0xE9967A
#define __salmon					0xFA8072
#define __lightsalmon				0xFFA07A
#define __orangered				0xFF4500
#define __darkorange				0xFF8C00
#define __orange					0xFFA500
#define __gold					0xFFD700
#define __darkgoldenrod			0xB8860B
#define __goldenrod				0xDAA520
#define __palegoldenrod			0xEEE8AA
#define __darkkhaki				0xBDB76B
#define __khaki					0xF0E68C
#define __olive					0x808000
#define __yellow					0xFFFF00
#define __yellowgreen				0x9ACD32
#define __darkolivegreen			0x556B2F
#define __olivedrab				0x6B8E23
#define __lawngreen				0x7CFC00
#define __chartreuse				0x7FFF00
#define __greenyellow				0xADFF2F
#define __darkgreen				0x006400
#define __green					0x008000
#define __forestgreen				0x228B22
#define __lime					0x00FF00
#define __limegreen				0x32CD32
#define __lightgreen				0x90EE90
#define __palegreen				0x98FB98
#define __darkseagreen			0x8FBC8F
#define __mediumspringgreen		0x00FA9A
#define __springgreen				0x00FF7F
#define __seagreen				0x2E8B57
#define __mediumaquamarine		0x66CDAA
#define __mediumseagreen			0x3CB371
#define __lightseagreen			0x20B2AA
#define __darkslategray			0x2F4F4F
#define __teal					0x008080
#define __darkcyan				0x008B8B
#define __aqua					0x00FFFF
#define __cyan					0x00FFFF
#define __lightcyan				0xE0FFFF
#define __darkturquoise			0x00CED1
#define __turquoise				0x40E0D0
#define __mediumturquoise			0x48D1CC
#define __paleturquoise			0xAFEEEE
#define __aquamarine				0x7FFFD4
#define __powderblue				0xB0E0E6
#define __cadetblue				0x5F9EA0
#define __steelblue				0x4682B4
#define __cornflowerblue			0x6495ED
#define __deepskyblue				0x00BFFF
#define __dodgerblue				0x1E90FF
#define __lightblue				0xADD8E6
#define __skyblue					0x87CEEB
#define __lightskyblue			0x87CEFA
#define __midnightblue			0x191970
#define __navy					0x000080
#define __darkblue				0x00008B
#define __mediumblue				0x0000CD
#define __blue					0x0000FF
#define __royalblue				0x4169E1
#define __blueviolet				0x8A2BE2
#define __indigo					0x4B0082
#define __darkslateblue			0x483D8B
#define __slateblue				0x6A5ACD
#define __mediumslateblue			0x7B68EE
#define __mediumpurple			0x9370DB
#define __darkmagenta				0x8B008B
#define __darkviolet				0x9400D3
#define __darkorchid				0x9932CC
#define __mediumorchid			0xBA55D3
#define __purple					0x800080
#define __thistle					0xD8BFD8
#define __plum					0xDDA0DD
#define __violet					0xEE82EE
#define __magentafuchsia			0xFF00FF
#define __orchid					0xDA70D6
#define __mediumvioletred			0xC71585
#define __palevioletred			0xDB7093
#define __deeppink				0xFF1493
#define __hotpink					0xFF69B4
#define __lightpink				0xFFB6C1
#define __pink					0xFFC0CB
#define __antiquewhite			0xFAEBD7
#define __beige					0xF5F5DC
#define __bisque					0xFFE4C4
#define __blanchedalmond			0xFFEBCD
#define __wheat					0xF5DEB3
#define __cornsilk				0xFFF8DC
#define __lemonchiffon			0xFFFACD
#define __lightgoldenrodyellow	0xFAFAD2
#define __lightyellow				0xFFFFE0
#define __saddlebrown				0x8B4513
#define __sienna					0xA0522D
#define __chocolate				0xD2691E
#define __peru					0xCD853F
#define __sandybrown				0xF4A460
#define __burlywood				0xDEB887
#define __tan						0xD2B48C
#define __rosybrown				0xBC8F8F
#define __moccasin				0xFFE4B5
#define __navajowhite				0xFFDEAD
#define __peachpuff				0xFFDAB9
#define __mistyrose				0xFFE4E1
#define __lavenderblush			0xFFF0F5
#define __linen					0xFAF0E6
#define __oldlace					0xFDF5E6
#define __papayawhip				0xFFEFD5
#define __seashell				0xFFF5EE
#define __mintcream				0xF5FFFA
#define __slategray				0x708090
#define __lightslategray			0x778899
#define __lightsteelblue			0xB0C4DE
#define __lavender				0xE6E6FA
#define __floralwhite				0xFFFAF0
#define __aliceblue				0xF0F8FF
#define __ghostwhite				0xF8F8FF
#define __honeydew				0xF0FFF0
#define __ivory					0xFFFFF0
#define __azure					0xF0FFFF
#define __snow					0xFFFAFA
#pragma endregion // RGB TABLE

enum __color_name
{
	_color_maroon = 0
	//, _color_darkred
	//, _color_brown
	, _color_firebrick
	, _color_crimson
	, _color_red
	, _color_tomato
	//, _color_coral
	//, _color_indianred
	//, _color_lightcoral
	//, _color_darksalmon
	, _color_salmon
	, _color_lightsalmon
	, _color_orangered
	//, _color_darkorange
	, _color_orange
	, _color_gold
	, _color_darkgoldenrod
	, _color_goldenrod
	//, _color_palegoldenrod
	, _color_darkkhaki
	, _color_khaki
	, _color_olive
	, _color_yellow
	, _color_yellowgreen
	//, _color_darkolivegreen
	//, _color_olivedrab
	//, _color_lawngreen
	//, _color_chartreuse
	, _color_greenyellow
	, _color_darkgreen
	, _color_green
	, _color_forestgreen
	, _color_lime
	, _color_limegreen
	//, _color_lightgreen
	, _color_palegreen
	, _color_darkseagreen
	, _color_mediumspringgreen
	//, _color_springgreen
	//, _color_seagreen
	, _color_mediumaquamarine
	//, _color_mediumseagreen
	, _color_lightseagreen
	, _color_darkslategray
	//, _color_teal
	, _color_darkcyan
	//, _color_aqua
	, _color_cyan
	//, _color_lightcyan
	, _color_darkturquoise
	//, _color_turquoise
	//, _color_mediumturquoise
	//, _color_paleturquoise
	, _color_aquamarine
	//, _color_powderblue
	, _color_cadetblue
	, _color_steelblue
	, _color_cornflowerblue
	, _color_deepskyblue
	, _color_dodgerblue
	//, _color_lightblue
	, _color_skyblue
	, _color_lightskyblue
	//, _color_midnightblue
	, _color_navy
	//, _color_darkblue
	, _color_mediumblue
	//, _color_blue
	, _color_royalblue
	, _color_blueviolet
	, _color_indigo
	, _color_darkslateblue
	, _color_slateblue
	//, _color_mediumslateblue
	, _color_mediumpurple
	, _color_darkmagenta
	, _color_darkviolet
	//, _color_darkorchid
	, _color_mediumorchid
	, _color_purple
	//, _color_thistle
	//, _color_plum
	//, _color_violet
	, _color_magentafuchsia
	//, _color_orchid
	//, _color_mediumvioletred
	//, _color_palevioletred
	, _color_deeppink
	, _color_hotpink
	//, _color_lightpink
	//, _color_pink
	//, _color_antiquewhite
	//, _color_beige
	//, _color_bisque
	//, _color_blanchedalmond
	//, _color_wheat
	//, _color_cornsilk
	//, _color_lemonchiffon
	//, _color_lightgoldenrodyellow
	//, _color_lightyellow
	//, _color_saddlebrown
	, _color_sienna
	, _color_chocolate
	, _color_peru
	//, _color_sandybrown
	//, _color_burlywood
	//, _color_tan
	, _color_rosybrown
	//, _color_moccasin
	//, _color_navajowhite
	//, _color_peachpuff
	//, _color_mistyrose
	//, _color_lavenderblush
	//, _color_linen
	//, _color_oldlace
	//, _color_papayawhip
	//, _color_seashell
	//, _color_mintcream
	, _color_slategray
	//, _color_lightslategray
	//, _color_lightsteelblue
	//, _color_lavender
	//, _color_floralwhite
	//, _color_aliceblue
	//, _color_ghostwhite
	//, _color_honeydew
	//, _color_ivory
	//, _color_azure
	//, _color_snow
};

static int __COLORS__[]{
	__maroon
	//, __darkred
	//, __brown
	, __firebrick
	, __crimson
	//, __red
	, __tomato
	//, __coral
	//, __indianred
	//, __lightcoral
	//, __darksalmon
	, __salmon
	, __lightsalmon
	, __orangered
	//, __darkorange
	, __orange
	, __gold
	, __darkgoldenrod
	, __goldenrod
	//, __palegoldenrod
	, __darkkhaki
	, __khaki
	, __olive
	, __yellow
	, __yellowgreen
	//, __darkolivegreen
	//, __olivedrab
	//, __lawngreen
	//, __chartreuse
	, __greenyellow
	, __darkgreen
	//, __green
	, __forestgreen
	, __lime
	, __limegreen
	//, __lightgreen
	, __palegreen
	, __darkseagreen
	, __mediumspringgreen
	//, __springgreen
	//, __seagreen
	, __mediumaquamarine
	//, __mediumseagreen
	, __lightseagreen
	, __darkslategray
	//, __teal
	, __darkcyan
	//, __aqua
	, __cyan
	//, __lightcyan
	, __darkturquoise
	//, __turquoise
	//, __mediumturquoise
	//, __paleturquoise
	, __aquamarine
	//, __powderblue
	, __cadetblue
	, __steelblue
	, __cornflowerblue
	, __deepskyblue
	, __dodgerblue
	//, __lightblue
	//, __skyblue
	, __lightskyblue
	//, __midnightblue
	, __navy
	//, __darkblue
	, __mediumblue
	//, __blue
	, __royalblue
	, __blueviolet
	, __indigo
	, __darkslateblue
	, __slateblue
	//, __mediumslateblue
	, __mediumpurple
	, __darkmagenta
	, __darkviolet
	//, __darkorchid
	, __mediumorchid
	, __purple
	//, __thistle
	//, __plum
	//, __violet
	, __magentafuchsia
	//, __orchid
	//, __mediumvioletred
	//, __palevioletred
	, __deeppink
	, __hotpink
	//, __lightpink
	//, __pink
	//, __antiquewhite
	//, __beige
	//, __bisque
	//, __blanchedalmond
	//, __wheat
	//, __cornsilk
	//, __lemonchiffon
	//, __lightgoldenrodyellow
	//, __lightyellow
	//, __saddlebrown
	, __sienna
	, __chocolate
	, __peru
	//, __sandybrown
	//, __burlywood
	//, __tan
	, __rosybrown
	//, __moccasin
	//, __navajowhite
	//, __peachpuff
	//, __mistyrose
	//, __lavenderblush
	//, __linen
	//, __oldlace
	//, __papayawhip
	//, __seashell
	//, __mintcream
	, __slategray
	//, __lightslategray
	//, __lightsteelblue
	//, __lavender
	//, __floralwhite
	//, __aliceblue
	//, __ghostwhite
	//, __honeydew
	//, __ivory
	//, __azure
	//, __snow
};

#define __COLOR_ARRAY_SIZE 1024
static int _color_array[__COLOR_ARRAY_SIZE];

inline __float3 color_int_2_float3(int _color)
{
	float r = (float)((_color >> 16) & 0xFF) / 255.f;
	float g = (float)((_color >> 8) & 0xFF) / 255.f;
	float b = (float)((_color >> 0) & 0xFF) / 255.f;
	return __float3(r, g, b);
}
inline __float3 color_clridx_2_float3(int _color_idx)
{
	return color_int_2_float3(__COLORS__[_color_idx]);
}

static void Fill_Color_Array(std::vector<std::pair<int, int>> color_idx_ctrs, int ary_size)
{
	__float3 prev_color(color_idx_ctrs[0].first);
	int prev_idx = 0;

	int end_color = color_idx_ctrs[color_idx_ctrs.size() - 1].first;
	int end_idx = ary_size - 1;
	color_idx_ctrs.push_back(std::pair<int, int>(end_color, end_idx));

	memset(_color_array, 0, sizeof(int) * ary_size);

	int idx_array = 0;
	for (int i = 0; i < (int)color_idx_ctrs.size(); i++)
	{
		__float3 next_color(color_idx_ctrs[i].first);
		int next_idx = color_idx_ctrs[i].second;

		while (idx_array <= next_idx && idx_array < ary_size)   //iAddrTable가 0부터 시작임으로 d4NextColor.w가 0이하일 경우는 패스 된다. 
		{
			double interpolate_ratio = 1.0;
			if (next_idx != prev_idx)
				interpolate_ratio = (double)(idx_array - prev_idx) / (double)(next_idx - prev_idx);

			int _R = (int)((interpolate_ratio * next_color.x + (1 - interpolate_ratio) * prev_color.x) * 255.0);
			int _G = (int)((interpolate_ratio * next_color.y + (1 - interpolate_ratio) * prev_color.y) * 255.0);
			int _B = (int)((interpolate_ratio * next_color.z + (1 - interpolate_ratio) * prev_color.z) * 255.0);

			_color_array[idx_array] = (_R << 16) | (_G << 8) | _B;

			idx_array++;
		}
		prev_color = next_color;
		prev_idx = next_idx;
	}
}
