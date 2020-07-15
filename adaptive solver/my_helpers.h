#pragma once

#include "InteropHeader.h"

const unsigned int g_edgeTable[16] = {
	0, 1, 2, 5,
	3, 6, 8, 4,
	4, 7, 6, 3,
	5, 2, 1, 0
};

template <typename T>
void extract_iso_contour(std::vector<__float2>& lines_list, const float iso_value, const T* pixs, const int w, const int h)
{
	auto __read_pixel = [&](const int x, int y)->float
	{
		if (x < 0 || x >= w || y < 0 || y >= h) return 0;
		return (float)pixs[x + y * w];
	};

	for (int i = 0; i < h - 1; i++)
	{
		for (int j = 0; j < w - 1; j++)
		{
			float value0 = __read_pixel(j, i);
			float value1 = __read_pixel(j + 1, i);
			float value2 = __read_pixel(j, i + 1);
			float value3 = __read_pixel(j + 1, i + 1);

			int tableIndex = 0;
			if (value0 < iso_value) tableIndex |= 0x1;
			if (value1 < iso_value) tableIndex |= 0x2;
			if (value2 < iso_value) tableIndex |= 0x4;
			if (value3 < iso_value) tableIndex |= 0x8;

#define EDGEPOS(A, VAL0, VAL1, ISO, S_POS, E_POS) { \
	float interval = (float)(VAL0 - VAL1); \
	float ratio = (float)(VAL0 - ISO) / interval; \
	A.x = S_POS.x * (1.f - ratio) + E_POS.x * ratio; \
	A.y = S_POS.y * (1.f - ratio) + E_POS.y * ratio; }

			__float2 pos0((float)j, (float)i);
			__float2 pos1((float)(j + 1), (float)i);
			__float2 pos2((float)j, (float)(i + 1));
			__float2 pos3((float)(j + 1), (float)(i + 1));
			// edge 0 : p0,p1
			// edge 1 : p0,p2
			// edge 2 : p1,p3
			// edge 3 : p2,p3
			__float2 pos_edges[4];
			switch (g_edgeTable[tableIndex])
			{
			case 1: // edge 0, 1
				EDGEPOS(pos_edges[0], value0, value1, iso_value, pos0, pos1);
				EDGEPOS(pos_edges[1], value0, value2, iso_value, pos0, pos2);
				lines_list.push_back(pos_edges[0]);
				lines_list.push_back(pos_edges[1]);
				break;
			case 2: // edge 0, 2
				EDGEPOS(pos_edges[0], value0, value1, iso_value, pos0, pos1);
				EDGEPOS(pos_edges[2], value1, value3, iso_value, pos1, pos3);
				lines_list.push_back(pos_edges[0]);
				lines_list.push_back(pos_edges[2]);
				break;
			case 3: // edge 1, 3
				EDGEPOS(pos_edges[1], value0, value2, iso_value, pos0, pos2);
				EDGEPOS(pos_edges[3], value2, value3, iso_value, pos2, pos3);
				lines_list.push_back(pos_edges[1]);
				lines_list.push_back(pos_edges[3]);
				break;
			case 4: // edge 2, 3
				EDGEPOS(pos_edges[2], value1, value3, iso_value, pos1, pos3);
				EDGEPOS(pos_edges[3], value2, value3, iso_value, pos2, pos3);
				lines_list.push_back(pos_edges[2]);
				lines_list.push_back(pos_edges[3]);
				break;
			case 5: // edge 1, 2
				EDGEPOS(pos_edges[1], value0, value2, iso_value, pos0, pos2);
				EDGEPOS(pos_edges[2], value1, value3, iso_value, pos1, pos3);
				lines_list.push_back(pos_edges[1]);
				lines_list.push_back(pos_edges[2]);
				break;
			case 6: // edge 0, 3
				EDGEPOS(pos_edges[0], value0, value1, iso_value, pos0, pos1);
				EDGEPOS(pos_edges[3], value2, value3, iso_value, pos2, pos3);
				lines_list.push_back(pos_edges[0]);
				lines_list.push_back(pos_edges[3]);
				break;
			case 7: // edge 0, 1, 2, 3
				EDGEPOS(pos_edges[0], value0, value1, iso_value, pos0, pos1);
				EDGEPOS(pos_edges[1], value0, value2, iso_value, pos0, pos2);
				EDGEPOS(pos_edges[2], value1, value3, iso_value, pos1, pos3);
				EDGEPOS(pos_edges[3], value2, value3, iso_value, pos2, pos3);
				lines_list.push_back(pos_edges[0]);
				lines_list.push_back(pos_edges[1]);
				lines_list.push_back(pos_edges[2]);
				lines_list.push_back(pos_edges[3]);
				break;
			case 8: // edge 0, 2, 1, 3
				EDGEPOS(pos_edges[0], value0, value1, iso_value, pos0, pos1);
				EDGEPOS(pos_edges[1], value0, value2, iso_value, pos0, pos2);
				EDGEPOS(pos_edges[2], value1, value3, iso_value, pos1, pos3);
				EDGEPOS(pos_edges[3], value2, value3, iso_value, pos2, pos3);
				lines_list.push_back(pos_edges[0]);
				lines_list.push_back(pos_edges[2]);
				lines_list.push_back(pos_edges[1]);
				lines_list.push_back(pos_edges[3]);
				break;
			default: continue;
			}
		}
	}
}