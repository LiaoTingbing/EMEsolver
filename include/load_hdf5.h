#pragma once

#include "common.h"

//	加载hdf5文件数据
void loadHdf5Data(map<string, cube>& dev,
	const string filePath,
	const field<string>& dataSetName);