
#include "include/common.h"


map<string, cube> dev;

void loadHdf5()
{
	// 读取x
	cout << "读取文件数据";
	string fileName = "E:/VisualStudioFiles/EMEsolver/lumerical/out.h5";

	int nmodes = 10;
	int ncells = 10;
	int ny = 51;
	int nz = 51;

	field<Cell> Cells(ncells);
	// 加载有效折射率

	mat realNeffs;
	mat imagNeffs;
	realNeffs.load(hdf5_name(fileName, "real_neffs"));
	imagNeffs.load(hdf5_name(fileName, "imag_neffs"));

	// 加载模场
	for (int indexCell = 0; indexCell < ncells;indexCell++)
	{
		// 加载折射率
		Cells(indexCell).neffs = realNeffs.col(indexCell)
			+ imagNeffs.col(indexCell) * IU;

		//	加载坐标轴
		Cells(indexCell).x.load(hdf5_name(fileName, "x"));
		Cells(indexCell).y.load(hdf5_name(fileName, "y"));
		Cells(indexCell).z.load(hdf5_name(fileName, "z"));
		Cells(indexCell).lambda = 1.55e-6;

		//	加载模场
		// Cells 初始化
		Cells(indexCell).modeFields.set_size(nmodes);
		for (int modeIndex = 0; modeIndex < nmodes;modeIndex++)
		{
			cube tmpData;
			string modeRealName = "cell_" + to_string(indexCell + 1) +
				"_real_eh" + to_string(modeIndex + 1);
			//cout << modeReal <<"\n";
			string modeImagName = "cell_" + to_string(indexCell + 1) +
				"_imag_eh" + to_string(modeIndex + 1);

			cube modeReal;
			modeReal.load(hdf5_name(fileName, modeRealName));
			cube modeImag;
			modeImag.load(hdf5_name(fileName, modeImagName));

			// 存储ex
			Cells(indexCell).modeFields(modeIndex).ex
				= modeReal.slice(0) + IU * modeImag.slice(0);
			// 存储ey
			Cells(indexCell).modeFields(modeIndex).ey
				= modeReal.slice(1) + IU * modeImag.slice(1);
			// 存储ez
			Cells(indexCell).modeFields(modeIndex).ez
				= modeReal.slice(2) + IU * modeImag.slice(2);
			// 存储hx
			Cells(indexCell).modeFields(modeIndex).hx
				= modeReal.slice(3) + IU * modeImag.slice(3);
			// 存储hy
			Cells(indexCell).modeFields(modeIndex).hy
				= modeReal.slice(4) + IU * modeImag.slice(4);
			// 存储hz
			Cells(indexCell).modeFields(modeIndex).hz
				= modeReal.slice(5) + IU * modeImag.slice(5);
		}
	}











}


int main() {

	loadHdf5();
	return 0;
}