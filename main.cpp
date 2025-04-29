
#include "include/common.h"


map<string, cube> dev;

void loadHdf5()
{
	// 读取x
	cout << "\t读取文件数据\n";
	string fileName = "E:/VisualStudioFiles/EMEsolver/lumerical/out.h5";

	int nmodes = 3;
	int ncells = 10;
	int ny = 51;
	int nz = 61;

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
			//cout << modeReal.n_rows << modeReal.n_cols
			//	<< modeReal.n_slices<<endl;

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

	// 计算归一化电场
	cout << "\t计算归一化电场\n";
	for (int cellIndex = 0;cellIndex < ncells;cellIndex++)
	{
		// 归一化
		for (int modeIndex = 0; modeIndex < nmodes; modeIndex++)
		{
			mat Px = real(
				Cells(cellIndex).modeFields(modeIndex).ey % Cells(cellIndex).modeFields(modeIndex).hz
				- Cells(cellIndex).modeFields(modeIndex).ez % Cells(cellIndex).modeFields(modeIndex).hy
			);
			mat  power = trapz(Cells(cellIndex).z,
				trapz(Cells(cellIndex).y, Px, 0), 1);
			Cells(cellIndex).modeFields(modeIndex).ex /= sqrt(power(0));
			Cells(cellIndex).modeFields(modeIndex).ey /= sqrt(power(0));
			Cells(cellIndex).modeFields(modeIndex).ez /= sqrt(power(0));
			Cells(cellIndex).modeFields(modeIndex).hx /= sqrt(power(0));
			Cells(cellIndex).modeFields(modeIndex).hy /= sqrt(power(0));
			Cells(cellIndex).modeFields(modeIndex).hz /= sqrt(power(0));
		}
	}

	// 计算传播散射矩阵
	cout << "\t计算传播散射矩阵\n";
	for (int cellIndex = 0; cellIndex < ncells; cellIndex++)
	{
		// 计算传播散射矩阵
		Cells(cellIndex).propagateSMAtrix.S11 = cx_mat(nmodes, nmodes);
		Cells(cellIndex).propagateSMAtrix.S12 =
			diagmat(exp(-IU * Cells(cellIndex).neffs * 2 * PI / Cells(cellIndex).lambda));
		Cells(cellIndex).propagateSMAtrix.S21 =
			diagmat(exp(+IU * Cells(cellIndex).neffs * 2 * PI / Cells(cellIndex).lambda));
		Cells(cellIndex).propagateSMAtrix.S22 = cx_mat(nmodes, nmodes);

	}

	// 计算重叠积分
	cout << "\t计算重叠积分\n";
	for (int cellIndex = 0; cellIndex < ncells; cellIndex++)
	{
		Cells(cellIndex).overlapNowPrevious.set_size(nmodes, nmodes);
		Cells(cellIndex).overlapNowNext.set_size(nmodes, nmodes);

		cx_mat Px;
		mat PxReal;
		for (int i = 0;i < nmodes;i++)
		{
			for (int j = 0;j < nmodes;j++)
			{
				if (cellIndex > 0)
				{
					Px = Cells(cellIndex).modeFields(i).ey % Cells(cellIndex - 1).modeFields(j).hz
						- Cells(cellIndex).modeFields(i).ez % Cells(cellIndex - 1).modeFields(j).hy;
					PxReal = trapz(Cells(cellIndex).z, trapz(Cells(cellIndex).y
						, real(Px), 0), 1);
					Cells(cellIndex).overlapNowPrevious(i, j) = PxReal(0);
				}
				if (cellIndex < ncells - 1)
				{
					Px = Cells(cellIndex).modeFields(i).ey % Cells(cellIndex + 1).modeFields(j).hz
						- Cells(cellIndex).modeFields(i).ez % Cells(cellIndex + 1).modeFields(j).hy;
					PxReal = trapz(Cells(cellIndex).z, trapz(Cells(cellIndex).y
						, real(Px), 0), 1);
					Cells(cellIndex).overlapNowNext(i, j) = PxReal(0);
				}
			}
		}

	}

	// 计算链接散射矩阵
	cout << "\t计算链接矩阵\n";









}


int main() {

	loadHdf5();
	return 0;
}