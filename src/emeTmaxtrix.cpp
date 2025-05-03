
#include "../include/emeTmaxtrix.h"


void emepropagateTmatrix()
{

	// 读取x
	cout << "\t读取文件数据\n";
	string fileName = "E:/VisualStudioFiles/EMEsolver/lumerical/out.h5";

	int nmodes = 6;
	int ncells = 40;
	int ny = 81;
	int nz = 81;
	double dx = 10e-6/ ncells;


	field<Cell> Cells(ncells);
	// 加载有效折射率

	mat realNeffs;
	mat imagNeffs;
	realNeffs.load(hdf5_name(fileName, "real_neffs"));
	imagNeffs.load(hdf5_name(fileName, "imag_neffs"));

	//realNeffs.print();
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


	// 计算重叠积分
	cout << "\t计算重叠积分\n";
	for (int cellIndex = 0; cellIndex < ncells; cellIndex++)
	{
		Cells(cellIndex).overlapNowPast.set_size(nmodes, nmodes);
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
					Cells(cellIndex).overlapNowPast(i, j) = PxReal(0);
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
		//cout << cellIndex + 1 << endl;

		//cout << "now next" << endl;
		//Cells(cellIndex).overlapNowNext.print();
		//cout << "now pre" << endl;

		//Cells(cellIndex).overlapNowPrevious.print();

	}
	// 计算传播矩阵
	cout << "\t计算传播矩阵\n";
	field<cx_mat> propTransfer(ncells);
	for (int cellIndex = 0; cellIndex < ncells; cellIndex++)
	{
		// 传播传输矩阵
		propTransfer(cellIndex).set_size(2 * nmodes, 2 * nmodes);

		propTransfer(cellIndex)(nmodes, 0, size(nmodes, nmodes))
			= cx_mat(nmodes, nmodes, fill::zeros);
		propTransfer(cellIndex)(0, nmodes, size(nmodes, nmodes))
			= cx_mat(nmodes, nmodes, fill::zeros);

		propTransfer(cellIndex)(0,0,size(nmodes,nmodes))=
			diagmat(exp(-IU * Cells(cellIndex).neffs * 2 * PI / Cells(cellIndex).lambda * dx / 2));
		propTransfer(cellIndex)(nmodes, nmodes, size(nmodes, nmodes)) =
			diagmat(exp(+IU * Cells(cellIndex).neffs * 2 * PI / Cells(cellIndex).lambda * dx / 2));
	}

	// 计算链接传输矩阵
	cout << "\t计算链接矩阵\n";
	//field<SMatrix> juncSMatrix(ncells - 1);
	field<cx_mat> juncTransfer(ncells - 1);
	for (int i = 0; i < juncTransfer.n_elem;i++)
	{
		// i + 1 
		mat O12 = Cells(i).overlapNowNext;
		mat O21 = Cells(i + 1).overlapNowPast;
		// 传输矩阵
		mat T11 = (O12.st() + O21) / 2;
		mat T12 = (O12.st() - O21) / 2;
		mat T21 = (O12.st() - O21) / 2;
		mat T22 = (O12.st() + O21) / 2;

		//	传输矩阵
		juncTransfer(i).set_size(2 * nmodes, 2 * nmodes);
		juncTransfer(i)(0, 0, size(nmodes, nmodes)) = T11+0.0*IU;
		juncTransfer(i)(0, nmodes, size(nmodes, nmodes)) = T12 + 0.0 * IU;
		juncTransfer(i)(nmodes, 0, size(nmodes, nmodes)) = T21 + 0.0 * IU;
		juncTransfer(i)(nmodes, nmodes, size(nmodes, nmodes)) = T22 + 0.0 * IU;
	}

	// 计算模式系数
		// 
	cx_mat aT(2 * nmodes, ncells + 2);
	aT(0) = 1.0 ;    //基模
	aT.col(1) = propTransfer(0) * aT.col(0); //第一个Cell中间值
	aT.col(2) = propTransfer(1) * juncTransfer(0) * propTransfer(0) * aT.col(1);
	for (int i = 0; i < ncells - 1; i++)
	{
		aT.col(i + 2) = propTransfer(i + 1) * juncTransfer(i) * propTransfer(i) * aT.col(i + 1);
	}
	aT.col(ncells + 1) = propTransfer(ncells - 1) * aT.col(ncells);//最后一个Cell中间值

	mat real_aT = real(aT);
	mat power_aT = real(aT % conj(aT));

	real_aT.print("real_aT");
	power_aT.print("power_aT");

 

	//	写入文件
	cout << "\t写入文件\n";
	string filePath = "E:/VisualStudioFiles/EMEsolver/lumerical/emeout.h5";
	remove(filePath.c_str());
	power_aT.save(hdf5_name(filePath, "power_aT", hdf5_opts::append));

 



}
