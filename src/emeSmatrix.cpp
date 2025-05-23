﻿
#include "../include/emeSmatrix.h"

void emepropagateSmatrix()
{

	// 读取x
	cout << "\t读取文件数据\n";
	string fileName = "E:/VisualStudioFiles/EMEsolver/lumerical/out.h5";

	int nmodes = 6;
	int ncells = 10;
	int ny = 81;
	int nz = 81;
	double dx = 10e-6 / ncells;
	double lambda = 1.55e-6;


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
		Cells(indexCell).lambda = lambda;

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
	// 计算传播散射矩阵
	cout << "\t计算传播散射矩阵\n";
	for (int cellIndex = 0; cellIndex < ncells; cellIndex++)
	{
		// 计算传播散射矩阵
		Cells(cellIndex).propagateSMAtrix.S11 = cx_mat(nmodes, nmodes);
		Cells(cellIndex).propagateSMAtrix.S12 =
			diagmat(exp(+IU * Cells(cellIndex).neffs * 2 * PI / Cells(cellIndex).lambda * dx / 2));
		Cells(cellIndex).propagateSMAtrix.S21 =
			diagmat(exp(+IU * Cells(cellIndex).neffs * 2 * PI / Cells(cellIndex).lambda * dx / 2));
		Cells(cellIndex).propagateSMAtrix.S22 = cx_mat(nmodes, nmodes);
	}

	// 计算链接散射矩阵
	cout << "\t计算链接矩阵\n";
	field<SMatrix> juncSMatrix(ncells - 1);
	for (int i = 0; i < juncSMatrix.n_elem;i++)
	{
		// i + 1 
		mat O12 = Cells(i).overlapNowNext;
		mat O21 = Cells(i + 1).overlapNowPast;

		mat T12 = 2.0 * inv(O21.st() + O12);
		mat R12 = 0.5 * (O21.st() - O12) * T12; //对的
		mat T21 = 2.0 * inv(O12.st() + O21);
		mat R21 = 0.5 * (O12.st() - O21) * T21;

		juncSMatrix(i).S11 = R12 + 0.0 * IU;
		juncSMatrix(i).S12 = T21 + 0.0 * IU;
		juncSMatrix(i).S21 = T12 + 0.0 * IU;
		juncSMatrix(i).S22 = R21 + 0.0 * IU;
		//cout << i + 1 << endl;
		//juncSMatrix(i).S12.print();
		//cout << -inv(T21)*R12 ;
	}

	// 计算左边界处散射矩阵
	SMatrix sLeft{ cx_mat(nmodes,nmodes) ,
	eye(nmodes, nmodes) + 0.0 * IU  ,
	eye(nmodes, nmodes) + 0.0 * IU ,
	cx_mat(nmodes,nmodes) };

	// 计算右边界处散射矩阵
	SMatrix sRight{ cx_mat(nmodes,nmodes) ,
			eye(nmodes, nmodes) + 0.0 * IU  ,
			eye(nmodes, nmodes) + 0.0 * IU ,
			cx_mat(nmodes,nmodes) };


	// 链接全局矩阵
	//	初始化
	SMatrix globalSMatrix{ cx_mat(nmodes,nmodes) ,
	eye(nmodes, nmodes) + 0.0 * IU  ,
	eye(nmodes, nmodes) + 0.0 * IU ,
	cx_mat(nmodes,nmodes) };

	// 计算模式系数
	cx_mat a(nmodes, ncells + 2);	//正向系数
	// 基模 index=1；前向
	a(0, 0) = 1;

	// 左边的链接矩阵
	SconnectRight(globalSMatrix, sLeft);
	// 传输到第一个Cell中点的模式系数
	SconnectRight(globalSMatrix, Cells(0).propagateSMAtrix);
	a.col(0 + 1) = globalSMatrix.S21 * a.col(0);

	//从左边开始传播
	for (int i = 0; i < a.n_cols - 3; i++)
	{
		//  -|-  计算位置循环
		SconnectRight(globalSMatrix, Cells(i).propagateSMAtrix);
		SconnectRight(globalSMatrix, juncSMatrix(i));
		SconnectRight(globalSMatrix, Cells(i + 1).propagateSMAtrix);
		a.col(i + 2) = globalSMatrix.S21 * a.col(0);
	}
	// 最后模式系数
	SconnectRight(globalSMatrix, Cells(Cells.n_elem - 1).propagateSMAtrix);

	//传输到最后端点
	//右边界的链接矩阵
	SconnectRight(globalSMatrix, sRight);
	a.col(a.n_cols - 1) = globalSMatrix.S21 * a.col(0);

	globalSMatrix.S11.print();

	mat real_a = real(a);
	mat imag_a = imag(a);
	mat abs_a = abs(a);
	//abs_a.print();

	mat power_a = real(a % conj(a));
	//abs_a.st().print();

	//	写入文件
	cout << "\t写入文件\n";
	string filePath = "E:/VisualStudioFiles/EMEsolver/lumerical/emeout.h5";
	remove(filePath.c_str()); //删除文件
	a.save(hdf5_name(filePath, "a", hdf5_opts::append));
	abs_a.save(hdf5_name(filePath, "abs_a", hdf5_opts::append));
	power_a.save(hdf5_name(filePath, "power_a", hdf5_opts::append));
}
