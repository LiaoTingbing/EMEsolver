#pragma once
 

#include<vector>
#include<map>
#include<cstring>
#include <iostream>
#define ARMA_USE_HDF5
#define ARMA_USE_SUPERLU
#include <armadillo>

using namespace arma;
using namespace std;

const double PI = 3.1415926;
const cx_double IU(0, 1);

struct ModeField
{
	//	x方向传播模场
	cx_mat ex;
	cx_mat ey;
	cx_mat ez;
	cx_mat hx;
	cx_mat hy;
	cx_mat hz;
};


struct SMatrix
{
	cx_mat S11;
	cx_mat S12;
	cx_mat S21;
	cx_mat S22;
};

struct Cell
{
	field<ModeField> modeFields;
	double lambda;
	cx_vec neffs;
	vec y;
	vec z;
	vec x;
	SMatrix propagateSMAtrix;
	mat overlapNowPast;
	mat overlapNowNext;

};

struct TMatrix
{
	cx_mat T11;
	cx_mat T12;
	cx_mat T21;
	cx_mat T22;
};