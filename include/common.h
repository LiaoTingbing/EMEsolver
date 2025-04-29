#pragma once
 

#include<map>
#include<cstring>
#include <iostream>
#define ARMA_USE_HDF5
#define ARMA_USE_SUPERLU
#include <armadillo>

using namespace arma;
using namespace std;


struct ModeField
{
	//	x方向传播模场
	cx_mat ex;
	cx_mat ey;
	cx_mat ez;
	cx_mat hx;
	cx_mat hy;
	cx_mat hz;
	cx_double neff;
	double lambda;
	vec y;
	vec z;
	vec x;
};

struct Cell
{
	field<ModeField> modeFields;
};
