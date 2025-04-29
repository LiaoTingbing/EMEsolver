
#include "include/common.h"


map<string, cube> dev;

void loadHdf5()
{
	// 读取x

	string filename = "E:/VisualStudioFiles/EMEsolver/lumerical/out.h5";
	cube tmp;
	tmp.load(hdf5_name(filename,"z"));
	cout << tmp.size();
	

}


int main() {
	
	loadHdf5();
	return 0;
}