#include "../include/load_hdf5.h"

void loadHdf5Data(map<string, cube>& dev,
	const string filePath,
	const field<string>& dataSetName)
 {
		cout << "加载数据	：\n";
		cube tmp;
		for (int i = 0; i < dataSetName.size(); i++) {
			tmp.load(hdf5_name(filePath, dataSetName(i)));
			cout << "\t" << dataSetName(i) << "\t\t" <<
				tmp.n_rows << " " << tmp.n_cols << " " << tmp.n_slices << endl;
			dev[dataSetName(i)] = tmp;
		}
}

