/*
MIT License

Copyright (c) 2022 Siyou Lin, Dong Xiao, Zuoqiang Shi, Bin Wang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "Octree.h"
#include <iostream>
#include <algorithm>
#include <cstring>
#include <string>
#include <fstream>
#include "Geometry.h"

#include <CLI11.hpp>
#include <cnpy.h>

int main(int argc, char** argv) {

    std::string ply_suffix(".ply");
	std::string normalized_npy_suffix = "_normalized.npy";
	std::string query_npy_suffix("_for_query.npy");
    
	std::string inFileName;
	std::string outFileName;
	int minDepth = 1;
	int maxDepth = 10;
    
    CLI::App app("PGRExportQuery");
    app.add_option("-i", inFileName, "input filename of xyz format")->required();
	app.add_option("-o", outFileName, "output filename with no suffix")->required();
	app.add_option("-m", minDepth, "");
	app.add_option("-d", maxDepth, "");
	
    CLI11_PARSE(app, argc, argv);

	if (maxDepth < minDepth) {
		cout << "[In PGRExportQuery] WARNING: minDepth "
			 << minDepth
			 << " smaller than maxDepth "
			 << maxDepth
			 << ", ignoring given minDepth\n";
	}
		
	Octree tree;
	tree.setTree(inFileName, maxDepth, minDepth);//1382_seahorse2_p

	//*** Nodes for query are from gridDataVector *** START ***
	unsigned long N_grid_pts = tree.gridDataVector.size();
	unsigned long N_sample_pts = tree.samplePoints.size();
	std::vector<float> grid_coords;
	
	for(int grid_idx=0; grid_idx<N_grid_pts; grid_idx++) {
		grid_coords.push_back( tree.gridDataVector[grid_idx]->coords[0] );
		grid_coords.push_back( tree.gridDataVector[grid_idx]->coords[1] );
		grid_coords.push_back( tree.gridDataVector[grid_idx]->coords[2] );
	}

	// exporting normalized point samples as npy
	std::vector<float> pts_normalized;
	for (int i=0; i<N_sample_pts; i++){
		pts_normalized.push_back( tree.samplePoints[i].x );
		pts_normalized.push_back( tree.samplePoints[i].y );
		pts_normalized.push_back( tree.samplePoints[i].z );
	}

	cnpy::npy_save(outFileName + normalized_npy_suffix, &pts_normalized[0], {N_sample_pts, 3}, "w");
	std::cout << "[In PGRExportQuery] Normalizing the point cloud. Result saved to " << outFileName + normalized_npy_suffix <<std::endl;
	cnpy::npy_save(outFileName + query_npy_suffix, &grid_coords[0], {N_grid_pts, 3}, "w");
	std::cout << "[In PGRExportQuery] Exporting points on octree for query. Result saved to " << outFileName + query_npy_suffix <<std::endl;
}
