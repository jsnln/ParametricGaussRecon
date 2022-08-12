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
#include <cstdlib>
#include <cstdio>
#include "Geometry.h"

#include "CLI11.hpp"

int main(int argc, char** argv) {

    std::string ply_suffix(".ply");
    std::string inFileName;
    std::string outFileName;
	std::string inGridValFileName;
	std::string inGridWidthFileName;
	float isovalue = 0.5;
	int minDepth = 1;
	int maxDepth = 10;
    
    CLI::App app("PGRLoadQuery");
    app.add_option("-i", inFileName, "input filename of xyz format")->required();
    app.add_option("-o", outFileName, "output filename of ply format")->required();
	app.add_option("--grid_val", inGridValFileName, "input filename of npy format for grid vals")->required();
	app.add_option("--grid_width", inGridWidthFileName, "input filename of npy format for grid widths")->required();
	app.add_option("-m", minDepth, "");
	app.add_option("-d", maxDepth, "");
	app.add_option("--isov", isovalue, "isovalue");

    CLI11_PARSE(app, argc, argv);

	if (maxDepth < minDepth) {
		cout << "[In PGRLoadQuery] WARNING: minDepth "
			 << minDepth
			 << " smaller than maxDepth "
			 << maxDepth
			 << ", ignoring given minDepth\n";
	}
	
	Octree tree;
	tree.setTree(inFileName, maxDepth, minDepth);//1382_seahorse2_p
	
	int N_grid = tree.gridDataVector.size();
	tree.loadImplicitFunctionFromNPY(inGridValFileName, N_grid);
	tree.loadGridWidthFromNPY(inGridWidthFileName, N_grid);
	
	std::cout << "[In PGRLoadQuery] Isovalue: " << isovalue << std::endl;

	CoredVectorMeshData mesh;
	tree.GetMCIsoTriangles(isovalue,  &mesh, 0, 1, 0, 0);
	char fileChar[255];
	strcpy(fileChar, (outFileName).c_str());
	tree.writePolygon2(&mesh, fileChar);
	std::cout << "[In PGRLoadQuery] Polygon Written to " << outFileName << std::endl;
}
