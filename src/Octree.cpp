/*
Copyright (c) 2018, Wenjia Lu, Zuoqiang Shi, Jian Sun and Bin Wang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of Tsinghua University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include "Octree.h"

#include <fstream>
#include "MarchingCubes.h"
#include "Geometry.h"
#include "ply.h"
#include <cnpy.h>

vector<WeightedPoint> Octree::samplePoints;

// constructor
Octree::Octree() {
	maxDepth = 0;
	bb.blx = 0; bb.bly = 0; bb.blz = 0;
	bb.xscale = 0; bb.yscale = 0; bb.zscale = 0;
	counterGrid = 0;
	counterIso = 0;
	scaleFactor = 1.1;
	samplePoints.clear();
	gridIdx = 0;
	//myTotalArea = 0;
}

bool Octree::setTree(const string filename, int d, int d_min) {
	if (d > DEPTH_LIMIT) {
		cout << "[In PGROctree] Default depth limit " << DEPTH_LIMIT << " exceeded, resetting d to " << DEPTH_LIMIT << " instead.\n";
		d = DEPTH_LIMIT;
	}
	maxDepth = d;
	if (!readFile(filename, d_min))
		return false;
	return true;
}

bool Octree::readFile(const string filename, int min_depth){
	//double start = Time();
	if (maxDepth <= 0){
		printf("[In PGROctree] Max Depth must be a positive number!\n");
		return false;
	}
	// get total number of samples and bounding box
	FILE* fp;
	fp = fopen(filename.data(), "r");
	if (!fp){
		printf("[In PGROctree] Cannot open file %s ... \n", filename);
		return false;
	}
	int minDepth = min_depth;
	unsigned long lineNo = 0;			// record number of lines, i.e., number of samples
	float x, y, z;
	float minX, maxX, minY, maxY, minZ, maxZ;
	samplePoints.clear();
	while (1){
		if (fscanf(fp, "%f %f %f", &x, &y, &z) != 3)
			break;
		
		if (!lineNo || x > maxX) maxX = x;
		if (!lineNo || x < minX) minX = x;
		if (!lineNo || y > maxY) maxY = y;
		if (!lineNo || y < minY) minY = y;
		if (!lineNo || z > maxZ) maxZ = z;
		if (!lineNo || z < minZ) minZ = z;

		samplePoints.push_back(WeightedPoint(x, y, z));
		lineNo++;

	}
	fclose(fp);
	//printf("%d samples input...\n", lineNo);
	if (lineNo == 0)
		return false;
	// set bounding box
	bb.blx = minX;
	bb.bly = minY;
	bb.blz = minZ;
	bb.xscale = maxX - minX;
	bb.yscale = maxY - minY;
	bb.zscale = maxZ - minZ;

	neighborKey.set(maxDepth);
	float center[3];
	float myCenter[3];
	OctNode* temp;
	float myWidth;
	maxScale = max(max(bb.xscale, bb.yscale), bb.zscale) * scaleFactor;
	
	//01-09
	center[0] = (maxX + minX - maxScale) / 2;
	center[1] = (maxY + minY - maxScale) / 2;
	center[2] = (maxZ + minZ - maxScale) / 2;
	
	// splat w.r.t. maxdepth-2
	int splatDepth = maxDepth - 2;
	if (splatDepth > 0){
		//cout << "Setting sample weights..." << endl;
		for (int i = 0; i < lineNo; i++){

			//01-09
			// normalizing sample points
			samplePoints[i].x = (samplePoints[i].x - center[0]) / maxScale;
			samplePoints[i].y = (samplePoints[i].y - center[1]) / maxScale;
			samplePoints[i].z = (samplePoints[i].z - center[2]) / maxScale;
	
			// setting tree nodes
			myCenter[0] = 0.5;
			myCenter[1] = 0.5;
			myCenter[2] = 0.5;
			myWidth = 1.0;

			temp = &root;
			int d = 0;
			float weight = 1.0;
			float positionArray[3] = { samplePoints[i].x, samplePoints[i].y, samplePoints[i].z };
			while (d < splatDepth){
				NonLinearUpdateWeightContribution(temp, positionArray, weight);
				if (!temp->children) temp->initChildren();
				int cIndex = Cube::CornerIndex(myCenter, positionArray);
				temp = &temp->children[cIndex];
				myWidth /= 2;
				if (cIndex & 1) myCenter[0] += myWidth / 2;
				else			 myCenter[0] -= myWidth / 2;
				if (cIndex & 2) myCenter[1] += myWidth / 2;
				else			 myCenter[1] -= myWidth / 2;
				if (cIndex & 4)  myCenter[2] += myWidth / 2;
				else			 myCenter[2] -= myWidth / 2;
				d++;
			}

			NonLinearUpdateWeightContribution(temp, positionArray, weight);
		}
	}
	//cout << "First: Leaves/Nodes: " << root.leaves() << "/" << root.nodes() << endl;

	// setting normals (smoothing?)
	normalLengths.clear();

	float normal_len = float(2 << maxDepth);

	for (int i = 0; i < lineNo; i++){
		NonLinearSplatOrientedPoint(samplePoints[i], normal_len, splatDepth, 1, minDepth, maxDepth, i);
	}
	//cout << "Second: Leaves/Nodes: " << root.leaves() << "/" << root.nodes() << endl;


	ClipTree();
	
	//cout << "\tSet And Clip: " << Time() - start << endl;

	//start = Time();
	maxDepth = root.maxDepth() + 1;

	minStep = 1.0 / pow(2, maxDepth - 1) * min(bb.xscale, min(bb.yscale, bb.zscale)) / max(bb.xscale, max(bb.yscale, bb.zscale));

	cout << "[In PGROctree] Max Depth: " << maxDepth - 1 << endl;

	OctNode* node = root.nextNode();
	while (node){
		int nodeSize = node->nodeIdx.size();
		if (nodeSize != 0){
			if (node->children){
				for (int i = 0; i < nodeSize; i++){
					node->centerAndWidth(center, myWidth);
					int cIndex = Cube::CornerIndex(center, &samplePoints[node->nodeIdx[i]]);
					node->children[cIndex].nodeIdx.push_back(node->nodeIdx[i]);
				}
				node->nodeIdx.clear();
			}
		}
		node = root.nextNode(node);
	}

	node = root.nextLeaf();
	while (node){
		if (node->nodeIdx.size() > 0){
			node->hasSample = true;
			temp = node;
			while (temp->parent && !temp->parent->hasSample){
				temp->parent->hasSample = true;
				temp = temp->parent;
			}
		}
		node = root.nextLeaf(node);
	}

	//cout << "\tFall Grid: " << Time() - start << endl;
	
	//start = Time();
	gridValueMap.clear();
	node = root.nextLeaf();
	float wholeGrid = 1.0f / pow(2, maxDepth + 1);

	while (node){
		for (int i = 0; i < Cube::CORNERS; i++){
			long long key = node->getCornerIndex(i, maxDepth);
			unordered_map<long long, GridData>::iterator iter = gridValueMap.find(key);
			//GridData* gd = &(gridValueMap[key]);
			if (iter == gridValueMap.end()){
				//cout << gd << endl;
				GridData gd;
				gd.value = 0;
				gd.maxDepth = node->Depth();
				gd.minDepth = node->Depth();
				gd.node = node;
				gd.key = key;
				gd.coords[0] = (key & ((1 << (DEPTH_LIMIT + 1)) - 1)) * (wholeGrid);
				gd.coords[1] = ((key >> DEPTH_LIMIT) & ((1 << (DEPTH_LIMIT + 1)) - 1)) * (wholeGrid);
				gd.coords[2] = ((key >> (DEPTH_LIMIT * 2))) * (wholeGrid);
				//gridValueMap[key] = gd;
				pair<unordered_map<long long, GridData>::iterator, bool> insertIter = gridValueMap.insert(pair<long long, GridData>(key, gd));
				node->cornerGrid[i] = &(insertIter.first->second);
			}
			else {
				node->cornerGrid[i] = &(iter->second);
			}

		}
		node = root.nextLeaf(node);
	}

	gridDataVector.resize(gridValueMap.size());
	//fstream out("result/normalVar.txt", ios::out);
	setGridNode(&root);
	for (unordered_map<long long, GridData>::iterator iter = gridValueMap.begin(); iter != gridValueMap.end(); iter++){

		for (int i = 0; i < Cube::NEIGHBORS; i++){
			if (iter->second.adjacent[i] == NULL && iter->second.adjacentKey[i] != -1){
				//iter->second.adjacent[i] = &(gridValueMap[iter->second.adjacentKey[i]]);
				iter->second.adjacent[i] = &(gridValueMap.find(iter->second.adjacentKey[i])->second);
				if (i % 2 == 0){
					iter->second.adjacent[i]->adjacent[i + 1] = &(iter->second);
				}
				else {
					iter->second.adjacent[i]->adjacent[i - 1] = &(iter->second);
				}
			}
		}
	}

	return true;
}

void Octree::setGridNode(OctNode* currentNode){
	if (currentNode->children == NULL){
		currentNode->gridStartIdx = gridIdx;
		for (int i = 0; i < Cube::CORNERS; i++){
			GridData* gd = currentNode->cornerGrid[i];
			if (gd == NULL){
				long long key = currentNode->getCornerIndex(i, maxDepth);
				gd = &(gridValueMap[key]);
				currentNode->cornerGrid[i] = gd;
			}
			if (gd->maxDepth < currentNode->Depth())
				gd->maxDepth = currentNode->Depth();
			if (gd->minDepth > currentNode->Depth())
				gd->minDepth = currentNode->Depth();
			for (int m = 0; m < Cube::NEIGHBORS; m++){
				if (Cube::CornerAdjacentMap[i][m] == -1)
					continue;
				if (gd->adjacentKey[m] == NULL)
					gd->adjacentKey[m] = currentNode->getCornerIndex(Cube::CornerAdjacentMap[i][m], maxDepth);
				else {
					long long adjacentKey = currentNode->getCornerIndex(Cube::CornerAdjacentMap[i][m], maxDepth);
					if (abs(gd->key - gd->adjacentKey[m]) > abs(gd->key - adjacentKey))
						gd->adjacentKey[m] = adjacentKey;
				}

			}
			if (gd->vectorSet)
				continue;
			gridDataVector[gridIdx++] = gd;
			gd->vectorSet = true;
		}
		currentNode->gridEndIdx = gridIdx;
		int gridCount = currentNode->gridEndIdx - currentNode->gridStartIdx;
		//fstream haha("result1/loc.txt", ios::app);
		for (int j = currentNode->gridStartIdx; j < currentNode->gridEndIdx; j++){
			currentNode->barycenter.x += gridDataVector[j]->coords[0] / gridCount;
			currentNode->barycenter.y += gridDataVector[j]->coords[1] / gridCount;
			currentNode->barycenter.z += gridDataVector[j]->coords[2] / gridCount;
		}
		//haha << currentNode->barycenter.x << " " << currentNode->barycenter.y << " " << currentNode->barycenter.z << endl;
		return;
	}
	for (int i = 0; i < Cube::CORNERS; i++){
		setGridNode(&currentNode->children[i]);
		currentNode->barycenter.x += currentNode->children[i].barycenter.x * (currentNode->children[i].gridEndIdx - currentNode->children[i].gridStartIdx);
		currentNode->barycenter.y += currentNode->children[i].barycenter.y * (currentNode->children[i].gridEndIdx - currentNode->children[i].gridStartIdx);
		currentNode->barycenter.z += currentNode->children[i].barycenter.z * (currentNode->children[i].gridEndIdx - currentNode->children[i].gridStartIdx);
	}
	currentNode->gridStartIdx = currentNode->children[0].gridStartIdx;
	currentNode->gridEndIdx = currentNode->children[Cube::CORNERS - 1].gridEndIdx;
	int gridCount2 = currentNode->gridEndIdx - currentNode->gridStartIdx;
	if (gridCount2 == 0) return;
	currentNode->barycenter.x /= gridCount2;
	currentNode->barycenter.y /= gridCount2;
	currentNode->barycenter.z /= gridCount2;
}

void Octree::initLeaf(){
	leafVector.clear();
	OctNode* leaf = root.nextLeaf();
	while (leaf){
		if (leaf->hasSample)
			leafVector.push_back(leaf);
		leaf = root.nextLeaf(leaf);
	}
}

void Octree::loadImplicitFunctionFromNPY(std::string npyFileName, int N_grid) {
	cnpy::NpyArray grid_vals_from_file = cnpy::npy_load(npyFileName);
    float* loaded_array = grid_vals_from_file.data<float>();
	
	for(int idx=0; idx<N_grid; idx++) {
		this->gridDataVector[idx]->value = *(loaded_array + idx);
	}
}

void Octree::loadGridWidthFromNPY(std::string npyFileName, int N_grid) {
	cnpy::NpyArray grid_widths_from_file = cnpy::npy_load(npyFileName);
    float* loaded_array = grid_widths_from_file.data<float>();
	
	for(int idx=0; idx<N_grid; idx++) {
		this->gridDataVector[idx]->width = *(loaded_array + idx);
	}
}


void Octree::SetIsoSurfaceCorners(const float& isovalue, const int& subdivisionDepth,
	const int& fullDepthIso){
	unordered_map< long long, float> values;
	float cornerValues[Cube::CORNERS];
	OctNode* temp;
	int leafCount = root.leaves();
	long long key;
	// sort tree nodes by depth and shift
	SortedNodes *sNodes = new SortedNodes();
	sNodes->set(root, 0);

	temp = root.nextNode();
	while (temp){
		temp->mcIdx = 0;
		temp = root.nextNode(temp);
	}

	// deal with nodes of subdivideDepth
	for (int i = 0; i < sNodes->nodeCount[subdivisionDepth]; i++){
		temp = sNodes->treeNodes[i];
		for (int c = 0; c < Cube::CORNERS; c++){
			long long cornerKey = temp->getCornerIndex(c, maxDepth);
			unordered_map<long long, GridData>::iterator iter = gridValueMap.find(cornerKey);
			GridData* gd = temp->cornerGrid[c];
			//if (&(iter->second) != gd){
			//	cout << "Error Noe Equal..." << endl;
			//}
			if (gd != NULL)
				cornerValues[c] = gd->value;
			else
				cerr << "[In PGROctree] Cannot find the specified value..." << endl;
			/*if (iter != gridValueMap.end())
				cornerValues[c] = iter->second.value;
				else
				cerr << "cannot find the specified value..." << endl;*/
		}
		// calculate 8 vertices MC value in each grid
		temp->mcIdx = MarchingCubes::GetIndex(cornerValues, isovalue);

		if (temp->parent){
			OctNode* parent = temp->parent;
			int c = int(temp - temp->parent->children);	// tell which child node it is
			// check if at parent node < isovalue
			int mcIdx = temp->mcIdx & (1 << MarchingCubes::cornerMap[c]);
			// if at the parent node < isovalue, pass the mcIndex to parent node and the parent node of the parent node
			if (mcIdx){
				parent->mcIdx |= mcIdx;
				while (1){
					if (parent->parent && (parent - parent->parent->children) == c){
						parent->parent->mcIdx |= mcIdx;
						parent = parent->parent;
					}
					else break;
				}
			}
		}
	}


	// deal with leaf nodes
	for (int i = sNodes->nodeCount[subdivisionDepth]; i < sNodes->nodeCount[subdivisionDepth + 1]; i++){
		temp = sNodes->treeNodes[i]->nextLeaf();
		while (temp){
			for (int c = 0; c < Cube::CORNERS; c++){
				GridData* gd = temp->cornerGrid[c];
				if (gd != NULL)
					cornerValues[c] = gd->value;
				else
					cerr << "[In PGROctree] cannot find the specified value..." << endl;		
			}


			temp->mcIdx = MarchingCubes::GetIndex(cornerValues, isovalue);
			if (temp->parent){
				OctNode* parent = temp->parent;
				int c = int(temp - temp->parent->children);
				int mcIdx = temp->mcIdx & (1 << MarchingCubes::cornerMap[c]);

				if (mcIdx){
					parent->mcIdx |= mcIdx;
					while (1){
						if (parent->parent && (parent - parent->parent->children) == c){
							parent->parent->mcIdx |= mcIdx;
							parent = parent->parent;
						}
						else break;
					}
				}
			}
			temp = sNodes->treeNodes[i]->nextLeaf(temp);
		}
	}
	delete sNodes;

}

int Octree::GetRootIndex(const OctNode* node, const int& edgeIndex, const int& maxDepth, RootInfo& ri){
	const OctNode *temp, *finest;
	int finestIndex;

	// ��node��edgeIndex���н���
	if (!(MarchingCubes::edgeMask[node->mcIdx] & (1 << edgeIndex)))
		return 0;

	int f1, f2;
	// ���edgeIndex���ڵ�������f1��f2
	Cube::FacesAdjacentToEdge(edgeIndex, f1, f2);

	finest = node;
	finestIndex = edgeIndex;
	if (node->Depth() < maxDepth){
		// �ҵ���node��f1���ڵĽڵ㣬�������ӽڵ㣨��ϸ�Ľڵ㣩
		temp = node->faceNeighbor(f1);
		if (temp && temp->children){
			finest = temp;
			finestIndex = Cube::FaceReflectEdgeIndex(edgeIndex, f1);
		}
		else {
			temp = node->faceNeighbor(f2);
			if (temp && temp->children){
				finest = temp;
				finestIndex = Cube::FaceReflectEdgeIndex(edgeIndex, f2);
			}
			else{
				temp = node->edgeNeighbor(edgeIndex);
				if (temp && temp->children){
					finest = temp;
					finestIndex = Cube::EdgeReflectEdgeIndex(edgeIndex);
				}
			}
		}
	}

	int c1, c2;
	Cube::EdgeCorners(finestIndex, c1, c2);
	// �����ϸһ��Ľڵ㻹û�е���ϸ����һ�㣬���еݹ�
	if (finest->children){
		if (GetRootIndex(&finest->children[c1], finestIndex, maxDepth, ri)) return 1;
		else if (GetRootIndex(&finest->children[c2], finestIndex, maxDepth, ri)) return 1;
		else return 0;
	}
	else {
		int o, i1, i2;
		Cube::FactorEdgeIndex(finestIndex, o, i1, i2);
		int d, off[3];
		finest->depthAndOffset(d, off);
		ri.node = finest;
		ri.edgeIndex = finestIndex;
		int offset, eIndex[2];
		offset = (1 << d) + off[o] - 1;
		switch (o){
		case 0:
			eIndex[0] = (off[1] + i1) << (maxDepth - d + 1);
			eIndex[1] = (off[2] + i2) << (maxDepth - d + 1);
			break;
		case 1:
			eIndex[0] = (off[0] + i1) << (maxDepth - d + 1);
			eIndex[1] = (off[2] + i2) << (maxDepth - d + 1);
			break;
		case 2:
			eIndex[0] = (off[0] + i1) << (maxDepth - d + 1);
			eIndex[1] = (off[1] + i2) << (maxDepth - d + 1);
			break;
		}
		ri.key = (long long)(o) | (long long)(eIndex[0]) << 5 | (long long)(eIndex[1]) << 25 | (long long)(offset) << 45;
		return 1;
	}
}

int Octree::GetRootIndex(const long long& key, unordered_map<long long, int>& boundaryRoots, unordered_map<long long, int>* interiorRoots, CoredPointIndex& index){
	unordered_map<long long, int>::iterator rootIter = boundaryRoots.find(key);
	// �����boundaryRoots���ҵ�key��Ӧ��root
	if (rootIter != boundaryRoots.end()){
		index.inCore = 1;
		index.index = rootIter->second;
		return 1;
	}
	else if (interiorRoots){
		rootIter = interiorRoots->find(key);
		if (rootIter != interiorRoots->end()){
			index.inCore = 0;
			index.index = rootIter->second;
			return 1;
		}
	}
	return 0;
}

int Octree::IsBoundaryEdge(const OctNode* node, const int& edgeIndex, const int& subdivideDepth){
	int dir, x, y;
	Cube::FactorEdgeIndex(edgeIndex, dir, x, y);
	return IsBoundaryEdge(node, dir, x, y, subdivideDepth);
}

int Octree::IsBoundaryEdge(const OctNode* node, const int& dir, const int& x, const int& y, const int& subdivideDepth){
	int d, o[3], idx1, idx2, mask;

	if (subdivideDepth < 0){ return 0; }
	if (node->Depth() <= subdivideDepth){ return 1; }
	node->depthAndOffset(d, o);

	switch (dir){
	case 0:
		idx1 = (int(o[1]) << 1) + (x << 1);
		idx2 = (int(o[2]) << 1) + (y << 1);
		break;
	case 1:
		idx1 = (int(o[0]) << 1) + (x << 1);
		idx2 = (int(o[2]) << 1) + (y << 1);
		break;
	case 2:
		idx1 = (int(o[0]) << 1) + (x << 1);
		idx2 = (int(o[1]) << 1) + (y << 1);
		break;
	}
	mask = 2 << (int(node->Depth()) - subdivideDepth);
	return !(idx1 % (mask)) || !(idx2 % (mask));
}

int Octree::IsBoundaryFace(const OctNode* node, const int& faceIndex, const int& subdivideDepth){
	int dir, offset, d, o[3], idx;

	if (subdivideDepth < 0){ return 0; }
	if (node->Depth() <= subdivideDepth){ return 1; }
	Cube::FactorFaceIndex(faceIndex, dir, offset);
	node->depthAndOffset(d, o);

	idx = (int(o[dir]) << 1) + (offset << 1);
	return !(idx % (2 << (int(node->Depth()) - subdivideDepth)));
}


int Octree::GetRoot(const RootInfo& ri, const float& isoValue, Point & position, unordered_map<long long, std::pair<float, Point> >& normalHash, const int& nonLinearFit){
	int c1, c2;	// edgeIndex��Ӧ������corner
	Cube::EdgeCorners(ri.edgeIndex, c1, c2);
	if (!MarchingCubes::HasEdgeRoots(ri.node->mcIdx, ri.edgeIndex)) return 0;

	float value1, value2;
	float width1, width2;
	GridData* gd1 = ri.node->cornerGrid[c1];
	if (gd1 == NULL)
		cerr << "[In PGROctree] Use Undefined Value: " << gd1->key;
	else
		//value1 = gridValueMap[key1].value;
		value1 = gd1->value;
		width1 = gd1->width;
	//if (gridValueMap.find(key2) == gridValueMap.end())
	GridData* gd2 = ri.node->cornerGrid[c2];
	if (gd2 == NULL)
		cerr << "[In PGROctree] Use Undefined Value: " << gd2->key;
	else
		//value2 = gridValueMap[key2].value;
		value2 = gd2->value;
		width2 = gd2->width;

	Point coords1, coords2;

	//coords1.x = gridValueMap[key1].coords[0];	coords1.y = gridValueMap[key1].coords[1];	coords1.z = gridValueMap[key1].coords[2];
	//coords2.x = gridValueMap[key2].coords[0];	coords2.y = gridValueMap[key2].coords[1];	coords2.z = gridValueMap[key2].coords[2];
	coords1.x = gd1->coords[0];
	coords1.y = gd1->coords[1];
	coords1.z = gd1->coords[2];
	
	coords2.x = gd2->coords[0];
	coords2.y = gd2->coords[1];
	coords2.z = gd2->coords[2];

	float isov = -isoValue;
	value1 *= -1;
	value2 *= -1;
	float ratio = (value1 - isov) * width1 / ((value1 - isov) * width1 - (value2 - isov) * width2);
	// float ratio = (value1 - isov) / (value1 - value2);
	//ratio = ratio * ratio;
	//cout << "(iso, v1, v2, w1, w2) = (" << isov << "," << value1 << ", " << value2 << ", " << width1 << "," << width2 << "), ratio = " << ratio << "\n";
	if (ratio > 1.0)
		ratio = 1.0;
	if (ratio < 0.0)
		ratio = 0.0;
	
	// if (ratio > 1.0 || ratio < 0)
	// 	cout << "[In PGROctree] Ratio error!" << endl;
	position.x = coords1.x - ratio * (coords1.x - coords2.x);
	position.y = coords1.y - ratio * (coords1.y - coords2.y);
	position.z = coords1.z - ratio * (coords1.z - coords2.z);

	return 1;
}

int Octree::GetRoot(const RootInfo& ri, const float& isoValue, const int& maxDepth, Point & position, unordered_map<long long, std::pair<float, Point > >& normals,
	Point* normal, const int& nonLinearFit){
	if (!MarchingCubes::HasRoots(ri.node->mcIdx)) return 0;
	return GetRoot(ri, isoValue, position, normals, nonLinearFit);
}



int Octree::SetMCRootPositions(OctNode* node, const int& sDepth, const float& isoValue,
	unordered_map<long long, int>& boundaryRoots, unordered_map<long long, int>* interiorRoots,
	unordered_map<long long, std::pair<float, Point > >& boundaryNormalHash, unordered_map<long long, std::pair<float, Point> >* interiorNormalHash,
	std::vector<Point >* interiorPositions,
	CoredVectorMeshData* mesh, const int& nonLinearFit){
	Point position;
	RootInfo ri;
	int count = 0;
	// �ж��Ƿ�Ҫ�Ե�ǰnodeִ��mc����
	if (!MarchingCubes::HasRoots(node->mcIdx)) return 0;

	// ����3�������ϵ�ÿ���߽��д���
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 2; k++){
				long long key;
				int eIndex = Cube::EdgeIndex(i, j, k);
				if (GetRootIndex(node, eIndex, maxDepth, ri)){
					key = ri.key;
					if (!interiorRoots || IsBoundaryEdge(node, i, j, k, sDepth)){
						if (boundaryRoots.find(key) == boundaryRoots.end()){
							GetRoot(ri, isoValue, maxDepth, position, boundaryNormalHash, NULL, nonLinearFit);
							mesh->inCorePoints.push_back(position);
							boundaryRoots[key] = int(mesh->inCorePoints.size()) - 1;
							count++;
						}
					}
					else {

					}
				}
			}
		}
	}
	return count;
}


void Octree::GetMCIsoTriangles(const float& isovalue, CoredVectorMeshData* mesh, const int& fullDepthIso,
	const int& nonLinearFit, bool addBarycenter, bool polygonMesh){

	OctNode* temp;

	unordered_map<long long, int> roots;	// roots: idx of incorePoints in mesh
	unordered_map<long long, pair<float, Point> > *normalHash = new unordered_map< long long, pair<float, Point>>();

	// set 8 values for the 8 vertices of a cube (checked)
	SetIsoSurfaceCorners(isovalue, 0, fullDepthIso);

	int leaves = root.leaves();

	// deal with all edges in leaf nodes, if one of them crosses the surface, add to roots
	temp = root.nextLeaf();
	while (temp){
		SetMCRootPositions(temp, 0, isovalue, roots, NULL, *normalHash, NULL, NULL, mesh, nonLinearFit);
		temp = root.nextLeaf(temp);
	}

	temp = root.nextLeaf();
	while (temp){
		GetMCIsoTriangles(temp, mesh, roots, NULL, NULL, 0, 0, addBarycenter, polygonMesh);
		temp = root.nextLeaf(temp);
	}
}

int Octree::GetRootPair(const RootInfo& ri, const int& maxDepth, RootInfo& pair){
	const OctNode* node = ri.node;
	int c1, c2, c;
	// ���ڵ����ڱ߶�Ӧ����������c1��c2
	Cube::EdgeCorners(ri.edgeIndex, c1, c2);
	while (node->parent){
		c = int(node - node->parent->children);
		if (c != c1 && c != c2) return 0;
		// ���ڵ������������û�и����и���������Ҫ�ر�Ĵ�����
		if (!MarchingCubes::HasEdgeRoots(node->parent->mcIdx, ri.edgeIndex)){
			// �����ǰ�ڵ���c1���غϣ�����c1��edgeIdx��root�Ѿ�֪���ˣ�������Ҫ����c2���д���
			if (c == c1) return GetRootIndex(&node->parent->children[c2], ri.edgeIndex, maxDepth, pair);
			else return GetRootIndex(&node->parent->children[c1], ri.edgeIndex, maxDepth, pair);
		}
		node = node->parent;
	}
	return 0;
}

void Octree::GetMCIsoEdges(OctNode* node, unordered_map<long long, int>& boundaryRoots, unordered_map<long long, int>* interiorRoots, const int& sDepth, std::vector<std::pair<long long, long long> >& edges){
	OctNode* temp;
	int count = 0, tris = 0;
	int isoTri[3 * MarchingCubes::MAX_TRIANGLES];//isoTri��3��5
	FaceEdgesFunction fef;
	unordered_map<long long, std::pair<RootInfo, int> >::iterator iter;
	unordered_map<long long, std::pair<RootInfo, int> > vertexCount;

	fef.edges = &edges;
	fef.maxDepth = maxDepth;
	fef.vertexCount = &vertexCount;

	// ����edgeMask����ߵ����ӹ�ϵ��count�����м���������
	count = MarchingCubes::AddTriangleIndices(node->mcIdx, isoTri);
	// ����ÿһ��cube���ڵ�6��face�����д���
	for (int fIndex = 0; fIndex < Cube::NEIGHBORS; fIndex++){
		// ��fIndex��Ե���
		int ref = Cube::FaceReflectFaceIndex(fIndex, fIndex);
		fef.fIndex = ref;
		// ��node����fIndex���ڵ���
		temp = node->faceNeighbor(fIndex);
		// temp�������node����fIndex���ڵ������
		// temp->children�����temp�ľ��ȱ�node��ϸ
		// IsBoundaryFace���ж�fIndex�Ƿ�Ϊ�߽��ϵ��棨���fIndex�Ǳ߽��ϵ��棬��ת�򲻳����ˣ�
		if (temp && temp->children && !IsBoundaryFace(node, fIndex, sDepth)){
			temp->processNodeFaces(temp, &fef, ref);
		}
		// �ӵ�ǰ���ȡnode
		else {
			RootInfo ri1, ri2;
			// �м�����������Ҫ����
			for (int j = 0; j < count; j++){
				// ÿ��������3����
				for (int k = 0; k < 3; k++){
					// �����ǰ���ڴ�������ᱻmarching cube���и�
					if (fIndex == Cube::FaceAdjacentToEdges(isoTri[j * 3 + k], isoTri[j * 3 + (k + 1) % 3]))	{
						if (GetRootIndex(node, isoTri[j * 3 + k], maxDepth, ri1) && GetRootIndex(node, isoTri[j * 3 + (k + 1) % 3], maxDepth, ri2)){
							// �߼����м�����������
							edges.push_back(pair<long long, long long>(ri1.key, ri2.key));
							// �㼯���м�����������
							iter = vertexCount.find(ri1.key);
							if (iter == vertexCount.end()){
								vertexCount[ri1.key].first = ri1;
								vertexCount[ri1.key].second = 0;
							}

							iter = vertexCount.find(ri2.key);
							if (iter == vertexCount.end()){
								vertexCount[ri2.key].first = ri2;
								vertexCount[ri2.key].second = 0;
							}
							// ͳ����ȡ�����
							vertexCount[ri1.key].second++;
							vertexCount[ri2.key].second--;
						}
						else
						{
							cerr << "Bad Edge 1: " << ri1.key << " " << ri2.key << endl;
						}
					}
				}
			}
		}
	}

	// ����ÿһ���߶����д���
	for (int i = 0; i < int(edges.size()); i++){
		// iter������vertex�е�һ�����������
		iter = vertexCount.find(edges[i].first);
		if (iter == vertexCount.end())
			cout << "[In PGROctree] Could not find vertex: " << edges[i].first << endl;
		// edges[i]�ĵ�һ���������ȳ��ȺͲ�Ϊ0������Ҫ��edges����������
		else if (vertexCount[edges[i].first].second) {
			RootInfo ri;
			// ��ȡ��Ӧ��rootinfo����ϸһ��ģ�
			if (GetRootPair(vertexCount[edges[i].first].first, maxDepth, ri)){
				iter = vertexCount.find(ri.key);
				if (iter == vertexCount.end())
					cout << "[In PGROctree] Vertex pair not in list" << endl;
				else{
					edges.push_back(pair<long long, long long>(ri.key, edges[i].first));
					vertexCount[ri.key].second++;
					vertexCount[edges[i].first].second--;
				}
			}
		}

		iter = vertexCount.find(edges[i].second);
		if (iter == vertexCount.end())
			cout << "[In PGROctree] Could not find vertex: " << edges[i].second << endl;
		else if (vertexCount[edges[i].second].second) {
			RootInfo ri;
			if (GetRootPair(vertexCount[edges[i].second].first, maxDepth, ri)){
				iter = vertexCount.find(ri.key);
				if (iter == vertexCount.end())
					cout << "[In PGROctree] Vertex pair not in list" << endl;
				else
				{
					edges.push_back(pair<long long, long long>(edges[i].second, ri.key));
					vertexCount[edges[i].second].second++;
					vertexCount[ri.key].second--;
				}
			}
		}
	}

}

int Octree::GetEdgeLoops(vector<pair<long long, long long>>& edges, vector<vector<pair<long long, long long>>>& loops){
	int loopSize = 0;
	long long frontIdx, backIdx;
	pair<long long, long long> e, temp;
	loops.clear();

	// ÿ��ѭ������һ��loop
	while (edges.size()){
		vector<pair<long long, long long>> front, back;
		// e�洢���ǵ�ǰ���ڴ����ı�
		e = edges[0];
		loops.resize(loopSize + 1);
		edges[0] = edges[edges.size() - 1];
		edges.pop_back();
		frontIdx = e.second;
		backIdx = e.first;
		for (int j = int(edges.size()) - 1; j >= 0; j--){
			// Ѱ����e������һ������ıߣ�����front����back����
			if (edges[j].first == frontIdx || edges[j].second == frontIdx){
				if (edges[j].first == frontIdx) temp = edges[j];
				else {
					temp.first = edges[j].second;
					temp.second = edges[j].first;
				}
				frontIdx = temp.second;
				front.push_back(temp);
				edges[j] = edges[edges.size() - 1];
				edges.pop_back();
				j = int(edges.size());
			}
			else if (edges[j].first == backIdx || edges[j].second == backIdx){
				if (edges[j].second == backIdx) temp = edges[j];
				else {
					temp.first = edges[j].second;
					temp.second = edges[j].first;
				}
				backIdx = temp.first;
				back.push_back(temp);
				edges[j] = edges[edges.size() - 1];
				edges.pop_back();
				j = int(edges.size());
			}
		}
		// ���γɵĻ�����loop
		for (int j = int(back.size()) - 1; j >= 0; j--){
			loops[loopSize].push_back(back[j]);
		}
		loops[loopSize].push_back(e);
		for (int j = 0; j < int(front.size()); j++){
			loops[loopSize].push_back(front[j]);
		}
		loopSize++;
	}
	return int(loops.size());
}

int Octree::GetMCIsoTriangles(OctNode* node, CoredVectorMeshData* mesh, unordered_map<long long, int>& boundaryRoots, unordered_map<long long, int>* interiorRoots, std::vector<Point>* interiorPositions, const int& offSet, const int& sDepth, bool addBarycenter, bool polygonMesh){
	int tris = 0;
	vector< pair< long long, long long> > edges;
	vector< vector< pair<long long, long long> > > edgeLoop;
	// ����Ӧ�ı߽�߶�����edges��
	GetMCIsoEdges(node, boundaryRoots, interiorRoots, sDepth, edges);
	// �������ɻ�
	GetEdgeLoops(edges, edgeLoop);
	for (int i = 0; i < int(edgeLoop.size()); i++){
		CoredPointIndex p;
		vector<CoredPointIndex> edgeIndices;
		for (int j = 0; j < int(edgeLoop[i].size()); j++){
			// ��ȡsecond��idx������edgeIndices
			if (!GetRootIndex(edgeLoop[i][j].first, boundaryRoots, interiorRoots, p)){
				cout << "Bad Point Index" << endl;
			}
			else
			{
				edgeIndices.push_back(p);
			}
		}
		tris += AddTriangles(mesh, edgeIndices, interiorPositions, offSet, addBarycenter, polygonMesh);
	}
	return tris;
}

int Octree::AddTriangles(CoredVectorMeshData* mesh, std::vector<CoredPointIndex>& edges, std::vector<Point >* interiorPositions, const int& offSet, bool addBarycenter, bool polygonMesh){
	if (polygonMesh){
		vector<CoredVertexIndex> vertices(edges.size());
		for (int i = 0; i < edges.size(); i++){
			vertices[i].idx = edges[i].index;
			vertices[i].inCore = edges[i].inCore;
		}
		mesh->addPolygon(vertices);
		return 1;
	}
	if (edges.size() > 3){
		bool isCoplanar = false;

		for (int i = 0; i < edges.size(); i++){
			for (int j = 0; j < i; j++){
				if ((i + 1) % edges.size() != j && (j + 1) % edges.size() != i){
					Point v1, v2;
					if (edges[i].inCore) {
						v1.x = mesh->inCorePoints[edges[i].index].x;
						v1.y = mesh->inCorePoints[edges[i].index].y;
						v1.z = mesh->inCorePoints[edges[i].index].z;
					}
					else {
						v1.x = (*interiorPositions)[edges[i].index - offSet].x;
						v1.y = (*interiorPositions)[edges[i].index - offSet].y;
						v1.z = (*interiorPositions)[edges[i].index - offSet].z;
					}
					if (edges[j].inCore){
						v2.x = mesh->inCorePoints[edges[j].index].x;
						v2.y = mesh->inCorePoints[edges[j].index].y;
						v2.z = mesh->inCorePoints[edges[j].index].z;
					}
					else {
						v2.x = (*interiorPositions)[edges[j].index - offSet].x;
						v2.y = (*interiorPositions)[edges[j].index - offSet].y;
						v2.z = (*interiorPositions)[edges[j].index - offSet].z;
					}
					// �����������һ��������ͬ������������㹲�㣿
					if (v1.x == v2.x || v1.y == v2.y || v1.z == v2.z)
						isCoplanar = true;
				}
			}
		}
		if (addBarycenter && isCoplanar){
			Point c;
			c.x = c.y = c.z = 0;
			for (int i = 0; i < int(edges.size()); i++){
				Point p;
				if (edges[i].inCore){
					p.x = mesh->inCorePoints[edges[i].index].x;
					p.y = mesh->inCorePoints[edges[i].index].y;
					p.z = mesh->inCorePoints[edges[i].index].z;
				}
				else {
					p.x = (*interiorPositions)[edges[i].index - offSet].x;
					p.y = (*interiorPositions)[edges[i].index - offSet].y;
					p.z = (*interiorPositions)[edges[i].index - offSet].z;
				}
				c.x += p.x; c.y += p.y; c.z += p.z;
			}
			c.x /= edges.size(); c.y /= edges.size(); c.z /= edges.size();
			int cIdx = mesh->addOutOfCorePoint(c);
			for (int i = 0; i < int(edges.size()); i++){
				vector<CoredVertexIndex> vertices(3);
				vertices[0].idx = edges[i].index;
				vertices[1].idx = edges[(i + 1) % edges.size()].index;
				vertices[2].idx = cIdx;
				vertices[0].inCore = edges[i].inCore;
				vertices[1].inCore = edges[(i + 1) % edges.size()].inCore;
				vertices[2].inCore = false;
				mesh->addPolygon(vertices);
			}
			return edges.size();
		}
		else {   // lwj
			Triangulation t;

			// Add the points to the triangulation
			for (int i = 0; i < int(edges.size()); i++){
				Point p;
				if (edges[i].inCore){
					p.x = mesh->inCorePoints[edges[i].index].x;
					p.y = mesh->inCorePoints[edges[i].index].y;
					p.z = mesh->inCorePoints[edges[i].index].z;
				}
				else {
					p.x = (*interiorPositions)[edges[i].index - offSet].x;
					p.y = (*interiorPositions)[edges[i].index - offSet].y;
					p.z = (*interiorPositions)[edges[i].index - offSet].z;
				}
				t.points.push_back(p);

			}

			// Create a fan triangulation
			for (int i = 1; i < int(edges.size()) - 1; i++){
				t.addTriangle(0, i, i + 1);
			}

			// Minimize
			while (1){
				int i;
				for (i = 0; i < int(t.edges.size()); i++){
					if (t.flipMinimize(i)) break;
				}
				if (i == t.edges.size())
					break;
			}

			// Add the triangles to the mesh
			for (int i = 0; i < int(t.triangles.size()); i++){
				vector<CoredVertexIndex> vertices(3);
				int idx[3];
				t.factor(i, idx[0], idx[1], idx[2]);
				for (int j = 0; j < 3; j++){
					vertices[j].idx = edges[idx[j]].index;
					vertices[j].inCore = edges[idx[j]].inCore;
				}
				mesh->addPolygon(vertices);
			}
		}
	}
	else if (edges.size() == 3){
		vector<CoredVertexIndex> vertices(3);
		for (int i = 0; i < 3; i++){
			vertices[i].idx = edges[i].index;
			vertices[i].inCore = edges[i].inCore;
		}
		mesh->addPolygon(vertices);
	}
	return int(edges.size()) - 2;
}

void FaceEdgesFunction::Function(const OctNode* node1){
	if (!node1->children && MarchingCubes::HasRoots(node1->mcIdx)){
		RootInfo ri1, ri2;
		unordered_map<long long, pair<RootInfo, int> >::iterator iter;
		int isoTri[3 * MarchingCubes::MAX_TRIANGLES];
		int count = MarchingCubes::AddTriangleIndices(node1->mcIdx, isoTri);

		for (int j = 0; j < count; j++){
			for (int k = 0; k < 3; k++){
				if (fIndex == Cube::FaceAdjacentToEdges(isoTri[j * 3 + k], isoTri[j * 3 + (k + 1) % 3])){
					if (Octree::GetRootIndex(node1, isoTri[j * 3 + k], maxDepth, ri1) && Octree::GetRootIndex(node1, isoTri[j * 3 + ((k + 1) % 3)], maxDepth, ri2)){
						edges->push_back(pair<long long, long long>(ri2.key, ri1.key));
						iter = vertexCount->find(ri1.key);
						if (iter == vertexCount->end()){
							(*vertexCount)[ri1.key].first = ri1;
							(*vertexCount)[ri1.key].second = 0;
						}
						iter = vertexCount->find(ri2.key);
						if (iter == vertexCount->end()){
							(*vertexCount)[ri2.key].first = ri2;
							(*vertexCount)[ri2.key].second = 0;
						}
						(*vertexCount)[ri1.key].second--;
						(*vertexCount)[ri2.key].second++;
					}
					else
					{
						cerr << "[In PGROctree] Bad Edge 1: " << ri1.key << " " << ri2.key << endl;
					}
				}
			}
		}

	}
}

void Octree::writePolygon2(CoredVectorMeshData* mesh, char* filename){
	Point translate;
	//01-09
	translate.x = (bb.blx + bb.xscale + bb.blx - maxScale) / 2.0;
	translate.y = (bb.bly + bb.yscale + bb.bly - maxScale) / 2.0;
	translate.z = (bb.blz + bb.zscale + bb.blz - maxScale) / 2.0;
	Point scale(bb.xscale*scaleFactor, bb.yscale*scaleFactor, bb.zscale*scaleFactor);
	PlyWritePolygons(filename, mesh, PLY_BINARY_NATIVE, translate, maxScale, NULL, 0);
}

void Octree::writePolygon(CoredVectorMeshData* mesh, string& filename){
	fstream out(filename, ios::out);
	if (!out.is_open()){
		cerr << "[In PGROctree] Cannot write into file: " << filename << endl;
		return;
	}

	mesh->resetIterator();
	int pointSize = mesh->inCorePoints.size();
	int polygonSize = mesh->polygonCount();
	// д���ļ�ͷ
	//out << "OFF" << endl;
	//out << pointSize << " " << polygonSize << " 0" << endl << endl;

	// д�������
	for (int i = 0; i < pointSize; i++){
		out << mesh->inCorePoints[i].x * maxScale + (bb.blx + bb.xscale + bb.blx - maxScale) / 2.0 << " "
			<< mesh->inCorePoints[i].y * maxScale + (bb.bly + bb.yscale + bb.bly - maxScale) / 2.0 << " "
			<< mesh->inCorePoints[i].z * maxScale + (bb.blz + bb.zscale + bb.blz - maxScale) / 2.0 << endl;
	}
	out.close();
}

float Octree::NonLinearUpdateWeightContribution(OctNode* node, const float* position, const float& weight){
	double x, dxdy, dx[3][3];	//dx�еĵ�һ��3��ά�ȣ��ڶ���3�ǽ���smooth�Ķ���ʽ����
	double width;
	float w;
	float center[3];
	node->centerAndWidth(center, w);
	width = w;
#if	NEW_SAMPLES_PER_NODE
	const double SAMPLE_SCALE = 1. / (0.125 * 0.125 + 0.75 * 0.75 + 0.125 * 0.125);
#endif
	for (int i = 0; i < 3; i++){
		x = (center[i] - position[i] - width) / width;
		dx[i][0] = 1.125 + 1.5 * x + 0.5 * x * x;	// B-spline
		x = (center[i] - position[i]) / width;
		dx[i][1] = 0.75 - x*x;
		dx[i][2] = 1.0 - dx[i][1] - dx[i][0];
#ifdef	NEW_SAMPLES_PER_NODE
		dx[i][0] *= SAMPLE_SCALE;
#endif
	}
	Neighbors& neighbors = neighborKey.setNeighbors(node);

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			dxdy = dx[0][i] * dx[1][j] * weight;
			for (int k = 0; k < 3; k++){
				if (neighbors.neighbors[i][j][k]){
					neighbors.neighbors[i][j][k]->centerWeightContribution += dxdy*dx[2][k];
				}
			}
		}
	}
	return 0;
}


void Octree::NonLinearSplatOrientedPoint(WeightedPoint& position, const float& normal_len, const int& splatDepth, const float& samplesPerNode, const int& minDepth, const int& maxDepth, const int& currentIdx){
	double dx;
	float nl;
	OctNode* temp;
	int cnt = 0;	//counter
	double width;
	float myCenter[3];
	float myWidth;
	myCenter[0] = myCenter[1] = myCenter[2] = 0.5;
	myWidth = 1.0;

	temp = &root;
	float pos[3] = { position.x, position.y, position.z };
	while (temp->Depth() < splatDepth){
		if (!temp->children){
			cout << "[In PGROctree] NonLinearSplatOrientedPoint error\n" << endl;
			return;
		}
		int cIndex = Cube::CornerIndex(myCenter, pos);
		temp = &temp->children[cIndex];
		myWidth /= 2;
		if (cIndex & 1)	myCenter[0] += myWidth / 2;
		else			myCenter[0] -= myWidth / 2;
		if (cIndex & 2)	myCenter[1] += myWidth / 2;
		else			myCenter[1] -= myWidth / 2;
		if (cIndex & 4)	myCenter[2] += myWidth / 2;
		else			myCenter[2] -= myWidth / 2;
	}

	float alpha, newDepth;
	NonLinearGetSampleDepthAndWeight(temp, &position, samplesPerNode, newDepth, alpha);

	if (newDepth < minDepth) newDepth = float(minDepth);
	if (newDepth > maxDepth) newDepth = float(maxDepth);

	int topDepth = int(ceil(newDepth));

	dx = 1.0 - (topDepth - newDepth);
	if (topDepth <= minDepth){
		topDepth = minDepth;
		dx = 1;
	}
	else if (topDepth > maxDepth){
		topDepth = maxDepth;
		dx = 1;
	}
	while (temp->Depth() > topDepth) temp = temp->parent;
	while (temp->Depth() < topDepth){
		if (!temp->children)
			temp->initChildren();
		int cIndex = Cube::CornerIndex(myCenter, pos);
		temp = &temp->children[cIndex];
		myWidth /= 2;
		if (cIndex & 1) myCenter[0] += myWidth / 2;
		else			myCenter[0] -= myWidth / 2;
		if (cIndex & 2) myCenter[1] += myWidth / 2;
		else			myCenter[1] -= myWidth / 2;
		if (cIndex & 4) myCenter[2] += myWidth / 2;
		else			myCenter[2] -= myWidth / 2;
	}
	width = 1.0 / (1 << temp->Depth());
	nl = normal_len * alpha / float(pow(width, 3))*dx;
	temp->nodeIdx.push_back(currentIdx);
	// editing
	NonLinearSplatOrientedPoint(temp, position, nl);
	if (fabs(1.0 - dx) > EPSILON){
		dx = float(1.0 - dx);
		temp = temp->parent;
		width = 1.0 / (1 << temp->Depth());

		nl = normal_len * alpha / pow(width, 3) * dx;
		NonLinearSplatOrientedPoint(temp, position, nl);
	}
}

void Octree::NonLinearGetSampleDepthAndWeight(OctNode* node, WeightedPoint* position, const float& samplesPerNode, float& depth, float& weight){
	OctNode* temp = node;
	float pos[3] = { position->x, position->y, position->z };
	weight = float(1.0) / NonLinearGetSampleWeight(temp, pos);
	position->weight = weight;

#if	NEW_SAMPLES_PER_NODE
	if (weight >= samplesPerNode) depth = float(temp->Depth() + log(weight / samplesPerNode) / log(double(1 << 2)));
#endif
	else{
		float oldAlpha, newAlpha;
		oldAlpha = newAlpha = weight;
#if NEW_SAMPLES_PER_NODE
		while (newAlpha < samplesPerNode && temp->parent){
#endif
			temp = temp->parent;
			oldAlpha = newAlpha;
			newAlpha = float(1.0) / NonLinearGetSampleWeight(temp, pos);
		}
#if NEW_SAMPLES_PER_NODE
		depth = float(temp->Depth() + log(newAlpha / samplesPerNode) / log(newAlpha / oldAlpha));
#endif
	}
	weight = float(pow(float(1 << 2), -double(depth)));
}

float Octree::NonLinearGetSampleWeight(OctNode* node, const float* position){
	float weight = 0;
	double x, dxdy, dx[3][3];
	Neighbors& neighbors = neighborKey.setNeighbors(node);

	double width;
	float center[3];
	float w;
	node->centerAndWidth(center, w);
	width = w;

	for (int i = 0; i < 3; i++){
		x = (center[i] - position[i] - width) / width;
		dx[i][0] = 1.125 + 1.5*x + 0.5*x*x;
		x = (center[i] - position[i]) / width;
		dx[i][1] = 0.75 - x*x;
		dx[i][2] = 1.0 - dx[i][1] - dx[i][0];
	}

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			dxdy = dx[0][i] * dx[1][j];
			for (int k = 0; k < 3; k++){
				if (neighbors.neighbors[i][j][k])
					weight += float(dxdy*dx[2][k] * neighbors.neighbors[i][j][k]->centerWeightContribution);
			}
		}
	}
	return 1.0 / weight;
}


int Octree::NonLinearSplatOrientedPoint(OctNode* node, const WeightedPoint& position, const float& normal_len){
	double x, dxdy, dxdydz, dx[3][3];
	Neighbors& neighbors = neighborKey.setNeighbors(node);

	double width;
	float center[3];
	float pos[3] = { position.x, position.y, position.z };
	float w;

	node->centerAndWidth(center, w);
	width = w;
	for (int i = 0; i < 3; i++){
		x = (center[i] - pos[i] - width) / width;
		dx[i][0] = 1.125 + 1.5*x + 0.5*x*x;
		x = (center[i] - pos[i]) / width;
		dx[i][1] = 0.75 - x*x;
		dx[i][2] = 1.0 - dx[i][1] - dx[i][0];
	}
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			dxdy = dx[0][i] * dx[1][j];
			for (int k = 0; k < 3; k++){
				if (neighbors.neighbors[i][j][k]){
					dxdydz = dxdy*dx[2][k];
					int idx = neighbors.neighbors[i][j][k]->normalIdx;
					if (idx < 0){
						float nl = 0;
						idx = neighbors.neighbors[i][j][k]->normalIdx = int(normalLengths.size());
						normalLengths.push_back(nl);
					}
					normalLengths[idx] += normal_len * dxdydz;
				}
			}
		}
	}
	return 0;
}

void Octree::ClipTree(void){
	OctNode* temp;
	temp = root.nextNode();
	while (temp){
		if (temp->children){
			int hasNormals = 0;
			for (int i = 0; i < Cube::CORNERS && !hasNormals; i++){
				hasNormals = HasNormals(&temp->children[i], EPSILON);
			}
			if (!hasNormals){
				temp->children = NULL;
			}
		}
		temp = root.nextNode(temp);
	}

}

int Octree::HasNormals(OctNode* node, const float& epsilon){
	int hasNormals = 0;
	if (node->normalIdx >= 0 && abs(normalLengths[node->normalIdx]) > epsilon){
		hasNormals = 1;
	}
	if (node->children) {
		for (int i = 0; i < Cube::CORNERS && !hasNormals; i++) {
			hasNormals |= HasNormals(&node->children[i], epsilon);
		}
	}
	return hasNormals;
}
