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

#include "Octnode.h"
#include <math.h>
#include <stdlib.h>

#define _SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS 1
#include <hash_map>

#include <iostream>
#include <algorithm>
#include "MarchingCubes.h"

using namespace std;
OctNode::OctNode(){

	for (int i = 0; i < Cube::CORNERS; i++)
		cornerGrid[i] = NULL;
	parent = NULL;	
	children = NULL;
	normalIdx = -1;
	mcIdx = -1;
	depth = offset[0] = offset[1] = offset[2] = 0;
	centerWeightContribution = 0;
	hasSample = false;
	gridStartIdx = -1;
	gridEndIdx = -1;
	barycenter.x = barycenter.y = barycenter.z = 0;
	nodeIdx.clear();
}

void OctNode::initChildren(){
	int d;
	int off[3];
	depthAndOffset( d, off);
	children = new OctNode[Cube::CORNERS];
	for( int i = 0; i < Cube::CORNERS; i++) {
		int currentOff[3];
		// ����i��ֵ�ֽ����ǰ�ӽڵ��λ��
		Cube::FactorCornerIndex( i, currentOff[0], currentOff[1], currentOff[2]);
		children[i].parent = this; children[i].children = NULL;
		int off2[3];

		off2[0] = (off[0] << 1) + currentOff[0];
		off2[1] = (off[1] << 1) + currentOff[1];
		off2[2] = (off[2] << 1) + currentOff[2];
		Index(d+1, off2, children[i].depth, children[i].offset);

	}

}

// ����λ�ƽ��н��룬Ϊ01000...��iλ��ʾ��i������֧�����ҷ�֧
void OctNode::depthAndOffset( int& depth, int off[3]) const{
//	cout << "depth and offset" << endl;
	depth = this->depth;
	off[0] = (int(offset[0])+1) & (~(1<<depth));
	off[1] = (int(offset[1])+1) & (~(1<<depth));
	off[2] = (int(offset[2])+1) & (~(1<<depth));
}

// ����ÿһ��node���б��룬Ϊ2^x+1...
void OctNode::Index( const int& depth, const int offset[3], int& d, int off[3]){
	d = depth;
	off[0] = ( (1<<depth) + offset[0] - 1);
	off[1] = ( (1<<depth) + offset[1] - 1);
	off[2] = ( (1<<depth) + offset[2] - 1);
}

long long OctNode::getCornerIndex( const int& childNo, const int& maxDepth){
	int d; 
	int off[3];
	depthAndOffset(d, off);
	int dir[3];
	Cube::FactorCornerIndex(childNo, dir[0], dir[1], dir[2]);
	for(int i = 0; i < 3; i++){
		off[i] = (off[i] + dir[i]) << (maxDepth - d + 1);
	}
	return (long long)off[0] | ((long long)off[1]) << DEPTH_LIMIT | ((long long)off[2]) << (DEPTH_LIMIT*2);
}

long long OctNode::getCornerIndex( const int& childNo, const int& maxDepth) const{
	int d; 
	int off[3];
	depthAndOffset(d, off);
	int dir[3];
	Cube::FactorCornerIndex(childNo, dir[0], dir[1], dir[2]);
	for(int i = 0; i < 3; i++){
		off[i] = (off[i] + dir[i]) << (maxDepth - d + 1);
	}
	return (long long)off[0] | ((long long)off[1]) << DEPTH_LIMIT | ((long long)off[2]) << (DEPTH_LIMIT*2);
}

int OctNode::leaves() {
	if( !children ) return 1;
	int c = 0;
	for( int i = 0; i < Cube::CORNERS; i++)
		c += children[i].leaves();
	return c;
}

int OctNode::nodes() {
	if( !children ) return 1;
	int c = 0;
	for( int i = 0; i < Cube::CORNERS; i++)
		c += children[i].nodes();
	return c + 1;
}

int OctNode::maxDepth(){
	if(!children) return 0;
	int c, d;
	for( int i = 0; i < Cube::CORNERS; i++){
		d = children[i].maxDepth();
		if( !i || d > c) c = d;
	}
	return c+1;
}

OctNode* OctNode::nextNode(OctNode* current /* = NULL */){
	if(!current) return this;
	else if( current->children ) {return &current->children[0];}
	else return nextBranch(current);
}

OctNode* OctNode::nextBranch(OctNode* current){
	if(!current->parent || current == this) return NULL;
	if( current - current->parent->children == Cube::CORNERS - 1) return nextBranch( current->parent);
	else return current + 1;
}

OctNode* OctNode::nextLeaf(OctNode* current){
	if(!current){	// ���ص�һ��Ҷ�ڵ�
		OctNode* temp = this;
		while(temp->children){
			temp = &temp->children[0];
		}
		return temp;
	}
	if(current->children) return current->nextLeaf();
	OctNode* temp = nextBranch(current);
	if(!temp) return NULL;
	else return temp->nextLeaf();
}

int OctNode::CompareForwardPointerDepths(const void* v1, const void* v2){
	const OctNode *n1, *n2;
	n1 = (*(const OctNode**)v1);
	n2 = (*(const OctNode**)v2);
	if( n1->depth != n2->depth) return n1->depth - n2->depth;
	while( n1->parent != n2->parent){
		n1 = n1->parent;
		n2 = n2->parent;
	}
	if( n1->offset[0] != n2->offset[0]) return n1->offset[0] - n2->offset[0];
	if( n1->offset[1] != n2->offset[1]) return n1->offset[1] - n2->offset[1];
	return n1->offset[2] - n2->offset[2];
}

int OctNode::Depth() const{
	return depth;
}

const OctNode* OctNode::__faceNeighbor(const int& dir,const int& off) const{
	if(!parent){return NULL;}
	int pIndex=int(this-parent->children);
	pIndex^=(1<<dir);
	if((pIndex & (1<<dir))==(off<<dir)){return &parent->children[pIndex];}
	//	if(!(((pIndex>>dir)^off)&1)){return &parent->children[pIndex];}
	else{
		const OctNode* temp=parent->__faceNeighbor(dir,off);
		if(!temp || !temp->children){return temp;}
		else{return &temp->children[pIndex];}
	}
}

OctNode* OctNode::__faceNeighbor(const int& dir, const int& off, const int& forceChildren){
	if(!parent) return NULL;
	int pIndex = int(this - parent->children);
	pIndex ^= (1<<dir);	// �ı�ָ�������ϵ���һλ
	// ������ڵ�ǰparent������ҵ�
	if ( (pIndex & (1 << dir)) == (off << dir)) return &parent->children[pIndex];
	// �Ҹ��ڵ���Ӧ�Ľ��
	else {
		OctNode* temp = parent->__faceNeighbor(dir, off, forceChildren);
		if(!temp) return NULL;
		if(!temp->children){
			if(forceChildren) temp->initChildren();
			else return temp;
		}
		return &temp->children[pIndex];
	}
}

const OctNode* OctNode::faceNeighbor(const int& faceIndex) const{
	return __faceNeighbor(faceIndex>>1,faceIndex&1);
}

OctNode* OctNode::faceNeighbor(const int& faceIndex, const int& forceChildren){
	return __faceNeighbor(faceIndex >> 1, faceIndex&1, forceChildren);
}

const OctNode* OctNode::__edgeNeighbor(const int& o,const int i[2],const int idx[2]) const{
	if(!parent){return NULL;}
	int pIndex=int(this-parent->children);
	int aIndex,x[3];

	Cube::FactorCornerIndex(pIndex,x[0],x[1],x[2]);
	aIndex=(~((i[0] ^ x[idx[0]]) | ((i[1] ^ x[idx[1]])<<1))) & 3;
	pIndex^=(7 ^ (1<<o));
	if(aIndex==1)	{	// I can get the neighbor from the parent's face adjacent neighbor
		const OctNode* temp=parent->__faceNeighbor(idx[0],i[0]);
		if(!temp || !temp->children){return NULL;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==2)	{	// I can get the neighbor from the parent's face adjacent neighbor
		const OctNode* temp=parent->__faceNeighbor(idx[1],i[1]);
		if(!temp || !temp->children){return NULL;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==0)	{	// I can get the neighbor from the parent
		return &parent->children[pIndex];
	}
	else if(aIndex==3)	{	// I can get the neighbor from the parent's edge adjacent neighbor
		const OctNode* temp=parent->__edgeNeighbor(o,i,idx);
		if(!temp || !temp->children){return temp;}
		else{return &temp->children[pIndex];}
	}
	else{return NULL;}
}

const OctNode* OctNode::edgeNeighbor(const int& edgeIndex) const{
	int idx[2],o,i[2];
	Cube::FactorEdgeIndex(edgeIndex,o,i[0],i[1]);
	switch(o){
	case 0:	idx[0]=1;	idx[1]=2;	break;
	case 1:	idx[0]=0;	idx[1]=2;	break;
	case 2:	idx[0]=0;	idx[1]=1;	break;
	};
	return __edgeNeighbor(o,i,idx);
}

SortedNodes::SortedNodes(){
	treeNodes = NULL;
	nodeCount = NULL;
	maxDepth = 0;
}

SortedNodes::~SortedNodes(){
	if( nodeCount ) delete[] nodeCount;
	if( treeNodes ) delete[] treeNodes;
	nodeCount = NULL;
	treeNodes = NULL;
}

void SortedNodes::set(OctNode& root, const int& setIndex){
	if(nodeCount) delete[] nodeCount;
	if(treeNodes) delete[] treeNodes;
	maxDepth = root.maxDepth() + 1;
	nodeCount = new int[maxDepth + 1];
	treeNodes = new OctNode*[root.nodes()];

	OctNode* temp = root.nextNode();
	int counter = 0;
	while( temp ){
		treeNodes[counter++] = temp;
		temp = root.nextNode(temp);
	}
	qsort( treeNodes, counter, sizeof( const OctNode*), OctNode::CompareForwardPointerDepths);
	for( int i = 0; i <= maxDepth; i++) nodeCount[i] = 0;
	for( int i = 0; i < counter; i++){
		nodeCount[treeNodes[i]->Depth() + 1]++;
	}
	for( int i = 1; i <= maxDepth; i++)
		nodeCount[i] += nodeCount[i-1];
}

void OctNode::processNodeFaces(OctNode* node, NodeAdjacencyFunction* F, const int& fIndex, const int& processCurrent /* = 1 */){
	if(processCurrent){ F->Function(node);}
	if(children){
		int c1, c2, c3, c4;
		Cube::FaceCorners(fIndex, c1, c2, c3, c4);
		__processNodeFaces(node, F, c1, c2, c3, c4);
	}
}

void OctNode::__processNodeFaces(OctNode* node, NodeAdjacencyFunction* F, const int& cIndex1, const int& cIndex2, const int& cIndex3, const int& cIndex4){
	F->Function(&children[cIndex1]);
	F->Function(&children[cIndex2]);
	F->Function(&children[cIndex3]);
	F->Function(&children[cIndex4]);
	if(children[cIndex1].children) {
		children[cIndex1].__processNodeFaces(node, F, cIndex1, cIndex2, cIndex3, cIndex4);
	}
	if(children[cIndex2].children){
		children[cIndex2].__processNodeFaces(node, F, cIndex1, cIndex2, cIndex3, cIndex4);
	}
	if(children[cIndex3].children){
		children[cIndex3].__processNodeFaces(node, F, cIndex1, cIndex2, cIndex3, cIndex4);
	}
	if(children[cIndex4].children){
		children[cIndex4].__processNodeFaces(node, F, cIndex1, cIndex2, cIndex3, cIndex4);
	}
}

void OctNode::centerAndWidth(float* c, float& width){
	int depth, off[3];
	depthAndOffset(depth, off);
	width = (float)1.0 / (1 << depth);
	for(int i = 0; i < 3; i++){
		c[i] = (0.5 + off[i]) * width;
	}

}

Neighbors::Neighbors(void){clear();}
void Neighbors::clear(void){
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				neighbors[i][j][k]=NULL;
			}
		}
	}
}

NeighborKey::NeighborKey(void){neighbors=NULL;}

NeighborKey::~NeighborKey(void){
	if(neighbors){delete[] neighbors;}
	neighbors=NULL;
}

void NeighborKey::set(const int& d){
	if(neighbors){delete[] neighbors;}
	neighbors=NULL;
	if(d<0){return;}
	neighbors=new Neighbors[d+1];
}

Neighbors& NeighborKey::setNeighbors(OctNode* node){
	//int d=node->depth();
	int d = node->Depth();
	if(node!=neighbors[d].neighbors[1][1][1]){
		neighbors[d].clear();

		if(!node->parent){neighbors[d].neighbors[1][1][1]=node;}
		else{
			int i,j,k,x1,y1,z1,x2,y2,z2;
			int idx=int(node-node->parent->children);
			Cube::FactorCornerIndex(  idx   ,x1,y1,z1);
			Cube::FactorCornerIndex((~idx)&7,x2,y2,z2);
			for(i=0;i<2;i++){
				for(j=0;j<2;j++){
					for(k=0;k<2;k++){
						neighbors[d].neighbors[x2+i][y2+j][z2+k]=&node->parent->children[Cube::CornerIndex(i,j,k)];
					}
				}
			}
			Neighbors& temp=setNeighbors(node->parent);

			// Set the neighbors from across the faces
			i=x1<<1;
			if(temp.neighbors[i][1][1]){
				if(!temp.neighbors[i][1][1]->children){temp.neighbors[i][1][1]->initChildren();}
				for(j=0;j<2;j++){for(k=0;k<2;k++){neighbors[d].neighbors[i][y2+j][z2+k]=&temp.neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];}}
			}
			j=y1<<1;
			if(temp.neighbors[1][j][1]){
				if(!temp.neighbors[1][j][1]->children){temp.neighbors[1][j][1]->initChildren();}
				for(i=0;i<2;i++){for(k=0;k<2;k++){neighbors[d].neighbors[x2+i][j][z2+k]=&temp.neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];}}
			}
			k=z1<<1;
			if(temp.neighbors[1][1][k]){
				if(!temp.neighbors[1][1][k]->children){temp.neighbors[1][1][k]->initChildren();}
				for(i=0;i<2;i++){for(j=0;j<2;j++){neighbors[d].neighbors[x2+i][y2+j][k]=&temp.neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];}}
			}

			// Set the neighbors from across the edges
			i=x1<<1;	j=y1<<1;
			if(temp.neighbors[i][j][1]){
				if(!temp.neighbors[i][j][1]->children){temp.neighbors[i][j][1]->initChildren();}
				for(k=0;k<2;k++){neighbors[d].neighbors[i][j][z2+k]=&temp.neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];}
			}
			i=x1<<1;	k=z1<<1;
			if(temp.neighbors[i][1][k]){
				if(!temp.neighbors[i][1][k]->children){temp.neighbors[i][1][k]->initChildren();}
				for(j=0;j<2;j++){neighbors[d].neighbors[i][y2+j][k]=&temp.neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];}
			}
			j=y1<<1;	k=z1<<1;
			if(temp.neighbors[1][j][k]){
				if(!temp.neighbors[1][j][k]->children){temp.neighbors[1][j][k]->initChildren();}
				for(i=0;i<2;i++){neighbors[d].neighbors[x2+i][j][k]=&temp.neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];}
			}

			// Set the neighbor from across the corner
			i=x1<<1;	j=y1<<1;	k=z1<<1;
			if(temp.neighbors[i][j][k]){
				if(!temp.neighbors[i][j][k]->children){temp.neighbors[i][j][k]->initChildren();}
				neighbors[d].neighbors[i][j][k]=&temp.neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
			}
		}
	}
	return neighbors[d];
}

Neighbors& NeighborKey::getNeighbors(OctNode* node){
	int d=node->Depth();
	if(node!=neighbors[d].neighbors[1][1][1]){
		neighbors[d].clear();

		if(!node->parent){neighbors[d].neighbors[1][1][1]=node;}
		else{
			int i,j,k,x1,y1,z1,x2,y2,z2;
			int idx=int(node-node->parent->children);
			Cube::FactorCornerIndex(  idx   ,x1,y1,z1);
			Cube::FactorCornerIndex((~idx)&7,x2,y2,z2);
			for(i=0;i<2;i++){
				for(j=0;j<2;j++){
					for(k=0;k<2;k++){
						neighbors[d].neighbors[x2+i][y2+j][z2+k]=&node->parent->children[Cube::CornerIndex(i,j,k)];
					}
				} 
			}
			Neighbors& temp=getNeighbors(node->parent);

			// Set the neighbors from across the faces
			i=x1<<1;
			if(temp.neighbors[i][1][1] && temp.neighbors[i][1][1]->children){
				for(j=0;j<2;j++){for(k=0;k<2;k++){neighbors[d].neighbors[i][y2+j][z2+k]=&temp.neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];}}
			}
			j=y1<<1;
			if(temp.neighbors[1][j][1] && temp.neighbors[1][j][1]->children){
				for(i=0;i<2;i++){for(k=0;k<2;k++){neighbors[d].neighbors[x2+i][j][z2+k]=&temp.neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];}}
			}
			k=z1<<1;
			if(temp.neighbors[1][1][k] && temp.neighbors[1][1][k]->children){
				for(i=0;i<2;i++){for(j=0;j<2;j++){neighbors[d].neighbors[x2+i][y2+j][k]=&temp.neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];}}
			}

			// Set the neighbors from across the edges
			i=x1<<1;	j=y1<<1;
			if(temp.neighbors[i][j][1] && temp.neighbors[i][j][1]->children){
				for(k=0;k<2;k++){neighbors[d].neighbors[i][j][z2+k]=&temp.neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];}
			}
			i=x1<<1;	k=z1<<1;
			if(temp.neighbors[i][1][k] && temp.neighbors[i][1][k]->children){
				for(j=0;j<2;j++){neighbors[d].neighbors[i][y2+j][k]=&temp.neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];}
			}
			j=y1<<1;	k=z1<<1;
			if(temp.neighbors[1][j][k] && temp.neighbors[1][j][k]->children){
				for(i=0;i<2;i++){neighbors[d].neighbors[x2+i][j][k]=&temp.neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];}
			}

			// Set the neighbor from across the corner
			i=x1<<1;	j=y1<<1;	k=z1<<1;
			if(temp.neighbors[i][j][k] && temp.neighbors[i][j][k]->children){
				neighbors[d].neighbors[i][j][k]=&temp.neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
			}
		}
	}
	return neighbors[node->Depth()];
}
