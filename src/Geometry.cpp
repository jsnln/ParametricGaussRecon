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

#include "Geometry.h"

float getDistance(const Point& p1, const Point& p2){
	return sqrt(getDistance2(p1, p2));
}

float getDistance2(const Point& p1, const Point& p2){
	return pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2);
}

float getDistance2(const float coords[3], const Point& p){
	return pow(coords[0] - p.x, 2) + pow( coords[1] - p.y, 2) + pow( coords[2] - p.z, 2);
}

float getDistance(const float coords[3], const Point& p){
	return sqrt(pow(coords[0] - p.x, 2) + pow(coords[1] - p.y, 2) + pow(coords[2] - p.z, 2));
}

float getDistance(const float p1[3], const float p2[3]){
	return sqrt( pow( p1[0] - p2[0], 2) + pow( p1[1] - p2[1], 2) + pow( p1[2] - p2[2], 2));
}

float getLength(const Point& p){
	return sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2));
}

float getSquareDistance(const Point& p1, const Point& p2){
	return pow( p1.x - p2.x, 2) + pow( p1.y - p2.y, 2) + pow( p1.z - p2.z, 2);
}

void CrossProduct(const Point& p1,const Point& p2,Point& p){
	p.x = p1.y*p2.z - p1.z*p2.y;
	p.y = -p1.x*p2.z + p1.z*p2.x;
	p.z = p1.x*p2.y - p1.y*p2.x;
}

TriangulationEdge::TriangulationEdge(void){pIndex[0]=pIndex[1]=tIndex[0]=tIndex[1]=-1;}
TriangulationTriangle::TriangulationTriangle(void){eIndex[0]=eIndex[1]=eIndex[2]=-1;}



long long Triangulation::EdgeIndex(const int& p1,const int& p2){
	if(p1>p2)	{return ((long long)(p1)<<32) | ((long long)(p2));}
	else		{return ((long long)(p2)<<32) | ((long long)(p1));}
}

int Triangulation::factor(const int& tIndex,int& p1,int& p2,int & p3){
	if(triangles[tIndex].eIndex[0]<0 || triangles[tIndex].eIndex[1]<0 || triangles[tIndex].eIndex[2]<0){return 0;}
	if(edges[triangles[tIndex].eIndex[0]].tIndex[0]==tIndex){p1=edges[triangles[tIndex].eIndex[0]].pIndex[0];}
	else													{p1=edges[triangles[tIndex].eIndex[0]].pIndex[1];}
	if(edges[triangles[tIndex].eIndex[1]].tIndex[0]==tIndex){p2=edges[triangles[tIndex].eIndex[1]].pIndex[0];}
	else													{p2=edges[triangles[tIndex].eIndex[1]].pIndex[1];}
	if(edges[triangles[tIndex].eIndex[2]].tIndex[0]==tIndex){p3=edges[triangles[tIndex].eIndex[2]].pIndex[0];}
	else													{p3=edges[triangles[tIndex].eIndex[2]].pIndex[1];}
	return 1;
}

double Triangulation::area(const int& p1,const int& p2,const int& p3){
	Point q1,q2,q;
	q1.x = points[p2].x - points[p1].x;
	q1.y = points[p2].y - points[p1].y;
	q1.z = points[p2].z - points[p1].z;
	q2.x = points[p3].x - points[p1].x;
	q2.y = points[p3].y - points[p1].y;
	q2.z = points[p3].z - points[p1].z;
	CrossProduct(q1,q2,q);
	return getLength(q);
}

double Triangulation::area(const int& tIndex){
	int p1,p2,p3;
	factor(tIndex,p1,p2,p3);
	return area(p1,p2,p3);
}

double Triangulation::area(void){
	double a=0;
	for(int i=0;i<int(triangles.size());i++){a+=area(i);}
	return a;
}

int Triangulation::addTriangle(const int& p1,const int& p2,const int& p3){
	unordered_map<long long,int>::iterator iter;
	int tIdx,eIdx,p[3];
	p[0]=p1;
	p[1]=p2;
	p[2]=p3;
	triangles.push_back(TriangulationTriangle());
	tIdx=int(triangles.size())-1;

	for(int i=0;i<3;i++)
	{
		long long e = EdgeIndex(p[i],p[(i+1)%3]);
		iter=edgeMap.find(e);
		if(iter==edgeMap.end())
		{
			TriangulationEdge edge;
			edge.pIndex[0]=p[i];
			edge.pIndex[1]=p[(i+1)%3];
			edges.push_back(edge);
			eIdx=int(edges.size())-1;
			edgeMap[e]=eIdx;
			edges[eIdx].tIndex[0]=tIdx;
		}
		else{
			eIdx=edgeMap[e];
			if(edges[eIdx].pIndex[0]==p[i]){
				if(edges[eIdx].tIndex[0]<0){edges[eIdx].tIndex[0]=tIdx;}
				else{printf("Edge Triangle in use 1\n");return 0;}
			}
			else{
				if(edges[eIdx].tIndex[1]<0){edges[eIdx].tIndex[1]=tIdx;}
				else{printf("Edge Triangle in use 2\n");return 0;}
			}

		}
		triangles[tIdx].eIndex[i]=eIdx;
	}
	return tIdx;
}


int Triangulation::flipMinimize(const int& eIndex){
	double oldArea,newArea;
	int oldP[3],oldQ[3],newP[3],newQ[3];
	TriangulationEdge newEdge;

	if(edges[eIndex].tIndex[0]<0 || edges[eIndex].tIndex[1]<0){return 0;}

	// ��ȡ��edge���ڵ����������μ�������
	if(!factor(edges[eIndex].tIndex[0],oldP[0],oldP[1],oldP[2])){return 0;}
	if(!factor(edges[eIndex].tIndex[1],oldQ[0],oldQ[1],oldQ[2])){return 0;}


	oldArea=area(oldP[0],oldP[1],oldP[2])+area(oldQ[0],oldQ[1],oldQ[2]);
	int idxP,idxQ;
	for(idxP=0;idxP<3;idxP++){
		int i;
		for(i=0;i<3;i++){if(oldP[idxP]==oldQ[i]){break;}}
		if(i==3){break;}
	}
	for(idxQ=0;idxQ<3;idxQ++){
		int i;
		for(i=0;i<3;i++){if(oldP[i]==oldQ[idxQ]){break;}}
		if(i==3){break;}
	}
	if(idxP==3 || idxQ==3){return 0;}
	newP[0]=oldP[idxP];
	newP[1]=oldP[(idxP+1)%3];
	newP[2]=oldQ[idxQ];
	newQ[0]=oldQ[idxQ];
	newQ[1]=oldP[(idxP+2)%3];
	newQ[2]=oldP[idxP];

	newArea=area(newP[0],newP[1],newP[2])+area(newQ[0],newQ[1],newQ[2]);
	if(oldArea<=newArea){return 0;}

	// Remove the entry in the hash_table for the old edge
	edgeMap.erase(EdgeIndex(edges[eIndex].pIndex[0],edges[eIndex].pIndex[1]));
	// Set the new edge so that the zero-side is newQ
	edges[eIndex].pIndex[0]=newP[0];
	edges[eIndex].pIndex[1]=newQ[0];
	// Insert the entry into the hash_table for the new edge
	edgeMap[EdgeIndex(newP[0],newQ[0])]=eIndex;
	// Update the triangle information
	for(int i=0;i<3;i++){
		int idx;
		idx=edgeMap[EdgeIndex(newQ[i],newQ[(i+1)%3])];
		triangles[edges[eIndex].tIndex[0]].eIndex[i]=idx;
		if(idx!=eIndex){
			if(edges[idx].tIndex[0]==edges[eIndex].tIndex[1]){edges[idx].tIndex[0]=edges[eIndex].tIndex[0];}
			if(edges[idx].tIndex[1]==edges[eIndex].tIndex[1]){edges[idx].tIndex[1]=edges[eIndex].tIndex[0];}
		}

		idx=edgeMap[EdgeIndex(newP[i],newP[(i+1)%3])];
		triangles[edges[eIndex].tIndex[1]].eIndex[i]=idx;
		if(idx!=eIndex){
			if(edges[idx].tIndex[0]==edges[eIndex].tIndex[0]){edges[idx].tIndex[0]=edges[eIndex].tIndex[1];}
			if(edges[idx].tIndex[1]==edges[eIndex].tIndex[0]){edges[idx].tIndex[1]=edges[eIndex].tIndex[1];}
		}
	}
	return 1;
}

void EdgeCollapse(const float& edgeRatio,std::vector<TriangleIndex>& triangles,std::vector< Point >& positions,std::vector< Point >* normals){
	int i,j,*remapTable,*pointCount,idx[3];
	Point p[3],q[2],c;
	double d[3],a;
	double Ratio=12.0/sqrt(3.0);	// (Sum of Squares Length / Area) for and equilateral triangle

	remapTable=new int[positions.size()];
	pointCount=new int[positions.size()];
	for(i=0;i<int(positions.size());i++){
		remapTable[i]=i;
		pointCount[i]=1;
	}
	for(i=int(triangles.size()-1);i>=0;i--){
		for(j=0;j<3;j++){
			idx[j]=triangles[i].idx[j];
			while(remapTable[idx[j]]<idx[j]){idx[j]=remapTable[idx[j]];}
		}
		if(idx[0]==idx[1] || idx[0]==idx[2] || idx[1]==idx[2]){
			triangles[i]=triangles[triangles.size()-1];
			triangles.pop_back();
			continue;
		}
		for( j = 0; j < 3; j++){
			p[j].x = positions[idx[j]].x / pointCount[idx[j]];
			p[j].y = positions[idx[j]].y / pointCount[idx[j]];
			p[j].z = positions[idx[j]].z / pointCount[idx[j]];
		}
		q[0].x = p[1].x - p[0].x;
		q[1].x = p[2].x - p[0].x;
		d[0] = getSquareDistance( p[0], p[1]);
		q[0].y = p[1].y - p[0].y;
		q[1].y = p[2].y - p[0].y;
		d[1] = getSquareDistance( p[1], p[2]);
		q[0].z = p[1].z - p[0].z;
		q[1].z = p[2].z - p[0].z;
		d[2] = getSquareDistance( p[2], p[0]);

		CrossProduct(q[0],q[1],c);
		a = getLength(c)/2;

		if((d[0]+d[1]+d[2])*edgeRatio > a*Ratio){
			// Find the smallest edge
			j=0;
			if(d[1]<d[j]){j=1;}
			if(d[2]<d[j]){j=2;}

			int idx1,idx2;
			if(idx[j]<idx[(j+1)%3]){
				idx1=idx[j];
				idx2=idx[(j+1)%3];
			}
			else{
				idx2=idx[j];
				idx1=idx[(j+1)%3];
			}
			positions[idx1].x += positions[idx2].x;
			positions[idx1].y += positions[idx2].y;
			positions[idx1].z += positions[idx2].z;
			if(normals){
				(*normals)[idx1].x += (*normals)[idx2].x;
				(*normals)[idx1].y += (*normals)[idx2].y;
				(*normals)[idx1].z += (*normals)[idx2].z;
			}
			pointCount[idx1]+=pointCount[idx2];
			remapTable[idx2]=idx1;
			triangles[i]=triangles[triangles.size()-1];
			triangles.pop_back();
		}
	}
	int pCount=0;
	for(i=0;i<int(positions.size());i++){
		//for(j=0;j<3;j++){positions[i].coords[j]/=pointCount[i];}
		positions[i].x /= pointCount[i];
		positions[i].y /= pointCount[i];
		positions[i].z /= pointCount[i];
		if(normals){
			float l = float(getLength((*normals)[i]));
			// for(j=0;j<3;j++){(*normals)[i].coords[j]/=l;}
			(*normals)[i].x /= l;
			(*normals)[i].y /= l;
			(*normals)[i].z /= l;
		}
		if(remapTable[i]==i){ // If vertex i is being used
			positions[pCount]=positions[i];
			if(normals){(*normals)[pCount]=(*normals)[i];}
			pointCount[i]=pCount;
			pCount++;
		}
	}
	positions.resize(pCount);
	for(i=int(triangles.size()-1);i>=0;i--){
		for(j=0;j<3;j++){
			idx[j]=triangles[i].idx[j];
			while(remapTable[idx[j]]<idx[j]){idx[j]=remapTable[idx[j]];}
			triangles[i].idx[j]=pointCount[idx[j]];
		}
		if(idx[0]==idx[1] || idx[0]==idx[2] || idx[1]==idx[2]){
			triangles[i]=triangles[triangles.size()-1];
			triangles.pop_back();
		}
	}

	delete[] pointCount;
	delete[] remapTable;
}

void TriangleCollapse(const float& edgeRatio,std::vector<TriangleIndex>& triangles,std::vector<Point>& positions,std::vector<Point>* normals){
	int i,j,*remapTable,*pointCount,idx[3];
	Point p[3],q[2],c;
	double d[3],a;
	double Ratio=12.0/sqrt(3.0);	// (Sum of Squares Length / Area) for and equilateral triangle

	remapTable=new int[positions.size()];
	pointCount=new int[positions.size()];
	for(i=0;i<int(positions.size());i++){
		remapTable[i]=i;
		pointCount[i]=1;
	}
	for(i=int(triangles.size()-1);i>=0;i--){
		for(j=0;j<3;j++){
			idx[j]=triangles[i].idx[j];
			while(remapTable[idx[j]]<idx[j]){idx[j]=remapTable[idx[j]];}
		}
		if(idx[0]==idx[1] || idx[0]==idx[2] || idx[1]==idx[2]){
			triangles[i]=triangles[triangles.size()-1];
			triangles.pop_back();
			continue;
		}
		for(j=0;j<3;j++){
			p[j].x = positions[idx[j]].x / pointCount[idx[j]];
			p[j].y =positions[idx[j]].y / pointCount[idx[j]];
			p[j].z =positions[idx[j]].z / pointCount[idx[j]];
		}

		q[0].x = p[1].x - p[0].x;
		q[1].x = p[2].x - p[0].x;
		d[0] = getSquareDistance(p[0], p[1]);
		q[0].y = p[1].y - p[0].y;
		q[1].y = p[2].y - p[0].y;
		d[1] = getSquareDistance(p[1], p[2]);
		q[0].z = p[1].z - p[0].z;
		q[1].z = p[2].z - p[0].z;
		d[2] = getSquareDistance(p[2], p[0]);
		CrossProduct(q[0],q[1],c);
		a = getLength(c)/2;

		if((d[0]+d[1]+d[2])*edgeRatio > a*Ratio){
			// Find the smallest edge
			j=0;
			if(d[1]<d[j]){j=1;}
			if(d[2]<d[j]){j=2;}

			int idx1,idx2,idx3;
			if(idx[0]<idx[1]){
				if(idx[0]<idx[2]){
					idx1=idx[0];
					idx2=idx[2];
					idx3=idx[1];
				}
				else{
					idx1=idx[2];
					idx2=idx[0];
					idx3=idx[1];
				}
			}
			else{
				if(idx[1]<idx[2]){
					idx1=idx[1];
					idx2=idx[2];
					idx3=idx[0];
				}
				else{
					idx1=idx[2];
					idx2=idx[1];
					idx3=idx[0];
				}
			}
			positions[idx1].x += positions[idx2].x + positions[idx3].x;
			positions[idx1].y += positions[idx2].y + positions[idx3].y;
			positions[idx1].z += positions[idx2].z + positions[idx3].z;
			if(normals){
				(*normals)[idx1].x += (*normals)[idx2].x + (*normals)[idx3].x;
				(*normals)[idx1].y += (*normals)[idx2].y + (*normals)[idx3].y;
				(*normals)[idx1].z += (*normals)[idx2].z + (*normals)[idx3].z;
			}
			pointCount[idx1]+=pointCount[idx2]+pointCount[idx3];
			remapTable[idx2]=idx1;
			remapTable[idx3]=idx1;
			triangles[i]=triangles[triangles.size()-1];
			triangles.pop_back();
		}
	}
	int pCount=0;
	for(i=0;i<int(positions.size());i++){
		//for(j=0;j<3;j++){positions[i].coords[j]/=pointCount[i];}
		positions[i].x /= pointCount[i];
		positions[i].y /= pointCount[i];
		positions[i].z /= pointCount[i];

		if(normals){
			float l=float(getLength((*normals)[i]));
			//for(j=0;j<3;j++){(*normals)[i].coords[j]/=l;}
			(*normals)[i].x /= l;
			(*normals)[i].y /= l;

		}
		if(remapTable[i]==i){ // If vertex i is being used
			positions[pCount]=positions[i];
			if(normals){(*normals)[pCount]=(*normals)[i];}
			pointCount[i]=pCount;
			pCount++;
		}
	}
	positions.resize(pCount);
	for(i=int(triangles.size()-1);i>=0;i--){
		for(j=0;j<3;j++){
			idx[j]=triangles[i].idx[j];
			while(remapTable[idx[j]]<idx[j]){idx[j]=remapTable[idx[j]];}
			triangles[i].idx[j]=pointCount[idx[j]];
		}
		if(idx[0]==idx[1] || idx[0]==idx[2] || idx[1]==idx[2]){
			triangles[i]=triangles[triangles.size()-1];
			triangles.pop_back();
		}
	}
	delete[] pointCount;
	delete[] remapTable;
}

