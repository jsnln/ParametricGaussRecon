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

#ifndef GEOMETRY_HEADER
#define GEOMETRY_HEADER
#include <math.h>
#include <vector>

#define _SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS 1
#include <unordered_map>

#include "BasicStructure.h"

using namespace std;
class Triangle{
public:
	double p[3][3];
	double Area(void) const{
		double v1[3],v2[3],v[3];
		for(int d=0;d<3;d++){
			v1[d]=p[1][d]-p[0][d];
			v2[d]=p[2][d]-p[0][d];
		}
		v[0]= v1[1]*v2[2]-v1[2]*v2[1];
		v[1]=-v1[0]*v2[2]+v1[2]*v2[0];
		v[2]= v1[0]*v2[1]-v1[1]*v2[0];
		return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/2;
	}
	double AspectRatio(void) const{
		double d=0;
		int i,j;
		for(i=0;i<3;i++){
			for(i=0;i<3;i++)
				for(j=0;j<3;j++){d+=(p[(i+1)%3][j]-p[i][j])*(p[(i+1)%3][j]-p[i][j]);}
		}
		return Area()/d;
	}

};

void CrossProduct(const Point& p1,const Point& p2,Point& p);
float getDistance(const Point& p1, const Point& p2);
float getDistance(const float coords[3], const Point& p);
float getDistance(const float p1[3], const float p2[3]);
float getDistance2(const Point& p1, const Point& p2);
float getDistance2(const float coords[3], const Point& p);

float getLength(const Point& p);
float getNormalLength(const Point& p);

float getLength(const Point& p);
float getSquareDistance( const Point& p1, const Point& p2);

class TriangleIndex{
public:
	int idx[3];
};

class TriangulationEdge
{
public:
	TriangulationEdge(void);
	int pIndex[2];
	int tIndex[2];
};

class TriangulationTriangle
{
public:
	TriangulationTriangle(void);
	int eIndex[3];
};

class Triangulation
{
public:

	std::vector<Point>		points;
	std::vector<TriangulationEdge>				edges;
	std::vector<TriangulationTriangle>			triangles;

	int factor(const int& tIndex,int& p1,int& p2,int& p3);
	double area(void);
	double area(const int& tIndex);
	double area(const int& p1,const int& p2,const int& p3);
	int flipMinimize(const int& eIndex);
	int addTriangle(const int& p1,const int& p2,const int& p3);

protected:
	unordered_map<long long,int> edgeMap;
	static long long EdgeIndex(const int& p1,const int& p2);
	double area(const Triangle& t);
};


void EdgeCollapse(const float& edgeRatio,std::vector<TriangleIndex>& triangles,std::vector< Point>& positions,std::vector<Point>* normals);

void TriangleCollapse(const float& edgeRatio,std::vector<TriangleIndex>& triangles,std::vector<Point>& positions,std::vector<Point>* normals);


#endif
