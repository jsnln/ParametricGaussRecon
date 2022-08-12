/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
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

#ifndef MARCHING_CUBES_HEADER
#define MARCHING_CUBES_HEADER

#include "Cube.h"
#include "Geometry.h"

class MarchingCubes{
	static double Interpolate(const double& v1,const double& v2);
	static void SetVertex(const int& e,const double values[Cube::CORNERS],const double& iso);
	static int GetFaceIndex(const double values[Cube::CORNERS],const double& iso,const int& faceIndex);

	static float Interpolate(const float& v1,const float& v2);
	static void SetVertex(const int& e,const float values[Cube::CORNERS],const float& iso);
	static int GetFaceIndex(const float values[Cube::CORNERS],const float& iso,const int& faceIndex);

	static int GetFaceIndex(const int& mcIndex,const int& faceIndex);
public:
	const static unsigned int MAX_TRIANGLES=5;
	static const int edgeMask[1<<Cube::CORNERS];
	static const int triangles[1<<Cube::CORNERS][3*MAX_TRIANGLES+1];
	static const int cornerMap[Cube::CORNERS];
	static double vertexList[Cube::EDGES][3];

	static int AddTriangleIndices(const int& mcIndex,int* triangles);

	static int GetIndex(const double values[Cube::CORNERS],const double& iso);
	static int IsAmbiguous(const double v[Cube::CORNERS],const double& isoValue,const int& faceIndex);
	static int HasRoots(const double v[Cube::CORNERS],const double& isoValue);
	static int HasRoots(const double v[Cube::CORNERS],const double& isoValue,const int& faceIndex);
	static int AddTriangles(const double v[Cube::CORNERS],const double& isoValue,Triangle* triangles);
	static int AddTriangleIndices(const double v[Cube::CORNERS],const double& isoValue,int* triangles);

	static int GetIndex(const float values[Cube::CORNERS],const float& iso);
	static int IsAmbiguous(const float v[Cube::CORNERS],const float& isoValue,const int& faceIndex);
	static int HasRoots(const float v[Cube::CORNERS],const float& isoValue);
	static int HasRoots(const float v[Cube::CORNERS],const float& isoValue,const int& faceIndex);
	static int AddTriangles(const float v[Cube::CORNERS],const float& isoValue,Triangle* triangles);
	static int AddTriangleIndices(const float v[Cube::CORNERS],const float& isoValue,int* triangles);

	static int IsAmbiguous(const int& mcIndex,const int& faceIndex);
	static int HasRoots(const int& mcIndex);
	static int HasFaceRoots(const int& mcIndex,const int& faceIndex);
	static int HasEdgeRoots(const int& mcIndex,const int& edgeIndex);
};

#endif 

