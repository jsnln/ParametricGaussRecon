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

#ifndef	CUBE_HEADER
#define CUBE_HEADER

#include "BasicStructure.h"

class Cube{
public:
	enum NEIGHBOR_DIRE{X_FRONT, X_BACK, Y_FRONT, Y_BACK, Z_FRONT, Z_BACK};
	const static int CORNERS = 8;
	const static int EDGES = 12;
	const static int NEIGHBORS = 6;
	const static int CornerAdjacentMap[CORNERS][NEIGHBORS];

	static int CornerIndex( const int& x, const int& y, const int& z);
	static int CornerIndex( const float* center, const float* position);
	static int CornerIndex( const float* center, const WeightedPoint* np);
	static void FactorCornerIndex( const int& idx,int& x,int& y,int& z);
	static void FactorFaceIndex		(const int& idx,int& x,int &y,int& z);
	static void FactorFaceIndex		(const int& idx,int& dir,int& offSet);
	static int  AntipodalCornerIndex	(const int& idx);
	static void EdgeCorners(const int& idx,int& c1,int &c2);
	static void FactorEdgeIndex		(const int& idx,int& orientation,int& i,int &j);
	static int  EdgeIndex			(const int& orientation,const int& i,const int& j);
	static void FacesAdjacentToEdge	(const int& eIndex,int& f1Index,int& f2Index);
	static int  FaceIndex			(const int& x,const int& y,const int& z);
	static int  FaceReflectEdgeIndex	(const int& idx,const int& faceIndex);
	static int	EdgeReflectEdgeIndex	(const int& edgeIndex);
	static int	FaceReflectFaceIndex	(const int& idx,const int& faceIndex);
	static int  FaceAdjacentToEdges	(const int& eIndex1,const int& eIndex2);
	static void FaceCorners(const int& idx,int& c1,int &c2,int& c3,int& c4);

};
#endif
