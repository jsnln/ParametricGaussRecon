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

#include "Cube.h"

const int Cube::CornerAdjacentMap[CORNERS][NEIGHBORS] = {
	{  1, -1,  2, -1,  4, -1},
	{ -1,  0,  3, -1,  5, -1},
	{  3, -1, -1,  0,  6, -1},
	{ -1,  2, -1,  1,  7, -1},
	{  5, -1,  6, -1, -1,  0},
	{ -1,  4,  7, -1, -1,  1},
	{  7, -1, -1,  4, -1,  2},
	{ -1,  6, -1,  5, -1,  3}
};

int Cube::CornerIndex(const int& x, const int& y, const int& z){
	return (z << 2) | (y << 1) | x;
}

int Cube::CornerIndex(const float* center, const float* position){
	int cIndex = 0;
	if( position[0] > center[0] ) cIndex |= 1;
	if( position[1] > center[1] ) cIndex |= 2;
	if( position[2] > center[2] ) cIndex |= 4;
	return cIndex;
}

int Cube::CornerIndex(const float* center, const WeightedPoint* np){
	int cIndex = 0;
	if(np->x > center[0]) cIndex |= 1;
	if(np->y > center[1]) cIndex |= 2;
	if(np->z > center[2]) cIndex |= 4;
	return cIndex;
}

void Cube::FactorCornerIndex(const int& idx,int& x,int& y,int& z){
	x=(idx>>0)%2;
	y=(idx>>1)%2;
	z=(idx>>2)%2;
}

int Cube::AntipodalCornerIndex(const int& idx){
	int x,y,z;
	FactorCornerIndex(idx,x,y,z);
	return CornerIndex((x+1)%2,(y+1)%2,(z+1)%2);
}

void Cube::FactorFaceIndex(const int& idx,int& x,int& y,int& z){
	x=y=z=0;
	switch(idx){
	case 0:		x=-1;	break;
	case 1:		x= 1;	break;
	case 2:		y=-1;	break;
	case 3:		y= 1;	break;
	case 4:		z=-1;	break;
	case 5:		z= 1;	break;
	};
}

void Cube::EdgeCorners(const int& idx,int& c1,int& c2){
	int orientation,i1,i2;
	FactorEdgeIndex(idx,orientation,i1,i2);
	switch(orientation){
	case 0:
		c1=CornerIndex(0,i1,i2);
		c2=CornerIndex(1,i1,i2);
		break;
	case 1:
		c1=CornerIndex(i1,0,i2);
		c2=CornerIndex(i1,1,i2);
		break;
	case 2:
		c1=CornerIndex(i1,i2,0);
		c2=CornerIndex(i1,i2,1);
		break;
	};
}

void Cube::FactorEdgeIndex(const int& idx,int& orientation,int& i,int &j){
	orientation=idx>>2;
	i=idx&1;
	j=(idx&2)>>1;
}

int Cube::EdgeIndex(const int& orientation,const int& i,const int& j){
	return (i | (j<<1))|(orientation<<2);
}

int Cube::FaceIndex(const int& x,const int& y,const int& z){
	if		(x<0)	{return  0;}
	else if	(x>0)	{return  1;}
	else if	(y<0)	{return  2;}
	else if	(y>0)	{return  3;}
	else if	(z<0)	{return  4;}
	else if	(z>0)	{return  5;}
	else			{return -1;}
}


void Cube::FacesAdjacentToEdge(const int& eIndex,int& f1Index,int& f2Index){
	int orientation,i1,i2;
	FactorEdgeIndex(eIndex,orientation,i1,i2);
	i1<<=1;
	i2<<=1;
	i1--;
	i2--;
	switch(orientation){
	case 0:
		f1Index=FaceIndex( 0,i1, 0);
		f2Index=FaceIndex( 0, 0,i2);
		break;
	case 1:
		f1Index=FaceIndex(i1, 0, 0);
		f2Index=FaceIndex( 0, 0,i2);
		break;
	case 2:
		f1Index=FaceIndex(i1, 0, 0);
		f2Index=FaceIndex( 0,i2, 0);
		break;
	};
}

int Cube::FaceReflectEdgeIndex(const int& idx,const int& faceIndex){
	int orientation=faceIndex/2;
	int o,i,j;
	FactorEdgeIndex(idx,o,i,j);
	if(o==orientation){return idx;}
	switch(orientation){
	case 0:	return EdgeIndex(o,(i+1)%2,j);
	case 1:
		switch(o){
		case 0:	return EdgeIndex(o,(i+1)%2,j);
		case 2:	return EdgeIndex(o,i,(j+1)%2);
		};
	case 2:	return EdgeIndex(o,i,(j+1)%2);
	};
	return -1;
}

int	Cube::EdgeReflectEdgeIndex(const int& edgeIndex){
	int o,i1,i2;
	FactorEdgeIndex(edgeIndex,o,i1,i2);
	return Cube::EdgeIndex(o,(i1+1)%2,(i2+1)%2);
}

int Cube::FaceReflectFaceIndex(const int& idx,const int& faceIndex){
	if(idx/2!=faceIndex/2){return idx;}
	else{
		if(idx%2)	{return idx-1;}
		else		{return idx+1;}
	}
}

// ��eIndex1��eIndex2�����ڵ���
int Cube::FaceAdjacentToEdges(const int& eIndex1,const int& eIndex2){
	int f1,f2,g1,g2;
	FacesAdjacentToEdge(eIndex1,f1,f2);
	FacesAdjacentToEdge(eIndex2,g1,g2);
	if(f1==g1 || f1==g2){return f1;}
	if(f2==g1 || f2==g2){return f2;}
	return -1;
}

void Cube::FaceCorners(const int& idx,int& c1,int& c2,int& c3,int& c4){
	int i=idx%2;
	switch(idx/2){
	case 0:
		c1=CornerIndex(i,0,0);
		c2=CornerIndex(i,1,0);
		c3=CornerIndex(i,0,1);
		c4=CornerIndex(i,1,1);
		return;
	case 1:
		c1=CornerIndex(0,i,0);
		c2=CornerIndex(1,i,0);
		c3=CornerIndex(0,i,1);
		c4=CornerIndex(1,i,1);
		return;
	case 2:
		c1=CornerIndex(0,0,i);
		c2=CornerIndex(1,0,i);
		c3=CornerIndex(0,1,i);
		c4=CornerIndex(1,1,i);
		return;
	}
}

void Cube::FactorFaceIndex(const int& idx,int& dir,int& offSet){
	dir  = idx>>1;
	offSet=idx &1;
}

