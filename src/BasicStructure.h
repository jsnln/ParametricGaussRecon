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

#ifndef	BASIC_STRUCTURE_HEADER
#define BASIC_STRUCTURE_HEADER

#define EPSILON 1e-6
#define NON_UNIT_NORM 0
#define NEW_SAMPLES_PER_NODE 1

//#include <limits>

struct Point{
	float x;
	float y;
	float z;
	Point(){
		x = 0; y = 0; z = 0;
	}
	Point(float xx, float yy, float zz){
		x = xx;
		y = yy;
		z = zz;
	}
};

struct WeightedPoint{
	float x;
	float y;
	float z;
	float weight;
	WeightedPoint(){
		x = y = z = 0;
	}
	WeightedPoint(float xx, float yy, float zz){
		x = xx; y = yy; z = zz;
	}
};

struct BoundingBox{
	float xscale;
	float yscale;
	float zscale;
	float blx;		// bottom-left x
	float bly;		// bottom-left y
	float blz;		// bottom-left z
};


#endif
