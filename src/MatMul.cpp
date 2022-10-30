/*
 * MatMul.cpp
 *
 *  Created on: Oct 28, 2022
 *      Author: cagnur
 */
#include <iostream>
#include <sys/time.h>
#include "MatMul.h"
using namespace std;


MatMul::MatMul(int dim) {
	row = dim;
	col = dim;
	point = dim;
	C=new float *[row];
	for(int i=0;i<row;i++)
	{
	  C[i]=new float[col];
	}
}

MatMul::~MatMul() {
    for(int i = 0; i < row; ++i)
        delete [] C[i];

    delete [] C;
}

void MatMul::reset(void){
	int i,j;
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			C[i][j] = 0.0;
		}
	}
}

void MatMul::ijk(float **firstMat, float **secondMat){
	for (int i = 0; i <row; i++)
		for(int j = 0; j< col; j++)
		{
			float sum = 0;
			for (int k = 0; k< point; k++ )
			{
				sum += firstMat[i][k]*secondMat[k][j];
			}
			C[i][j]= sum;
		}
}

void MatMul::jik(float **firstMat, float **secondMat ){
	for (int j = 0; j <row; j++)
		for(int i = 0; i< col; i++)
		{
			float sum = 0;
			for (int k = 0; k< point; k++ )
			{
				sum += firstMat[i][k]*secondMat[k][j];
			}
			C[i][j]= sum;
		}
}

void MatMul::jki(float **firstMat, float **secondMat ){
	float r;
	for (int j = 0; j <col; j++)
		for(int k = 0; k< point; k++)
		{
			r = secondMat[k][j];
			for (int i = 0; i< col; i++ )
			{
				C[i][j] += firstMat[i][k]*r;
			}
		}
}

void MatMul::kji(float **firstMat, float **secondMat ){
	float r;
	for (int k = 0; k <col; k++)
		for(int j = 0; j< point; j++)
		{
			r = secondMat[k][j];
			for (int i = 0; i< col; i++ )
			{
				C[i][j] += firstMat[i][k]*r;
			}
		}
}

void MatMul::kij(float **firstMat, float **secondMat ){
	float r;
	for (int k = 0; k<point; k++)
		for(int i = 0; i< row; i++)
		{
			r = firstMat[i][k];
			for (int j = 0; j< row; j++ )
			{
				C[i][j] += r * secondMat[k][j];
			}
		}
}

void MatMul::ikj(float **firstMat, float **secondMat ){
	for (int i = 0; i<point; i++)
		for(int k = 0; k< row; k++)
		{
			float r = firstMat[i][k];
			for (int j = 0; j< row; j++ )
			{
				C[i][j] += r * secondMat[k][j];
			}
		}
}

void MatMul::bijk(float **firstMat, float **secondMat){
	int i, j, k, kk, jj;
	float sum;
	int n = col; /* or row*/
	for (jj=0; jj<n; jj+=bsize) {
		for (i=0; i<n; i++)
			for (j=jj; j < min(jj+bsize,n); j++)
				C[i][j] = 0.0;
		for (kk=0; kk<n; kk+=bsize) {
			for (i=0; i<n; i++) {
				for (j=jj; j < min(jj+bsize,n); j++) {
					sum = 0.0;
					for (k=kk; k < min(kk+bsize,n); k++) {
						sum += firstMat[i][k] * secondMat[k][j];
					}
					C[i][j] += sum;
				}
			}
		}
	}
}

void MatMul::possibleIndexingTimeAnalysis(float **firstMat, float **secondMat){
	this->reset();

    struct timeval t;
    double time1, time2, wct; /* unit is second */

    cout<< "For bijk multiplication wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->bijk(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "For ijk multiplication wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->ijk(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "For jik multiplication wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->jik(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "For jki multiplication wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->jki(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "For kji multiplication wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->kji(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "For kij multiplication wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->kij(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "For ikj multiplication wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->ikj(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;
}

void MatMul::unrollingTwo(float **firstMat, float **secondMat){
	int unroll = 2;
	for (int i = 0; i<point; i++)
		for(int k = 0; k< row; k++)
		{
			float r = firstMat[i][k];
			for (int j = 0; j< row;j+=unroll)
			{
				C[i][j] += r * secondMat[k][j];
				C[i][j+1] += r * secondMat[k][j+1];
			}
		}
}

void MatMul::unrollingFour(float **firstMat, float **secondMat){
	int unroll = 4;
	for (int i = 0; i<point; i++)
		for(int k = 0; k< row; k++)
		{
			float r = firstMat[i][k];
			for (int j = 0; j< row;j+=unroll)
			{
				C[i][j] += r * secondMat[k][j];
				C[i][j+1] += r * secondMat[k][j+1];
				C[i][j+2] += r * secondMat[k][j+2];
				C[i][j+3] += r * secondMat[k][j+3];
			}
		}
}

void MatMul::unrollingEight(float **firstMat, float **secondMat){
	int unroll = 8;
	for (int i = 0; i<point; i++)
		for(int k = 0; k< row; k++)
		{
			float r = firstMat[i][k];
			for (int j = 0; j< row; j+=unroll)
			{
				C[i][j] += r * secondMat[k][j];
				C[i][j+1] += r * secondMat[k][j+1];
				C[i][j+2] += r * secondMat[k][j+2];
				C[i][j+3] += r * secondMat[k][j+3];
				C[i][j+4] += r * secondMat[k][j+4];
				C[i][j+5] += r * secondMat[k][j+5];
				C[i][j+6] += r * secondMat[k][j+6];
				C[i][j+7] += r * secondMat[k][j+7];
			}
		}
}

void MatMul::unrollingSixteen(float **firstMat, float **secondMat){
	int unroll = 16;
	for (int i = 0; i<point; i++)
		for(int k = 0; k< row; k++)
		{
			float r = firstMat[i][k];
			for (int j = 0; j< row; j+=unroll)
			{
				C[i][j] += r * secondMat[k][j];
				C[i][j+1] += r * secondMat[k][j+1];
				C[i][j+2] += r * secondMat[k][j+2];
				C[i][j+3] += r * secondMat[k][j+3];
				C[i][j+4] += r * secondMat[k][j+4];
				C[i][j+5] += r * secondMat[k][j+5];
				C[i][j+6] += r * secondMat[k][j+6];
				C[i][j+7] += r * secondMat[k][j+7];
				C[i][j+8] += r * secondMat[k][j+8];
				C[i][j+9] += r * secondMat[k][j+9];
				C[i][j+10] += r * secondMat[k][j+10];
				C[i][j+11] += r * secondMat[k][j+11];
				C[i][j+12] += r * secondMat[k][j+12];
				C[i][j+13] += r * secondMat[k][j+13];
				C[i][j+14] += r * secondMat[k][j+14];
				C[i][j+15] += r * secondMat[k][j+15];
			}
		}
}

void MatMul::onlyUnrollingTimeAnalysis(float **firstMat, float **secondMat){
	this->reset();

	cout << "Unrolling use only ikj indexing"<<endl;
    struct timeval t;
    double time1, time2, wct; /* unit is second */

    cout<< "Loop Unrolling-2, wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->unrollingTwo(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "Loop Unrolling-4, wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->unrollingFour(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "Loop Unrolling-8, wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->unrollingEight(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "Loop Unrolling-16, wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->unrollingSixteen(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;
}

void MatMul::unrollingFourWithFission(float **firstMat, float **secondMat){
	int unroll = 4;
	for (int i = 0; i<point; i++)
		for(int k = 0; k< row; k++)
		{
			float r = firstMat[i][k];
			for (int j = 0; j< row;j+=unroll)
			{
				C[i][j] += r * secondMat[k][j];
				C[i][j+1] += r * secondMat[k][j+1];
			}
			for(int j = 0; j< row;j+=unroll)
			{
				C[i][j+2] += r * secondMat[k][j+2];
				C[i][j+3] += r * secondMat[k][j+3];
			}
		}
}

void MatMul::unrollingEightWithFission(float **firstMat, float **secondMat){
	int unroll = 8;
	for (int i = 0; i<point; i++)
		for(int k = 0; k< row; k++)
		{
			float r = firstMat[i][k];
			for (int j = 0; j< row; j+=unroll)
			{
				C[i][j] += r * secondMat[k][j];
				C[i][j+1] += r * secondMat[k][j+1];
				C[i][j+2] += r * secondMat[k][j+2];
				C[i][j+3] += r * secondMat[k][j+3];
			}
			for (int j = 0; j< row; j+=unroll)
			{
				C[i][j+4] += r * secondMat[k][j+4];
				C[i][j+5] += r * secondMat[k][j+5];
				C[i][j+6] += r * secondMat[k][j+6];
				C[i][j+7] += r * secondMat[k][j+7];
			}
		}
}

void MatMul::unrollingSixteenWithFission(float **firstMat, float **secondMat){
	int unroll = 16;
	for (int i = 0; i<point; i++)
		for(int k = 0; k< row; k++)
		{
			float r = firstMat[i][k];
			for (int j = 0; j< row; j+=unroll)
			{
				C[i][j] += r * secondMat[k][j];
				C[i][j+1] += r * secondMat[k][j+1];
				C[i][j+2] += r * secondMat[k][j+2];
				C[i][j+3] += r * secondMat[k][j+3];
				C[i][j+4] += r * secondMat[k][j+4];
				C[i][j+5] += r * secondMat[k][j+5];
				C[i][j+6] += r * secondMat[k][j+6];
				C[i][j+7] += r * secondMat[k][j+7];
			}
			for (int j = 0; j< row; j+=unroll)
			{
				C[i][j+8] += r * secondMat[k][j+8];
				C[i][j+9] += r * secondMat[k][j+9];
				C[i][j+10] += r * secondMat[k][j+10];
				C[i][j+11] += r * secondMat[k][j+11];
				C[i][j+12] += r * secondMat[k][j+12];
				C[i][j+13] += r * secondMat[k][j+13];
				C[i][j+14] += r * secondMat[k][j+14];
				C[i][j+15] += r * secondMat[k][j+15];
			}
		}
}

void MatMul::unrollingWithFissionTimeAnalysis(float **firstMat, float **secondMat){
	this->reset();

	cout << "Unrolling use only ikj indexing"<<endl;
    struct timeval t;
    double time1, time2, wct; /* unit is second */

    cout<< "Loop Unrolling-4 with fission, wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->unrollingFourWithFission(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "Loop Unrolling-8 with fission, wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->unrollingEightWithFission(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "Loop Unrolling-16 with fission, wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->unrollingSixteenWithFission(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;
}

void MatMul::printMat(void)
{
	int i,j;
	cout << "Matrix is the following:"<< endl;
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			cout << C[i][j] << " ";
		}
		cout<<' '<<endl;
	}
}

void MatMul::changeBlockSize(int blockSize){
	bsize =blockSize;
	cout << "New block size is "<<bsize<<endl;
}

void MatMul::blockSizeTimeAnalysis(float **firstMat, float **secondMat, int start, int stop, int stepSize){
    struct timeval t;
    double time1, time2, wct; /* unit is second */
	for(int b= start; b<stop; b*=stepSize){
		this->reset();
		this->changeBlockSize(b);

	    cout<<"Wall time for bijk: ";
	    gettimeofday(&t, NULL);
	    time1= t.tv_sec + 1.0e-6*t.tv_usec;
	    this->bijk(firstMat, secondMat);
	    gettimeofday(&t, NULL);
	    time2= t.tv_sec + 1.0e-6*t.tv_usec;
	    wct = time2-time1;
	    cout << wct<< endl;
	}
}

void MatMul::bikj(float **firstMat, float **secondMat){
	int i, j, k, kk, ii;
	float r;
	int n = col; /* or row*/
	for (ii=0; ii<row; ii+=bsize) {
		for (kk=0; kk<row; kk+=bsize) {
			for (i=ii; i < min(ii+bsize,n); i++) {
				for (k=kk; k < min(kk+bsize,n); k++) {
					r = firstMat[i][k];
					for(j = 0; j <n; j++)
						C[i][j] += r * secondMat[k][j];
					}
				}
			}
		}
}

void MatMul::bijk_updated(float **firstMat, float **secondMat){
	int i, j, k, kk, jj;
	float sum;
	int n = col; /* or row*/
	for (jj=0; jj<n; jj+=bsize) {
		for (kk=0; kk<n; kk+=bsize) {
			for (i=0; i<n; i++) {
				for (j=jj; j < min(jj+bsize,n); j++) {
					sum = 0.0;
					for (k=kk; k < min(kk+bsize,n); k++) {
						sum += firstMat[i][k] * secondMat[k][j];
					}
					C[i][j] += sum;
				}
			}
		}
	}
}

void MatMul::loopBlockingTimeAnalysis(float **firstMat, float **secondMat){
	this->reset();

    struct timeval t;
    double time1, time2, wct; /* unit is second */

    cout<< "Wall time for ijk is: ";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->ijk(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "Wall time for bijk_updated is: ";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->bijk_updated(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "Wall time for ikj is: ";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->ikj(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "Wall time for bikj is: ";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->bikj(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;
}

void MatMul::bikjUnrolledFour(float **firstMat, float **secondMat){
	int i, j, k, kk, ii;
	float r;
	int unroll = 4;
	int n = col; /* or row*/
	for (ii=0; ii<row; ii+=bsize) {
		for (kk=0; kk<row; kk+=bsize) {
			for (i=ii; i < min(ii+bsize,n); i++) {
				for (k=kk; k < min(kk+bsize,n); k++) {
					r = firstMat[i][k];
					for(j = 0; j <n; j+=unroll){
						C[i][j] += r * secondMat[k][j];
						C[i][j+1] += r * secondMat[k][j+1];
						C[i][j+2] += r * secondMat[k][j+2];
						C[i][j+3] += r * secondMat[k][j+3];
					}
				}
			}
		}
	}
}

void MatMul::bikjUnrolledEight(float **firstMat, float **secondMat){
	int i, j, k, kk, ii;
	float r;
	int unroll = 8;
	int n = col; /* or row*/
	for (ii=0; ii<row; ii+=bsize) {
		for (kk=0; kk<row; kk+=bsize) {
			for (i=ii; i < min(ii+bsize,n); i++) {
				for (k=kk; k < min(kk+bsize,n); k++) {
					r = firstMat[i][k];
					for(j = 0; j <n; j+=unroll){
						C[i][j] += r * secondMat[k][j];
						C[i][j+1] += r * secondMat[k][j+1];
						C[i][j+2] += r * secondMat[k][j+2];
						C[i][j+3] += r * secondMat[k][j+3];
						C[i][j+4] += r * secondMat[k][j+4];
						C[i][j+5] += r * secondMat[k][j+5];
						C[i][j+6] += r * secondMat[k][j+6];
						C[i][j+7] += r * secondMat[k][j+7];
					}
				}
			}
		}
	}
}

void MatMul::bikjUnrolledSixteen(float **firstMat, float **secondMat){
	int i, j, k, kk, ii;
	float r;
	int unroll = 16;
	int n = col; /* or row*/
	for (ii=0; ii<row; ii+=bsize) {
		for (kk=0; kk<row; kk+=bsize) {
			for (i=ii; i < min(ii+bsize,n); i++) {
				for (k=kk; k < min(kk+bsize,n); k++) {
					r = firstMat[i][k];
					for(j = 0; j <n; j+=unroll){
						C[i][j] += r * secondMat[k][j];
						C[i][j+1] += r * secondMat[k][j+1];
						C[i][j+2] += r * secondMat[k][j+2];
						C[i][j+3] += r * secondMat[k][j+3];
						C[i][j+4] += r * secondMat[k][j+4];
						C[i][j+5] += r * secondMat[k][j+5];
						C[i][j+6] += r * secondMat[k][j+6];
						C[i][j+7] += r * secondMat[k][j+7];
						C[i][j+8] += r * secondMat[k][j+8];
						C[i][j+9] += r * secondMat[k][j+9];
						C[i][j+10] += r * secondMat[k][j+10];
						C[i][j+11] += r * secondMat[k][j+11];
						C[i][j+12] += r * secondMat[k][j+12];
						C[i][j+13] += r * secondMat[k][j+13];
						C[i][j+14] += r * secondMat[k][j+14];
						C[i][j+15] += r * secondMat[k][j+15];
					}
				}
			}
		}
	}
}

void MatMul::bikjUnrolledTimeAnalysis(float **firstMat, float **secondMat){
	this->reset();

    struct timeval t;
    double time1, time2, wct; /* unit is second */

    cout<< "bikj + Loop Unrolling-4 , Wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->bikjUnrolledFour(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "bikj + Loop Unrolling-8 , Wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->bikjUnrolledEight(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;

    this->reset();

    cout<< "bikj + Loop Unrolling-16 , Wall time is:";
    gettimeofday(&t, NULL);
    time1= t.tv_sec + 1.0e-6*t.tv_usec;
    this->bikjUnrolledSixteen(firstMat, secondMat);
    gettimeofday(&t, NULL);
    time2= t.tv_sec + 1.0e-6*t.tv_usec;
    wct = time2-time1;
    cout << wct<< endl;
}



