/*
 * MatMul.h
 *
 *  Created on: Oct 28, 2022
 *      Author: cagnur
 */

#ifndef MATMUL_H_
#define MATMUL_H_

class MatMul {

public:
	int bsize=100;
	float **C;
	int col;
	int row;
	int point;
	MatMul(int dim);
	virtual ~MatMul();
	/* Experiment-1*/
	void reset(void);
	void ijk(float **firstMat, float **secondMat);
	void jik(float **firstMat, float **secondMat);
	void jki(float **firstMat, float **secondMat);
	void kji(float **firstMat, float **secondMat);
	void kij(float **firstMat, float **secondMat);
	void ikj(float **firstMat, float **secondMat);
	void bijk(float **firstMat, float **secondMat);
	void possibleIndexingTimeAnalysis(float **firstMat, float **secondMat);
	/* Experiment-2*/
	void unrollingTwo(float **firstMat, float **secondMat);
	void unrollingFour(float **firstMat, float **secondMat);
	void unrollingEight(float **firstMat, float **secondMat);
	void unrollingSixteen(float **firstMat, float **secondMat);
	void onlyUnrollingTimeAnalysis(float **firstMat, float **secondMat);
	/*Experiment-3*/
	void unrollingFourWithFission(float **firstMat, float **secondMat);
	void unrollingEightWithFission(float **firstMat, float **secondMat);
	void unrollingSixteenWithFission(float **firstMat, float **secondMat);
	void unrollingWithFissionTimeAnalysis(float **firstMat, float **secondMat);
	/* Experiment-4*/
	void changeBlockSize(int blockSize);
	void blockSizeTimeAnalysis(float **firstMat, float **secondMat, int start, int stop, int stepSize);
	/* Experiment-5*/
	void printMat(void);
	void bikj(float **firstMat, float **secondMat);
	void bijk_updated(float **firstMat, float **secondMat);
	void loopBlockingTimeAnalysis(float **firstMat, float **secondMat);
	/* Experiment-6*/
	void bikjUnrolledFour(float **firstMat, float **secondMat);
	void bikjUnrolledEight(float **firstMat, float **secondMat);
	void bikjUnrolledSixteen(float **firstMat, float **secondMat);
	void bikjUnrolledTimeAnalysis(float **firstMat, float **secondMat);

};

#endif /* MATMUL_H_ */
