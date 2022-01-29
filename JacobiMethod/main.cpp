#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

ifstream fileOut;
ofstream file;
static double *matrix;
static int numOfVariables;
static int numOfEquations;
static int numOfIter = 0;

void getMatrix(int numOfVar, int numOfEq, double *matrix/*out*/) {
	for (int i = 0; i < numOfEq; i++) {
		for (int j = 0; j < numOfVar; j++) {
			fileOut >> matrix[i * numOfVar + j];
		}
	}
}

void getMartixSize(int *numOfVar/*out*/, int *numOfEq/*out*/) {
	fileOut >> *numOfEq >> *numOfVar;
	return;
}

void PrintMatrix(double *matrix, int numOfVariables, int numOfEquations) {
	file.precision(5);

	for (int i = 0; i < numOfEquations; i++) {
		for (int j = 0; j < numOfVariables; j++) {
			file << matrix[i * numOfVariables + j] << "     ";
		}
		file << endl;
	}
	file << "--------------------------------------";
}

void MultMatrixAB(double *matrixA, double *matrixB) {
	

	for (int i = 0; i < numOfEquations; i++) {
		for (int j = 0; j < numOfVariables; j++) {
			double sum = 0;
			for (int k = 0; k < numOfVariables; k++) {
				sum += matrixA[i * numOfVariables + k] * matrixB[k * numOfVariables + j];
			}
			matrixA[i * numOfVariables + j] = sum;
		}
	}
}

double normMatrixWithoutDiagElements(double *A, int numOfVar, int numOfEq) {
	double sum = 0;
	for (int i = 0; i < numOfEq; i++) {
		for (int j = 0; j < numOfVar; j++) {
			if (i != j) {
				double tmpValue = 0;
				tmpValue = A[i * numOfVar + j];
				tmpValue = powl(tmpValue, 2);
				sum += tmpValue;
			}
		}
	}
	return sqrt(sum);
}

void CopyMatrix(double *to, double *from, int numOfRows, int numOfColumns, int fromColumn, int fromRow) {
	for (int i = fromRow; i < numOfRows; i++) {
		for (int j = fromColumn; j < numOfColumns; j++) {
			to[i * numOfColumns + j] = from[i * numOfColumns + j];
		}
	}
}

void CalculateRotationMatrixParameters(double *matrix, double *c /*cos(a)*/, double *s /*sin(a)*/, int i, int j) {
	double p = 2 * matrix[i * numOfVariables + j];
	double q = matrix[i *numOfVariables + i] - matrix[j *numOfVariables + j];
	double d = sqrt(powl(p, 2) + powl(q, 2));
	int signPQ = (int)((p * q) / abs(p * q));

	if (q != 0) {
		double r = abs(q) / (2 * d);
		*c = sqrt(0.5 + r);
		*s = sqrt(0.5 - r) * signPQ;
	}
	else {
		*c = *s = 1 / sqrt(2);
	}
	return;
}

void InitByZero(double *data, int size) {
	for (int i = 0; i < size; i++) {
		data[i] = 0;
	}
}

double *GetIdentityMatrix(int numOfColumns, int numOfRows) {
	double *identityMatrix = new double[numOfColumns * numOfRows];


	InitByZero(identityMatrix, numOfColumns * numOfRows);
	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			if (i == j) {
				identityMatrix[i * numOfColumns + j] = 1;
			}
		}
	}
	return identityMatrix;
}

void ChooseSpecialElementFromSymmetricMatrix(double *matrix, int *iRow, int *jColumn) {
	double max = 0;

	for (int i = 1; i < numOfEquations; i++) {
		for (int j = 0; j < i; j++) {
			if (matrix[i * numOfVariables + j] > max) {
				*iRow = i;
				*jColumn = j;
				max = matrix[i * numOfVariables + j];
			}
		}
	}
	return;
}

void InitRotationMatrix(double *identityMatrix, int i, int j, double c, double s) {
	identityMatrix[i * numOfVariables + i] = identityMatrix[j * numOfVariables + j] = c;
	identityMatrix[i * numOfVariables + j] = -s;
	identityMatrix[j * numOfVariables + i] = s;
}

void JacobiMethod(double error, double **diagMatrix /* out */, double **eigenVectors /* out */) {
	// A = G^(-1) * B * G
	double *A = new double [numOfVariables * numOfEquations];
	double *B = new double [numOfVariables * numOfEquations];
	double *MatrixOfEigenvectors = GetIdentityMatrix(numOfVariables, numOfEquations);
	double *TmpMatrixOfEigenvectors = GetIdentityMatrix(numOfVariables, numOfEquations);
	double normB = normMatrixWithoutDiagElements(B, numOfVariables, numOfEquations);
	double prevNorm = normB + 1;
	bool isFirstIteration = true;

	CopyMatrix(A, matrix, numOfEquations, numOfVariables, 0, 0);
	CopyMatrix(B, A, numOfEquations, numOfVariables, 0, 0);

	while (normB > error) {
		if (!isFirstIteration) {
			CopyMatrix(A, B, numOfEquations, numOfVariables, 0, 0);
			prevNorm = normB;
		}

		int i, j;
		ChooseSpecialElementFromSymmetricMatrix(A, &i, &j);

		double c, s;
		CalculateRotationMatrixParameters(A, &c, &s, i, j);

		if (isFirstIteration) {
			InitRotationMatrix(MatrixOfEigenvectors, i, j, c, s);
			isFirstIteration = false;
		}
		else if (!isFirstIteration) {
			InitRotationMatrix(TmpMatrixOfEigenvectors, i, j, c, s);
			MultMatrixAB(MatrixOfEigenvectors, TmpMatrixOfEigenvectors);
		}

		B[i * numOfVariables + i] = powl(c, 2) * A[i * numOfVariables + i] + powl(s, 2) * A[j * numOfVariables + j] + 2 * c * s * A[i * numOfVariables + j];
		B[j * numOfVariables + j] = powl(s, 2) * A[i * numOfVariables + i] + powl(c, 2) * A[j * numOfVariables + j] - 2 * c * s * A[i * numOfVariables + j];

		B[i * numOfVariables + j] = B[j * numOfVariables + i] = (powl(c, 2) - powl(s, 2)) * A[i * numOfVariables + j] + c * s * (A[j * numOfVariables + j] - A[i * numOfVariables + i]);
		B[i * numOfVariables + j] = B[j * numOfVariables + i] = 0;

		for (int m = 0; m < numOfEquations; m++) {
			if (m != i && m != j) {
				B[i * numOfVariables + m] = B[m * numOfVariables + i] = c * A[m * numOfVariables + i] + s * A[m * numOfVariables + j];
				B[j * numOfVariables + m] = B[m * numOfVariables + j] = - s * A[m * numOfVariables + i] + c * A[m * numOfVariables + j];
			}
		}
		numOfIter++;
		normB = normMatrixWithoutDiagElements(B, numOfVariables, numOfEquations);
		if (!isFirstIteration && prevNorm == normB) {
			break;
		}
	}
	*diagMatrix = B;
	*eigenVectors = MatrixOfEigenvectors;

	delete[] A;
	delete[] TmpMatrixOfEigenvectors;
}

int main(void) {
	fileOut.open("matr1.txt");

	getMartixSize(&numOfVariables, &numOfEquations);

	matrix = new double[numOfVariables * numOfEquations];
	getMatrix(numOfVariables, numOfEquations, matrix);

	fileOut.close();

	//--------------------------------------------------
	double error = 0.0001;
	double *diagMatrix = NULL, *eigenVectors = NULL;
	JacobiMethod(error, &diagMatrix, &eigenVectors);


	file.open("result.txt");

	PrintMatrix(diagMatrix, numOfVariables, numOfEquations);
	PrintMatrix(eigenVectors, numOfVariables, numOfEquations);

	file.close();
	return 0;
}