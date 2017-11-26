#include <iostream>
#include <fstream>
#include <time.h>
#include "image_util.h"


double alpha;
double alpha_b;
double tmp;
int block_width;
int block_high;
int neuronAmount;
int N;
int x_w;
int x_h;
double L;
double e;
//double E;

MatrixOfImage *matrixOfImage;
RGB **rgbInfo;
double **X;
double **Y;
double **X_B;
double **delta_X;
double **tmpMatrixPN;
double **tmpMatrixNP;
double *tmpArray;
double **W;
double **W_B;
double **W_B_T;


void calculateW(int index);

double **createMatrix(int i, int j) {
    double **array = new double *[i]; // две строки в массиве
    for (int count = 0; count < i; count++)
        array[count] = new double[j];

    for (int index_i = 0; index_i < i; index_i++)
        for (int index_j = 0; index_j < j; index_j++) array[index_i][index_j] = 0;
    return array;
}

void deleteMatrix(double **array, int size) {
    for (int i = 0; i < size; i++) {
        delete array[i];
    }

    delete array;
}


void normalize_array() {
    int index = 0;
    int sub_index = 0;
    for (int block_number_i_index = 0; block_number_i_index < x_h; block_number_i_index++) {
        for (int block_number_j_index = 0; block_number_j_index < x_w; block_number_j_index++) {
            for (int i = block_number_i_index * block_high; i < block_number_i_index * block_high + block_high; i++) {
                for (int j = block_number_j_index * block_width;
                     j < block_number_j_index * block_width + block_width; j++) {
                    X[index][sub_index] = (2 * (double) rgbInfo[i][j].red / 255) - 1;
                    sub_index++;
                    X[index][sub_index] = (2 * (double) rgbInfo[i][j].green / 255) - 1;
                    sub_index++;
                    X[index][sub_index] = (2 * (double) rgbInfo[i][j].blue / 255) - 1;
                    sub_index++;
                }
            }
            sub_index = 0;
            index++;
        }
    }

}

void restore_image() {
    int index = 0;
    int sub_index = 0;
    RGB **rgbInfo = new RGB *[matrixOfImage->height];
    for (int i = 0; i < matrixOfImage->height; i++) rgbInfo[i] = new RGB[matrixOfImage->width];

    for (int block_number_i_index = 0; block_number_i_index < x_h; block_number_i_index++) {
        for (int block_number_j_index = 0; block_number_j_index < x_w; block_number_j_index++) {
            for (int i = block_number_i_index * block_high; i < block_number_i_index * block_high + block_high; i++) {
                for (int j = block_number_j_index * block_width;
                     j < block_number_j_index * block_width + block_width; j++) {
                    int tmp = 255 * (X_B[index][sub_index] + 1) / 2;
                    if (tmp > 255) tmp = 255;
                    if (tmp < 0) tmp = 0;
                    rgbInfo[i][j].red = tmp;
                    sub_index++;
                    tmp = 255 * (X_B[index][sub_index] + 1) / 2;
                    if (tmp > 255) tmp = 255;
                    if (tmp < 0) tmp = 0;
                    rgbInfo[i][j].green = tmp;
                    sub_index++;
                    tmp = 255 * (X_B[index][sub_index] + 1) / 2;
                    if (tmp > 255) tmp = 255;
                    if (tmp < 0) tmp = 0;
                    rgbInfo[i][j].blue = tmp;
                    sub_index++;
                }
            }
            sub_index = 0;
            index++;
        }
    }
    matrixOfImage->matrixOfPixels = rgbInfo;
    formImage(matrixOfImage);
}

double calculateError() {
    double E = 0;
    for (int index = 0; index < x_w * x_h; index++) {
        //  double E_q = 0;
        for (int i = 0; i < N; i++) E += delta_X[index][i] * delta_X[index][i];
        //   E = E + E_q;
    }
    std::cout << "Error: " << E << std::endl;
    return E;
}

void calculateY(int index) {
    //Y[index] = X[index]*W(2)
    for (int j = 0; j < neuronAmount; j++) Y[index][j] = 0;
    for (int j = 0; j < neuronAmount; j++) {
        for (int i = 0; i < N; i++) {
            Y[index][j] += X[index][i] * W[i][j];
        }
    }
}

void calculateXB(int index) {
    //X_B[index] = Y[index] * W_B(3)
    for (int j = 0; j < N; j++) X_B[index][j] = 0;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < neuronAmount; i++) {
            X_B[index][j] += Y[index][i] * W_B[i][j];
        }
    }
}

void calculateDeltaX(int index) {
    //delta_X[index] = X_b[index] - X[index] (4)
    for (int j = 0; j < N; j++) {
        delta_X[index][j] = X_B[index][j] - X[index][j];
    }
}


void calculateWB(int index) {
    //W_B = W_B - alpha * Y_B[index] * delta_x[index](5)
    //calculate alpha_B
    tmp = N * N;
    for (int i = 0; i < neuronAmount; i++) {
        tmp += Y[index][i] * Y[index][i];
    }
    alpha_b = 1 / tmp;
    //tmpMatrixPN = alpha_B * (Y)_T[index] * delta_X[index]
    for (int i = 0; i < neuronAmount; i++) {
        for (int j = 0; j < N; j++) {
            tmpMatrixPN[i][j] = alpha_b * Y[index][i] * delta_X[index][j];
        }
    }

    //calculate new W_b
    for (int i = 0; i < neuronAmount; i++)
        for (int j = 0; j < N; j++) {
            W_B[i][j] = W_B[i][j] - tmpMatrixPN[i][j];
        }
}

void calculateW(int index) {
    //W = W - alpha * (X[index])_T * delta_X[i] * (W_B)_T(6)
    //(W_B)_T
    for (int i = 0; i < neuronAmount; i++)
        for (int j = 0; j < N; j++) {
            W_B_T[j][i] = W_B[i][j];
        }
    // delta_X[i] * (W_B)_T
    for (int i = 0; i < neuronAmount; i++) tmpArray[i] = 0;
    for (int j = 0; j < neuronAmount; j++)
        for (int i = 0; i < N; i++) {
            tmpArray[j] += delta_X[index][i] * W_B_T[i][j];
        }
    //calculate alpha
    tmp = N * N;
    for (int i = 0; i < N; i++) {
        tmp += X[index][i] * X[index][i];
    }
    alpha = 1 / tmp;

    //(X[index])_T * delta_X[i] * (W_B)_T
    for (int i = 0; i < N; i++)
        for (int j = 0; j < neuronAmount; j++) {
            tmpMatrixNP[i][j] = alpha * X[index][i] * tmpArray[j];
        }
    //calculate new W
    for (int i = 0; i < N; i++)
        for (int j = 0; j < neuronAmount; j++) {
            W[i][j] = W[i][j] - tmpMatrixNP[i][j];
        }
}

void freeResources() {
    deleteMatrix(W, N);
    deleteMatrix(W_B, neuronAmount);
    deleteMatrix(W_B_T, N);
    deleteMatrix(X, x_w * x_h);
    deleteMatrix(Y, x_w * x_h);
    deleteMatrix(X_B, x_w * x_h);
    deleteMatrix(delta_X, x_w * x_h);
    deleteMatrix(tmpMatrixPN, neuronAmount);
    deleteMatrix(tmpMatrixNP, N);
    for (int i = 0; i < matrixOfImage->height; i++) delete rgbInfo[i];
    delete rgbInfo;
    delete matrixOfImage;
    delete tmpArray;
}


int main(int argc, char *argv[]) {
    srand(time(0));
    std::string path = "test2.bmp";
    double MAX_VALUE = 1;

    matrixOfImage = getMatrixOfImage((char *) path.c_str());
    rgbInfo = matrixOfImage->matrixOfPixels;
    block_width;
    block_high;
    neuronAmount = 0;
    std::cout << "input block width: ";
    std::cin >> block_width;
    std::cout << std::endl << "input block high: ";
    std::cin >> block_high;
    std::cout << std::endl << "input neuron amount: ";
    std::cin >> neuronAmount;

    N = block_width * block_high * 3;

    //---init X----------------------------------------------------------------------------
    x_w = matrixOfImage->width / block_width;
    x_h = matrixOfImage->height / block_high;
    L = x_w * x_h;

    X = createMatrix(x_w * x_h, N);

    normalize_array();

    //-----------create W & W_B-----------------------------------------------------------
    W = createMatrix(N, neuronAmount);
    W_B = createMatrix(neuronAmount, N);
    W_B_T = createMatrix(N, neuronAmount);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < neuronAmount; j++) {
            W[i][j] = ((double) (rand()) / RAND_MAX) * MAX_VALUE;
        }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < neuronAmount; j++) {
            W_B[j][i] = W[i][j];
        }
    //-------------------------------------------------------------------------------------
    Y = createMatrix(x_w * x_h, neuronAmount);
    X_B = createMatrix(x_w * x_h, N);
    delta_X = createMatrix(x_w * x_h, N);
    tmpMatrixPN = createMatrix(neuronAmount, N);
    tmpMatrixNP = createMatrix(N, neuronAmount);
    alpha;
    alpha_b;
    tmp;
    tmpArray = new double[neuronAmount];


    for (int i = 0; i < neuronAmount; i++) tmpArray[i] = 0;
    std::cout << std::endl << "input max error: ";
    std::cin >> e;

    //--------------------------------------------------------------------------------------
    int iteration = 0;
    do {
        iteration++;
        std::cout << "iteration: " << iteration << " ";
        for (int index = 0; index < x_w * x_h; index++) {
            calculateY(index);
            calculateXB(index);
            calculateDeltaX(index);
            calculateWB(index);
            calculateW(index);
        }
    } while (calculateError() > e);

    double Z = (N * L) / ((N + L) * neuronAmount + 2);
    std::cout << "Z: " << Z << std::endl;


    restore_image();
    freeResources();
    return 0;
}


