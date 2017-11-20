#include <iostream>
#include <fstream>
#include <time.h>
#include <assert.h>
#include "image_util.h"

double **createMatrix(int i, int j) {
    double **array = new double *[i]; // две строки в массиве
    for (int count = 0; count < i; count++)
        array[count] = new double[j];

    for (int index_i = 0; index_i < i; index_i++)
        for (int index_j = 0; index_j < j; index_j++) array[index_i][index_j] = 0;
    return array;
}

//void freeMatrix(double** M, int i, int j){
//
//}

void normalize_array(double **X, int x_h, int x_w, int block_high, int block_width, RGB **rgbInfo) {
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


int main(int argc, char *argv[]) {
    srand(time(0));
    std::string path = "Lenna.bmp";
//    BITMAPINFOHEADER fileInfoHeader = getImageInfo(name);
//    RGBQUAD **rgbInfo = getRGB(fileInfoHeader);
    //RGB **rgbInfo = getMatrixOfPixels((char *) path.c_str());
    MatrixOfImage *matrixOfImage = getMatrixOfImage((char *) path.c_str());
    RGB **rgbInfo = matrixOfImage->matrixOfPixels;
    int block_width;
    int block_high;
    int neuronAmount = 0;
    std::cout << "input block width: ";
    std::cin >> block_width;
    std::cout << std::endl << "input block high: ";
    std::cin >> block_high;
    std::cout << std::endl << "input neuron amount: ";
    std::cin >> neuronAmount;

    int N = block_width * block_high * 3;

    //---init X----------------------------------------------------------------------------
    int x_w = matrixOfImage->width / block_width;
    int x_h = matrixOfImage->height / block_high;
    double L = x_w * x_h;
    //double X[x_w * x_h][N];

    double **X = createMatrix(x_w * x_h, N);

    normalize_array(X, x_h, x_w, block_high, block_width, rgbInfo);

    //------------------------------------------------------------------------------------
    //-----------create W & W_B-----------------------------------------------------------
    double **W = createMatrix(N, neuronAmount);
    double **W_B = createMatrix(neuronAmount, N);
    double **W_B_T = createMatrix(N, neuronAmount);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < neuronAmount; j++) {
            W[i][j] = ((double) (rand()) / RAND_MAX) * 0.1;
        }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < neuronAmount; j++) {
            W_B[j][i] = W[i][j];
        }
    //-------------------------------------------------------------------------------------
    double **Y = createMatrix(x_w * x_h, neuronAmount);
    double **X_B = createMatrix(x_w * x_h, N);
    double **delta_X = createMatrix(x_w * x_h, N);
    double **tmpMatrixPN = createMatrix(neuronAmount, N);
    double **tmpMatrixNP = createMatrix(N, neuronAmount);
    double alpha;
    double alpha_b;
    double tmp;
    double *tmpArray = new double[neuronAmount];
    for (int i = 0; i < neuronAmount; i++) tmpArray[i] = 0;
    double e = 1*neuronAmount;
    double E = e + 1;
    //--------------------------------------------------------------------------------------
    int iteration = 0;
    while (E >= e) {
        std::cout << "iteration: " << iteration << " ";
        E = 0;
        for (int index = 0; index < x_w * x_h; index++) {
            //Y[index] = X[index]*W(2)
            for (int j = 0; j < neuronAmount; j++) Y[index][j] = 0;
            for (int j = 0; j < neuronAmount; j++) {
                for (int i = 0; i < N; i++) {
                    Y[index][j] += X[index][i] * W[i][j];
                }
            }

            //X_B[index] = Y[index] * W_B(3)
            for (int j = 0; j < N; j++) X_B[index][j] = 0;
            for (int j = 0; j < N; j++) {
                for (int i = 0; i < neuronAmount; i++) {
                    X_B[index][j] += Y[index][i] * W_B[i][j];
                }
            }

            //delta_X[index] = X_b[index] - X[index] (4)
            for (int j = 0; j < N; j++) {
                delta_X[index][j] = X_B[index][j] - X[index][j];
            }

            //W_B = W_B - alpha * Y_B[index] * delta_x[index](5)
            //calculate alpha_B
            tmp = N*N;
            for (int i = 0; i < neuronAmount; i++) {
                tmp += Y[index][i] * Y[index][i];
            }
            alpha_b = 1 / tmp;
            //      alpha_b = 0.0379105;
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
            tmp = N*N;
            for (int i = 0; i < N; i++) {
                tmp += X[index][i] * X[index][i];
            }
            alpha = 1 / tmp;
            //     alpha = 0.0379105;

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
            //calculate Error
            for (int i = 0; i < N; i++) E += delta_X[index][i] * delta_X[index][i];
        }
        std::cout << "Error: " << E << std::endl;
        iteration++;
    }

    double Z = (N * L) / ((N + L) * neuronAmount + 2);
    std::cout << "Z: " << Z << std::endl;

    int index = 0;
    int sub_index = 0;
    rgbInfo = new RGB *[matrixOfImage->height];
    for(int i = 0; i < matrixOfImage->height; i++) rgbInfo[i] = new RGB[matrixOfImage -> width];

    for (int block_number_i_index = 0; block_number_i_index < x_h; block_number_i_index++) {
        for (int block_number_j_index = 0; block_number_j_index < x_w; block_number_j_index++) {
            for (int i = block_number_i_index * block_high; i < block_number_i_index * block_high + block_high; i++) {
                for (int j = block_number_j_index * block_width;
                     j < block_number_j_index * block_width + block_width; j++) {
                    rgbInfo[i][j].red = 255 * (X_B[index][sub_index] + 1) / 2;
                    sub_index++;
                    rgbInfo[i][j].green = 255 * (X_B[index][sub_index] + 1) / 2;
                    sub_index++;
                    rgbInfo[i][j].blue = 255 * (X_B[index][sub_index] + 1) / 2;
                    sub_index++;
                }
            }
            sub_index = 0;
            index++;
        }
    }

    for (unsigned int i = 0; i < matrixOfImage->height; i++) {
        for (unsigned int j = 0; j < matrixOfImage->width; j++) {
            std::cout << +rgbInfo[i][j].red << " "
                      << +rgbInfo[i][j].green << " "
                      << +rgbInfo[i][j].blue << " "
                      //     << +rgbInfo[i][j].rgbReserved
                      << std::endl;
        }
        std::cout << std::endl;
    }
    return 0;
}


//for (unsigned int i = 0; i < matrixOfImage->height; i++) {
//for (unsigned int j = 0; j < matrixOfImage->width; j++) {
//std::cout << +rgbInfo[i][j].red << " "
//<< +rgbInfo[i][j].green << " "
//<< +rgbInfo[i][j].blue << " "
////     << +rgbInfo[i][j].rgbReserved
//<< std::endl;
//}
//std::cout << std::endl;
//}