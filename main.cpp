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

void normalize_array(double **X, int x_h, int x_w, int block_high, int block_width, RGBQUAD **rgbInfo) {
    int index = 0;
    int sub_index = 0;
    for (int block_number_i_index = 0; block_number_i_index < x_h; block_number_i_index++) {
        for (int block_number_j_index = 0; block_number_j_index < x_w; block_number_j_index++) {
            for (int i = block_number_i_index * block_high; i < block_number_i_index * block_high + block_high; i++) {
                for (int j = block_number_j_index * block_width;
                     j < block_number_j_index * block_width + block_width; j++) {
                    X[index][sub_index] = (2 * (double) rgbInfo[i][j].rgbRed / 255) - 1;
                    sub_index++;
                    X[index][sub_index] = (2 * (double) rgbInfo[i][j].rgbGreen / 255) - 1;
                    sub_index++;
                    X[index][sub_index] = (2 * (double) rgbInfo[i][j].rgbBlue / 255) - 1;
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
    char *name = "Lenna.bmp";
    BITMAPINFOHEADER fileInfoHeader = getImageInfo(name);
    RGBQUAD **rgbInfo = getRGB(fileInfoHeader);
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
    int x_w = fileInfoHeader.biWidth / block_width;
    int x_h = fileInfoHeader.biHeight / block_high;
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
    double e = 0.3*neuronAmount;
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

//            std::cout << std::endl << "-----------------------------------" << std::endl << index << std::endl;
//            for (int i = 0; i < neuronAmount; i++)
//                for (int j = 0; j < N; j++) {
//                    std::cout << tmpMatrixPN[index][j] << " ";
//                }

//            //alpha_b * tmpMatrixPN
//        for(int i = 0; i<neuronAmount; i++)
//            for(int j = 0; j<N; j++){
//                tmpMatrixPN[i][j] = tmpMatrixPN[i][j]*alpha_b;
//            }

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
    std::cout << std::endl << "-----------------------------" << std::endl;
    for (int i = 0; i < N; i++) std::cout << X[0][i] << " ";
    std::cout << std::endl << "-----------------------------" << std::endl;
    for (int i = 0; i < N; i++) std::cout << X_B[0][i] << " ";
    return 0;
}

// вывод
// std::cout<<fileInfoHeader.biHeight<<" "<<fileInfoHeader.biWidth;
//    for (unsigned int i = 0; i < fileInfoHeader.biHeight; i++) {
//        for (unsigned int j = 0; j < fileInfoHeader.biWidth; j++) {
//            std::cout << +rgbInfo[i][j].rgbRed << " "
//                      << +rgbInfo[i][j].rgbGreen << " "
//                      << +rgbInfo[i][j].rgbBlue << " "
//                      //     << +rgbInfo[i][j].rgbReserved
//                      << std::endl;
//        }
//        std::cout << std::endl;
//    }