#include <malloc.h>
#include "image_util.h"

BmpImageInfo bmpImageInfo;

MatrixOfImage* getMatrixOfImage(char *imagePath) {
    MatrixOfImage *matrixOfImage = (MatrixOfImage *) malloc(sizeof(MatrixOfImage));
    matrixOfImage->matrixOfPixels = getMatrixOfPixels(imagePath);
    matrixOfImage->height = bmpImageInfo.height;
    matrixOfImage->width = bmpImageInfo.width;
    return matrixOfImage;
}

RGB **getMatrixOfPixels(char *imagePath) {
    FILE *readImage = fopen(imagePath, "rb");
    readBmpHeader(readImage);
    bmpImageInfo = readBmpImageInfo(readImage);
    readBmpImagePalette(readImage, bmpImageInfo.colorsCount);
    RGB **matrixOfPixels = readMatrixOfPixels(readImage, bmpImageInfo.height, bmpImageInfo.width);
    fclose(readImage);
    return matrixOfPixels;
}

BmpHeader readBmpHeader(FILE *readImage) {
    BmpHeader bmpHeader;
    if (readImage != NULL) {
        fread(&bmpHeader, 1, sizeof(BmpHeader), readImage);
    }
    return bmpHeader;
}

BmpImageInfo readBmpImageInfo(FILE *readImage) {
    BmpImageInfo bmpImageInfo;
    if (readImage != NULL) {
        fread(&bmpImageInfo, 1, sizeof(BmpImageInfo), readImage);
    }
    return bmpImageInfo;
}

RGB* readBmpImagePalette(FILE *readImage, unsigned int colorsCount) {
    RGB* palette = (RGB *) malloc(sizeof(RGB) * colorsCount);
    fread(palette, sizeof(RGB), colorsCount, readImage);
    return palette;
}

RGB** createMatrixOfPixels(unsigned int imageHeight, unsigned int imageWidth) {
    RGB **matrixOfPixels = (RGB **) malloc(sizeof(RGB*) * imageHeight);
    for (int heightIndex = 0; heightIndex < imageHeight; heightIndex++) {
        matrixOfPixels[heightIndex] = (RGB *) malloc(sizeof(RGB) * imageWidth);
    }
    return matrixOfPixels;
}

RGB** readMatrixOfPixels(FILE *readImage, unsigned int imageHeight, unsigned int imageWidth) {
    RGB **matrixOfPixels = createMatrixOfPixels(imageWidth, imageHeight);
    RGB *pixel = (RGB *) malloc(sizeof(RGB*));
    for (int indexHeight = 0; indexHeight < imageHeight; indexHeight++) {
        for (int indexWidth = 0; indexWidth < imageWidth; indexWidth++) {
            fread(pixel, 1, sizeof(RGB), readImage);
            matrixOfPixels[indexHeight][indexWidth].red = pixel->red;
            matrixOfPixels[indexHeight][indexWidth].green = pixel->green;
            matrixOfPixels[indexHeight][indexWidth].blue = pixel->blue;
        }
    }
    return matrixOfPixels;
}
