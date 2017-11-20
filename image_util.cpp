#include <malloc.h>
#include "image_util.h"

BmpHeader bmpHeader;
BmpImageInfo bmpImageInfo;
RGB *palette;

MatrixOfImage* getMatrixOfImage(char *imagePath) {
    MatrixOfImage *matrixOfImage = (MatrixOfImage *) malloc(sizeof(MatrixOfImage));
    matrixOfImage->matrixOfPixels = getMatrixOfPixels(imagePath);
    matrixOfImage->height = bmpImageInfo.height;
    matrixOfImage->width = bmpImageInfo.width;
    return matrixOfImage;
}

RGB **getMatrixOfPixels(char *imagePath) {
    FILE *readImage = fopen(imagePath, "rb");
    bmpHeader = readBmpHeader(readImage);
    bmpImageInfo = readBmpImageInfo(readImage);
    fseek(readImage, bmpHeader.bfOffBits, SEEK_SET);
//    readBmpImagePalette(readImage, bmpImageInfo.colorsCount);
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

void formImage(MatrixOfImage *matrixOfImage) {
    FILE *file = fopen("result.bmp", "wb");

    fwrite(&bmpHeader,1,sizeof(bmpHeader),file);
    fwrite(&bmpImageInfo, 1, sizeof(bmpImageInfo), file);
    fseek(file, bmpHeader.bfOffBits, SEEK_SET);

//    fwrite(&palette, 1, sizeof(palette), file);

    for (int heightIndex = 0; heightIndex < matrixOfImage->height; heightIndex++) {
        for (int widthIndex = 0; widthIndex < matrixOfImage->width; widthIndex++) {

            unsigned char red = matrixOfImage->matrixOfPixels[heightIndex][widthIndex].red;
            unsigned char green = matrixOfImage->matrixOfPixels[heightIndex][widthIndex].green;
            unsigned char blue = matrixOfImage->matrixOfPixels[heightIndex][widthIndex].blue;
            fwrite(&red, 1, 1, file);
            fwrite(&green, 1, 1, file);
            fwrite(&blue, 1, 1, file);
        }
    }
    fclose(file);

}