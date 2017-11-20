#ifndef COMPRESSIONOFGRAPHICINFORMATIONBYNN_IMAGEPARSER_H
#define COMPRESSIONOFGRAPHICINFORMATIONBYNN_IMAGEPARSER_H

#include <stdio.h>

#pragma pack(push, 2)
typedef struct
{
    unsigned short int type;
    unsigned int size;
    unsigned short int reserved1;
    unsigned short int reserved2;
    unsigned int offset;
} BmpHeader;

typedef struct
{
    int headerSize;
    int width;
    int height;
    unsigned short planeCount;
    unsigned short bitsPerPixel;
    unsigned int compression;
    unsigned int imageSize;
    int xResolution;
    int yResolution;
    unsigned int colorsCount;
    unsigned int importantColors;

} BmpImageInfo;

typedef struct {
    unsigned char red;
    unsigned char green;
    unsigned char blue;
    unsigned char junk;
} RGB;
#pragma pack(pop)

typedef struct {
    RGB **matrixOfPixels;
    unsigned int height;
    unsigned int width;
} MatrixOfImage;

MatrixOfImage* getMatrixOfImage(char *imagePath);
RGB **getMatrixOfPixels(char *imagePath);
BmpHeader readBmpHeader(FILE *readImage);
BmpImageInfo readBmpImageInfo(FILE *readImage);
RGB* readBmpImagePalette(FILE *readImage, unsigned int colorsCount);
RGB** createMatrixOfPixels(unsigned int imageWidth, unsigned int imageHeight);
RGB** readMatrixOfPixels(FILE *readImage, unsigned int imageWidth, unsigned int imageHeight);
#endif //COMPRESSIONOFGRAPHICINFORMATIONBYNN_IMAGEPARSER_H