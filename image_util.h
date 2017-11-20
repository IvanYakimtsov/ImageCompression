#ifndef COMPRESSIONOFGRAPHICINFORMATIONBYNN_IMAGEPARSER_H
#define COMPRESSIONOFGRAPHICINFORMATIONBYNN_IMAGEPARSER_H

#include <stdio.h>

#pragma pack(push, 1)
typedef struct
{
    unsigned char	b1,b2;
    unsigned long	bfSize;
    unsigned short	bfReserved1;
    unsigned short	bfReserved2;
    unsigned long	bfOffBits;
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
//    unsigned char reserved;
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
void formImage(MatrixOfImage *matrixOfImage);
#endif //COMPRESSIONOFGRAPHICINFORMATIONBYNN_IMAGEPARSER_H