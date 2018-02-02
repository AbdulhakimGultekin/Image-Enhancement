/*
**********************************************************************************************************
An example to dictate how to use the "Image" class with PGM format images. This program reads an image and
makes a histogram equalization process on the image. Run the code and see the results.
**********************************************************************************************************  
*/

#include "Image.h"

int main(){

// initializing objects that belong to input images
Image firstInputImage;

// input file name
char firstInputImageFileName[] = "Exp.pgm";

// initializing objects that belong to output images
Image firstOutputImage;

// output file name
char fourier_with_scaling[] = "HistEq.pgm";

firstInputImage.readImage(firstInputImageFileName);
firstOutputImage = firstInputImage.HistogramEqualization();
firstOutputImage.writeImage(fourier_with_scaling, true);

return 0;
}