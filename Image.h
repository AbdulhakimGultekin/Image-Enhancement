/*
***********************************************************************************************************************
This library has been developed for "ELM 463 / Introduction to Image Processing" course of "Gebze Technical University"
and contains codes which are equivalent to some basic image enhancement functions. Only PGM format images can be
processed with this library.

Author		:	Abdülhakim Gültekin
Date		:	31.01.2018
Institute	:	Gebze Technical University
Course		:	ELM 463 / Introduction to Image Processing
Instructor	:	Assist. Prof. A. Köksal Hocaoğlu
***********************************************************************************************************************
*/

#ifndef IMAGE_H
#define IMAGE_H

class Image{
	public:
		// constructors and destructor
		Image(int = 0, int = 0);									// default constructor
		Image(const Image &);										// copy constructor
		~Image();													// destructor
		// create image
		void createImage(int, int);									// create an image, parameters all set
		void initImage(float init = 0.0);							// initiate the pixel value of an image, the default is 0.0
		// get and set functions
		int getRow() const;											// get row # / the height of the image
		int getCol() const;											// get col # / the width of the image
		float getMaximum() const;									// get the maximum pixel value
		float getMinimum() const;									// get the mininum pixel value
		float getPix(int rows, int cols);							// get pixel value at (rows, cols)
		Image getImage() const;										// get the image

		void setRow(int);											// set row number
		void setCol(int);											// set column number
		void setPix(int rows, int cols, float value);				// set pixel value at (rows, cols)
		void setImage(Image &);										// set the image
		// operator overloading functions
		float & operator()(int, int c = 0) const;					// operator overloading (i,j), when c = 0, a column vector
		const Image operator=(const Image &);						// overloading = operator 
		// read and write functions
		void readImage(char *fname);								// read the PGM format image file
		void writeImage(char *fname, bool flag = false);			// write the PGM image to the file
	
		// MEMBER FUNCTIONS FOR IMAGE ENHANCEMENT
		Image thresholdImage(float threshold = 127.0, float lowValue = 0.0, float highValue = 255.0);
		Image negativeImage();
		Image smoothImage(int dim);
		Image Gaussian(int dim);
		Image Sobel();
		Image HistogramEqualization();
		Image FourierTransform(bool flag = true);
		Image BasicThreshold();
		Image OtsuThreshold();
		Image uint2double();
		Image HoughImage();
		// END OF MEMBER FUNCTIONS
	private:
		int nrows;													// number of rows / height
		int ncols;													// number of columns / width
		int maximum;												// the maximum pixel value
		float *image;												// image buffer
};

#endif