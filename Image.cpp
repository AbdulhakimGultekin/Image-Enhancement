#include "Image.h"

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

// default constructor
Image::Image(int nRows, int nCols){
    if (nRows < 0 || nCols < 0){
        cout << "Image: Index out of range.\n";
        exit(1);
    }
    image = NULL;
    nrows = nRows;
    ncols = nCols;
    maximum = 255;
    if(nrows > 0 && ncols > 0)
        createImage(nrows, ncols);
}

// copy constructor
Image::Image(const Image &img) {
    image = NULL;
    nrows = img.getRow();
    ncols = img.getCol();
    createImage(nrows, ncols);

    for (int rows = 0; rows < nrows; rows++)
        for (int cols = 0; cols < ncols; cols++)
            image[rows * ncols + cols] = img(rows, cols);
}

// destructor
Image::~Image(){
    if(image)
        delete [] image;
}

// Allocates memory and initialize the image.
void Image::createImage(int nRows, int nCols){
    if(image)
        delete [] image;
    nrows = nRows;
    ncols = nCols;
    image = (float *)new float[nrows * ncols];
    if(!image){
        cout << "CREATEIMAGE: Out of memory!";
        exit(1);
    }
    initImage();			
}

// Initialize image.
void Image::initImage(float initVal){
    for(int numofPix = 0; numofPix < nrows * ncols; numofPix++)
        image[numofPix] = initVal;
}

// Returns the number of rows.
int Image::getRow() const {
    return nrows;
}

// Returns the number of cols.
int Image::getCol() const {
    return ncols;
}

// Returns the maximum number.
float Image::getMaximum() const {
    float max = -10000;
    for(int rows = 0; rows < nrows; rows++){
        for(int cols = 0; cols < ncols; cols++){
            if (max < image[rows * ncols + cols])
                max = image[rows * ncols + cols];
        }
    }
    return max;
}

// Returns the minimum number.
float Image::getMinimum() const {
    float min = 10000;
    for(int rows = 0; rows < nrows; rows++){
        for(int cols = 0; cols < ncols; cols++){
            if (min > image[rows * ncols + cols])
                min = image[rows * ncols + cols];
        }
    }
    return min;
}

// Returns the pixel value.
float Image::getPix(int Rows, int Cols){
    return image[Rows * ncols + Cols];
}

// Returns the image using image buffer.
Image Image::getImage() const {
    Image temp;
    temp.createImage(nrows, ncols);
    for (int rows = 0; rows < nrows; rows++)
        for (int cols = 0; cols < ncols; cols++)
            temp(rows, cols) = image[rows * ncols + cols];

    return temp;
}

// Sets the row number.
void Image::setRow(int numberOfRows){
    nrows = numberOfRows;
}

// Sets the column number.
void Image::setCol(int numberOfColumns){
    ncols = numberOfColumns;
}

// Sets the pixel value.
void Image::setPix(int rows, int cols, float value){
    image[rows * ncols + cols] = value;
}

// Sets the image.
void Image::setImage(Image &img){
    int rows, cols;

    for (rows = 0; rows < nrows; rows++)
        for (cols = 0; cols < ncols; cols++)
            image[rows * ncols + cols] = img(rows, cols);
}

float & Image::operator()(int rows, int cols) const {
    return image[rows * ncols + cols];
}

const Image Image::operator=(const Image& img) {
    int rows, cols;

    if (this == &img)
        return *this;

    nrows = img.getRow();
    ncols = img.getCol();
    createImage(nrows, ncols);

    for (rows = 0; rows < nrows; rows++)
		for (cols = 0; cols < ncols; cols++)
            (*this)(rows, cols) = img(rows, cols);
	
    return *this;
}

/*
Read image from a file.
*/
	void Image::readImage(char *fname) {
	ifstream ifp;
	char dummy[80];
	unsigned char *img;
	int rows, cols;
	int nRows, nCols, nt, maxi;

	ifp.open(fname, ios::in | ios::binary);

	if (!ifp) {
		cout << "readImage: Can't read image: " << fname << endl;
		exit(1);
	}

	// identify image format
	ifp.getline(dummy, 80, '\n');

	if (dummy[0] == 'P' && dummy[1] == '5')
		;
	else {
		cout << "readImage: Can't identify image format." << endl;
		exit(1);
	}

	// skip the comments
	ifp.getline(dummy, 80, '\n');

	while (dummy[0] == '#') {
		ifp.getline(dummy, 80, '\n');
	}

	// read the row number and column number
	sscanf(dummy, "%d %d", &nCols, &nRows);

	// read the maximum pixel value
	ifp.getline(dummy, 80, '\n');
	sscanf(dummy, "%d", &maxi);
	if (maxi > 255) {
		cout << "Don't know what to do: maximum value is over 255.\n";
		exit(1);
	}

	if (image != NULL)
		delete [] image;

	nrows = nRows;
	ncols = nCols;
	maximum = 255;

	// read the image data
	img = (unsigned char *) new unsigned char [nRows * nCols];
	if (!img) {
		cout << "READIMAGE: Out of memory.\n";
		exit(1);
	}
	image = (float *) new float [nRows * nCols];
	if (!image) {
		cout << "READIMAGE: Out of memory.\n";
		exit(1);
	}

    ifp.read((char *)img, (nRows * nCols * sizeof(unsigned char)));

    for (rows = 0; rows < nRows; rows++)
        for (cols = 0; cols < nCols; cols++)
            image[rows * nCols + cols] = (float) img[rows * nCols + cols];

    ifp.close();

    delete [] img;
}


/*
Write image buffer to a file.
*/
void Image::writeImage(char *fname, bool flag) {
	ofstream ofp;
	int i, j;
	int nRows, nCols, nt;
	unsigned char *img;

	ofp.open(fname, ios::out | ios::binary);

	if (!ofp) {
		cout << "writeImage: Can't write image: " << fname << endl;
		exit(1);
	}


	ofp << "P5" << endl;
	ofp << ncols << " " << nrows << endl;


	ofp << 255 << endl;

	// convert the image data type back to unsigned char
	img = (unsigned char *) new unsigned char [nrows * ncols];
	if (!img) {
		cout << "WRITEIMAGE: Out of memory.\n";
		exit(1);
	}

	float maxi = getMaximum();
	float mini = getMinimum();


    for (i = 0; i< nrows; i++){
        for (j = 0; j < ncols; j++){
            // rescale if the flag is set
		    if ((maxi != mini) && flag == true)
			    img[i * ncols + j] = (unsigned char)  ((image[i * ncols + j]-mini)/(float)(maxi-mini)*255.0);
		    // any intensity that is larger than the maximum would be set as maximum
		    else if (image[i * ncols + j] > 255)
			    img[i * ncols + j] = 255;
		    else if (image[i * ncols + j] < 0)
			    img[i * ncols + j] = 0;
		    else
			img[i * ncols + j] = (unsigned char)  image[i * ncols + j];
        }
    }
	
    ofp.write((char *)img, (nrows * ncols * sizeof(unsigned char)));

    ofp.close();
    delete [] img;
}

// MEMBER FUNCTIONS FOR IMAGE ENHANCEMENT

Image Image::thresholdImage(float thresholdValue, float lowValue, float highValue) {
	Image temp;
	int rows, cols;

	temp.createImage(nrows, ncols);   // temp is a gray-scale image
	for (rows = 0; rows < nrows; rows++)
		for (cols = 0; cols < ncols; cols++)
			if (image[rows * ncols + cols] <= thresholdValue)
				temp(rows, cols) = lowValue;
			else
				temp(rows, cols) = highValue;
	return temp;
}

// "negativeImage" function makes the image negatived.
Image Image::negativeImage(){
    int L = 256;
    Image temp;
    temp.createImage(nrows, ncols);
    // So as to get negative of an image used method is s = L - 1 - r.
    for (int rows = 0; rows < nrows; rows++)
        for (int cols = 0; cols < ncols; cols++)
            temp.setPix(rows, cols, L - 1 - image[rows * ncols + cols]);  // it's a way to set

    return temp;
}
// "smoothImage" function can be used to blur any given image by using dimXdim average filter.
Image Image::smoothImage(int dim){
	Image temp;
	temp.createImage(nrows, ncols);
	float result = 0;
	float temp_val = 0;
	float var;

	// Blurring the image with an averaging mask by regarding neighborhood specifications.
	for (int rows = 0; rows < nrows; rows++){
        for (int cols = 0; cols < ncols; cols++){
            for (int i = -1 * ((dim - 1) / 2); i < ((dim - 1) / 2 + 1); i++){
                for(int k = -1 * ((dim - 1) / 2); k < ((dim - 1) / 2 + 1); k++){
        			if((rows + i < 0) || (cols + k < 0) || (rows + i > nrows - 1) || (cols + k > ncols - 1))
        				var = 0;
        			else
        				var = this->getPix(rows + i, cols + k);

                    temp_val = temp_val + var;
        	    }
            }
			// Result is the average value for every pixels (x,y).
            result = (1.0 / (dim * dim)) * temp_val;
            temp.setPix(rows, cols, result);
            temp_val = 0;
        }
    }
    return temp;
}

// "Gaussian" function can be used for smoothing purposes.
Image Image::Gaussian(int dim){
	Image temp;
	temp.createImage(nrows, ncols);
	Image h(dim,dim);
	float temp_val = 0;
	float var;
	float sigma = 3.0;
    float result;

    // Generating Gaussian spatial filter
    for (int rows = 0; rows < h.getRow(); rows++)
        for (int cols = 0; cols < h.getCol(); cols++){
            h.setPix(rows, cols, exp(-1 * (rows * rows + cols * cols) / (2 * sigma * sigma)));
        }

    for (int rows = 0; rows < nrows; rows++){
        for (int cols = 0; cols < ncols; cols++){
            for (int i = -1 * ((dim - 1) / 2); i < ((dim - 1) / 2 + 1); i++){
        	    for(int k = -1 * ((dim - 1) / 2); k < ((dim - 1) / 2 + 1); k++){
        		    if((rows + i < 0) || (cols + k < 0) || (rows + i > nrows - 1) || (cols + k > ncols - 1))
        			    var = 0;
                    else
        			    var = this->getPix(rows + i,cols + k) * h.getPix(i + ((dim - 1) / 2), k + ((dim - 1) / 2));
					
                    temp_val = temp_val + var;
                }
            }
            result = (1.0 / (dim * dim)) * temp_val;
            temp.setPix(rows, cols, result);
            temp_val = 0;
        }
    }
    return temp;
}

// "Sobel" function has been generated to find edges of any given image by using Sobel operators.
Image Image::Sobel(){
	Image temp1(this->getRow(), this->getCol());
	Image temp2(this->getRow(), this->getCol());
	Image temp_ret(this->getRow(), this->getCol());
	float temp_val1 = 0;
	float temp_val2 = 0;
	float var1, var2;
	const int width = 3;
	const int height = 3;
	// Sobel operators (filters/masks) generation
	int sobely[width][height] = {{-1,0,1},{-2,0,2},{-1,0,1}};
	int sobelx[width][height] = {{-1,-2,-1},{0,0,0},{1,2,1}};


	Image sobelX(3,3);
	Image sobelY(3,3);

	for (int rows = 0; rows < sobelX.getRow(); rows++)
        for (int cols = 0; cols < sobelY.getCol(); cols++){
            sobelX.setPix(rows, cols, sobelx[rows][cols]);
		    sobelY.setPix(rows, cols, sobely[rows][cols]);
        }

	for (int rows = 0; rows < this->getRow(); rows++)
        for (int cols = 0; cols < this->getCol(); cols++){
            for (int i = -1; i < 2; i++){
                for(int k = -1; k < 2; k++){
        		    if((rows + i < 0) || (cols + k < 0) || (rows + i > nrows - 1) || (cols + k > ncols - 1 )){
        			    var1 = 0;
        			    var2 = 0;
                    }
                    else{
        			    var1 = this->getPix(rows + i, cols + k) * sobelX.getPix(i + 1, k + 1);
        			    var2 = this->getPix(rows + i, cols + k) * sobelY.getPix(i + 1, k + 1);
                    }

                    temp_val1 = temp_val1 + var1;
                    temp_val2 = temp_val2 + var2;
                }
            }

            temp1.setPix(rows, cols, temp_val1);
            temp2.setPix(rows, cols, temp_val2);
            temp_val1 = 0;
            temp_val2 = 0;
        }
	// sqrt(Gx.^2 + Gy.^2) operation is being handled
    for (int rows = 0; rows < this->getRow(); rows++)
        for (int cols = 0; cols < this->getCol(); cols++)
            temp_ret.setPix(rows, cols, sqrt(pow(temp1.getPix(rows, cols), 2) + pow(temp2.getPix(rows, cols), 2)));

    return temp_ret;
}

/* "HistogramEqualization" function has been generated to use in histogram equalization operations. This operation makes a continuous signal has a uniform
   histogram. But in discrete cases it won't be as uniform as continuous cases are. Histogram equalization is useful, when an image has a gaussian histogram,
   to make the image more sharp.
*/
Image Image::HistogramEqualization(){
    Image temp(nrows, ncols);
    int histogram[256];
    int cumhistogram[256];
	int Sk[256];
	double PrRk[256];

    // initialize all intensity values with 0
    for(int i = 0; i < 256; i++)
        histogram[i] = 0;

    // calculate the number of pixels for each intensity values
    for(int rows = 0; rows < this->nrows; rows++)
        for(int cols = 0; cols < this->ncols; cols++)
            histogram[(int)this->getPix(rows, cols)]++;

    // set the cumulative histogram
    cumhistogram[0] = histogram[0];
    for(int i = 1; i < 256; i++)
        cumhistogram[i] = histogram[i] + cumhistogram[i-1];

    // calculate the size of image
    int size = this->getRow() * this->getCol();
    float alpha = 255.0 / size; // alpha = (L-1)/M*N

    // calculate the probability of each intensity
    for(int i = 0; i < 256; i++)
        PrRk[i] = (double)histogram[i] / size;

    // scale the histogram
    for(int i = 0; i < 256; i++)
        Sk[i] = round((double)cumhistogram[i] * alpha); // round to the nearest integer value

    // set the new values of corresponding pixel values
    for(int rows = 0; rows < this->nrows; rows++)
        for(int cols = 0; cols < this->ncols; cols++)
            temp.setPix(rows, cols, Sk[(int)this->getPix(rows, cols)]);

    return temp;
}

// "FourierTransform" function performs a transformation to Fourier space for a given image. This function will be controlled.
Image Image::FourierTransform(bool flag){
	Image temp(nrows, ncols);

	// creating a structure type which, for all samples, has two variables that are real and imaginary parts of a complex number
	struct NewComplex{
		float real;
		float imag;
	};

	NewComplex complexImage[nrows * ncols];

	// initializing complex image array
    for(int rows = 0; rows < nrows; rows++){
        for(int cols = 0; cols < ncols; cols++){
            complexImage[rows * ncols + cols].real = 0;
            complexImage[rows * ncols + cols].imag = 0;
	    }
    }
	
    float result = 0;
    float pi = 3.1419265;
    float real = 0;
    float imag = 0;


    // 1-D fourier transform for rows
    for(int rows = 0; rows < nrows; rows++){
        for(int v = 0; v < ncols; v++){
            for(int cols = 0; cols < ncols; cols++){
			    real += this->getPix(rows, cols) * cos((2 * pi * v * cols) / ncols);
			    imag += this->getPix(rows, cols) * -1 * sin((2 * pi * v * cols) / ncols);
            }
            complexImage[rows * ncols + v].real = real;
            complexImage[rows * ncols + v].imag = imag;
            real = 0;
            imag = 0;
        }
    }

    // then 1-D fourier transform for columns
    for(int v = 0; v < ncols; v++){
        for(int u = 0; u < nrows; u++){
            for(int rows = 0; rows < nrows; rows++){
                // there is a minus in the middle of the equation since j * j = -1 equality occurs there
                real += complexImage[rows * ncols + u].real * cos((2* pi * u * rows) / nrows) - complexImage[rows * ncols + u].imag * -1 * sin((2* pi * u * rows) / nrows);
                imag += complexImage[rows * ncols + u].real * -1 * sin((2* pi * u * rows) / nrows) + complexImage[rows * ncols + u].imag * cos((2* pi * u * rows) / nrows);
            }
            result = sqrt(pow(real,2) + pow(imag,2));
            temp.setPix(u, v, result); // set new pixel values on temporal image object
            real = 0;
            imag = 0;
        }
    }

    // if flag is true, logarithmic scaling will be applied on the spectrum
    if(flag == true){
        for(int rows = 0; rows < nrows; rows++)
            for(int cols = 0; cols < ncols; cols++)
                temp.setPix(rows, cols, log10(1 + abs(temp.getPix(rows, cols))));
    }

    return temp;
}

// "BasicThreshold" function performs a simple thresholding process on a given image.
Image Image::BasicThreshold(){
	Image temp(nrows,ncols);
	float thVal = 0.0;
	float max = 0.0;

	for(int row = 0; row < nrows; row++)
		for(int col = 0; col < ncols; col++){
			if(this->getPix(row, col) > max)
				max = this->getPix(row,col);
		}

	thVal = (float)0.33 * max;
	for(int row = 0; row < nrows; row++)
		for(int col = 0; col < ncols; col++){
			if(this->getPix(row, col) < thVal)
				temp.setPix(row, col, 0);
			else
				temp.setPix(row, col, 255);
		}

    return temp;
}

// "OtsuThreshold" performs a thresholding regarding Otsu's histogram processes as the name indicates. Will be reviewed. Not usable for now.
Image Image::OtsuThreshold(){
    // initializing an image object with sizes of interested image
    Image temp(nrows, ncols);
    // initializing all arrays needed to calculate max. between group variance value
    int histogram[256];
    float q1[256];
    float av1[256];
    float av2[256];
    float betVar[256];
    float glMean = 0.0;
	float max = 0;
    int th;

    // initialize all intensity values with 0
    for(int i = 0; i < 256; i++)
        histogram[i] = 0;

    // calculate the number of pixels for each intensity values
    for(int rows = 0; rows < this->nrows; rows++)
        for(int cols = 0; cols < this->ncols; cols++)
            histogram[(int)this->getPix(rows, cols)]++;

    // calculate the size of image
    int size = this->getRow() * this->getCol();

    // calculate the probability of each intensity, normalized histogram
    float PrRk[256];
	for(int i = 0; i < 256; i++)
		PrRk[i] = 0;

    for(int i = 0; i < 256; i++){
        PrRk[i] = (float)histogram[i] / size;
        //cout << "PrRk[" << i << "] = " << PrRk[i] << endl;
    }
    //  global mean -glMean- calculation
    for(int i = 0; i < 256; i++)
        glMean += i * PrRk[i];

    // initializing arrays
    for(int i = 0; i < 256; i++){
        q1[i] = 0;
        av1[i] = 0;
        av2[i] = 0;
        betVar[i] = 0;
    }

    // calculations for between class variance -betVar-
    for(int i = 1; i < 256; i++){
        q1[i] = q1[i-1] + PrRk[i];
        //cout << "q1[" << i << "] = " << q1[i] << endl;
        av1[i] = ((q1[i-1] * av1[i-1]) + (i * PrRk[i])) / (q1[i] + 0.000001);
        av2[i] = (glMean - (q1[i] * av1[i])) / (1 - q1[i]);
        //cout << "av2[" << i << "] = " << av2[i] << endl;
        betVar[i] = q1[i] * (1 - q1[i]) * pow((av1[i] - av2[i]), 2);
        //cout << "betVar[" << i << "] = " << betVar[i] << endl;
		if(betVar[i] > max){
            max = betVar[i];
            th = i;
		}
    }
    //cout << "th val = " << th << endl;

    // thresholding the image by regarding the pixel intensity value founded at previous step
    for(int rows = 0; rows < nrows; rows++)
        for(int cols = 0; cols < ncols; cols++){
            if(this->getPix(rows, cols) > th)
                temp.setPix(rows, cols, 255);
            else
                temp.setPix(rows, cols, 0);
    }

    return temp;
}

// "uint2double" function takes an image which is in 'uint8' type and converts it to 'double'.
Image Image::uint2double(){
    Image temp(nrows, ncols);

    for(int rows = 0; rows < nrows; rows++)
        for(int cols; cols < ncols; cols++)
            temp.setPix(rows, cols, this->getPix(rows, cols) * (1.0 / 255));

    return temp;
}

// This is to generate a specific image.
Image Image::HoughImage(){
    int X = 101, Y = 101;
    Image temp(X, Y);

    for(int rows = 0; rows < X; rows++)
        for(int cols = 0; cols < Y; cols++)
            if(rows == 0 && cols == 0)
                temp.setPix(rows, cols, 255);
            else if(rows == 0 && cols == Y -1)
                temp.setPix(rows, cols, 255);
            else if(rows == (X - 1) / 2 && cols == (Y - 1) / 2)
                temp.setPix(rows, cols, 255);
            else if(rows == X - 1 && cols == 0)
                temp.setPix(rows, cols, 255);
            else if(rows == X - 1 && cols == Y - 1)
                temp.setPix(rows, cols, 255);
            else
                temp.setPix(rows, cols, 0);

    return temp;
}
// END OF MEMBER FUNCTIONS