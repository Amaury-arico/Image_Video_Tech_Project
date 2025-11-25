// Image and Video Technology, Colas Schretter <colas.schretter@vub.be>
// This example program compares the C syntax for linear and multidimensional arrays
// Compilation: g++ -Wall -Wextra -pedantic -o ivt ivt_exercises.cpp


#include <iostream>
using std::cout;
using std::endl;

// storing images in memory with a one-dimensional array
float *create_image(const int height, const int width) {
    // dynamic memory allocation of a one-dimensional array
    float *image = new float[height * width];

    // fill in some values
    for(int y = 0; y < height; y++)
        for(int x = 0; x < 256; x++)
            image[y * width + x] = x + y * width;

    // return the pointer to the first image element
    return image;
}

// storing images in memory with a two-dimensional array
float (*create_image(const int height))[256] {
    // dynamic memory allocation of a two-dimensional array
    float (*image)[256] = new float[height][256];

    // fill in some values
    for(int y = 0; y < height; y++)
        for(int x = 0; x < 256; x++)
            image[y][x] = x + y * 256;

    // return the pointer to the first image element
    return image;
}

int main() {
    // create a two-dimensional array of 320x256 pixels
    float (*image2D)[256] = create_image(320);

    // create a one-dimensional array containing 81920 (320x256) pixels
    float *image1D = create_image(320,256);

    // access to a specific pixel using the 2D and 1D arrays
    const int y = 7, x = 5;
    cout << "the value of element at row number " << y << " and column number " << x << " is " << image2D[y][x] << endl;
    cout << "the value of element at row number " << y << " and column number " << x << " is " << image1D[y * 256 + x] << endl;

    // free the dynamic memory allocations
    delete [] image2D;
    delete [] image1D;

	return EXIT_SUCCESS;
}
