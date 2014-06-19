#include "qdbmp/qdbmp.h"
#include "../Cube2Cyl.h"

#include <iostream>

using namespace std;

const char* cubeNames[CUBE_FACE_NUM] =
{
    "TOP.bmp",
    "LEFT.bmp",
    "FRONT.bmp",
    "RIGHT.bmp",
    "BACK.bmp",
    "DOWN.bmp"
};

int main()
{
    unsigned int i = 0;
    unsigned int j = 0;

    unsigned char rr;
    unsigned char gg;
    unsigned char bb;

    BMP *bmpCube[CUBE_FACE_NUM];

    // Read the 6 input images
    for (i = 0; i < CUBE_FACE_NUM; ++i)
    {
        bmpCube[i] = BMP_ReadFile(cubeNames[i]);

        if (BMP_GetError() != BMP_OK)
        {
            return 1;
        }
    }

	// Get image's dimensions
	int width  = BMP_GetWidth( bmpCube[0]);
	int height = BMP_GetHeight(bmpCube[0]);
    int depth  = BMP_GetDepth( bmpCube[0]);

    // The input images must be square
    if (width != height)
    {
        return 1;
    }

    // Create a instance of Cube2Cyl algorithm
    Cube2Cyl algo;

    /*
       Initialise the algorithm:
         the width of each input is 640 pixel,
         the vertical view portion is PI (180 degrees),
         the horizontal view portion is 2*PI (360 degress).

       In this case, the output image size will be calculated accordingly.
       There is another more detailed init function you can play with.
     */
    algo.init(width, M_PI, 2.0*M_PI);
    // Generate the mapping from panorama to cubic
    algo.genMap();

    // Access the dimension of the panorama image
    unsigned int panoWidth  = algo.pxPanoSizeH;
    unsigned int panoHeight = algo.pxPanoSizeV;

    // Create the panorama output image
    BMP *output = BMP_Create(panoWidth, panoHeight, depth);

    const CUBE_COORD* coord = NULL;

    // Map the pixels from the panorama back to the source image
    for (i = 0; i < panoWidth; ++i)
    {
        for (j = 0; j < panoHeight; ++j)
        {
            // Get the corresponding position of (i, j)
            coord = algo.getCoord(i, j);

            // Access the pixel
            BMP_GetPixelRGB(bmpCube[coord->face], coord->x, coord->y, &rr, &gg, &bb);

            // Write the pixel to the panorama
            BMP_SetPixelRGB(output, i, j, rr, gg, bb);
        }
    }

    // Write the output file
    BMP_WriteFile(output, "PANO.bmp");

    // Release memory
    BMP_Free(output);

    for (i = 0; i < CUBE_FACE_NUM; ++i)
    {
        BMP_Free(bmpCube[i]);
    }

    return 0;
}
