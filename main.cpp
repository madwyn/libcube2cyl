#include "lib/qdbmp/qdbmp.h"
#include "Cube2Cyl.h"

using namespace std;

char* cubeNames[CUBE_FACE_NUM] =
{
    "top.bmp",
    "left.bmp",
    "front.bmp",
    "right.bmp",
    "back.bmp",
    "down.bmp"
};

int main()
{
    unsigned int i = 0;
    unsigned int j = 0;

    int xx = 0;
    int yy = 0;

    unsigned char rr;
    unsigned char gg;
    unsigned char bb;

    BMP *bmpCube[CUBE_FACE_NUM];

    // read the 6 images
    for (i = 0; i < CUBE_FACE_NUM; ++i)
    {
        bmpCube[i] = BMP_ReadFile(cubeNames[i]);

        if (BMP_GetError() != BMP_OK)
        {
            return 1;
        }
    }

	/* Get image's dimensions */
	int width  = BMP_GetWidth( bmpCube[0]);
	int height = BMP_GetHeight(bmpCube[0]);
    int depth  = BMP_GetDepth( bmpCube[0]);

    // the input images must be square
    if (width != height)
    {
        return 1;
    }

    // map the 6 images to 3D space
    Cube2Cyl algo;

    algo.init(width, M_PI, 2.0*M_PI);

    unsigned int panoWidth  = algo.pxPanoSizeH;
    unsigned int panoHeight = algo.pxPanoSizeV;

    // create panorama image
    BMP *output = BMP_Create(panoWidth, panoHeight, depth);

    // process
    for (i = 0; i < panoWidth; ++i)
    {
        for (j = 0; j < panoHeight; ++j)
        {
            algo.calXY(i, j, xx, yy);

            BMP_GetPixelRGB(bmpCube[algo.cubeFaceId], xx, yy, &rr, &gg, &bb);

            BMP_SetPixelRGB(output, i, j, rr, gg, bb);
        }
    }

    BMP_WriteFile(output, "pano.bmp");

    BMP_Free(output);

    for (i = 0; i < CUBE_FACE_NUM; ++i)
    {
        BMP_Free(bmpCube[i]);
    }

    return 0;
}
