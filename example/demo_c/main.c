<<<<<<< HEAD:example/demo_c/main.c
#include <cube_to_cyl.h>
#include <qdbmp.h>

int main() {
    static const char* cubeNames[CC_FACE_NUM] = {
        "TOP.bmp",
        "LEFT.bmp",
        "FRONT.bmp",
        "RIGHT.bmp",
        "BACK.bmp",
        "DOWN.bmp"
    };

    struct cc_context ctx;

    unsigned int i = 0;
    unsigned int j = 0;

    unsigned char rr;
    unsigned char gg;
    unsigned char bb;

    BMP *bmpCube[CC_FACE_NUM];

	unsigned int   width  = 0;
    unsigned int   height = 0;
    unsigned short depth  = 0;

    BMP *output = NULL;

    unsigned int pano_width  = 0;
    unsigned int pano_height = 0;

    const struct cc_coord* coord = NULL;

    // Read the 6 input images
    for (i = 0; i < CC_FACE_NUM; ++i) {
        bmpCube[i] = BMP_ReadFile(cubeNames[i]);

        if (BMP_GetError() != BMP_OK) {
            return 1;
        }
    }

	// Get image's dimensions
	width  = (unsigned int)BMP_GetWidth( bmpCube[0]);
	height = (unsigned int)BMP_GetHeight(bmpCube[0]);
    depth  = BMP_GetDepth( bmpCube[0]);

    // The input images must be square
    if (width != height) {
        return 1;
    }

    /*
       Initialise the algorithm:
         the width of each input is 640 pixel,
         the vertical view portion is PI (180 degrees),
         the horizontal view portion is 2*PI (360 degress).

       In this case, the output image size will be calculated accordingly.
       There is another more detailed init function you can play with.
     */
    cc_init(&ctx, width, M_PI*2.0, M_PI);

    // Generate the mapping from panorama to cubic
    cc_gen_map(&ctx);

    // Access the dimension of the panorama image
    pano_width  = ctx.px_pano_h;
    pano_height = ctx.px_pano_v;

    // Create the panorama output image
    output = BMP_Create(pano_width, pano_height, depth);

    // Map the pixels from the panorama back to the source image
    for (i = 0; i < pano_width; ++i) {
        for (j = 0; j < pano_height; ++j) {
            // Get the corresponding position of (i, j)
            coord = cc_get_coord(&ctx, i, j);

            // Access the pixel
            BMP_GetPixelRGB(bmpCube[coord->face], (unsigned long)coord->x, (unsigned long)coord->y, &rr, &gg, &bb);

            // Write the pixel to the panorama
            BMP_SetPixelRGB(output, i, j, rr, gg, bb);
        }
    }

    // Write the output file
    BMP_WriteFile(output, "PANO.bmp");

    // Release memory
    BMP_Free(output);

    for (i = 0; i < CC_FACE_NUM; ++i) {
        BMP_Free(bmpCube[i]);
    }

    cc_close(&ctx);

    return 0;
}
=======
#include "../cube_to_cyl.h"
#include "qdbmp/qdbmp.h"

int main() {
    static const char* cubeNames[CC_FACE_NUM] = {
        "TOP.bmp",
        "LEFT.bmp",
        "FRONT.bmp",
        "RIGHT.bmp",
        "BACK.bmp",
        "DOWN.bmp"
    };

    struct cc_context ctx;

    unsigned int i = 0;
    unsigned int j = 0;

    unsigned char rr;
    unsigned char gg;
    unsigned char bb;

    BMP *bmpCube[CC_FACE_NUM];

    int width  = 0;
    int height = 0;
    int depth  = 0;

    BMP *output = NULL;

    unsigned int pano_width  = 0;
    unsigned int pano_height = 0;

    const struct cc_coord* coord = NULL;

    // Read the 6 input images
    for (i = 0; i < CC_FACE_NUM; ++i) {
        bmpCube[i] = BMP_ReadFile(cubeNames[i]);

        if (BMP_GetError() != BMP_OK) {
            return 1;
        }
    }

    // Get image's dimensions
    width  = BMP_GetWidth( bmpCube[0]);
    height = BMP_GetHeight(bmpCube[0]);
    depth  = BMP_GetDepth( bmpCube[0]);

    // The input images must be square
    if (width != height) {
        return 1;
    }

    /*
       Initialise the algorithm:
         the width of each input is 640 pixel,
         the vertical view portion is PI (180 degrees),
         the horizontal view portion is 2*PI (360 degress).

       In this case, the output image size will be calculated accordingly.
       There is another more detailed init function you can play with.
     */
    cc_init(&ctx, width, M_PI*2.0, M_PI);

    // Generate the mapping from panorama to cubic
    cc_gen_map(&ctx);

    // Access the dimension of the panorama image
    pano_width  = ctx.px_pano_h;
    pano_height = ctx.px_pano_v;

    // Create the panorama output image
    output = BMP_Create(pano_width, pano_height, depth);

    // Map the pixels from the panorama back to the source image
    for (i = 0; i < pano_width; ++i) {
        for (j = 0; j < pano_height; ++j) {
            // Get the corresponding position of (i, j)
            coord = cc_get_coord(&ctx, i, j);

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

    for (i = 0; i < CC_FACE_NUM; ++i) {
        BMP_Free(bmpCube[i]);
    }

    cc_close(&ctx);

    return 0;
}
>>>>>>> d0f7c516e3da66bb1fe9abf63324c5b574abe774:demo/main.c
