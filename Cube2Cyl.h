#ifndef CUBE2CYL_H_INCLUDED
#define CUBE2CYL_H_INCLUDED

/*  Cube2Cyl v1.0.0 - 2012-10-10
 *
 *  Cube2Cyl is a cubic projection to cylindrical projection conversion lib.
 *
 *  Homepage: http://www.wenyanan.com/cube2cyl/
 *  Please check the web page for further information, upgrades and bug fixes.
 *
 *  Copyright (c) 2012 Yanan Wen (WenYananWork@gmail.com)
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 *
 ******************************************************************************/

#include <math.h>

#ifndef M_PI
#define M_PI        3.14159265358979323846
#define M_PI_2      1.57079632679489661923
#define M_PI_4      0.78539816339744830962
#endif

enum CUBE_FACES
{
    CUBE_TOP = 0,
    CUBE_LEFT,
    CUBE_FRONT,
    CUBE_RIGHT,
    CUBE_BACK,
    CUBE_DOWN,
    CUBE_FACE_NUM
};

class Cube2Cyl
{
public:

    //-------- input information
    unsigned int pxCamV;    /**< The vertical pixel of a camera */
    unsigned int pxCamH;    /**< The horizontal pixel of a camera */

    double diameter;        /**< The diameter of the sphere */

    double radPanoV;        /**< The vertical view portion */
    double radPanoH;        /**< The horizontal view portion */

    //-------- output information
    unsigned int pxPanoSizeV;   /**< The vertical pixels of the panorama */
    unsigned int pxPanoSizeH;   /**< The horizontal pixels of the panorama */

    unsigned int cubeFaceId;    /**< The cube face to be read */

    void init(unsigned int pxInW, double radInV, double radInH);
    void init(unsigned int pxPanoH, unsigned int pxPanoV, unsigned int pxInW, double radInV, double radInH);

    inline void calXY(const int& i, const int& j, int& xx, int& yy);

    Cube2Cyl(void);
    ~Cube2Cyl(void);

private:

    inline void calXYZ(const int& i, const int& j, double& x, double& y, double& z);

    inline void calNormXY(const int& i, const int& j, double& x, double& y);
    inline void calThetaAndPhi(double& x, double& y, double& theta, double& phi);
    inline void calXyzFromThetaPhi(double& theta, double& phi, double& x, double& y, double& z);

    inline void calCubeFace(const double& theta, const double& phi);

    inline void locateTop(   const double& x, const double& y, const double& z, int& xx, int& yy);
    inline void locateDown(  const double& x, const double& y, const double& z, int& xx, int& yy);
    inline void locateFront( const double& x, const double& y, const double& z, int& xx, int& yy);
    inline void locateBack(  const double& x, const double& y, const double& z, int& xx, int& yy);
    inline void locateLeft(  const double& x, const double& y, const double& z, int& xx, int& yy);
    inline void locateRight( const double& x, const double& y, const double& z, int& xx, int& yy);

    // the helper functions
    inline bool cmpDoubleEqual(        const double &a, const double &b, const double &epsilon);
    inline bool cmpDoubleSmaller(      const double &a, const double &b, const double &epsilon);
    inline bool cmpDoubleEqualSmaller( const double &a, const double &b, const double &epsilon);
    inline bool Cube2CylcmpDoubleLager(const double &a, const double &b, const double &epsilon);
    inline bool cmpDoubleEqualLager(   const double &a, const double &b, const double &epsilon);

    inline void rotRad(double rad, double& x, double& y, double& temp);
    inline void transDis(double dis, double& x, double& y);


    //-------- The temp variables
    double normTheta;   /**< The normalised theta */
    double resCal;      /**< The resolution used for calculation */
    double normFactorX; /**< The normalisation factor for x */
    double normFactorY; /**< The normalisation factor for y */

    double sizeRatio;   /**< The size ratio of the mapped x and the actual diameter */
    double mappedX;     /**< The x coordinate mapped on the cube face */
    double mappedY;     /**< The y coordinate mapped on the cube face */

    double tX;           /**< x coordinate in 3D space */
    double tY;           /**< y coordinate in 3D space */
    double tZ;           /**< z coordinate in 3D space */

    double tTheta;       /**< The radian horizontally */
    double tPhi;         /**< The radian vertically */

    double phiThreshold;    /**< The threshold of phi, it seprates the top, middle and down
     of the cube, which are the edges of the top and down surface */
};

/** \brief  The initialise function of Cube2Cyl. This one takes fewer parameters, and will
 *          generate the width and height of the panorama automatically based on the input
 *
 * \param pxInW unsigned int    The width of the input image
 * \param radInV double         The radian of the view portion vertically, range [0.01, PI]
 * \param radInH double         The radian of the view portion horizontally, range [0.01, 2*PI]
 * \return void
 *
 */
void Cube2Cyl::init(unsigned int pxInW, double radInV, double radInH)
{
    unsigned int pxPanoH = (unsigned int)((double)(radInH/M_PI_2) * pxInW);
    unsigned int pxPanoV = (unsigned int)((double)(radInV/M_PI_2) * pxInW);
    init(pxPanoH, pxPanoV, pxInW, radInV, radInH);
}

/** \brief  The initialise function of Cube2Cyl.
 *
 * \param pxPanoH unsigned int  The desired panorama width
 * \param pxPanoV unsigned int  The desired panorama height
 * \param pxInW unsigned int    The width of the input image
 * \param radInV double         The radian of the view portion vertically, range [0.01, PI]
 * \param radInH double         The radian of the view portion horizontally, range [0.01, 2*PI]
 * \return void
 *
 */
void Cube2Cyl::init(unsigned int pxPanoH, unsigned int pxPanoV, unsigned int pxInW, double radInV, double radInH)
{
    // check parameters
    if (   (pxInW == 0)
        || (pxPanoH == 0)
        || (pxPanoV == 0))
    {
        return;
    }

    // check the view portion
    if (   ((radInV<0.01) || (radInV>M_PI))
        || ((radInH<0.01) || (radInH>(M_PI*2.0))))
    {
        return;
    }

    pxCamV = pxInW;
    pxCamH = pxInW;

    radPanoV = radInV;
    radPanoH = radInH;

    pxPanoSizeH = pxPanoH;
    pxPanoSizeV = pxPanoV;

    diameter = ((double)pxInW) / 2.0;

    // the actual calculation resolution is 10 times bigger than the texture resolution
    resCal = M_PI_4 / (double)pxCamH / 10.0;

    // the normalisation factors
    normFactorX = radPanoH/(M_PI*2);
    normFactorY = radPanoV/M_PI;
}

/** \brief Rotate the point for a given radian
 *
 * \param rad double
 * \param x double&
 * \param y double&
 * \param temp double&
 * \return void
 *
 */
inline void Cube2Cyl::rotRad(double rad, double& x, double& y, double& temp)
{
    temp = x;
    x = x*cos(rad) - y*sin(rad);
    y = temp*sin(rad) + y*cos(rad);
}

/** \brief Translate the point for a given distance in both x and y coordinates
 *
 * \param dis double
 * \param x double&
 * \param y double&
 * \return void
 *
 */
inline void Cube2Cyl::transDis(double dis, double& x, double& y)
{
    x += dis;
    y += dis;
}

/** \brief Calculate the x and y coordinates of given i and j
 *
 * \param i const int&  The coordinate along the width axis
 * \param j const int&  The coordinate along the height axis
 * \param xx int&       The x coordinate on the cube face
 * \param yy int&       The y coordinate on the cube face
 * \return void
 *
 */
inline void Cube2Cyl::calXY(const int& i, const int& j, int& xx, int& yy)
{
    calXYZ(i, j, tX, tY, tZ);

    switch (cubeFaceId)
    {
    case CUBE_TOP:
        {
            locateTop(tX, tY, tZ, xx, yy);
            break;
        }
    case CUBE_DOWN:
        {
            locateDown(tX, tY, tZ, xx, yy);
            break;
        }
    case CUBE_LEFT:
        {
            locateLeft(tX, tY, tZ, xx, yy);
            break;
        }
    case CUBE_RIGHT:
        {
            locateRight(tX, tY, tZ, xx, yy);
            break;
        }
    case CUBE_FRONT:
        {
            locateFront(tX, tY, tZ, xx, yy);
            break;
        }
    case CUBE_BACK:
        {
            locateBack(tX, tY, tZ, xx, yy);
            break;
        }
    default:
        {
            break;
        }
    }
}

/** \brief Locate the point in the top image
 *
 * \param x const double&
 * \param y const double&
 * \param z const double&
 * \param xx int&
 * \param yy int&
 * \return void
 *
 */
inline void Cube2Cyl::locateTop(const double& x, const double& y, const double& z, int& xx, int& yy)
{
    sizeRatio = diameter / y;

    mappedX = sizeRatio * z;
    mappedY = sizeRatio * x;

    // rotate pi
    rotRad(M_PI, mappedX, mappedY, sizeRatio);

    // translate
    transDis(diameter, mappedX, mappedY);

    xx = (int) mappedX;
    yy = (int) mappedY;
}

/** \brief Locate the point in the down image
 *
 * \param x const double&
 * \param y const double&
 * \param z const double&
 * \param xx int&
 * \param yy int&
 * \return void
 *
 */
inline void Cube2Cyl::locateDown(const double& x, const double& y, const double& z, int& xx, int& yy)
{
    sizeRatio = diameter / y;

    mappedX = sizeRatio * x;
    mappedY = sizeRatio * z;

    // rotate -pi
    rotRad(-M_PI_2, mappedX, mappedY, sizeRatio);

    // translate
    transDis(diameter, mappedX, mappedY);

    xx = (int) mappedX;
    yy = (int) mappedY;
}

/** \brief Locate the point in the front image
 *
 * \param x const double&
 * \param y const double&
 * \param z const double&
 * \param xx int&
 * \param yy int&
 * \return void
 *
 */
inline void Cube2Cyl::locateFront(const double& x, const double& y, const double& z, int& xx, int& yy)
{
    sizeRatio = diameter / x;

    mappedX = sizeRatio * z;
    mappedY = sizeRatio * y;

    // translate
    transDis(diameter, mappedX, mappedY);

    xx = (int) mappedX;
    yy = (int) mappedY;
}

/** \brief Locate the point in the back image
 *
 * \param x const double&
 * \param y const double&
 * \param z const double&
 * \param xx int&
 * \param yy int&
 * \return void
 *
 */
inline void Cube2Cyl::locateBack(const double& x, const double& y, const double& z, int& xx, int& yy)
{
    sizeRatio = diameter / x;

    mappedX = sizeRatio * y;
    mappedY = sizeRatio * z;

    // rotate -pi
    rotRad(-M_PI_2, mappedX, mappedY, sizeRatio);

    // translate
    transDis(diameter, mappedX, mappedY);

    xx = (int) mappedX;
    yy = (int) mappedY;
}


/** \brief Locate the point in the left image
 *
 * \param x const double&
 * \param y const double&
 * \param z const double&
 * \param xx int&
 * \param yy int&
 * \return void
 *
 */
inline void Cube2Cyl::locateLeft(const double& x, const double& y, const double& z, int& xx, int& yy)
{
    sizeRatio = diameter / z;

    mappedX = sizeRatio * x;
    mappedY = sizeRatio * y;

    // rotate pi
    rotRad(M_PI, mappedX, mappedY, sizeRatio);

    // translate
    transDis(diameter, mappedX, mappedY);

    xx = (int) mappedX;
    yy = (int) mappedY;
}

/** \brief Locate the point in the right image
 *
 * \param x const double&
 * \param y const double&
 * \param z const double&
 * \param xx int&
 * \param yy int&
 * \return void
 *
 */
inline void Cube2Cyl::locateRight(const double& x, const double& y, const double& z, int& xx, int& yy)
{
    sizeRatio = diameter / z;

    mappedX = sizeRatio * y;
    mappedY = sizeRatio * x;

    // rotate pi/2
    rotRad(M_PI_2, mappedX, mappedY, sizeRatio);

    // translate
    transDis(diameter, mappedX, mappedY);

    xx = (int) mappedX;
    yy = (int) mappedY;
}

/** \brief Calculate the face which the point is on
 *
 * \param theta const double&
 * \param phi const double&
 * \return void
 *
 */
inline void Cube2Cyl::calCubeFace(const double& theta, const double& phi)
{
    // Looking at the cube from top down
    // FRONT zone
    if (   cmpDoubleEqualLager(theta, -M_PI_4, resCal)
        && cmpDoubleSmaller(theta, M_PI_4, resCal))
    {
        cubeFaceId = CUBE_FRONT;
        normTheta  = theta;
    }
    // LEFT zone
    else if (   cmpDoubleEqualLager(theta, -(M_PI_2+M_PI_4), resCal)
             && cmpDoubleSmaller(theta, -M_PI_4, resCal))
    {
        cubeFaceId = CUBE_LEFT;
        normTheta  = theta + M_PI_2;
    }
    // RIGHT zone
    else if (   cmpDoubleEqualLager(theta, M_PI_4, resCal)
             && cmpDoubleSmaller(theta, M_PI_2+M_PI_4, resCal))
    {
        cubeFaceId = CUBE_RIGHT;
        normTheta  = theta - M_PI_2;
    }
    else
    {
        cubeFaceId = CUBE_BACK;

        if (theta > 0.0)
        {
            normTheta = theta - M_PI;
        }
        else
        {
            normTheta = theta + M_PI;
        }
    }

    // find out which segment the line strikes to
    phiThreshold = atan2(1.0, 1.0/cos(normTheta));

    // in the top segment
    if (phi > phiThreshold)
    {
        cubeFaceId = CUBE_DOWN;
    }
    // in the bottom segment
    else if (phi < -phiThreshold)
    {
        cubeFaceId = CUBE_TOP;
    }
    // in the middle segment
    else
    {
        ;
    }
}

/** \brief Calculate the 3D space coordinates of given 2D coordinates
 *
 * \param i const int&  The coordinate along the x-axis
 * \param j const int&  The coordinate along the y-axis
 * \param x double&
 * \param y double&
 * \param z double&
 * \return void
 *
 */
inline void Cube2Cyl::calXYZ(const int& i, const int& j, double& x, double& y, double& z)
{
    calNormXY(i, j, x, y);
    calThetaAndPhi(x, y, tTheta, tPhi);
    calXyzFromThetaPhi(tTheta, tPhi, x, y, z);

    calCubeFace(tTheta, tPhi);
}

inline void Cube2Cyl::calNormXY(const int& i, const int& j, double& x, double& y)
{
    x = ((2.0*i)/pxPanoSizeH - 1.0) * normFactorX;
    y = ((2.0*j)/pxPanoSizeV - 1.0) * normFactorY;
    // y = 1.0 - (2.0*j)/pxPanoSizeV;
}

inline void Cube2Cyl::calThetaAndPhi(double& x, double& y, double& theta, double& phi)
{
    theta = x * M_PI;
    phi   = y * M_PI_2;
    // for spherical vertical distortion
    // phi = asin(y);
}

inline void Cube2Cyl::calXyzFromThetaPhi(double& theta, double& phi, double& x, double& y, double& z)
{
    x = cos(phi) * cos(theta);
    y = sin(phi);
    z = cos(phi) * sin(theta);
}


/** \brief Test if a == b
 *
 * \param a const double&
 * \param b const double&
 * \param epsilon const double& The precision
 * \return bool
 *
 */
inline bool Cube2Cyl::cmpDoubleEqual(const double &a, const double &b, const double &epsilon)
{
    return (fabs(a - b) < epsilon);
}

/** \brief Test if a < b
 *
 * \param a const double&
 * \param b const double&
 * \param epsilon const double&
 * \return bool
 *
 */
inline bool Cube2Cyl::cmpDoubleSmaller(const double &a, const double &b, const double &epsilon)
{
    return ((a - b) < 0) && (!cmpDoubleEqual(a, b, epsilon));
}

/** \brief Test if a <= b
 *
 * \param a const double&
 * \param b const double&
 * \param epsilon const double&
 * \return bool
 *
 */
inline bool Cube2Cyl::cmpDoubleEqualSmaller(const double &a, const double &b, const double &epsilon)
{
    return ((a - b) < 0) || Cube2Cyl::cmpDoubleEqual(a, b, epsilon);
}


/** \brief Test if a > b
 *
 * \param a const double&
 * \param b const double&
 * \param epsilon const double&
 * \return bool
 *
 */
inline bool Cube2Cyl::Cube2CylcmpDoubleLager(const double &a, const double &b, const double &epsilon)
{
    return ((a - b) > 0) && (!cmpDoubleEqual(a, b, epsilon));
}

/** \brief Test if a >= b
 *
 * \param a const double&
 * \param b const double&
 * \param epsilon const double&
 * \return bool
 *
 */
inline bool Cube2Cyl::cmpDoubleEqualLager(const double &a, const double &b, const double &epsilon)
{
    return ((a - b) > 0) || Cube2Cyl::cmpDoubleEqual(a, b, epsilon);
}


/** \brief The default constructor
 *
 * \param void
 *
 */
Cube2Cyl::Cube2Cyl(void)
:   pxCamV(0),
    pxCamH(0),
    diameter(0.0),
    radPanoV(0.0),
    radPanoH(0.0),
    pxPanoSizeV(0),
    pxPanoSizeH(0),
    cubeFaceId(CUBE_FRONT),
    normTheta(0.0),
    resCal(0.001),
    normFactorX(1.0),
    normFactorY(1.0),
    sizeRatio(1.0),
    mappedX(0),
    mappedY(0),
    tX(0.0),
    tY(0.0),
    tZ(0.0),
    tTheta(0.0),
    tPhi(0.0),
    phiThreshold(0.0)
{

}

/** \brief The default destructor
 *
 * \param void
 *
 */
Cube2Cyl::~Cube2Cyl(void)
{

}

#endif // CUBE2CYL_H_INCLUDED
