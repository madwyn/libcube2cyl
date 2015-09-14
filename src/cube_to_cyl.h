#ifndef CUBE_TO_CYL_H_INCLUDED
#define CUBE_TO_CYL_H_INCLUDED

#include <stdlib.h>     // required for malloc/free
#include <math.h>

#ifndef M_PI
#define M_PI        3.14159265358979323846
#define M_PI_2      1.57079632679489661923
#define M_PI_4      0.78539816339744830962
#endif

// sphere radius
#define S_RADIUS    1.0

#ifdef _MSC_VER
    #ifndef inline
        #define inline __inline
    #endif
#endif

// defines the faces of a cube
enum cc_face {
    CC_TOP = 0,
    CC_LEFT,
    CC_FRONT,
    CC_RIGHT,
    CC_BACK,
    CC_DOWN,
    CC_FACE_NUM
};

// the cubic image coordinates
struct cc_coord {
    enum cc_face face;    // the face of the cube
    double       x;       // the x coordinate
    double       y;       // the y coordinate
};

struct cc_context {
    //-------- input information
    unsigned int px_cam_h;  /**< The horizontal pixel of a camera */
    unsigned int px_cam_v;  /**< The vertical pixel of a camera */

    double radius;         /**< The radius of the sphere for projection */

    double rd_pano_h;       /**< The horizontal view portion */
    double rd_pano_v;       /**< The vertical view portion */

    //-------- output information
    unsigned int px_pano_h; /**< The horizontal pixels of the panorama */
    unsigned int px_pano_v; /**< The vertical pixels of the panorama */

    //-------- to access the pixel
    enum cc_face cube_face; /**< The cube face to be read */
    double       cube_x;    /**< The x coordinate mapped on the cube face */
    double       cube_y;    /**< The y coordinate mapped on the cube face */

    //-------- The temp variables
    double theta_norm;      /**< The normalised theta */
    double res_cal;         /**< The resolution used for calculation */
    double norm_factor_x;   /**< The normalisation factor for x */
    double norm_factor_y;   /**< The normalisation factor for y */

    double size_ratio;      /**< The size ratio of the mapped x and the actual radius */

    double t_x;             /**< x coordinate in 3D space */
    double t_y;             /**< y coordinate in 3D space */
    double t_z;             /**< z coordinate in 3D space */

    double t_theta;         /**< The radian horizontally */
    double t_phi;           /**< The radian vertically */

    double phi_threshold;   /**< The threshold of phi, it separates the top, middle and down
                                 of the cube, which are the edges of the top and down surface */

    struct cc_coord *map;   /**< The map from panorama coordinates to cubic coordinates */
};

/** \brief Test if a == b
 *
 * \param a       const double*
 * \param b       const double*
 * \param epsilon const double* The precision
 * \return int
 *
 */
static inline int
__cc_db_equal(const double *a, const double *b, const double *epsilon) {
    return (fabs(*a - *b) < *epsilon);
}


/** \brief Test if a < b
 *
 * \param a const double*
 * \param b const double*
 * \param epsilon const double*
 * \return int
 *
 */
static inline int
__cc_db_smaller(const double *a, const double *b, const double *epsilon) {
    return ((*a - *b) < 0.0) && (!__cc_db_equal(a, b, epsilon));
}

/** \brief Test if a >= b
 *
 * \param a const double*
 * \param b const double*
 * \param epsilon const double*
 * \return int
 *
 */
static inline int
__cc_db_larger(const double *a, const double *b, const double *epsilon) {
    return ((*a - *b) > 0.0) && (!__cc_db_equal(a, b, epsilon));
}

/** \brief rd in range [small, large)
 *
 * \param rd    const double*
 * \param small const double*
 * \param large const double*
 * \param res   const double*
 * \return int
 *
 */
static inline int
__cc_db_in_range_c_o(const double rd, const double small, const double large, const double res) {
   return    !__cc_db_smaller(&rd, &small, &res)    /* rd >= small */
          &&  __cc_db_smaller(&rd, &large, &res);
}

/** \brief rd in range (small, large]
 *
 * \param rd    const double*
 * \param small const double*
 * \param large const double*
 * \param res   const double*
 * \return int
 *
 */
static inline int
__cc_db_in_range_o_c(const double rd, const double small, const double large, const double res) {
   return     __cc_db_larger(&rd, &small, &res)
          && !__cc_db_larger(&rd, &large, &res);    /* rd <= small */
}

/** \brief
 *
 * \param ctx       struct cc_context*    The context for calculation
 * \param px_in_w   const unsigned int    The input square image width/height
 * \param rd_view_h const       double    The radian of the panorama view portion horizontally, range [0.01, 2*PI]
 * \param rd_view_v const       double    The radian of the panorama view portion vertically, range [0.01, PI]
 * \param px_pano_h const unsigned int    The desired panorama size horizontally
 * \param px_pano_v const unsigned int    The desired panorama size vertically
 * \return void
 *
 */
static inline void
cc_init_full(struct  cc_context *ctx,       const unsigned int px_in_w,
             const       double  rd_view_h, const       double rd_view_v,
             const unsigned int  px_pano_h, const unsigned int px_pano_v) {

    // check parameters
    if (   (px_in_w   == 0)
        || (px_pano_h == 0)
        || (px_pano_v == 0)) {
        return;
    }

    // check the view portion
    if (   !__cc_db_in_range_o_c(rd_view_v, 0.001, M_PI,     0.001)
        || !__cc_db_in_range_o_c(rd_view_h, 0.001, M_PI*2.0, 0.001)) {
        return;
    }

    ctx->px_cam_h  = px_in_w;
    ctx->px_cam_v  = px_in_w;

    ctx->rd_pano_h = rd_view_h;
    ctx->rd_pano_v = rd_view_v;

    ctx->px_pano_h = px_pano_h;
    ctx->px_pano_v = px_pano_v;

    ctx->radius = ((double)px_in_w) / 2.0;

    // the actual calculation resolution is 10 times bigger than the texture resolution
    ctx->res_cal = M_PI_4 / (double)px_in_w / 10.0;

    // the normalisation factors
    ctx->norm_factor_x = rd_view_h / (M_PI * 2);
    ctx->norm_factor_y = rd_view_v /  M_PI;

    ctx->size_ratio = 0.0;

    ctx->t_x = 0.0;
    ctx->t_y = 0.0;
    ctx->t_z = 0.0;

    ctx->t_theta = 0.0;
    ctx->t_phi   = 0.0;

    ctx->phi_threshold = 0.0;

    ctx->map = NULL;
}

/** \brief Initialise the context
 *
 * \param ctx       struct cc_context*    The context for calculation
 * \param px_in_w   const unsigned int    The input square image width/height
 * \param rd_view_h const       double    The radian of the panorama view portion horizontally, range (0.001, 2*PI]
 * \param rd_view_v const       double    The radian of the panorama view portion vertically,   range (0.001, PI]
 * \return void
 *
 */
static inline void
cc_init(struct  cc_context *ctx,       const unsigned int px_in_w,
        const       double  rd_view_h, const       double rd_view_v) {

    unsigned int px_pano_h = (unsigned int)((rd_view_h/M_PI_2) * (double)px_in_w);
    unsigned int px_pano_v = (unsigned int)((rd_view_v/M_PI_2) * (double)px_in_w);
    cc_init_full(ctx, px_in_w, rd_view_h, rd_view_v, px_pano_h, px_pano_v);
}


/** \brief Visit the coordinate
 *
 * \param *ctx const struct   cc_context
 * \param x    const unsigned int
 * \param y    const unsigned int
 * \return const struct cc_coord*
 *
 */
static inline const struct cc_coord*
cc_get_coord(const struct cc_context *ctx, const unsigned int x, const unsigned int y) {
    return &(ctx->map[x*ctx->px_pano_v + y]);
}

/** \brief Release the context
 *
 * \param ctx struct cc_context*
 * \return void
 *
 */
static inline void
cc_close(struct cc_context *ctx) {
    free(ctx->map);
}

static inline void
__cc_cal_norm_xy(struct cc_context *ctx, const int x, const int y) {
    ctx->t_x = ((2.0 * x) / ctx->px_pano_h - 1.0) * ctx->norm_factor_x;
    ctx->t_y = ((2.0 * y) / ctx->px_pano_v - 1.0) * ctx->norm_factor_y;
    // ctx->t_y = 1.0 - (2.0*j)/ctx->px_pano_v;
}

static inline void
__cc_cal_theta_phi(struct cc_context *ctx) {
    ctx->t_theta = ctx->t_x * M_PI;
    ctx->t_phi   = ctx->t_y * M_PI_2;
    // for spherical vertical distortion
//    ctx->t_phi   = asin(ctx->t_y);
}

static inline void
__cc_cal_xyz_from_theta_phi(struct cc_context *ctx) {
    ctx->t_x = cos(ctx->t_phi) * cos(ctx->t_theta);
    ctx->t_y = sin(ctx->t_phi);
    ctx->t_z = cos(ctx->t_phi) * sin(ctx->t_theta);
}

static inline void
__cc_cal_cube_face(struct cc_context *ctx) {
    // Looking at the cube from top down
    // FRONT zone
    if (__cc_db_in_range_c_o(ctx->t_theta, -M_PI_4, M_PI_4, ctx->res_cal)) {
        ctx->cube_face  = CC_FRONT;
        ctx->theta_norm = ctx->t_theta;
    }
    // LEFT zone
    else if (__cc_db_in_range_c_o(ctx->t_theta, -(M_PI_2 + M_PI_4), -M_PI_4, ctx->res_cal)) {
        ctx->cube_face  = CC_LEFT;
        ctx->theta_norm = ctx->t_theta + M_PI_2;
    }
    // RIGHT zone
    else if (__cc_db_in_range_c_o(ctx->t_theta, M_PI_4, M_PI_2 + M_PI_4, ctx->res_cal)) {
        ctx->cube_face  = CC_RIGHT;
        ctx->theta_norm = ctx->t_theta - M_PI_2;
    }
    else {
        ctx->cube_face  = CC_BACK;
        ctx->theta_norm = ctx->t_theta + ((ctx->t_theta > 0.0) ? -M_PI : M_PI);
    }

    // find out which segment the line strikes to
    ctx->phi_threshold = atan2(S_RADIUS, S_RADIUS / cos(ctx->theta_norm));

    // in the top segment
    if (ctx->t_phi > ctx->phi_threshold) {
        ctx->cube_face = CC_DOWN;
    }
    // in the bottom segment
    else if (ctx->t_phi < -ctx->phi_threshold) {
        ctx->cube_face = CC_TOP;
    }
    // in the middle segment
    else {
        ;
    }
}

/** \brief Calculate the point(x, y) mapped on the sphere
 *
 * \param ctx struct cc_context*
 * \param x   const unsigned int
 * \param y   const unsigned int
 * \return void
 *
 */
static inline void
__cc_cal_sphere_xyz(struct cc_context *ctx, const unsigned int x, const unsigned int y) {
    // normalise the x, y
    __cc_cal_norm_xy(ctx, x, y);

    // calculate theta and phi
    __cc_cal_theta_phi(ctx);

    // calculate the sphere xyz from theta and phi
    __cc_cal_xyz_from_theta_phi(ctx);

    // calculate the cube face
    __cc_cal_cube_face(ctx);
}

/** \brief Rotate the point for a given radian
 *
 * \param x double*
 * \param y double*
 * \param rad const double
 * \param temp double*
 * \return void
 *
 */
static inline void
__cc_rot_rad(double *x, double *y, const double rad, double *temp) {
    *temp = *x;
    *x = *x    * cos(rad) - *y * sin(rad);
    *y = *temp * sin(rad) + *y * cos(rad);
}

/** \brief Translate the point for a given distance in both x and y coordinates
 *
 * \param x double*
 * \param y double*
 * \param dis const double*
 * \return void
 *
 */
static inline void
__cc_trans_dis(double *x, double *y, const double *dis) {
    *x += *dis;
    *y += *dis;
}

/** \brief Locate the projected coordinates from sphere to 2D plane.
 *
 * \param ctx  struct cc_context*
 * \param axis const  double*       The coordinate on the axis
 * \param px   const  double*
 * \param py   const  double*
 * \param rad  const  double
 * \return void
 *
 */
static inline void
__cc_loc(struct cc_context *ctx, const double axis, const double px, const double py, const double rad) {
    ctx->size_ratio = ctx->radius / axis;

    ctx->cube_x = ctx->size_ratio * px;
    ctx->cube_y = ctx->size_ratio * py;

    // rotate pi
    __cc_rot_rad(&(ctx->cube_x), &(ctx->cube_y), rad, &(ctx->size_ratio));

    // translate
    __cc_trans_dis(&(ctx->cube_x), &(ctx->cube_y), &(ctx->radius));
}

/** \brief Calculate the projected coordinates on the cube of x y
 *
 * \param ctx struct cc_context*
 * \param x   const unsignedint
 * \param y   const unsignedint
 * \return void
 *
 */
static inline void
__cc_cal_cube_xy(struct cc_context *ctx, const unsigned int x, const unsigned int y) {
    // calculate the sphere coordinates
    __cc_cal_sphere_xyz(ctx, x, y);

    switch (ctx->cube_face)
    {
        case CC_TOP:
        {
            __cc_loc(ctx, ctx->t_y, ctx->t_z, ctx->t_x, M_PI);
            break;
        }
        case CC_DOWN:
        {
            __cc_loc(ctx, ctx->t_y, ctx->t_x, ctx->t_z, -M_PI_2);
            break;
        }
        case CC_LEFT:
        {
            __cc_loc(ctx, ctx->t_z, ctx->t_x, ctx->t_y, M_PI);
            break;
        }
        case CC_RIGHT:
        {
            __cc_loc(ctx, ctx->t_z, ctx->t_y, ctx->t_x, M_PI_2);
            break;
        }
        case CC_FRONT:
        {
            __cc_loc(ctx, ctx->t_x, ctx->t_z, ctx->t_y, 0.0);
            break;
        }
        case CC_BACK:
        {
            __cc_loc(ctx, ctx->t_x, ctx->t_y, ctx->t_z, -M_PI_2);
            break;
        }
        default:
        {
            break;
        }
    }
}

/** \brief Generate the mapping from panorama to cubic images
 *
 * \param ctx struct cc_context*
 * \return void
 *
 */
static inline void
cc_gen_map(struct cc_context *ctx) {
    unsigned int x, y, pos;

    if (NULL != ctx->map) {
        free(ctx->map);
    }

    ctx->map = (struct cc_coord*)malloc(ctx->px_pano_v * ctx->px_pano_h * sizeof(struct cc_coord));

    for (x = 0, y = 0, pos = 0; x < ctx->px_pano_h; ++x) {
        for (y = 0; y < ctx->px_pano_v; ++y) {
            __cc_cal_cube_xy(ctx, x, y);

            ctx->map[pos  ].face = ctx->cube_face;
            ctx->map[pos  ].x    = ctx->cube_x;
            ctx->map[pos++].y    = ctx->cube_y;
        }
    }
}

#endif // CUBE_TO_CYL_H_INCLUDED
