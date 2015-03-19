project(libcube2cyl)

set(libcube2cyl_HDRS
    Cube2Cyl.hpp
    cube_to_cyl.h
)

set(INS_INCLUDE_DIR
    "/libcube2cyl"
)

install(FILES ${libcube2cyl_HDRS} DESTINATION "include${INS_INCLUDE_DIR}" COMPONENT headers)
