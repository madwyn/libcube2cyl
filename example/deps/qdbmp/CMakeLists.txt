project(qdbmp)

set(qdbmp_SRCS
    qdbmp.c
)

set(qdbmp_HDRS
    qdbmp.h
)

add_library(qdbmp ${qdbmp_SRCS})

set(INS_INCLUDE_DIR
    "/qdbmp"
)

set(INS_TARGETS
    qdbmp
)

# Install
install(TARGETS ${INS_TARGETS}
    RUNTIME DESTINATION bin COMPONENT applications
    LIBRARY DESTINATION lib COMPONENT libraries
    ARCHIVE DESTINATION lib COMPONENT libraries
)

install(FILES ${qdbmp_HDRS} DESTINATION "include${INS_INCLUDE_DIR}" COMPONENT headers)
