project(demo_cpp)

include_directories(
    ${qdbmp_SOURCE_DIR}
    ${libcube2cyl_SOURCE_DIR}
)

set(demo_cpp_SRCS
    main.cpp
)

add_executable(demo_cpp
    ${demo_cpp_SRCS}
)

target_link_libraries(demo_cpp
    qdbmp
)

set(INS_TARGETS
    demo_cpp
)

# Install
install(TARGETS ${INS_TARGETS}
    RUNTIME DESTINATION bin COMPONENT applications
    LIBRARY DESTINATION lib COMPONENT libraries
    ARCHIVE DESTINATION lib COMPONENT libraries
)
