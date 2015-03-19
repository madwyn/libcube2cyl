project(demo_c)

include_directories(
    ${qdbmp_SOURCE_DIR}
    ${libcube2cyl_SOURCE_DIR}
)

set(demo_c_SRCS
    main.c
)

add_executable(demo_c
    ${demo_c_SRCS}
)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_link_libraries(demo_c
        m
        qdbmp
    )
else()
    target_link_libraries(demo_c
        qdbmp
    )
endif()

set(INS_TARGETS
    demo_c
)

# Install
install(TARGETS ${INS_TARGETS}
    RUNTIME DESTINATION bin COMPONENT applications
    LIBRARY DESTINATION lib COMPONENT libraries
    ARCHIVE DESTINATION lib COMPONENT libraries
)
