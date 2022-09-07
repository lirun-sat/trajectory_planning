

set( GK_PUBLIC_HEADER_DIR ${PROJECT_SOURCE_DIR}/include)

add_library (opengjk 
STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/openGJK.c
    ${GK_PUBLIC_HEADER_DIR}/openGJK.h
    )

target_include_directories(opengjk 
PUBLIC 
	# ${PROJECT_SOURCE_DIR}/include
	${CMAKE_SOURCE_DIR}/include
	${CMAKE_BINARY_DIR}
)

target_link_libraries(opengjk PRIVATE m)

set (LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
