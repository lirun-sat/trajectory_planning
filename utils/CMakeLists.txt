# 查找当前目录下的所有源文件
# 并将名称保存到 DIR_LIB_SRCS 变量

aux_source_directory(. DIR_LIB_SRCS)

# 生成链接库

add_library (utils ${DIR_LIB_SRCS})

target_include_directories(utils PUBLIC ${PROJECT_SOURCE_DIR}/include)

set (LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
