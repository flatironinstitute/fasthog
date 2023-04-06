cmake_minimum_required(VERSION 3.0)
project(
  hog
  LANGUAGES C
  )

set(HOG_INCLUDES
  ${PROJECT_SOURCE_DIR}/include
  )

file(GLOB HOG_SOURCES "src/*.c")
add_library(hog SHARED ${HOG_SOURCES})
target_include_directories(hog PUBLIC ${HOG_INCLUDES})
target_compile_options(hog PRIVATE)
