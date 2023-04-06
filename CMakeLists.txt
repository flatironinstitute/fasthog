cmake_minimum_required(VERSION 3.0)
project(
  hog
  LANGUAGES CXX
  )

set(HOG_INCLUDES
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/extern
)

set(CMAKE_CXX_STANDARD 17)

file(GLOB HOG_SOURCES "src/*.cpp")
add_library(hog SHARED ${HOG_SOURCES})
target_include_directories(hog PUBLIC ${HOG_INCLUDES})
target_compile_options(hog PRIVATE)
