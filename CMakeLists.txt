cmake_minimum_required(VERSION 3.0)
project(
  fasthog
  LANGUAGES CXX
  )

set(HOG_INCLUDES
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/extern
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mfma")

file(GLOB HOG_SOURCES "src/*.cpp")
add_library(fasthog SHARED ${HOG_SOURCES})
target_include_directories(fasthog PUBLIC ${HOG_INCLUDES})
target_compile_options(fasthog PRIVATE)

install(TARGETS fasthog)
