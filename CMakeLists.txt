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

file(GLOB HOG_SOURCES "src/hog.cpp")
if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
  add_library(hog_sse4 OBJECT src/magnitude_orientation.cpp)
  target_compile_options(hog_sse4 PRIVATE -msse4 -DVEC_T=Vec2d)
  target_include_directories(hog_sse4 PRIVATE ${HOG_INCLUDES})

  add_library(hog_avx2 OBJECT src/magnitude_orientation.cpp)
  target_compile_options(hog_avx2 PUBLIC -mavx2 -mfma -DVEC_T=Vec4d)
  target_include_directories(hog_avx2 PRIVATE ${HOG_INCLUDES})

  add_library(hog_avx512 OBJECT src/magnitude_orientation.cpp)
  target_compile_options(hog_avx512 PUBLIC -mavx512f -mfma -DVEC_T=Vec8d)
  target_include_directories(hog_avx512 PRIVATE ${HOG_INCLUDES})

  add_library(fasthog SHARED ${HOG_SOURCES} $<TARGET_OBJECTS:hog_sse4> $<TARGET_OBJECTS:hog_avx2> $<TARGET_OBJECTS:hog_avx512>)
else()
  add_library(fasthog SHARED ${HOG_SOURCES})
endif()

target_include_directories(fasthog PUBLIC ${HOG_INCLUDES})
target_compile_options(fasthog PRIVATE)

install(TARGETS fasthog)
