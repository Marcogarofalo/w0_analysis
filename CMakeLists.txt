cmake_minimum_required(VERSION 3.13)

project(w0 CXX C)

IF(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all")
  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment")
  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -static-libasan -g   ") # -fsanitize=float-divide-by-zero
  set(LDFLAGS "${LDFLAGS}  -fsanitize=address")
  set(CMAKE_LDFLAGS_FLAGS " ${CMAKE_LDFLAGS_FLAGS} ${LDFLAGS} ")
ENDIF()

# ############################################
# load analysis_program
# ############################################
find_package(analysis_program)
include_directories(${ANALYSIS_PROGRAM_INCLUDE_DIR})
MESSAGE("ANALYSIS_PROGRAM_INCLUDE_DIR   ${ANALYSIS_PROGRAM_INCLUDE_DIR}")
MESSAGE("ANALYSIS_PROGRAM_LIBRARIES: ${ANALYSIS_PROGRAM_LIBRARIES}")
MESSAGE("ANALYSIS_PROGRAM_LIBRARY_DIR: ${ANALYSIS_PROGRAM_LIBRARY_DIR}")

IF(NOT analysis_program_FOUND)
  MESSAGE("find_package faild to find yaml-cpp: you need to do it manually")
  MESSAGE("hint:  analysis_program_FOUND_DIR or CMAKE_PREFIX_PATH need to be global path or relative to the source dir")
  MESSAGE(FATAL_ERROR)
ENDIF()

find_package(Eigen3 REQUIRED) # we need eigen
MESSAGE("IF the analysis_program was compiled with OpenMP then you also need to load it")
MESSAGE("uncomment the lines after this message")

# FIND_PACKAGE( OpenMP)
# if(OPENMP_FOUND)
# message("OPENMP FOUND")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()
# include_directories(SYSTEM ${OpenMP_INCLUDE_PATH})
# ############################################
# add a local library for this project
# ############################################
add_library(
  functions_w0 STATIC
  functions_w0.hpp
  functions_w0.cpp
)
target_link_libraries(functions_w0 PUBLIC ${ANALYSIS_PROGRAM_LIBRARIES})
target_include_directories(functions_w0 PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(w0 w0.cpp)
target_link_libraries(w0 PUBLIC functions_w0)
IF(WITH_ARB)
  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -DWITH_ARB   ")
  target_link_libraries(w0 PUBLIC flint-arb flint)
ENDIF()
target_link_libraries(w0 PUBLIC ${ANALYSIS_PROGRAM_LIBRARIES})

add_executable(w0_rew w0_rew.cpp)
target_link_libraries(w0_rew PUBLIC functions_w0)
IF(WITH_ARB)
  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -DWITH_ARB   ")
  target_link_libraries(w0_rew PUBLIC flint-arb flint)
ENDIF()
target_link_libraries(w0_rew PUBLIC ${ANALYSIS_PROGRAM_LIBRARIES})

macro(add_target_with_libs TARGETNAME)
  # Make the executable
  add_executable(${TARGETNAME} ${ARGN})

  target_link_libraries(${TARGETNAME} PUBLIC functions_w0)
  IF(WITH_ARB)
    set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -DWITH_ARB   ")
    target_link_libraries(${TARGETNAME} PUBLIC flint-arb flint)
  ENDIF()
  target_link_libraries(${TARGETNAME} PUBLIC ${ANALYSIS_PROGRAM_LIBRARIES})

endmacro()

add_target_with_libs(pion_rew pion_rew.cpp)
add_target_with_libs(pion pion.cpp)
add_target_with_libs(fit_all_rew_der fit_all_rew_der.cpp)

add_target_with_libs(create_fake_jack create_fake_jack.cpp)
add_target_with_libs(fit_all_beta_rew_der fit_all_beta_rew_der.cpp)

add_target_with_libs(pion_rew_ratio pion_rew_ratio.cpp)

add_target_with_libs(fit_all_beta_w0rew_der fit_all_beta_w0rew_der.cpp)
