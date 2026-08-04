cmake_minimum_required(VERSION 3.18)

project(OpenPPFD
  DESCRIPTION "A radiative-transfer PPFD / spectrum solver for plant factory lighting"
  LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ── Windows (MSVC) 対応 ──────────────────────────────────────────────
# libm は存在しない (MATH_LIB を空にする)。ソースは UTF-8 (日本語コメント)。
# M_PI 等は _USE_MATH_DEFINES が必要。
if(MSVC)
    add_compile_options(/utf-8)
    add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS NOMINMAX)
    set(MATH_LIB "")
else()
    set(MATH_LIB m)
endif()

# OpenMP は任意 (無ければシリアルでビルド : macOS の AppleClang 等)
find_package(OpenMP QUIET)

# output folder
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(oppfd
    src/main.c
    src/input_data.c
    src/spectrum.c
    src/photometry.c
    src/geometry.c
    src/formfactor.c
    src/direct.c
    src/radiosity.c
    src/solve.c
    src/output.c
    src/utils.c
)

target_include_directories(oppfd PRIVATE ${PROJECT_SOURCE_DIR}/include)

if(OpenMP_C_FOUND)
    target_link_libraries(oppfd PRIVATE OpenMP::OpenMP_C)
    message(STATUS "OpenMP が見つかりました。並列ビルドを有効にします。")
else()
    message(STATUS "OpenMP が見つかりません。シリアルでビルドします。")
endif()

target_link_libraries(oppfd PRIVATE ${MATH_LIB})
