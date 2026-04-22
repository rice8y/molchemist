cmake_minimum_required(VERSION 3.9)

project(molchemist_smiles_plugin LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-fno-exceptions)

set(COORDGEN_RIGOROUS_BUILD OFF CACHE BOOL "" FORCE)
set(COORDGEN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(COORDGEN_BUILD_EXAMPLE OFF CACHE BOOL "" FORCE)
set(COORDGEN_USE_MAEPARSER OFF CACHE BOOL "" FORCE)
set(COORDGEN_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/vendor/coordgenlibs ${CMAKE_CURRENT_BINARY_DIR}/vendor/coordgenlibs)

add_library(${PROJECT_NAME}-object OBJECT src/coordgen_smiles.cpp)
target_include_directories(${PROJECT_NAME}-object
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor/coordgenlibs
)
target_link_libraries(${PROJECT_NAME}-object PRIVATE coordgen)

set(EMCC_FLAGS
  -Oz
  --no-entry
  -fno-exceptions
  -sFILESYSTEM=0
  -sASSERTIONS=0
  -sEXPORT_KEEPALIVE=1
  -sINITIAL_MEMORY=128mb
  -sSTACK_SIZE=1mb
  -sALLOW_MEMORY_GROWTH=0
  -Wl,--allow-undefined
)

add_custom_target(${PROJECT_NAME}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/dist
  COMMAND emcc
    -o ${CMAKE_CURRENT_BINARY_DIR}/dist/${PROJECT_NAME}.wasm
    ${EMCC_FLAGS}
    $<TARGET_OBJECTS:${PROJECT_NAME}-object>
    $<TARGET_FILE:coordgen>
  COMMAND wasi-stub
    ${CMAKE_CURRENT_BINARY_DIR}/dist/${PROJECT_NAME}.wasm
    --stub-function env:__syscall_unlinkat,env:__syscall_faccessat
    -o ${CMAKE_CURRENT_BINARY_DIR}/dist/${PROJECT_NAME}.stub.wasm
  COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_BINARY_DIR}/dist/${PROJECT_NAME}.stub.wasm
    ${CMAKE_CURRENT_BINARY_DIR}/dist/${PROJECT_NAME}.wasm
  DEPENDS
    ${PROJECT_NAME}-object
)
