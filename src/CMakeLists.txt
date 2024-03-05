Primary navigation
Homepage

 Search or go to…
Project
P
PI_Luca_Killian

Pinned
Issues
0
Merge requests
0

Manage

Plan

Code
Merge requests
0
Repository
Branches
Commits
Tags
Repository graph
Compare revisions
Snippets
Locked files

Build

Secure

Deploy

Operate

Monitor

Analyze

Settings
Santoro Luca
PI_Luca_Killian
Repository
Killian's avatar
test paraview
Kollasch Killian authored 5 minutes ago
ab7aa6ce
 Code owners
Assign users and groups as approvers for specific file changes. Learn more.
pi_luca_killian
CMakeLists.txt
CMakeLists.txt
1.47 KiB
CMAKE_MINIMUM_REQUIRED(VERSION 3.18)
PROJECT(SOLVER CXX)
INCLUDE(CMakePrintHelpers)
# ------------------------------------------------------------------------------
# Find libraries and setup compiler
# ------------------------------------------------------------------------------
# change "Debug" => "Release" for a faster executable (but less runtime checks)!
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
ENDIF()
# enable C++17
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
IF(APPLE)
    # on macOS, do not give priority to frameworks/apps
    SET(CMAKE_FIND_APPBUNDLE LAST)
    SET(CMAKE_FIND_FRAMEWORK LAST)
ENDIF()
# increase warning level
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-Wall" WALL)
IF(WALL)
    ADD_COMPILE_OPTIONS(-Wall)
ENDIF()
# find OpenMP
FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    MESSAGE (STATUS "OpenMP enabled (CXX_FLAGS=${OpenMP_CXX_FLAGS})")
ELSE()
    MESSAGE (STATUS "OpenMP disabled.")
ENDIF()
# Find zlib for compressed output files
FIND_PACKAGE(ZLIB)
# ------------------------------------------------------------------------------

SET(SRCS
    main.cpp
    initialize.cpp
    generate_particle.cpp
    export.cpp
    # Add more files as needed
)
ADD_EXECUTABLE(solver ${SRCS})
IF(OpenMP_CXX_FOUND)
    TARGET_LINK_LIBRARIES(solver PUBLIC OpenMP::OpenMP_CXX)
ENDIF()
IF(ZLIB_FOUND)
    TARGET_COMPILE_DEFINITIONS(solver PUBLIC -DUSE_ZLIB)
    TARGET_LINK_LIBRARIES(solver PUBLIC ZLIB::ZLIB)
ENDIF()
