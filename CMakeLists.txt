cmake_minimum_required(VERSION 3.16)
project(STDAtmo VERSION 1.0.0)

# Auto-find Qt - Add PrintSupport for QCustomPlot
find_package(Qt6 REQUIRED COMPONENTS Widgets Core Gui PrintSupport)

# CRITICAL FIX: Apply big-obj flag BEFORE adding sources
if(MINGW)
    add_compile_options(-Wa,-mbig-obj)
    add_compile_options(-O1)  # Lower optimization helps
endif()



set(SOURCES
    main.cpp
    stdatmo.cpp
    libs/qcustomplot/qcustomplot.cpp  # Direct path, no variable
)

set(HEADERS
    stdatmo.h
    units.h
    libs/qcustomplot/qcustomplot.h
)

set(UI_FILES
    stdatmo.ui
)

set(RESOURCE_FILES
    resources.qrc
)

# Enable Qt's auto features
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

add_executable(STDAtmo WIN32
    ${SOURCES}
    ${UI_FILES}
    ${HEADERS}
    ${RESOURCE_FILES}
)

# Link Qt libraries - ADD PrintSupport
target_link_libraries(STDAtmo
    Qt6::Widgets 
    Qt6::Core
    Qt6::Gui
    Qt6::PrintSupport  # REQUIRED for QCustomPlot
)

# Include directories
target_include_directories(STDAtmo PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}  # For your headers
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/qcustomplot
)

# Set C++ standard
set_target_properties(STDAtmo PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

# Windows-specific settings
if(WIN32)
    # For MSVC only, not MinGW
    if(MSVC)
        add_compile_options("/bigobj")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS")
    endif()
    
    # Copy Qt DLLs automatically (Qt6 feature)
    set(QT6_COPY_DLLS ON)
endif()
