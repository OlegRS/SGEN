# Minimum required version of CMake
cmake_minimum_required(VERSION 3.14)

# Project name and language
project(SGEN LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify the include directories
include_directories(include)

# Collect source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Create a static library
add_library(SGEN STATIC ${SOURCES})
set_target_properties(SGEN PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Add dependencies (optional)
find_package(Armadillo REQUIRED)
target_include_directories(SGEN PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(SGEN PUBLIC ${ARMADILLO_LIBRARIES})

# # Add an executable for the main .cpp file
# add_executable(optimised_moments main/optimised_moments.cpp)

# # Link the library to the executable
# target_link_libraries(optimised_moments PRIVATE SGEN)

# # Specify include directories for the executable
# target_include_directories(optimised_moments PRIVATE include)

# Add Python bindings
# Specify pybind11 directory if not found automatically
# set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "../extern/pybind11/build/")
find_package(pybind11 REQUIRED)

# Collect all binding files
set(BINDING_SOURCES
  bindings/bindings.cpp
  bindings/compartments/Compartment_binding.cpp
  bindings/compartments/Soma_binding.cpp
  bindings/compartments/Spine_binding.cpp
  bindings/compartments/Dendritic_segment_binding.cpp
  bindings/Neuron_binding.cpp
  bindings/engines/Analytic_engine_binding.cpp
  bindings/engines/Morphologic_engine_binding.cpp
  bindings/engines/Gillespie_engine.cpp
)

# Add Python bindings module
pybind11_add_module(_SGEN_Py ${BINDING_SOURCES})

# Link core library to the Python bindings
target_link_libraries(_SGEN_Py PRIVATE SGEN)

# Include directories for Python bindings
target_include_directories(_SGEN_Py PRIVATE ${CMAKE_SOURCE_DIR}/include)



# # Install the shared object file to the site-packages directory
# # Find Python's site-packages directory
# execute_process(COMMAND ${Python3_EXECUTABLE} -c
#     "import site; print(site.getsitepackages()[0])"
#     OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)

# # Define the install target
# install(TARGETS SGEN_Py
#     LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES})
