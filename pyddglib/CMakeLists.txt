set(PYDDG_SOURCES "src/pyddg.cpp")

pybind11_add_module(pyddglib MODULE ${PYDDG_SOURCES})
#add_library(pyddglib MODULE ${PYDDG_SOURCES})
target_link_libraries(pyddglib PRIVATE ddg_solver)
#target_link_libraries(pyddglib PRIVATE pybind11::module)
