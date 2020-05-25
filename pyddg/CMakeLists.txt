set(PYDDG_SOURCES "src/pyddg.cpp")

pybind11_add_module(pyddg MODULE ${PYDDG_SOURCES})
target_link_libraries(pyddg PRIVATE ddg_solver)
