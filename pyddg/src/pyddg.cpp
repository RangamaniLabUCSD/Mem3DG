
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ddgsolver/ddgsolver.h"

namespace ddgsolver {
namespace py = pybind11;

// Initialize the `pyddg` module
PYBIND11_MODULE(pyddg, pyddg) {
  pyddg.doc() = "Python wrapper around the DDG solver C++ library.";

  pyddg.def("driver", &driver, py::arg("filename"),
            R"delim(
                Run the driver.

                Args:
                    filename (:py:class:`str`): Filename to read.
            
                Returns:
                    :py:class:`int`: success.
            )delim");
};
} // namespace ddgsolver
