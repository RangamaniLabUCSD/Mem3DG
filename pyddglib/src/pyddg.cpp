
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ddgsolver/ddgsolver.h"

namespace ddgsolver {
namespace py = pybind11;

// Initialize the `pyddg` module
PYBIND11_MODULE(pyddglib, pyddg) {
  pyddg.doc() = "Python wrapper around the DDG solver C++ library.";

  pyddg.def("driver", &driver, " a driver function",
    py::arg("run") = "visualization", py::arg("option") = "sphere",
    py::arg("nSub"), py::arg("H0"), py::arg("Vt"), 
    py::arg("h"), py::arg("T"), py::arg("eps"));
  //pyddg.def("driver", &driver," run the driver", py::arg("filename"));
            //R"delim(
            //    Run the driver.

            //    Args:
            //        filename (:py:class:`str`): Filename to read.
            //
            //    Returns:
            //        :py:class:`int`: success.
            //)delim");
};
} // namespace ddgsolver
