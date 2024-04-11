#
# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
#
# Copyright 2024- The Mem3DG Authors
# and the project initiators Cuncheng Zhu, Christopher T. Lee, and
# Padmini Rangamani.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Please help us support Mem3DG development by citing the research
# papers on the package. Check out https://github.com/RangamaniLabUCSD/Mem3DG/
# for more information.

find_package(Python QUIET COMPONENTS Interpreter)
if(Python_Interpreter_FOUND)
  execute_process(
    COMMAND
      ${Python_EXECUTABLE} -c
      "from importlib import util;print(util.find_spec('pytest') is not None)"
    OUTPUT_VARIABLE PYTEST_IMPORT_STATUS
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif(Python_Interpreter_FOUND)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(pytest DEFAULT_MSG PYTEST_IMPORT_STATUS)
