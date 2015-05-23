set(TARGET mhdc3dl)
set(PY_MOD_TARGET mhdc3dl_python)

find_package (PythonLibs)
if(PYTHONLIBS_FOUND)
  include_directories(${PYTHON_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${PYTHON_LIBRARIES})
endif(PYTHONLIBS_FOUND)

find_package (Cat)
if (CAT_FOUND)
  include_directories(${CAT_INCLUDE_DIR})
  set(LIBS ${LIBS} ${CAT_LIBRARIES})
endif (CAT_FOUND)

find_package (FFTW)
if (FFTW_FOUND)
  include_directories(${FFTW_INCLUDE_DIR})
  set(LIBS ${LIBS} ${FFTW_LIBRARIES})
endif (FFTW_FOUND)

find_package (Goops)
if (GOOPS_FOUND)
  include_directories(${GOOPS_INCLUDE_DIR})
  set(LIBS ${LIBS} ${GOOPS_LIBRARIES})
endif (GOOPS_FOUND)

find_package (Godess)
if (GODESS_FOUND)
  include_directories(${GODESS_INCLUDE_DIR})
  set(LIBS ${LIBS} ${GODESS_LIBRARIES})
endif (GODESS_FOUND)

find_package (Lass)
if (LASS_FOUND)
  include_directories(${LASS_INCLUDE_DIR})
  set(LIBS ${LIBS} ${LASS_LIBRARIES})
endif (LASS_FOUND)

message("LASS: " ${LIBS})

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")

include_directories(${SRC_DIR})

set(HEADERS

"${SRC_DIR}/basic.h"
"${SRC_DIR}/block_vector.h"
"${SRC_DIR}/gen_random.h"
"${SRC_DIR}/globals.h"
"${SRC_DIR}/input.h"
"${SRC_DIR}/io.h"
"${SRC_DIR}/linops.h"
"${SRC_DIR}/lss.h"
"${SRC_DIR}/mhdc3dl.h"
"${SRC_DIR}/spectral.h"
"${SRC_DIR}/sss.h"

)

set(TEMPLATE_IMPLEMENTATIONS

"${SRC_DIR}/io.C"

)

set(SOURCES

"${SRC_DIR}/basic.C"
"${SRC_DIR}/gen_random.C"
"${SRC_DIR}/globals.C"
"${SRC_DIR}/input.C"
"${SRC_DIR}/linops.C"
"${SRC_DIR}/lss.C"
"${SRC_DIR}/spectral.C"
"${SRC_DIR}/sss.C"

)

set(MAIN_SRC

"${SRC_DIR}/mhdc3dl.C"

)

set(PY_MOD_SRC

"${SRC_DIR}/mhdc3dl_python.C"

)

add_executable(${TARGET} ${SOURCES} ${MAIN_SRC})
target_link_libraries(${TARGET} ${LIBS} )

add_library(${PY_MOD_TARGET} MODULE ${SOURCES} ${PY_MOD_SRC})
set_target_properties(${PY_MOD_TARGET} PROPERTIES PREFIX "")
target_link_libraries(${PY_MOD_TARGET} ${LIBS} )

install (TARGETS ${TARGET} DESTINATION bin)
install (TARGETS ${PY_MOD_TARGET} DESTINATION lib)


