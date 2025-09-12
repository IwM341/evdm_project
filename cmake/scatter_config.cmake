
set(P_SOURCES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/evdm_python)

set(SUPPORTED_ALGOLS 0 1 2)
set(SUPPORTED_TYPES float double)
set(SUPPORTED_NMK_TYPES "size_t" "std::vector<size_t>")

set(SUPPORTED_GRIDTYPES "GridCUU" "GridCVV")
set(SUPPORTED_THERMGENS "Naive" "Full" "Soft8" "NoTherm" "Soft8_Treshold")

set(COUNTER 0)

set(COMBINED_FILE ${P_SOURCES_DIR}/scatter_impl/scatter_impl_defs.hpp)

file(WRITE ${COMBINED_FILE} "// Combined template specializations\n")
file(APPEND ${COMBINED_FILE} "// Generated automatically - do not edit\n\n")
file(APPEND ${COMBINED_FILE} "#pragma once\n\n")
file(APPEND ${COMBINED_FILE} "#include <evdm/core/core_dynamics_scatter.hpp>\n\n")
file(APPEND ${COMBINED_FILE} "#include <evdm/utils/prng.hpp>\n\n")
file(APPEND ${COMBINED_FILE} "#include <evdm/measure.hpp>\n\n")

foreach(ALGOL_N IN LISTS SUPPORTED_ALGOLS)
foreach(VAL_T IN LISTS SUPPORTED_TYPES)
foreach(BODY_VAL_T IN LISTS SUPPORTED_TYPES)
foreach(GRID_VAL_T IN LISTS SUPPORTED_TYPES)
foreach(GRID_TYPE IN LISTS SUPPORTED_GRIDTYPES)
foreach(NMK_T IN LISTS SUPPORTED_NMK_TYPES)
foreach(THERMGEN_TYPE IN LISTS SUPPORTED_THERMGENS)
    configure_file(
        ${P_SOURCES_DIR}/scatter_impl/scatter_impl.cpp.in
        ${P_SOURCES_DIR}/scatter_impl/cppimpl/scatter_impl_${COUNTER}.cpp
        @ONLY
    )

    file(READ ${P_SOURCES_DIR}/scatter_impl/scatter_impl_defs.hpp.in TEMPLATE_CONTENT)
    string(REPLACE "@ALGOL_N@" "${ALGOL_N}" MCONTENT "${TEMPLATE_CONTENT}")
    string(REPLACE "@VAL_T@" "${VAL_T}" MCONTENT "${MCONTENT}")
    string(REPLACE "@BODY_VAL_T@" "${BODY_VAL_T}" MCONTENT "${MCONTENT}")
    string(REPLACE "@GRID_VAL_T@" "${GRID_VAL_T}" MCONTENT "${MCONTENT}")
    string(REPLACE "@GRID_TYPE@" "${GRID_TYPE}" MCONTENT "${MCONTENT}")
    string(REPLACE "@NMK_T@" "${NMK_T}" MCONTENT "${MCONTENT}")
    string(REPLACE "@THERMGEN_TYPE@" "${THERMGEN_TYPE}" MCONTENT "${MCONTENT}")
    

    file(APPEND ${COMBINED_FILE} "${MCONTENT}")
    file(APPEND ${COMBINED_FILE} "\n// ======================\n\n")
math(EXPR COUNTER "${COUNTER} + 1")
endforeach()
endforeach()
endforeach()
endforeach()
endforeach()
endforeach()
endforeach()


