option(BODY_MODEL_USE_FLOAT "BODY_MODEL_USE_FLOAT" ON)
if(BODY_MODEL_USE_FLOAT)
    add_compile_definitions(BODY_MODEL_USE_FLOAT)
endif()

option(BODY_MODEL_USE_DOUBLE "BODY_MODEL_USE_DOUBLE" ON)
if(BODY_MODEL_USE_DOUBLE)
    add_compile_definitions(BODY_MODEL_USE_DOUBLE)
endif()


option(GRID_EL_USE_FLOAT "GRID_EL_USE_FLOAT" ON)
if(GRID_EL_USE_FLOAT)
    add_compile_definitions(GRID_EL_USE_FLOAT)
endif()

option(GRID_EL_USE_DOUBLE "GRID_EL_USE_DOUBLE" ON)
if(GRID_EL_USE_DOUBLE)
    add_compile_definitions(GRID_EL_USE_DOUBLE)
endif()

option(GRID_EL_USE_CUU "GRID_EL_USE_CUU" ON)
if(GRID_EL_USE_CUU)
    add_compile_definitions(GRID_EL_USE_CUU)
endif()

option(GRID_EL_USE_CVV "GRID_EL_USE_CVV" ON)
if(GRID_EL_USE_CVV)
    add_compile_definitions(GRID_EL_USE_CVV)
endif()

option(DISTRIB_USE_FLOAT "DISTRIB_USE_FLOAT" ON)
option(DISTRIB_USE_DOUBLE "DISTRIB_USE_DOUBLE" ON)
if(DISTRIB_USE_FLOAT)
    add_compile_definitions(DISTRIB_USE_FLOAT)
endif()
if(DISTRIB_USE_DOUBLE)
    add_compile_definitions(DISTRIB_USE_DOUBLE)
endif()