add_library(libpoly OBJECT callback.cc gurobi-Dmax.cc)

add_executable(master-g master-gurobi.cc)
target_compile_options(master-g PRIVATE -g)
target_link_libraries(master-g ${GUROBI_CXX_LIBRARY} ${GUROBI_CC_LIBRARY} libpoly)

add_executable(master-oa master-oa.cc Dmax.cc)
target_compile_options(master-oa PRIVATE -g)

add_executable(robj reobjective.cc)
target_compile_options(robj PRIVATE -g)