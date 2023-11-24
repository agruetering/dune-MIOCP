set(modules DuneMIOCPMacros.cmake FINDGUROBI.cmake)

install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
