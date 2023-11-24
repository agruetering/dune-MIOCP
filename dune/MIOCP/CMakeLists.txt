#install headers
install(FILES probleminterface.hh 
	      Parameter.hh
	      heatdriver.hh 
	      adjointdriver.hh
	      linearheatfem.hh
	      adjointheatfem.hh
	      hfunctions.hh
	      activeset.hh
	      outerapprox.hh
	      linsolver.hh
	      gurobi.hh
	      callback.hh
	      switchconstr.hh
	      Dmax.hh
 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/MIOCP)

