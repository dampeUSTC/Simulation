#include <boost/python.hpp>
#include "DmpSimAlg.h"

BOOST_PYTHON_MODULE(libDmpSimAlg){
  using namespace boost::python;

  class_<DmpSimAlg,boost::noncopyable,bases<DmpVAlg> >("DmpSimAlg",init<>())
     .def("SaveTrackVertex",    &DmpSimAlg::SaveTrackVertex)
     .def("SetGPSType",    &DmpSimAlg::SetGPSType)
     .def("ResetGPS",    &DmpSimAlg::ResetGPS)
   ;
}

