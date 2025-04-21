#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Partition.hpp"

namespace py = pybind11;

// custom type caster for big_float
namespace pybind11 {namespace detail {
	template <> struct type_caster<big_float>{
		public:
			PYBIND11_TYPE_CASTER(big_float, _("big_float"));

			// Python -> C++ conversion
			bool load(handle src, bool) {
				PyObject* src_ptr = src.ptr();
				PyObject* float_ptr = PyNumber_Float(src_ptr);
				if (!float_ptr) {
					return false;
				}
				
				double d = PyFloat_AsDouble(float_ptr);
				Py_DECREF(float_ptr);
				value = big_float(d);
				return true;
			}

			// C++ -> Python conversion
			static handle cast(big_float src, return_value_policy, handle){
				try {
					// Get the string representation of big_float
					std::string str = src.str(50);
					if (str.empty()) {
						str = "0.0";
					}

					// Import mpmath module
					py::module mpmath = py::module::import("mpmath");
					
					// Convert string to mpmath.mpf
					py::object mpf = mpmath.attr("mpf")(str);
					
					return mpf.release();  // Return the mpmath object
				}
				catch(const std::exception& e) {
					throw std::runtime_error("Failed to convert big_float to mpmath: " + std::string(e.what()));
				}
			}
	};
}}	// namespace pybind11::detail

PYBIND11_MODULE(partition_module, m){
	// module docstring
	m.doc() = R"pbdoc(
		Partition module
		-----------------------
		A module for statistical mechanical calculations of partition functions.

		This module provides the Partition class which computes various
		thermodynamic quantities for a physical system.
		-----------------------
	)pbdoc";

	// Define the Partition class
	py::class_<Partition>(m, "Partition", R"pbdoc(
		A class to compute partition functions and thermodynamic quantities.
	)pbdoc")
	// Constructor
	.def(py::init<int, int, double, double, double, double, double, double, double, double, std::string>(),
		 py::arg("Nm"),
		 py::arg("Nmm_max"),
		 py::arg("T"),
		 py::arg("pw"),
		 py::arg("EMM_A"),
		 py::arg("EMM_N"),
		 py::arg("EMS_A"),
		 py::arg("EMS_N"),
		 py::arg("alpha"),
		 py::arg("beta"),
		 py::arg("dfile"),
		 "Initialize the Partition object with given parameters.")

	// read-write properties
	.def_readwrite("Nm",              &Partition::Nm,        "Number of monomers")
	.def_readwrite("Nmm_max",         &Partition::Nmm_max,   "Maximum number of monomer-monomer contacts")
	.def_readwrite("T",               &Partition::T,         "Temperature")
	.def_readwrite("z",               &Partition::z,         "Coordination number")
	.def_readwrite("pw",              &Partition::pw,        "Probability weight")
	.def_readwrite("EMM_A",           &Partition::EMM_A,     "Energy parameter EMM_A")
	.def_readwrite("EMM_N",           &Partition::EMM_N,     "Energy parameter EMM_N")
	.def_readwrite("EMS_A",           &Partition::EMS_A,     "Energy parameter EMS_A")
	.def_readwrite("EMS_N",           &Partition::EMS_N,     "Energy parameter EMS_N")
	.def_readwrite("alpha",           &Partition::alpha,     "Beta distribution alpha parameter")
	.def_readwrite("beta",            &Partition::beta,      "Beta distribution beta parameter")
	.def_readwrite("loc",             &Partition::loc,       "Beta distribution loc parameter")
	.def_readwrite("scale",           &Partition::scale,     "Beta distribution scale parameter")
	.def_readwrite("dumpfile",        &Partition::dumpfile,  "File to dump information to")
	.def_readwrite("Z",               &Partition::Z,         "Partition function value")
	.def_readwrite("free_energy",     &Partition::free_energy,      "Free energy of the summation")
	.def_readwrite("Cv",              &Partition::Cv,               "Heat capacity")
	.def_readwrite("ave_E",           &Partition::ave_E,            "Average energy")
	.def_readwrite("ave_Nmm",         &Partition::ave_Nmm,          "Average Nmm")
	.def_readwrite("ave_Nmm_a",       &Partition::ave_Nmm_a,        "Average aligned Nmm")
	.def_readwrite("ave_Nms_a",       &Partition::ave_Nms_a,        "Average aligned Nms")
	.def_readwrite("weights_store",   &Partition::weights_store,    "Weights from the boltzmann summation" )
	.def_readwrite("energy_store",    &Partition::energy_store,     "Energies from the boltzmann summation")
	.def_readwrite("energy2_store",   &Partition::energy2_store,    "Energy^2 from the boltzmann summation")
	.def_readwrite("Nmm_store",       &Partition::Nmm_store,        "Nmm from the boltzmann summation"     )
	.def_readwrite("Nmm_a_store",     &Partition::Nmm_a_store,      "Nmm_a from the boltzmann summation"   )
	.def_readwrite("Nms_a_store",     &Partition::Nms_a_store,      "Nms_a from the boltzmann summation"   )
	.def_readwrite("factorial_cache", &Partition::factorial_cache,  "store of pre-computed factorials")
	.def_readwrite("pw_power_cache",  &Partition::pw_power_cache,   "store of pre-computed powers of pw")
	.def_readwrite("ipw_power_cache", &Partition::ipw_power_cache,  "store of pre-computed powers of 1-pw")
	.def_readwrite("contact_constraint_constant", &Partition::contact_constraint_constant, "Discrete constraint on the number of contacts")

	// Methods
	.def("Nms_tot", &Partition::Nms_tot, py::arg("Nmm"),
			"Calculate total Nms for given Nmm")
	.def("p_chain", &Partition::p_chain, py::arg("Nmm"),
			"Calculate chain probability for given Nmm")
	.def("init_factorial_cache", &Partition::init_factorial_cache,
			"Pre-compute factorials")
	.def("init_pw_power_cache",  &Partition::init_pw_power_cache,
			"Pre-compute powers of pw")
	.def("init_ipw_power_cache", &Partition::init_ipw_power_cache,
			"Pre-compute powers of ipw")
	.def("comb_mm", &Partition::comb_mm, py::arg("Nmm"), py::arg("Nmm_a"),
			"Calculate the combinatorial factor of Nmm.")
	.def("comb_ms", &Partition::comb_ms, py::arg("Nmm"), py::arg("Nms_a"),
			"Calculate the combinatorial factor of Nms.")
	.def("Energy",  &Partition::Energy, py::arg("Nmm"), py::arg("Nmm_a"), py::arg("Nms_a"),
			"Calculate the energy of a certain configuration.")
	.def("boltzmann", &Partition::boltzmann, py::arg("Nmm"), py::arg("Nmm_a"), py::arg("Nms_a"),
			"Calculate the boltzmann weight.")
	.def("get_partition_weights", &Partition::get_partition_weights,
			"Calculate and store partition weights")
	.def("get_partition_weights_quick", &Partition::get_partition_weights_quick,
			"Calculate and store partition weights")
	.def("get_free_energy", &Partition::get_free_energy,
			"Calculate the free energy")
	.def("get_partition", &Partition::get_partition,
			"Calculate partition function")
	.def("get_average_energy", &Partition::get_average_energy,
			"Calculate average energy")
	.def("get_average_Nmm", &Partition::get_average_Nmm,
			"Calculate average Nmm")
	.def("get_average_Nmm_a", &Partition::get_average_Nmm_a,
			"Calculate average aligned Nmm")
	.def("get_average_Nms_a", &Partition::get_average_Nms_a,
			"Calculate average aligned Nms")
	.def("get_Cv", &Partition::get_Cv,
			"Calculate heat capacity")
	.def("write", &Partition::write,
			"Write results to file")

	// Add string representation
	.def("__repr__",
			[](const Partition &p) {
				return "Partition(Nm=" + std::to_string(p.Nm) + 
					", T=" + std::to_string(p.T) + 
					", EMM_A=" + std::to_string(p.EMM_A) + 
					", EMM_N=" + std::to_string(p.EMM_N) + ")";
			});

	// Add version information
	m.attr("__version__") = "1.0.0";

}
