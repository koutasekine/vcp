#pragma once

#ifndef VCP_MAKE_MARRIX_LDBASE_ON_INTERVAL_HPP
#define VCP_MAKE_MARRIX_LDBASE_ON_INTERVAL_HPP

namespace vcp {
#if (defined(RDOUBLE_HPP) || (defined(DD_HPP) && defined(RDD_HPP))) && defined(RMPFR_HPP) && defined(MPFR_HPP)
	template <int _N, typename _T = kv::dd >
	class make_Matrix_ldbase_on_interval : public Legendre_Bases_Generator< kv::interval< kv::mpfr< _N > >, double, vcp::mats< double > > {
		int Div_number;
		vcp::matrix< kv::interval< _T >, vcp::imats< _T > > p_phi;

	public:	
		void make(const int mesh_order = 9) {
			if ((*this).mode <= 2) {
				std::cout << ">> ERROR : make : This function only use the Approximation mode (Mode > 2)" << std::endl;
				exit(1);
			}
			kv::mpfr< _N > mesh_size = pow(kv::mpfr< _N >(2), -mesh_order);
			(*this).Div_number = std::pow(2, mesh_order);

			int mmax = (*this).phi.size();
			vcp::matrix< kv::interval< kv::mpfr< _N >>, vcp::imats< kv::mpfr< _N > > > p_phi_T;
			
			p_phi_T.zeros(mmax, (*this).Div_number);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < mmax; i++) {
				kv::interval< kv::mpfr< _N >> interval_j;
				for (int j = 0; j < (*this).Div_number; j++) {
					interval_j.lower() = j * mesh_size;
					interval_j.upper() = (j + 1) * mesh_size;
					LegendreBaseFunctions::LegendreFunc((*this).phi[i], interval_j, p_phi_T(i, j));
				}
			}
			convert(p_phi_T, (*this).p_phi);
		}
		void save(const char* name = "") {
			std::string filename = name;
			filename += "Matrix_ldbase_on_interval_LegenderOrder";
			filename += std::to_string((*this).Order_of_Base);
			filename += "_Div";
			filename += std::to_string((*this).Div_number);
			vcp::save((*this).p_phi, filename);
		}
	};

	template < int _N >
	void make_interval_value_of_ldbase_demo(int Legendre_Order, int mesh_order,  const char* name = "") {
		vcp::make_Matrix_ldbase_on_interval< _N > A;
		A.setting(Legendre_Order, 1, 1, 1, 3);
		A.make(mesh_order);
		A.save(name);
	}
#endif
}
#endif // VCP_MAKE_MARRIX_LDBASE_ON_INTERVAL_HPP