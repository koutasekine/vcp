#include <iostream>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/pdblas.hpp>
#include <vcp/pidblas.hpp>

#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

#include <vcp/ldbase.hpp>
#include <vcp/ldbase_assist.hpp>


#include <cmath>

int main(void){
	int Legendre_Order = 50;
	for (int mesh_size = 1; mesh_size <= 5; mesh_size++){
		std::cout << ">> Making :: mesh_size = " << mesh_size << std::endl;
		vcp::make_interval_value_of_ldbase_demo< 1000 >(Legendre_Order, mesh_size, "Data_Matrix_ldbase_on_interval/");
	}
}
