#include <iostream>

#include <omp.h>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/imats.hpp>
#include <vcp/pdblas.hpp>
#include <vcp/pidblas.hpp>

#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

int main(void) {
	int n = 100;
	vcp::matrix< double > Ad;
	vcp::matrix< kv::dd > Add;
	vcp::matrix< kv::mpfr< 300 > > Am300;
	vcp::matrix< kv::mpfr< 500 > > Am500;

	vcp::matrix< kv::interval< double >, vcp::pidblas > Aid;
	vcp::matrix< kv::interval< kv::dd >, vcp::imats< kv::dd > > Aidd;
	vcp::matrix< kv::interval< kv::mpfr< 300 > >, vcp::imats< kv::mpfr< 300 > > > Aim300;
	vcp::matrix< kv::interval< kv::mpfr< 500 > >, vcp::imats< kv::mpfr< 500 > > > Aim500;

	Ad.rand(n);
	vcp::convert(Ad, Add);
	vcp::convert(Add, Am300);
	std::cout << Am300(0,0) << std::endl;
	vcp::convert(Am300, Am500);
	vcp::convert(Am500, Ad);
	std::cout << Ad(0, 0) << std::endl;

	vcp::convert(Am500, Aim300);
	std::cout << Aim300(0, 0) << std::endl;
	vcp::convert(Aim300, Add);
	std::cout << Ad(0, 0) << std::endl;
	vcp::convert(Aim300, Aidd);
	std::cout << Aidd(0, 0) << std::endl;
	vcp::convert(Aidd, Am500);
	std::cout << Am500(0, 0) << std::endl;
	return 0;
}