#include <kv/interval.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>

int main(void) {
	vcp::matrix< kv::interval< kv::dd >, vcp::imats< kv::dd > > A, B;
	A.rand(15);
	std::cout << A << std::endl;
	B = A.submatrix({}, {});
	//std::cout << B << std::endl;
	B = A.submatrix({}, {4});
	//std::cout << B << std::endl;
	B = A.submatrix({}, {4,14});
	//std::cout << B << std::endl;
	B = A.submatrix({}, {4,4,14});
	//std::cout << B << std::endl;
	B = A.submatrix({14}, {});
	//std::cout << B << std::endl;
	B = A.submatrix({2}, {5});
	//std::cout << B << std::endl;
	B = A.submatrix({ 2 }, { 5,6 });
	//std::cout << B << std::endl;
	B = A.submatrix({ 2 }, { 5,3,11 });
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 4 }, {});
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 14 }, {14});
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 4 }, { 5,7 });
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 3, 8 }, { });
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 3, 10 }, {3});
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 3, 10 }, { 3,5 });
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 3, 10 }, { 3,5,10 });
	//std::cout << B << std::endl;

	std::cout << A.submatrix({1}, {})*A.submatrix({}, {1}) << std::endl;
	std::cout << B[1, 2] << std::endl;
}