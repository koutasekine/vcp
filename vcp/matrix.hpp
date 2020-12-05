#pragma once

#ifndef VCP_MATRIX_HPP
#define VCP_MATRIX_HPP

#include <vcp/mats.hpp>

namespace vcp {
	template <typename _T, class _P = mats< _T >> class matrix : protected _P {
	public:
		matrix< _T, _P >() {
			(*this).row = 0;
			(*this).column = 0;
			(*this).n = 0;
			(*this).type = 'N';
		}
		~matrix< _T, _P >() = default;
		matrix< _T, _P >(const matrix< _T, _P >&) = default;
		matrix< _T, _P >(matrix< _T, _P >&&) = default;
		matrix< _T, _P >& operator=(const matrix< _T, _P >& A) = default;
		matrix< _T, _P >& operator=(matrix< _T, _P >&& A) = default;

		_T& operator () (const int i) {
			return (*this).v[i];
		}
		_T operator () (const int i) const {
			return (*this).v[i];
		}
		_T& operator () (const int i, const int j) {
			if ((*this).type == 'R') {
				return (*this).v[j];
			}
			else {
				return (*this).v[i + (*this).row*j];
			}
		}
		_T operator () (const int i, const int j)const {
			if ((*this).type == 'R') {
				return (*this).v[j];
			}
			else {
				return (*this).v[i + (*this).row*j];
			}
		}
	/*
		_T& operator [] (const int i) {
			return (*this).v[i];
		}
		_T operator [] (const int i) const {
			return (*this).v[i];
		}
		*/
		
		matrix< _T, _P > submatrix(const std::initializer_list<int>& list1, const std::initializer_list<int>& list2) const {
			matrix< _T, _P > A;
			(*this).submat(A, list1, list2);
			return std::move(A);
		}

		int elementsize()const { return (*this).n; }
		int columnsize()const { return (*this).column; }
		int rowsize()const { return (*this).row; }
		char matstype()const {
			return (*this).type;
		}
		std::vector< _T > vecpointer()const {
			return (*this).v;
		}

		void eye(const int r) { _P::eye(r); }
		void ones(const int i) { _P::ones(i); }
		void ones(const int r, const int c) { _P::ones(r, c); }
		void zeros(const int i) { _P::zeros(i); }
		void zeros(const int r, const int c) { _P::zeros(r, c); }
		void rand(const int i) { _P::rand(i); }
		void rand(const int r, const int c) { _P::rand(r, c); }
		void resize(const int i, const int j) { _P::resize(i, j); }
		void clear() { _P::clear(); }

		//***************** Operator Overload *****************//
		friend matrix< _T, _P > operator+(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			C = A;
			C.addmm(B);
			return std::move(C);
		}
		friend matrix< _T, _P > operator+(matrix< _T, _P >&& A, const matrix< _T, _P >& B) {
			A.addmm(B);
			return std::move(A);
		}
		friend matrix< _T, _P > operator+(const matrix< _T, _P >& A, matrix< _T, _P >&& B) {
			B.addmm(A);
			return std::move(B);
		}
		friend matrix< _T, _P > operator+(matrix< _T, _P >&& A, matrix< _T, _P >&& B) {
			A.addmm(B);
			return std::move(A);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator+(const _Tm a, const matrix< _T, _P >& B) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.addsm(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator+(const _Tm a, matrix< _T, _P >&& B) {
			_T Ta = _T(a);
			B.addsm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator+(const matrix< _T, _P >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.addms(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator+(matrix< _T, _P >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.addms(Ta);
			return std::move(B);
		}
		friend matrix< _T, _P > operator+(const matrix< _T, _P >& A) {
			//		A.plusm();
			return A;
		}
		friend matrix< _T, _P > operator+(matrix< _T, _P >&& A) {
			//		A.plusm();
			return std::move(A);
		}

		friend matrix< _T, _P >& operator+=(matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			A.addmm(B);
			return A;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P >& >::type operator+=(matrix< _T, _P >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.addms(Ta);
			return A;
		}

		friend matrix< _T, _P > operator-(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			C = A;
			C.subsmmA(B);
			return std::move(C);
		}
		friend matrix< _T, _P > operator-(matrix< _T, _P >&& A, const matrix< _T, _P >& B) {
			A.subsmmA(B);
			return std::move(A);
		}
		friend matrix< _T, _P > operator-(const matrix< _T, _P >& A, matrix< _T, _P >&& B) {
			B.subsmmB(A);
			return std::move(B);
		}
		friend matrix< _T, _P > operator-(matrix< _T, _P >&& A, matrix< _T, _P >&& B) {
			A.subsmmA(B);
			return std::move(A);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator-(const _Tm a, const matrix< _T, _P >& B) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.subssm(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator-(const _Tm a, matrix< _T, _P >&& B) {
			_T Ta = _T(a);
			B.subssm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator-(const matrix< _T, _P >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.subsms(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator-(matrix< _T, _P >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.subsms(Ta);
			return std::move(B);
		}
		friend matrix< _T, _P > operator-(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.minusm();
			return std::move(C);
		}
		friend matrix< _T, _P > operator-(matrix< _T, _P >&& A) {
			A.minusm();
			return std::move(A);
		}

		friend matrix< _T, _P >& operator-=(matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			A.subsmmA(B);
			return A;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P >& >::type operator-=(matrix< _T, _P >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.subsms(Ta);
			return A;
		}

		friend matrix< _T, _P > operator*(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			A.mulmm(B, C);
			return std::move(C);
		}
		friend matrix< _T, _P > operator*(matrix< _T, _P >&& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			A.mulmm(B, C);
			A = std::move(C);
			return std::move(A);
		}
		friend matrix< _T, _P > operator*(const matrix< _T, _P >& A, matrix< _T, _P >&& B) {
			matrix< _T, _P > C;
			A.mulmm(B, C);
			B = std::move(C);
			return std::move(B);
		}
		friend matrix< _T, _P > operator*(matrix< _T, _P >&& A, matrix< _T, _P >&& B) {
			matrix< _T, _P > C;
			A.mulmm(B, C);
			A = std::move(C);
			return std::move(A);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator*(const _Tm a, const matrix< _T, _P >& B) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.mulsm(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator*(const _Tm a, matrix< _T, _P >&& B) {
			_T Ta = _T(a);
			B.mulsm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator*(const matrix< _T, _P >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.mulms(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator*(matrix< _T, _P >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.mulms(Ta);
			return std::move(B);
		}

		friend matrix< _T, _P >& operator*=(matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			A = A * B;
			return A;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P >& >::type operator*=(matrix< _T, _P >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.mulms(Ta);
			return A;
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator/(const _Tm a, const matrix< _T, _P >& B) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.divsm(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator/(const _Tm a, matrix< _T, _P >&& B) {
			_T Ta = _T(a);
			B.divsm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator/(const matrix< _T, _P >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.divms(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator/(matrix< _T, _P >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.divms(Ta);
			return std::move(B);
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P >& >::type operator/=(matrix< _T, _P >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.divms(Ta);
			return A;
		}

		friend matrix< _T, _P > pow(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			C = A;
			C.powmmA(B);
			return std::move(C);
		}
		friend matrix< _T, _P > pow(matrix< _T, _P >&& A, const matrix< _T, _P >& B) {
			A.powmmA(B);
			return std::move(A);
		}
		friend matrix< _T, _P > pow(const matrix< _T, _P > A, matrix< _T, _P >&& B) {
			B.powmmB(A);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type pow(const matrix< _T, _P >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.powms(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type pow(matrix< _T, _P >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.powms(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type pow(const _Tm a, const matrix< _T, _P >& B) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.powsm(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type pow(const _Tm a, matrix< _T, _P >&& B) {
			_T Ta = _T(a);
			B.powsm(Ta);
			return std::move(B);
		}


		//***************** Special Matrix Arithemetic *****************//
		friend matrix< _T, _P > ltransmul(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			A.mulltmm(C);
			return std::move(C);
		}
		friend matrix< _T, _P > ltransmul(matrix< _T, _P >&& A) {
			matrix< _T, _P > C;
			A.mulltmm(C);
			A = std::move(C);
			return std::move(A);
		}

		//***************** Math functions *****************//
		friend matrix< _T, _P > abs(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.abs();
			return std::move(C);
		}
		friend matrix< _T, _P > abs(matrix< _T, _P >&& A) {
			A.abs();
			return std::move(A);
		}
		friend matrix< _T, _P > sqrt(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.sqrt();
			return std::move(C);
		}
		friend matrix< _T, _P > sqrt(matrix< _T, _P >&& A) {
			A.sqrt();
			return std::move(A);
		}
		friend matrix< _T, _P > sin(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.sin();
			return std::move(C);
		}
		friend matrix< _T, _P > sin(matrix< _T, _P >&& A) {
			A.sin();
			return std::move(A);
		}
		friend matrix< _T, _P > cos(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.cos();
			return std::move(C);
		}
		friend matrix< _T, _P > cos(matrix< _T, _P >&& A) {
			A.cos();
			return std::move(A);
		}
		friend matrix< _T, _P > exp(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.exp();
			return std::move(C);
		}
		friend matrix< _T, _P > exp(matrix< _T, _P >&& A) {
			A.exp();
			return std::move(A);
		}
		friend matrix< _T, _P > log(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.log();
			return std::move(C);
		}
		friend matrix< _T, _P > log(matrix< _T, _P >&& A) {
			A.log();
			return std::move(A);
		}

		//************* matlab like functions *************//	
		friend matrix< _T, _P > sum(const matrix< _T, _P >& A) {
			matrix< _T, _P > c;
			A.sum(c);
			return std::move(c);
		}
		friend matrix< _T, _P > diag(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			A.diag(C);
			return std::move(C);
		}
		friend matrix< _T, _P > transpose(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			A.transpose(C);
			return std::move(C);
		}
		friend matrix< _T, _P > max(const matrix< _T, _P >& A) {
			matrix< _T, _P > c;
			A.max(c);
			return std::move(c);
		}
		friend matrix< _T, _P > min(const matrix< _T, _P >& A) {
			matrix< _T, _P > c;
			A.min(c);
			return std::move(c);
		}
		friend matrix< _T, _P > normone(const matrix< _T, _P >& A) {
			matrix< _T, _P > c;
			A.normone(c);
			return std::move(c);
		}
		friend matrix< _T, _P > norminf(const matrix< _T, _P >& A) {
			matrix< _T, _P > c;
			A.norminf(c);
			return std::move(c);
		}

		friend matrix< _T, _P > tril(const matrix< _T, _P >& A) {
			matrix< _T, _P > c = A;
			c.tril();
			return std::move(c);
		}
		friend matrix< _T, _P > tril(matrix< _T, _P >&& A) {
			A.tril();
			return std::move(A);
		}
		friend matrix< _T, _P > triu(const matrix< _T, _P >& A) {
			matrix< _T, _P > c = A;
			c.triu();
			return std::move(c);
		}
		friend matrix< _T, _P > triu(matrix< _T, _P >&& A) {
			A.triu();
			return std::move(A);
		}

		friend matrix< _T, _P > SymTridiagonalization(const matrix< _T, _P >& A) {
			matrix< _T, _P > AA;
			AA = A;
			AA.SymTridiagonalization();
			return AA;
		}
		friend matrix< _T, _P > SymTridiagonalization(matrix< _T, _P >&& A) {
			A.SymTridiagonalization();
			return A;
		}

		friend matrix< _T, _P > HessenbergTrans(const matrix< _T, _P >& A) {
			matrix< _T, _P > AA;
			AA = A;
			AA.HessenbergTrans();
			return AA;
		}
		friend matrix< _T, _P > HessenbergTrans(matrix< _T, _P >&& A) {
			A.HessenbergTrans();
			return A;
		}

		//  L = tril(ipiv*A) + eye
		//  U = triu(ipiv*A) + eye
		friend void lu(matrix< _T, _P >& A, matrix< int >& ipiv) {
			A.ludecomposition(ipiv);
		}
		// [A, ipiv] = lu(A)
		//  L = tril(ipiv*LU) + eye
		//  U = triu(ipiv*LU) + eye
		friend void lu(const matrix< _T, _P >& A, matrix< _T, _P >& LU, matrix< int >& ipiv) {
			LU = A;
			LU.ludecomposition(ipiv);
		}
		// ipiv * A = L * U
		friend void lu(const matrix< _T, _P >& A, matrix< _T, _P >& L, matrix< _T, _P >& U, matrix< int >& ipiv) {
			U = A;
			U.ludecomposition(ipiv);
			U.LUtoLandU(L, ipiv);
		}

		friend void Trilu(const matrix< _T, _P >& A, matrix< _T, _P >& LU, matrix< int >& ipiv) {
			LU = A;
			LU.TriLudecomposition(ipiv);
		}

		// [Q, R] = qr(A), A=Q*R
		friend void qr(const matrix< _T, _P >& A, matrix< _T, _P >& Q, matrix< _T, _P >& R) {
			R = A;
			R.Householder_qrdecomposition(Q);
		}
		// [Q, R] = qr(A), A=Q*R Hessenberg
		friend void qr_Hessenberg(const matrix< _T, _P >& A, matrix< _T, _P >& Q, matrix< _T, _P >& R) {
			R = A;
			R.Householder_qrdecomposition_Hessenberg(Q);
		}
		// [Q, R] = qr(A), A=Q*R Tridiagonal
		friend void qr_Tridiag(const matrix< _T, _P >& A, matrix< _T, _P >& Q, matrix< _T, _P >& R) {
			R = A;
			R.Householder_qrdecomposition_Tridiag(Q);
		}

		// x = A\b
		friend void lss(const matrix< _T, _P >& A, const matrix< _T, _P >& b, matrix< _T, _P >& x) {
			matrix< _T, _P > AA = A;
			matrix< _T, _P > bb = b;
			AA.linearsolve(bb, x);
		}
		friend matrix< _T, _P > lss(const matrix< _T, _P >& A, const matrix< _T, _P >& b) {
			matrix< _T, _P > AA = A;
			matrix< _T, _P > bb = b;
			matrix< _T, _P > x;
			AA.linearsolve(bb, x);
			return std::move(x);
		}
		friend matrix< _T, _P > lss(matrix< _T, _P >&& A, const matrix< _T, _P >& b) {
			matrix< _T, _P > bb = b;
			matrix< _T, _P > x;
			A.linearsolve(bb, x);
			return std::move(x);
		}
		friend matrix< _T, _P > lss(const matrix< _T, _P >& A, matrix< _T, _P >&& b) {
			matrix< _T, _P > AA = A;
			matrix< _T, _P > x;
			AA.linearsolve(b, x);
			return std::move(x);
		}
		friend matrix< _T, _P > lss(matrix< _T, _P >&& A, matrix< _T, _P >&& b) {
			matrix< _T, _P > x;
			A.linearsolve(b, x);
			return std::move(x);
		}
		// R = inv(A)
		friend matrix< _T, _P > inv(const matrix< _T, _P >& A) {
			matrix< _T, _P > R = A;
			R.inv();
			return std::move(R);
		}
		friend matrix< _T, _P > inv(matrix< _T, _P >&& A) {
			A.inv();
			return std::move(A);
		}
		friend matrix< _T, _P > Cholesky(const matrix< _T, _P >& A) {
			matrix< _T, _P > C = A;
			C.Cholesky();
			return std::move(C);
		}
		friend matrix< _T, _P > Cholesky(matrix< _T, _P >&& A) {
			A.Cholesky();
			return std::move(A);
		}
		friend void eigsym(const matrix< _T, _P >& A, matrix< _T, _P >& E, int itep = 1) {
			E = A;
			E.eigsym(itep);
		}
		friend void eigsym(const matrix< _T, _P >& A, matrix< _T, _P >& E, matrix< _T, _P >& V, int itep = 1) {
			E = A;
			E.eigsym(V, itep);
		}
		friend void eigsymge(const matrix< _T, _P >& A, const matrix< _T, _P >& B, matrix< _T, _P >& E, int itep = 1) {
			E = A;
			matrix< _T, _P > BD = B;
			E.eigsymge(BD, itep);
		}
		friend void eigsymge(const matrix< _T, _P >& A, const matrix< _T, _P >& B, matrix< _T, _P >& E, matrix< _T, _P >& V, int itep = 1) {
			E = A;
			matrix< _T, _P > BD = B;
			E.eigsymge(BD, V, itep);
		}

		//Matlab C = [A,B];
		friend matrix< _T, _P > horzcat(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			A.horzcat(B, C);
			return std::move(C);
		}
		template<typename... Args> friend matrix< _T, _P > horzcat(const matrix< _T, _P >& A, const matrix< _T, _P >& B, const Args&... args) {
			matrix< _T, _P > C;
			A.horzcat(B, C);
			return std::move(horzcat(C, args...));
		}
		
		//Matlab C = [A;B];
		friend matrix< _T, _P > vercat(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			A.vercat(B, C);
			return std::move(C);
		}
		template<typename... Args> friend matrix< _T, _P > vercat(const matrix< _T, _P >& A, const matrix< _T, _P >& B, const Args&... args) {
			matrix< _T, _P > C;
			A.vercat(B, C);
			return std::move(vercat(C, args...));
		}

		friend int length(matrix< _T, _P >& A) {
			return A.length();
		}

		//**************** display function ***************//
		friend std::ostream& operator<<(std::ostream& os, matrix< _T, _P >& A) {
			return A.display(os);
		}
		friend std::ostream& operator<<(std::ostream& os, matrix< _T, _P >&& A) {
			return A.display(os);
		}
	};
}
#endif // VCP_MATRIX_HPP
