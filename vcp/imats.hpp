// VCP Library
// http ://verified.computation.jp
//   
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
// Copyright(c) 2017, Kouta Sekine <k.sekine@computation.jp>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met :
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and / or other materials provided with the distribution.
// * Neither the name of the Kouta Sekine nor the names of its contributors
//   may be used to endorse or promote products derived from this software
//   without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED.IN NO EVENT SHALL KOUTA SEKINE BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

#ifndef VCP_IMATS_HPP
#define VCP_IMATS_HPP

#include <vcp/mats.hpp>

namespace vcp {
	template <typename _T, class _P = vcp::mats< _T > > class imats : public mats< kv::interval< _T > > {
	protected:
		void imid( _P& A )const{
			A.zeros((*this).row, (*this).column);
			for (int i = 0; i < A.row; i++) {
				for (int j = 0; j < A.column; j++) {
					A.v[i + (*this).row*j] = mid( (*this).v[i + (*this).row*j] );
				}
			}
		}

		//C = IA * B
		virtual void mul_im_m(const _P& B, vcp::imats< _T, _P >& c)const {
			if ((*this).type == 'S' && (B.type == 'C' || B.type == 'R' || B.type == 'M')) {
				c.row = B.row;
				c.column = B.column;
				c.n = c.row * c.column;
				c.v.resize(c.n);
				c.type = B.type;
				for (int i = 0; i < B.n; i++) {
					c.v[i] = B.v[i]*(*this).v[0];
				}
				return;
			}
			else if (((*this).type == 'C' || (*this).type == 'R' || (*this).type == 'M') && B.type == 'S') {
				c = *this;
				for (int i = 0; i < (*this).n; i++) {
					c.v[i] *= B.v[0];
				}
				return;
			}
			c.zeros((*this).row, B.column);
			if ((*this).type == 'S' && B.type == 'S') {
				c.type = 'S';
				c.v[0] = (*this).v[0] * B.v[0];
				return;
			}
			else if ((*this).type == 'R' && B.type == 'C' && (*this).column == B.row) {
				c.type = 'S';
				c.v[0] = (*this).v[0] * B.v[0];
				for (int i = 1; i < (*this).n; i++) {
					c.v[0] = (*this).v[i] * B.v[i] + c.v[0];
				}
				return;
			}
			else if ((*this).type == 'C' && B.type == 'R') {
				c.type = 'M';
				int k = 0;
				for (int i = 0; i < (*this).n; i++) {
					for (int j = 0; j < B.n; j++) {
						c.v[k] = (*this).v[i] * B.v[j];
						k++;
					}
				}
				return;
			}
			else if ((*this).type == 'M' && B.type == 'C' && (*this).column == B.row) {
				c.type = 'C';
				for (int i = 0; i < (*this).row; i++) {
					c.v[i] = (*this).v[i] * B.v[0];
					for (int j = 1; j < (*this).column; j++) {
						c.v[i] = (*this).v[i + (*this).row*j] * B.v[j] + c.v[i];
					}
				}
				return;
			}
			else if ((*this).type == 'R' && B.type == 'M' && (*this).column == B.row) {
				c.type = 'R';
				for (int i = 0; i < B.column; i++) {
					c.v[i] = (*this).v[0] * B.v[B.row*i];
					for (int j = 1; j < B.row; j++) {
						c.v[i] = (*this).v[j] * B.v[j + B.row*i] + c.v[i];
					}
				}
				return;
			}
			else if ((*this).type == 'M' && B.type == 'M' && (*this).column == B.row) {
				c.type = 'M';

				for (int i = 0; i < c.n; i++) {
					c.v[i] = kv::interval< _T >(0);
				}

				for (int k = 0; k < B.column; k++) {
					for (int j = 0; j < B.row; j++) {
						for (int i = 0; i < (*this).row; i++) {
							c.v[i + c.row*k] += (*this).v[i + (*this).row*j] * B.v[j + B.row*k];
						}
					}
				}
				return;
			}
			else {
				std::cout << "*:error " << (*this).type << " , " << B.type << std::endl;
				exit(1);
			}
		}	
		//C = A * IB
		virtual void mul_m_im(const _P& B, vcp::imats< _T, _P >& c) const {
			if (B.type == 'S' && ((*this).type == 'C' || (*this).type == 'R' || (*this).type == 'M')) {
				c.row = (*this).row;
				c.column = (*this).column;
				c.n = c.row * c.column;
				c.v.resize(c.n);
				c.type = (*this).type;
				for (int i = 0; i < (*this).n; i++) {
					c.v[i] = (*this).v[i] * B.v[0];
				}
				return;
			}
			else if ((B.type == 'C' || B.type == 'R' || B.type == 'M') && (*this).type == 'S') {
				c.row = B.row;
				c.column = B.column;
				c.n = c.row * c.column;
				c.v.resize(c.n);
				c.type = B.type;
				for (int i = 0; i < B.n; i++) {
					c.v[i] = (*this).v[0] * B.v[i];
				}
				return;
			}
			c.zeros(B.row, (*this).column);
			if (B.type == 'S' && (*this).type == 'S') {
				c.type = 'S';
				c.v[0] = B.v[0] * (*this).v[0];
				return;
			}
			else if (B.type == 'R' && (*this).type == 'C' && B.column == (*this).row) {
				c.type = 'S';
				c.v[0] = B.v[0] * (*this).v[0];
				for (int i = 1; i < B.n; i++) {
					c.v[0] = B.v[i] * (*this).v[i] + c.v[0];
				}
				return;
			}
			else if (B.type == 'C' && (*this).type == 'R') {
				c.type = 'M';
				int k = 0;
				for (int i = 0; i < B.n; i++) {
					for (int j = 0; j < (*this).n; j++) {
						c.v[k] = B.v[i] * (*this).v[j];
						k++;
					}
				}
				return;
			}
			else if (B.type == 'M' && (*this).type == 'C' && B.column == (*this).row) {
				c.type = 'C';
				for (int i = 0; i < B.row; i++) {
					c.v[i] = B.v[i] * (*this).v[0];
					for (int j = 1; j < B.column; j++) {
						c.v[i] = B.v[i + B.row*j] * (*this).v[j] + c.v[i];
					}
				}
				return;
			}
			else if (B.type == 'R' && (*this).type == 'M' && B.column == (*this).row) {
				c.type = 'R';
				for (int i = 0; i < (*this).column; i++) {
					c.v[i] = B.v[0] * (*this).v[(*this).row*i];
					for (int j = 1; j < (*this).row; j++) {
						c.v[i] = B.v[j] * (*this).v[j + (*this).row*i] + c.v[i];
					}
				}
				return;
			}
			else if (B.type == 'M' && (*this).type == 'M' && B.column == (*this).row) {
				c.type = 'M';
				for (int i = 0; i < c.n; i++) {
					c.v[i] = kv::interval< _T >(0);
				}
				for (int k = 0; k < (*this).column; k++) {
					for (int j = 0; j < (*this).row; j++) {
						for (int i = 0; i < B.row; i++) {
							c.v[i + c.row*k] += B.v[i + B.row*j] * (*this).v[j + (*this).row*k];
						}
					}
				}
				return;
			}
			else {
				std::cout << "*:error " << B.type << " , " << (*this).type << std::endl;
				exit(1);
			}
		}
		//IA = IA + B
		virtual void add_im_m(const _P& B) {
			if ((*this).row != B.row && (*this).column != B.column) {
				std::cout << "+:error " << std::endl;
				exit(1);
			}
			if ((*this).type == 'S') {
				(*this).v[0] += B.v[0];
				return;
			}
			else {
				for (int i = 0; i < (*this).n; i++) {
					(*this).v[i] += B.v[i];
				}
				return;
			}
		}

		//IC = transpose(A)*A with verification 
		virtual void vmulmm(const _P& C ){
			(*this).zeros(C.column);
			if (C.type == 'S') {
				using std::pow;
				(*this).type = 'S';
				(*this).v[0] = pow(kv::interval< _T >(C.v[0]), 2);
				return;
			}
			else if (C.type == 'C') {
				using std::pow;
				(*this).type = 'S';
				(*this).v[0] = kv::interval< _T >(0);
				for (int i = 0; i < C.n; i++) {
					(*this).v[0] += pow(kv::interval< _T >(C.v[i]), 2);
				}
			}
			else if (C.type == 'R') {
				(*this).type = 'M';
				for (int i = 0; i < C.column; i++) {
					for (int j = i; j < C.column; j++) {
						(*this).v[i + C.column*j] = kv::interval< _T >(C.v[i]) * C.v[j];
					}
				}
				for (int i = 0; i < C.column; i++) {
					for (int j = i + 1; j < C.column; j++) {
						(*this).v[j + C.column*i] = (*this).v[i + C.column*j];
					}
				}
			}
			else if (C.type == 'M') {
				(*this).type = 'M';
				for (int i = 0; i < C.column; i++) {
					for (int j = i; j < C.column; j++) {
						for (int k = 0; k < C.row; k++) {
							(*this).v[i + C.column*j] += kv::interval< _T >(C.v[k + C.row*i]) * C.v[k + C.row*j];
						}
					}
				}
				for (int i = 0; i < C.column; i++) {
					for (int j = i + 1; j < C.column; j++) {
						(*this).v[j + C.column*i] = (*this).v[i + C.column*j];
					}
				}
			}
			else {
				std::cout << "*:error " << C.type << std::endl;
				exit(1);
			}
		}

	public:
		void linearsolve(imats< _T, _P >& b, imats< _T, _P >& x) {
			if ((*this).row != (*this).column || b.row != (*this).row) {
				std::cout << "error : linearsolve row != column || b.row != row" << std::endl;
				exit(1);
			}
			_P mx, mb, R;
			(*this).imid(R);
			b.imid(mb);
			R.inv();				
			R.mulmm(mb, mx);
			// x = A*mx
			(*this).mul_im_m(mx, x);
			// x = x - b, (x = A*mx - b)
			x.subsmmA(b);
			// x = |A*mx - b| 
			x.abs();
			// b = R*A
			(*this).mul_m_im( R, b);
			_T one = _T(1);
			// b = RA - I
			for (int i = 0; i < (*this).row; i++) {
				b.v[i + (*this).row*i] -= one;
			}
			imats< _T, _P > G, T;
			// G = || RA - I ||_inf
			b.norminf(G);
			if (G.v[0].upper() >= one) {
				std::cout << "error : linearsolve verification is failed ||RA - I|| <= " << G.v[0].upper() << std::endl;
				exit(1);
			}
			// T = || A*mx - b ||_inf
			x.norminf(T);
			// b = | A*mx - b |
			b.abs();
			// G = || A*mx - b ||_inf/(1 - || RA - I ||_inf)
			G.v[0] = T.v[0] / (one - G.v[0]);
			// T = || A*mx - b ||_inf/(1 - || RA - I ||_inf)*e
			T.zeros(b.row, b.column);
			for (int i = 0; i < b.row; i++) {
				for (int j = 0; j < b.column; j++){
					T.v[i + b.row*j] = G.v[0];
				}
			}
			// G = |RA-I|*||Ax - b||/(1 - ||RA-I||)*e
			b.mulmm(T, G);
			// x = |A*mx - b| + |RA-I|*||Ax - b||/(1 - ||RA-I||)*e
			x.addmm(G);
			// x = infsup(-x.upper, x.upper)
			for (int i = 0; i < x.row; i++) {
				for (int j = 0; j < x.column; j++) {
					x.v[i+x.row*j].lower() = -x.v[i + x.row*j].upper();
				}
			}
			// x = x + mx;
			x.add_im_m(mx);	
		}

		void eigsym(int itep = 1) {
			_P mA, V;
			imats< _T, _P > C;
			// mA = mid(A);
			(*this).imid(mA);
			// mA <- eigehvalue, V <- eigenvector
			mA.eigsym(V);
			// C = A*V;
			(*this).mul_im_m(V, C);
			// C = A*V - V*Lambda;
			for (int i = 0; i < (*this).row; i++) {
				for (int j = 0; j < (*this).row; j++) {
					C.v[i + (*this).row*j] -= V.v[i + (*this).row*j]*kv::interval< _T >(mA.v[j + (*this).row*j]);
				}
			}
			imats< _T, _P > Ginf, G, Tinf, T;
			// ||A*V - V*Lambda||_inf <= Tinf
			C.norminf(Tinf);
			// ||A*V - V*Lambda||_1 <= T
			C.normone(T);
			// T = sqrt(Tone*Tinf); ||A*V - V*Lambda||_2 <= sqrt(||A*V - V*Lambda||_1 * ||A*V - V*Lambda||_inf)
			T.v[0] = sqrt(T.v[0] *Tinf.v[0]);
			int thisrow = (*this).row;
			// (*this) = transpose(V)*V, C is Transpose matrix.
			(*this).vmulmm(V);
			_T one = _T(1);
			// C = transpose(V)*V - I, C is Transpose matrix.
			for (int i = 0; i < (*this).row; i++) {
				(*this).v[i + (*this).row*i] -= one;
			}
			// ||transpose(V)*V - I||_inf <= Ginf
			(*this).norminf(Ginf);
			// ||transpose(V)*V - I||_1 <= G
			(*this).normone(G);
			// G = sqrt(Gone*Ginf); ||transpose(V)*V - I||_2 <= sqrt(||transpose(V)*V - I||_1 * ||transpose(V)*V - I||_inf)
			G.v[0] = sqrt(G.v[0] *Ginf.v[0]);
			if (G.v[0].upper() >= one) {
				std::cout << "error : eigsym verification is failed || I - X^T X|| <= " << G.v[0].upper() << std::endl;
				exit(1);
			}
			G.v[0] = T.v[0] / (one - G.v[0]);
			T.zeros(thisrow, 1);
			for (int i = 0; i < thisrow; i++) {
				T.v[i].upper() = G.v[0].upper();
				T.v[i].lower() = -T.v[i].upper();
			}
			(*this).zeros(mA.row, mA.column);
			for (int i = 0; i < thisrow; i++) {
				(*this).v[i + thisrow*i] = mA.v[i + thisrow*i] + T.v[i];
			}

		}

		void eigsymge(imats< _T, _P >& B, int itep = 1) {
			int thisrow = (*this).row;
			imats< _T, _P > Binv_norm;
			{
				imats< _T, _P > BB, DD;
				BB = B;
				B.eigsym();
				B.diag(DD);
				DD.min(Binv_norm);
				if (Binv_norm.v[0].lower() <= 0) {
					std::cout << "error : eigsymge verification is failed. B is not positive " << Binv_norm.v[0].lower() << std::endl;
					exit(1);
				}
				B = BB;
			}
			
			_P mA, V;
			imats< _T, _P > C;
			{
				_P mB;
				// mA = mid(A);
				(*this).imid(mA);
				B.imid(mB);
				// mA <- eigehvalue, V <- eigenvector
				mA.eigsymge(mB, V);
			}
			// C = A*V;
			(*this).mul_im_m(V, C);
			// (*this) = B*V;
			B.mul_im_m(V, (*this));

			// C = A*V - B*V*Lambda;
			for (int i = 0; i < thisrow; i++) {
				for (int j = 0; j < thisrow; j++) {
					C.v[i + thisrow*j] -= (*this).v[i + thisrow*j] * kv::interval< _T >(mA.v[j + thisrow*j]);
				}
			}

			imats< _T, _P > Ginf, G, Tinf, T;
			// ||A*V - V*Lambda||_inf <= Tinf
			C.norminf(Tinf);
			// ||A*V - V*Lambda||_1 <= T
			C.normone(T);
			// T = sqrt(Tone*Tinf); ||A*V - B*V*Lambda||_2 <= sqrt(||A*V - B*V*Lambda||_1 * ||A*V - B*V*Lambda||_inf)
			T.v[0] = sqrt(T.v[0] * Tinf.v[0]);
			
			// C = B*V
			B.mul_im_m(V, (*this));
			{
				_P Vt;
				V.transpose(Vt);
				// (*this) = transpose(V)*B*V
				(*this).mul_m_im(Vt, C);
				(*this) = C;
			}
			_T one = _T(1);
			// C = transpose(V)*B*V - I, C is Transpose matrix.
			for (int i = 0; i < (*this).row; i++) {
				(*this).v[i + (*this).row*i] -= one;
			}
			// ||transpose(V)*V - I||_inf <= Ginf
			(*this).norminf(Ginf);
			// ||transpose(V)*V - I||_1 <= G
			(*this).normone(G);
			// G = sqrt(Gone*Ginf); ||transpose(V)*V - I||_2 <= sqrt(||transpose(V)*V - I||_1 * ||transpose(V)*V - I||_inf)
			G.v[0] = sqrt(G.v[0] * Ginf.v[0]);
			if (G.v[0].upper() >= one) {
				std::cout << "error : eigsym verification is failed || I - X^T X|| <= " << G.v[0].upper() << std::endl;
				exit(1);
			}

			G.v[0] = T.v[0] / ((one - G.v[0])*Binv_norm.v[0]);
			T.zeros(thisrow, 1);
			for (int i = 0; i < thisrow; i++) {
				T.v[i].upper() = G.v[0].upper();
				T.v[i].lower() = -T.v[i].upper();
			}
			(*this).zeros(mA.row, mA.column);
			for (int i = 0; i < thisrow; i++) {
				(*this).v[i + thisrow*i] = mA.v[i + thisrow*i] + T.v[i];
			}
		}

	};
}
#endif // VCP_IMATS_HPP