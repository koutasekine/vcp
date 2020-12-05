#pragma once

#ifndef VCP_PIDBLAS_HPP
#define VCP_PIDBLAS_HPP

#include <vcp/pdblas.hpp>
#include <vcp/imats.hpp>

namespace vcp {
	namespace pidblas_assist{
		void midrad(const vcp::mats< kv::interval< double > >& A, vcp::pdblas& B, vcp::pdblas& C) {
			B.zeros(A.row, A.column);
			C.zeros(A.row, A.column);
			for (int i = 0; i < A.row; i++) {
				for (int j = 0; j < A.column; j++) {
					midrad(A.v[i + A.row*j], B.v[i + A.row*j], C.v[i + A.row*j]);
				}
			}
		}
		void midrad(const vcp::imats< double, vcp::pdblas >& A, vcp::pdblas& B, vcp::pdblas& C) {
			B.zeros(A.row, A.column);
			C.zeros(A.row, A.column);
			for (int i = 0; i < A.row; i++) {
				for (int j = 0; j < A.column; j++) {
					midrad(A.v[i + A.row*j], B.v[i + A.row*j], C.v[i + A.row*j]);
				}
			}
		}
	}

	class pidblas : public vcp::imats< double, vcp::pdblas > { 
	protected:	
		//C = IA * B
		void mul_im_m(const vcp::pdblas& B, vcp::imats< double, vcp::pdblas >& c) const override {
			if ((*this).type == 'S' && (B.type == 'C' || B.type == 'R' || B.type == 'M')) {
				c.row = B.row;
				c.column = B.column;
				c.n = c.row * c.column;
				c.v.resize(c.n);
				c.type = B.type;
				for (int i = 0; i < B.n; i++) {
					c.v[i] = B.v[i] * (*this).v[0];
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
			vcp::pdblas lC, uC;
			{
				vcp::pdblas mA, rA, temp;
				pidblas_assist::midrad((*this), mA, rA);

				//DOWN
				kv::hwround::rounddown();
				// lC =  mA*B
				mA.mulmm(B, lC);

				//UP
				kv::hwround::roundup();
				// uC =  mA*B
				mA.mulmm(B, uC);
				//mA = B
				mA = B;
				//mA = |B|
				mA.abs();
				// temp = |B|*rA
				rA.mulmm(mA, temp);
				uC.addmm(temp);
				kv::hwround::rounddown();
				uC.subsmmA(temp);
			}
			kv::hwround::roundnear();
			c.zeros(uC.row, uC.column);
			for (int i = 0; i < c.n; i++) {
				c.v[i].lower() = lC.v[i];
				c.v[i].upper() = uC.v[i];
			}
		}
		//C = B * IA
		void mul_m_im(const vcp::pdblas& B, vcp::imats< double, vcp::pdblas >& c) const override {
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
			vcp::pdblas lC, uC;
			{
				vcp::pdblas mA, rA, temp;
				pidblas_assist::midrad((*this), mA, rA);

				//DOWN
				kv::hwround::rounddown();
				// lC =  B*mA
				B.mulmm(mA, lC);

				//UP
				kv::hwround::roundup();
				// uC =  B*mA
				B.mulmm(mA, uC);
				//mA = B
				mA = B;
				//mA = |B|
				mA.abs();
				// temp = |B|*rA
				mA.mulmm(rA, temp);
				uC.addmm(temp);
				kv::hwround::rounddown();
				uC.subsmmA(temp);
			}
			kv::hwround::roundnear();
			c.zeros(uC.row, uC.column);
			for (int i = 0; i < c.n; i++) {
				c.v[i].lower() = lC.v[i];
				c.v[i].upper() = uC.v[i];
			}
		}
		//IC = transpose(A)*A with verification 
		void vmulmm(const vcp::pdblas& C) override {
			(*this).zeros(C.column);
			vcp::pdblas tmp;
			kv::hwround::roundup();
			C.mulltmm(tmp);
			if (tmp.type == 'S') {
				(*this).v[0].upper() = tmp.v[0];
			}
			else if(tmp.type == 'M'){
				for (int i = 0; i < tmp.row; i++) {
					for (int j = 0; j < tmp.column; j++) {
						(*this).v[i + (*this).row*j].upper() = tmp.v[i + tmp.row*j];
					}
				}
			}
			else {
				kv::hwround::roundnear();
				std::cout << "Error : type miss" << std::endl;
				exit(1);
			}

			kv::hwround::rounddown();
			C.mulltmm(tmp);
			if (tmp.type == 'S') {
				(*this).v[0].lower() = tmp.v[0];
			}
			else if (tmp.type == 'M') {
				for (int i = 0; i < tmp.row; i++) {
					for (int j = 0; j < tmp.column; j++) {
						(*this).v[i + (*this).row*j].lower() = tmp.v[i + tmp.row*j];
					}
				}
			}
			else {
				kv::hwround::roundnear();
				std::cout << "Error : type miss" << std::endl;
				exit(1);
			}
			kv::hwround::roundnear();
		}
	public:
		//C = IA * IB
		void mulmm(const vcp::mats< kv::interval< double > >& B, vcp::mats< kv::interval< double > >& c)const override {
			kv::hwround::roundnear();
			if (type == 'S' && (B.type == 'C' || B.type == 'R' || B.type == 'M')) {
				c = B;
				for (int i = 0; i < B.n; i++) {
					c.v[i] *= v[0];
				}
				return;
			}
			else if ((type == 'C' || type == 'R' || type == 'M') && B.type == 'S') {
				c = *this;
				for (int i = 0; i < n; i++) {
					c.v[i] *= B.v[0];
				}
				return;
			}
			
			vcp::pdblas lC, uC;
			{
				vcp::pdblas mA, rA, mB, rB, temp;
				pidblas_assist::midrad((*this), mA, rA);
				pidblas_assist::midrad(B, mB, rB);
			
				//DOWN
				kv::hwround::rounddown();
				// lC =  mA*mB
				mA.mulmm(mB, lC);
				
				//UP
				kv::hwround::roundup();
				// uC =  mA*mB
				mA.mulmm(mB, uC);
				//mA = |mA|
				mA.abs();
				//mB = |mB|
				mB.abs();
				// mA = |mA| + rA
				mA.addmm(rA);
				// temp = (|mA|+rA)*rB
				mA.mulmm(rB, temp);
				// mA = rA*|mB|
				rA.mulmm(mB, mA);
				// temp = (|mA|+rA)*rB + rA*|mB|
				temp.addmm(mA);
				uC.addmm(temp);
				kv::hwround::rounddown();
				lC.subsmmA(temp);
			}
			kv::hwround::roundnear();
			c.zeros(uC.row, uC.column);
			for (int i = 0; i < c.n; i++) {
				c.v[i].lower() = lC.v[i];
				c.v[i].upper() = uC.v[i];
			}
		}
		// C = transpose(A)*A : multiplication left side transpose
		void mulltmm(vcp::mats< kv::interval< double > >& c)const override{
			c.zeros((*this).column, (*this).column);
			if ((*this).type == 'S') {
				using std::pow;
				c.v[0] = (*this).v[0] * (*this).v[0];
				return;
			}
			else if ((*this).type == 'M' || (*this).type == 'R' || (*this).type == 'C') {
				vcp::pdblas lC, uC;
				{
					vcp::pdblas mA, mAt, rA, temp;
					pidblas_assist::midrad((*this), mA, rA);
					//DOWN
					kv::hwround::rounddown();
					// lC =  transpose(mA)*mA
					mA.mulltmm(lC);

					//UP
					kv::hwround::roundup();
					// uC =  transpose(mA)*mA
					mA.mulltmm(uC);
					//mA = |mA|
					mA.abs();
					//mAt = transpose(|mA|)
					mA.transpose(mAt);
					// temp = transpose(|mA|)*rA
					mAt.mulmm(rA, temp);
					// mAt = transpose(rA)*rA
					rA.mulltmm(mAt);
					for (int i = 0; i < temp.row; i++) {
						for (int j = i; j < temp.column; j++) {
							temp.v[i + temp.row * j] += temp.v[j + temp.row * i] + mAt.v[i + temp.row * j];
						}
					}
					for (int i = 0; i < temp.row; i++) {
						for (int j = i+1; j < temp.column; j++) {
							temp.v[j + temp.row * i] = temp.v[i + temp.row * j];
						}
					}

					// temp = (|mA|+rA)*rB + rA*|mB|
					uC.addmm(temp);
					kv::hwround::rounddown();
					lC.subsmmA(temp);
				}
				kv::hwround::roundnear();
				c.zeros(uC.row, uC.column);
				for (int i = 0; i < c.n; i++) {
					c.v[i].lower() = lC.v[i];
					c.v[i].upper() = uC.v[i];
				}

				for (int i = 0; i < (*this).column; i++) {
					for (int j = i + 1; j < (*this).column; j++) {
						c.v[j + (*this).column*i] = c.v[i + (*this).column*j];
					}
				}
			}
			else {
				std::cout << "*:error " << (*this).type << std::endl;
				exit(1);
			}
		}
	};
}
#endif // VCP_PIDBLAS_HPP