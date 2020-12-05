#pragma once

#ifndef VCP_CONVERTER_HPP
#define VCP_CONVERTER_HPP

namespace vcp {
#if defined(RDOUBLE_HPP) && defined(RMPFR_HPP) && defined(MPFR_HPP)
	template <int N> void convert(const kv::mpfr<N>& x, double& y, int rnd = 0)
	{
		mp_rnd_t mode;

		if (rnd == 1) {
			mode = MPFR_RNDU;
		}
		else if (rnd == -1) {
			mode = MPFR_RNDD;
		}
		else {
			mode = MPFR_RNDN;
		}

		y = mpfr_get_d(x.a, mode);
	}
	template <int N> void convert(const double& x, kv::mpfr<N>& y) {
		y = x;
	}
#endif

#if defined(DD_HPP) && defined(RDD_HPP) && defined(RMPFR_HPP) && defined(MPFR_HPP)
	template <int N> void convert(const kv::dd& x, kv::mpfr<N>& y, int rnd = 0)
	{
		mp_rnd_t mode;
		if (rnd == 1) {
			mode = MPFR_RNDU;
		}
		else if (rnd == -1) {
			mode = MPFR_RNDD;
		}
		else {
			mode = MPFR_RNDN;
		}

		// y = x.a;
		mpfr_set_d(y.a, x.a1, mode);
		mpfr_add_d(y.a, y.a, x.a2, mode);
	}
	template <int N> void convert(const kv::mpfr<N>& x, kv::dd& y, int rnd = 0)
	{
		kv::mpfr<N> mtmp;
		double dtmp1, dtmp2;

		mp_rnd_t mode;

		if (rnd == 1) {
			mode = MPFR_RNDU;
		}
		else if (rnd == -1) {
			mode = MPFR_RNDD;
		}
		else {
			mode = MPFR_RNDN;
		}

		convert(x, dtmp1, 0);

		// mtmp = x - dtmp1;
		// theoretically no rounding error.
		// use rounded subtraction just to be safe.
		mpfr_sub_d(mtmp.a, x.a, dtmp1, mode);

		dtmp2 = mpfr_get_d(mtmp.a, mode);

		kv::dd::twosum(dtmp1, dtmp2, y.a1, y.a2);
	}
#endif

#if defined(MPFR_HPP) && defined(RMPFR_HPP)
	template <int N, int M> void convert(const kv::mpfr<N>& x, kv::mpfr<M>& y, int rnd = 0)
	{
		mp_rnd_t mode;

		if (rnd == 1) {
			mode = MPFR_RNDU;
		}
		else if (rnd == -1) {
			mode = MPFR_RNDD;
		}
		else {
			mode = MPFR_RNDN;
		}

		mpfr_set(y.a, x.a, mode);
	}
#endif

#if defined(DD_HPP) && defined(RDD_HPP) && defined(RDOUBLE_HPP)
	void convert(const kv::dd& x, double& y, int rnd = 0)
	{
		if (rnd == 1) {
			kv::rop<double>::begin();
			y = kv::rop<double>::add_up(x.a1, x.a2);
			kv::rop<double>::end();
		}
		else if (rnd == -1) {
			kv::rop<double>::begin();
			y = kv::rop<double>::add_down(x.a1, x.a2);
			kv::rop<double>::end();
		}
		else {
			y = x.a1 + x.a2;
		}
	}
	void convert(const double& x, kv::dd& y) {
		y = x;
	}
#endif

#if defined(INTERVAL_HPP) && defined(DD_HPP) && defined(RDD_HPP) && defined(RDOUBLE_HPP)
	void convert(const kv::interval< kv::dd >& x, kv::interval<double>& y)
	{
		convert(x.lower(), y.lower(), -1);
		convert(x.upper(), y.upper(), 1);
	}
	void convert(const kv::interval<double>& x, kv::interval< kv::dd >& y)
	{
		y.lower() = x.lower();
		y.upper() = x.upper();
	}
#endif

#if defined(INTERVAL_HPP) && defined(MPFR_HPP) && defined(RMPFR_HPP) && defined(RDOUBLE_HPP)
	template <int N> void convert(const kv::interval< kv::mpfr<N> >& x, kv::interval<double>& y)
	{
		convert(x.lower(), y.lower(), -1);
		convert(x.upper(), y.upper(), 1);
	}
	template <int N> void convert(const kv::interval<double>& x, kv::interval< kv::mpfr<N> >& y)
	{
		y.lower() = x.lower();
		y.upper() = x.upper();
	}
#endif

#if defined(INTERVAL_HPP) && defined(MPFR_HPP) && defined(RMPFR_HPP) && defined(DD_HPP) && defined(RDD_HPP)
	template <int N> void convert(const kv::interval< kv::mpfr<N> >& x, kv::interval< kv::dd >& y)
	{
		convert(x.lower(), y.lower(), -1);
		convert(x.upper(), y.upper(), 1);
	}
	template <int N> void convert(const kv::interval< kv::dd >& x, kv::interval< kv::mpfr<N> >& y)
	{
		convert(x.lower(), y.lower(), -1);
		convert(x.upper(), y.upper(), 1);
	}
#endif

#if defined(INTERVAL_HPP) && defined(MPFR_HPP) && defined(RMPFR_HPP)
	template <int N, int M> void convert(const kv::interval< kv::mpfr<N> >& x, kv::interval< kv::mpfr<M> >& y)
	{
		convert(x.lower(), y.lower(), -1);
		convert(x.upper(), y.upper(), 1);
	}
#endif

#if defined(INTERVAL_HPP) && defined(RDOUBLE_HPP)
	template<typename _T> void convert(const kv::interval< _T >& x, double& y) {
		convert(mid(x), y);
	}
	template<typename _T> void convert(const double& x, kv::interval< _T >& y) {
		kv::interval< double > yy = kv::interval< double >(x);
		convert(yy, y);
	}
#endif
#if defined(INTERVAL_HPP) && defined(DD_HPP) && defined(RDD_HPP)
	template<typename _T> void convert(const kv::interval< _T >& x, kv::dd& y) {
		convert(mid(x), y);
	}
	template<typename _T> void convert(const kv::dd& x, kv::interval< _T >& y) {
		kv::interval< kv::dd > yy = kv::interval< kv::dd >(x);
		convert(yy, y);
	}
#endif
#if defined(INTERVAL_HPP) && defined(MPFR_HPP) && defined(RMPFR_HPP)
	template<typename _T, int N> void convert(const kv::interval< _T >& x, kv::mpfr< N >& y) {
		convert(mid(x), y);
	}
	template<typename _T, int N> void convert(const kv::mpfr< N >& x, kv::interval< _T >& y) {
		kv::interval< kv::mpfr< N > > yy = kv::interval< kv::mpfr< N > >(x);
		convert(yy, y);
	}
#endif

	template<typename _T> void convert(const _T& x, _T& y) {
		y = x;
	}
}
#endif //VCP_CONVERTER_HPP