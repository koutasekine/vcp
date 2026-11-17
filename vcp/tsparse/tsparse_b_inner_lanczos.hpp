// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_B_INNER_LANCZOS_HPP
#define VCP_TSPARSE_B_INNER_LANCZOS_HPP

// EIG-T1 T-6: the B-inner-product Lanczos core (b_inner_lanczos_eigs and its
// helpers b_inner_product / b_norm / b_lanczos_result) was removed.  The
// generalized symmetric problem without shift is solved by the thick-restart
// kernel on L^{-1} P A P^T L^{-T} (B = P^T L L^T P, sparse Cholesky); see
// vcp/spmats_base/spmats_eigs.hpp (t2_generalized_chol_).  This header is kept
// empty because vcp/tsparse/tsparse.hpp includes it.

#endif // VCP_TSPARSE_B_INNER_LANCZOS_HPP
