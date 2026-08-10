// vcp/bfem/constants/dict/registry_3d.hpp
//
// CONST-B2a: generated 3D class registry (dictionary first layer).  This file
// is written by sandbox/probes/constb2a_registry_gen.cpp (mode
// generate); gate G-C5 checks that the committed file and a
// regeneration agree bit for bit.  Do not edit by hand.
//
// Entry semantics: Chat^2 := C_d(K_rep)^2 / h^2(K_rep) as a
// 12-significant-digit UPWARD decimal string; consumption is
// C_d(K)^2 <= Chat^2 * h^2(K).  Keys are the canonical similarity
// keys of dict_entry.hpp as integer fractions.  rep-vertices are
// the integer-rescaled representatives (addendum section 7).
// Diagnostic components (mu, ch2) of every entry are recorded in
// the generation log referenced per entry, not here.
//
// Lexical regime B: decimal strings are authorized INSIDE the
// marker block only (cm1_ch_tests G7).

#ifndef VCP_BFEM_CONSTANTS_DICT_REGISTRY_3D_HPP
#define VCP_BFEM_CONSTANTS_DICT_REGISTRY_3D_HPP

#include <vcp/bfem/constants/dict/dict_entry.hpp>

namespace vcp {
namespace bfem {
namespace constants {

// missing (class, d) list -- ruled local scope (addendum sections 1, 8):
//   every class below is generated for d = 0..2; d = 3..8 is NOT generated locally
//   (ruled local scope; B-2b' remote scope), and d = 9 is
//   PERMANENTLY missing (engine factorial cap 2d+D <= 20; R18 amended to d <= 8):
//   class 3d-c00 missing d=3..8 and d=9
//   class 3d-c01 missing d=3..8 and d=9
//   class 3d-c02 missing d=3..8 and d=9
//   class 3d-c03 missing d=3..8 and d=9
//   class 3d-c04 missing d=3..8 and d=9
//   class 3d-c05 missing d=3..8 and d=9
//   class 3d-c06 missing d=3..8 and d=9
//   class 3d-c07 missing d=3..8 and d=9
//   class 3d-c08 missing d=3..8 and d=9
//   class 3d-c09 missing d=3..8 and d=9
//   class 3d-c10 missing d=3..8 and d=9
//   class 3d-c11 missing d=3..8 and d=9
//   class 3d-c12 missing d=3..8 and d=9
//   class 3d-c13 missing d=3..8 and d=9
//   class 3d-c14 missing d=3..8 and d=9
//   class 3d-c15 missing d=3..8 and d=9
//   class 3d-c16 missing d=3..8 and d=9
//   class 3d-c17 missing d=3..8 and d=9
//   class 3d-c18 missing d=3..8 and d=9
//   class 3d-c19 missing d=3..8 and d=9
//   class 3d-c20 missing d=3..8 and d=9
//   class 3d-c21 missing d=3..8 and d=9
//   class 3d-c22 missing d=3..8 and d=9
//   class 3d-c23 missing d=3..8 and d=9
//   class 3d-c24 missing d=3..8 and d=9
//   class 3d-c25 missing d=3..8 and d=9
//   class 3d-c26 missing d=3..8 and d=9
//   class 3d-c27 missing d=3..8 and d=9
//   class 3d-c28 missing d=3..8 and d=9
//   class 3d-c29 missing d=3..8 and d=9
//   class 3d-c30 missing d=3..8 and d=9
//   class 3d-c31 missing d=3..8 and d=9
//   class 3d-c32 missing d=3..8 and d=9
//   class 3d-c33 missing d=3..8 and d=9
//   class 3d-c34 missing d=3..8 and d=9
//   class 3d-c35 missing d=3..8 and d=9
//   class 3d-c36 missing d=3..8 and d=9
//   class 3d-c37 missing d=3..8 and d=9
//   class 3d-c38 missing d=3..8 and d=9
//   class 3d-c39 missing d=3..8 and d=9
//   class 3d-c40 missing d=3..8 and d=9
//   class 3d-c41 missing d=3..8 and d=9
//   class 3d-c42 missing d=3..8 and d=9
//   class 3d-c43 missing d=3..8 and d=9
//   class 3d-c44 missing d=3..8 and d=9
//   class 3d-c45 missing d=3..8 and d=9
//   class 3d-c46 missing d=3..8 and d=9
//   class 3d-c47 missing d=3..8 and d=9
//   class 3d-c48 missing d=3..8 and d=9
//   class 3d-c49 missing d=3..8 and d=9
//   class 3d-c50 missing d=3..8 and d=9

// VCP_CONSTANTS_TABLE_BEGIN  (authorized decimal strings; generated table)

namespace detail {

inline const l2_projection_registry_entry_3d* l2_projection_registry_3d_entries(int& count) {
    static const l2_projection_registry_entry_3d entries[] = {
        // class 3d-c00 key=[1/1 1/1 1/1 1/1 1/1 1/1]
        //   rep-vertices: (1,0,0) (0,1,0) (0,0,1) (1,1,1)  prov=seed:T2
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c00-d0, alpha-achieved=0.671 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          0, "0.0675446739209" },
        // class 3d-c00 key=[1/1 1/1 1/1 1/1 1/1 1/1]
        //   rep-vertices: (1,0,0) (0,1,0) (0,0,1) (1,1,1)  prov=seed:T2
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c00-d1, alpha-achieved=1.56 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          1, "0.0445350335374" },
        // class 3d-c00 key=[1/1 1/1 1/1 1/1 1/1 1/1]
        //   rep-vertices: (1,0,0) (0,1,0) (0,0,1) (1,1,1)  prov=seed:T2
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c00-d2, alpha-achieved=2.83 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          2, "0.0367222558805" },
        // class 3d-c01 key=[1/1 1/1 1/1 1/1 1/1 2/1]
        //   rep-vertices: (1,1,0) (1,0,1) (2,1,1) (1,1,2)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c01-d0, alpha-achieved=0.69 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 1LL, 1LL, 2LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          0, "0.0553559422224" },
        // class 3d-c01 key=[1/1 1/1 1/1 1/1 1/1 2/1]
        //   rep-vertices: (1,1,0) (1,0,1) (2,1,1) (1,1,2)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c01-d1, alpha-achieved=1.71 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 1LL, 1LL, 2LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          1, "0.0358131546554" },
        // class 3d-c01 key=[1/1 1/1 1/1 1/1 1/1 2/1]
        //   rep-vertices: (1,1,0) (1,0,1) (2,1,1) (1,1,2)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c01-d2, alpha-achieved=3.22 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 1LL, 1LL, 2LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          2, "0.0296241912944" },
        // class 3d-c02 key=[1/1 1/1 1/1 16/9 16/9 32/9]
        //   rep-vertices: (1,2,2) (0,4,0) (0,0,4) (0,4,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c02-d0, alpha-achieved=0.739 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 16LL, 16LL, 32LL },
          { 1LL, 1LL, 1LL, 9LL, 9LL, 9LL },
          0, "0.0532094653345" },
        // class 3d-c02 key=[1/1 1/1 1/1 16/9 16/9 32/9]
        //   rep-vertices: (1,2,2) (0,4,0) (0,0,4) (0,4,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c02-d1, alpha-achieved=1.85 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 16LL, 16LL, 32LL },
          { 1LL, 1LL, 1LL, 9LL, 9LL, 9LL },
          1, "0.0348485216551" },
        // class 3d-c02 key=[1/1 1/1 1/1 16/9 16/9 32/9]
        //   rep-vertices: (1,2,2) (0,4,0) (0,0,4) (0,4,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c02-d2, alpha-achieved=3.44 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 16LL, 16LL, 32LL },
          { 1LL, 1LL, 1LL, 9LL, 9LL, 9LL },
          2, "0.0291773100175" },
        // class 3d-c03 key=[1/1 1/1 1/1 2/1 2/1 2/1]
        //   rep-vertices: (0,0,0) (1,0,0) (0,1,0) (0,0,1)  prov=seed:T1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c03-d0, alpha-achieved=0.737 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 2LL, 2LL, 2LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          0, "0.0639498452441" },
        // class 3d-c03 key=[1/1 1/1 1/1 2/1 2/1 2/1]
        //   rep-vertices: (0,0,0) (1,0,0) (0,1,0) (0,0,1)  prov=seed:T1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c03-d1, alpha-achieved=1.77 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 2LL, 2LL, 2LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          1, "0.0424665512209" },
        // class 3d-c03 key=[1/1 1/1 1/1 2/1 2/1 2/1]
        //   rep-vertices: (0,0,0) (1,0,0) (0,1,0) (0,0,1)  prov=seed:T1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c03-d2, alpha-achieved=3.21 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 2LL, 2LL, 2LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          2, "0.0355866107061" },
        // class 3d-c04 key=[1/1 1/1 1/1 8/3 8/3 8/3]
        //   rep-vertices: (2,0,0) (0,2,0) (0,0,2) (1,1,1)  prov=seed:T5
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c04-d0, alpha-achieved=0.756 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 8LL, 8LL, 8LL },
          { 1LL, 1LL, 1LL, 3LL, 3LL, 3LL },
          0, "0.0630196245480" },
        // class 3d-c04 key=[1/1 1/1 1/1 8/3 8/3 8/3]
        //   rep-vertices: (2,0,0) (0,2,0) (0,0,2) (1,1,1)  prov=seed:T5
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c04-d1, alpha-achieved=1.82 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 8LL, 8LL, 8LL },
          { 1LL, 1LL, 1LL, 3LL, 3LL, 3LL },
          1, "0.0420605630435" },
        // class 3d-c04 key=[1/1 1/1 1/1 8/3 8/3 8/3]
        //   rep-vertices: (2,0,0) (0,2,0) (0,0,2) (1,1,1)  prov=seed:T5
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c04-d2, alpha-achieved=3.28 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 8LL, 8LL, 8LL },
          { 1LL, 1LL, 1LL, 3LL, 3LL, 3LL },
          2, "0.0354141272280" },
        // class 3d-c05 key=[1/1 1/1 1/1 32/11 32/11 32/11]
        //   rep-vertices: (4,0,0) (0,4,0) (0,0,4) (1,1,1)  prov=seed:T4
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c05-d0, alpha-achieved=0.761 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 32LL, 32LL, 32LL },
          { 1LL, 1LL, 1LL, 11LL, 11LL, 11LL },
          0, "0.0627649611432" },
        // class 3d-c05 key=[1/1 1/1 1/1 32/11 32/11 32/11]
        //   rep-vertices: (4,0,0) (0,4,0) (0,0,4) (1,1,1)  prov=seed:T4
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c05-d1, alpha-achieved=1.83 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 32LL, 32LL, 32LL },
          { 1LL, 1LL, 1LL, 11LL, 11LL, 11LL },
          1, "0.0419522809243" },
        // class 3d-c05 key=[1/1 1/1 1/1 32/11 32/11 32/11]
        //   rep-vertices: (4,0,0) (0,4,0) (0,0,4) (1,1,1)  prov=seed:T4
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c05-d2, alpha-achieved=3.29 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 32LL, 32LL, 32LL },
          { 1LL, 1LL, 1LL, 11LL, 11LL, 11LL },
          2, "0.0353710885162" },
        // class 3d-c06 key=[1/1 1/1 1/1 128/43 128/43 128/43]
        //   rep-vertices: (8,0,0) (0,8,0) (0,0,8) (3,3,3)  prov=seed:T5>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c06-d0, alpha-achieved=0.763 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 128LL, 128LL, 128LL },
          { 1LL, 1LL, 1LL, 43LL, 43LL, 43LL },
          0, "0.0626989630047" },
        // class 3d-c06 key=[1/1 1/1 1/1 128/43 128/43 128/43]
        //   rep-vertices: (8,0,0) (0,8,0) (0,0,8) (3,3,3)  prov=seed:T5>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c06-d1, alpha-achieved=1.83 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 128LL, 128LL, 128LL },
          { 1LL, 1LL, 1LL, 43LL, 43LL, 43LL },
          1, "0.0419245947266" },
        // class 3d-c06 key=[1/1 1/1 1/1 128/43 128/43 128/43]
        //   rep-vertices: (8,0,0) (0,8,0) (0,0,8) (3,3,3)  prov=seed:T5>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c06-d2, alpha-achieved=3.3 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 128LL, 128LL, 128LL },
          { 1LL, 1LL, 1LL, 43LL, 43LL, 43LL },
          2, "0.0353602764271" },
        // class 3d-c07 key=[1/1 1/1 1/1 512/171 512/171 512/171]
        //   rep-vertices: (16,0,0) (0,16,0) (0,0,16) (5,5,5)  prov=seed:T4>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c07-d0, alpha-achieved=0.763 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 512LL, 512LL, 512LL },
          { 1LL, 1LL, 1LL, 171LL, 171LL, 171LL },
          0, "0.0626823041136" },
        // class 3d-c07 key=[1/1 1/1 1/1 512/171 512/171 512/171]
        //   rep-vertices: (16,0,0) (0,16,0) (0,0,16) (5,5,5)  prov=seed:T4>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c07-d1, alpha-achieved=1.84 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 512LL, 512LL, 512LL },
          { 1LL, 1LL, 1LL, 171LL, 171LL, 171LL },
          1, "0.0419176305087" },
        // class 3d-c07 key=[1/1 1/1 1/1 512/171 512/171 512/171]
        //   rep-vertices: (16,0,0) (0,16,0) (0,0,16) (5,5,5)  prov=seed:T4>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c07-d2, alpha-achieved=3.3 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 1LL, 512LL, 512LL, 512LL },
          { 1LL, 1LL, 1LL, 171LL, 171LL, 171LL },
          2, "0.0353575691394" },
        // class 3d-c08 key=[1/1 1/1 9/5 16/5 16/5 16/5]
        //   rep-vertices: (5,3,4) (4,0,4) (8,4,4) (4,4,8)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c08-d0, alpha-achieved=0.715 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 9LL, 16LL, 16LL, 16LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          0, "0.0650820060024" },
        // class 3d-c08 key=[1/1 1/1 9/5 16/5 16/5 16/5]
        //   rep-vertices: (5,3,4) (4,0,4) (8,4,4) (4,4,8)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c08-d1, alpha-achieved=1.78 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 9LL, 16LL, 16LL, 16LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          1, "0.0424040794495" },
        // class 3d-c08 key=[1/1 1/1 9/5 16/5 16/5 16/5]
        //   rep-vertices: (5,3,4) (4,0,4) (8,4,4) (4,4,8)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c08-d2, alpha-achieved=3.19 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 9LL, 16LL, 16LL, 16LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          2, "0.0356449356585" },
        // class 3d-c09 key=[1/1 1/1 17/9 16/9 32/9 16/3]
        //   rep-vertices: (4,0,0) (1,2,2) (0,0,4) (0,4,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c09-d0, alpha-achieved=0.814 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 17LL, 16LL, 32LL, 16LL },
          { 1LL, 1LL, 9LL, 9LL, 9LL, 3LL },
          0, "0.0604496333200" },
        // class 3d-c09 key=[1/1 1/1 17/9 16/9 32/9 16/3]
        //   rep-vertices: (4,0,0) (1,2,2) (0,0,4) (0,4,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c09-d1, alpha-achieved=2.02 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 17LL, 16LL, 32LL, 16LL },
          { 1LL, 1LL, 9LL, 9LL, 9LL, 3LL },
          1, "0.0405942666493" },
        // class 3d-c09 key=[1/1 1/1 17/9 16/9 32/9 16/3]
        //   rep-vertices: (4,0,0) (1,2,2) (0,0,4) (0,4,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c09-d2, alpha-achieved=3.45 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 17LL, 16LL, 32LL, 16LL },
          { 1LL, 1LL, 9LL, 9LL, 9LL, 3LL },
          2, "0.0350077970905" },
        // class 3d-c10 key=[1/1 1/1 17/9 32/9 32/9 32/9]
        //   rep-vertices: (4,0,0) (0,4,0) (0,0,4) (1,2,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c10-d0, alpha-achieved=0.721 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 17LL, 32LL, 32LL, 32LL },
          { 1LL, 1LL, 9LL, 9LL, 9LL, 9LL },
          0, "0.0647773900334" },
        // class 3d-c10 key=[1/1 1/1 17/9 32/9 32/9 32/9]
        //   rep-vertices: (4,0,0) (0,4,0) (0,0,4) (1,2,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c10-d1, alpha-achieved=1.79 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 17LL, 32LL, 32LL, 32LL },
          { 1LL, 1LL, 9LL, 9LL, 9LL, 9LL },
          1, "0.0422774732119" },
        // class 3d-c10 key=[1/1 1/1 17/9 32/9 32/9 32/9]
        //   rep-vertices: (4,0,0) (0,4,0) (0,0,4) (1,2,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c10-d2, alpha-achieved=3.22 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 17LL, 32LL, 32LL, 32LL },
          { 1LL, 1LL, 9LL, 9LL, 9LL, 9LL },
          2, "0.0355647057588" },
        // class 3d-c11 key=[1/1 1/1 67/35 128/35 128/35 128/35]
        //   rep-vertices: (7,7,5) (12,4,4) (4,12,4) (4,4,12)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c11-d0, alpha-achieved=0.722 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 67LL, 128LL, 128LL, 128LL },
          { 1LL, 1LL, 35LL, 35LL, 35LL, 35LL },
          0, "0.0646989784707" },
        // class 3d-c11 key=[1/1 1/1 67/35 128/35 128/35 128/35]
        //   rep-vertices: (7,7,5) (12,4,4) (4,12,4) (4,4,12)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c11-d1, alpha-achieved=1.8 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 67LL, 128LL, 128LL, 128LL },
          { 1LL, 1LL, 35LL, 35LL, 35LL, 35LL },
          1, "0.0422453569362" },
        // class 3d-c11 key=[1/1 1/1 67/35 128/35 128/35 128/35]
        //   rep-vertices: (7,7,5) (12,4,4) (4,12,4) (4,4,12)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c11-d2, alpha-achieved=3.22 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 67LL, 128LL, 128LL, 128LL },
          { 1LL, 1LL, 35LL, 35LL, 35LL, 35LL },
          2, "0.0355453568883" },
        // class 3d-c12 key=[1/1 1/1 267/139 512/139 512/139 512/139]
        //   rep-vertices: (11,11,7) (20,4,4) (4,20,4) (4,4,20)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c12-d0, alpha-achieved=0.723 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 267LL, 512LL, 512LL, 512LL },
          { 1LL, 1LL, 139LL, 139LL, 139LL, 139LL },
          0, "0.0646792235484" },
        // class 3d-c12 key=[1/1 1/1 267/139 512/139 512/139 512/139]
        //   rep-vertices: (11,11,7) (20,4,4) (4,20,4) (4,4,20)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c12-d1, alpha-achieved=1.8 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 267LL, 512LL, 512LL, 512LL },
          { 1LL, 1LL, 139LL, 139LL, 139LL, 139LL },
          1, "0.0422372953064" },
        // class 3d-c12 key=[1/1 1/1 267/139 512/139 512/139 512/139]
        //   rep-vertices: (11,11,7) (20,4,4) (4,20,4) (4,4,20)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c12-d2, alpha-achieved=3.23 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 267LL, 512LL, 512LL, 512LL },
          { 1LL, 1LL, 139LL, 139LL, 139LL, 139LL },
          2, "0.0355405672737" },
        // class 3d-c13 key=[1/1 1/1 2/1 2/1 1/1 3/1]
        //   rep-vertices: (0,0,0) (1,0,0) (1,1,0) (1,1,1)  prov=seed:Kuhn0
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c13-d0, alpha-achieved=0.768 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 2LL, 2LL, 1LL, 3LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          0, "0.0624662949174" },
        // class 3d-c13 key=[1/1 1/1 2/1 2/1 1/1 3/1]
        //   rep-vertices: (0,0,0) (1,0,0) (1,1,0) (1,1,1)  prov=seed:Kuhn0
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c13-d1, alpha-achieved=1.88 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 2LL, 2LL, 1LL, 3LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          1, "0.0415495813406" },
        // class 3d-c13 key=[1/1 1/1 2/1 2/1 1/1 3/1]
        //   rep-vertices: (0,0,0) (1,0,0) (1,1,0) (1,1,1)  prov=seed:Kuhn0
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c13-d2, alpha-achieved=3.33 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 2LL, 2LL, 1LL, 3LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          2, "0.0352754854574" },
        // class 3d-c14 key=[1/1 1/1 7/3 8/3 8/3 16/3]
        //   rep-vertices: (3,2,1) (4,0,0) (4,4,0) (4,4,4)  prov=seed:Kuhn0>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c14-d0, alpha-achieved=0.756 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 7LL, 8LL, 8LL, 16LL },
          { 1LL, 1LL, 3LL, 3LL, 3LL, 3LL },
          0, "0.0564632696277" },
        // class 3d-c14 key=[1/1 1/1 7/3 8/3 8/3 16/3]
        //   rep-vertices: (3,2,1) (4,0,0) (4,4,0) (4,4,4)  prov=seed:Kuhn0>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c14-d1, alpha-achieved=1.95 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 7LL, 8LL, 8LL, 16LL },
          { 1LL, 1LL, 3LL, 3LL, 3LL, 3LL },
          1, "0.0367634601106" },
        // class 3d-c14 key=[1/1 1/1 7/3 8/3 8/3 16/3]
        //   rep-vertices: (3,2,1) (4,0,0) (4,4,0) (4,4,4)  prov=seed:Kuhn0>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c14-d2, alpha-achieved=3.46 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 7LL, 8LL, 8LL, 16LL },
          { 1LL, 1LL, 3LL, 3LL, 3LL, 3LL },
          2, "0.0313255739818" },
        // class 3d-c15 key=[1/1 1/1 13/5 12/5 12/5 32/5]
        //   rep-vertices: (4,2,3) (4,0,4) (6,2,2) (2,2,6)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c15-d0, alpha-achieved=0.434 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 13LL, 12LL, 12LL, 32LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          0, "0.0420256141872" },
        // class 3d-c15 key=[1/1 1/1 13/5 12/5 12/5 32/5]
        //   rep-vertices: (4,2,3) (4,0,4) (6,2,2) (2,2,6)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c15-d1, alpha-achieved=1.23 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 13LL, 12LL, 12LL, 32LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          1, "0.0230387657186" },
        // class 3d-c15 key=[1/1 1/1 13/5 12/5 12/5 32/5]
        //   rep-vertices: (4,2,3) (4,0,4) (6,2,2) (2,2,6)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c15-d2, alpha-achieved=2.39 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 13LL, 12LL, 12LL, 32LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          2, "0.0180314135000" },
        // class 3d-c16 key=[1/1 1/1 8/3 8/3 1/1 11/3]
        //   rep-vertices: (2,2,0) (2,0,2) (3,1,1) (1,1,3)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c16-d0, alpha-achieved=0.461 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 8LL, 8LL, 1LL, 11LL },
          { 1LL, 1LL, 3LL, 3LL, 1LL, 3LL },
          0, "0.0521371539472" },
        // class 3d-c16 key=[1/1 1/1 8/3 8/3 1/1 11/3]
        //   rep-vertices: (2,2,0) (2,0,2) (3,1,1) (1,1,3)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c16-d1, alpha-achieved=1.29 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 8LL, 8LL, 1LL, 11LL },
          { 1LL, 1LL, 3LL, 3LL, 1LL, 3LL },
          1, "0.0292349118464" },
        // class 3d-c16 key=[1/1 1/1 8/3 8/3 1/1 11/3]
        //   rep-vertices: (2,2,0) (2,0,2) (3,1,1) (1,1,3)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c16-d2, alpha-achieved=2.39 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 8LL, 8LL, 1LL, 11LL },
          { 1LL, 1LL, 3LL, 3LL, 1LL, 3LL },
          2, "0.0233274447952" },
        // class 3d-c17 key=[1/1 1/1 51/19 44/19 44/19 128/19]
        //   rep-vertices: (7,3,5) (8,0,8) (10,2,2) (2,2,10)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c17-d0, alpha-achieved=0.42 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 51LL, 44LL, 44LL, 128LL },
          { 1LL, 1LL, 19LL, 19LL, 19LL, 19LL },
          0, "0.0408333371180" },
        // class 3d-c17 key=[1/1 1/1 51/19 44/19 44/19 128/19]
        //   rep-vertices: (7,3,5) (8,0,8) (10,2,2) (2,2,10)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c17-d1, alpha-achieved=1.2 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 51LL, 44LL, 44LL, 128LL },
          { 1LL, 1LL, 19LL, 19LL, 19LL, 19LL },
          1, "0.0221295021625" },
        // class 3d-c17 key=[1/1 1/1 51/19 44/19 44/19 128/19]
        //   rep-vertices: (7,3,5) (8,0,8) (10,2,2) (2,2,10)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c17-d2, alpha-achieved=2.36 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 51LL, 44LL, 44LL, 128LL },
          { 1LL, 1LL, 19LL, 19LL, 19LL, 19LL },
          2, "0.0171916588839" },
        // class 3d-c18 key=[1/1 1/1 32/11 32/11 1/1 43/11]
        //   rep-vertices: (4,4,0) (4,0,4) (5,1,1) (1,1,5)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c18-d0, alpha-achieved=0.446 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 32LL, 32LL, 1LL, 43LL },
          { 1LL, 1LL, 11LL, 11LL, 1LL, 11LL },
          0, "0.0517884360775" },
        // class 3d-c18 key=[1/1 1/1 32/11 32/11 1/1 43/11]
        //   rep-vertices: (4,4,0) (4,0,4) (5,1,1) (1,1,5)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c18-d1, alpha-achieved=1.25 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 32LL, 32LL, 1LL, 43LL },
          { 1LL, 1LL, 11LL, 11LL, 1LL, 11LL },
          1, "0.0287600481178" },
        // class 3d-c18 key=[1/1 1/1 32/11 32/11 1/1 43/11]
        //   rep-vertices: (4,4,0) (4,0,4) (5,1,1) (1,1,5)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c18-d2, alpha-achieved=2.31 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 32LL, 32LL, 1LL, 43LL },
          { 1LL, 1LL, 11LL, 11LL, 1LL, 11LL },
          2, "0.0228915328049" },
        // class 3d-c19 key=[1/1 1/1 3/1 2/1 2/1 2/1]
        //   rep-vertices: (1,0,0) (0,1,0) (0,0,1) (0,1,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c19-d0, alpha-achieved=0.759 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 3LL, 2LL, 2LL, 2LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          0, "0.0628787005653" },
        // class 3d-c19 key=[1/1 1/1 3/1 2/1 2/1 2/1]
        //   rep-vertices: (1,0,0) (0,1,0) (0,0,1) (0,1,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c19-d1, alpha-achieved=1.88 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 3LL, 2LL, 2LL, 2LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          1, "0.0415884692727" },
        // class 3d-c19 key=[1/1 1/1 3/1 2/1 2/1 2/1]
        //   rep-vertices: (1,0,0) (0,1,0) (0,0,1) (0,1,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c19-d2, alpha-achieved=3.21 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 3LL, 2LL, 2LL, 2LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          2, "0.0355733393130" },
        // class 3d-c20 key=[1/1 1/1 3/1 2/1 2/1 6/1]
        //   rep-vertices: (1,0,0) (1,1,0) (1,1,1) (2,2,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c20-d0, alpha-achieved=0.825 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 3LL, 2LL, 2LL, 6LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          0, "0.0566887212144" },
        // class 3d-c20 key=[1/1 1/1 3/1 2/1 2/1 6/1]
        //   rep-vertices: (1,0,0) (1,1,0) (1,1,1) (2,2,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c20-d1, alpha-achieved=2.04 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 3LL, 2LL, 2LL, 6LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          1, "0.0382036689946" },
        // class 3d-c20 key=[1/1 1/1 3/1 2/1 2/1 6/1]
        //   rep-vertices: (1,0,0) (1,1,0) (1,1,1) (2,2,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c20-d2, alpha-achieved=3.55 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 3LL, 2LL, 2LL, 6LL },
          { 1LL, 1LL, 1LL, 1LL, 1LL, 1LL },
          2, "0.0328470152870" },
        // class 3d-c21 key=[1/1 1/1 11/3 8/3 8/3 8/3]
        //   rep-vertices: (2,2,0) (3,1,1) (1,3,1) (1,1,3)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c21-d0, alpha-achieved=0.762 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 11LL, 8LL, 8LL, 8LL },
          { 1LL, 1LL, 3LL, 3LL, 3LL, 3LL },
          0, "0.0665233561093" },
        // class 3d-c21 key=[1/1 1/1 11/3 8/3 8/3 8/3]
        //   rep-vertices: (2,2,0) (3,1,1) (1,3,1) (1,1,3)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c21-d1, alpha-achieved=1.89 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 11LL, 8LL, 8LL, 8LL },
          { 1LL, 1LL, 3LL, 3LL, 3LL, 3LL },
          1, "0.0439804697226" },
        // class 3d-c21 key=[1/1 1/1 11/3 8/3 8/3 8/3]
        //   rep-vertices: (2,2,0) (3,1,1) (1,3,1) (1,1,3)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c21-d2, alpha-achieved=3.21 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 11LL, 8LL, 8LL, 8LL },
          { 1LL, 1LL, 3LL, 3LL, 3LL, 3LL },
          2, "0.0377520618582" },
        // class 3d-c22 key=[1/1 1/1 43/11 32/11 32/11 32/11]
        //   rep-vertices: (4,4,0) (5,1,1) (1,5,1) (1,1,5)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c22-d0, alpha-achieved=0.763 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 43LL, 32LL, 32LL, 32LL },
          { 1LL, 1LL, 11LL, 11LL, 11LL, 11LL },
          0, "0.0675725685465" },
        // class 3d-c22 key=[1/1 1/1 43/11 32/11 32/11 32/11]
        //   rep-vertices: (4,4,0) (5,1,1) (1,5,1) (1,1,5)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c22-d1, alpha-achieved=1.9 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 43LL, 32LL, 32LL, 32LL },
          { 1LL, 1LL, 11LL, 11LL, 11LL, 11LL },
          1, "0.0446571019241" },
        // class 3d-c22 key=[1/1 1/1 43/11 32/11 32/11 32/11]
        //   rep-vertices: (4,4,0) (5,1,1) (1,5,1) (1,1,5)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c22-d2, alpha-achieved=3.2 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 43LL, 32LL, 32LL, 32LL },
          { 1LL, 1LL, 11LL, 11LL, 11LL, 11LL },
          2, "0.0383589754272" },
        // class 3d-c23 key=[1/1 1/1 21/5 16/5 16/5 32/5]
        //   rep-vertices: (4,0,0) (4,4,0) (4,4,4) (5,4,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c23-d0, alpha-achieved=0.801 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 21LL, 16LL, 16LL, 32LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          0, "0.0641920830110" },
        // class 3d-c23 key=[1/1 1/1 21/5 16/5 16/5 32/5]
        //   rep-vertices: (4,0,0) (4,4,0) (4,4,4) (5,4,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c23-d1, alpha-achieved=2.03 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 21LL, 16LL, 16LL, 32LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          1, "0.0426298325755" },
        // class 3d-c23 key=[1/1 1/1 21/5 16/5 16/5 32/5]
        //   rep-vertices: (4,0,0) (4,4,0) (4,4,4) (5,4,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c23-d2, alpha-achieved=3.4 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 21LL, 16LL, 16LL, 32LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          2, "0.0369316086370" },
        // class 3d-c24 key=[1/1 1/1 29/5 16/5 32/5 48/5]
        //   rep-vertices: (5,4,2) (4,4,0) (4,4,4) (8,8,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c24-d0, alpha-achieved=0.753 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 29LL, 16LL, 32LL, 48LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          0, "0.0653372450722" },
        // class 3d-c24 key=[1/1 1/1 29/5 16/5 32/5 48/5]
        //   rep-vertices: (5,4,2) (4,4,0) (4,4,4) (8,8,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c24-d1, alpha-achieved=1.95 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 29LL, 16LL, 32LL, 48LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          1, "0.0425062628827" },
        // class 3d-c24 key=[1/1 1/1 29/5 16/5 32/5 48/5]
        //   rep-vertices: (5,4,2) (4,4,0) (4,4,4) (8,8,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c24-d2, alpha-achieved=3.34 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 1LL, 29LL, 16LL, 32LL, 48LL },
          { 1LL, 1LL, 5LL, 5LL, 5LL, 5LL },
          2, "0.0364906759281" },
        // class 3d-c25 key=[1/1 35/27 35/27 16/9 16/9 128/27]
        //   rep-vertices: (8,8,0) (12,4,4) (4,12,4) (7,7,5)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c25-d0, alpha-achieved=0.74 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 35LL, 35LL, 16LL, 16LL, 128LL },
          { 1LL, 27LL, 27LL, 9LL, 9LL, 27LL },
          0, "0.0505180645364" },
        // class 3d-c25 key=[1/1 35/27 35/27 16/9 16/9 128/27]
        //   rep-vertices: (8,8,0) (12,4,4) (4,12,4) (7,7,5)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c25-d1, alpha-achieved=1.87 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 35LL, 35LL, 16LL, 16LL, 128LL },
          { 1LL, 27LL, 27LL, 9LL, 9LL, 27LL },
          1, "0.0329835205660" },
        // class 3d-c25 key=[1/1 35/27 35/27 16/9 16/9 128/27]
        //   rep-vertices: (8,8,0) (12,4,4) (4,12,4) (7,7,5)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c25-d2, alpha-achieved=3.5 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 35LL, 35LL, 16LL, 16LL, 128LL },
          { 1LL, 27LL, 27LL, 9LL, 9LL, 27LL },
          2, "0.0276237848182" },
        // class 3d-c26 key=[1/1 35/27 67/27 16/9 176/27 128/27]
        //   rep-vertices: (8,8,0) (7,7,5) (4,12,4) (4,4,12)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c26-d0, alpha-achieved=0.689 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 35LL, 67LL, 16LL, 176LL, 128LL },
          { 1LL, 27LL, 27LL, 9LL, 27LL, 27LL },
          0, "0.0578414342760" },
        // class 3d-c26 key=[1/1 35/27 67/27 16/9 176/27 128/27]
        //   rep-vertices: (8,8,0) (7,7,5) (4,12,4) (4,4,12)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c26-d1, alpha-achieved=1.8 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 35LL, 67LL, 16LL, 176LL, 128LL },
          { 1LL, 27LL, 27LL, 9LL, 27LL, 27LL },
          1, "0.0367259427381" },
        // class 3d-c26 key=[1/1 35/27 67/27 16/9 176/27 128/27]
        //   rep-vertices: (8,8,0) (7,7,5) (4,12,4) (4,4,12)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c26-d2, alpha-achieved=3.25 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 35LL, 67LL, 16LL, 176LL, 128LL },
          { 1LL, 27LL, 27LL, 9LL, 27LL, 27LL },
          2, "0.0308405902046" },
        // class 3d-c27 key=[1/1 139/99 139/99 16/9 16/9 512/99]
        //   rep-vertices: (16,16,0) (20,4,4) (4,20,4) (11,11,7)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c27-d0, alpha-achieved=0.742 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 139LL, 139LL, 16LL, 16LL, 512LL },
          { 1LL, 99LL, 99LL, 9LL, 9LL, 99LL },
          0, "0.0497599246054" },
        // class 3d-c27 key=[1/1 139/99 139/99 16/9 16/9 512/99]
        //   rep-vertices: (16,16,0) (20,4,4) (4,20,4) (11,11,7)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c27-d1, alpha-achieved=1.88 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 139LL, 139LL, 16LL, 16LL, 512LL },
          { 1LL, 99LL, 99LL, 9LL, 9LL, 99LL },
          1, "0.0324893859185" },
        // class 3d-c27 key=[1/1 139/99 139/99 16/9 16/9 512/99]
        //   rep-vertices: (16,16,0) (20,4,4) (4,20,4) (11,11,7)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c27-d2, alpha-achieved=3.5 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 139LL, 139LL, 16LL, 16LL, 512LL },
          { 1LL, 99LL, 99LL, 9LL, 9LL, 99LL },
          2, "0.0272549443641" },
        // class 3d-c28 key=[1/1 139/99 89/33 16/9 688/99 512/99]
        //   rep-vertices: (16,16,0) (11,11,7) (4,20,4) (4,4,20)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c28-d0, alpha-achieved=0.686 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 139LL, 89LL, 16LL, 688LL, 512LL },
          { 1LL, 99LL, 33LL, 9LL, 99LL, 99LL },
          0, "0.0584961939272" },
        // class 3d-c28 key=[1/1 139/99 89/33 16/9 688/99 512/99]
        //   rep-vertices: (16,16,0) (11,11,7) (4,20,4) (4,4,20)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c28-d1, alpha-achieved=1.79 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 139LL, 89LL, 16LL, 688LL, 512LL },
          { 1LL, 99LL, 33LL, 9LL, 99LL, 99LL },
          1, "0.0370952419644" },
        // class 3d-c28 key=[1/1 139/99 89/33 16/9 688/99 512/99]
        //   rep-vertices: (16,16,0) (11,11,7) (4,20,4) (4,4,20)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c28-d2, alpha-achieved=3.24 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 139LL, 89LL, 16LL, 688LL, 512LL },
          { 1LL, 99LL, 33LL, 9LL, 99LL, 99LL },
          2, "0.0311573579654" },
        // class 3d-c29 key=[1/1 83/59 219/59 176/59 256/59 560/59]
        //   rep-vertices: (16,0,0) (4,4,4) (11,7,3) (4,20,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c29-d0, alpha-achieved=0.762 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 83LL, 219LL, 176LL, 256LL, 560LL },
          { 1LL, 59LL, 59LL, 59LL, 59LL, 59LL },
          0, "0.0543753538744" },
        // class 3d-c29 key=[1/1 83/59 219/59 176/59 256/59 560/59]
        //   rep-vertices: (16,0,0) (4,4,4) (11,7,3) (4,20,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c29-d1, alpha-achieved=1.95 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 83LL, 219LL, 176LL, 256LL, 560LL },
          { 1LL, 59LL, 59LL, 59LL, 59LL, 59LL },
          1, "0.0355546812822" },
        // class 3d-c29 key=[1/1 83/59 219/59 176/59 256/59 560/59]
        //   rep-vertices: (16,0,0) (4,4,4) (11,7,3) (4,20,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c29-d2, alpha-achieved=3.48 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 83LL, 219LL, 176LL, 256LL, 560LL },
          { 1LL, 59LL, 59LL, 59LL, 59LL, 59LL },
          2, "0.0302791612532" },
        // class 3d-c30 key=[1/1 91/59 219/59 256/59 256/59 512/59]
        //   rep-vertices: (11,7,3) (4,4,4) (20,4,4) (4,20,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c30-d0, alpha-achieved=0.5 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 91LL, 219LL, 256LL, 256LL, 512LL },
          { 1LL, 59LL, 59LL, 59LL, 59LL, 59LL },
          0, "0.0471953399278" },
        // class 3d-c30 key=[1/1 91/59 219/59 256/59 256/59 512/59]
        //   rep-vertices: (11,7,3) (4,4,4) (20,4,4) (4,20,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c30-d1, alpha-achieved=1.38 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 91LL, 219LL, 256LL, 256LL, 512LL },
          { 1LL, 59LL, 59LL, 59LL, 59LL, 59LL },
          1, "0.0271417753789" },
        // class 3d-c30 key=[1/1 91/59 219/59 256/59 256/59 512/59]
        //   rep-vertices: (11,7,3) (4,4,4) (20,4,4) (4,20,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c30-d2, alpha-achieved=2.57 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 91LL, 219LL, 256LL, 256LL, 512LL },
          { 1LL, 59LL, 59LL, 59LL, 59LL, 59LL },
          2, "0.0218596016896" },
        // class 3d-c31 key=[1/1 19/12 32/3 9/4 35/3 17/4]
        //   rep-vertices: (8,0,0) (0,8,0) (3,5,1) (2,10,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c31-d0, alpha-achieved=0.61 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 19LL, 32LL, 9LL, 35LL, 17LL },
          { 1LL, 12LL, 3LL, 4LL, 3LL, 4LL },
          0, "0.0604980535090" },
        // class 3d-c31 key=[1/1 19/12 32/3 9/4 35/3 17/4]
        //   rep-vertices: (8,0,0) (0,8,0) (3,5,1) (2,10,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c31-d1, alpha-achieved=1.66 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 19LL, 32LL, 9LL, 35LL, 17LL },
          { 1LL, 12LL, 3LL, 4LL, 3LL, 4LL },
          1, "0.0367796043252" },
        // class 3d-c31 key=[1/1 19/12 32/3 9/4 35/3 17/4]
        //   rep-vertices: (8,0,0) (0,8,0) (3,5,1) (2,10,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c31-d2, alpha-achieved=3.05 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 19LL, 32LL, 9LL, 35LL, 17LL },
          { 1LL, 12LL, 3LL, 4LL, 3LL, 4LL },
          2, "0.0304589702622" },
        // class 3d-c32 key=[1/1 19/11 27/11 4/1 64/11 12/11]
        //   rep-vertices: (3,5,1) (0,8,0) (2,2,2) (2,10,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c32-d0, alpha-achieved=0.392 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 19LL, 27LL, 4LL, 64LL, 12LL },
          { 1LL, 11LL, 11LL, 1LL, 11LL, 11LL },
          0, "0.0457052741972" },
        // class 3d-c32 key=[1/1 19/11 27/11 4/1 64/11 12/11]
        //   rep-vertices: (3,5,1) (0,8,0) (2,2,2) (2,10,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c32-d1, alpha-achieved=1.13 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 19LL, 27LL, 4LL, 64LL, 12LL },
          { 1LL, 11LL, 11LL, 1LL, 11LL, 11LL },
          1, "0.0242735695654" },
        // class 3d-c32 key=[1/1 19/11 27/11 4/1 64/11 12/11]
        //   rep-vertices: (3,5,1) (0,8,0) (2,2,2) (2,10,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c32-d2, alpha-achieved=2.16 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 19LL, 27LL, 4LL, 64LL, 12LL },
          { 1LL, 11LL, 11LL, 1LL, 11LL, 11LL },
          2, "0.0188151386119" },
        // class 3d-c33 key=[1/1 19/11 51/11 4/1 4/1 128/11]
        //   rep-vertices: (8,0,0) (0,8,0) (2,2,2) (3,5,1)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c33-d0, alpha-achieved=0.435 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 19LL, 51LL, 4LL, 4LL, 128LL },
          { 1LL, 11LL, 11LL, 1LL, 1LL, 11LL },
          0, "0.0417124107601" },
        // class 3d-c33 key=[1/1 19/11 51/11 4/1 4/1 128/11]
        //   rep-vertices: (8,0,0) (0,8,0) (2,2,2) (3,5,1)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c33-d1, alpha-achieved=1.24 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 19LL, 51LL, 4LL, 4LL, 128LL },
          { 1LL, 11LL, 11LL, 1LL, 1LL, 11LL },
          1, "0.0228535459746" },
        // class 3d-c33 key=[1/1 19/11 51/11 4/1 4/1 128/11]
        //   rep-vertices: (8,0,0) (0,8,0) (2,2,2) (3,5,1)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c33-d2, alpha-achieved=2.43 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 19LL, 51LL, 4LL, 4LL, 128LL },
          { 1LL, 11LL, 11LL, 1LL, 1LL, 11LL },
          2, "0.0178550960548" },
        // class 3d-c34 key=[1/1 83/48 11/3 91/48 16/3 59/48]
        //   rep-vertices: (16,0,0) (4,4,4) (20,4,4) (11,7,3)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c34-d0, alpha-achieved=0.573 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 83LL, 11LL, 91LL, 16LL, 59LL },
          { 1LL, 48LL, 3LL, 48LL, 3LL, 48LL },
          0, "0.0522538896685" },
        // class 3d-c34 key=[1/1 83/48 11/3 91/48 16/3 59/48]
        //   rep-vertices: (16,0,0) (4,4,4) (20,4,4) (11,7,3)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c34-d1, alpha-achieved=1.57 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 83LL, 11LL, 91LL, 16LL, 59LL },
          { 1LL, 48LL, 3LL, 48LL, 3LL, 48LL },
          1, "0.0311545926142" },
        // class 3d-c34 key=[1/1 83/48 11/3 91/48 16/3 59/48]
        //   rep-vertices: (16,0,0) (4,4,4) (20,4,4) (11,7,3)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c34-d2, alpha-achieved=3 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 83LL, 11LL, 91LL, 16LL, 59LL },
          { 1LL, 48LL, 3LL, 48LL, 3LL, 48LL },
          2, "0.0253995505726" },
        // class 3d-c35 key=[1/1 83/48 35/3 91/48 32/3 73/16]
        //   rep-vertices: (16,0,0) (11,7,3) (20,4,4) (4,20,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c35-d0, alpha-achieved=0.625 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 83LL, 35LL, 91LL, 32LL, 73LL },
          { 1LL, 48LL, 3LL, 48LL, 3LL, 16LL },
          0, "0.0619529046025" },
        // class 3d-c35 key=[1/1 83/48 35/3 91/48 32/3 73/16]
        //   rep-vertices: (16,0,0) (11,7,3) (20,4,4) (4,20,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c35-d1, alpha-achieved=1.68 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 83LL, 35LL, 91LL, 32LL, 73LL },
          { 1LL, 48LL, 3LL, 48LL, 3LL, 16LL },
          1, "0.0380092372712" },
        // class 3d-c35 key=[1/1 83/48 35/3 91/48 32/3 73/16]
        //   rep-vertices: (16,0,0) (11,7,3) (20,4,4) (4,20,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c35-d2, alpha-achieved=3.07 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 83LL, 35LL, 91LL, 32LL, 73LL },
          { 1LL, 48LL, 3LL, 48LL, 3LL, 16LL },
          2, "0.0315787937761" },
        // class 3d-c36 key=[1/1 9/5 9/5 16/5 16/5 32/5]
        //   rep-vertices: (4,4,0) (5,3,4) (8,4,4) (4,4,8)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c36-d0, alpha-achieved=0.73 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 9LL, 9LL, 16LL, 16LL, 32LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          0, "0.0535846815783" },
        // class 3d-c36 key=[1/1 9/5 9/5 16/5 16/5 32/5]
        //   rep-vertices: (4,4,0) (5,3,4) (8,4,4) (4,4,8)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c36-d1, alpha-achieved=1.85 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 9LL, 9LL, 16LL, 16LL, 32LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          1, "0.0348323163201" },
        // class 3d-c36 key=[1/1 9/5 9/5 16/5 16/5 32/5]
        //   rep-vertices: (4,4,0) (5,3,4) (8,4,4) (4,4,8)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c36-d2, alpha-achieved=3.48 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 9LL, 9LL, 16LL, 16LL, 32LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          2, "0.0291116337613" },
        // class 3d-c37 key=[1/1 44/19 128/19 51/19 51/19 172/19]
        //   rep-vertices: (8,8,0) (7,3,5) (10,2,2) (2,2,10)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c37-d0, alpha-achieved=0.621 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 44LL, 128LL, 51LL, 51LL, 172LL },
          { 1LL, 19LL, 19LL, 19LL, 19LL, 19LL },
          0, "0.0542008053349" },
        // class 3d-c37 key=[1/1 44/19 128/19 51/19 51/19 172/19]
        //   rep-vertices: (8,8,0) (7,3,5) (10,2,2) (2,2,10)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c37-d1, alpha-achieved=1.67 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 44LL, 128LL, 51LL, 51LL, 172LL },
          { 1LL, 19LL, 19LL, 19LL, 19LL, 19LL },
          1, "0.0331725677428" },
        // class 3d-c37 key=[1/1 44/19 128/19 51/19 51/19 172/19]
        //   rep-vertices: (8,8,0) (7,3,5) (10,2,2) (2,2,10)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c37-d2, alpha-achieved=3.14 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 44LL, 128LL, 51LL, 51LL, 172LL },
          { 1LL, 19LL, 19LL, 19LL, 19LL, 19LL },
          2, "0.0273829927676" },
        // class 3d-c38 key=[1/1 7/3 7/3 8/3 16/3 8/1]
        //   rep-vertices: (0,0,0) (3,2,1) (4,4,0) (4,4,4)  prov=seed:Kuhn0>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c38-d0, alpha-achieved=0.822 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 7LL, 7LL, 8LL, 16LL, 8LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 1LL },
          0, "0.0601233117732" },
        // class 3d-c38 key=[1/1 7/3 7/3 8/3 16/3 8/1]
        //   rep-vertices: (0,0,0) (3,2,1) (4,4,0) (4,4,4)  prov=seed:Kuhn0>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c38-d1, alpha-achieved=2.03 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 7LL, 7LL, 8LL, 16LL, 8LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 1LL },
          1, "0.0405057066259" },
        // class 3d-c38 key=[1/1 7/3 7/3 8/3 16/3 8/1]
        //   rep-vertices: (0,0,0) (3,2,1) (4,4,0) (4,4,4)  prov=seed:Kuhn0>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c38-d2, alpha-achieved=3.48 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 7LL, 7LL, 8LL, 16LL, 8LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 1LL },
          2, "0.0349182785291" },
        // class 3d-c39 key=[1/1 12/5 32/5 13/5 13/5 44/5]
        //   rep-vertices: (4,4,0) (4,2,3) (6,2,2) (2,2,6)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c39-d0, alpha-achieved=0.626 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 12LL, 32LL, 13LL, 13LL, 44LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          0, "0.0539276712087" },
        // class 3d-c39 key=[1/1 12/5 32/5 13/5 13/5 44/5]
        //   rep-vertices: (4,4,0) (4,2,3) (6,2,2) (2,2,6)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c39-d1, alpha-achieved=1.68 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 12LL, 32LL, 13LL, 13LL, 44LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          1, "0.0331018892393" },
        // class 3d-c39 key=[1/1 12/5 32/5 13/5 13/5 44/5]
        //   rep-vertices: (4,4,0) (4,2,3) (6,2,2) (2,2,6)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c39-d2, alpha-achieved=3.16 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 12LL, 32LL, 13LL, 13LL, 44LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          2, "0.0273297924390" },
        // class 3d-c40 key=[1/1 27/11 51/11 64/11 4/1 140/11]
        //   rep-vertices: (8,0,0) (3,5,1) (2,2,2) (2,10,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c40-d0, alpha-achieved=0.658 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 27LL, 51LL, 64LL, 4LL, 140LL },
          { 1LL, 11LL, 11LL, 11LL, 1LL, 11LL },
          0, "0.0494634852568" },
        // class 3d-c40 key=[1/1 27/11 51/11 64/11 4/1 140/11]
        //   rep-vertices: (8,0,0) (3,5,1) (2,2,2) (2,10,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c40-d1, alpha-achieved=1.76 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 27LL, 51LL, 64LL, 4LL, 140LL },
          { 1LL, 11LL, 11LL, 11LL, 1LL, 11LL },
          1, "0.0307703286026" },
        // class 3d-c40 key=[1/1 27/11 51/11 64/11 4/1 140/11]
        //   rep-vertices: (8,0,0) (3,5,1) (2,2,2) (2,10,2)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c40-d2, alpha-achieved=3.43 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 27LL, 51LL, 64LL, 4LL, 140LL },
          { 1LL, 11LL, 11LL, 11LL, 1LL, 11LL },
          2, "0.0253631241844" },
        // class 3d-c41 key=[1/1 49/17 49/17 256/51 256/51 512/51]
        //   rep-vertices: (0,0,0) (16,0,0) (0,16,0) (5,5,1)  prov=seed:T3>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c41-d0, alpha-achieved=0.724 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 49LL, 49LL, 256LL, 256LL, 512LL },
          { 1LL, 17LL, 17LL, 51LL, 51LL, 51LL },
          0, "0.0538399257073" },
        // class 3d-c41 key=[1/1 49/17 49/17 256/51 256/51 512/51]
        //   rep-vertices: (0,0,0) (16,0,0) (0,16,0) (5,5,1)  prov=seed:T3>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c41-d1, alpha-achieved=1.85 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 49LL, 49LL, 256LL, 256LL, 512LL },
          { 1LL, 17LL, 17LL, 51LL, 51LL, 51LL },
          1, "0.0348571063798" },
        // class 3d-c41 key=[1/1 49/17 49/17 256/51 256/51 512/51]
        //   rep-vertices: (0,0,0) (16,0,0) (0,16,0) (5,5,1)  prov=seed:T3>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c41-d2, alpha-achieved=3.48 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 49LL, 49LL, 256LL, 256LL, 512LL },
          { 1LL, 17LL, 17LL, 51LL, 51LL, 51LL },
          2, "0.0290982518895" },
        // class 3d-c42 key=[1/1 16/5 48/5 21/5 29/5 96/5]
        //   rep-vertices: (4,0,0) (4,4,0) (5,4,2) (8,8,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c42-d0, alpha-achieved=0.838 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 16LL, 48LL, 21LL, 29LL, 96LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          0, "0.0561897062205" },
        // class 3d-c42 key=[1/1 16/5 48/5 21/5 29/5 96/5]
        //   rep-vertices: (4,0,0) (4,4,0) (5,4,2) (8,8,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c42-d1, alpha-achieved=2.09 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 16LL, 48LL, 21LL, 29LL, 96LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          1, "0.0378635574728" },
        // class 3d-c42 key=[1/1 16/5 48/5 21/5 29/5 96/5]
        //   rep-vertices: (4,0,0) (4,4,0) (5,4,2) (8,8,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c42-d2, alpha-achieved=3.6 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 16LL, 48LL, 21LL, 29LL, 96LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          2, "0.0327392840555" },
        // class 3d-c43 key=[1/1 11/3 11/3 16/3 16/3 32/3]
        //   rep-vertices: (0,0,0) (4,0,0) (0,4,0) (1,1,1)  prov=seed:T3
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c43-d0, alpha-achieved=0.705 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 11LL, 11LL, 16LL, 16LL, 32LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 3LL },
          0, "0.0546669252993" },
        // class 3d-c43 key=[1/1 11/3 11/3 16/3 16/3 32/3]
        //   rep-vertices: (0,0,0) (4,0,0) (0,4,0) (1,1,1)  prov=seed:T3
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c43-d1, alpha-achieved=1.8 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 11LL, 11LL, 16LL, 16LL, 32LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 3LL },
          1, "0.0351401475060" },
        // class 3d-c43 key=[1/1 11/3 11/3 16/3 16/3 32/3]
        //   rep-vertices: (0,0,0) (4,0,0) (0,4,0) (1,1,1)  prov=seed:T3
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c43-d2, alpha-achieved=3.43 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 11LL, 11LL, 16LL, 16LL, 32LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 3LL },
          2, "0.0292004381892" },
        // class 3d-c44 key=[1/1 11/3 32/3 16/3 35/3 11/3]
        //   rep-vertices: (4,0,0) (0,4,0) (1,1,1) (1,5,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c44-d0, alpha-achieved=0.519 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 11LL, 32LL, 16LL, 35LL, 11LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 3LL },
          0, "0.0575045971561" },
        // class 3d-c44 key=[1/1 11/3 32/3 16/3 35/3 11/3]
        //   rep-vertices: (4,0,0) (0,4,0) (1,1,1) (1,5,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c44-d1, alpha-achieved=1.45 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 11LL, 32LL, 16LL, 35LL, 11LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 3LL },
          1, "0.0332008551213" },
        // class 3d-c44 key=[1/1 11/3 32/3 16/3 35/3 11/3]
        //   rep-vertices: (4,0,0) (0,4,0) (1,1,1) (1,5,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c44-d2, alpha-achieved=2.72 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 11LL, 32LL, 16LL, 35LL, 11LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 3LL },
          2, "0.0268656135074" },
        // class 3d-c45 key=[1/1 11/3 35/3 16/3 32/3 16/3]
        //   rep-vertices: (4,0,0) (1,1,1) (5,1,1) (1,5,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c45-d0, alpha-achieved=0.599 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 11LL, 35LL, 16LL, 32LL, 16LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 3LL },
          0, "0.0627840463054" },
        // class 3d-c45 key=[1/1 11/3 35/3 16/3 32/3 16/3]
        //   rep-vertices: (4,0,0) (1,1,1) (5,1,1) (1,5,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c45-d1, alpha-achieved=1.61 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 11LL, 35LL, 16LL, 32LL, 16LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 3LL },
          1, "0.0381279428252" },
        // class 3d-c45 key=[1/1 11/3 35/3 16/3 32/3 16/3]
        //   rep-vertices: (4,0,0) (1,1,1) (5,1,1) (1,5,1)  prov=bey-layer>bey1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c45-d2, alpha-achieved=2.94 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 11LL, 35LL, 16LL, 32LL, 16LL },
          { 1LL, 3LL, 3LL, 3LL, 3LL, 3LL },
          2, "0.0315184817333" },
        // class 3d-c46 key=[1/1 21/5 29/5 32/5 32/5 96/5]
        //   rep-vertices: (4,0,0) (5,4,2) (4,4,4) (8,8,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c46-d0, alpha-achieved=0.75 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 21LL, 29LL, 32LL, 32LL, 96LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          0, "0.0503447928355" },
        // class 3d-c46 key=[1/1 21/5 29/5 32/5 32/5 96/5]
        //   rep-vertices: (4,0,0) (5,4,2) (4,4,4) (8,8,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c46-d1, alpha-achieved=1.92 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 21LL, 29LL, 32LL, 32LL, 96LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          1, "0.0328125247457" },
        // class 3d-c46 key=[1/1 21/5 29/5 32/5 32/5 96/5]
        //   rep-vertices: (4,0,0) (5,4,2) (4,4,4) (8,8,4)  prov=bey-layer>bey1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c46-d2, alpha-achieved=3.56 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 21LL, 29LL, 32LL, 32LL, 96LL },
          { 1LL, 5LL, 5LL, 5LL, 5LL, 5LL },
          2, "0.0276316511245" },
        // class 3d-c47 key=[1/1 48/11 16/1 51/11 147/11 256/11]
        //   rep-vertices: (0,0,0) (5,5,1) (0,16,0) (4,4,4)  prov=seed:T3>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c47-d0, alpha-achieved=0.714 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 48LL, 16LL, 51LL, 147LL, 256LL },
          { 1LL, 11LL, 1LL, 11LL, 11LL, 11LL },
          0, "0.0640124620544" },
        // class 3d-c47 key=[1/1 48/11 16/1 51/11 147/11 256/11]
        //   rep-vertices: (0,0,0) (5,5,1) (0,16,0) (4,4,4)  prov=seed:T3>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c47-d1, alpha-achieved=1.87 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 48LL, 16LL, 51LL, 147LL, 256LL },
          { 1LL, 11LL, 1LL, 11LL, 11LL, 11LL },
          1, "0.0408995472784" },
        // class 3d-c47 key=[1/1 48/11 16/1 51/11 147/11 256/11]
        //   rep-vertices: (0,0,0) (5,5,1) (0,16,0) (4,4,4)  prov=seed:T3>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c47-d2, alpha-achieved=3.27 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 48LL, 16LL, 51LL, 147LL, 256LL },
          { 1LL, 11LL, 1LL, 11LL, 11LL, 11LL },
          2, "0.0348312544997" },
        // class 3d-c48 key=[1/1 147/11 147/11 16/1 16/1 512/11]
        //   rep-vertices: (5,5,1) (16,0,0) (0,16,0) (4,4,4)  prov=seed:T3>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c48-d0, alpha-achieved=0.721 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 147LL, 147LL, 16LL, 16LL, 512LL },
          { 1LL, 11LL, 11LL, 1LL, 1LL, 11LL },
          0, "0.0506159148982" },
        // class 3d-c48 key=[1/1 147/11 147/11 16/1 16/1 512/11]
        //   rep-vertices: (5,5,1) (16,0,0) (0,16,0) (4,4,4)  prov=seed:T3>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c48-d1, alpha-achieved=1.87 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 147LL, 147LL, 16LL, 16LL, 512LL },
          { 1LL, 11LL, 11LL, 1LL, 1LL, 11LL },
          1, "0.0325210950833" },
        // class 3d-c48 key=[1/1 147/11 147/11 16/1 16/1 512/11]
        //   rep-vertices: (5,5,1) (16,0,0) (0,16,0) (4,4,4)  prov=seed:T3>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c48-d2, alpha-achieved=3.55 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 147LL, 147LL, 16LL, 16LL, 512LL },
          { 1LL, 11LL, 11LL, 1LL, 1LL, 11LL },
          2, "0.0271682189361" },
        // class 3d-c49 key=[1/1 43/3 43/3 16/1 16/1 128/3]
        //   rep-vertices: (3,3,3) (0,8,0) (0,0,8) (4,4,4)  prov=seed:T5>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c49-d0, alpha-achieved=0.701 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 43LL, 43LL, 16LL, 16LL, 128LL },
          { 1LL, 3LL, 3LL, 1LL, 1LL, 3LL },
          0, "0.0521297470838" },
        // class 3d-c49 key=[1/1 43/3 43/3 16/1 16/1 128/3]
        //   rep-vertices: (3,3,3) (0,8,0) (0,0,8) (4,4,4)  prov=seed:T5>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c49-d1, alpha-achieved=1.82 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 43LL, 43LL, 16LL, 16LL, 128LL },
          { 1LL, 3LL, 3LL, 1LL, 1LL, 3LL },
          1, "0.0333033004821" },
        // class 3d-c49 key=[1/1 43/3 43/3 16/1 16/1 128/3]
        //   rep-vertices: (3,3,3) (0,8,0) (0,0,8) (4,4,4)  prov=seed:T5>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c49-d2, alpha-achieved=3.51 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 43LL, 43LL, 16LL, 16LL, 128LL },
          { 1LL, 3LL, 3LL, 1LL, 1LL, 3LL },
          2, "0.0275924511958" },
        // class 3d-c50 key=[1/1 57/1 57/1 176/3 176/3 512/3]
        //   rep-vertices: (5,5,5) (0,16,0) (0,0,16) (4,4,4)  prov=seed:T4>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c50-d0, alpha-achieved=0.701 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 57LL, 57LL, 176LL, 176LL, 512LL },
          { 1LL, 1LL, 1LL, 3LL, 3LL, 3LL },
          0, "0.0514194424654" },
        // class 3d-c50 key=[1/1 57/1 57/1 176/3 176/3 512/3]
        //   rep-vertices: (5,5,5) (0,16,0) (0,0,16) (4,4,4)  prov=seed:T4>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c50-d1, alpha-achieved=1.82 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 57LL, 57LL, 176LL, 176LL, 512LL },
          { 1LL, 1LL, 1LL, 3LL, 3LL, 3LL },
          1, "0.0328160026169" },
        // class 3d-c50 key=[1/1 57/1 57/1 176/3 176/3 512/3]
        //   rep-vertices: (5,5,5) (0,16,0) (0,0,16) (4,4,4)  prov=seed:T4>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=2,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-3d-c50-d2, alpha-achieved=3.51 (alpha not met, cap-curtailed at L=2; B-2b' L=3 sharpening planned)
        { { 1LL, 57LL, 57LL, 176LL, 176LL, 512LL },
          { 1LL, 1LL, 1LL, 3LL, 3LL, 3LL },
          2, "0.0272302982412" },
    };
    count = static_cast<int>(sizeof entries / sizeof entries[0]);
    return entries;
}

} // namespace detail

// VCP_CONSTANTS_TABLE_END

inline const l2_projection_registry_entry_3d*
l2_projection_registry_3d_table(int& count) {
    return detail::l2_projection_registry_3d_entries(count);
}

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_DICT_REGISTRY_3D_HPP
