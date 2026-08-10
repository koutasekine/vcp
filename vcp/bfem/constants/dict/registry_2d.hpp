// vcp/bfem/constants/dict/registry_2d.hpp
//
// CONST-B2a: generated 2D class registry (dictionary first layer).  This file
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

#ifndef VCP_BFEM_CONSTANTS_DICT_REGISTRY_2D_HPP
#define VCP_BFEM_CONSTANTS_DICT_REGISTRY_2D_HPP

#include <vcp/bfem/constants/dict/dict_entry.hpp>

namespace vcp {
namespace bfem {
namespace constants {

// missing (class, d) list -- ruled local scope (addendum sections 1, 8):
//   every class below is generated for d = 0..5; d = 6..9 is NOT generated locally
//   (alpha unmet within the job cap at L = 3; B-2b' remote sharpening scope):
//   class 2d-c00 missing d=6..9
//   class 2d-c01 missing d=6..9
//   class 2d-c02 missing d=6..9
//   class 2d-c03 missing d=6..9
//   class 2d-c04 missing d=6..9
//   class 2d-c05 missing d=6..9
//   class 2d-c06 missing d=6..9
//   class 2d-c07 missing d=6..9
//   class 2d-c08 missing d=6..9
//   class 2d-c09 missing d=6..9
//   class 2d-c10 missing d=6..9
//   class 2d-c11 missing d=6..9
//   class 2d-c12 missing d=6..9
//   class 2d-c13 missing d=6..9
//   class 2d-c14 missing d=6..9
//   class 2d-c15 missing d=6..9
//   class 2d-c16 missing d=6..9
//   class 2d-c17 missing d=6..9
//   class 2d-c18 missing d=6..9
//   class 2d-c19 missing d=6..9
//   class 2d-c20 missing d=6..9
//   class 2d-c21 missing d=6..9
//   class 2d-c22 missing d=6..9
//   class 2d-c23 missing d=6..9
//   class 2d-c24 missing d=6..9
//   class 2d-c25 missing d=6..9
//   class 2d-c26 missing d=6..9
//   class 2d-c27 missing d=6..9
//   class 2d-c28 missing d=6..9
//   class 2d-c29 missing d=6..9
//   class 2d-c30 missing d=6..9
//   class 2d-c31 missing d=6..9
//   class 2d-c32 missing d=6..9

// VCP_CONSTANTS_TABLE_BEGIN  (authorized decimal strings; generated table)

namespace detail {

inline const l2_projection_registry_entry_2d* l2_projection_registry_2d_entries(int& count) {
    static const l2_projection_registry_entry_2d entries[] = {
        // class 2d-c00 key=[1/1 1/1 2/1]
        //   rep-vertices: (0,0) (1,0) (0,1)  prov=seed:K1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c00-d0, alpha-achieved=0.011 (alpha ok, beta ok)
        { { 1LL, 1LL, 2LL },
          { 1LL, 1LL, 1LL },
          0, "0.0516586587537" },
        // class 2d-c00 key=[1/1 1/1 2/1]
        //   rep-vertices: (0,0) (1,0) (0,1)  prov=seed:K1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c00-d1, alpha-achieved=0.0355 (alpha ok, beta ok)
        { { 1LL, 1LL, 2LL },
          { 1LL, 1LL, 1LL },
          1, "0.0163360924264" },
        // class 2d-c00 key=[1/1 1/1 2/1]
        //   rep-vertices: (0,0) (1,0) (0,1)  prov=seed:K1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c00-d2, alpha-achieved=0.0714 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 2LL },
          { 1LL, 1LL, 1LL },
          2, "0.00840003917742" },
        // class 2d-c00 key=[1/1 1/1 2/1]
        //   rep-vertices: (0,0) (1,0) (0,1)  prov=seed:K1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c00-d3, alpha-achieved=0.113 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 2LL },
          { 1LL, 1LL, 1LL },
          3, "0.00552866969022" },
        // class 2d-c00 key=[1/1 1/1 2/1]
        //   rep-vertices: (0,0) (1,0) (0,1)  prov=seed:K1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c00-d4, alpha-achieved=0.158 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 2LL },
          { 1LL, 1LL, 1LL },
          4, "0.00410576026661" },
        // class 2d-c00 key=[1/1 1/1 2/1]
        //   rep-vertices: (0,0) (1,0) (0,1)  prov=seed:K1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c00-d5, alpha-achieved=0.206 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 2LL },
          { 1LL, 1LL, 1LL },
          5, "0.00328136596448" },
        // class 2d-c01 key=[1/1 1/1 36/13]
        //   rep-vertices: (12,0) (12,6) (10,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c01-d0, alpha-achieved=0.0116 (alpha ok, beta ok)
        { { 1LL, 1LL, 36LL },
          { 1LL, 1LL, 13LL },
          0, "0.0486453134125" },
        // class 2d-c01 key=[1/1 1/1 36/13]
        //   rep-vertices: (12,0) (12,6) (10,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c01-d1, alpha-achieved=0.0378 (alpha ok, beta ok)
        { { 1LL, 1LL, 36LL },
          { 1LL, 1LL, 13LL },
          1, "0.0153599173628" },
        // class 2d-c01 key=[1/1 1/1 36/13]
        //   rep-vertices: (12,0) (12,6) (10,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c01-d2, alpha-achieved=0.0786 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 36LL },
          { 1LL, 1LL, 13LL },
          2, "0.00768613587611" },
        // class 2d-c01 key=[1/1 1/1 36/13]
        //   rep-vertices: (12,0) (12,6) (10,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c01-d3, alpha-achieved=0.124 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 36LL },
          { 1LL, 1LL, 13LL },
          3, "0.00506836272850" },
        // class 2d-c01 key=[1/1 1/1 36/13]
        //   rep-vertices: (12,0) (12,6) (10,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c01-d4, alpha-achieved=0.172 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 36LL },
          { 1LL, 1LL, 13LL },
          4, "0.00380616691999" },
        // class 2d-c01 key=[1/1 1/1 36/13]
        //   rep-vertices: (12,0) (12,6) (10,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c01-d5, alpha-achieved=0.221 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 36LL },
          { 1LL, 1LL, 13LL },
          5, "0.00308981525765" },
        // class 2d-c02 key=[1/1 1/1 16/5]
        //   rep-vertices: (4,0) (2,1) (0,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c02-d0, alpha-achieved=0.012 (alpha ok, beta ok)
        { { 1LL, 1LL, 16LL },
          { 1LL, 1LL, 5LL },
          0, "0.0472868321768" },
        // class 2d-c02 key=[1/1 1/1 16/5]
        //   rep-vertices: (4,0) (2,1) (0,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c02-d1, alpha-achieved=0.039 (alpha ok, beta ok)
        { { 1LL, 1LL, 16LL },
          { 1LL, 1LL, 5LL },
          1, "0.0149233805926" },
        // class 2d-c02 key=[1/1 1/1 16/5]
        //   rep-vertices: (4,0) (2,1) (0,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c02-d2, alpha-achieved=0.0824 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 16LL },
          { 1LL, 1LL, 5LL },
          2, "0.00735669445392" },
        // class 2d-c02 key=[1/1 1/1 16/5]
        //   rep-vertices: (4,0) (2,1) (0,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c02-d3, alpha-achieved=0.13 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 16LL },
          { 1LL, 1LL, 5LL },
          3, "0.00485166098959" },
        // class 2d-c02 key=[1/1 1/1 16/5]
        //   rep-vertices: (4,0) (2,1) (0,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c02-d4, alpha-achieved=0.18 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 16LL },
          { 1LL, 1LL, 5LL },
          4, "0.00366850463374" },
        // class 2d-c02 key=[1/1 1/1 16/5]
        //   rep-vertices: (4,0) (2,1) (0,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c02-d5, alpha-achieved=0.229 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 16LL },
          { 1LL, 1LL, 5LL },
          5, "0.00300705888286" },
        // class 2d-c03 key=[1/1 1/1 18/5]
        //   rep-vertices: (1,1) (3,0) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c03-d0, alpha-achieved=0.0123 (alpha ok, beta ok)
        { { 1LL, 1LL, 18LL },
          { 1LL, 1LL, 5LL },
          0, "0.0461719782692" },
        // class 2d-c03 key=[1/1 1/1 18/5]
        //   rep-vertices: (1,1) (3,0) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c03-d1, alpha-achieved=0.04 (alpha ok, beta ok)
        { { 1LL, 1LL, 18LL },
          { 1LL, 1LL, 5LL },
          1, "0.0145616468566" },
        // class 2d-c03 key=[1/1 1/1 18/5]
        //   rep-vertices: (1,1) (3,0) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c03-d2, alpha-achieved=0.0858 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 18LL },
          { 1LL, 1LL, 5LL },
          2, "0.00708398159035" },
        // class 2d-c03 key=[1/1 1/1 18/5]
        //   rep-vertices: (1,1) (3,0) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c03-d3, alpha-achieved=0.136 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 18LL },
          { 1LL, 1LL, 5LL },
          3, "0.00467233591485" },
        // class 2d-c03 key=[1/1 1/1 18/5]
        //   rep-vertices: (1,1) (3,0) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c03-d4, alpha-achieved=0.187 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 18LL },
          { 1LL, 1LL, 5LL },
          4, "0.00355649392688" },
        // class 2d-c03 key=[1/1 1/1 18/5]
        //   rep-vertices: (1,1) (3,0) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c03-d5, alpha-achieved=0.235 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 18LL },
          { 1LL, 1LL, 5LL },
          5, "0.00294179696077" },
        // class 2d-c04 key=[1/1 1/1 144/37]
        //   rep-vertices: (12,0) (6,1) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c04-d0, alpha-achieved=0.0125 (alpha ok, beta ok)
        { { 1LL, 1LL, 144LL },
          { 1LL, 1LL, 37LL },
          0, "0.0454303644012" },
        // class 2d-c04 key=[1/1 1/1 144/37]
        //   rep-vertices: (12,0) (6,1) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c04-d1, alpha-achieved=0.0407 (alpha ok, beta ok)
        { { 1LL, 1LL, 144LL },
          { 1LL, 1LL, 37LL },
          1, "0.0143183624634" },
        // class 2d-c04 key=[1/1 1/1 144/37]
        //   rep-vertices: (12,0) (6,1) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c04-d2, alpha-achieved=0.0883 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 144LL },
          { 1LL, 1LL, 37LL },
          2, "0.00690272447246" },
        // class 2d-c04 key=[1/1 1/1 144/37]
        //   rep-vertices: (12,0) (6,1) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c04-d3, alpha-achieved=0.14 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 144LL },
          { 1LL, 1LL, 37LL },
          3, "0.00455440788917" },
        // class 2d-c04 key=[1/1 1/1 144/37]
        //   rep-vertices: (12,0) (6,1) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c04-d4, alpha-achieved=0.191 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 144LL },
          { 1LL, 1LL, 37LL },
          4, "0.00348386411423" },
        // class 2d-c04 key=[1/1 1/1 144/37]
        //   rep-vertices: (12,0) (6,1) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c04-d5, alpha-achieved=0.239 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 1LL, 144LL },
          { 1LL, 1LL, 37LL },
          5, "0.00290019553303" },
        // class 2d-c05 key=[1/1 32/29 117/29]
        //   rep-vertices: (18,9) (20,4) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c05-d0, alpha-achieved=0.0124 (alpha ok, beta ok)
        { { 1LL, 32LL, 117LL },
          { 1LL, 29LL, 29LL },
          0, "0.0455873240106" },
        // class 2d-c05 key=[1/1 32/29 117/29]
        //   rep-vertices: (18,9) (20,4) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c05-d1, alpha-achieved=0.0405 (alpha ok, beta ok)
        { { 1LL, 32LL, 117LL },
          { 1LL, 29LL, 29LL },
          1, "0.0143695942688" },
        // class 2d-c05 key=[1/1 32/29 117/29]
        //   rep-vertices: (18,9) (20,4) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c05-d2, alpha-achieved=0.0877 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 32LL, 117LL },
          { 1LL, 29LL, 29LL },
          2, "0.00694447428806" },
        // class 2d-c05 key=[1/1 32/29 117/29]
        //   rep-vertices: (18,9) (20,4) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c05-d3, alpha-achieved=0.139 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 32LL, 117LL },
          { 1LL, 29LL, 29LL },
          3, "0.00458048164807" },
        // class 2d-c05 key=[1/1 32/29 117/29]
        //   rep-vertices: (18,9) (20,4) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c05-d4, alpha-achieved=0.19 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 32LL, 117LL },
          { 1LL, 29LL, 29LL },
          4, "0.00349979352271" },
        // class 2d-c05 key=[1/1 32/29 117/29]
        //   rep-vertices: (18,9) (20,4) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c05-d5, alpha-achieved=0.238 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 32LL, 117LL },
          { 1LL, 29LL, 29LL },
          5, "0.00290946355156" },
        // class 2d-c06 key=[1/1 16/13 45/13]
        //   rep-vertices: (10,3) (12,6) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c06-d0, alpha-achieved=0.0119 (alpha ok, beta ok)
        { { 1LL, 16LL, 45LL },
          { 1LL, 13LL, 13LL },
          0, "0.0477055540517" },
        // class 2d-c06 key=[1/1 16/13 45/13]
        //   rep-vertices: (10,3) (12,6) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c06-d1, alpha-achieved=0.0386 (alpha ok, beta ok)
        { { 1LL, 16LL, 45LL },
          { 1LL, 13LL, 13LL },
          1, "0.0150571184265" },
        // class 2d-c06 key=[1/1 16/13 45/13]
        //   rep-vertices: (10,3) (12,6) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c06-d2, alpha-achieved=0.081 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 16LL, 45LL },
          { 1LL, 13LL, 13LL },
          2, "0.00746995543720" },
        // class 2d-c06 key=[1/1 16/13 45/13]
        //   rep-vertices: (10,3) (12,6) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c06-d3, alpha-achieved=0.128 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 16LL, 45LL },
          { 1LL, 13LL, 13LL },
          3, "0.00492226371718" },
        // class 2d-c06 key=[1/1 16/13 45/13]
        //   rep-vertices: (10,3) (12,6) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c06-d4, alpha-achieved=0.178 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 16LL, 45LL },
          { 1LL, 13LL, 13LL },
          4, "0.00371329181288" },
        // class 2d-c06 key=[1/1 16/13 45/13]
        //   rep-vertices: (10,3) (12,6) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c06-d5, alpha-achieved=0.226 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 16LL, 45LL },
          { 1LL, 13LL, 13LL },
          5, "0.00303534553723" },
        // class 2d-c07 key=[1/1 5/4 5/4]
        //   rep-vertices: (4,0) (4,2) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c07-d0, alpha-achieved=0.00942 (alpha ok, beta ok)
        { { 1LL, 5LL, 5LL },
          { 1LL, 4LL, 4LL },
          0, "0.0600019854060" },
        // class 2d-c07 key=[1/1 5/4 5/4]
        //   rep-vertices: (4,0) (4,2) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c07-d1, alpha-achieved=0.0298 (alpha ok, beta ok)
        { { 1LL, 5LL, 5LL },
          { 1LL, 4LL, 4LL },
          1, "0.0193773153121" },
        // class 2d-c07 key=[1/1 5/4 5/4]
        //   rep-vertices: (4,0) (4,2) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c07-d2, alpha-achieved=0.0584 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 5LL },
          { 1LL, 4LL, 4LL },
          2, "0.0101431038379" },
        // class 2d-c07 key=[1/1 5/4 5/4]
        //   rep-vertices: (4,0) (4,2) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c07-d3, alpha-achieved=0.0923 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 5LL },
          { 1LL, 4LL, 4LL },
          3, "0.00662568109817" },
        // class 2d-c07 key=[1/1 5/4 5/4]
        //   rep-vertices: (4,0) (4,2) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c07-d4, alpha-achieved=0.131 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 5LL },
          { 1LL, 4LL, 4LL },
          4, "0.00484831760645" },
        // class 2d-c07 key=[1/1 5/4 5/4]
        //   rep-vertices: (4,0) (4,2) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c07-d5, alpha-achieved=0.172 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 5LL },
          { 1LL, 4LL, 4LL },
          5, "0.00381837039697" },
        // class 2d-c08 key=[1/1 5/4 13/4]
        //   rep-vertices: (6,3) (6,1) (8,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c08-d0, alpha-achieved=0.0117 (alpha ok, beta ok)
        { { 1LL, 5LL, 13LL },
          { 1LL, 4LL, 4LL },
          0, "0.0484224515179" },
        // class 2d-c08 key=[1/1 5/4 13/4]
        //   rep-vertices: (6,3) (6,1) (8,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c08-d1, alpha-achieved=0.038 (alpha ok, beta ok)
        { { 1LL, 5LL, 13LL },
          { 1LL, 4LL, 4LL },
          1, "0.0152874129757" },
        // class 2d-c08 key=[1/1 5/4 13/4]
        //   rep-vertices: (6,3) (6,1) (8,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c08-d2, alpha-achieved=0.079 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 13LL },
          { 1LL, 4LL, 4LL },
          2, "0.00764398505029" },
        // class 2d-c08 key=[1/1 5/4 13/4]
        //   rep-vertices: (6,3) (6,1) (8,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c08-d3, alpha-achieved=0.125 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 13LL },
          { 1LL, 4LL, 4LL },
          3, "0.00503645838703" },
        // class 2d-c08 key=[1/1 5/4 13/4]
        //   rep-vertices: (6,3) (6,1) (8,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c08-d4, alpha-achieved=0.174 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 13LL },
          { 1LL, 4LL, 4LL },
          4, "0.00378614100436" },
        // class 2d-c08 key=[1/1 5/4 13/4]
        //   rep-vertices: (6,3) (6,1) (8,0)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c08-d5, alpha-achieved=0.222 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 13LL },
          { 1LL, 4LL, 4LL },
          5, "0.00307981744440" },
        // class 2d-c09 key=[1/1 13/10 9/2]
        //   rep-vertices: (0,0) (6,3) (3,2)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c09-d0, alpha-achieved=0.0125 (alpha ok, beta ok)
        { { 1LL, 13LL, 9LL },
          { 1LL, 10LL, 2LL },
          0, "0.0455087912810" },
        // class 2d-c09 key=[1/1 13/10 9/2]
        //   rep-vertices: (0,0) (6,3) (3,2)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c09-d1, alpha-achieved=0.0406 (alpha ok, beta ok)
        { { 1LL, 13LL, 9LL },
          { 1LL, 10LL, 2LL },
          1, "0.0143409002394" },
        // class 2d-c09 key=[1/1 13/10 9/2]
        //   rep-vertices: (0,0) (6,3) (3,2)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c09-d2, alpha-achieved=0.0877 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 9LL },
          { 1LL, 10LL, 2LL },
          2, "0.00694498568648" },
        // class 2d-c09 key=[1/1 13/10 9/2]
        //   rep-vertices: (0,0) (6,3) (3,2)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c09-d3, alpha-achieved=0.139 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 9LL },
          { 1LL, 10LL, 2LL },
          3, "0.00457569093429" },
        // class 2d-c09 key=[1/1 13/10 9/2]
        //   rep-vertices: (0,0) (6,3) (3,2)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c09-d4, alpha-achieved=0.191 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 9LL },
          { 1LL, 10LL, 2LL },
          4, "0.00349666826153" },
        // class 2d-c09 key=[1/1 13/10 9/2]
        //   rep-vertices: (0,0) (6,3) (3,2)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c09-d5, alpha-achieved=0.238 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 9LL },
          { 1LL, 10LL, 2LL },
          5, "0.00290873205003" },
        // class 2d-c10 key=[1/1 29/20 9/4]
        //   rep-vertices: (22,7) (24,12) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c10-d0, alpha-achieved=0.0105 (alpha ok, beta ok)
        { { 1LL, 29LL, 9LL },
          { 1LL, 20LL, 4LL },
          0, "0.0537232614461" },
        // class 2d-c10 key=[1/1 29/20 9/4]
        //   rep-vertices: (22,7) (24,12) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c10-d1, alpha-achieved=0.0341 (alpha ok, beta ok)
        { { 1LL, 29LL, 9LL },
          { 1LL, 20LL, 4LL },
          1, "0.0169989169730" },
        // class 2d-c10 key=[1/1 29/20 9/4]
        //   rep-vertices: (22,7) (24,12) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c10-d2, alpha-achieved=0.0675 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 29LL, 9LL },
          { 1LL, 20LL, 4LL },
          2, "0.00885883705898" },
        // class 2d-c10 key=[1/1 29/20 9/4]
        //   rep-vertices: (22,7) (24,12) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c10-d3, alpha-achieved=0.106 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 29LL, 9LL },
          { 1LL, 20LL, 4LL },
          3, "0.00581846658626" },
        // class 2d-c10 key=[1/1 29/20 9/4]
        //   rep-vertices: (22,7) (24,12) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c10-d4, alpha-achieved=0.15 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 29LL, 9LL },
          { 1LL, 20LL, 4LL },
          4, "0.00430310392798" },
        // class 2d-c10 key=[1/1 29/20 9/4]
        //   rep-vertices: (22,7) (24,12) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c10-d5, alpha-achieved=0.195 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 29LL, 9LL },
          { 1LL, 20LL, 4LL },
          5, "0.00342558248547" },
        // class 2d-c11 key=[1/1 8/5 9/5]
        //   rep-vertices: (4,1) (6,0) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c11-d0, alpha-achieved=0.00956 (alpha ok, beta ok)
        { { 1LL, 8LL, 9LL },
          { 1LL, 5LL, 5LL },
          0, "0.0591332154490" },
        // class 2d-c11 key=[1/1 8/5 9/5]
        //   rep-vertices: (4,1) (6,0) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c11-d1, alpha-achieved=0.0308 (alpha ok, beta ok)
        { { 1LL, 8LL, 9LL },
          { 1LL, 5LL, 5LL },
          1, "0.0187633116270" },
        // class 2d-c11 key=[1/1 8/5 9/5]
        //   rep-vertices: (4,1) (6,0) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c11-d2, alpha-achieved=0.0599 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 8LL, 9LL },
          { 1LL, 5LL, 5LL },
          2, "0.00990387295188" },
        // class 2d-c11 key=[1/1 8/5 9/5]
        //   rep-vertices: (4,1) (6,0) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c11-d3, alpha-achieved=0.0946 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 8LL, 9LL },
          { 1LL, 5LL, 5LL },
          3, "0.00647949350041" },
        // class 2d-c11 key=[1/1 8/5 9/5]
        //   rep-vertices: (4,1) (6,0) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c11-d4, alpha-achieved=0.133 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 8LL, 9LL },
          { 1LL, 5LL, 5LL },
          4, "0.00475963350388" },
        // class 2d-c11 key=[1/1 8/5 9/5]
        //   rep-vertices: (4,1) (6,0) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c11-d5, alpha-achieved=0.175 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 8LL, 9LL },
          { 1LL, 5LL, 5LL },
          5, "0.00376151822297" },
        // class 2d-c12 key=[1/1 61/37 180/37]
        //   rep-vertices: (24,0) (18,5) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c12-d0, alpha-achieved=0.0122 (alpha ok, beta ok)
        { { 1LL, 61LL, 180LL },
          { 1LL, 37LL, 37LL },
          0, "0.0464918505974" },
        // class 2d-c12 key=[1/1 61/37 180/37]
        //   rep-vertices: (24,0) (18,5) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c12-d1, alpha-achieved=0.0397 (alpha ok, beta ok)
        { { 1LL, 61LL, 180LL },
          { 1LL, 37LL, 37LL },
          1, "0.0146559513942" },
        // class 2d-c12 key=[1/1 61/37 180/37]
        //   rep-vertices: (24,0) (18,5) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c12-d2, alpha-achieved=0.084 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 61LL, 180LL },
          { 1LL, 37LL, 37LL },
          2, "0.00722367370806" },
        // class 2d-c12 key=[1/1 61/37 180/37]
        //   rep-vertices: (24,0) (18,5) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c12-d3, alpha-achieved=0.134 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 61LL, 180LL },
          { 1LL, 37LL, 37LL },
          3, "0.00474487358400" },
        // class 2d-c12 key=[1/1 61/37 180/37]
        //   rep-vertices: (24,0) (18,5) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c12-d4, alpha-achieved=0.184 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 61LL, 180LL },
          { 1LL, 37LL, 37LL },
          4, "0.00360241059215" },
        // class 2d-c12 key=[1/1 61/37 180/37]
        //   rep-vertices: (24,0) (18,5) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c12-d5, alpha-achieved=0.232 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 61LL, 180LL },
          { 1LL, 37LL, 37LL },
          5, "0.00297518283283" },
        // class 2d-c13 key=[1/1 53/29 144/29]
        //   rep-vertices: (24,0) (24,12) (22,7)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c13-d0, alpha-achieved=0.012 (alpha ok, beta ok)
        { { 1LL, 53LL, 144LL },
          { 1LL, 29LL, 29LL },
          0, "0.0471969905966" },
        // class 2d-c13 key=[1/1 53/29 144/29]
        //   rep-vertices: (24,0) (24,12) (22,7)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c13-d1, alpha-achieved=0.0391 (alpha ok, beta ok)
        { { 1LL, 53LL, 144LL },
          { 1LL, 29LL, 29LL },
          1, "0.0148819544712" },
        // class 2d-c13 key=[1/1 53/29 144/29]
        //   rep-vertices: (24,0) (24,12) (22,7)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c13-d2, alpha-achieved=0.0818 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 53LL, 144LL },
          { 1LL, 29LL, 29LL },
          2, "0.00740769550076" },
        // class 2d-c13 key=[1/1 53/29 144/29]
        //   rep-vertices: (24,0) (24,12) (22,7)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c13-d3, alpha-achieved=0.13 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 53LL, 144LL },
          { 1LL, 29LL, 29LL },
          3, "0.00486035947378" },
        // class 2d-c13 key=[1/1 53/29 144/29]
        //   rep-vertices: (24,0) (24,12) (22,7)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c13-d4, alpha-achieved=0.18 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 53LL, 144LL },
          { 1LL, 29LL, 29LL },
          4, "0.00367615926391" },
        // class 2d-c13 key=[1/1 53/29 144/29]
        //   rep-vertices: (24,0) (24,12) (22,7)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c13-d5, alpha-achieved=0.227 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 53LL, 144LL },
          { 1LL, 29LL, 29LL },
          5, "0.00302272224861" },
        // class 2d-c14 key=[1/1 25/13 72/13]
        //   rep-vertices: (6,6) (0,0) (4,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c14-d0, alpha-achieved=0.0122 (alpha ok, beta ok)
        { { 1LL, 25LL, 72LL },
          { 1LL, 13LL, 13LL },
          0, "0.0464762464763" },
        // class 2d-c14 key=[1/1 25/13 72/13]
        //   rep-vertices: (6,6) (0,0) (4,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c14-d1, alpha-achieved=0.0398 (alpha ok, beta ok)
        { { 1LL, 25LL, 72LL },
          { 1LL, 13LL, 13LL },
          1, "0.0146438074195" },
        // class 2d-c14 key=[1/1 25/13 72/13]
        //   rep-vertices: (6,6) (0,0) (4,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c14-d2, alpha-achieved=0.0837 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 25LL, 72LL },
          { 1LL, 13LL, 13LL },
          2, "0.00724939419360" },
        // class 2d-c14 key=[1/1 25/13 72/13]
        //   rep-vertices: (6,6) (0,0) (4,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c14-d3, alpha-achieved=0.134 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 25LL, 72LL },
          { 1LL, 13LL, 13LL },
          3, "0.00475316668433" },
        // class 2d-c14 key=[1/1 25/13 72/13]
        //   rep-vertices: (6,6) (0,0) (4,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c14-d4, alpha-achieved=0.184 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 25LL, 72LL },
          { 1LL, 13LL, 13LL },
          4, "0.00360863653479" },
        // class 2d-c14 key=[1/1 25/13 72/13]
        //   rep-vertices: (6,6) (0,0) (4,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c14-d5, alpha-achieved=0.231 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 25LL, 72LL },
          { 1LL, 13LL, 13LL },
          5, "0.00298293670427" },
        // class 2d-c15 key=[1/1 2/1 5/1]
        //   rep-vertices: (0,0) (2,1) (1,1)  prov=seed:earclip:E2
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c15-d0, alpha-achieved=0.0118 (alpha ok, beta ok)
        { { 1LL, 2LL, 5LL },
          { 1LL, 1LL, 1LL },
          0, "0.0480257591155" },
        // class 2d-c15 key=[1/1 2/1 5/1]
        //   rep-vertices: (0,0) (2,1) (1,1)  prov=seed:earclip:E2
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c15-d1, alpha-achieved=0.0384 (alpha ok, beta ok)
        { { 1LL, 2LL, 5LL },
          { 1LL, 1LL, 1LL },
          1, "0.0151475275038" },
        // class 2d-c15 key=[1/1 2/1 5/1]
        //   rep-vertices: (0,0) (2,1) (1,1)  prov=seed:earclip:E2
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c15-d2, alpha-achieved=0.0794 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 2LL, 5LL },
          { 1LL, 1LL, 1LL },
          2, "0.00761288445269" },
        // class 2d-c15 key=[1/1 2/1 5/1]
        //   rep-vertices: (0,0) (2,1) (1,1)  prov=seed:earclip:E2
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c15-d3, alpha-achieved=0.126 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 2LL, 5LL },
          { 1LL, 1LL, 1LL },
          3, "0.00499175821642" },
        // class 2d-c15 key=[1/1 2/1 5/1]
        //   rep-vertices: (0,0) (2,1) (1,1)  prov=seed:earclip:E2
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c15-d4, alpha-achieved=0.175 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 2LL, 5LL },
          { 1LL, 1LL, 1LL },
          4, "0.00376119096324" },
        // class 2d-c15 key=[1/1 2/1 5/1]
        //   rep-vertices: (0,0) (2,1) (1,1)  prov=seed:earclip:E2
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c15-d5, alpha-achieved=0.222 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 2LL, 5LL },
          { 1LL, 1LL, 1LL },
          5, "0.00307825901426" },
        // class 2d-c16 key=[1/1 17/8 45/8]
        //   rep-vertices: (0,0) (4,1) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c16-d0, alpha-achieved=0.012 (alpha ok, beta ok)
        { { 1LL, 17LL, 45LL },
          { 1LL, 8LL, 8LL },
          0, "0.0472968655413" },
        // class 2d-c16 key=[1/1 17/8 45/8]
        //   rep-vertices: (0,0) (4,1) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c16-d1, alpha-achieved=0.039 (alpha ok, beta ok)
        { { 1LL, 17LL, 45LL },
          { 1LL, 8LL, 8LL },
          1, "0.0149074475234" },
        // class 2d-c16 key=[1/1 17/8 45/8]
        //   rep-vertices: (0,0) (4,1) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c16-d2, alpha-achieved=0.0812 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 17LL, 45LL },
          { 1LL, 8LL, 8LL },
          2, "0.00745549763517" },
        // class 2d-c16 key=[1/1 17/8 45/8]
        //   rep-vertices: (0,0) (4,1) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c16-d3, alpha-achieved=0.129 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 17LL, 45LL },
          { 1LL, 8LL, 8LL },
          3, "0.00488370500500" },
        // class 2d-c16 key=[1/1 17/8 45/8]
        //   rep-vertices: (0,0) (4,1) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c16-d4, alpha-achieved=0.179 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 17LL, 45LL },
          { 1LL, 8LL, 8LL },
          4, "0.00369234321618" },
        // class 2d-c16 key=[1/1 17/8 45/8]
        //   rep-vertices: (0,0) (4,1) (6,3)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c16-d5, alpha-achieved=0.226 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 17LL, 45LL },
          { 1LL, 8LL, 8LL },
          5, "0.00303736372501" },
        // class 2d-c17 key=[1/1 9/4 13/4]
        //   rep-vertices: (6,6) (4,3) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c17-d0, alpha-achieved=0.0102 (alpha ok, beta ok)
        { { 1LL, 9LL, 13LL },
          { 1LL, 4LL, 4LL },
          0, "0.0557137683986" },
        // class 2d-c17 key=[1/1 9/4 13/4]
        //   rep-vertices: (6,6) (4,3) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c17-d1, alpha-achieved=0.0329 (alpha ok, beta ok)
        { { 1LL, 9LL, 13LL },
          { 1LL, 4LL, 4LL },
          1, "0.0176040765152" },
        // class 2d-c17 key=[1/1 9/4 13/4]
        //   rep-vertices: (6,6) (4,3) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c17-d2, alpha-achieved=0.0645 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 9LL, 13LL },
          { 1LL, 4LL, 4LL },
          2, "0.00923656749823" },
        // class 2d-c17 key=[1/1 9/4 13/4]
        //   rep-vertices: (6,6) (4,3) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c17-d3, alpha-achieved=0.102 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 9LL, 13LL },
          { 1LL, 4LL, 4LL },
          3, "0.00606126815765" },
        // class 2d-c17 key=[1/1 9/4 13/4]
        //   rep-vertices: (6,6) (4,3) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c17-d4, alpha-achieved=0.143 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 9LL, 13LL },
          { 1LL, 4LL, 4LL },
          4, "0.00447944211044" },
        // class 2d-c17 key=[1/1 9/4 13/4]
        //   rep-vertices: (6,6) (4,3) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c17-d5, alpha-achieved=0.186 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 9LL, 13LL },
          { 1LL, 4LL, 4LL },
          5, "0.00356777618054" },
        // class 2d-c18 key=[1/1 37/16 45/16]
        //   rep-vertices: (18,5) (18,9) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c18-d0, alpha-achieved=0.00955 (alpha ok, beta ok)
        { { 1LL, 37LL, 45LL },
          { 1LL, 16LL, 16LL },
          0, "0.0591643870436" },
        // class 2d-c18 key=[1/1 37/16 45/16]
        //   rep-vertices: (18,5) (18,9) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c18-d1, alpha-achieved=0.0309 (alpha ok, beta ok)
        { { 1LL, 37LL, 45LL },
          { 1LL, 16LL, 16LL },
          1, "0.0186950340109" },
        // class 2d-c18 key=[1/1 37/16 45/16]
        //   rep-vertices: (18,5) (18,9) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c18-d2, alpha-achieved=0.0601 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 37LL, 45LL },
          { 1LL, 16LL, 16LL },
          2, "0.00986855249445" },
        // class 2d-c18 key=[1/1 37/16 45/16]
        //   rep-vertices: (18,5) (18,9) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c18-d3, alpha-achieved=0.0948 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 37LL, 45LL },
          { 1LL, 16LL, 16LL },
          3, "0.00646822126247" },
        // class 2d-c18 key=[1/1 37/16 45/16]
        //   rep-vertices: (18,5) (18,9) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c18-d4, alpha-achieved=0.133 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 37LL, 45LL },
          { 1LL, 16LL, 16LL },
          4, "0.00476254919971" },
        // class 2d-c18 key=[1/1 37/16 45/16]
        //   rep-vertices: (18,5) (18,9) (12,6)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c18-d5, alpha-achieved=0.174 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 37LL, 45LL },
          { 1LL, 16LL, 16LL },
          5, "0.00377547678151" },
        // class 2d-c19 key=[1/1 5/2 9/2]
        //   rep-vertices: (0,0) (1,1) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c19-d0, alpha-achieved=0.0108 (alpha ok, beta ok)
        { { 1LL, 5LL, 9LL },
          { 1LL, 2LL, 2LL },
          0, "0.0523092894690" },
        // class 2d-c19 key=[1/1 5/2 9/2]
        //   rep-vertices: (0,0) (1,1) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c19-d1, alpha-achieved=0.0351 (alpha ok, beta ok)
        { { 1LL, 5LL, 9LL },
          { 1LL, 2LL, 2LL },
          1, "0.0165167394085" },
        // class 2d-c19 key=[1/1 5/2 9/2]
        //   rep-vertices: (0,0) (1,1) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c19-d2, alpha-achieved=0.07 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 9LL },
          { 1LL, 2LL, 2LL },
          2, "0.00856125224720" },
        // class 2d-c19 key=[1/1 5/2 9/2]
        //   rep-vertices: (0,0) (1,1) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c19-d3, alpha-achieved=0.111 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 9LL },
          { 1LL, 2LL, 2LL },
          3, "0.00561572766624" },
        // class 2d-c19 key=[1/1 5/2 9/2]
        //   rep-vertices: (0,0) (1,1) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c19-d4, alpha-achieved=0.155 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 9LL },
          { 1LL, 2LL, 2LL },
          4, "0.00417697185880" },
        // class 2d-c19 key=[1/1 5/2 9/2]
        //   rep-vertices: (0,0) (1,1) (0,3)  prov=seed:K1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c19-d5, alpha-achieved=0.2 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 9LL },
          { 1LL, 2LL, 2LL },
          5, "0.00335893335813" },
        // class 2d-c20 key=[1/1 13/5 16/5]
        //   rep-vertices: (8,0) (8,4) (6,3)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c20-d0, alpha-achieved=0.00953 (alpha ok, beta ok)
        { { 1LL, 13LL, 16LL },
          { 1LL, 5LL, 5LL },
          0, "0.0593013627054" },
        // class 2d-c20 key=[1/1 13/5 16/5]
        //   rep-vertices: (8,0) (8,4) (6,3)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c20-d1, alpha-achieved=0.0308 (alpha ok, beta ok)
        { { 1LL, 13LL, 16LL },
          { 1LL, 5LL, 5LL },
          1, "0.0187310120140" },
        // class 2d-c20 key=[1/1 13/5 16/5]
        //   rep-vertices: (8,0) (8,4) (6,3)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c20-d2, alpha-achieved=0.06 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 16LL },
          { 1LL, 5LL, 5LL },
          2, "0.00988740606883" },
        // class 2d-c20 key=[1/1 13/5 16/5]
        //   rep-vertices: (8,0) (8,4) (6,3)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c20-d3, alpha-achieved=0.0945 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 16LL },
          { 1LL, 5LL, 5LL },
          3, "0.00648202283873" },
        // class 2d-c20 key=[1/1 13/5 16/5]
        //   rep-vertices: (8,0) (8,4) (6,3)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c20-d4, alpha-achieved=0.133 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 16LL },
          { 1LL, 5LL, 5LL },
          4, "0.00477482024225" },
        // class 2d-c20 key=[1/1 13/5 16/5]
        //   rep-vertices: (8,0) (8,4) (6,3)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c20-d5, alpha-achieved=0.173 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 16LL },
          { 1LL, 5LL, 5LL },
          5, "0.00378761707286" },
        // class 2d-c21 key=[1/1 13/5 4/1]
        //   rep-vertices: (8,0) (6,3) (4,2)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c21-d0, alpha-achieved=0.0103 (alpha ok, beta ok)
        { { 1LL, 13LL, 4LL },
          { 1LL, 5LL, 1LL },
          0, "0.0549009300019" },
        // class 2d-c21 key=[1/1 13/5 4/1]
        //   rep-vertices: (8,0) (6,3) (4,2)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c21-d1, alpha-achieved=0.0334 (alpha ok, beta ok)
        { { 1LL, 13LL, 4LL },
          { 1LL, 5LL, 1LL },
          1, "0.0173407284185" },
        // class 2d-c21 key=[1/1 13/5 4/1]
        //   rep-vertices: (8,0) (6,3) (4,2)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c21-d2, alpha-achieved=0.0657 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 4LL },
          { 1LL, 5LL, 1LL },
          2, "0.00907745949763" },
        // class 2d-c21 key=[1/1 13/5 4/1]
        //   rep-vertices: (8,0) (6,3) (4,2)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c21-d3, alpha-achieved=0.104 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 4LL },
          { 1LL, 5LL, 1LL },
          3, "0.00595619214848" },
        // class 2d-c21 key=[1/1 13/5 4/1]
        //   rep-vertices: (8,0) (6,3) (4,2)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c21-d4, alpha-achieved=0.145 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 4LL },
          { 1LL, 5LL, 1LL },
          4, "0.00440992167383" },
        // class 2d-c21 key=[1/1 13/5 4/1]
        //   rep-vertices: (8,0) (6,3) (4,2)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c21-d5, alpha-achieved=0.189 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 4LL },
          { 1LL, 5LL, 1LL },
          5, "0.00352268994364" },
        // class 2d-c22 key=[1/1 53/20 117/20]
        //   rep-vertices: (24,0) (22,7) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c22-d0, alpha-achieved=0.0115 (alpha ok, beta ok)
        { { 1LL, 53LL, 117LL },
          { 1LL, 20LL, 20LL },
          0, "0.0494259091069" },
        // class 2d-c22 key=[1/1 53/20 117/20]
        //   rep-vertices: (24,0) (22,7) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c22-d1, alpha-achieved=0.0373 (alpha ok, beta ok)
        { { 1LL, 53LL, 117LL },
          { 1LL, 20LL, 20LL },
          1, "0.0155883595746" },
        // class 2d-c22 key=[1/1 53/20 117/20]
        //   rep-vertices: (24,0) (22,7) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c22-d2, alpha-achieved=0.0757 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 53LL, 117LL },
          { 1LL, 20LL, 20LL },
          2, "0.00795588837409" },
        // class 2d-c22 key=[1/1 53/20 117/20]
        //   rep-vertices: (24,0) (22,7) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c22-d3, alpha-achieved=0.12 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 53LL, 117LL },
          { 1LL, 20LL, 20LL },
          3, "0.00520764673488" },
        // class 2d-c22 key=[1/1 53/20 117/20]
        //   rep-vertices: (24,0) (22,7) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c22-d4, alpha-achieved=0.167 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 53LL, 117LL },
          { 1LL, 20LL, 20LL },
          4, "0.00390493639569" },
        // class 2d-c22 key=[1/1 53/20 117/20]
        //   rep-vertices: (24,0) (22,7) (18,9)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c22-d5, alpha-achieved=0.214 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 53LL, 117LL },
          { 1LL, 20LL, 20LL },
          5, "0.00317969096369" },
        // class 2d-c23 key=[1/1 17/5 36/5]
        //   rep-vertices: (0,0) (6,0) (4,1)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c23-d0, alpha-achieved=0.0113 (alpha ok, beta ok)
        { { 1LL, 17LL, 36LL },
          { 1LL, 5LL, 5LL },
          0, "0.0500017458304" },
        // class 2d-c23 key=[1/1 17/5 36/5]
        //   rep-vertices: (0,0) (6,0) (4,1)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c23-d1, alpha-achieved=0.0368 (alpha ok, beta ok)
        { { 1LL, 17LL, 36LL },
          { 1LL, 5LL, 5LL },
          1, "0.0157626004505" },
        // class 2d-c23 key=[1/1 17/5 36/5]
        //   rep-vertices: (0,0) (6,0) (4,1)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c23-d2, alpha-achieved=0.0743 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 17LL, 36LL },
          { 1LL, 5LL, 5LL },
          2, "0.00809691135742" },
        // class 2d-c23 key=[1/1 17/5 36/5]
        //   rep-vertices: (0,0) (6,0) (4,1)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c23-d3, alpha-achieved=0.118 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 17LL, 36LL },
          { 1LL, 5LL, 5LL },
          3, "0.00529271243766" },
        // class 2d-c23 key=[1/1 17/5 36/5]
        //   rep-vertices: (0,0) (6,0) (4,1)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c23-d4, alpha-achieved=0.165 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 17LL, 36LL },
          { 1LL, 5LL, 5LL },
          4, "0.00396327305802" },
        // class 2d-c23 key=[1/1 17/5 36/5]
        //   rep-vertices: (0,0) (6,0) (4,1)  prov=seed:earclip:E1>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c23-d5, alpha-achieved=0.21 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 17LL, 36LL },
          { 1LL, 5LL, 5LL },
          5, "0.00322576176757" },
        // class 2d-c24 key=[1/1 61/16 117/16]
        //   rep-vertices: (24,0) (18,9) (18,5)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c24-d0, alpha-achieved=0.011 (alpha ok, beta ok)
        { { 1LL, 61LL, 117LL },
          { 1LL, 16LL, 16LL },
          0, "0.0515624074297" },
        // class 2d-c24 key=[1/1 61/16 117/16]
        //   rep-vertices: (24,0) (18,9) (18,5)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c24-d1, alpha-achieved=0.0357 (alpha ok, beta ok)
        { { 1LL, 61LL, 117LL },
          { 1LL, 16LL, 16LL },
          1, "0.0162616591596" },
        // class 2d-c24 key=[1/1 61/16 117/16]
        //   rep-vertices: (24,0) (18,9) (18,5)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c24-d2, alpha-achieved=0.0712 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 61LL, 117LL },
          { 1LL, 16LL, 16LL },
          2, "0.00842363586962" },
        // class 2d-c24 key=[1/1 61/16 117/16]
        //   rep-vertices: (24,0) (18,9) (18,5)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c24-d3, alpha-achieved=0.113 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 61LL, 117LL },
          { 1LL, 16LL, 16LL },
          3, "0.00550995704393" },
        // class 2d-c24 key=[1/1 61/16 117/16]
        //   rep-vertices: (24,0) (18,9) (18,5)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c24-d4, alpha-achieved=0.158 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 61LL, 117LL },
          { 1LL, 16LL, 16LL },
          4, "0.00410971071212" },
        // class 2d-c24 key=[1/1 61/16 117/16]
        //   rep-vertices: (24,0) (18,9) (18,5)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c24-d5, alpha-achieved=0.202 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 61LL, 117LL },
          { 1LL, 16LL, 16LL },
          5, "0.00332643869763" },
        // class 2d-c25 key=[1/1 4/1 5/1]
        //   rep-vertices: (0,0) (2,0) (2,1)  prov=seed:earclip:E1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c25-d0, alpha-achieved=0.00941 (alpha ok, beta ok)
        { { 1LL, 4LL, 5LL },
          { 1LL, 1LL, 1LL },
          0, "0.0600620216432" },
        // class 2d-c25 key=[1/1 4/1 5/1]
        //   rep-vertices: (0,0) (2,0) (2,1)  prov=seed:earclip:E1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c25-d1, alpha-achieved=0.0304 (alpha ok, beta ok)
        { { 1LL, 4LL, 5LL },
          { 1LL, 1LL, 1LL },
          1, "0.0189590446020" },
        // class 2d-c25 key=[1/1 4/1 5/1]
        //   rep-vertices: (0,0) (2,0) (2,1)  prov=seed:earclip:E1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c25-d2, alpha-achieved=0.0592 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 4LL, 5LL },
          { 1LL, 1LL, 1LL },
          2, "0.0100118897400" },
        // class 2d-c25 key=[1/1 4/1 5/1]
        //   rep-vertices: (0,0) (2,0) (2,1)  prov=seed:earclip:E1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c25-d3, alpha-achieved=0.0932 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 4LL, 5LL },
          { 1LL, 1LL, 1LL },
          3, "0.00656515824377" },
        // class 2d-c25 key=[1/1 4/1 5/1]
        //   rep-vertices: (0,0) (2,0) (2,1)  prov=seed:earclip:E1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c25-d4, alpha-achieved=0.131 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 4LL, 5LL },
          { 1LL, 1LL, 1LL },
          4, "0.00484010915599" },
        // class 2d-c25 key=[1/1 4/1 5/1]
        //   rep-vertices: (0,0) (2,0) (2,1)  prov=seed:earclip:E1
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c25-d5, alpha-achieved=0.17 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 4LL, 5LL },
          { 1LL, 1LL, 1LL },
          5, "0.00384478967942" },
        // class 2d-c26 key=[1/1 5/1 8/1]
        //   rep-vertices: (2,2) (0,0) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c26-d0, alpha-achieved=0.0103 (alpha ok, beta ok)
        { { 1LL, 5LL, 8LL },
          { 1LL, 1LL, 1LL },
          0, "0.0549907118829" },
        // class 2d-c26 key=[1/1 5/1 8/1]
        //   rep-vertices: (2,2) (0,0) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c26-d1, alpha-achieved=0.0333 (alpha ok, beta ok)
        { { 1LL, 5LL, 8LL },
          { 1LL, 1LL, 1LL },
          1, "0.0173516672927" },
        // class 2d-c26 key=[1/1 5/1 8/1]
        //   rep-vertices: (2,2) (0,0) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c26-d2, alpha-achieved=0.0656 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 8LL },
          { 1LL, 1LL, 1LL },
          2, "0.00909645432955" },
        // class 2d-c26 key=[1/1 5/1 8/1]
        //   rep-vertices: (2,2) (0,0) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c26-d3, alpha-achieved=0.104 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 8LL },
          { 1LL, 1LL, 1LL },
          3, "0.00595780523304" },
        // class 2d-c26 key=[1/1 5/1 8/1]
        //   rep-vertices: (2,2) (0,0) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c26-d4, alpha-achieved=0.145 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 8LL },
          { 1LL, 1LL, 1LL },
          4, "0.00441781008645" },
        // class 2d-c26 key=[1/1 5/1 8/1]
        //   rep-vertices: (2,2) (0,0) (2,1)  prov=nvb-closure
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c26-d5, alpha-achieved=0.188 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 5LL, 8LL },
          { 1LL, 1LL, 1LL },
          5, "0.00354474535486" },
        // class 2d-c27 key=[1/1 29/5 36/5]
        //   rep-vertices: (18,9) (18,3) (20,4)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c27-d0, alpha-achieved=0.00928 (alpha ok, beta ok)
        { { 1LL, 29LL, 36LL },
          { 1LL, 5LL, 5LL },
          0, "0.0608869248374" },
        // class 2d-c27 key=[1/1 29/5 36/5]
        //   rep-vertices: (18,9) (18,3) (20,4)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c27-d1, alpha-achieved=0.03 (alpha ok, beta ok)
        { { 1LL, 29LL, 36LL },
          { 1LL, 5LL, 5LL },
          1, "0.0192146506801" },
        // class 2d-c27 key=[1/1 29/5 36/5]
        //   rep-vertices: (18,9) (18,3) (20,4)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c27-d2, alpha-achieved=0.0584 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 29LL, 36LL },
          { 1LL, 5LL, 5LL },
          2, "0.0101528898702" },
        // class 2d-c27 key=[1/1 29/5 36/5]
        //   rep-vertices: (18,9) (18,3) (20,4)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c27-d3, alpha-achieved=0.0918 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 29LL, 36LL },
          { 1LL, 5LL, 5LL },
          3, "0.00665734902334" },
        // class 2d-c27 key=[1/1 29/5 36/5]
        //   rep-vertices: (18,9) (18,3) (20,4)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c27-d4, alpha-achieved=0.129 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 29LL, 36LL },
          { 1LL, 5LL, 5LL },
          4, "0.00490855847642" },
        // class 2d-c27 key=[1/1 29/5 36/5]
        //   rep-vertices: (18,9) (18,3) (20,4)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c27-d5, alpha-achieved=0.168 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 29LL, 36LL },
          { 1LL, 5LL, 5LL },
          5, "0.00390073144629" },
        // class 2d-c28 key=[1/1 25/4 45/4]
        //   rep-vertices: (4,3) (0,0) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c28-d0, alpha-achieved=0.0107 (alpha ok, beta ok)
        { { 1LL, 25LL, 45LL },
          { 1LL, 4LL, 4LL },
          0, "0.0528696212392" },
        // class 2d-c28 key=[1/1 25/4 45/4]
        //   rep-vertices: (4,3) (0,0) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c28-d1, alpha-achieved=0.0348 (alpha ok, beta ok)
        { { 1LL, 25LL, 45LL },
          { 1LL, 4LL, 4LL },
          1, "0.0166631168572" },
        // class 2d-c28 key=[1/1 25/4 45/4]
        //   rep-vertices: (4,3) (0,0) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c28-d2, alpha-achieved=0.0688 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 25LL, 45LL },
          { 1LL, 4LL, 4LL },
          2, "0.00869529225622" },
        // class 2d-c28 key=[1/1 25/4 45/4]
        //   rep-vertices: (4,3) (0,0) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c28-d3, alpha-achieved=0.109 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 25LL, 45LL },
          { 1LL, 4LL, 4LL },
          3, "0.00568220194424" },
        // class 2d-c28 key=[1/1 25/4 45/4]
        //   rep-vertices: (4,3) (0,0) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c28-d4, alpha-achieved=0.153 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 25LL, 45LL },
          { 1LL, 4LL, 4LL },
          4, "0.00422817567604" },
        // class 2d-c28 key=[1/1 25/4 45/4]
        //   rep-vertices: (4,3) (0,0) (6,3)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c28-d5, alpha-achieved=0.196 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 25LL, 45LL },
          { 1LL, 4LL, 4LL },
          5, "0.00341811053477" },
        // class 2d-c29 key=[1/1 32/5 9/1]
        //   rep-vertices: (20,4) (18,3) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c29-d0, alpha-achieved=0.00975 (alpha ok, beta ok)
        { { 1LL, 32LL, 9LL },
          { 1LL, 5LL, 1LL },
          0, "0.0579799485519" },
        // class 2d-c29 key=[1/1 32/5 9/1]
        //   rep-vertices: (20,4) (18,3) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c29-d1, alpha-achieved=0.0316 (alpha ok, beta ok)
        { { 1LL, 32LL, 9LL },
          { 1LL, 5LL, 1LL },
          1, "0.0182971100167" },
        // class 2d-c29 key=[1/1 32/5 9/1]
        //   rep-vertices: (20,4) (18,3) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c29-d2, alpha-achieved=0.0616 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 32LL, 9LL },
          { 1LL, 5LL, 1LL },
          2, "0.00964442368366" },
        // class 2d-c29 key=[1/1 32/5 9/1]
        //   rep-vertices: (20,4) (18,3) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c29-d3, alpha-achieved=0.0972 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 32LL, 9LL },
          { 1LL, 5LL, 1LL },
          3, "0.00632148041900" },
        // class 2d-c29 key=[1/1 32/5 9/1]
        //   rep-vertices: (20,4) (18,3) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c29-d4, alpha-achieved=0.136 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 32LL, 9LL },
          { 1LL, 5LL, 1LL },
          4, "0.00467258320836" },
        // class 2d-c29 key=[1/1 32/5 9/1]
        //   rep-vertices: (20,4) (18,3) (24,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c29-d5, alpha-achieved=0.177 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 32LL, 9LL },
          { 1LL, 5LL, 1LL },
          5, "0.00373002427415" },
        // class 2d-c30 key=[1/1 9/1 10/1]
        //   rep-vertices: (3,2) (6,3) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c30-d0, alpha-achieved=0.00875 (alpha ok, beta ok)
        { { 1LL, 9LL, 10LL },
          { 1LL, 1LL, 1LL },
          0, "0.0645259117571" },
        // class 2d-c30 key=[1/1 9/1 10/1]
        //   rep-vertices: (3,2) (6,3) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c30-d1, alpha-achieved=0.0283 (alpha ok, beta ok)
        { { 1LL, 9LL, 10LL },
          { 1LL, 1LL, 1LL },
          1, "0.0203523325422" },
        // class 2d-c30 key=[1/1 9/1 10/1]
        //   rep-vertices: (3,2) (6,3) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c30-d2, alpha-achieved=0.0549 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 9LL, 10LL },
          { 1LL, 1LL, 1LL },
          2, "0.0107651373721" },
        // class 2d-c30 key=[1/1 9/1 10/1]
        //   rep-vertices: (3,2) (6,3) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c30-d3, alpha-achieved=0.0862 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 9LL, 10LL },
          { 1LL, 1LL, 1LL },
          3, "0.00705632802042" },
        // class 2d-c30 key=[1/1 9/1 10/1]
        //   rep-vertices: (3,2) (6,3) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c30-d4, alpha-achieved=0.121 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 9LL, 10LL },
          { 1LL, 1LL, 1LL },
          4, "0.00519394613772" },
        // class 2d-c30 key=[1/1 9/1 10/1]
        //   rep-vertices: (3,2) (6,3) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c30-d5, alpha-achieved=0.157 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 9LL, 10LL },
          { 1LL, 1LL, 1LL },
          5, "0.00411736443914" },
        // class 2d-c31 key=[1/1 37/4 45/4]
        //   rep-vertices: (6,1) (6,3) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c31-d0, alpha-achieved=0.00912 (alpha ok, beta ok)
        { { 1LL, 37LL, 45LL },
          { 1LL, 4LL, 4LL },
          0, "0.0619850563240" },
        // class 2d-c31 key=[1/1 37/4 45/4]
        //   rep-vertices: (6,1) (6,3) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c31-d1, alpha-achieved=0.0295 (alpha ok, beta ok)
        { { 1LL, 37LL, 45LL },
          { 1LL, 4LL, 4LL },
          1, "0.0195573283870" },
        // class 2d-c31 key=[1/1 37/4 45/4]
        //   rep-vertices: (6,1) (6,3) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c31-d2, alpha-achieved=0.0572 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 37LL, 45LL },
          { 1LL, 4LL, 4LL },
          2, "0.0103407429380" },
        // class 2d-c31 key=[1/1 37/4 45/4]
        //   rep-vertices: (6,1) (6,3) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c31-d3, alpha-achieved=0.09 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 37LL, 45LL },
          { 1LL, 4LL, 4LL },
          3, "0.00678015959240" },
        // class 2d-c31 key=[1/1 37/4 45/4]
        //   rep-vertices: (6,1) (6,3) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c31-d4, alpha-achieved=0.126 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 37LL, 45LL },
          { 1LL, 4LL, 4LL },
          4, "0.00499789192923" },
        // class 2d-c31 key=[1/1 37/4 45/4]
        //   rep-vertices: (6,1) (6,3) (0,0)  prov=nvb-closure>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c31-d5, alpha-achieved=0.164 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 37LL, 45LL },
          { 1LL, 4LL, 4LL },
          5, "0.00397138071378" },
        // class 2d-c32 key=[1/1 13/1 18/1]
        //   rep-vertices: (0,0) (3,2) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c32-d0, alpha-achieved=0.0096 (alpha ok, beta ok)
        { { 1LL, 13LL, 18LL },
          { 1LL, 1LL, 1LL },
          0, "0.0588676218232" },
        // class 2d-c32 key=[1/1 13/1 18/1]
        //   rep-vertices: (0,0) (3,2) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c32-d1, alpha-achieved=0.0311 (alpha ok, beta ok)
        { { 1LL, 13LL, 18LL },
          { 1LL, 1LL, 1LL },
          1, "0.0185692752407" },
        // class 2d-c32 key=[1/1 13/1 18/1]
        //   rep-vertices: (0,0) (3,2) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c32-d2, alpha-achieved=0.0606 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 18LL },
          { 1LL, 1LL, 1LL },
          2, "0.00980339963465" },
        // class 2d-c32 key=[1/1 13/1 18/1]
        //   rep-vertices: (0,0) (3,2) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c32-d3, alpha-achieved=0.0955 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 18LL },
          { 1LL, 1LL, 1LL },
          3, "0.00642329235261" },
        // class 2d-c32 key=[1/1 13/1 18/1]
        //   rep-vertices: (0,0) (3,2) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c32-d4, alpha-achieved=0.134 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 18LL },
          { 1LL, 1LL, 1LL },
          4, "0.00474529271456" },
        // class 2d-c32 key=[1/1 13/1 18/1]
        //   rep-vertices: (0,0) (3,2) (3,3)  prov=seed:earclip:E2>alf
        // source: CONST-B2a gen on SekineMainCorei7-11700 2026-08-10, gates G-C3/G-C5, L=3,
        //   log-ref: constb2a_gen_SekineMainCorei7-11700.txt job-2d-c32-d5, alpha-achieved=0.173 (alpha ok; beta not met at L=3 -- B-2b' sharpening planned)
        { { 1LL, 13LL, 18LL },
          { 1LL, 1LL, 1LL },
          5, "0.00379001450031" },
    };
    count = static_cast<int>(sizeof entries / sizeof entries[0]);
    return entries;
}

} // namespace detail

// VCP_CONSTANTS_TABLE_END

inline const l2_projection_registry_entry_2d*
l2_projection_registry_2d_table(int& count) {
    return detail::l2_projection_registry_2d_entries(count);
}

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_DICT_REGISTRY_2D_HPP
