#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cell.h"
#include "delaunay.h"
#include "hall_symbol.h"
#include "mathfunc.h"
#include "niggli.h"
#include "magnetic_primitive.h"
#include "magnetic_symmetry.h"
#include "mspacegroup.h"
#include "magnetic_cell.h"
#include "pointgroup.h"
#include "primitive.h"
#include "spacegroup.h"
#include "spg_database.h"
#include "symmetry.h"
#include "debug.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define REDUCE_RATE 0.95
#define NUM_ATTEMPT 100
#define INT_PREC 0.1
#define ZERO_PREC 1e-10

static double change_of_basis_monocli[36][3][3] = {
  {{ 1, 0, 0 }, /* b  first turn; two axes are flipped in second turn */
   { 0, 1, 0 },
   { 0, 0, 1 }},
  {{ 0, 0, 1 }, /* b */
   { 0,-1, 0 },
   { 1, 0, 0 }},
  {{ 0, 0, 1 }, /* a */
   { 1, 0, 0 },
   { 0, 1, 0 }},
  {{ 1, 0, 0 }, /* c */
   { 0, 0, 1 },
   { 0,-1, 0 }},
  {{ 0, 1, 0 }, /* c */
   { 0, 0, 1 },
   { 1, 0, 0 }},
  {{ 0,-1, 0 }, /* a */
   { 1, 0, 0 },
   { 0, 0, 1 }},
  {{-1, 0, 1 }, /* b */
   { 0, 1, 0 },
   {-1, 0, 0 }},
  {{ 1, 0,-1 }, /* b */
   { 0,-1, 0 },
   { 0, 0,-1 }},
  {{ 0, 1,-1 }, /* a */
   { 1, 0, 0 },
   { 0, 0,-1 }},
  {{-1,-1, 0 }, /* c */
   { 0, 0, 1 },
   {-1, 0, 0 }},
  {{ 1,-1, 0 }, /* c */
   { 0, 0, 1 },
   { 0,-1, 0 }},
  {{ 0, 1, 1 }, /* a */
   { 1, 0, 0 },
   { 0, 1, 0 }},
  {{ 0, 0,-1 }, /* b */
   { 0, 1, 0 },
   { 1, 0,-1 }},
  {{-1, 0, 0 }, /* b */
   { 0,-1, 0 },
   {-1, 0, 1 }},
  {{ 0,-1, 0 }, /* a */
   { 1, 0, 0 },
   { 0,-1, 1 }},
  {{ 0, 1, 0 }, /* c */
   { 0, 0, 1 },
   { 1, 1, 0 }},
  {{-1, 0, 0 }, /* c */
   { 0, 0, 1 },
   {-1, 1, 0 }},
  {{ 0, 0,-1 }, /* a */
   { 1, 0, 0 },
   { 0,-1,-1 }},
  {{ 1, 0, 0 }, /* b  two axes are flipped to look for non-acute axes */
   { 0,-1, 0 },
   { 0, 0,-1 }},
  {{ 0, 0,-1 }, /* b */
   { 0, 1, 0 },
   { 1, 0, 0 }},
  {{ 0, 0, 1 }, /* a */
   {-1, 0, 0 },
   { 0,-1, 0 }},
  {{-1, 0, 0 }, /* c */
   { 0, 0,-1 },
   { 0,-1, 0 }},
  {{ 0, 1, 0 }, /* c */
   { 0, 0,-1 },
   {-1, 0, 0 }},
  {{ 0, 1, 0 }, /* a */
   {-1, 0, 0 },
   { 0, 0, 1 }},
  {{-1, 0,-1 }, /* b */
   { 0,-1, 0 },
   {-1, 0, 0 }},
  {{ 1, 0, 1 }, /* b */
   { 0, 1, 0 },
   { 0, 0, 1 }},
  {{ 0,-1,-1 }, /* a */
   {-1, 0, 0 },
   { 0, 0,-1 }},
  {{ 1,-1, 0 }, /* c */
   { 0, 0,-1 },
   { 1, 0, 0 }},
  {{-1,-1, 0 }, /* c */
   { 0, 0,-1 },
   { 0,-1, 0 }},
  {{ 0,-1, 1 }, /* a */
   {-1, 0, 0 },
   { 0,-1, 0 }},
  {{ 0, 0, 1 }, /* b */
   { 0,-1, 0 },
   { 1, 0, 1 }},
  {{-1, 0, 0 }, /* b */
   { 0, 1, 0 },
   {-1, 0,-1 }},
  {{ 0, 1, 0 }, /* a */
   {-1, 0, 0 },
   { 0, 1, 1 }},
  {{ 0, 1, 0 }, /* c */
   { 0, 0,-1 },
   {-1, 1, 0 }},
  {{ 1, 0, 0 }, /* c */
   { 0, 0,-1 },
   { 1, 1, 0 }},
  {{ 0, 0,-1 }, /* a */
   {-1, 0, 0 },
   { 0, 1,-1 }}};

static Centering change_of_centering_monocli[36] = {
  C_FACE, /* first turn */
  A_FACE,
  B_FACE,
  B_FACE,
  A_FACE,
  C_FACE,
  A_FACE,
  C_FACE,
  C_FACE,
  A_FACE,
  B_FACE,
  B_FACE,
  BODY,
  BODY,
  BODY,
  BODY,
  BODY,
  BODY,
  C_FACE, /* second turn */
  A_FACE,
  B_FACE,
  B_FACE,
  A_FACE,
  C_FACE,
  A_FACE,
  C_FACE,
  C_FACE,
  A_FACE,
  B_FACE,
  B_FACE,
  BODY,
  BODY,
  BODY,
  BODY,
  BODY,
  BODY};

static int change_of_unique_axis_monocli[36] = {
  1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0,
  1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0};

static double change_of_basis_ortho[6][3][3] = {{{ 1, 0, 0 },
                                                 { 0, 1, 0 },
                                                 { 0, 0, 1 }},
                                                {{ 0, 0, 1 },
                                                 { 1, 0, 0 },
                                                 { 0, 1, 0 }},
                                                {{ 0, 1, 0 },
                                                 { 0, 0, 1 },
                                                 { 1, 0, 0 }},
                                                {{ 0, 1, 0 },
                                                 { 1, 0, 0 },
                                                 { 0, 0,-1 }},
                                                {{ 1, 0, 0 },
                                                 { 0, 0, 1 },
                                                 { 0,-1, 0 }},
                                                {{ 0, 0, 1 },
                                                 { 0, 1, 0 },
                                                 {-1, 0, 0 }}};

static Centering change_of_centering_ortho[6] = {C_FACE,
                                                 B_FACE,
                                                 A_FACE,
                                                 C_FACE,
                                                 B_FACE,
                                                 A_FACE};
static int change_of_unique_axis_ortho[6] = {2, 1, 0, 2, 1, 0};

/* n_l : the index of L(g) in N_\epsilon(g) of SPG No.16-74 */
/* See ITA: Affine normalizer or highest symmetry Euclidean normalizer */
/* Previous implementation below was not correct for 67, 68, 73, 74, */
/* Cmma, Ccca, Ibca, Imma */
/* 6 / ((Number of hall symbols of each spg-type) / x) */
/* where x=1 and x=2 with without and with centring. */
/* 6, 2, 2, 6, 2, */
/* 2, 6, 6, 6, 2, 1, 2, 1, 1, 1, */
/* 1, 2, 1, 2, 2, 1, 2, 1, 1, 1, */
/* 1, 2, 2, 2, 2, 1, 6, 6, 2, 2, */
/* 1, 1, 1, 1, 2, 2, 1, 2, 2, 1, */
/* 3, 1, 1, 1, 2, 2, 1, 1, 6, 6, */
/* 6, 2, 3, 1 */
static int num_axis_choices_ortho[59] = {
  6, 2, 2, 6, 2,  /* 16-20 */
  2, 6, 6, 6, 2, 1, 2, 1, 1, 1,  /* 21-30 */
  1, 2, 1, 2, 2, 1, 2, 1, 1, 1,  /* 31-40 */
  1, 2, 2, 2, 2, 1, 6, 6, 2, 2,  /* 41-50 */
  1, 1, 1, 1, 2, 2, 1, 2, 2, 1,  /* 51-60 */
  3, 1, 1, 1, 2, 2, 2, 2, 6, 6,  /* 61-70 */
  6, 2, 6, 2};  /* 71-74 */

/* hR_to_hP is the inverse of rhombo_obverse and used until the */
/* commit  e0398ba9. But now this is included in */
/* change_of_basis_rhombo in a way of */
/* dot(change_of_basis_rhombo_hex,hR_to_hP). */
/* static double hR_to_hP[3][3] = {{ 1, 0, 1 }, */
/*                                 {-1, 1, 1 }, */
/*                                 { 0,-1, 1 }}; */

static double change_of_basis_rhombo[6][3][3] = {{{ 1, 0, 0 },
                                                  { 0, 1, 0 },
                                                  { 0, 0, 1 }},
                                                 {{ 0, 0, 1 },
                                                  { 1, 0, 0 },
                                                  { 0, 1, 0 }},
                                                 {{ 0, 1, 0 },
                                                  { 0, 0, 1 },
                                                  { 1, 0, 0 }},
                                                 {{ 0, 0,-1 },
                                                  { 0,-1, 0 },
                                                  {-1, 0, 0 }},
                                                 {{ 0,-1, 0 },
                                                  {-1, 0, 0 },
                                                  { 0, 0,-1 }},
                                                 {{-1, 0, 0 },
                                                  { 0, 0,-1 },
                                                  { 0,-1, 0 }}};

static double change_of_basis_rhombo_hex[6][3][3] = {{{ 1,  0,  1},
                                                      {-1,  1,  1},
                                                      { 0, -1,  1}},
                                                     {{ 0, -1,  1},
                                                      { 1,  0,  1},
                                                      {-1,  1,  1}},
                                                     {{-1,  1,  1},
                                                      { 0, -1,  1},
                                                      { 1,  0,  1}},
                                                     {{ 0,  1, -1},
                                                      { 1, -1, -1},
                                                      {-1,  0, -1}},
                                                     {{ 1, -1, -1},
                                                      {-1,  0, -1},
                                                      { 0,  1, -1}},
                                                     {{-1,  0, -1},
                                                      { 0,  1, -1},
                                                      { 1, -1, -1}}};

static double change_of_basis_C4[8][3][3] = {{{ 1, 0, 0 },
                                              { 0, 1, 0 },
                                              { 0, 0, 1 }},
                                             {{ 0,-1, 0 },
                                              { 1, 0, 0 },
                                              { 0, 0, 1 }},
                                             {{-1, 0, 0 },
                                              { 0,-1, 0 },
                                              { 0, 0, 1 }},
                                             {{ 0, 1, 0 },
                                              {-1, 0, 0 },
                                              { 0, 0, 1 }},
                                             {{ 0, 1, 0 },
                                              { 1, 0, 0 },
                                              { 0, 0,-1 }},
                                             {{-1, 0, 0 },
                                              { 0, 1, 0 },
                                              { 0, 0,-1 }},
                                             {{ 0,-1, 0 },
                                              {-1, 0, 0 },
                                              { 0, 0,-1 }},
                                             {{ 1, 0, 0 },
                                              { 0,-1, 0 },
                                              { 0, 0,-1 }}};

static double change_of_basis_C6[12][3][3] = {{{ 1, 0, 0 }, /* 1 */
                                               { 0, 1, 0 },
                                               { 0, 0, 1 }},
                                              {{ 1,-1, 0 }, /* 2 */
                                               { 1, 0, 0 },
                                               { 0, 0, 1 }},
                                              {{ 0,-1, 0 }, /* 3 */
                                               { 1,-1, 0 },
                                               { 0, 0, 1 }},
                                              {{-1, 0, 0 }, /* 4 */
                                               { 0,-1, 0 },
                                               { 0, 0, 1 }},
                                              {{-1, 1, 0 }, /* 5 */
                                               {-1, 0, 0 },
                                               { 0, 0, 1 }},
                                              {{ 0, 1, 0 }, /* 6 */
                                               {-1, 1, 0 },
                                               { 0, 0, 1 }},
                                              {{ 0, 1, 0 }, /* 7 */
                                               { 1, 0, 0 },
                                               { 0, 0,-1 }},
                                              {{-1, 1, 0 }, /* 8 */
                                               { 0, 1, 0 },
                                               { 0, 0,-1 }},
                                              {{-1, 0, 0 }, /* 9 */
                                               {-1, 1, 0 },
                                               { 0, 0,-1 }},
                                              {{ 0,-1, 0 }, /* 10 */
                                               {-1, 0, 0 },
                                               { 0, 0,-1 }},
                                              {{ 1,-1, 0 }, /* 11 */
                                               { 0,-1, 0 },
                                               { 0, 0,-1 }},
                                              {{ 1, 0, 0 }, /* 12 */
                                               { 1,-1, 0 },
                                               { 0, 0,-1 }}};

/* Removed after commit 67997654 */
/* static double change_of_basis_501[3][3] = {{ 0, 0, 1}, */
/*                                            { 0,-1, 0}, */
/*                                            { 1, 0, 0}}; */

/* It is BSC standard conventions to mapping from spacegroup name to hall numbers,
 which is not the one used in spglib */

static int spacegroup_to_hall_number[230] = {
     1,   2,   3,   6,   9,  18,  21,  30,  39,  57,
    60,  63,  72,  81,  90, 108, 109, 112, 115, 116,
   119, 122, 123, 124, 125, 128, 134, 137, 143, 149,
   155, 161, 164, 170, 173, 176, 182, 185, 191, 197,
   203, 209, 212, 215, 218, 221, 227, 229, 230, 234,
   239, 245, 251, 257, 263, 266, 269, 275, 279, 284,
   290, 292, 298, 304, 310, 313, 316, 323, 334, 336,
   337, 338, 341, 343, 349, 350, 351, 352, 353, 354,
   355, 356, 357, 358, 360, 362, 363, 365, 366, 367,
   368, 369, 370, 371, 372, 373, 374, 375, 376, 377,
   378, 379, 380, 381, 382, 383, 384, 385, 386, 387,
   388, 389, 390, 391, 392, 393, 394, 395, 396, 397,
   398, 399, 400, 401, 403, 405, 406, 407, 409, 411,
   412, 413, 415, 417, 418, 419, 421, 423, 424, 425,
   427, 429, 430, 431, 432, 433, 435, 436, 438, 439,
   440, 441, 442, 443, 444, 446, 447, 448, 449, 450,
   452, 454, 455, 456, 457, 458, 460, 462, 463, 464,
   465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
   475, 476, 477, 478, 479, 480, 481, 482, 483, 484,
   485, 486, 487, 488, 489, 490, 491, 492, 493, 494,
   496, 497, 499, 500, 501, 502, 503, 504, 505, 506,
   507, 508, 509, 510, 511, 512, 513, 514, 515, 516,
   517, 519, 520, 522, 523, 524, 526, 528, 529, 530,
};

static double identity[3][3] = {{ 1, 0, 0 },
                                { 0, 1, 0 },
                                { 0, 0, 1 }};
static double monocli_i2c[3][3] = {{ 1, 0,-1 },
                                   { 0, 1, 0 },
                                   { 1, 0, 0 }};
static double monocli_a2c[3][3] = {{ 0, 0, 1 },
                                   { 0,-1, 0 },
                                   { 1, 0, 0 }};
static double rhombo_obverse[3][3] = {{ 2./3,-1./3,-1./3 },
                                      { 1./3, 1./3,-2./3 },
                                      { 1./3, 1./3, 1./3 }};
static double rhomb_reverse[3][3] = {{ 1./3,-2./3, 1./3 },
                                     { 2./3,-1./3,-1./3 },
                                     { 1./3, 1./3, 1./3 }};
static double a2c[3][3] = {{ 0, 0, 1 },
                           { 1, 0, 0 },
                           { 0, 1, 0 }};
static double b2c[3][3] = {{ 0, 1, 0 },
                           { 0, 0, 1 },
                           { 1, 0, 0 }};

static double A_mat[3][3] = {{    1,    0,    0},
                             {    0, 1./2,-1./2},
                             {    0, 1./2, 1./2}};
static double C_mat[3][3] = {{ 1./2, 1./2,    0},
                             {-1./2, 1./2,    0},
                             {    0,    0,    1}};
static double R_mat[3][3] = {{ 2./3,-1./3,-1./3 },
                             { 1./3, 1./3,-2./3 },
                             { 1./3, 1./3, 1./3 }};
static double I_mat[3][3] = {{-1./2, 1./2, 1./2 },
                             { 1./2,-1./2, 1./2 },
                             { 1./2, 1./2,-1./2 }};
static double F_mat[3][3] = {{    0, 1./2, 1./2 },
                             { 1./2,    0, 1./2 },
                             { 1./2, 1./2,    0 }};

static Spacegroup * search_spacegroup_with_symmetry(const Primitive * primitive,
                                                    const int candidates[],
                                                    const int num_candidates,
                                                    const Symmetry *symmetry,
                                                    const double symprec,
                                                    const double angle_tolerance);
static Spacegroup * get_spacegroup(const int hall_number,
                                   const double origin_shift[3],
                                   SPGCONST double conv_lattice[3][3]);
static int iterative_search_hall_number(double origin_shift[3],
                                        double conv_lattice[3][3],
                                        const int candidates[],
                                        const int num_candidates,
                                        const Primitive * primitive,
                                        const Symmetry * symmetry,
                                        const double symprec,
                                        const double angle_tolerance);
static int change_basis_tricli(int tmat_int[3][3],
                               SPGCONST double conv_lattice[3][3],
                               SPGCONST double primitive_lattice[3][3],
                               const double symprec);
static int change_basis_monocli(int tmat_int[3][3],
                                SPGCONST double conv_lattice[3][3],
                                SPGCONST double primitive_lattice[3][3],
                                const double symprec);
static Symmetry *
get_initial_conventional_symmetry(const Centering centering,
                                  SPGCONST double tmat[3][3],
                                  const Symmetry * symmetry);
static int search_hall_number(double origin_shift[3],
                              double conv_lattice[3][3],
                              const int candidates[],
                              const int num_candidates,
                              const Primitive *primitive,
                              const Symmetry * symmetry,
                              const double symprec);
static int match_hall_symbol_db(double origin_shift[3],
                                double conv_lattice[3][3],
                                SPGCONST double (*orig_lattice)[3],
                                const int hall_number,
                                const int pointgroup_number,
                                const Holohedry holohedry,
                                const Centering centering,
                                const Symmetry *symmetry,
                                const double symprec);
static int match_hall_symbol_db_monocli(double origin_shift[3],
                                        double conv_lattice[3][3],
                                        SPGCONST double (*orig_lattice)[3],
                                        const int hall_number,
                                        const int space_group_number,
                                        const Centering centering,
                                        const Symmetry *conv_symmetry,
                                        const double symprec);
static int match_hall_symbol_db_monocli_in_loop(double origin_shift[3],
                                                double conv_lattice[3][3],
                                                double norms_squared[2],
                                                const int i,
                                                SPGCONST double (*orig_lattice)[3],
                                                const int check_norms,
                                                const int hall_number,
                                                const Centering centering,
                                                const Symmetry *conv_symmetry,
                                                const double symprec);
static int match_hall_symbol_db_ortho(double origin_shift[3],
                                      double conv_lattice[3][3],
                                      SPGCONST double (*orig_lattice)[3],
                                      const int hall_number,
                                      const Centering centering,
                                      const Symmetry *symmetry,
                                      const int num_free_axes,
                                      const double symprec);
static int match_hall_symbol_db_ortho_in_loop(double origin_shift[3],
                                              double lattice[3][3],
                                              SPGCONST double (*orig_lattice)[3],
                                              const int i,
                                              const int hall_number,
                                              const Centering centering,
                                              const Symmetry *symmetry,
                                              const int num_free_axes,
                                              const double symprec);
static int match_hall_symbol_db_cubic(double origin_shift[3],
                                      double conv_lattice[3][3],
                                      SPGCONST double (*orig_lattice)[3],
                                      const int hall_number,
                                      const Centering centering,
                                      const Symmetry *conv_symmetry,
                                      const double symprec);
static int match_hall_symbol_db_cubic_in_loop(double origin_shift[3],
                                              double conv_lattice[3][3],
                                              SPGCONST double (*orig_lattice)[3],
                                              const int i,
                                              const int hall_number,
                                              const Centering centering,
                                              const Symmetry *conv_symmetry,
                                              const double symprec);
static int match_hall_symbol_db_rhombo(double origin_shift[3],
                                       double conv_lattice[3][3],
                                       SPGCONST double (*orig_lattice)[3],
                                       const int hall_number,
                                       const Symmetry *conv_symmetry,
                                       const double symprec);
static int match_hall_symbol_db_others(double origin_shift[3],
                                       double conv_lattice[3][3],
                                       SPGCONST double (*orig_lattice)[3],
                                       const int hall_number,
                                       const Centering centering,
                                       const Holohedry holohedry,
                                       const Symmetry *conv_symmetry,
                                       const double symprec);
static int match_hall_symbol_db_change_of_basis_loop(
  double origin_shift[3],
  double conv_lattice[3][3],
  SPGCONST double (*orig_lattice)[3],
  SPGCONST double (*change_of_basis)[3][3],
  const int num_change_of_basis,
  const int hall_number,
  const Centering centering,
  const Symmetry *conv_symmetry,
  const double symprec);
static Symmetry * get_conventional_symmetry(SPGCONST double tmat[3][3],
                                            const Centering centering,
                                            const Symmetry *primitive_sym);
static Centering get_centering(double correction_mat[3][3],
                               SPGCONST int tmat[3][3],
                               const Laue laue);
static Centering get_base_center(SPGCONST int tmat[3][3]);
static int get_centering_shifts(double shift[3][3],
                                const Centering centering);
static int is_equivalent_lattice(double tmat[3][3],
                                 const int allow_flip,
                                 SPGCONST double lattice[3][3],
                                 SPGCONST double orig_lattice[3][3],
                                 const double symprec);




/* modified part start */

int magnetic_type_check( const Magnetic_Symmetry * symmetry, const double symprec);
static int check_vector(const int j,
                        const int k,
                        const double vectors[][3],
                        SPGCONST int rot[3][3],
                        SPGCONST double lattice[3][3],
                        const int T,
                        const double symprec);

static int mat_check_identity_vector_d3(const double a[3], const double b[3], const double symprec);

static Magnetic_Symmetry * get_magnetic_operations(const Symmetry *sym_nonspin,
                                                   const Magnetic_Cell *cell,
                                                   const double symprec);

Magnetic_Symmetry * sym_get_magnetic_operation(const Magnetic_Cell * primitive,
                                               const double symprec,
                                               const double angle_tolerance);

static Symmetry * nonmagnetic_symmetries_for_standard_cell(const Magnetic_Symmetry * symmetry, const int msg_type);


Spacegroup * spa_search_mspacegroup(const Magnetic_Primitive * primitive,
                                    const int hall_number,
                                    const double symprec,
                                    const double angle_tolerance)
{
  Spacegroup *spacegroup;
  Magnetic_Symmetry *symmetry;
  Symmetry *symmetry_new;
  Primitive *primitive_new;
  int mtype;  
  int candidate[1];

  debug_print("search_spacegroup (tolerance = %f):\n", symprec);

  symmetry = NULL;
  spacegroup = NULL;
  symmetry_new = NULL;
  primitive_new = NULL;

  if ((symmetry =
       sym_get_magnetic_operation(primitive->cell, symprec, angle_tolerance)) == NULL) {
    return NULL;
  }
  
  mtype = magnetic_type_check( symmetry, symprec);
  if ((symmetry_new =
       nonmagnetic_symmetries_for_standard_cell(symmetry,mtype)) == NULL) {
    return NULL;
  }

  if ((primitive_new =
       Magnetic_Primitive_to_Primitive(primitive)) == NULL) {
    return NULL;
  }

  if (hall_number > 0) {
    candidate[0] = hall_number;
  }
  



  if (mtype > 0 && mtype < 5)
  {


  if (hall_number) {
    spacegroup = search_spacegroup_with_symmetry(primitive_new,
                                                 candidate,
                                                 1,
                                                 symmetry_new,
                                                 symprec,
                                                 angle_tolerance);
  } else {
    spacegroup = search_spacegroup_with_symmetry(primitive_new,
                                                 spacegroup_to_hall_number,
                                                 230,
                                                 symmetry_new,
                                                 symprec,
                                                 angle_tolerance);
  }
  
  sym_free_magnetic_symmetry(symmetry);
  sym_free_symmetry(symmetry_new);
  prm_free_primitive(primitive_new);
  symmetry = NULL;
  symmetry_new = NULL;
  primitive_new =  NULL;
  } else {
    spacegroup = NULL; 
  }

  return spacegroup;
}

Magnetic_Primitive * standard_primitive_cell(const Magnetic_Primitive * primitive, SPGCONST Spacegroup * spacegroup)
{
  Magnetic_Primitive *primitive_new;
  Magnetic_Cell *conv_prim;
  Centering centering;
  int i,j;
  double origin_shift_prim[3];
  double Cmat_inv[3][3], Cmat[3][3];
  double Pmat[3][3];
  double trans_mat[3][3],trans_mat_inv[3][3];
  double Pmat_inv[3][3] = {{1.0,0.0,0.0},
                           {0.0,1.0,0.0},
                           {0.0,0.0,1.0}};

  centering = -10;

  if ((primitive_new = magn_prm_alloc_primitive(primitive->size))==NULL)
  {
    return NULL;
  }

  centering = spgdb_get_spacegroup_type(spacegroup->hall_number).centering;

  switch (centering) {
  case PRIMITIVE:
    break;
  case A_FACE:
    mat_copy_matrix_d3(Pmat_inv, A_mat);
    break;
  case C_FACE:
    mat_copy_matrix_d3(Pmat_inv, C_mat);
    break;
  case FACE:
    mat_copy_matrix_d3(Pmat_inv, F_mat);
    break;
  case BODY:
    mat_copy_matrix_d3(Pmat_inv, I_mat);
    break;
  case R_CENTER:
    mat_copy_matrix_d3(Pmat_inv, R_mat);
    break;
  default:
    primitive_new = NULL;
    goto ret;
  }

  if (! mat_inverse_matrix_d3(Cmat_inv,spacegroup->bravais_lattice,0))
  {
    return NULL;
  }


  if (!mat_inverse_matrix_d3(Pmat,Pmat_inv,0))
  {
    goto ret;
  }

  mat_multiply_matrix_d3(Cmat,Cmat_inv,primitive->cell->lattice);

  mat_multiply_matrix_d3(trans_mat,Pmat,Cmat);


    if ((conv_prim = magnetic_cel_alloc_cell(primitive->cell->size)) == NULL) {
    return NULL;
  }
 


 mat_multiply_matrix_vector_d3(origin_shift_prim,Pmat,spacegroup->origin_shift);
 
  if (!mat_inverse_matrix_d3(trans_mat_inv,trans_mat,0))
  {
    goto ret;
  }
 
 mat_multiply_matrix_d3(conv_prim->lattice,primitive->cell->lattice,trans_mat_inv);

  for (i = 0; i < primitive->cell->size; i++) {
    conv_prim->types[i] = primitive->cell->types[i];
    mat_copy_vector_d3(conv_prim->magmom[i],primitive->cell->magmom[i]);
    mat_multiply_matrix_vector_d3(conv_prim->position[i],
                                  trans_mat,
                                  primitive->cell->position[i]);
    for (j = 0; j < 3; j++) {
      conv_prim->position[i][j] += origin_shift_prim[j];

      conv_prim->position[i][j] = mat_Dmod1(conv_prim->position[i][j]);
    }
  }


  if ((primitive_new->orig_lattice =
       (double (*)[3]) malloc(sizeof(double[3]) * 3)) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    return NULL;
  }

  primitive_new->size = primitive->size;
  primitive_new->angle_tolerance = primitive->angle_tolerance;
  primitive_new->cell = conv_prim;

  mat_copy_matrix_d3(primitive_new->orig_lattice,primitive->orig_lattice);

  for (i=0;i<primitive->size;i++)
  {
    primitive_new->mapping_table[i] = primitive->mapping_table[i];
  }

  ret:
    return primitive_new;
}


/* Return NULL if failed */
Magnetic_Symmetry * sym_get_magnetic_operation(const Magnetic_Cell * primitive,
                                               const double symprec,
                                               const double angle_tolerance)
{
  Symmetry *symmetry;
  Magnetic_Symmetry *symmetry_new;
  Cell *cell;
  
  cell = NULL;
  symmetry = NULL;
  symmetry_new = NULL;

  cell = Magnetic_Cell_to_Cell(primitive);
  symmetry = sym_get_operation(cell,symprec,angle_tolerance);

  symmetry_new = get_magnetic_operations(symmetry, primitive, symprec);

  sym_free_symmetry(symmetry);
  cel_free_cell(cell);
  symmetry = NULL;
  cell = NULL;

  return symmetry_new;
}

static Magnetic_Symmetry * get_magnetic_operations(const Symmetry *sym_nonspin,
                                                   const Magnetic_Cell *cell,
                                                   const double symprec)
{
  Magnetic_Symmetry *symmetry;
  int i, j, k, T, num_sym, is_found;
  int *time_reversal;
  double pos[3];
  MatINT * rotations;
  VecDBL * trans;

  /* double operations with TR = 1 and TR = -1 */
  if ((time_reversal = (int *) malloc(sizeof(int)*sym_nonspin->size*2))==NULL)
  {
    return NULL;
  }

  rotations = mat_alloc_MatINT(sym_nonspin->size*2);
  trans = mat_alloc_VecDBL(sym_nonspin->size*2);
  num_sym = 0;
  is_found = 0;
  
  for (T = 1;T > -2;T = T-2)
  {
    for (i = 0; i < sym_nonspin->size; i++) {   
      for (j = 0; j < cell->size; j++) {
        mat_multiply_matrix_vector_id3(pos, sym_nonspin->rot[i], cell->position[j]);
        for (k = 0; k < 3; k++) {
          pos[k] += sym_nonspin->trans[i][k];
        }
        
        for (k = 0; k < cell->size; k++) {
          if (cel_is_overlap_with_same_type(cell->position[k],
                                            pos,
                                            cell->types[k],
                                            cell->types[j],
                                            cell->lattice,
                                            symprec)) {

            
            is_found = check_vector(j, k, cell->magmom, sym_nonspin->rot[i],
                                      cell->lattice, T, symprec);
            /* Break because cel_is_overlap_with_same_type == true */
            /* for only one atom. */
            break;
          }
        }
        if (! is_found) {
          break;
        }
      }
      if (is_found) {
        mat_copy_matrix_i3(rotations->mat[num_sym], sym_nonspin->rot[i]);
        mat_copy_vector_d3(trans->vec[num_sym], sym_nonspin->trans[i]);
        time_reversal[num_sym]=T;
        num_sym++;
      }
    }
  }



  symmetry = sym_alloc_magnetic_symmetry(num_sym);

  for (i = 0; i < num_sym; i++) {
    mat_copy_matrix_i3(symmetry->rot[i], rotations->mat[i]);
    mat_copy_vector_d3(symmetry->trans[i], trans->vec[i]);
    symmetry->time_reversal[i] = time_reversal[i];
  }
  

  free(time_reversal);

  mat_free_MatINT(rotations);
  rotations = NULL;
  mat_free_VecDBL(trans);
  trans = NULL;
  return symmetry;
}

/* Work in Cartesian coordinates. */
static int check_vector(const int j,
                        const int k,
                        const double vectors[][3],
                        SPGCONST int rot[3][3],
                        SPGCONST double lattice[3][3],
                        const int T,
                        const double symprec)
{
  int i, detR;
  double vec_jp[3], diff[3];
  double inv_lat[3][3], rot_cart[3][3];


  /* is_magnetic: Non-collinear magnetic moment m' = |detR|Rm */
  /* !is_magnetic: Usual vector: v' = Rv */

  mat_inverse_matrix_d3(inv_lat, lattice, 0);
  mat_multiply_matrix_id3(rot_cart, rot, inv_lat);
  mat_multiply_matrix_d3(rot_cart, lattice, rot_cart);


  /* v_j' = R v_j */
  mat_multiply_matrix_vector_d3(vec_jp, rot_cart, vectors[j]);

  detR = mat_get_determinant_i3(rot);
  for (i = 0; i < 3; i++) {
    diff[i] = mat_Dabs(T*detR * vec_jp[i] - vectors[k][i]);
  }

  if (diff[0] < symprec && diff[1] < symprec && diff[2] < symprec) {
    return 1;
  } else {
    return 0;
  }
}

int magnetic_type_check( const Magnetic_Symmetry * symmetry, const double symprec)
{
  /* only have T=1 type 1 
   (T=1) == (T=-1) type 1 
   pure translation operation with T = -1 type 4 
   else type 3 */


  int i,j;
  int rot_eq_num,trans_eq_num;
  double (*D_trans)[3],(*DT_trans)[3];
  int (*D_rot)[3][3],(*DT_rot)[3][3];

  if (symmetry->time_reversal[symmetry->size-1]  == 1)
  {
    return 1;
  }else
  {
    if ((symmetry->size % 2 ==1) || symmetry->time_reversal[0]==symmetry->time_reversal[symmetry->size/2])
    {
      return 0;
    }
    


    if ((D_trans = (double (*)[3]) malloc(sizeof(double[3])*symmetry->size/2 )) == NULL)
      {
        return 0;
      }
    if ((DT_trans = (double (*)[3]) malloc(sizeof(double[3])*symmetry->size/2 )) == NULL)
      {
        return 0;
      }
    if ((D_rot = (int (*)[3][3]) malloc(sizeof(int[3][3])*symmetry->size/2 )) == NULL)
      {
        return 0;
      }
      }
    if ((DT_rot = (int (*)[3][3]) malloc(sizeof(int[3][3])*symmetry->size/2 )) == NULL)
      {
        return 0;
      }


    for (i=0;i<symmetry->size/2;i++)
    {
      mat_copy_matrix_i3(D_rot[i],symmetry->rot[i]);
      mat_copy_matrix_i3(DT_rot[i],symmetry->rot[i+symmetry->size/2]);
      mat_copy_vector_d3(D_trans[i],symmetry->trans[i]);
      mat_copy_vector_d3(DT_trans[i],symmetry->trans[i+symmetry->size/2]);
    }

    /* compare T = 1 and T =-1 */

    rot_eq_num = 0;
    trans_eq_num = 0;

    for (i=0;i<symmetry->size/2;i++)
    {
      for (j=0;j<symmetry->size/2;j++)
      {
        if (mat_check_identity_matrix_i3(D_rot[i],DT_rot[j]))
        {
          rot_eq_num++;
          if (mat_check_identity_vector_d3(D_trans[i],DT_trans[j],symprec))
          {
            trans_eq_num++;
          }
        }
      }
    }

    free(D_trans);
    free(D_rot);
    free(DT_trans);
    free(DT_rot);

    if (rot_eq_num==symmetry->size/2 && trans_eq_num==symmetry->size/2 )
    {
      return 2;
    }else if (rot_eq_num==symmetry->size/2 && trans_eq_num!=symmetry->size/2 )
    {
      return 4;
    }else if (rot_eq_num > symmetry->size/2 && trans_eq_num > symmetry->size/2)
    {
      return 0;
    }else{
      return 3;
    }

    return 0;

  }

static Symmetry * nonmagnetic_symmetries_for_standard_cell(const Magnetic_Symmetry * symmetry, const int msg_type)
{
  Symmetry * symmetry_new;
  int i;
  
  symmetry_new = NULL;
  if (msg_type == 1 || msg_type == 3 )
  {
    if ((symmetry_new =
       magnetic_symmetry_to_symmetry(symmetry)) == NULL) {
    return NULL;
    }

  }else if (msg_type == 4 || msg_type ==2)
  {
    if ((symmetry_new = sym_alloc_symmetry(symmetry->size/2))==NULL){
      return NULL;
    }

    for (i = 0; i < symmetry->size/2;i++ )
    {
      mat_copy_matrix_i3(symmetry_new->rot[i],symmetry->rot[i]);
      mat_copy_vector_d3(symmetry_new->trans[i],symmetry->trans[i]);
    }



  }

  return symmetry_new;
}

static int mat_check_identity_vector_d3(const double a[3], const double b[3], const double symprec)
{
  if ( mat_Dabs( a[0] - b[0] ) > symprec ||
       mat_Dabs( a[1] - b[1] ) > symprec ||
       mat_Dabs( a[2] - b[2] ) > symprec )  
  {
    return 0;
  }
  else {
    return 1;
  }
}
/* Return NULL if failed */
static Spacegroup * search_spacegroup_with_symmetry(const Primitive * primitive,
                                                    const int candidates[],
                                                    const int num_candidates,
                                                    const Symmetry *symmetry,
                                                    const double symprec,
                                                    const double angle_tolerance)
{ 
  int hall_number;
  double conv_lattice[3][3];
  double origin_shift[3];
  Spacegroup *spacegroup;
  PointSymmetry pointsym;

  debug_print("search_spacegroup (tolerance = %f):\n", symprec);

  spacegroup = NULL;
  origin_shift[0] = 0;
  origin_shift[1] = 0;
  origin_shift[2] = 0;

  pointsym = ptg_get_pointsymmetry(symmetry->rot, symmetry->size);
  if (pointsym.size < symmetry->size) {
    warning_print("spglib: Point symmetry of primitive cell is broken. ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
    return NULL;
  }

  hall_number = iterative_search_hall_number(origin_shift,
                                             conv_lattice,
                                             candidates,
                                             num_candidates,
                                             primitive,
                                             symmetry,
                                             symprec,
                                             angle_tolerance);
  if (hall_number == 0) {
    return NULL;
  }



  spacegroup = get_spacegroup(hall_number, origin_shift, conv_lattice);

  return spacegroup;
}

/* modified part end*/

/* Return spacegroup.number = 0 if failed */
static Spacegroup * get_spacegroup(const int hall_number,
                                   const double origin_shift[3],
                                   SPGCONST double conv_lattice[3][3])
{
  Spacegroup *spacegroup;
  SpacegroupType spacegroup_type;

  spacegroup = NULL;

  if ((spacegroup = (Spacegroup*) malloc(sizeof(Spacegroup))) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    return NULL;
  }

  if (0 < hall_number && hall_number < 531) {
    spacegroup_type = spgdb_get_spacegroup_type(hall_number);
    mat_copy_matrix_d3(spacegroup->bravais_lattice, conv_lattice);
    mat_copy_vector_d3(spacegroup->origin_shift, origin_shift);
    spacegroup->number = spacegroup_type.number;
    spacegroup->hall_number = hall_number;
    spacegroup->pointgroup_number = spacegroup_type.pointgroup_number;
    memcpy(spacegroup->schoenflies, spacegroup_type.schoenflies, 7);
    memcpy(spacegroup->hall_symbol, spacegroup_type.hall_symbol, 17);
    memcpy(spacegroup->international, spacegroup_type.international, 32);
    memcpy(spacegroup->international_long, spacegroup_type.international_full, 20);
    memcpy(spacegroup->international_short,
           spacegroup_type.international_short, 11);
    memcpy(spacegroup->choice, spacegroup_type.choice, 6);
  }

  return spacegroup;
}

/* Return 0 if failed */
static int iterative_search_hall_number(double origin_shift[3],
                                        double conv_lattice[3][3],
                                        const int candidates[],
                                        const int num_candidates,
                                        const Primitive * primitive,
                                        const Symmetry * symmetry,
                                        const double symprec,
                                        const double angle_tolerance)
{
  int attempt, hall_number;
  double tolerance;
  Symmetry * sym_reduced;

  debug_print("iterative_search_hall_number:\n");

  hall_number = 0;
  sym_reduced = NULL;

  hall_number = search_hall_number(origin_shift,
                                   conv_lattice,
                                   candidates,
                                   num_candidates,
                                   primitive,
                                   symmetry,
                                   symprec);

  if (hall_number > 0) {
    goto ret;
  }

  tolerance = symprec;
  for (attempt = 0; attempt < NUM_ATTEMPT; attempt++) {

    warning_print("spglib: Attempt %d tolerance = %e failed",
                  attempt, tolerance);
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);

    tolerance *= REDUCE_RATE;
    /* need modified for magnetic */
    sym_reduced = sym_reduce_operation(primitive->cell,
                                       symmetry,
                                       tolerance,
                                       angle_tolerance);
    hall_number = search_hall_number(origin_shift,
                                     conv_lattice,
                                     candidates,
                                     num_candidates,
                                     primitive,
                                     sym_reduced,
                                     symprec);
    /* need modified for mgnetic */
    sym_free_symmetry(sym_reduced);
    sym_reduced = NULL;
    if (hall_number > 0) {
      break;
    }
  }

 ret:
  return hall_number;
}

/* Return 0 if failed */
static int search_hall_number(double origin_shift[3],
                              double conv_lattice[3][3],
                              const int candidates[],
                              const int num_candidates,
                              const Primitive * primitive,
                              const Symmetry * symmetry,
                              const double symprec)
{
  int i, hall_number;
  Centering centering;
  Pointgroup pointgroup;
  Symmetry * conv_symmetry;
  int tmat_int[3][3];
  double correction_mat[3][3], tmat[3][3], conv_lattice_tmp[3][3];

  debug_print("search_hall_number:\n");

  hall_number = 0;
  conv_symmetry = NULL;

  /* primitive->cell->lattice is set right handed. */
  /* tmat_int is set right handed. */

  /* For rhombohedral system, basis vectors for hP */
  /* (hexagonal lattice basis) are provided by tmat_int, */
  /* but may be either obverse or reverse setting. */

  /*no need to information about structure it is for initial transformation matrix*/
  pointgroup = ptg_get_transformation_matrix(tmat_int,
                                             symmetry->rot,
                                             symmetry->size);

  debug_print("[line %d, %s]\n", __LINE__, __FILE__);
  debug_print("initial tranformation matrix\n");
  debug_print_matrix_i3(tmat_int);

  if (pointgroup.number == 0) {
    goto err;
  }

  /* For LAUE1 and LAUE2M, update tmat_int by making smallest lattice */
  if (pointgroup.laue == LAUE1 || pointgroup.laue == LAUE2M) {
    mat_multiply_matrix_di3(
      conv_lattice_tmp, primitive->cell->lattice, tmat_int);

    if (pointgroup.laue == LAUE1) {
      if (! change_basis_tricli(tmat_int,
                                conv_lattice_tmp,
                                primitive->cell->lattice,
                                symprec)) {
        goto err;
      }
    }

    /* The angle between non-unique axes may be acute or non-acute. */
    if (pointgroup.laue == LAUE2M) {
      if (! change_basis_monocli(tmat_int,
                                 conv_lattice_tmp,
                                 primitive->cell->lattice,
                                 symprec)) {
        goto err;
      }
    }
  }

  /* For rhombohedral system, correction matrix to hR */
  /* (a=b=c primitive lattice basis) is returned judging */
  /* either obverse or reverse setting. centering=R_CENTER. */

  /* seems has nothing to do with information about structure */
  if ((centering = get_centering(correction_mat,
                                 tmat_int,
                                 pointgroup.laue)) == CENTERING_ERROR) {
    goto err;
  }

  mat_multiply_matrix_id3(tmat, tmat_int, correction_mat);
  mat_multiply_matrix_d3(conv_lattice, primitive->cell->lattice, tmat);

  debug_print("[line %d, %s]\n", __LINE__, __FILE__);
  debug_print("tranformation matrix\n");
  debug_print_matrix_d3(tmat);

  /* For rhombohedral system, symmetry for a=b=c primitive lattice */
  /* basis is returned although centering == R_CENTER. */

  /* seems transform form of symmetry no need information about structure*/
  if ((conv_symmetry = get_initial_conventional_symmetry(centering,
                                                         tmat,
                                                         symmetry)) == NULL) {
    goto err;
  }

  for (i = 0; i < num_candidates; i++) {
    /* Check if hall_number is that of rhombohedral. */
    if (match_hall_symbol_db(origin_shift,
                             conv_lattice, /* <-- modified only matched */
                             primitive->orig_lattice,
                             candidates[i],
                             pointgroup.number,
                             pointgroup.holohedry,
                             centering,
                             conv_symmetry,
                             symprec)) {

      debug_print("[line %d, %s]\n", __LINE__, __FILE__);
      debug_print("origin shift\n");
      debug_print_vector_d3(origin_shift);

      hall_number = candidates[i];
      break;
    }
  }

  sym_free_symmetry(conv_symmetry);
  conv_symmetry = NULL;

  return hall_number;

 err:
  return 0;
}

/* Triclinic: Niggli cell reduction */
/* Return 0 if failed */
static int change_basis_tricli(int tmat_int[3][3],
                               SPGCONST double conv_lattice[3][3],
                               SPGCONST double primitive_lattice[3][3],
                               const double symprec)
{
  int i, j;
  double niggli_cell[9];
  double smallest_lattice[3][3], inv_lattice[3][3], tmat[3][3];

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      niggli_cell[i * 3 + j] = conv_lattice[i][j];
    }
  }

  if (! niggli_reduce(niggli_cell, symprec * symprec)) {
    return 0;
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      smallest_lattice[i][j] = niggli_cell[i * 3 + j];
    }
  }
  if (mat_get_determinant_d3(smallest_lattice) < 0) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        smallest_lattice[i][j] = -smallest_lattice[i][j];
      }
    }
  }
  mat_inverse_matrix_d3(inv_lattice, primitive_lattice, 0);
  mat_multiply_matrix_d3(tmat, inv_lattice, smallest_lattice);
  mat_cast_matrix_3d_to_3i(tmat_int, tmat);

  return 1;
}

/* Monoclinic: choose shortest a, c lattice vectors (|a| < |c|) */
/* Return 0 if failed */
static int change_basis_monocli(int tmat_int[3][3],
                                SPGCONST double conv_lattice[3][3],
                                SPGCONST double primitive_lattice[3][3],
                                const double symprec)
{
  double smallest_lattice[3][3], inv_lattice[3][3], tmat[3][3];

  if (! del_delaunay_reduce_2D(smallest_lattice,
                               conv_lattice,
                               1, /* unique axis of b */
                               symprec)) {
    return 0;
  }

  mat_inverse_matrix_d3(inv_lattice, primitive_lattice, 0);
  mat_multiply_matrix_d3(tmat, inv_lattice, smallest_lattice);
  mat_cast_matrix_3d_to_3i(tmat_int, tmat);
  return 1;
}

/* Return NULL if failed */

/*depend on centering type*/
static Symmetry *
get_initial_conventional_symmetry(const Centering centering,
                                  SPGCONST double tmat[3][3],
                                  const Symmetry * symmetry)
{
  Symmetry * conv_symmetry;

  debug_print("get_initial_conventional_symmetry\n");

  conv_symmetry = NULL;

  if (centering == R_CENTER) {
    /* hP for rhombohedral */

    /*what is different from get_initial_conventional_symmetry*/
    conv_symmetry = get_conventional_symmetry(tmat,
                                              PRIMITIVE,
                                              symmetry);
  } else {
    conv_symmetry = get_conventional_symmetry(tmat,
                                              centering,
                                              symmetry);
  }

  return conv_symmetry;
}

/* Return 0 if failed */
static int match_hall_symbol_db(double origin_shift[3],
                                double conv_lattice[3][3],
                                SPGCONST double (*orig_lattice)[3],
                                const int hall_number,
                                const int pointgroup_number,
                                const Holohedry holohedry,
                                const Centering centering,
                                const Symmetry *symmetry,
                                const double symprec)
{
  int is_found;
  SpacegroupType spacegroup_type;
  Symmetry * changed_symmetry;
  double changed_lattice[3][3], inv_lattice[3][3], tmat[3][3];

  changed_symmetry = NULL;

  spacegroup_type = spgdb_get_spacegroup_type(hall_number);

  if (pointgroup_number != spacegroup_type.pointgroup_number) {
    goto err;
  }

  switch (holohedry) {
  case MONOCLI:
    if (match_hall_symbol_db_monocli(origin_shift,
                                     conv_lattice,
                                     orig_lattice,
                                     hall_number,
                                     spacegroup_type.number,
                                     centering,
                                     symmetry,
                                     symprec)) {return 1;}
    break;

  case ORTHO:
    if (num_axis_choices_ortho[spacegroup_type.number - 16] == 6) {
      if (match_hall_symbol_db_ortho(origin_shift,
                                     conv_lattice,
                                     orig_lattice,
                                     hall_number,
                                     centering,
                                     symmetry,
                                     6,
                                     symprec)) {return 1;}
      break;
    }

    if (num_axis_choices_ortho[spacegroup_type.number - 16] == 3) {
      if (match_hall_symbol_db_ortho(origin_shift,
                                     conv_lattice,
                                     orig_lattice,
                                     hall_number,
                                     centering,
                                     symmetry,
                                     3,
                                     symprec)) {return 1;}
      break;
    }

    /* Switching two axes */
    /* Two steps: */
    /*   1. Finding principal axis for the representative hall symbol */
    /*      of the specified hall symbol without checking basis vector */
    /*      lengths preference. */
    /*   2. Finding transformation matrix and origin shift for the */
    /*      specified hall symbol. */
    if (num_axis_choices_ortho[spacegroup_type.number - 16] == 2) {
      mat_copy_matrix_d3(changed_lattice, conv_lattice);
      if (! match_hall_symbol_db_ortho(
            origin_shift,
            changed_lattice,
            orig_lattice,
            spacegroup_to_hall_number[spacegroup_type.number - 1],
            centering,
            symmetry,
            0,
            symprec)) {break;}
      mat_inverse_matrix_d3(inv_lattice, conv_lattice, 0);
      mat_multiply_matrix_d3(tmat, inv_lattice, changed_lattice);

      if ((changed_symmetry = get_conventional_symmetry(tmat,
                                                        PRIMITIVE,
                                                        symmetry)) == NULL) {
        goto err;
      }

      is_found = match_hall_symbol_db_ortho(origin_shift,
                                            changed_lattice,
                                            orig_lattice,
                                            hall_number,
                                            centering,
                                            changed_symmetry,
                                            2,
                                            symprec);
      sym_free_symmetry(changed_symmetry);
      changed_symmetry = NULL;
      if (is_found) {
        mat_copy_matrix_d3(conv_lattice, changed_lattice);
        return 1;
      }
      break;
    }

    if (num_axis_choices_ortho[spacegroup_type.number - 16] == 1) {
      if (match_hall_symbol_db_ortho(origin_shift,
                                     conv_lattice,
                                     orig_lattice,
                                     hall_number,
                                     centering,
                                     symmetry,
                                     1,
                                     symprec)) {return 1;}
      break;
    }

    break;

  case CUBIC:
    if (match_hall_symbol_db_cubic(origin_shift,
                                   conv_lattice,
                                   orig_lattice,
                                   hall_number,
                                   centering,
                                   symmetry,
                                   symprec)) {return 1;}
    break;

  case TRIGO:
    if ((centering == R_CENTER) && (
          hall_number == 433 ||
          hall_number == 434 ||
          hall_number == 436 ||
          hall_number == 437 ||
          hall_number == 444 ||
          hall_number == 445 ||
          hall_number == 450 ||
          hall_number == 451 ||
          hall_number == 452 ||
          hall_number == 453 ||
          hall_number == 458 ||
          hall_number == 459 ||
          hall_number == 460 ||
          hall_number == 461)) {
      /* Rhombohedral. symmetry is for a=b=c basis. */
      if (match_hall_symbol_db_rhombo(origin_shift,
                                      conv_lattice,
                                      orig_lattice,
                                      hall_number,
                                      symmetry,
                                      symprec)) {return 1;}
      break;
    }
    /* Do not break for other trigonal cases */
  default: /* HEXA, TETRA, TRICLI and rest of TRIGO */
    if (match_hall_symbol_db_others(origin_shift,
                                    conv_lattice,
                                    orig_lattice,
                                    hall_number,
                                    centering,
                                    holohedry,
                                    symmetry,
                                    symprec)) {return 1;}
    break;
  }

 err:
  return 0;
}

/* Return 0 if failed */
static int match_hall_symbol_db_monocli(double origin_shift[3],
                                        double conv_lattice[3][3],
                                        SPGCONST double (*orig_lattice)[3],
                                        const int hall_number,
                                        const int space_group_number,
                                        const Centering centering,
                                        const Symmetry *conv_symmetry,
                                        const double symprec)
{
  int i, check_norms, i_shortest, is_found_any;
  int is_found[36];
  double shortest_norm_sum, norm_sum;
  double norms_squared[36][2];
  double all_conv_lattices[36][3][3], all_origin_shifts[36][3];

  /* E. Parthe and L. M. Gelato */
  /* "The best unit cell for monoclinic structure..." (1983) */
  if (space_group_number == 3 ||
      space_group_number == 4 ||
      space_group_number == 6 ||
      space_group_number == 10 ||
      space_group_number == 11) {
    /* |a| < |c| for unique axis b. (This is as written in the paper 1983.) */
    /* |a| < |b| for unique axis c. (This is spgilb definition.) */
    /* |b| < |c| for unique axis a. (This is spgilb definition.) */
    check_norms = 1;
  } else {
    check_norms = 0;
  }

  for (i = 0; i < 36; i++) {
    mat_copy_matrix_d3(all_conv_lattices[i], conv_lattice);
    /* all_origin_shifts[i] is possibly overwritten. */
    /* is_found == 0: Not found */
    /* is_found == 1: Found. conv_lattice may not be similar to orig_lattice */
    /* is_found == 2: Found. conv_lattice is similar to orig_lattice */
    is_found[i] = match_hall_symbol_db_monocli_in_loop(all_origin_shifts[i],
                                                       all_conv_lattices[i],
                                                       norms_squared[i],
                                                       i,
                                                       orig_lattice,
                                                       check_norms,
                                                       hall_number,
                                                       centering,
                                                       conv_symmetry,
                                                       symprec);
  }

  is_found_any = 0;
  for (i = 0; i < 36; i++) {
    if (is_found[i]) {
      i_shortest = i;
      shortest_norm_sum = sqrt(norms_squared[i][0]) + sqrt(norms_squared[i][1]);
      is_found_any = 1;
      break;
    }
  }

  if (! is_found_any) {
    goto err;
  }

  /* Find shortest non-unique axes lengths */
  for (i = 0; i < 36; i++) {
    if (is_found[i]) {
      norm_sum = sqrt(norms_squared[i][0]) + sqrt(norms_squared[i][1]);
      if (shortest_norm_sum > norm_sum) {
        shortest_norm_sum = norm_sum;
      }
    }
  }

  /* Prefers is_found == 2, i.e., similar to orig_lattice */
  i_shortest = -1;
  for (i = 0; i < 36; i++) {
    if (is_found[i]) {
      norm_sum = sqrt(norms_squared[i][0]) + sqrt(norms_squared[i][1]);
      if (mat_Dabs(norm_sum - shortest_norm_sum) < symprec) {
        if (is_found[i] == 2) {
          i_shortest = i;
          break;
        }
        if (i_shortest < 0) {
          i_shortest = i;
        }
      }
    }
  }

  mat_copy_vector_d3(origin_shift, all_origin_shifts[i_shortest]);
  mat_copy_matrix_d3(conv_lattice, all_conv_lattices[i_shortest]);
  return 1;

 err:
  return 0;
}

static int match_hall_symbol_db_monocli_in_loop(double origin_shift[3],
                                                double conv_lattice[3][3],
                                                double norms_squared[2],
                                                const int i,
                                                SPGCONST double (*orig_lattice)[3],
                                                const int check_norms,
                                                const int hall_number,
                                                const Centering centering,
                                                const Symmetry *conv_symmetry,
                                                const double symprec)
{
  int j, k, l, is_found, retval, unique_axis;
  double vec[2][3];
  double tmat[3][3], change_of_basis[3][3];
  Centering changed_centering;
  Symmetry * changed_symmetry;

  /* centring type should be P or C */
  if (centering == C_FACE) {
    changed_centering = change_of_centering_monocli[i];
  } else { /* suppose PRIMITIVE */
    changed_centering = centering;
  }

  mat_copy_matrix_d3(change_of_basis, change_of_basis_monocli[i]);
  mat_multiply_matrix_d3(conv_lattice, conv_lattice, change_of_basis);
  unique_axis = change_of_unique_axis_monocli[i];

  /* Make non-acute and length preference */
  /* norms_squared[0] and norms_squared[1] are the norms of the two */
  /* non unique axes among a,b,c. */
  l = 0;
  for (j = 0; j < 3; j++) {
    if (j == unique_axis) {continue;}
    for (k = 0; k < 3; k++) {
      vec[l][k] = conv_lattice[k][j];
    }
    norms_squared[l] = mat_norm_squared_d3(vec[l]);
    l++;
  }

  /* discard if principal angle is acute. */
  if (vec[0][0] * vec[1][0] +
      vec[0][1] * vec[1][1] +
      vec[0][2] * vec[1][2] > ZERO_PREC) {
    goto cont;
  }

  /* Choose |a| < |b| < |c| among two non-principles axes */
  /* if there are freedom. */
  if (check_norms) {
    if (norms_squared[0] > norms_squared[1] + ZERO_PREC) {
      goto cont;
    }
  }

  /* When orig_lattice is given not NULL, try to find similar (a, b, c) */
  /* choises to the input (a, b, c) by flipping a, b, c axes. */
  /* Here flipping means a -> -a, and so on. */
  /* Note that flipped (a,b,c) that match those of input should not */
  /* change centring for monoclinic case. */
  retval = 1;
  if (orig_lattice != NULL) {
    if (mat_get_determinant_d3(orig_lattice) > symprec) {
      /* This (mode=1) effectively checks C2 rotation along unique axis. */
      if (is_equivalent_lattice(
            tmat, 1, conv_lattice, orig_lattice, symprec)) {
        if (tmat[(unique_axis + 1) % 3][(unique_axis + 1) % 3]
            * tmat[(unique_axis + 2) % 3][(unique_axis + 2) % 3] > ZERO_PREC) {
          mat_multiply_matrix_d3(conv_lattice, conv_lattice, tmat);
          mat_multiply_matrix_d3(change_of_basis, change_of_basis, tmat);
          retval = 2;
        }
      }
    }
  }

  if ((changed_symmetry = get_conventional_symmetry(change_of_basis,
                                                    PRIMITIVE,
                                                    conv_symmetry)) == NULL) {
    goto cont;
  }

  is_found = hal_match_hall_symbol_db(origin_shift,
                                      conv_lattice,
                                      hall_number,
                                      changed_centering,
                                      changed_symmetry,
                                      symprec);
  sym_free_symmetry(changed_symmetry);
  changed_symmetry = NULL;

  if (is_found) {
    return retval;
  }

cont:
  return 0;
}

/* Return 0 if failed */
static int match_hall_symbol_db_ortho(double origin_shift[3],
                                      double conv_lattice[3][3],
                                      SPGCONST double (*orig_lattice)[3],
                                      const int hall_number,
                                      const Centering centering,
                                      const Symmetry *conv_symmetry,
                                      const int num_free_axes,
                                      const double symprec)
{
  int i;

  /* Try to find the best (a, b, c) by flipping axes. */
  if (orig_lattice != NULL) {
    if (mat_get_determinant_d3(orig_lattice) > symprec) {
      for (i = 0; i < 6; i++) {
        if (match_hall_symbol_db_ortho_in_loop(origin_shift,
                                               conv_lattice,
                                               orig_lattice,
                                               i,
                                               hall_number,
                                               centering,
                                               conv_symmetry,
                                               num_free_axes,
                                               symprec)) {
          return 1;
        }
      }
    }
  }

  /* If flipping didn't work, usual search is made. */
  for (i = 0; i < 6; i++) {
    if (match_hall_symbol_db_ortho_in_loop(origin_shift,
                                           conv_lattice,
                                           NULL,
                                           i,
                                           hall_number,
                                           centering,
                                           conv_symmetry,
                                           num_free_axes,
                                           symprec)) {
      return 1;
    }
  }

  return 0;
}

static int match_hall_symbol_db_ortho_in_loop(double origin_shift[3],
                                              double conv_lattice[3][3],
                                              SPGCONST double (*orig_lattice)[3],
                                              const int i,
                                              const int hall_number,
                                              const Centering centering,
                                              const Symmetry *symmetry,
                                              const int num_free_axes,
                                              const double symprec)
{
  int j, k, l, is_found;
  double vec[3], norms[3];
  Centering changed_centering;
  Symmetry * changed_symmetry;
  double changed_lattice[3][3], tmat[3][3], change_of_basis[3][3];

  changed_symmetry = NULL;

  if (centering == C_FACE) {
    changed_centering = change_of_centering_ortho[i];
  } else {
    changed_centering = centering;
  }

  mat_multiply_matrix_d3(changed_lattice,
                         conv_lattice,
                         change_of_basis_ortho[i]);

  /* When orig_lattice is given not NULL, try to find similar (a, b, c) */
  /* choises to the input (a, b, c) by flipping a, b, c axes. */
  /* Here flipping means a -> -a, and so on. */
  /* Note that flipping of axes doesn't change centring. */
  mat_copy_matrix_d3(change_of_basis, change_of_basis_ortho[i]);
  if (orig_lattice != NULL) {
    if (is_equivalent_lattice(
          tmat, 1, changed_lattice, orig_lattice, symprec)) {
      mat_multiply_matrix_d3(changed_lattice, changed_lattice, tmat);
      mat_multiply_matrix_d3(change_of_basis, change_of_basis, tmat);
    } else {
      goto cont;  /* This is necessary to run through all */
                  /* change_of_basis_ortho. */
    }
  }

  if (num_free_axes == 2) {
    l = 0;
    for (j = 0; j < 3; j++) {
      if (j == change_of_unique_axis_ortho[i]) {continue;}
      for (k = 0; k < 3; k++) {vec[k] = changed_lattice[k][j];}
      norms[l] = mat_norm_squared_d3(vec);
      l++;
    }
    if (norms[0] > norms[1] + ZERO_PREC) {goto cont;}
  }

  if (num_free_axes == 3) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {vec[k] = changed_lattice[k][j];}
      norms[j] = mat_norm_squared_d3(vec);
    }
    if ((norms[0] > norms[1] + ZERO_PREC) ||
        (norms[0] > norms[2] + ZERO_PREC)) {goto cont;}
  }

  if (num_free_axes == 6) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {vec[k] = changed_lattice[k][j];}
      norms[j] = mat_norm_squared_d3(vec);
    }
    if ((norms[0] > norms[1] + ZERO_PREC) ||
        (norms[1] > norms[2] + ZERO_PREC)) {goto cont;}
  }

  if ((changed_symmetry = get_conventional_symmetry(change_of_basis,
                                                    PRIMITIVE,
                                                    symmetry)) == NULL) {
    goto cont;
  }

  is_found = hal_match_hall_symbol_db(origin_shift,
                                      changed_lattice,
                                      hall_number,
                                      changed_centering,
                                      changed_symmetry,
                                      symprec);
  sym_free_symmetry(changed_symmetry);
  changed_symmetry = NULL;
  if (is_found) {
    mat_copy_matrix_d3(conv_lattice, changed_lattice);
    return 1;
  }

 cont:
  return 0;
}

static int match_hall_symbol_db_cubic(double origin_shift[3],
                                      double conv_lattice[3][3],
                                      SPGCONST double (*orig_lattice)[3],
                                      const int hall_number,
                                      const Centering centering,
                                      const Symmetry *conv_symmetry,
                                      const double symprec)
{
  int i;

  /* Special treatment for No. 205 (501) is included in this change of */
  /* basis. To see old code, commit hash 67997654 and change_of_basis_501. */
  if (orig_lattice != NULL) {
    if (mat_get_determinant_d3(orig_lattice) > symprec) {
      for (i = 0; i < 6; i++) {
        if (match_hall_symbol_db_cubic_in_loop(origin_shift,
                                               conv_lattice,
                                               orig_lattice,
                                               i,
                                               hall_number,
                                               centering,
                                               conv_symmetry,
                                               symprec)) {
          return 1;
        }
      }
    }
  }

  for (i = 0; i < 6; i++) {
    if (match_hall_symbol_db_cubic_in_loop(origin_shift,
                                           conv_lattice,
                                           NULL,
                                           i,
                                           hall_number,
                                           centering,
                                           conv_symmetry,
                                           symprec)) {
      return 1;
    }
  }
  return 0;
}

static int match_hall_symbol_db_cubic_in_loop(double origin_shift[3],
                                              double conv_lattice[3][3],
                                              SPGCONST double (*orig_lattice)[3],
                                              const int i,
                                              const int hall_number,
                                              const Centering centering,
                                              const Symmetry *conv_symmetry,
                                              const double symprec)
{
  int is_found;
  double changed_lattice[3][3], tmat[3][3], change_of_basis[3][3];
  Symmetry *changed_symmetry;

  changed_symmetry = NULL;

  mat_copy_matrix_d3(change_of_basis, change_of_basis_ortho[i]);
  mat_multiply_matrix_d3(changed_lattice, conv_lattice, change_of_basis);

  if (orig_lattice != NULL) {
    if (is_equivalent_lattice(
          tmat, 1, changed_lattice, orig_lattice, symprec)) {
      mat_multiply_matrix_d3(changed_lattice, changed_lattice, tmat);
      mat_multiply_matrix_d3(change_of_basis, change_of_basis, tmat);
    } else {
      goto cont;  /* This is necessary to run through all */
                  /* change_of_basis_ortho. */
    }
  }

  if ((changed_symmetry = get_conventional_symmetry(change_of_basis,
                                                    PRIMITIVE,
                                                    conv_symmetry)) == NULL) {
    goto cont;
  }

  is_found = hal_match_hall_symbol_db(origin_shift,
                                      changed_lattice,
                                      hall_number,
                                      centering,
                                      changed_symmetry,
                                      symprec);

  sym_free_symmetry(changed_symmetry);
  changed_symmetry = NULL;

  if (is_found) {
    mat_copy_matrix_d3(conv_lattice, changed_lattice);
    return 1;
  }

cont:
  return 0;
}

static int match_hall_symbol_db_rhombo(double origin_shift[3],
                                       double conv_lattice[3][3],
                                       SPGCONST double (*orig_lattice)[3],
                                       const int hall_number,
                                       const Symmetry *conv_symmetry,
                                       const double symprec)
{
  int is_found;

  if (hall_number == 433 ||
      hall_number == 436 ||
      hall_number == 444 ||
      hall_number == 450 ||
      hall_number == 452 ||
      hall_number == 458 ||
      hall_number == 460) {  /* Hexagonal lattice */
    is_found = match_hall_symbol_db_change_of_basis_loop(
      origin_shift,
      conv_lattice,
      orig_lattice,
      change_of_basis_rhombo_hex,
      6,
      hall_number,
      R_CENTER,
      conv_symmetry,
      symprec);
    if (is_found) {
      /* mat_copy_matrix_d3(conv_lattice, changed_lattice); */
      return 1;
    } else {
      return 0;
    }
  } else {  /* Rhombohedral (a=b=c) lattice */
    return match_hall_symbol_db_change_of_basis_loop(
      origin_shift,
      conv_lattice,
      orig_lattice,
      change_of_basis_rhombo,
      6,
      hall_number,
      PRIMITIVE,
      conv_symmetry,
      symprec);
  }

  return 0;
}

/* HEXA, TETRA, TRICLI and TRIGO but not Rhombohedral a=b=c */
static int match_hall_symbol_db_others(double origin_shift[3],
                                       double conv_lattice[3][3],
                                       SPGCONST double (*orig_lattice)[3],
                                       const int hall_number,
                                       const Centering centering,
                                       const Holohedry holohedry,
                                       const Symmetry *conv_symmetry,
                                       const double symprec)
{
  /* TRICLI: No check. */
  if (holohedry == TRICLI) {
    return hal_match_hall_symbol_db(origin_shift,
                                    conv_lattice,
                                    hall_number,
                                    centering,
                                    conv_symmetry,
                                    symprec);
  }

  /* TETRA: Check over 4 fold rotations with/without flipping c-axis */
  if (holohedry == TETRA) {
    return match_hall_symbol_db_change_of_basis_loop(
      origin_shift,
      conv_lattice,
      orig_lattice,
      change_of_basis_C4,
      8,
      hall_number,
      centering,
      conv_symmetry,
      symprec);
  }

  /* HEXA and TRIGO but not Rhombohedral: */
  /*    Check over 6 fold rotations with/without flipping c-axis */
  return match_hall_symbol_db_change_of_basis_loop(
    origin_shift,
    conv_lattice,
    orig_lattice,
    change_of_basis_C6,
    12,
    hall_number,
    centering,
    conv_symmetry,
    symprec);
}

static int match_hall_symbol_db_change_of_basis_loop(
  double origin_shift[3],
  double conv_lattice[3][3],
  SPGCONST double (*orig_lattice)[3],
  SPGCONST double (*change_of_basis)[3][3],
  const int num_change_of_basis,
  const int hall_number,
  const Centering centering,
  const Symmetry *conv_symmetry,
  const double symprec)
{
  int i, is_found;
  double changed_lattice[3][3], tmat[3][3];
  Symmetry *changed_symmetry;
  Centering centering_for_symmetry;

  changed_symmetry = NULL;

  if (centering == R_CENTER) {
    centering_for_symmetry = R_CENTER;
  } else {
    centering_for_symmetry = PRIMITIVE;
  }

  /* Check of similarity of basis vectors to those of input */
  if (orig_lattice != NULL) {
    if (mat_get_determinant_d3(orig_lattice) > symprec) {
      for (i = 0; i < num_change_of_basis; i++) {
        if ((changed_symmetry = get_conventional_symmetry(
               change_of_basis[i], centering_for_symmetry, conv_symmetry))
            == NULL) {
          continue;
        }
        mat_multiply_matrix_d3(
          changed_lattice, conv_lattice, change_of_basis[i]);
        if (is_equivalent_lattice(
              tmat, 0, changed_lattice, orig_lattice, symprec)) {
          /* Here R_CENTER means hP (hexagonal) setting. */
          is_found = hal_match_hall_symbol_db(origin_shift,
                                              changed_lattice,
                                              hall_number,
                                              centering,
                                              changed_symmetry,
                                              symprec);
        } else {
          is_found = 0;
        }
        sym_free_symmetry(changed_symmetry);
        changed_symmetry = NULL;
        if (is_found) {
          mat_copy_matrix_d3(conv_lattice, changed_lattice);
          return 1;
        }
      }
    }
  }

  /* No check of similarity of basis vectors to those of input */
  for (i = 0; i < num_change_of_basis; i++) {
    if ((changed_symmetry = get_conventional_symmetry(
           change_of_basis[i], centering_for_symmetry, conv_symmetry)) == NULL) {
      continue;
    }
    mat_multiply_matrix_d3(changed_lattice, conv_lattice, change_of_basis[i]);
    is_found = hal_match_hall_symbol_db(origin_shift,
                                        changed_lattice,
                                        hall_number,
                                        centering,
                                        changed_symmetry,
                                        symprec);
    sym_free_symmetry(changed_symmetry);
    changed_symmetry = NULL;
    if (is_found) {
      mat_copy_matrix_d3(conv_lattice, changed_lattice);
      return 1;
    }
  }

  return 0;
}

/* Return NULL if failed */
static Symmetry * get_conventional_symmetry(SPGCONST double tmat[3][3],
                                            const Centering centering,
                                            const Symmetry *primitive_sym)
{
  int i, j, k, multi, size;
  double inv_tmat[3][3], shift[3][3];
  double symmetry_rot_d3[3][3], primitive_sym_rot_d3[3][3];
  Symmetry *symmetry;

  symmetry = NULL;

  size = primitive_sym->size;

  switch (centering) {
  case FACE:
    if ((symmetry = sym_alloc_symmetry(size * 4)) == NULL) {
      return NULL;
    }
    break;
  case R_CENTER:
    if ((symmetry = sym_alloc_symmetry(size * 3)) == NULL) {
      return NULL;
    }
    break;
  case BODY:
  case A_FACE:
  case B_FACE:
  case C_FACE:
    if ((symmetry = sym_alloc_symmetry(size * 2)) == NULL) {
      return NULL;
    }
    break;
  default:
    if ((symmetry = sym_alloc_symmetry(size)) == NULL) {
      return NULL;
    }
    break;
  }

  for (i = 0; i < size; i++) {
    mat_cast_matrix_3i_to_3d(primitive_sym_rot_d3, primitive_sym->rot[i]);

    /* C*S*C^-1: recover conventional cell symmetry operation */
    mat_get_similar_matrix_d3(symmetry_rot_d3,
                              primitive_sym_rot_d3,
                              tmat,
                              0);
    mat_cast_matrix_3d_to_3i(symmetry->rot[i], symmetry_rot_d3);

    /* translation in conventional cell: C = B^-1*P */
    mat_inverse_matrix_d3(inv_tmat, tmat, 0);
    mat_multiply_matrix_vector_d3(symmetry->trans[i],
                                  inv_tmat,
                                  primitive_sym->trans[i]);
  }

  multi = 1;

  if (centering != PRIMITIVE) {
    multi = get_centering_shifts(shift, centering);
    for (i = 0; i < multi - 1; i++) {
      for (j = 0; j < size; j++) {
        mat_copy_matrix_i3(symmetry->rot[(i+1) * size + j],
                           symmetry->rot[j]);
        for (k = 0; k < 3; k++) {
          symmetry->trans[(i+1) * size + j][k] =
            symmetry->trans[j][k] + shift[i][k];
        }
      }
    }
  }

  for (i = 0; i < multi; i++) {
    for (j = 0; j < size; j++) {
      for (k = 0; k < 3; k++) {
        symmetry->trans[i * size + j][k] =
          mat_Dmod1(symmetry->trans[i * size + j][k]);
      }
    }
  }

  return symmetry;
}

/* Return CENTERING_ERROR if failed */
static Centering get_centering(double correction_mat[3][3],
                               SPGCONST int tmat[3][3],
                               const Laue laue)
{
  int det;
  double trans_corr_mat[3][3];
  Centering centering;

  mat_copy_matrix_d3(correction_mat, identity);
  det = abs(mat_get_determinant_i3(tmat));
  debug_print("laue class: %d\n", laue);
  debug_print("multiplicity: %d\n", det);

  switch (det) {

  case 1:
    centering = PRIMITIVE;
    break;

  case 2:
    centering = get_base_center(tmat);
    if (centering == A_FACE) {
      if (laue == LAUE2M) {
        debug_print("Monocli A to C\n");
        mat_copy_matrix_d3(correction_mat, monocli_a2c);
      } else {
        mat_copy_matrix_d3(correction_mat, a2c);
      }
      centering = C_FACE;
    }
    if (centering == B_FACE) {
      mat_copy_matrix_d3(correction_mat, b2c);
      centering = C_FACE;
    }
    if (laue == LAUE2M && centering == BODY) {
      debug_print("Monocli I to C\n");
      mat_copy_matrix_d3(correction_mat, monocli_i2c);
      centering = C_FACE;
    }
    break;

  case 3:
    /* hP (a=b) but not hR (a=b=c) */
    centering = R_CENTER;
    mat_multiply_matrix_id3(trans_corr_mat, tmat, rhombo_obverse);
    if (mat_is_int_matrix(trans_corr_mat, INT_PREC)) {
      mat_copy_matrix_d3(correction_mat, rhombo_obverse);
      debug_print("R-center observe setting\n");
      debug_print_matrix_d3(trans_corr_mat);
    }
    mat_multiply_matrix_id3(trans_corr_mat, tmat, rhomb_reverse);
    if (mat_is_int_matrix(trans_corr_mat, INT_PREC)) {
      mat_copy_matrix_d3(correction_mat, rhomb_reverse);
      debug_print("R-center reverse setting\n");
      debug_print_matrix_d3(trans_corr_mat);
    }
    break;

  case 4:
    centering = FACE;
    break;

  default:
    centering = CENTERING_ERROR;
    break;
  }

  return centering;
}

static Centering get_base_center(SPGCONST int tmat[3][3])
{
  int i;
  Centering centering = PRIMITIVE;

  debug_print("lat_get_base_center\n");

  /* C center */
  for (i = 0; i < 3; i++) {
    if (tmat[i][0] == 0 &&
        tmat[i][1] == 0 &&
        abs(tmat[i][2]) == 1) {
      centering = C_FACE;
      goto end;
    }
  }

  /* A center */
  for (i = 0; i < 3; i++) {
    if (abs(tmat[i][0]) == 1 &&
        tmat[i][1] == 0 &&
        tmat[i][2] == 0) {
      centering = A_FACE;
      goto end;
    }
  }

  /* B center */
  for (i = 0; i < 3; i++) {
    if (tmat[i][0] == 0 &&
        abs(tmat[i][1]) == 1 &&
        tmat[i][2] == 0) {
      centering = B_FACE;
      goto end;
    }
  }

  /* body center */
  if ((abs(tmat[0][0]) +
       abs(tmat[0][1]) +
       abs(tmat[0][2]) == 2) &&
      (abs(tmat[1][0]) +
       abs(tmat[1][1]) +
       abs(tmat[1][2]) == 2) &&
      (abs(tmat[2][0]) +
       abs(tmat[2][1]) +
       abs(tmat[2][2]) == 2)) {
    centering = BODY;
    goto end;
  }

  /* This should not happen. */
  warning_print("spglib: No centring was found (line %d, %s).\n", __LINE__, __FILE__);
  return PRIMITIVE;

 end:
  return centering;
}

static int get_centering_shifts(double shift[3][3],
                                const Centering centering)
{
  int i, j, multi;

  multi = 1;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      shift[i][j] = 0;
    }
  }

  if (centering != PRIMITIVE) {
    if (centering != FACE && centering != R_CENTER) {
      for (i = 0; i < 3; i++) { shift[0][i] = 0.5; } /* BODY */
      if (centering == A_FACE) { shift[0][0] = 0; }
      if (centering == B_FACE) { shift[0][1] = 0; }
      if (centering == C_FACE) { shift[0][2] = 0; }
      multi = 2;
    }

    if (centering == R_CENTER) {
      shift[0][0] = 2. / 3;
      shift[0][1] = 1. / 3;
      shift[0][2] = 1. / 3;
      shift[1][0] = 1. / 3;
      shift[1][1] = 2. / 3;
      shift[1][2] = 2. / 3;
      multi = 3;
    }

    if (centering == FACE) {
      shift[0][0] = 0;
      shift[0][1] = 0.5;
      shift[0][2] = 0.5;
      shift[1][0] = 0.5;
      shift[1][1] = 0;
      shift[1][2] = 0.5;
      shift[2][0] = 0.5;
      shift[2][1] = 0.5;
      shift[2][2] = 0;
      multi = 4;
    }
  }

  return multi;
}

/* allow_flip: Flip axes when abs(P)=I but P!=I. */
/* lattice is overwritten. */
static int is_equivalent_lattice(double tmat[3][3],
                                 const int mode,
                                 SPGCONST double lattice[3][3],
                                 SPGCONST double orig_lattice[3][3],
                                 const double symprec)
{
  int i, j;
  double inv_lat[3][3], tmat_abs[3][3], metric[3][3], orig_metric[3][3];
  int tmat_int[3][3];

  if (mat_Dabs(mat_get_determinant_d3(lattice) -
               mat_get_determinant_d3(orig_lattice)) > symprec) {
    goto fail;
  }

  if (!mat_inverse_matrix_d3(inv_lat, lattice, symprec)) {
    goto fail;
  }

  mat_multiply_matrix_d3(tmat, inv_lat, orig_lattice);

  switch (mode) {
  case 0:  /* Check identity of all elements */
    if (mat_check_identity_matrix_d3(identity, tmat, symprec)) {
      return 1;
    }
    break;

  case 1:  /* Check identity of all elements allowing axes flips */
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        tmat_abs[i][j] = mat_Dabs(tmat[i][j]);
      }
    }

    if (mat_check_identity_matrix_d3(identity, tmat_abs, symprec)) {
      return 1;
    }
    break;

  case 2:  /* Check metric tensors */
    mat_cast_matrix_3d_to_3i(tmat_int, tmat);
    if (! mat_check_identity_matrix_id3(tmat_int, tmat, symprec)) {
      break;
    }
    if (mat_get_determinant_i3(tmat_int) != 1) {
      break;
    }
    mat_get_metric(orig_metric, orig_lattice);
    mat_get_metric(metric, lattice);
    if (mat_check_identity_matrix_d3(orig_metric, metric, symprec)) {
      return 1;
    }
    break;
  }

fail:
  return 0;
}
