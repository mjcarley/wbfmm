#include <glib.h>

#include "wbfmm-private.h"

gint _wbfmm_shift_angles[] = {
     36,  -9,   5,  18,    31,  -9,   5,  17,    25,  -9,   5,  16,
     24,  -9,   5,  15,    23,  -9,   5,  16,    17,  -9,   5,  17,
     12,  -9,   5,  18,    37,  -9,   4,  17,    33,  -9,   4,  14,
     26,  -9,   4,  13,    24,  -9,   4,  12,    22,  -9,   4,  13,
     15,  -9,   4,  14,    11,  -9,   4,  17,    39,  -9,   2,  16,
     34,  -9,   2,  13,    27,  -9,   2,  10,    24,  -9,   2,   9,
     21,  -9,   2,  10,    14,  -9,   2,  13,     9,  -9,   2,  16,
     40,  -9,  -1,  15,    35,  -9,  -1,  12,    28,  -9,  -1,   9,
     24,  -9,  -1,   8,    20,  -9,  -1,   9,    13,  -9,  -1,  12,
      8,  -9,  -1,  15,    39,  -9,  -2,  16,    34,  -9,  -2,  13,
     27,  -9,  -2,  10,    24,  -9,  -2,   9,    21,  -9,  -2,  10,
     14,  -9,  -2,  13,     9,  -9,  -2,  16,    37,  -9,  -4,  17,
     33,  -9,  -4,  14,    26,  -9,  -4,  13,    24,  -9,  -4,  12,
     22,  -9,  -4,  13,    15,  -9,  -4,  14,    11,  -9,  -4,  17,
     36,  -9,  -5,  18,    31,  -9,  -5,  17,    25,  -9,  -5,  16,
     24,  -9,  -5,  15,    23,  -9,  -5,  16,    17,  -9,  -5,  17,
     12,  -9,  -5,  18,    37,  -9,   6,  17,    33,  -9,   6,  14,
     26,  -9,   6,  13,    24,  -9,   6,  12,    22,  -9,   6,  13,
     15,  -9,   6,  14,    11,  -9,   6,  17,    41,  -9,   5,  14,
     36,  -9,   5,  11,    29,  -9,   5,   8,    24,  -9,   5,   7,
     19,  -9,   5,   8,    12,  -9,   5,  11,     7,  -9,   5,  14,
     42,  -9,   3,  13,    38,  -9,   3,   8,    30,  -9,   3,   6,
     24,  -9,   3,   5,    18,  -9,   3,   6,    10,  -9,   3,   8,
      6,  -9,   3,  13,    44,  -9,  -1,  12,    40,  -9,  -1,   7,
     32,  -9,  -1,   5,    24,  -9,  -1,   4,    16,  -9,  -1,   5,
      8,  -9,  -1,   7,     4,  -9,  -1,  12,    42,  -9,  -3,  13,
     38,  -9,  -3,   8,    30,  -9,  -3,   6,    24,  -9,  -3,   5,
     18,  -9,  -3,   6,    10,  -9,  -3,   8,     6,  -9,  -3,  13,
     41,  -9,  -5,  14,    36,  -9,  -5,  11,    29,  -9,  -5,   8,
     24,  -9,  -5,   7,    19,  -9,  -5,   8,    12,  -9,  -5,  11,
      7,  -9,  -5,  14,    37,  -9,  -6,  17,    33,  -9,  -6,  14,
     26,  -9,  -6,  13,    24,  -9,  -6,  12,    22,  -9,  -6,  13,
     15,  -9,  -6,  14,    11,  -9,  -6,  17,    39,  -9,   8,  16,
     34,  -9,   8,  13,    27,  -9,   8,  10,    24,  -9,   8,   9,
     21,  -9,   8,  10,    14,  -9,   8,  13,     9,  -9,   8,  16,
     42,  -9,   7,  13,    38,  -9,   7,   8,    30,  -9,   7,   6,
     24,  -9,   7,   5,    18,  -9,   7,   6,    10,  -9,   7,   8,
      6,  -9,   7,  13,    46,  -9,   5,  10,    43,  -9,   5,   6,
     36,  -9,   5,   3,    24,  -9,   5,   2,    12,  -9,   5,   3,
      5,  -9,   5,   6,     2,  -9,   5,  10,    47,  -9,  -1,   9,
     45,  -9,  -1,   5,    40,  -9,  -1,   2,    24,  -9,  -1,   1,
      8,  -9,  -1,   2,     3,  -9,  -1,   5,     1,  -9,  -1,   9,
     46,  -9,  -5,  10,    43,  -9,  -5,   6,    36,  -9,  -5,   3,
     24,  -9,  -5,   2,    12,  -9,  -5,   3,     5,  -9,  -5,   6,
      2,  -9,  -5,  10,    42,  -9,  -7,  13,    38,  -9,  -7,   8,
     30,  -9,  -7,   6,    24,  -9,  -7,   5,    18,  -9,  -7,   6,
     10,  -9,  -7,   8,     6,  -9,  -7,  13,    39,  -9,  -8,  16,
     34,  -9,  -8,  13,    27,  -9,  -8,  10,    24,  -9,  -8,   9,
     21,  -9,  -8,  10,    14,  -9,  -8,  13,     9,  -9,  -8,  16,
     40, -13,   9,  15,    35, -13,   9,  12,    28, -13,   9,   9,
     24, -13,   9,   8,    20, -13,   9,   9,    13, -13,   9,  12,
      8, -13,   9,  15,    44, -13,   9,  12,    40, -13,   9,   7,
     32, -13,   9,   5,    24, -13,   9,   4,    16, -13,   9,   5,
      8, -13,   9,   7,     4, -13,   9,  12,    47, -13,   9,   9,
     45, -13,   9,   5,    40, -13,   9,   2,    24, -13,   9,   1,
      8, -13,   9,   2,     3, -13,   9,   5,     1, -13,   9,   9,
     48,  -1,  -1,   8,    48,  -1,  -1,   4,    48,  -1,  -1,   1,
     63,  64,  64,   0,     0,  -1,  -1,   1,     0,  -1,  -1,   4,
      0,  -1,  -1,   8,    47,   5,  -9,   9,    45,   5,  -9,   5,
     40,   5,  -9,   2,    24,   5,  -9,   1,     8,   5,  -9,   2,
      3,   5,  -9,   5,     1,   5,  -9,   9,    44,   5,  -9,  12,
     40,   5,  -9,   7,    32,   5,  -9,   5,    24,   5,  -9,   4,
     16,   5,  -9,   5,     8,   5,  -9,   7,     4,   5,  -9,  12,
     40,   5,  -9,  15,    35,   5,  -9,  12,    28,   5,  -9,   9,
     24,   5,  -9,   8,    20,   5,  -9,   9,    13,   5,  -9,  12,
      8,   5,  -9,  15,    39,  -9,  10,  16,    34,  -9,  10,  13,
     27,  -9,  10,  10,    24,  -9,  10,   9,    21,  -9,  10,  10,
     14,  -9,  10,  13,     9,  -9,  10,  16,    42,  -9,  11,  13,
     38,  -9,  11,   8,    30,  -9,  11,   6,    24,  -9,  11,   5,
     18,  -9,  11,   6,    10,  -9,  11,   8,     6,  -9,  11,  13,
     46,  -9,  13,  10,    43,  -9,  13,   6,    36,  -9,  13,   3,
     24,  -9,  13,   2,    12,  -9,  13,   3,     5,  -9,  13,   6,
      2,  -9,  13,  10,    47,  -9,  17,   9,    45,  -9,  17,   5,
     40,  -9,  17,   2,    24,  -9,  17,   1,     8,  -9,  17,   2,
      3,  -9,  17,   5,     1,  -9,  17,   9,    46,  -9, -13,  10,
     43,  -9, -13,   6,    36,  -9, -13,   3,    24,  -9, -13,   2,
     12,  -9, -13,   3,     5,  -9, -13,   6,     2,  -9, -13,  10,
     42,  -9, -11,  13,    38,  -9, -11,   8,    30,  -9, -11,   6,
     24,  -9, -11,   5,    18,  -9, -11,   6,    10,  -9, -11,   8,
      6,  -9, -11,  13,    39,  -9, -10,  16,    34,  -9, -10,  13,
     27,  -9, -10,  10,    24,  -9, -10,   9,    21,  -9, -10,  10,
     14,  -9, -10,  13,     9,  -9, -10,  16,    37,  -9,  12,  17,
     33,  -9,  12,  14,    26,  -9,  12,  13,    24,  -9,  12,  12,
     22,  -9,  12,  13,    15,  -9,  12,  14,    11,  -9,  12,  17,
     41,  -9,  13,  14,    36,  -9,  13,  11,    29,  -9,  13,   8,
     24,  -9,  13,   7,    19,  -9,  13,   8,    12,  -9,  13,  11,
      7,  -9,  13,  14,    42,  -9,  15,  13,    38,  -9,  15,   8,
     30,  -9,  15,   6,    24,  -9,  15,   5,    18,  -9,  15,   6,
     10,  -9,  15,   8,     6,  -9,  15,  13,    44,  -9,  17,  12,
     40,  -9,  17,   7,    32,  -9,  17,   5,    24,  -9,  17,   4,
     16,  -9,  17,   5,     8,  -9,  17,   7,     4,  -9,  17,  12,
     42,  -9, -15,  13,    38,  -9, -15,   8,    30,  -9, -15,   6,
     24,  -9, -15,   5,    18,  -9, -15,   6,    10,  -9, -15,   8,
      6,  -9, -15,  13,    41,  -9, -13,  14,    36,  -9, -13,  11,
     29,  -9, -13,   8,    24,  -9, -13,   7,    19,  -9, -13,   8,
     12,  -9, -13,  11,     7,  -9, -13,  14,    37,  -9, -12,  17,
     33,  -9, -12,  14,    26,  -9, -12,  13,    24,  -9, -12,  12,
     22,  -9, -12,  13,    15,  -9, -12,  14,    11,  -9, -12,  17,
     36,  -9,  13,  18,    31,  -9,  13,  17,    25,  -9,  13,  16,
     24,  -9,  13,  15,    23,  -9,  13,  16,    17,  -9,  13,  17,
     12,  -9,  13,  18,    37,  -9,  14,  17,    33,  -9,  14,  14,
     26,  -9,  14,  13,    24,  -9,  14,  12,    22,  -9,  14,  13,
     15,  -9,  14,  14,    11,  -9,  14,  17,    39,  -9,  16,  16,
     34,  -9,  16,  13,    27,  -9,  16,  10,    24,  -9,  16,   9,
     21,  -9,  16,  10,    14,  -9,  16,  13,     9,  -9,  16,  16,
     40,  -9,  17,  15,    35,  -9,  17,  12,    28,  -9,  17,   9,
     24,  -9,  17,   8,    20,  -9,  17,   9,    13,  -9,  17,  12,
      8,  -9,  17,  15,    39,  -9, -16,  16,    34,  -9, -16,  13,
     27,  -9, -16,  10,    24,  -9, -16,   9,    21,  -9, -16,  10,
     14,  -9, -16,  13,     9,  -9, -16,  16,    37,  -9, -14,  17,
     33,  -9, -14,  14,    26,  -9, -14,  13,    24,  -9, -14,  12,
     22,  -9, -14,  13,    15,  -9, -14,  14,    11,  -9, -14,  17,
     36,  -9, -13,  18,    31,  -9, -13,  17,    25,  -9, -13,  16,
     24,  -9, -13,  15,    23,  -9, -13,  16,    17,  -9, -13,  17,
     12,  -9, -13,  18} ;
