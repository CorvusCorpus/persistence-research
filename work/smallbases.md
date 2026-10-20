# Multiplicative persistence in small bases: complete antecedent trees

For a base b and a nonzero final digit d, the tree T_d consists of the digit products whose base-b
digit-product orbit ends at d. Each node row lists EVERY integer y (all sizes) that is a product of
primes below b, has no 0 digit, and has digit product equal to the node; when every node of T_d has
been treated, T_d is complete and every positive integer whose orbit ends at d has persistence at
most height(T_d) + 1.

These statements concern nonzero final digits only. As in base 10, numbers whose orbit reaches 0
(where the longest persistences live) are not covered.

## Base 3

| final digit | depth | node (base 10) | node (base 3) | largest digit set | proved length bound | all preimages |
|---|---|---|---|---|---|---|
| 1 | 0 | 1 | 1 | 0 | 44 | 4 (11_3) |
| 1 | 1 | 4 | 11 | 2 | 93 | 8 (22_3) |
| 1 | 2 | 8 | 22 | 3 | 94 | none |
| 2 | 0 | 2 | 2 | 1 | 94 | 16 (121_3) |
| 2 | 1 | 16 | 121 | 4 | 93 | none |

- **Theorem (base 3, final digit 1).** The tree has exactly 3 node(s) and height 2; every positive integer whose base-3 digit-product orbit ends at 1 has persistence at most 3, and this is attained.
- **Theorem (base 3, final digit 2).** The tree has exactly 2 node(s) and height 1; every positive integer whose base-3 digit-product orbit ends at 2 has persistence at most 2, and this is attained.

## Base 4

| final digit | depth | node (base 10) | node (base 4) | largest digit set | proved length bound | all preimages |
|---|---|---|---|---|---|---|
| 1 | 0 | 1 | 1 | 0 | 36 | none |
| 2 | 0 | 2 | 2 | 1 | 76 | 6 (12_4), 9 (21_4) |
| 2 | 1 | 6 | 12 | 2 | 75 | 27 (123_4), 54 (312_4) |
| 2 | 1 | 9 | 21 | 2 | 75 | none |
| 2 | 2 | 27 | 123 | 3 | 76 | none |
| 2 | 2 | 54 | 312 | 4 | 75 | none |
| 3 | 0 | 3 | 3 | 1 | 76 | none |

- **Theorem (base 4, final digit 1).** The tree has exactly 1 node(s) and height 0; every positive integer whose base-4 digit-product orbit ends at 1 has persistence at most 1, and this is attained.
- **Theorem (base 4, final digit 2).** The tree has exactly 5 node(s) and height 2; every positive integer whose base-4 digit-product orbit ends at 2 has persistence at most 3, and this is attained.
- **Theorem (base 4, final digit 3).** The tree has exactly 1 node(s) and height 0; every positive integer whose base-4 digit-product orbit ends at 3 has persistence at most 1, and this is attained.

## Base 6

| final digit | depth | node (base 10) | node (base 6) | largest digit set | proved length bound | all preimages |
|---|---|---|---|---|---|---|
| 1 | 0 | 1 | 1 | 0 | 29 | none |
| 2 | 0 | 2 | 2 | 1 | 60 | 8 (12_6) |
| 2 | 1 | 8 | 12 | 3 | 59 | 16 (24_6), 512 (2212_6) |
| 2 | 2 | 16 | 24 | 4 | 59 | 64 (144_6), 1600 (11224_6) |
| 2 | 2 | 512 | 2212 | 9 | 64 | none |
| 2 | 3 | 64 | 144 | 6 | 61 | 10000 (114144_6), 64000 (1212144_6), 80000 (1414212_6), 390625 (12212241_6) |
| 2 | 3 | 1600 | 11224 | 8 | 61 | none |
| 2 | 4 | 10000 | 114144 | 8 | 61 | none |
| 2 | 4 | 64000 | 1212144 | 12 | 78 | none |
| 2 | 4 | 80000 | 1414212 | 11 | 71 | none |
| 2 | 4 | 390625 | 12212241 | 8 | 77 | none |
| 3 | 0 | 3 | 3 | 1 | 60 | 9 (13_6), 45 (113_6) |
| 3 | 1 | 9 | 13 | 2 | 59 | none |
| 3 | 1 | 45 | 113 | 3 | 60 | none |
| 4 | 0 | 4 | 4 | 2 | 60 | 10 (14_6), 25 (41_6), 50 (122_6), 80 (212_6) |
| 4 | 1 | 10 | 14 | 2 | 59 | 32 (52_6) |
| 4 | 1 | 25 | 41 | 2 | 59 | none |
| 4 | 1 | 50 | 122 | 3 | 60 | none |
| 4 | 1 | 80 | 212 | 5 | 60 | 10240 (115224_6) |
| 4 | 2 | 32 | 52 | 5 | 60 | 100 (244_6), 160 (424_6), 32768 (411412_6) |
| 4 | 2 | 10240 | 115224 | 12 | 78 | none |
| 4 | 3 | 100 | 244 | 4 | 60 | none |
| 4 | 3 | 160 | 424 | 6 | 60 | 640 (2544_6), 1250 (5442_6), 3125 (22245_6), 3200 (22452_6) |
| 4 | 3 | 32768 | 411412 | 15 | 78 | none |
| 4 | 4 | 640 | 2544 | 8 | 62 | none |
| 4 | 4 | 1250 | 5442 | 5 | 60 | none |
| 4 | 4 | 3125 | 22245 | 5 | 60 | none |
| 4 | 4 | 3200 | 22452 | 9 | 65 | none |
| 5 | 0 | 5 | 5 | 1 | 60 | none |

- **Theorem (base 6, final digit 1).** The tree has exactly 1 node(s) and height 0; every positive integer whose base-6 digit-product orbit ends at 1 has persistence at most 1, and this is attained.
- **Theorem (base 6, final digit 2).** The tree has exactly 10 node(s) and height 4; every positive integer whose base-6 digit-product orbit ends at 2 has persistence at most 5, and this is attained.
- **Theorem (base 6, final digit 3).** The tree has exactly 3 node(s) and height 1; every positive integer whose base-6 digit-product orbit ends at 3 has persistence at most 2, and this is attained.
- **Theorem (base 6, final digit 4).** The tree has exactly 14 node(s) and height 4; every positive integer whose base-6 digit-product orbit ends at 4 has persistence at most 5, and this is attained.
- **Theorem (base 6, final digit 5).** The tree has exactly 1 node(s) and height 0; every positive integer whose base-6 digit-product orbit ends at 5 has persistence at most 1, and this is attained.

