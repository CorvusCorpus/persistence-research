#!/bin/bash
# tsieve6 against gsieve on the same main-path region: every DONE counter, then the S2/R3/R4/S4/MEMBER dump sets.
# usage: val_tsieve6.sh PARAMFILE G|E M N [N0] [threads] [GSIEVE_EXE] [TSIEVE_EXE]
P=$1; D=$2; M=$3; N=$4; N0=${5:--1}; T=${6:-2}; G=${7:-../gsieve/gsieve_k13}; X=${8:-./tsieve6}
S=${TMPDIR:-/tmp}/vt6.$$; mkdir -p $S
n0arg=(); [ "$N0" -ge 0 ] && n0arg=(n0=$N0)
$G main $D $N 0 $((N + 1)) dump "${n0arg[@]}" threads=$T > $S/g.out 2>/dev/null
$X $P $D $M $N "${n0arg[@]}" threads=$T dump > $S/t.out
g=$(grep ^DONE $S/g.out); t=$(grep ^DONE $S/t.out)
ok=1
for k in points s1 s2 s3 s4 gmp full members; do
    gv=$(echo "$g" | grep -o " $k=[0-9]*" | cut -d= -f2); tv=$(echo "$t" | grep -o " $k=[0-9]*" | cut -d= -f2)
    [ -n "$gv" ] && [ "$gv" = "$tv" ] || { ok=0; echo "  $k: gsieve $gv tsieve6 $tv"; }
done
for k in S2 R3 R4 S4 MEMBER; do
    cmp -s <(grep "^$k " $S/g.out | sort) <(grep "^$k " $S/t.out | sort) || { ok=0; echo "  $k lines differ"; }
done
echo "$D N=$N N0=$N0 M=$M: $(echo "$t" | grep -o 'points=.*members=[0-9]*') -> $([ $ok = 1 ] && echo IDENTICAL || echo DIFFERENT)"
rm -rf $S; [ $ok = 1 ]
