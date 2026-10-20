#!/usr/bin/env python3
"""Exact, independent check of MEMBER lines of the chain sieve (a MEMBER line is a new chain number).

For each line: recompute x = 2^a 3^b 7^c (form 0) or 3^a 5^b 7^c (form 1) with Python integers, its digit count, the
chain condition (no 0, not both a 5 and an even digit) and its multiplicative persistence orbit x -> p(x) -> ... A
7-smooth x >= 10 is the digit product of some n (spell x as a product of digits 2..9), so a chain number of persistence k
gives a number of persistence k + 1: persistence(x) >= 11 means a persistence-12 number exists.
usage: check_member.py FILE_OR_LINE [...]      (log files, or lines like 'MEMBER form=0 a=.. b=.. c=..')
"""
import os, sys
sys.set_int_max_str_digits(0)

def digit_product(n):
    s = str(n)
    r = 1
    for d in '23456789': r *= int(d) ** s.count(d)
    return 0 if '0' in s else r

def short(n):
    s = str(n)
    return s if len(s) <= 60 else f"{s[:25]}...{s[-25:]} ({len(s)} digits)"

def check(form, a, b, c):
    x = (2 ** a * 3 ** b if form == 0 else 3 ** a * 5 ** b) * 7 ** c
    s = str(x)
    chain = x >= 10 and '0' not in s and not ('5' in s and any(d in s for d in '2468'))
    orbit = [x]
    while orbit[-1] >= 10: orbit.append(digit_product(orbit[-1]))
    pers = len(orbit) - 1
    print(f"form={form} a={a} b={b} c={c}: x has {len(s)} digits; chain number: {chain}; persistence of x: {pers}")
    for i, y in enumerate(orbit): print(f"  step {i}: {short(y)}")
    if chain and pers >= 11: print("  !!! persistence(x) >= 11: x is the digit product of a number with persistence >= 12 !!!")
    return chain, pers

lines = []
for arg in sys.argv[1:]:
    lines += open(arg).read().splitlines() if os.path.exists(arg) else [arg]
found = 0
for line in lines:
    if not line.startswith('MEMBER'): continue
    d = dict(kv.split('=') for kv in line.split()[1:])
    check(int(d['form']), int(d['a']), int(d['b']), int(d['c'])); found += 1
print(f"{found} MEMBER line(s) checked")
