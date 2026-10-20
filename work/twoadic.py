"""Exact maximal 2-adic valuation of y over all integers y whose decimal digits are a fixed
multiset D of non-1 digits plus any number of 1s (re-implementation of the recursion in
Fonga, arXiv 2608.27802, Thm 4.1; the only input is the elementary fact that
v2(T + c*10^t) = v2(T) whenever t + v2(c) > v2(T)).

9y = -1 + sum_i 9(d_i-1) 10^{a_i} + 10^L ; v2(y) = v2(9y)."""
from functools import lru_cache

def v2(n):
    n = abs(n)
    if n == 0: return 10**9
    return (n & -n).bit_length() - 1

def tmax(bag):
    bag = tuple(sorted(bag))
    @lru_cache(maxsize=None)
    def V(rem, a, T):
        best = v2(T)
        lim = v2(T)                         # only positions t <= v2(T) can change the valuation
        if rem:
            for i, d in enumerate(rem):
                if i and rem[i] == rem[i-1]: continue
                nrem = rem[:i] + rem[i+1:]
                for t in range(a+1, lim+1):
                    best = max(best, V(nrem, t, T + 9*(d-1)*10**t))
        else:
            for L in range(a+1, lim+1):
                best = max(best, v2(T + 10**L))
        return best
    best = 0                                # units digit 1 or odd -> valuation 0
    for i, d in enumerate(bag):
        if d % 2 == 0 and (i == 0 or bag[i] != bag[i-1]):
            best = max(best, V(bag[:i] + bag[i+1:], 0, 9*(d-1) - 1))
    return best

if __name__ == "__main__":
    import itertools
    for b, want in (((4,4,7),7), ((2,7,8),9), ((2,2,4,7),15), ((2,2,2,2,7),13)):
        print(b, tmax(b), "published", want)
    # brute force cross-check: all arrangements with up to 9 extra ones
    def brute(bag, extra):
        best = 0
        for n1 in range(extra+1):
            digs = list(bag) + [1]*n1
            for perm in set(itertools.permutations(digs)):
                y = int(''.join(map(str, perm)))
                best = max(best, v2(y))
        return best
    for bag in ((4,4,7),(2,7,8),(2,2,9),(8,9),(2,4,3,3),(2,6,6),(2,2,2,3,3)):
        print(bag, "recursion", tmax(bag), "brute(<= +6 ones)", brute(bag, 6 if len(bag) <= 3 else 4))
