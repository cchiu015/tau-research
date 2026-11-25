pari.allocatemem(2*10^9)
pari.default("parisizemax", 4*10^9)

_TNF = {}

def generate_F_d_minus_1(d):
    """
    Generates the homogeneous polynomial F_{d-1}(X, Y) in ZZ[X, Y]
    of degree m = (d - 1) // 2, satisfying tau(p^(d-1)) = F_{d-1}(p^11, tau(p)^2).
    """
    d = ZZ(d)
    if d <= 5 or d % 2 == 0:
        raise ValueError("Expected odd d > 5 so deg(F_{d-1}) >= 3.")
    m = (d - 1) // 2
    R.<X, Y> = PolynomialRing(ZZ)
    return sum(((-1)^j) * binomial(2*m - j, j) * (X^j) * (Y^(m - j)) for j in range(m + 1))

def get_tnf(d, certify=False):
    key = (d, certify)
    if key not in _TNF:
        F = generate_F_d_minus_1(d)
        S.<x> = PolynomialRing(ZZ)
        f = S(F(x, 1))
        _TNF[key] = pari.thueinit(f, 1 if certify else 0)
    return _TNF[key]

def solve_thue(ell, d, certify=False):
    tnf = get_tnf(d, certify)
    raw = pari.thue(tnf, ell)
    all_sols = [(ZZ(s[0]), ZZ(s[1])) for s in raw]

    tau_sols = []
    for X_val, Y_val in all_sols:
        if X_val > 1 and Y_val >= 0:
            p_cand, is_11th = X_val.nth_root(11, truncate_mode=True)
            if is_11th and p_cand.is_prime() and Y_val.is_square():
                tau_p = delta_qexp(p_cand + 1)[p_cand]
                if tau_p^2 == Y_val:
                    tau_sols.append((p_cand, tau_p))
    return all_sols, tau_sols

l_and_d = [(-1319,659),(1381,23),(-1609,23),(1789,179),(1801,53),
    (1933,7),(1933,23),(1933,967)]
# order in terms of size of d
ordered_pairs = sorted(l_and_d, key=lambda x: x[1])
# solve thue eqns for ordered pairs
for pair in ordered_pairs[:5]:
    print(pair)
    print(solve_thue(pair[0],pair[1],True))