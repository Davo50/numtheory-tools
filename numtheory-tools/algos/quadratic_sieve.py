from sage.all import is_prime, Matrix, GF, vector, ZZ, gcd


def legendre_symbol(a, p):
    return pow(a, (p - 1) // 2, p)


def try_factor(Q, FB):
    res = []
    for p in FB:
        c = 0
        while Q % p == 0:
            c += 1
            Q //= p
        res.append(c)
    if Q == 1:
        return True, res
    return False, None


def quadratic_sieve(N, B):
    FB = []
    for i in range(2, B + 1):
        if is_prime(i):
            if legendre_symbol(N, i) == 1:
                FB.append(i)

    sqrt_N = int(N**0.5)
    i = -1
    M = []
    T = []
    while True:
        v = i + sqrt_N
        Q = v ** 2 - N
        done, res = try_factor(abs(Q), FB)
        if done:
            if Q < 0:
                res = vector(ZZ, [1] + res)
            else:
                res = vector(ZZ, [0] + res)
            M.append(res)
            T.append(v)
        i += 1

        if len(M) > len(FB) * 1.2:
            MM = Matrix(GF(2), M)
            K = MM.left_kernel_matrix()
            for r in K:
                r = list(r)
                if r.count(1) == 2:
                    i1 = r.index(1)
                    i2 = r.index(1, i1 + 1)

                    print(T[i1], T[i2])


                    ll = T[i1] * T[i2]
                    rr = 1
                    for p, a in zip([-1] + FB, M[i1] + M[i2]):
                        rr *= pow(p, a // 2)
                    d = gcd(rr - ll, N)
                    if 1 < d < N:
                        return d, N // d
