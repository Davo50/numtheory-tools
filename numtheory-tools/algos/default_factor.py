from sympy.abc import x
from sympy import poly, sqrt, lcm_list, gcd, isprime
from sympy.ntheory.primetest import is_square


def fermat_factorization(N):
    steps = []
    x = int(N**0.5) + 1
    while True:
        y = x ** 2 - N
        steps.append((x, x**2, y))
        if is_square(y):
            y = sqrt(y)
            return (x - y, x + y), steps
        x += 1


def p1_pollard(N, B):
    steps = []
    k = int(lcm_list(range(2, B + 1)))
    steps.append(k)
    for a in range(2, N):
        a_ = pow(a, k, N)
        d = gcd(a_ - 1, N)
        steps.append((a, a_, d))
        if 1 < d < N:
            return (d, N // d), steps


def rho_pollard(N, f, x0):
    steps = []
    x1 = f(x0)
    x2 = f(f(x0))
    d = gcd(x1 - x2, N)
    steps.append((x1, x2, d))
    while d == 1:
        x1 = f(x1)
        x2 = f(f(x2))
        d = gcd(x1 - x2, N)
        steps.append((x1, x2, d))
    return (d, N // d), steps


def try_double(P, a, N):
    d = gcd(P[1], N)
    if d != 1:
        return (), d
    lam = ((3 * P[0] + a) * pow(2 * P[1], -1, N)) % N
    x3 = (pow(lam, 2, N) - P[0] - P[0]) % N
    return (x3, (-(lam * (x3 - P[0]) + P[1])) % N), 1


def try_add(P, Q, N):
    d = gcd(P[0] - Q[0], N)
    if d != 1:
        return (), d
    lam = ((P[1] - Q[1]) * pow(P[0] - Q[0], -1, N)) % N
    x3 = (pow(lam, 2, N) - P[0] - Q[0]) % N
    return (x3, (-(lam * (x3 - P[0]) + P[1])) % N), 1


def try_scalar_mult(k, P, a, N):
    l = k.bit_length()
    Q = P
    for i in range(l - 2, -1, -1):
        Q, v = try_double(Q, a, N)
        if v != 1:
            return v
        if (k & (1 << (i))):
            Q, v = try_add(P, Q, N)
            if v != 1:
                return v
    return 1


def lenstra(N, B, x, y, a):
    b = y**2 - x**3 - a*x
    d = gcd(-16 * (4 * a ** 3 + 27 * b ** 2), N)
    if 1 < d < N:
        return (d, N // d)
    k = 1
    P = (x, y)
    for p in range(2, B + 1):
        if isprime(p):
            k *= p
            d = try_scalar_mult(k, P, a, N)
            if d != 1:
                return (d, N // d)



# print(rho_pollard(1359331, lambda v: poly("x**2 + 5")(v) % 1359331, 1))
