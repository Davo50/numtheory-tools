"""
Калькулятор числовых свойств для зачёта.

Для целого n (или пары a, n) выдаёт всё, что нужно:
— разложение n на простые
— функция Эйлера φ(n), функция Кармайкла λ(n)
— количество и список делителей
— первообразные корни (если n простое)
— порядок элемента a mod n
— символы Лежандра/Якоби (a/n)
— B-гладкость для разных B
"""
from sympy import factorint, isprime, totient, reduced_totient, divisors, jacobi_symbol, gcd as sympy_gcd, primerange
from functools import reduce
from math import lcm


def _lcm_list(xs):
    return reduce(lcm, xs, 1)


def number_properties(n):
    """Собирает ВСЕ базовые свойства числа n в один словарь."""
    n = int(n)
    if n < 2:
        return {'error': 'n должно быть ≥ 2'}

    fact = factorint(n)
    phi = int(totient(n))
    lam = int(reduced_totient(n))  # Carmichael
    divs = [int(d) for d in divisors(n)]
    is_prime = isprime(n)

    # B-гладкость для разных B
    largest_prime = max(fact.keys())
    smoothness = {}
    for B in [5, 11, 13, 19, 31, 50, 100]:
        smoothness[B] = largest_prime <= B

    result = {
        'n': n,
        'is_prime': is_prime,
        'factorization': {int(p): int(e) for p, e in fact.items()},
        'factorization_str': ' · '.join(f"{p}^{e}" if e > 1 else str(p)
                                        for p, e in sorted(fact.items())),
        'phi': phi,
        'lambda_carmichael': lam,
        'num_divisors': len(divs),
        'divisors': divs if len(divs) <= 30 else divs[:15] + ['…'] + divs[-5:],
        'divisors_total': len(divs),
        'largest_prime_factor': int(largest_prime),
        'smoothness': smoothness,
    }

    # Первообразные корни — только для простого n
    if is_prime:
        result['primitive_roots'] = _find_primitive_roots(n, limit=10)

    return result


def _find_primitive_roots(p, limit=10):
    """Первые limit первообразных корней по простому модулю p.
    Первообразный корень = элемент порядка p−1."""
    if not isprime(p):
        return None
    n = p - 1
    # Для каждого простого делителя q числа n: элемент g не должен давать g^(n/q) = 1
    q_factors = list(factorint(n).keys())
    roots = []
    for g in range(2, p):
        if all(pow(g, n // q, p) != 1 for q in q_factors):
            roots.append(g)
            if len(roots) >= limit:
                break
    return roots


def element_order(a, n):
    """Порядок элемента a по модулю n.
    Возвращает словарь с результатом и шагами вычисления."""
    a, n = int(a), int(n)
    steps = []

    if int(sympy_gcd(a, n)) != 1:
        return {
            'order': None,
            'steps': [{'k': 'Ошибка',
                       'description': f"НОД({a}, {n}) = {int(sympy_gcd(a, n))} ≠ 1 — "
                                      f"элемент {a} не обратим по модулю {n}, "
                                      f"порядок не определён."}]
        }

    # Порядок делит λ(n) (функция Кармайкла). Перебираем делители λ(n)
    # в порядке возрастания — первый, при котором a^d ≡ 1, и есть порядок.
    lam = int(reduced_totient(n))
    divs = sorted([int(d) for d in divisors(lam)])

    steps.append({
        'k': 'init',
        'description': f"Ищем порядок a = {a} по модулю n = {n}. "
                       f"По теореме Лагранжа порядок делит λ(n) = {lam}."
    })
    steps.append({
        'k': 'divisors',
        'description': f"Делители λ(n) = {lam}: {divs}"
    })

    for d in divs:
        val = pow(a, d, n)
        steps.append({
            'k': d,
            'description': f"{a}^{d} mod {n} = {val}{'   ← порядок!' if val == 1 else ''}"
        })
        if val == 1:
            return {'order': d, 'steps': steps}

    return {'order': None, 'steps': steps}


def modular_sqrt(a, p):
    """
    Квадратный корень из a по простому модулю p (алгоритм Тонелли-Шенкса).
    Возвращает одно из решений (второе = −x mod p), или None.
    """
    a = int(a) % int(p)
    p = int(p)
    if a == 0:
        return 0
    if p == 2:
        return a
    # Эйлеров критерий
    if pow(a, (p - 1) // 2, p) != 1:
        return None

    # p ≡ 3 mod 4 — простой случай
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)

    # Общий случай Тонелли-Шенкса
    s, q = 0, p - 1
    while q % 2 == 0:
        q //= 2
        s += 1

    # Находим z — невычет
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1

    M = s
    c = pow(z, q, p)
    t = pow(a, q, p)
    R = pow(a, (q + 1) // 2, p)

    while True:
        if t == 1:
            return R
        i = 0
        temp = t
        while temp != 1 and i < M:
            temp = pow(temp, 2, p)
            i += 1
        b = pow(c, 1 << (M - i - 1), p)
        M = i
        c = pow(b, 2, p)
        t = (t * c) % p
        R = (R * b) % p


def legendre_symbol_info(a, p):
    """Символ Лежандра (a/p) со смыслом: 1 = квадратичный вычет, -1 = невычет."""
    a, p = int(a), int(p)
    if not isprime(p):
        return {'error': f'p = {p} должно быть простым для символа Лежандра'}

    a_mod = a % p
    if a_mod == 0:
        return {'value': 0, 'meaning': f"{a} ≡ 0 (mod {p}), символ = 0"}

    # Вычисляем через Эйлеров критерий: a^((p-1)/2) mod p
    power = pow(a_mod, (p - 1) // 2, p)
    value = 1 if power == 1 else (-1 if power == p - 1 else power)

    result = {
        'a': a,
        'p': p,
        'a_mod_p': a_mod,
        'value': value,
        'power_computation': f"{a_mod}^({p}−1)/2 = {a_mod}^{(p-1)//2} mod {p} = {power}"
    }

    if value == 1:
        # Найдём корни
        sqrt_val = modular_sqrt(a_mod, p)
        result['meaning'] = f"{a} — квадратичный ВЫЧЕТ по модулю {p}"
        if sqrt_val is not None:
            result['sqrt'] = [sqrt_val, (p - sqrt_val) % p]
    elif value == -1:
        result['meaning'] = f"{a} — квадратичный НЕВЫЧЕТ по модулю {p}"

    return result


def jacobi_info(a, n):
    """Символ Якоби (a/n) для нечётного n."""
    a, n = int(a), int(n)
    if n % 2 == 0 or n < 1:
        return {'error': 'n должно быть положительным нечётным'}

    jac = int(jacobi_symbol(a, n))

    return {
        'a': a,
        'n': n,
        'value': jac,
        'factorization_of_n': {int(p): int(e) for p, e in factorint(n).items()},
        'meaning': ("1 — возможно квадратичный вычет" if jac == 1
                    else ("−1 — точно квадратичный НЕвычет" if jac == -1
                          else "0 — a и n не взаимно просты"))
    }


def check_b_smooth(n, B):
    """Проверяет, является ли n B-гладким, и если да — раскладывает."""
    n_orig, n = int(n), int(abs(int(n)))
    B = int(B)

    primes = list(primerange(2, B + 1))
    exps = {}
    x = n
    for p in primes:
        e = 0
        while x % p == 0:
            x //= p
            e += 1
        if e > 0:
            exps[p] = e

    is_smooth = (x == 1)
    return {
        'n': n_orig,
        'B': B,
        'primes_up_to_B': primes,
        'is_B_smooth': is_smooth,
        'factorization': exps,
        'remainder': x,  # что осталось после деления — должно быть 1 если гладкое
        'factorization_str': (' · '.join(f"{p}^{e}" if e > 1 else str(p)
                                         for p, e in sorted(exps.items())) +
                              (f" · {x}" if x > 1 else '')) if exps or x > 1 else '1'
    }


if __name__ == '__main__':
    # Self-test
    print("=== n = 180 ===")
    r = number_properties(180)
    print(f"  факторизация: {r['factorization_str']}")
    print(f"  φ(180) = {r['phi']}")
    print(f"  λ(180) = {r['lambda_carmichael']}")
    print(f"  делителей: {r['num_divisors']}")

    print("\n=== n = 181 (простое) ===")
    r = number_properties(181)
    print(f"  is_prime: {r['is_prime']}")
    print(f"  φ(181) = {r['phi']}")
    print(f"  первообразные корни: {r['primitive_roots']}")

    print("\n=== порядок 2 mod 181 ===")
    r = element_order(2, 181)
    print(f"  ord = {r['order']}")

    print("\n=== Символ Лежандра (5 / 181) ===")
    r = legendre_symbol_info(5, 181)
    print(f"  {r['value']} — {r['meaning']}")
    if 'sqrt' in r:
        print(f"  корни: {r['sqrt']}")

    print("\n=== 7350 является 11-гладким? ===")
    r = check_b_smooth(7350, 11)
    print(f"  {r['is_B_smooth']}  разложение: {r['factorization_str']}")
