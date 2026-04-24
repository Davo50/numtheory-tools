"""
Эллиптические кривые над F_p по лекции 8 Полякова М.В.

Работаем с кривой в форме Вейерштрасса:
    E : Y² = X³ + a·X + b   над F_p

Точки задаём парой (x, y); бесконечность обозначаем как None (нейтральный O).

Формулы (Л8 стр. 40, 42):
  Сложение P + Q = (x₃, y₃), P ≠ ±Q:
    λ = (y_2 − y_1)/(x_2 − x_1)
    x₃ = λ² − x_1 − x_2
    y₃ = λ·(x_1 − x_3) − y_1
  Удвоение 2P:
    λ = (3·x_1² + a)/(2·y_1)
    x₃ = λ² − 2·x_1
    y₃ = λ·(x_1 − x_3) − y_1
  P + (−P) = O (где −P = (x_1, −y_1)).

Реализованные операции:
  — ec_add(P, Q, a, b, p)     — сложение точек
  — ec_double(P, a, b, p)     — удвоение
  — ec_scalar_mult(k, P, …)   — скалярное умножение (бинарный метод)
  — ec_order(P, a, b, p)      — порядок точки
  — ec_bsgs(P, Q, a, b, p, n) — ECDLP: найти k такое, что Q = k·P (BSGS)
"""

from sympy import sqrt as sympy_sqrt, isprime, gcd as sympy_gcd


# ────────────────────────────────────────────────────────────
# Низкоуровневая арифметика
# ────────────────────────────────────────────────────────────

def _modinv(x, p):
    """Обратный по простому модулю p."""
    return pow(int(x) % p, -1, p)


def _is_on_curve(P, a, b, p):
    """P лежит на кривой y² ≡ x³ + ax + b (mod p)?"""
    if P is None:
        return True
    x, y = P
    return (y * y - x * x * x - a * x - b) % p == 0


def ec_add(P, Q, a, b, p):
    """Сложение двух точек. Возвращает P+Q (или None = бесконечность)."""
    if P is None:
        return Q
    if Q is None:
        return P

    x1, y1 = P
    x2, y2 = Q

    if x1 % p == x2 % p:
        if (y1 + y2) % p == 0:
            return None  # P + (−P) = O
        # P == Q → удвоение
        return ec_double(P, a, b, p)

    # общий случай
    lam = ((y2 - y1) * _modinv(x2 - x1, p)) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def ec_double(P, a, b, p):
    """Удвоение точки. Возвращает 2P."""
    if P is None:
        return None
    x1, y1 = P
    if y1 % p == 0:
        return None  # касательная вертикальная → O

    lam = ((3 * x1 * x1 + a) * _modinv(2 * y1, p)) % p
    x3 = (lam * lam - 2 * x1) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def ec_scalar_mult(k, P, a, b, p, verbose=False):
    """
    k·P через бинарный метод (double-and-add).
    Если verbose=True, возвращает (результат, шаги), иначе просто результат.
    """
    steps = []
    k = int(k)
    if k == 0 or P is None:
        if verbose:
            steps.append({'iteration': 'Результат',
                          'description': f"0·P = O" if k == 0 else "k·O = O"})
            return None, steps
        return None

    negate = (k < 0)
    k = abs(k)

    bits = bin(k)[2:]
    if verbose:
        steps.append({
            'iteration': 0,
            'phase': 'init',
            'description': f"Вычисляем k·P при k = {k} (= {bits}₂), P = {P}."
        })

    R = None  # накапливаем результат
    for i, bit in enumerate(bits):
        # удвоение
        R = ec_double(R, a, b, p)
        desc = f"шаг {i+1}/{len(bits)}, бит '{bit}': удвоение → R = {R}"
        if bit == '1':
            R = ec_add(R, P, a, b, p)
            desc += f";   бит = 1 ⇒ R ← R + P = {R}"
        if verbose:
            steps.append({'iteration': i + 1, 'description': desc})

    if negate and R is not None:
        R = (R[0], (-R[1]) % p)

    if verbose:
        steps.append({
            'iteration': 'Результат',
            'description': f"k·P = {R}"
        })
        return R, steps
    return R


# ────────────────────────────────────────────────────────────
# Высокоуровневые задачи: порядок точки и ECDLP
# ────────────────────────────────────────────────────────────

def ec_order(P, a, b, p, limit=None):
    """
    Находит порядок точки P на кривой простым перебором до limit
    (если не указано — ищем до 4·√p + 1, что заведомо больше |E(F_p)|).
    """
    if P is None:
        return 1
    if limit is None:
        # Граница Хассе: |E(F_p)| ≤ p + 1 + 2·√p
        limit = p + 1 + 2 * int(p ** 0.5) + 2
    R = P
    for k in range(1, limit + 1):
        if R is None:
            return k
        R = ec_add(R, P, a, b, p)
    return None


def ec_bsgs(P, Q, a, b, p, n=None):
    """
    Baby-Step Giant-Step для ECDLP: найти k такое, что Q = k·P (mod p).
    n — порядок точки P (если известен). Если не указан, пытаемся найти.
    """
    steps = []

    if not _is_on_curve(P, a, b, p):
        steps.append({'iteration': 'Ошибка',
                      'description': f'Точка P = {P} не лежит на кривой y² = x³ + {a}x + {b} (mod {p})'})
        return None, steps
    if Q is not None and not _is_on_curve(Q, a, b, p):
        steps.append({'iteration': 'Ошибка',
                      'description': f'Точка Q = {Q} не лежит на кривой'})
        return None, steps

    if n is None:
        n = ec_order(P, a, b, p)
        if n is None:
            steps.append({'iteration': 'Ошибка',
                          'description': f'Не удалось определить порядок P'})
            return None, steps

    s = int(n ** 0.5)
    if s * s < n:
        s += 1

    steps.append({
        'iteration': 1,
        'phase': 'init',
        'description': f"ECDLP на кривой y² = x³ + {a}x + {b} (mod {p}), "
                       f"P = {P}, Q = {Q}, ord(P) = n = {n}. Ищем k: Q = k·P. "
                       f"s = ⌈√n⌉ = {s}."
    })

    # Baby steps: {Q + r·P : r = 0..s−1} — храним значение → r
    baby = {}
    cur = Q
    baby_list = []
    for r in range(s):
        key = None if cur is None else (cur[0] % p, cur[1] % p)
        if key not in baby:
            baby[key] = r
        baby_list.append((cur, r))
        cur = ec_add(cur, P, a, b, p)

    steps.append({
        'iteration': 2,
        'phase': 'small',
        'description': f"Малый шаг: S = {{(Q + r·P, r) : r = 0…{s-1}}} = "
                       f"{{{', '.join(f'({v}, {r})' for v, r in baby_list)}}}"
    })

    # Giant steps: {t·(s·P) : t = 0..s} — ищем совпадение с baby
    sP = ec_scalar_mult(s, P, a, b, p)
    steps.append({
        'iteration': 3,
        'description': f"s·P = {sP}"
    })

    giant_list = []
    cur = None  # 0·(sP) = O
    found = None
    for t in range(s + 1):
        key = None if cur is None else (cur[0] % p, cur[1] % p)
        giant_list.append((cur, t * s))
        if key in baby:
            r_match = baby[key]
            found = (r_match, t, cur)
        cur = ec_add(cur, sP, a, b, p)

    steps.append({
        'iteration': 4,
        'phase': 'big',
        'description': f"Большой шаг: T = {{(t·(s·P), ts) : t = 0…{s}}} = "
                       f"{{{', '.join(f'({v}, {ts})' for v, ts in giant_list)}}}"
    })

    if found is None:
        steps.append({
            'iteration': 'Ошибка',
            'description': "Совпадение не найдено. Проверьте, что Q ∈ ⟨P⟩."
        })
        return None, steps

    r_match, t, common = found
    k = (t * s - r_match) % n

    # Проверка
    kP = ec_scalar_mult(k, P, a, b, p)
    ok = (kP == Q)
    steps.append({
        'iteration': 'Результат',
        'description': f"Совпадение: (Q + {r_match}·P) = ({t}·s)·P = {common}. "
                       f"Значит k = t·s − r = {t}·{s} − {r_match} = {k}.   "
                       f"Проверка: {k}·P = {kP} {'✓' if ok else '✗'}"
    })
    return k, steps


# ────────────────────────────────────────────────────────────
# Self-test
# ────────────────────────────────────────────────────────────
if __name__ == '__main__':
    # Пример: кривая y² = x³ + 2x + 3 mod 97
    a, b, p = 2, 3, 97
    P = (3, 6)  # проверим, что на кривой: 36 ≡ 27 + 6 + 3 = 36 ✓
    print(f"=== Кривая y² = x³ + {a}x + {b} (mod {p}) ===")
    print(f"P = {P}, на кривой: {_is_on_curve(P, a, b, p)}")

    print("\n--- Арифметика ---")
    P2 = ec_double(P, a, b, p)
    print(f"2P = {P2}")
    P3 = ec_add(P, P2, a, b, p)
    print(f"3P = {P3}")
    P5 = ec_scalar_mult(5, P, a, b, p)
    print(f"5P = {P5}")

    ord_P = ec_order(P, a, b, p)
    print(f"\nord(P) = {ord_P}")

    print("\n--- ECDLP (BSGS) ---")
    for k_true in [3, 7, 17, 31]:
        Q = ec_scalar_mult(k_true, P, a, b, p)
        k_found, _ = ec_bsgs(P, Q, a, b, p, n=ord_P)
        # k_found может отличаться на кратное ord_P
        ok = (ec_scalar_mult(k_found, P, a, b, p) == Q) if k_found is not None else False
        print(f"  Q = {k_true}·P = {Q},  BSGS → k = {k_found}  {'OK' if ok else 'FAIL'}")
