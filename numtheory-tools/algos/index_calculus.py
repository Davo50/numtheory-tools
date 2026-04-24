"""
Метод исчисления индексов (Index Calculus) — алгоритм Адлемана
по Л5 стр. 24.

Идея. Для фиксированной факторной базы ФБ = {q_1, ..., q_k} (малых простых)
собираем соотношения вида g^{k_j} ≡ q_1^{a_{j,1}} · ... · q_k^{a_{j,k}} (mod p).
Каждое такое соотношение даёт уравнение в Z/(p-1)Z:
    k_j ≡ Σ a_{j,i} · log_g q_i   (mod p-1).
Набрав > k таких уравнений (система не должна быть вырожденной),
решаем её и получаем log_g q_i для всех i.

Затем для цели y: ищем случайное r такое, что y · g^r является B-гладким.
Разложим: y · g^r = Π q_i^{β_i} (mod p). Тогда
    log_g y + r ≡ Σ β_i · log_g q_i   (mod p-1)
    ⇒ x = log_g y ≡ Σ β_i · log_g q_i − r   (mod p-1).

Пример на Л5 стр. 25–27: G = Z*_229, g = 6, y = 13. Ответ: x = 117.
"""
from sympy import factorint, isprime, primerange, gcd as sympy_gcd
import random


def _factor_smooth(n, base):
    """
    Пытается разложить n по простым из base = [q_1, ..., q_k].
    Возвращает список степеней [e_1, ..., e_k], если n B-гладкое, иначе None.
    """
    exps = [0] * len(base)
    x = int(n)
    for i, q in enumerate(base):
        while x % q == 0:
            x //= q
            exps[i] += 1
    if x == 1:
        return exps
    return None


def _solve_linear_system_mod(A, b, m):
    """
    Решает A · z ≡ b (mod m) по каждому простому делителю m и собирает
    решение по китайской теореме об остатках. Это надёжнее, чем гаусс по
    составному модулю (где далеко не все ненулевые обратимы).

    Возвращает [z_1, ..., z_k] mod m, или None при несовместности.
    """
    if not A or not A[0]:
        return None
    k = len(A[0])

    # Разложим m на простые степени
    factors = factorint(m)   # {p_i: e_i}

    # Для каждой степени p^e решаем СЛАУ mod p^e, используя p-адический подъём
    partial_solutions = []
    for pi, ei in factors.items():
        mod = pi ** ei
        sol = _solve_linear_system_mod_prime_power(A, b, pi, ei)
        if sol is None:
            return None
        partial_solutions.append((sol, mod))

    # КТО покомпонентно
    result = []
    for j in range(k):
        vals = [(sol[j], mod) for sol, mod in partial_solutions]
        # Простейший CRT
        X, M = 0, 1
        for a, mi in vals:
            # X ≡ X (mod M), X ≡ a (mod mi) → новое X
            diff = (a - X) % mi
            inv = pow(M % mi, -1, mi)
            t = (diff * inv) % mi
            X = X + M * t
            M = M * mi
            X %= M
        result.append(X % m)
    return result


def _solve_linear_system_mod_prime_power(A, b, p, e):
    """
    СЛАУ mod p^e методом Гаусса. В кольце Z/p^eZ элемент обратим ⇔
    не делится на p. Выбираем главный элемент с минимальной степенью p.
    Для простоты: предполагаем e ≤ ... (работаем через p-adic, но для простых p²
    как в примере 228=4·3·19 — e ≤ 2 — этого достаточно).
    """
    mod = p ** e
    n = len(A)
    k = len(A[0])
    M = [[row[j] % mod for j in range(k)] + [b[i] % mod] for i, row in enumerate(A)]

    row = 0
    for col in range(k):
        # Найдём строку, где M[r][col] минимально делится на p (в идеале — не делится)
        best_row, best_val_p = -1, e + 1
        for r in range(row, n):
            v = M[r][col] % mod
            if v == 0:
                continue
            # Степень p в v
            vp = 0
            tmp = v
            while tmp % p == 0 and vp < e:
                tmp //= p
                vp += 1
            if vp < best_val_p:
                best_val_p, best_row = vp, r
                if vp == 0:
                    break

        if best_row == -1 or best_val_p >= e:
            # столбец не содержит обратимых → не можем выделить — пропускаем
            continue
        M[row], M[best_row] = M[best_row], M[row]

        # Если pivot делится на p^vp, нужно работать в Z/p^{e-vp}Z для этой строки.
        # Чтобы не усложнять, ограничимся случаем vp = 0 (pivot обратим).
        if best_val_p > 0:
            # пропускаем, для наших задач обычно найдётся обратимый pivot
            continue

        inv = pow(M[row][col], -1, mod)
        M[row] = [(v * inv) % mod for v in M[row]]
        for r in range(n):
            if r != row and M[r][col] % mod != 0:
                factor = M[r][col] % mod
                M[r] = [(M[r][j] - factor * M[row][j]) % mod for j in range(k + 1)]
        row += 1
        if row == n:
            break

    # Извлекаем решение, отдельно отслеживая свободные переменные
    z = [0] * k
    for i in range(min(k, row)):
        for j in range(k):
            if M[i][j] % mod == 1:
                ok = all(M[r][j] % mod == (1 if r == i else 0) for r in range(row))
                if ok:
                    z[j] = M[i][k] % mod
                    break

    return z


def index_calculus(g, p, x=None, h=None, base_bound=None,
                   max_relations_extra=3, max_tries=50000, seed=None):
    """
    Алгоритм Адлемана для вычисления x = log_g h (mod p).

    Вход:
      g, p — параметры задачи (p простое, g — образующий F*_p)
      x    — (опционально) задать степень; тогда h = g^x mod p
      h    — (опционально) значение
      base_bound — граница B для факторной базы; если None, выбирается
                   разумный маленький бoundary (например 13 для p≈200).
      max_relations_extra — сколько лишних соотношений собирать сверх |ФБ|
      max_tries — максимальное число случайных степеней k_j

    Выход: (x, steps).
    """
    if seed is not None:
        random.seed(seed)

    steps = []
    g, p = int(g), int(p)

    if not isprime(p):
        steps.append({'iteration': 'Ошибка',
                      'description': f'p = {p} не простое'})
        return None, steps

    # Разбор h/x
    if h is None and x is not None:
        h = pow(g, int(x), p)
        steps.append({
            'iteration': 0,
            'phase': 'init',
            'description': f"Вычислили h = g^x mod p = {g}^{x} mod {p} = {h}"
        })
    elif h is None:
        steps.append({'iteration': 'Ошибка',
                      'description': 'Задайте h или x_known'})
        return None, steps
    h = int(h)

    n = p - 1

    # Выбираем факторную базу
    if base_bound is None:
        # эвристика: ФБ из первых нескольких простых
        # для p ~ 200-1000 хватит B = 11-13
        if p < 50:
            base_bound = 5
        elif p < 300:
            base_bound = 11
        elif p < 10000:
            base_bound = 19
        else:
            base_bound = 31
    FB = list(primerange(2, base_bound + 1))
    k = len(FB)

    steps.append({
        'iteration': 1,
        'phase': 'init',
        'description': f"Дано: G = F*_{p}, g = {g}, y = {h}, n = p − 1 = {n}. "
                       f"Факторная база ℱ𝓑 = {{{', '.join(str(q) for q in FB)}}} "
                       f"(границa B = {base_bound})."
    })

    # ─── Этап 1: собираем соотношения g^{k_j} ≡ Π q_i^{a_{j,i}} (mod p) ───
    need = k + max_relations_extra
    A = []
    b = []
    used_k = set()
    tries = 0
    relations = []  # (k_j, [exponents])
    steps.append({
        'iteration': 2,
        'phase': 'phase',
        'description': f"Этап 1: ищем ≥ {need} соотношений "
                       f"g^{{k_j}} ≡ q_1^{{a_1}} · … · q_{{{k}}}^{{a_{{{k}}}}} (mod p)"
    })

    while len(relations) < need and tries < max_tries:
        tries += 1
        kj = random.randint(1, n - 1)
        if kj in used_k:
            continue
        val = pow(g, kj, p)
        exps = _factor_smooth(val, FB)
        if exps is None:
            continue
        used_k.add(kj)
        relations.append((kj, exps))
        factoring = " · ".join(f"{q}^{e}" for q, e in zip(FB, exps) if e > 0) or "1"
        steps.append({
            'iteration': f"rel {len(relations)}",
            'phase': 'relation',
            'description': f"g^{kj} = {g}^{kj} mod {p} = {val} = {factoring}"
        })

    if len(relations) < k:
        steps.append({
            'iteration': 'Ошибка',
            'description': f"Собрано только {len(relations)} соотношений из "
                           f"{k} необходимых за {tries} попыток — попробуйте "
                           f"больше max_tries или меньший base_bound."
        })
        return None, steps

    # Формируем СЛАУ: A[j] · z = b[j] mod n, где z_i = log_g q_i
    A = [exps for kj, exps in relations]
    b = [kj for kj, exps in relations]

    steps.append({
        'iteration': 3,
        'phase': 'phase',
        'description': f"Этап 2: решаем СЛАУ (по модулю n = {n}) "
                       f"на неизвестные log_g q_i для q_i ∈ ℱ𝓑."
    })

    # Решаем. Модуль n, возможно составной.
    z = _solve_linear_system_mod(A, b, n)
    if z is None:
        steps.append({
            'iteration': 'Ошибка',
            'description': f"Не удалось решить СЛАУ по mod {n} (нет обратимых "
                           f"главных элементов). Попробуйте собрать больше "
                           f"соотношений или другую факторную базу."
        })
        return None, steps

    # Проверим и покажем log_g q_i
    for qi, zi in zip(FB, z):
        check = pow(g, zi, p)
        ok = (check == qi)
        steps.append({
            'iteration': f"log {qi}",
            'phase': 'discrete_log',
            'description': f"log_g {qi} = {zi}   "
                           f"(проверка: {g}^{zi} mod {p} = {check} "
                           f"{'✓' if ok else '✗'})"
        })

    # ─── Этап 3: находим r такое, что y · g^r гладкое по ФБ ───
    steps.append({
        'iteration': 4,
        'phase': 'phase',
        'description': f"Этап 3: ищем случайное r такое, что y · g^r "
                       f"разлагается по ℱ𝓑."
    })

    tries2 = 0
    found_r, found_exps = None, None
    while tries2 < max_tries:
        tries2 += 1
        r = random.randint(0, n - 1)
        val = (h * pow(g, r, p)) % p
        exps = _factor_smooth(val, FB)
        if exps is not None:
            found_r, found_exps = r, exps
            factoring = " · ".join(f"{q}^{e}" for q, e in zip(FB, exps) if e > 0) or "1"
            steps.append({
                'iteration': f"r найдено",
                'phase': 'relation',
                'description': f"r = {r}:   y · g^r = {h}·{g}^{r} mod {p} = {val} = {factoring}"
            })
            break

    if found_r is None:
        steps.append({
            'iteration': 'Ошибка',
            'description': "Не удалось найти подходящее r за лимит попыток."
        })
        return None, steps

    # x = Σ β_i · log_g q_i − r  (mod n)
    x_found = sum(beta * zi for beta, zi in zip(found_exps, z)) - found_r
    x_found %= n

    # Проверка: g^x ?= h
    check_val = pow(g, x_found, p)
    ok = (check_val == h)
    steps.append({
        'iteration': 'Результат',
        'description': f"x = (Σ β_i · log_g q_i) − r = "
                       f"({' + '.join(f'{beta}·{zi}' for beta, zi in zip(found_exps, z))}) − "
                       f"{found_r} = {x_found} (mod {n}).   "
                       f"Проверка: g^x = {g}^{x_found} mod {p} = {check_val} "
                       f"{'✓' if ok else '✗'}"
    })

    return x_found, steps


if __name__ == '__main__':
    # Эталон Л5 стр. 25–27: G = Z*_229, g = 6, y = 13, ожидаем x = 117
    print("=== Index Calculus: пример Л5 стр. 25–27 (g=6, p=229, y=13) ===")
    random.seed(42)
    x, steps = index_calculus(6, 229, h=13, base_bound=11, seed=42)
    print(f"Найден x = {x},  6^{x} mod 229 = {pow(6, x, 229)} (ожидали 13)")
    assert pow(6, x, 229) == 13
    print("OK")

    print("\n=== Ещё одна проверка ===")
    random.seed(1)
    x, _ = index_calculus(2, 181, h=62)
    print(f"log_2 62 mod 181 = {x},  2^{x} mod 181 = {pow(2, x, 181)} (ожидали 62)")
    assert pow(2, x, 181) == 62
    print("OK")
