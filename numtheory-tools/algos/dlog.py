"""
Алгоритмы дискретного логарифмирования по лекции 5 Полякова М.В.

Решают: найти x такой, что g^x ≡ h (mod p), где порядок g равен n.

Реализации буквально следуют формулировкам препода:
- BSGS (Гельфонда-Шэнкса), Л5 стр. 4
- ρ-алгоритм Полларда для ДЛ, Л5 стр. 18–20
- Алгоритм Полига-Хеллмана, Л5 стр. 9–10
"""
from sympy import factorint, isprime, gcd as sympy_gcd


# ────────────────────────────────────────────────────────────
# Вспомогательные штуки
# ────────────────────────────────────────────────────────────

def _order_of_g(g, p):
    """
    Порядок элемента g в F*_p. Работает, когда p — простое.
    Если g — первообразный корень, то порядок = p-1.
    Иначе ищем минимальный делитель d числа p-1 такой, что g^d ≡ 1.
    """
    if not isprime(p):
        # Для составных p порядок группы сложнее — считаем по Эйлеру;
        # для учебных задач всегда p простое.
        return p - 1
    n = p - 1
    # Разложение n — делители по убыванию (для поиска min d с g^d=1)
    factors = factorint(n)
    order = n
    for q in factors:
        while order % q == 0 and pow(g, order // q, p) == 1:
            order //= q
    return order


def _resolve_h(g, p, x=None, h=None):
    """
    Входные варианты:
      — даны g, p, x → вычисляем h = g^x mod p (формат препода)
      — даны g, p, h → работаем с h напрямую
    Возвращает (h, подсказка_строкой_или_None).
    """
    if h is None and x is None:
        raise ValueError("Нужно задать либо h, либо x")
    if h is None:
        h = pow(int(g), int(x), int(p))
        return h, f"Вычислили h = g^x mod p = {g}^{x} mod {p} = {h}"
    return int(h), None


# ────────────────────────────────────────────────────────────
# 1. BSGS (Гельфонда-Шэнкса)
# ────────────────────────────────────────────────────────────

def bsgs(g, p, x=None, h=None):
    """
    Алгоритм Гельфонда-Шэнкса по слайду Л5 стр. 4:

      Дано: |G| = n, y = g^x (mod p)
      Выход: x = log_g y

      1) s = ⌈√n⌉
      2) Малый шаг: S = {(y·g^r, r)} для r = 0..s-1
      3) Большой шаг: T = {(g^(ts), ts)} для t = 0..s
      4) Найти совпадение y·g^r = g^(ts) → x = ts - r

    Внутри используются обозначения препода: y = h.
    """
    steps = []
    g, p = int(g), int(p)
    h, hint = _resolve_h(g, p, x=x, h=h)
    if hint:
        steps.append({'iteration': 0, 'phase': 'init', 'description': hint})

    n = _order_of_g(g, p)
    # s = ⌈√n⌉
    s = int(n ** 0.5)
    if s * s < n:
        s += 1

    steps.append({
        'iteration': 1,
        'phase': 'init',
        'description': f"Дано: G = ⟨g⟩ ≤ F*_{p}, g = {g}, y = {h}, "
                       f"|G| = n = {n}. Ищем x = log_g y."
    })
    steps.append({
        'iteration': 2,
        'description': f"Вычисляем s = ⌈√n⌉ = ⌈√{n}⌉ = {s}"
    })

    # Малый шаг: S = {(y·g^r mod p, r) : r=0..s-1}
    S = {}   # значение → r (для поиска совпадения)
    S_list = []
    val = h % p
    for r in range(s):
        if val not in S:
            S[val] = r  # сохраняем первый попавшийся r
        S_list.append((val, r))
        val = (val * g) % p

    steps.append({
        'iteration': 3,
        'phase': 'small',
        'description': f"Малый шаг: S = {{(y·g^r mod p, r) : r = 0…{s-1}}} = "
                       f"{{{', '.join(f'({v}, {r})' for v, r in S_list)}}}"
    })

    # Большой шаг: T = {(g^(ts) mod p, ts) : t=0..s}
    gs = pow(g, s, p)          # = g^s mod p
    T_list = []
    cur = 1
    found = None
    for t in range(s + 1):
        ts = t * s
        T_list.append((cur, ts))
        if cur in S:
            r_match = S[cur]
            found = (r_match, ts, cur)
            # НЕ выходим, соберём всю таблицу для наглядности —
            # в примере препода видно всю T целиком
        cur = (cur * gs) % p

    steps.append({
        'iteration': 4,
        'phase': 'big',
        'description': f"Большой шаг: T = {{(g^(ts) mod p, ts) : t = 0…{s}}} = "
                       f"{{{', '.join(f'({v}, {ts})' for v, ts in T_list)}}}"
    })

    if found is None:
        steps.append({
            'iteration': 'Ошибка',
            'description': "Совпадение y·g^r = g^(ts) не найдено. "
                           "Проверь, что g — образующий подгруппы нужного порядка."
        })
        return None, steps

    r_match, ts, common = found
    x_found = (ts - r_match) % n

    steps.append({
        'iteration': 5,
        'description': f"Общий элемент: y·g^{r_match} = g^{ts} = {common} (mod {p})"
    })
    steps.append({
        'iteration': 'Результат',
        'description': f"x = ts − r = {ts} − {r_match} = {x_found}   "
                       f"(проверка: g^x = {g}^{x_found} mod {p} = {pow(g, x_found, p)})"
    })
    return x_found, steps


# ────────────────────────────────────────────────────────────
# 2. ρ-алгоритм Полларда для ДЛ
# ────────────────────────────────────────────────────────────

def rho_pollard_dlog(g, p, x=None, h=None):
    """
    ρ-алгоритм Полларда по псевдокоду Л5 стр. 20 и таблице-примеру стр. 21.

    Инвариант: z_i = g^{α_i} · y^{β_i} (mod p).
    Разбиение G = A ∪ B ∪ C по остатку z mod 3 — согласовано с примером
    со стр. 21 (z₀=1, α₀=0, β₀=0, z₁=228=y, α₁=0, β₁=1):
        z ≡ 1 (mod 3) → A:  f(z) = y·z,       β += 1, α не меняется
        z ≡ 0 (mod 3) → B:  f(z) = z²,        α *= 2, β *= 2
        z ≡ 2 (mod 3) → C:  f(z) = g·z,       α += 1, β не меняется

    (На слайде препода ветки "α+1 / β+1" и функции f в них подписаны
    слегка непоследовательно — реализация выбрана так, чтобы z_i
    совпадало с таблицей примера и сохранялся инвариант.)

    Итоговая формула: x_log = (α_{2i} − α_i) · (β_i − β_{2i})⁻¹ (mod n).
    """
    steps = []
    g, p = int(g), int(p)
    h, hint = _resolve_h(g, p, x=x, h=h)
    if hint:
        steps.append({'iteration': 0, 'phase': 'init', 'description': hint})

    n = _order_of_g(g, p)
    steps.append({
        'iteration': 1,
        'phase': 'init',
        'description': f"Дано: g = {g}, y = {h}, p = {p}, |G| = n = {n}. "
                       f"Разбиение G = A ∪ B ∪ C по остатку z mod 3."
    })

    def step_f(z, alpha, beta):
        """Один шаг: (z, α, β) → (f(z), α', β')."""
        r = z % 3
        if r == 1:
            # A: f(z) = y·z, β += 1
            new_z = (h * z) % p
            return new_z, alpha, (beta + 1) % n, 'A'
        elif r == 0:
            # B: f(z) = z², α↦2α, β↦2β
            new_z = (z * z) % p
            return new_z, (2 * alpha) % n, (2 * beta) % n, 'B'
        else:  # r == 2
            # C: f(z) = g·z, α += 1
            new_z = (g * z) % p
            return new_z, (alpha + 1) % n, beta, 'C'

    # Инициализация по псевдокоду: z₀ = 1, α₀ = 0, β₀ = 0
    z_i, alpha_i, beta_i = 1, 0, 0
    z_2i, alpha_2i, beta_2i = 1, 0, 0

    steps.append({
        'iteration': 2,
        'description': f"Инициализация: z₀ = 1, α₀ = 0, β₀ = 0"
    })

    # Таблица-строки (как на слайде стр. 21): i | z_i | α_i | β_i | z_{2i} | α_{2i} | β_{2i}
    table_header = {
        'iteration': 'table',
        'phase': 'table_header',
        'description': "i | z_i | α_i | β_i | z_{2i} | α_{2i} | β_{2i}"
    }
    steps.append(table_header)

    max_iter = 10 * n  # практический лимит
    for i in range(1, max_iter + 1):
        # Один шаг для медленного указателя
        z_i, alpha_i, beta_i, tag_i = step_f(z_i, alpha_i, beta_i)
        # Два шага для быстрого
        z_2i, alpha_2i, beta_2i, _ = step_f(z_2i, alpha_2i, beta_2i)
        z_2i, alpha_2i, beta_2i, _ = step_f(z_2i, alpha_2i, beta_2i)

        steps.append({
            'iteration': i,
            'phase': 'table_row',
            'description': f"{i} | {z_i} | {alpha_i} | {beta_i} | "
                           f"{z_2i} | {alpha_2i} | {beta_2i}"
        })

        if z_i == z_2i:
            # Нашли коллизию. Формула: g^(α_i) · y^(β_i) = g^(α_{2i}) · y^(β_{2i})
            # ⇒ y^(β_i − β_{2i}) = g^(α_{2i} − α_i)
            # ⇒ (β_i − β_{2i}) · x = (α_{2i} − α_i)  (mod n)
            delta_a = (alpha_2i - alpha_i) % n
            delta_b = (beta_i - beta_2i) % n
            steps.append({
                'iteration': 'collision',
                'description': f"Коллизия на шаге i = {i}: z_i = z_{{2i}} = {z_i}. "
                               f"Получаем (β_i − β_{{2i}})·x ≡ (α_{{2i}} − α_i) (mod n), "
                               f"то есть {delta_b}·x ≡ {delta_a} (mod {n})."
            })
            if delta_b == 0:
                if delta_a == 0:
                    steps.append({
                        'iteration': 'Ошибка',
                        'description': "0 ≡ 0 — тривиальное равенство, алгоритм не дал "
                                       "информации. Попробуй другое z₀ или другую функцию f."
                    })
                else:
                    steps.append({
                        'iteration': 'Ошибка',
                        'description': f"β_i − β_{{2i}} = 0, но α_{{2i}} − α_i = {delta_a} ≠ 0. "
                                       "Решения нет, алгоритм сработал неудачно — "
                                       "попробуй другое z₀."
                    })
                return None, steps

            # Решаем delta_b · x ≡ delta_a (mod n)
            d = int(sympy_gcd(delta_b, n))
            if delta_a % d != 0:
                steps.append({
                    'iteration': 'Ошибка',
                    'description': f"НОД({delta_b}, {n}) = {d} не делит {delta_a} → "
                                   "решений у сравнения нет."
                })
                return None, steps

            # Общее решение: x ≡ (delta_a/d) · (delta_b/d)^(−1) (mod n/d)
            nn = n // d
            a_ = delta_a // d
            b_ = delta_b // d
            inv = pow(b_, -1, nn)
            x_base = (a_ * inv) % nn

            # Перебираем d кандидатов по mod n
            candidates = [(x_base + k * nn) % n for k in range(d)]
            for cand in candidates:
                if pow(g, cand, p) == h % p:
                    steps.append({
                        'iteration': 'Результат',
                        'description': f"x ≡ (α_{{2i}} − α_i) · (β_i − β_{{2i}})⁻¹ "
                                       f"= {delta_a} · {delta_b}⁻¹ (mod {n}) = {cand}.  "
                                       f"Проверка: g^x = {g}^{cand} mod {p} = {pow(g, cand, p)}."
                    })
                    return cand, steps
            steps.append({
                'iteration': 'Ошибка',
                'description': "Ни один из кандидатов не подошёл при проверке g^x = h."
            })
            return None, steps

    steps.append({
        'iteration': 'Ошибка',
        'description': f"Коллизия не найдена за {max_iter} шагов."
    })
    return None, steps


# ────────────────────────────────────────────────────────────
# 3. Полиг-Хеллман
# ────────────────────────────────────────────────────────────

def _crt(congruences):
    """
    Китайская теорема об остатках.
    Вход: [(a1, m1), (a2, m2), ...]; модули попарно взаимно просты.
    Выход: (x, M), где x ≡ ai (mod mi), M = произведение модулей.
    """
    x, M = 0, 1
    for a, m in congruences:
        # x ≡ x (mod M), x ≡ a (mod m); ищем новое x mod (M*m)
        # новое x = x + M * k, где M*k ≡ (a − x) (mod m)
        diff = (a - x) % m
        inv = pow(M % m, -1, m)
        k = (diff * inv) % m
        x = x + M * k
        M = M * m
        x %= M
    return x, M


def pohlig_hellman(g, p, x=None, h=None):
    """
    Алгоритм Полига-Хеллмана по Л5 стр. 9–16.

    Работает для F*_p, p простое. Порядок n = p − 1 = ∏ p_i^{α_i}.
    Для каждого простого p_i:
      1) считаем таблицу r_{p_i, j} = g^{j·n/p_i} (mod p), j = 0..p_i-1
      2) поочерёдно находим x_0, x_1, …, x_{α_i − 1}
         — разряды x в p_i-ичном представлении (mod p_i^{α_i})
    Затем по КТО собираем общее x (mod n).

    Формат вывода буквально соответствует примеру стр. 13–16
    (log_2 62 mod 181, ответ x = 100).
    """
    steps = []
    g, p = int(g), int(p)
    h, hint = _resolve_h(g, p, x=x, h=h)
    if hint:
        steps.append({'iteration': 0, 'phase': 'init', 'description': hint})

    if not isprime(p):
        steps.append({
            'iteration': 'Ошибка',
            'description': f"p = {p} не простое. Версия Полига-Хеллмана реализована для F*_p."
        })
        return None, steps

    n = p - 1
    factors = factorint(n)   # {p_i: α_i}

    steps.append({
        'iteration': 1,
        'phase': 'init',
        'description': f"1. Раскладываем q = p − 1 = {n} = "
                       + " · ".join(f"{pi}^{ai}" if ai > 1 else f"{pi}"
                                    for pi, ai in factors.items())
    })

    congruences = []  # для КТО: (x mod p_i^α_i, p_i^α_i)

    step_counter = 2
    for pi, alpha in factors.items():
        # 2.1 Таблица r_{p_i, j}
        exp = n // pi
        table = []
        for j in range(pi):
            rj = pow(g, j * exp, p)
            table.append(rj)
        table_str = ",  ".join(
            f"r_{{{pi},{j}}} = {g}^({exp}·{j}) mod {p} = {table[j]}"
            for j in range(pi)
        )
        steps.append({
            'iteration': step_counter,
            'phase': 'table',
            'description': f"Для p_i = {pi}: таблица r_{{p_i, j}} = "
                           f"g^{{j·(p−1)/p_i}} (mod p):  {table_str}"
        })
        step_counter += 1

        # 2.2 Ищем цифры x_0, x_1, ..., x_{α-1} в p_i-ичной системе
        digits = []    # x = x_0 + x_1·pi + x_2·pi² + ...
        # Текущий "хвост" hg^{-accum}, где accum = x_0 + x_1·pi + ...
        # Для j-й цифры: (h · g^{-accum})^{(p-1)/pi^{j+1}} = g^{x_j·(p-1)/pi}

        accum = 0   # Σ x_k · pi^k
        pi_pow = 1  # pi^j; начинаем с j=0 → pi^0 = 1

        for j in range(alpha):
            # Основание в сравнении
            # для j=0: h^{(p-1)/pi} == r_{pi, x_0}
            # для j>0: (h · g^{-accum})^{(p-1)/pi^{j+1}} == r_{pi, x_j}
            pi_pow_next = pi_pow * pi  # pi^{j+1}
            exp_j = n // pi_pow_next   # (p-1)/pi^{j+1}

            # база: h * g^{-accum}
            if accum == 0:
                base = h % p
                base_str = f"y"
            else:
                ginv_accum = pow(g, -accum, p)
                base = (h * ginv_accum) % p
                base_str = f"(y · g^(−{accum}))"

            target = pow(base, exp_j, p)

            # Ищем в таблице
            xj = None
            for candidate, val in enumerate(table):
                if val == target:
                    xj = candidate
                    break

            if xj is None:
                steps.append({
                    'iteration': 'Ошибка',
                    'description': f"Не нашли x_{j} в таблице для p_i = {pi}: target = {target}"
                })
                return None, steps

            steps.append({
                'iteration': step_counter,
                'phase': 'digit',
                'description': f"Ищем x_{j} (цифра x в {pi}-ичной записи, "
                               f"разряд {j}): {base_str}^({n}/{pi_pow_next}) "
                               f"= {base_str}^{exp_j} mod {p} = {target} "
                               f"= r_{{{pi},{xj}}} ⇒ x_{j} = {xj}"
            })
            step_counter += 1

            digits.append(xj)
            accum += xj * pi_pow
            pi_pow = pi_pow_next

        x_mod = accum % (pi ** alpha)
        congruences.append((x_mod, pi ** alpha))
        poly_str = " + ".join(
            f"{d}·{pi}^{k}" if k > 0 else f"{d}"
            for k, d in enumerate(digits)
        )
        steps.append({
            'iteration': step_counter,
            'phase': 'partial',
            'description': f"Итого x ≡ {poly_str} = {x_mod} (mod {pi}^{alpha} = {pi**alpha})"
        })
        step_counter += 1

    # КТО
    x_final, M = _crt(congruences)
    system_str = "; ".join(f"x ≡ {a} (mod {m})" for a, m in congruences)
    steps.append({
        'iteration': step_counter,
        'phase': 'crt',
        'description': f"Система: {system_str}. По китайской теореме об остатках: x = {x_final}."
    })
    steps.append({
        'iteration': 'Результат',
        'description': f"x = {x_final}.  Проверка: g^x = {g}^{x_final} mod {p} "
                       f"= {pow(g, x_final, p)}  (ожидалось y = {h})"
    })
    return x_final, steps


# ────────────────────────────────────────────────────────────
# Self-test на примерах препода
# ────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("=== BSGS: пример Л5 стр. 7 ===")
    print("G = F*_31, g = 3, y = 24 ⇒ ожидаем x = 13")
    x_, _ = bsgs(3, 31, h=24)
    print(f"Получили x = {x_}  {'OK' if x_ == 13 else 'FAIL'}\n")

    print("=== Полиг-Хеллман: пример Л5 стр. 13–16 ===")
    print("2^x ≡ 62 (mod 181) ⇒ ожидаем x = 100")
    x_, _ = pohlig_hellman(2, 181, h=62)
    print(f"Получили x = {x_}  {'OK' if x_ == 100 else 'FAIL'}\n")

    print("=== ρ-Полларда: проверка корректности (не на примере препода) ===")
    # Пример простой: p=23, g=5 (первообразный корень), h=g^10 mod 23 = ?
    p_test = 23
    g_test = 5
    for x_true in [1, 2, 7, 10, 15, 19]:
        h_test = pow(g_test, x_true, p_test)
        x_, _ = rho_pollard_dlog(g_test, p_test, h=h_test)
        ok = (x_ is not None and pow(g_test, x_, p_test) == h_test)
        print(f"  log_{g_test} {h_test} (mod {p_test}): expect {x_true}, "
              f"got {x_}  {'OK' if ok else 'FAIL'}")
