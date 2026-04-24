"""
Алгоритмы Евклида по лекции 1 Полякова М.В.
— Классический (в 10-ричной записи, через деление с остатком)
— Бинарный (по псевдокоду Л1 стр.18)
— НОД нескольких чисел
— Расширенный (коэффициенты Безу) по Л1 стр.19–20
— Разложение в цепную дробь по Л1 стр.22–28
"""
from sympy import gcd as sympy_gcd


def euclid_classic(a, b):
    """
    Классический алгоритм Евклида по схеме со слайда:
        a = b·q₁ + r₁        (a, b) = (b, r₁)
        b = r₁·q₂ + r₂       (b, r₁) = (r₁, r₂)
        r₁ = r₂·q₃ + r₃      ...
        ...
        r_{n−1} = r_n · q_{n+1}   (последний r_n — это НОД)

    Возвращает (НОД, список шагов).
    """
    steps = []

    # Всегда работаем с неотрицательными; большее — первым
    a, b = abs(int(a)), abs(int(b))
    if a < b:
        a, b = b, a
        steps.append({
            'iteration': 0,
            'phase': 'init',
            'description': f"Поменяли местами: теперь a = {a}, b = {b} (a ≥ b)"
        })

    if b == 0:
        steps.append({
            'iteration': 'Результат',
            'description': f"b = 0, поэтому НОД(a, b) = a = {a}"
        })
        return a, steps

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Дано: a = {a}, b = {b}. Ищем НОД(a, b)."
    })

    # Основной цикл. Ведём пару (x, y), где x — делимое, y — делитель.
    x, y = a, b
    i = 0
    while y != 0:
        i += 1
        q = x // y
        r = x % y
        steps.append({
            'iteration': i,
            'description': f"{x} = {y} · {q} + {r}   →   НОД({x}, {y}) = НОД({y}, {r})"
        })
        x, y = y, r

    steps.append({
        'iteration': 'Результат',
        'description': f"Последний ненулевой остаток = {x}. Значит НОД(a, b) = {x}."
    })
    return x, steps


def euclid_multi(numbers):
    """
    НОД нескольких чисел через цепочку:
    НОД(a, b, c, ...) = НОД(НОД(a, b), c, ...)
    """
    steps = []
    nums = [abs(int(n)) for n in numbers if int(n) != 0]
    if not nums:
        steps.append({'iteration': 'Ошибка',
                      'description': 'Все числа нулевые'})
        return None, steps

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Дано чисел: {len(nums)}. Считаем НОД последовательно, "
                       f"используя НОД(a₁, …, aₖ) = НОД(НОД(a₁, …, a_{{k−1}}), aₖ)."
    })

    cur = nums[0]
    steps.append({
        'iteration': 1,
        'description': f"Стартовое значение: g = a₁ = {cur}"
    })

    for idx, n in enumerate(nums[1:], start=2):
        prev = cur
        g, sub_steps = euclid_classic(cur, n)
        cur = g
        # Вложенно выводим шаги текущего парного Евклида
        steps.append({
            'iteration': f"{idx}",
            'phase': 'pair',
            'description': f"Шаг {idx - 1}: считаем НОД(g, a_{idx}) = НОД({prev}, {n}):"
        })
        for s in sub_steps:
            s2 = dict(s)
            s2['iteration'] = f"{idx}.{s['iteration']}"
            steps.append(s2)
        steps.append({
            'iteration': f"{idx}.итог",
            'description': f"g := НОД({prev}, {n}) = {cur}"
        })

    steps.append({
        'iteration': 'Результат',
        'description': f"НОД всех {len(nums)} чисел = {cur}"
    })
    return cur, steps


def euclid_binary(a, b):
    """
    Бинарный алгоритм Евклида по псевдокоду из Л1 (стр. 18):
        1: пока t > 0:
        2:    t ← |a − b|
        3:    пока a mod 2 = 0:
        4:       a ← a/2
        5:    если a ≥ b:
        6:       a ← t; goto 1
        7:    иначе:
        8:       b ← t; goto 1
        9: вернуть a

    Замечание: препод даёт упрощённую версию без вынесения общей степени 2.
    Реализация следует слайду буквально.
    """
    steps = []
    a, b = abs(int(a)), abs(int(b))
    a_orig, b_orig = a, b  # запомним для итоговой сверки

    if a == 0:
        steps.append({'iteration': 'Результат',
                      'description': f"a = 0, НОД = b = {b}"})
        return b, steps
    if b == 0:
        steps.append({'iteration': 'Результат',
                      'description': f"b = 0, НОД = a = {a}"})
        return a, steps

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Вход: a = {a}, b = {b}. Следуем псевдокоду из лекции."
    })

    # t инициализируем так, чтобы войти в цикл
    t = abs(a - b)
    i = 0
    while t > 0:
        i += 1
        t = abs(a - b)
        steps.append({
            'iteration': f"{i}.t",
            'description': f"t = |a − b| = |{a} − {b}| = {t}"
        })
        if t == 0:
            # Из цикла выйдем по условию while, но зафиксируем:
            steps.append({
                'iteration': f"{i}.eq",
                'description': f"a = b ⇒ выход из внешнего цикла"
            })
            break

        # Сокращение: пока a чётное, делим на 2
        if a % 2 == 0:
            reductions = []
            while a % 2 == 0:
                reductions.append(a // 2)
                a //= 2
            steps.append({
                'iteration': f"{i}.half",
                'description': f"a — чётное, делим на 2 пока можно: "
                               f"a → {' → '.join(str(x) for x in reductions)}"
            })

        # Сравнение и присваивание t
        if a >= b:
            steps.append({
                'iteration': f"{i}.cmp",
                'description': f"a = {a} ≥ b = {b} ⇒ a ← t = {t}"
            })
            a = t
        else:
            steps.append({
                'iteration': f"{i}.cmp",
                'description': f"a = {a} < b = {b} ⇒ b ← t = {t}"
            })
            b = t

        if i > 10000:
            steps.append({'iteration': 'Ошибка',
                          'description': 'Превышено число итераций'})
            return None, steps

    result = a
    steps.append({
        'iteration': 'Результат',
        'description': f"Возвращаем a = {result}. Значит НОД = {result} "
                       f"(по псевдокоду из Л1 стр. 18)."
    })
    # Сверка с истинным НОД — псевдокод со слайда может давать НЕВЕРНЫЙ ответ
    # на паре чётных чисел, т.к. в нём делится только a на 2 (общая степень 2
    # не выносится). Честно об этом предупреждаем.
    true_gcd = int(sympy_gcd(a_orig, b_orig))
    if result != true_gcd:
        steps.append({
            'iteration': '⚠',
            'description': f"Внимание: истинный НОД({a_orig}, {b_orig}) = {true_gcd}, "
                           f"но псевдокод препода даёт {result}. "
                           f"Причина: в слайде делится только a на 2, "
                           f"поэтому общая степень 2 у a и b может теряться. "
                           f"Если на контрольной требуется верный ответ — "
                           f"используй классический (10-ричный) Евклид."
        })
    return result, steps


# ────────────────────────────────────────────────────────────
# Расширенный алгоритм Евклида (коэффициенты Безу) — Л1 стр. 19–20
# ────────────────────────────────────────────────────────────

def euclid_extended(a, b):
    """
    Расширенный алгоритм Евклида: находит d = НОД(a, b) и такие x, y,
    что a·x + b·y = d.

    Следуем схеме со слайда Л1 стр. 19: ведём тройки (xᵢ, yᵢ, zᵢ) такие,
    что a·xᵢ + b·yᵢ = zᵢ. Начальные значения:
        (x₋₁, y₋₁, z₋₁) = (1, 0, a),
        (x₀,  y₀,  z₀ ) = (0, 1, b).
    Рекуррентность:
        (xᵢ, yᵢ, zᵢ) = (xᵢ₋₂, yᵢ₋₂, zᵢ₋₂) − qᵢ·(xᵢ₋₁, yᵢ₋₁, zᵢ₋₁),
        где qᵢ = ⌊zᵢ₋₂ / zᵢ₋₁⌋.
    Процесс останавливается, когда очередное zᵢ становится равно 0.
    Тогда предыдущее zᵢ₋₁ = НОД(a, b), а xᵢ₋₁, yᵢ₋₁ — коэффициенты Безу.
    """
    steps = []
    a_orig, b_orig = int(a), int(b)
    a, b = abs(int(a)), abs(int(b))

    if a == 0 and b == 0:
        steps.append({'iteration': 'Ошибка',
                      'description': 'a и b не могут быть оба равны 0'})
        return None, steps

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Вход: a = {a}, b = {b}. Ищем d = НОД(a,b) "
                       f"и x, y такие, что a·x + b·y = d."
    })

    # Инициализация троек
    x_prev, y_prev, z_prev = 1, 0, a
    x_cur,  y_cur,  z_cur  = 0, 1, b

    steps.append({
        'iteration': 'i=−1',
        'phase': 'row',
        'description': f"(x₋₁, y₋₁, z₋₁) = (1, 0, {a})"
    })
    steps.append({
        'iteration': 'i=0',
        'phase': 'row',
        'description': f"(x₀,  y₀,  z₀ ) = (0, 1, {b})"
    })

    i = 0
    while z_cur != 0:
        i += 1
        q = z_prev // z_cur
        x_new = x_prev - q * x_cur
        y_new = y_prev - q * y_cur
        z_new = z_prev - q * z_cur

        steps.append({
            'iteration': f"i={i}",
            'phase': 'row',
            'description': f"q_{i} = ⌊{z_prev}/{z_cur}⌋ = {q};   "
                           f"(x_{i}, y_{i}, z_{i}) = ({x_prev}, {y_prev}, {z_prev}) − "
                           f"{q}·({x_cur}, {y_cur}, {z_cur}) = ({x_new}, {y_new}, {z_new})"
        })

        x_prev, y_prev, z_prev = x_cur, y_cur, z_cur
        x_cur, y_cur, z_cur = x_new, y_new, z_new

        if i > 10000:
            steps.append({'iteration': 'Ошибка',
                          'description': 'Слишком много итераций'})
            return None, steps

    # Последняя ненулевая z — это НОД, соответствующие x, y — коэффициенты Безу
    d, x, y = z_prev, x_prev, y_prev

    # Если исходные a или b были отрицательны — возвращаем знак коэффициентам
    if a_orig < 0:
        x = -x
    if b_orig < 0:
        y = -y

    check = a_orig * x + b_orig * y
    steps.append({
        'iteration': 'Результат',
        'description': f"d = НОД({a_orig}, {b_orig}) = {d};   "
                       f"x = {x}, y = {y}.   "
                       f"Проверка: {a_orig}·({x}) + {b_orig}·({y}) = {check} "
                       f"{'✓' if check == d else '✗'}"
    })
    return (d, x, y), steps


# ────────────────────────────────────────────────────────────
# Цепные (непрерывные) дроби — Л1 стр. 22–28
# ────────────────────────────────────────────────────────────

def continued_fraction_rational(numerator, denominator):
    """
    Разложение рационального числа num/den в конечную цепную дробь
    num/den = q₁ + 1/(q₂ + 1/(q₃ + ...))

    Используем алгоритм Евклида: частные q_i — это частные делений.
    Одновременно строим таблицу подходящих дробей P_s/Q_s по рекуррентности
    со слайда Л1 стр. 28:
        P_s = q_s · P_{s-1} + P_{s-2},   P_{-1} = 1, P_0 = q_1
        Q_s = q_s · Q_{s-1} + Q_{s-2},   Q_{-1} = 0, Q_0 = 1

    Возвращает ((quotients, convergents), steps), где
        quotients   — [q₁, q₂, …, q_n]
        convergents — [(P_0, Q_0), (P_1, Q_1), …, (P_n, Q_n)]
    """
    steps = []
    num, den = int(numerator), int(denominator)
    if den == 0:
        steps.append({'iteration': 'Ошибка',
                      'description': 'Знаменатель не может быть 0'})
        return None, steps

    sign = -1 if (num < 0) ^ (den < 0) else 1
    num, den = abs(num), abs(den)
    orig_num, orig_den = num, den

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Раскладываем дробь {'-' if sign < 0 else ''}{orig_num}/{orig_den} "
                       f"в цепную дробь через алгоритм Евклида."
    })

    # Собираем частные q_s
    quotients = []
    a, b = num, den
    i = 0
    while b != 0:
        i += 1
        q = a // b
        r = a % b
        quotients.append(q)
        steps.append({
            'iteration': f"q_{i}",
            'phase': 'step',
            'description': f"{a} = {b}·{q} + {r}   ⇒   q_{i} = {q}"
        })
        a, b = b, r

    # Таблица подходящих дробей. По слайду P_{-1}=1, P_0=q_1, Q_{-1}=0, Q_0=1,
    # далее P_s = q_{s+1}·P_{s-1} + P_{s-2}.
    # Для простоты занумеруем подходящие с 0 и будем использовать квадратный
    # рекуррент с двумя начальными:
    P_minus1, Q_minus1 = 1, 0
    P0, Q0 = quotients[0], 1
    convergents = [(P0, Q0)]
    table_rows = [f"P₀/Q₀ = {P0}/{Q0}  (= q₁ = {quotients[0]})"]

    P_prev2, Q_prev2 = P_minus1, Q_minus1
    P_prev, Q_prev = P0, Q0
    for s in range(1, len(quotients)):
        q = quotients[s]
        P_s = q * P_prev + P_prev2
        Q_s = q * Q_prev + Q_prev2
        convergents.append((P_s, Q_s))
        table_rows.append(f"P_{s}/Q_{s} = {q}·{P_prev} + {P_prev2} / {q}·{Q_prev} + {Q_prev2} "
                          f"= {P_s}/{Q_s}")
        P_prev2, Q_prev2 = P_prev, Q_prev
        P_prev, Q_prev = P_s, Q_s

    steps.append({
        'iteration': 'CF',
        'phase': 'result',
        'description': f"Цепная дробь: [{quotients[0]}; " +
                       (", ".join(str(q) for q in quotients[1:]) if len(quotients) > 1 else "") +
                       f"] = {'+'.join(['q_'+str(i+1) for i in range(len(quotients))])} вложенно"
    })

    for row in table_rows:
        steps.append({
            'iteration': 'conv',
            'phase': 'convergent',
            'description': row
        })

    # Проверка: последняя подходящая дробь должна равняться исходной
    last_P, last_Q = convergents[-1]
    check_ok = (last_P * orig_den == last_Q * orig_num)
    steps.append({
        'iteration': 'Результат',
        'description': f"{'-' if sign < 0 else ''}{orig_num}/{orig_den} = "
                       f"[{quotients[0]}; {', '.join(str(q) for q in quotients[1:])}]. "
                       f"Последняя подходящая дробь = {last_P}/{last_Q} "
                       f"{'✓' if check_ok else '✗'}"
    })

    return {'quotients': quotients, 'convergents': convergents,
            'sign': sign, 'original': (orig_num, orig_den)}, steps


# Быстрый self-test при прямом запуске
if __name__ == '__main__':
    for x, y in [(1071, 462), (252, 198), (17, 5), (100, 75)]:
        r, _ = euclid_classic(x, y)
        r2, _ = euclid_binary(x, y)
        expect = int(sympy_gcd(x, y))
        print(f"НОД({x},{y}) = classic {r}, binary {r2}, sympy {expect}  "
              f"{'OK' if r == expect == r2 else 'FAIL'}")
    r, _ = euclid_multi([120, 36, 60, 84, 24])
    print(f"НОД(120,36,60,84,24) = {r} (ожидается 12)")

    print("\n=== Расширенный Евклид ===")
    for a, b in [(1071, 462), (99, 78), (17, 5)]:
        (d, x, y), _ = euclid_extended(a, b)
        ok = (a * x + b * y == d) and (d == int(sympy_gcd(a, b)))
        print(f"  extended({a}, {b}) = d={d}, x={x}, y={y}  "
              f"({a}·{x} + {b}·{y} = {a*x+b*y})  {'OK' if ok else 'FAIL'}")

    print("\n=== Цепные дроби ===")
    data, _ = continued_fraction_rational(1071, 462)
    print(f"  1071/462 = [{data['quotients']}]")
    print(f"  Подходящие: {data['convergents']}")
    # 1071/462 = 2 + 147/462; 462/147=3 ост 21; 147/21=7 ост 0
    # Ожидаем [2, 3, 7]  (поскольку 147 = 3·49 не совсем, но как Евклид)
    # На самом деле 1071 = 462·2 + 147; 462 = 147·3 + 21; 147 = 21·7 + 0
    # → [2; 3, 7]
    assert data['quotients'] == [2, 3, 7], f"Ожидалось [2,3,7], получили {data['quotients']}"
    print("  ✓ OK")
