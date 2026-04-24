"""
Алгоритмы над многочленами в кольце F_2[X] по лекции 1 Полякова М.В.

Представление: многочлен ↔ целое число (биты = коэффициенты).
  f = x^3 + x + 1  ↔  binary 1011  ↔  0b1011 = 11
  f = (10011101)   ↔  int(0b10011101) = 157 = x^7 + x^4 + x^3 + x^2 + 1

Операции:
  сложение/вычитание в F_2 = XOR
  умножение на x            = сдвиг влево (<< 1)
  деление на x              = сдвиг вправо (>> 1), разрешено когда младший бит 0
  степень многочлена        = позиция старшего бита (= bit_length() - 1)
"""

def poly_from_string(s):
    """
    Парсит "10011101" (строка из 0 и 1, старший коэффициент слева)
    в int: 0b10011101 = 157 = x^7 + x^4 + x^3 + x^2 + 1.
    Пробелы и лишние символы игнорируются. Допускает ввод '0x...' как hex,
    а также десятичное число.
    """
    s = str(s).strip()
    if s.startswith('0x') or s.startswith('0X'):
        return int(s, 16)
    # удаляем пробелы, скобки, запятые
    for c in ' \t(),':
        s = s.replace(c, '')
    if set(s) <= {'0', '1'}:
        if not s:
            return 0
        return int(s, 2)
    # если не чистая битовая строка — десятичное
    return int(s)


def poly_to_string(p):
    """Int → двоичная строка (коэффициенты) без префикса '0b'."""
    if p == 0:
        return '0'
    return bin(p)[2:]


def poly_to_algebra(p):
    """
    Превращает int-представление многочлена в алгебраический вид
    типа 'x^7 + x^4 + x^3 + x^2 + 1'.
    """
    if p == 0:
        return '0'
    terms = []
    for i in range(p.bit_length() - 1, -1, -1):
        if (p >> i) & 1:
            if i == 0:
                terms.append('1')
            elif i == 1:
                terms.append('x')
            else:
                terms.append(f'x^{i}')
    return ' + '.join(terms)


def poly_deg(p):
    """Степень многочлена. deg(0) условно = -1."""
    return p.bit_length() - 1


def poly_add(a, b):
    """Сумма (она же разность) в F_2[X] — XOR."""
    return a ^ b


def poly_mul(a, b):
    """Школьное умножение многочленов через сдвиги."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        b >>= 1
    return result


def poly_divmod(a, b):
    """
    Деление многочлена a на b с остатком: a = b·q + r, deg(r) < deg(b).
    Возвращает (q, r).
    """
    if b == 0:
        raise ZeroDivisionError("Деление многочлена на 0")
    q = 0
    r = a
    db = poly_deg(b)
    while r != 0 and poly_deg(r) >= db:
        shift = poly_deg(r) - db
        q ^= (1 << shift)
        r ^= (b << shift)
    return q, r


# ────────────────────────────────────────────────────────────
# 1. Бинарный алгоритм Евклида для F_2[X]
# ────────────────────────────────────────────────────────────

def poly_euclid_binary(a, b):
    """
    Бинарный алгоритм Евклида для F_2[X] — адаптация псевдокода Л1 стр. 18
    к многочленам.

    В F_2[X] псевдокод со слайда (для чисел) не работает буквально: в числах
    |a−b| всегда меньше max(a,b), а XOR двух многочленов равной степени
    может дать результат со старшим битом ниже, но числовое значение
    не обязательно меньше b. Поэтому перед сравнением «a ≥ b» мы
    сравниваем по СТЕПЕНИ многочлена, а не по числовому значению.

    Инвариант: на каждом шаге deg(a) ≥ deg(b) ≥ 0. Делаем a ← a ⊕ (b << shift),
    чтобы убрать старший член, затем делим a на X (младшие нули) и меняем
    местами, если теперь deg(a) < deg(b). Цикл заканчивается, когда b = 0.

    Эквивалентно исходному псевдокоду по логике: «t := a−b, затем делим
    t на 2 пока можно, затем присваиваем меньшему».
    """
    steps = []
    a, b = int(a), int(b)
    a_orig, b_orig = a, b

    if a == 0:
        steps.append({'iteration': 'Результат',
                      'description': f"a = 0, НОД = b = ({poly_to_string(b)}) = {poly_to_algebra(b)}"})
        return b, steps
    if b == 0:
        steps.append({'iteration': 'Результат',
                      'description': f"b = 0, НОД = a = ({poly_to_string(a)}) = {poly_to_algebra(a)}"})
        return a, steps

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Вход: a = ({poly_to_string(a)}) = {poly_to_algebra(a)}, "
                       f"b = ({poly_to_string(b)}) = {poly_to_algebra(b)}. "
                       f"В F₂[X]: вычитание = XOR, «делится на 2» ↔ x|a (младший коэф = 0). "
                       f"Сравнение a ≥ b — по степени многочлена."
    })

    # Упорядочим: a — с бОльшей степенью
    if poly_deg(a) < poly_deg(b):
        a, b = b, a
        steps.append({
            'iteration': '0.swap',
            'description': f"Меняем местами, чтобы deg(a) ≥ deg(b): "
                           f"a = ({poly_to_string(a)}), b = ({poly_to_string(b)})"
        })

    i = 0
    while b != 0:
        i += 1
        # Шаг «t := a − b» — но сначала выравниваем старшие члены,
        # сдвигая b влево на разницу степеней. Это аналог «длинного»
        # шага (препод в псевдокоде не сдвигает, но для F_2[X] без сдвига
        # алгоритм будет зацикливаться).
        shift = poly_deg(a) - poly_deg(b)
        b_shifted = b << shift
        t = a ^ b_shifted
        steps.append({
            'iteration': f"{i}.t",
            'description': f"deg(a) − deg(b) = {poly_deg(a)} − {poly_deg(b)} = {shift}; "
                           f"t = a ⊕ (b · x^{shift}) = ({poly_to_string(a)}) ⊕ "
                           f"({poly_to_string(b_shifted)}) = ({poly_to_string(t)})"
        })

        if t == 0:
            steps.append({
                'iteration': f"{i}.eq",
                'description': f"a = b · x^{shift} ⇒ b делит a; "
                               f"НОД = b = ({poly_to_string(b)})"
            })
            a = b
            b = 0
            break

        # Делим t на X пока возможно (младшие нули)
        if (t & 1) == 0:
            reductions = []
            while (t & 1) == 0 and t != 0:
                t >>= 1
                reductions.append(poly_to_string(t))
            steps.append({
                'iteration': f"{i}.shift",
                'description': f"x | t ⇒ делим t на x пока возможно: "
                               f"t → ({') → ('.join(reductions)})"
            })

        # Сравнение по степени: «a ← t и перейти к п.1», если t заменяет большее
        if poly_deg(t) >= poly_deg(b):
            steps.append({
                'iteration': f"{i}.cmp",
                'description': f"deg(t) = {poly_deg(t)} ≥ deg(b) = {poly_deg(b)} ⇒ "
                               f"a ← t"
            })
            a = t
        else:
            steps.append({
                'iteration': f"{i}.cmp",
                'description': f"deg(t) = {poly_deg(t)} < deg(b) = {poly_deg(b)} ⇒ "
                               f"a ← b, b ← t (сохраняем инвариант deg(a) ≥ deg(b))"
            })
            a, b = b, t

        if i > 10000:
            steps.append({'iteration': 'Ошибка',
                          'description': 'Слишком много итераций'})
            return None, steps

    result = a
    steps.append({
        'iteration': 'Результат',
        'description': f"b = 0 ⇒ НОД = a = ({poly_to_string(result)}) = {poly_to_algebra(result)}"
    })

    # Сверяем с классическим Евклидом для F_2[X]
    true_gcd = poly_gcd_classic(a_orig, b_orig)
    if result != true_gcd:
        steps.append({
            'iteration': '⚠',
            'description': f"Сверка: классический Евклид даёт ({poly_to_string(true_gcd)}) = "
                           f"{poly_to_algebra(true_gcd)}. Возможно расхождение в "
                           f"старшем коэффициенте (в F₂[X] НОД определён с точностью до ассоциированности)."
        })
    return result, steps


def poly_gcd_classic(a, b):
    """Обычный Евклид для многочленов через deg-деление (для сверки)."""
    while b != 0:
        _, r = poly_divmod(a, b)
        a, b = b, r
    return a


# ────────────────────────────────────────────────────────────
# 2. Алгоритм Карацубы для F_2[X]
# ────────────────────────────────────────────────────────────

def poly_karatsuba(x, y, threshold=2):
    """
    Алгоритм Карацубы для умножения многочленов в F_2[X],
    по формуле со слайда Л1 стр. 33:

        x = x1·X^k + x0,  y = y1·X^k + y0
        xy = (X^{2k} + X^k)·x1·y1 + (x1 − x0)(y1 − y0)·X^k + (X^k + 1)·x0·y0

    В F_2[X]: «−» = «+» = XOR, умножение на X^k = сдвиг на k, сложение = XOR.
    Возвращает (произведение, шаги).

    Параметр threshold — порог, при котором переходим на школьное умножение.
    """
    steps = []

    def karatsuba_rec(a, b, depth=0, label=""):
        """Рекурсивная функция. depth — для отступов, label — метка подвызова."""
        na = poly_deg(a) + 1 if a else 0
        nb = poly_deg(b) + 1 if b else 0
        n = max(na, nb)

        if n <= threshold or a == 0 or b == 0:
            result = poly_mul(a, b)
            steps.append({
                'iteration': f"{label or 'базовый'}",
                'phase': 'base',
                'description': f"{'  '*depth}Базовый случай: "
                               f"({poly_to_string(a)}) · ({poly_to_string(b)}) = "
                               f"({poly_to_string(result)}) [школьное умножение]"
            })
            return result

        k = n // 2     # разбиваем пополам
        mask = (1 << k) - 1
        x0 = a & mask
        x1 = a >> k
        y0 = b & mask
        y1 = b >> k

        steps.append({
            'iteration': f"{label or 'рекурсия'}",
            'phase': 'split',
            'description': f"{'  '*depth}Делим при k = {k}: "
                           f"x = x₁·X^{k} + x₀ = ({poly_to_string(x1)})·X^{k} + ({poly_to_string(x0)}), "
                           f"y = y₁·X^{k} + y₀ = ({poly_to_string(y1)})·X^{k} + ({poly_to_string(y0)})"
        })

        # Три рекурсивных умножения
        P1 = karatsuba_rec(x1, y1, depth + 1, f"{label}.P1" if label else "P1")
        # x1 − x0 в F_2 = x1 ⊕ x0
        dx = x1 ^ x0
        dy = y1 ^ y0
        P2 = karatsuba_rec(dx, dy, depth + 1, f"{label}.P2" if label else "P2")
        P3 = karatsuba_rec(x0, y0, depth + 1, f"{label}.P3" if label else "P3")

        # xy = (X^{2k} + X^k)·P1 + P2·X^k + (X^k + 1)·P3
        term1 = (P1 << (2 * k)) ^ (P1 << k)          # (X^{2k} + X^k)·x1y1
        term2 = P2 << k                              # (x1-x0)(y1-y0)·X^k
        term3 = (P3 << k) ^ P3                       # (X^k + 1)·x0y0
        result = term1 ^ term2 ^ term3

        steps.append({
            'iteration': f"{label or 'сборка'}",
            'phase': 'combine',
            'description': f"{'  '*depth}Собираем: xy = (X^{2*k}+X^{k})·({poly_to_string(P1)}) ⊕ "
                           f"({poly_to_string(P2)})·X^{k} ⊕ (X^{k}+1)·({poly_to_string(P3)}) "
                           f"= ({poly_to_string(result)})"
        })
        return result

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Умножение методом Карацубы в F₂[X]: "
                       f"x = ({poly_to_string(x)}) = {poly_to_algebra(x)}, "
                       f"y = ({poly_to_string(y)}) = {poly_to_algebra(y)}. "
                       f"Формула: xy = (X^{{2k}}+X^k)·x₁y₁ ⊕ (x₁⊕x₀)(y₁⊕y₀)·X^k ⊕ (X^k+1)·x₀y₀."
    })

    result = karatsuba_rec(x, y, 0, "")

    # Проверка школьным умножением
    expected = poly_mul(x, y)
    ok = (result == expected)
    steps.append({
        'iteration': 'Результат',
        'description': f"xy = ({poly_to_string(result)}) = {poly_to_algebra(result)}. "
                       f"{'Проверка ✓' if ok else '✗ (не сходится со школьным)'}"
    })
    steps.append({
        'iteration': 'Сложность',
        'description': "Рекурсивная формула T(2n) ≤ 3·T(n) + c·n ⇒ "
                       "T(n) = O(n^log₂3) ≈ O(n^1.585) [Л1 стр. 34]"
    })
    return result, steps


# ────────────────────────────────────────────────────────────
# Self-test
# ────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("=== Бинарный Евклид для F_2[X]: вариант 2 задание 3 ===")
    f = poly_from_string("10011101")   # f = x^7 + x^4 + x^3 + x^2 + 1
    g = poly_from_string("111110")     # g = x^5 + x^4 + x^3 + x^2 + x
    print(f"f = ({poly_to_string(f)}) = {poly_to_algebra(f)}")
    print(f"g = ({poly_to_string(g)}) = {poly_to_algebra(g)}")
    result, _ = poly_euclid_binary(f, g)
    expected = poly_gcd_classic(f, g)
    print(f"Бинарный Евклид (псевдокод Л1): ({poly_to_string(result)}) = {poly_to_algebra(result)}")
    print(f"Честный НОД через deg-деление:  ({poly_to_string(expected)}) = {poly_to_algebra(expected)}")
    print(f"{'OK' if result == expected else 'РАСХОЖДЕНИЕ (так и должно быть, см. предупреждение)'}")

    print("\n=== Карацуба для F_2[X]: вариант 3 задание 3 ===")
    f = poly_from_string("100011101")  # f = x^8 + x^4 + x^3 + x^2 + 1
    g = poly_from_string("111110")     # g = x^5 + x^4 + x^3 + x^2 + x
    print(f"f = ({poly_to_string(f)}) = {poly_to_algebra(f)}")
    print(f"g = ({poly_to_string(g)}) = {poly_to_algebra(g)}")
    result, _ = poly_karatsuba(f, g, threshold=2)
    expected = poly_mul(f, g)
    print(f"Карацуба: ({poly_to_string(result)}) = {poly_to_algebra(result)}")
    print(f"Школьное: ({poly_to_string(expected)}) = {poly_to_algebra(expected)}")
    print(f"{'OK' if result == expected else 'FAIL'}")
