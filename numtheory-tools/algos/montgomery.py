"""
Алгоритм Монтгомери по лекции 1 Полякова М.В., стр. 36–40.

Определение (стр. 36):
  x ⊗ y = x·y·R⁻¹ (mod N),
где R > N, (R, N) = 1, и известно представление R·R' − N·N' = 1, 0 < R' < N.

Функция φ_R : Z_N → Z_N, φ_R : x ↦ x·R⁻¹ (mod N).

Алгоритм вычисления φ_R(x) (Л1 стр. 37):
  1: m = x·N' mod R
  2: t = (x + m·N) / R
  3: Если 0 ≤ t < N: вернуть t
  4: Иначе:           вернуть t − N.
"""

from sympy import gcd as sympy_gcd


def _extended_gcd(a, b):
    """Расширенный Евклид: возвращает (g, x, y): a·x + b·y = g."""
    if b == 0:
        return a, 1, 0
    g, x1, y1 = _extended_gcd(b, a % b)
    return g, y1, x1 - (a // b) * y1


def montgomery_reduce(x, N, R=None):
    """
    Применяет φ_R(x) = x·R⁻¹ mod N по алгоритму со слайда стр. 37.
    Если R не задан — выбирается минимальная степень 2, большая N
    (стандартный выбор, т.к. деление на R = сдвиг).

    Возвращает (результат, шаги, параметры_R_N'_R').
    """
    steps = []
    N = int(N)
    if R is None:
        # Берём R как наименьшую степень 2, большую N — удобно для сдвигов
        R = 1
        while R <= N:
            R <<= 1

    if int(sympy_gcd(R, N)) != 1:
        steps.append({'iteration': 'Ошибка',
                      'description': f"Требуется НОД(R, N) = 1, а получили "
                                     f"НОД({R}, {N}) = {int(sympy_gcd(R, N))}"})
        return None, steps, None

    # Расширенный Евклид: R·R' − N·N' = 1  ⇒  R·R' + N·(−N') = 1
    # То есть если Euclid даёт R·u + N·v = 1, то R' = u mod N, N' = (−v) mod R
    g, u, v = _extended_gcd(R, N)
    # u = R⁻¹ mod N → R' = u mod N
    # v = N⁻¹ mod R с противоположным знаком → N' = (−v) mod R
    R_prime = u % N
    N_prime = (-v) % R

    # Проверка: R·R' − N·N' = 1
    check = R * R_prime - N * N_prime
    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Параметры: N = {N}, R = {R} (R > N, НОД(R, N) = 1). "
                       f"Ищем R' и N' такие, что R·R' − N·N' = 1."
    })
    steps.append({
        'iteration': 1,
        'description': f"Расширенный алгоритм Евклида для (R, N) = ({R}, {N}) даёт: "
                       f"R' = {R_prime}, N' = {N_prime}. "
                       f"Проверка: R·R' − N·N' = {R}·{R_prime} − {N}·{N_prime} = {check} "
                       f"{'✓' if check == 1 else '✗'}"
    })

    # Алгоритм: 1) m = x·N' mod R; 2) t = (x + m·N)/R; 3) возврат t или t−N
    if not (0 <= x < R * N):
        steps.append({
            'iteration': '⚠',
            'description': f"Замечание: алгоритм корректен для 0 ≤ x < R·N = {R*N}. "
                           f"Для x = {x} сначала приведём: x = x mod (R·N) = {x % (R*N)}."
        })
        x = x % (R * N)

    m = (x * N_prime) % R
    steps.append({
        'iteration': 2,
        'description': f"Шаг 1: m = x · N' mod R = {x} · {N_prime} mod {R} = {m}"
    })

    t_num = x + m * N
    t = t_num // R
    steps.append({
        'iteration': 3,
        'description': f"Шаг 2: t = (x + m·N) / R = ({x} + {m}·{N}) / {R} = "
                       f"{t_num} / {R} = {t}"
    })

    if 0 <= t < N:
        steps.append({
            'iteration': 4,
            'description': f"Шаг 3: 0 ≤ t < N ({t} < {N}) ⇒ возвращаем t = {t}"
        })
        result = t
    else:
        result = t - N
        steps.append({
            'iteration': 4,
            'description': f"Шаг 3: t ≥ N ({t} ≥ {N}) ⇒ возвращаем t − N = {t} − {N} = {result}"
        })

    # Проверка: φ_R(x) должно равняться x·R⁻¹ mod N
    expected = (x * R_prime) % N
    steps.append({
        'iteration': 'Результат',
        'description': f"φ_R(x) = {result}.  "
                       f"Проверка: x·R⁻¹ mod N = {x}·{R_prime} mod {N} = {expected} "
                       f"{'✓' if expected == result else '✗'}"
    })
    return result, steps, {'R': R, 'R_prime': R_prime, 'N_prime': N_prime}


def montgomery_multiply(x, y, N, R=None):
    """
    Умножение по Монтгомери: x ⊗ y = x·y·R⁻¹ mod N (определение стр. 36).
    По алгоритму «Умножение по модулю» (стр. 40):
      1: z = φ_R(x·y) = x·y·R⁻¹ mod N
      2: возвращаем x·y mod N = φ_R⁻¹(z) = z·R mod N
    На практике «⊗» — это именно операция z из шага 1; для удобства
    возвращаем и её, и обычное x·y mod N.

    Задача из Варианта 5: N = 299, x = 155, y = 117.
    """
    steps = []
    x, y, N = int(x), int(y), int(N)

    if R is None:
        R = 1
        while R <= N:
            R <<= 1

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Умножение по Монтгомери: x ⊗ y = x·y·R⁻¹ mod N, "
                       f"где x = {x}, y = {y}, N = {N}, R = {R}."
    })

    # Сначала посчитаем x·y (это подставим в φ_R)
    xy = x * y
    steps.append({
        'iteration': 1,
        'description': f"Вычисляем x·y = {x}·{y} = {xy}. Теперь применим φ_R к xy."
    })

    # φ_R(xy) = xy·R⁻¹ mod N
    reduce_result, reduce_steps, params = montgomery_reduce(xy, N, R)
    if reduce_result is None:
        steps.extend(reduce_steps)
        return None, steps

    # Сместим нумерацию шагов для читаемости
    for s in reduce_steps:
        s_copy = dict(s)
        if isinstance(s_copy.get('iteration'), int):
            s_copy['iteration'] = f"φ.{s_copy['iteration']}"
        else:
            s_copy['iteration'] = f"φ.{s_copy['iteration']}"
        steps.append(s_copy)

    # Окончательный ответ: обычное x·y mod N — для справки
    xy_mod_N = xy % N
    steps.append({
        'iteration': 'Итог ⊗',
        'description': f"x ⊗ y = x·y·R⁻¹ mod N = {reduce_result}"
    })
    steps.append({
        'iteration': 'Итог ·',
        'description': f"Для сравнения: обычное x·y mod N = {xy} mod {N} = {xy_mod_N}"
    })
    return reduce_result, steps


# ────────────────────────────────────────────────────────────
# Self-test
# ────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("=== Вариант 5 задание 3: N = 299, элементы 155 и 117 ===\n")
    result, steps = montgomery_multiply(155, 117, 299)
    for s in steps:
        print(f"  [{s['iteration']}] {s['description']}")
    print(f"\n⇒ 155 ⊗ 117 = {result} (в представлении Монтгомери)")
    print(f"  Обычное 155·117 mod 299 = {(155 * 117) % 299}")

    print("\n=== Независимая проверка ===")
    # Возьмём пример попроще: N=7, R=8, x=3, y=5
    print("N=7, x=3, y=5, R=8: ожидаем 3·5·8⁻¹ mod 7")
    # 8 mod 7 = 1, значит 8⁻¹ mod 7 = 1; 3·5·1 mod 7 = 15 mod 7 = 1
    r, _ = montgomery_multiply(3, 5, 7, 8)
    print(f"  результат {r}, ожидали 1: {'OK' if r == 1 else 'FAIL'}")
