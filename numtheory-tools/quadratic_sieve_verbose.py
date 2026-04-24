from math import gcd, isqrt
from sympy import isprime


def _legendre_symbol_mod_prime(a, p):
    if p == 2:
        return 1
    a %= p
    if a == 0:
        return 0
    return pow(a, (p - 1) // 2, p)


def _trial_factor_over_base(value, factor_base):
    exponents = []
    rest = value
    for prime in factor_base:
        count = 0
        while rest % prime == 0:
            rest //= prime
            count += 1
        exponents.append(count)
    return rest == 1, exponents


def _gf2_left_kernel_basis(parity_rows):
    # parity_rows: m x n, ищем c (длины m), где c^T * A = 0
    if not parity_rows:
        return []

    m = len(parity_rows)
    n = len(parity_rows[0])
    at = [[parity_rows[row][col] for row in range(m)] for col in range(n)]  # n x m

    row = 0
    pivot_cols = []
    for col in range(m):
        pivot = None
        for r in range(row, n):
            if at[r][col]:
                pivot = r
                break
        if pivot is None:
            continue

        at[row], at[pivot] = at[pivot], at[row]

        for r in range(n):
            if r != row and at[r][col]:
                at[r] = [x ^ y for x, y in zip(at[r], at[row])]

        pivot_cols.append(col)
        row += 1
        if row == n:
            break

    pivot_set = set(pivot_cols)
    free_cols = [c for c in range(m) if c not in pivot_set]

    if not free_cols:
        return []

    # Нужно понять, какая строка соответствует какому pivot-столбцу
    pivot_row_of_col = {}
    for r in range(min(n, len(pivot_cols))):
        pivot_col = None
        for c in range(m):
            if at[r][c]:
                pivot_col = c
                break
        if pivot_col is not None:
            pivot_row_of_col[pivot_col] = r

    basis = []
    for free in free_cols:
        vec = [0] * m
        vec[free] = 1
        for pivot_col in pivot_cols:
            pr = pivot_row_of_col.get(pivot_col)
            if pr is None:
                continue
            vec[pivot_col] = at[pr][free]
        basis.append(vec)

    return basis


def quadratic_sieve_verbose(N, B, max_iterations=20000):
    steps = []

    if N < 2:
        steps.append({"iteration": "Ошибка", "description": "N должно быть >= 2"})
        return None, steps

    if B < 2:
        steps.append({"iteration": "Ошибка", "description": "B должно быть >= 2"})
        return None, steps

    factor_base = []
    steps.append(
        {
            "iteration": "FB.0",
            "phase": "init",
            "description": f"Строим фактор-базу: простые p <= {B}, для которых N является квадратичным вычетом mod p.",
        }
    )

    for p in range(2, B + 1):
        if not isprime(p):
            continue
        symbol = _legendre_symbol_mod_prime(N, p)
        if symbol == 1 or (N % p == 0):
            factor_base.append(p)

    if not factor_base:
        steps.append(
            {
                "iteration": "Ошибка",
                "description": "Фактор-база пуста, попробуйте увеличить B.",
            }
        )
        return None, steps

    steps.append(
        {
            "iteration": "FB.final",
            "phase": "factor_base",
            "description": f"Фактор-база готова: FB = {[-1] + factor_base}.",
        }
    )

    sqrt_n = isqrt(N)
    target_relations = len(factor_base) + 3
    steps.append(
        {
            "iteration": "INIT",
            "phase": "init",
            "description": f"sqrt(N) = {sqrt_n}, размер фактор-базы = {len(factor_base)}, целевое число отношений = {target_relations}.",
        }
    )

    relations = []
    parity_rows = []
    relation_values = []
    offset = -1

    for iteration in range(1, max_iterations + 1):
        v = sqrt_n + offset
        q = v * v - N
        abs_q = abs(q)

        smooth, exponents = _trial_factor_over_base(abs_q, factor_base)
        if smooth:
            sign_bit = 1 if q < 0 else 0
            full_vec = [sign_bit] + exponents
            parity_vec = [x % 2 for x in full_vec]

            relations.append(full_vec)
            parity_rows.append(parity_vec)
            relation_values.append(v)

            factor_str = " * ".join(
                f"{p}^{a}" for p, a in zip(factor_base, exponents) if a > 0
            )
            if not factor_str:
                factor_str = "1"

            vector_description = f"v(sign, p1..pk) = {full_vec}"

            steps.append(
                {
                    "iteration": iteration,
                    "phase": "relation",
                    "description": (
                        f"offset={offset}, v={v}, Q=v^2-N={q}. "
                        f"|Q| гладкое: {abs_q} = {factor_str}. "
                        f"Вектор разложения: {vector_description}. "
                        f"Добавили отношение #{len(relations)}."
                    ),
                }
            )

        offset += 1

        if len(relations) < target_relations:
            continue

        matrix_lines = []
        for row_idx, row in enumerate(parity_rows, start=1):
            matrix_lines.append(f"r{row_idx} = {row}")
        matrix_block = "\n".join(matrix_lines) if matrix_lines else "(пусто)"

        steps.append(
            {
                "iteration": f"LIN.{iteration}",
                "phase": "linear",
                "description": (
                    f"Собрано отношений: {len(relations)}. Запускаем линейную алгебру над GF(2), "
                    "ищем зависимости между строками."
                ),
            }
        )
        steps.append(
            {
                "iteration": f"LIN.{iteration}.matrix",
                "phase": "linear",
                "description": (
                    f"Двоичная матрица показателей mod 2 ({len(parity_rows)}x{len(parity_rows[0])}):\n"
                    f"{matrix_block}"
                ),
            }
        )

        dependencies = _gf2_left_kernel_basis(parity_rows)
        if not dependencies:
            steps.append(
                {
                    "iteration": f"LIN.{iteration}.0",
                    "phase": "linear",
                    "description": "Ядро пустое на текущем наборе отношений, продолжаем сбор.",
                }
            )
            continue

        steps.append(
            {
                "iteration": f"LIN.{iteration}.k",
                "phase": "linear",
                "description": f"Найдено зависимостей: {len(dependencies)}.",
            }
        )

        for dep_idx, dep in enumerate(dependencies, start=1):
            selected = [idx for idx, bit in enumerate(dep) if bit == 1]
            if len(selected) < 2:
                continue

            left = 1
            exp_sums = [0] * (len(factor_base) + 1)
            for idx in selected:
                left = (left * relation_values[idx]) % N
                exp_sums = [a + b for a, b in zip(exp_sums, relations[idx])]

            right = 1
            for prime, exp in zip([-1] + factor_base, exp_sums):
                right *= pow(prime, exp // 2)
            right %= N

            g1 = gcd(abs(right - left), N)
            steps.append(
                {
                    "iteration": f"DEP.{iteration}.{dep_idx}",
                    "phase": "dependency",
                    "description": (
                        f"Зависимость использует отношения {selected}. "
                        f"x = prod(v_i) mod N = {left}, y = sqrt(prod(Q_i)) mod N = {right}, "
                        f"gcd(|y-x|, N) = {g1}."
                    ),
                }
            )

            if 1 < g1 < N:
                steps.append(
                    {
                        "iteration": "Результат",
                        "description": f"Найден нетривиальный делитель: N = {g1} * {N // g1}.",
                    }
                )
                return (int(g1), int(N // g1)), steps

    steps.append(
        {
            "iteration": "Ошибка",
            "description": (
                f"Не удалось найти делитель за {max_iterations} итераций. "
                "Попробуйте увеличить B или max_iterations."
            ),
        }
    )
    return None, steps
