"""
Быстрое возведение в степень по модулю (бинарный метод) по Л1 стр. 32.

Идея: x^e = Π x^{e_i · 2^i}, где e = Σ e_i · 2^i — двоичная запись e.
Сложность по числу умножений: не более L₂(e) = длина e в битах (= ⌊log₂ e⌋ + 1).

Алгоритм:
  Проходим биты e от старшего к младшему. На каждом шаге:
    — возводим текущее значение в квадрат (всегда)
    — если текущий бит = 1, дополнительно умножаем на x
  В конце получаем x^e mod N.
"""


def fast_pow_mod(x, e, N):
    """
    Бинарный метод возведения в степень: x^e mod N.

    Возвращает (результат, шаги).
    """
    steps = []
    x, e, N = int(x), int(e), int(N)

    if N < 1:
        steps.append({'iteration': 'Ошибка', 'description': 'N должно быть ≥ 1'})
        return None, steps
    if e < 0:
        steps.append({'iteration': 'Ошибка',
                      'description': 'Отрицательные степени не поддерживаются '
                                     '(для них нужен обратный по модулю)'})
        return None, steps

    if e == 0:
        steps.append({
            'iteration': 'Результат',
            'description': f"e = 0 ⇒ x^0 = 1 mod {N} = {1 % N if N > 1 else 0}"
        })
        return 1 % N if N > 1 else 0, steps

    # Двоичное представление степени (старший бит слева)
    bits = bin(e)[2:]
    L = len(bits)

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Вычисляем x^e mod N при x = {x}, e = {e}, N = {N}. "
                       f"Двоичная запись e = {bits}₂ (длина L₂(e) = {L})."
    })

    # Проходим биты слева направо
    result = 1
    for i, bit in enumerate(bits):
        # 1) квадрат
        prev = result
        result = (result * result) % N
        phase_desc = f"шаг {i+1}/{L}, бит '{bit}': квадрат: {prev}² mod {N} = {result}"

        # 2) если бит 1 — умножить на x
        if bit == '1':
            prev = result
            result = (result * x) % N
            phase_desc += f";  бит = 1 ⇒ умножаем на x: {prev}·{x} mod {N} = {result}"
        else:
            phase_desc += ";  бит = 0, умножение пропускаем"

        steps.append({
            'iteration': i + 1,
            'description': phase_desc
        })

    steps.append({
        'iteration': 'Результат',
        'description': f"x^e mod N = {x}^{e} mod {N} = {result}.  "
                       f"Потребовалось {L} возведений в квадрат и "
                       f"{bits.count('1')} дополнительных умножений — "
                       f"всего не более 2·L₂(e) = {2*L} операций."
    })

    # Сверка
    expected = pow(x, e, N)
    assert result == expected, f"Несовпадение: {result} vs {expected}"

    return result, steps


if __name__ == '__main__':
    print("=== Быстрое возведение в степень ===")
    for x, e, N in [(3, 13, 31), (2, 100, 181), (5, 17, 23), (7, 0, 13)]:
        r, _ = fast_pow_mod(x, e, N)
        expected = pow(x, e, N)
        print(f"  {x}^{e} mod {N} = {r}  (ожидали {expected})  "
              f"{'OK' if r == expected else 'FAIL'}")
