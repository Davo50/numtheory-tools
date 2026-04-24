from flask import Flask, render_template, request, jsonify
from sympy import sqrt, lcm, gcd, isprime, sympify, Symbol
from sympy.ntheory.primetest import is_square
from quadratic_sieve_verbose import quadratic_sieve_verbose
from algos.euclid import (euclid_classic, euclid_binary, euclid_multi,
                          euclid_extended, continued_fraction_rational)
from algos.dlog import bsgs, rho_pollard_dlog, pohlig_hellman
from algos.poly_f2 import (poly_from_string, poly_to_string, poly_to_algebra,
                           poly_euclid_binary, poly_karatsuba)
from algos.montgomery import montgomery_multiply, montgomery_reduce
from algos.fastpow import fast_pow_mod
from algos.index_calculus import index_calculus
from algos.elliptic import (ec_add, ec_double, ec_scalar_mult, ec_order,
                            ec_bsgs, _is_on_curve)

app = Flask(__name__)


def fermat_factorization(N):
    steps = []
    x = int(N ** 0.5) + 1
    iteration = 0
    while True:
        iteration += 1
        x_sq = x ** 2
        y_sq = x_sq - N
        is_sq = is_square(y_sq)
        step = {
            'iteration': iteration,
            'description': f"x = {x}, x² = {x_sq}, y² = x² − N = {y_sq}, "
                           f"√y² ∈ ℤ? {'Да ✓' if is_sq else 'Нет ✗'}"
        }
        steps.append(step)
        if is_sq:
            y = int(sqrt(y_sq))
            factor1 = x - y
            factor2 = x + y
            steps.append({
                'iteration': 'Результат',
                'description': f"y = {y}, N = (x − y)(x + y) = ({x} − {y})({x} + {y}) = {factor1} × {factor2}"
            })
            return (int(factor1), int(factor2)), steps
        x += 1
        if iteration > 10000:
            steps.append({'iteration': 'Ошибка', 'description': 'Превышено максимальное число итераций'})
            return None, steps


def p1_pollard(N, B):
    steps = []
    k = 1
    for i in range(2, B + 1):
        k = int(lcm(k, i))

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Вычисляем k = LCM(1, 2, ..., {B}) = {k}"
    })

    iteration = 0
    for a in range(2, N):
        iteration += 1
        a_k = pow(a, k, N)
        d = int(gcd(a_k - 1, N))

        step = {
            'iteration': iteration,
            'phase': 'search',
            'description': f"a = {a}: a^k mod N = {a}^{k} mod {N} = {a_k}, "
                           f"GCD({a_k} − 1, {N}) = GCD({a_k - 1}, {N}) = {d}"
        }

        if 1 < d < N:
            step['description'] += f" → Найден делитель! ✓"
            steps.append(step)
            steps.append({
                'iteration': 'Результат',
                'description': f"N = {d} × {N // d}"
            })
            return (d, N // d), steps
        elif d == 1:
            step['description'] += " → d = 1, пробуем следующее a"
        else:
            step['description'] += f" → d = N, пробуем следующее a"

        steps.append(step)

        if iteration > 1000:
            steps.append({'iteration': 'Ошибка', 'description': 'Превышено максимальное число итераций'})
            return None, steps

    return None, steps


def rho_pollard(N, f_expr_str, x0):
    """
    f_expr_str — строка вида "x**2 + 5", парсится через sympy.
    Итеративная функция: f(v) = f_expr(v) mod N
    """
    steps = []

    x_sym = Symbol('x')
    try:
        f_expr = sympify(f_expr_str)
    except Exception as e:
        steps.append({'iteration': 'Ошибка', 'description': f'Не удалось распознать многочлен: {e}'})
        return None, steps

    def f(v):
        return int(f_expr.subs(x_sym, v)) % N

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"f(x) = ({f_expr_str}) mod {N}, x₀ = {x0}"
    })

    x1 = f(x0)
    x2 = f(f(x0))
    d = int(gcd(abs(x1 - x2), N))
    iteration = 1

    steps.append({
        'iteration': iteration,
        'description': f"x₁ = f({x0}) = {x1}, "
                       f"x₂ = f(f({x0})) = {x2}, "
                       f"GCD(|{x1} − {x2}|, {N}) = GCD({abs(x1 - x2)}, {N}) = {d}"
    })

    while d == 1:
        iteration += 1
        prev_x1 = x1
        prev_x2 = x2
        x1 = f(x1)
        x2 = f(f(x2))
        d = int(gcd(abs(x1 - x2), N))

        steps.append({
            'iteration': iteration,
            'description': f"x_{iteration} = f({prev_x1}) = {x1}, "
                           f"x_{2 * iteration} = f(f({prev_x2})) = {x2}, "
                           f"GCD(|{x1} − {x2}|, {N}) = GCD({abs(x1 - x2)}, {N}) = {d}"
        })

        if iteration > 10000:
            steps.append({'iteration': 'Ошибка', 'description': 'Превышено максимальное число итераций'})
            return None, steps

    if d == N:
        steps.append({
            'iteration': 'Ошибка',
            'description': f'GCD = N = {N}, метод не сработал с данными параметрами. Попробуйте другой многочлен или x₀.'
        })
        return None, steps

    steps.append({
        'iteration': 'Результат',
        'description': f"GCD = {d} ≠ 1 → Найден делитель! N = {d} × {N // d}"
    })
    return (int(d), int(N // d)), steps


def try_double(P, a, N, steps, step_num):
    d = int(gcd(P[1], N))
    if d != 1:
        steps.append({
            'iteration': step_num,
            'operation': 'double',
            'description': f"Удвоение P = {P}: GCD({P[1]}, {N}) = {d} ≠ 1 → Найден делитель!"
        })
        return (), d
    inv_2y = pow(2 * P[1], -1, N)
    lam = ((3 * P[0] * P[0] + a) * inv_2y) % N
    x3 = (lam * lam - P[0] - P[0]) % N
    y3 = (-(lam * (x3 - P[0]) + P[1])) % N
    steps.append({
        'iteration': step_num,
        'operation': 'double',
        'description': f"Удвоение: 2·{P} → λ = (3·{P[0]}² + {a})·(2·{P[1]})⁻¹ mod {N} = {lam}, "
                       f"Результат: ({x3}, {y3})"
    })
    return (x3, y3), 1


def try_add(P, Q, N, steps, step_num):
    d = int(gcd(abs(P[0] - Q[0]), N))
    if d != 1:
        steps.append({
            'iteration': step_num,
            'operation': 'add',
            'description': f"Сложение P={P}, Q={Q}: GCD(|{P[0]}−{Q[0]}|, {N}) = {d} ≠ 1 → Найден делитель!"
        })
        return (), d
    inv = pow(P[0] - Q[0], -1, N)
    lam = ((P[1] - Q[1]) * inv) % N
    x3 = (lam * lam - P[0] - Q[0]) % N
    y3 = (-(lam * (x3 - P[0]) + P[1])) % N
    steps.append({
        'iteration': step_num,
        'operation': 'add',
        'description': f"Сложение: {P} + {Q} → λ = ({P[1]}−{Q[1]})·({P[0]}−{Q[0]})⁻¹ mod {N} = {lam}, "
                       f"Результат: ({x3}, {y3})"
    })
    return (x3, y3), 1


def try_scalar_mult(k, P, a, N, steps, base_step):
    l = k.bit_length()
    Q = P
    sub_step = 0
    steps.append({
        'iteration': f"{base_step}.init",
        'operation': 'scalar_init',
        'description': f"Скалярное умножение: {k}·P, bin(k) = {bin(k)}, длина = {l} бит"
    })
    for i in range(l - 2, -1, -1):
        sub_step += 1
        bit = (k >> i) & 1
        steps.append({
            'iteration': f"{base_step}.{sub_step}",
            'operation': 'bit',
            'description': f"Бит {l - 1 - i} (позиция {i}): {bit}"
        })
        Q, v = try_double(Q, a, N, steps, f"{base_step}.{sub_step}d")
        if v != 1:
            return v
        if bit:
            Q, v = try_add(P, Q, N, steps, f"{base_step}.{sub_step}a")
            if v != 1:
                return v
    return 1


def lenstra(N, B, x, y, a):
    steps = []

    b = (y ** 2 - x ** 3 - a * x) % N
    disc = -16 * (4 * a ** 3 + 27 * b ** 2)
    d = int(gcd(disc, N))

    steps.append({
        'iteration': 0,
        'phase': 'init',
        'description': f"Кривая: a = {a}, b = y² − x³ − ax = {y}² − {x}³ − {a}·{x} = {b}, "
                       f"Точка P = ({x}, {y})"
    })
    steps.append({
        'iteration': '0.5',
        'phase': 'discriminant',
        'description': f"Δ = −16(4a³ + 27b²) = {disc}, GCD(Δ, N) = {d}"
    })

    if 1 < d < N:
        steps.append({
            'iteration': 'Результат',
            'description': f"GCD дискриминанта и N дал делитель! N = {d} × {N // d}"
        })
        return (d, N // d), steps

    k = 1
    P = (x % N, y % N)
    prime_step = 0
    for p in range(2, B + 1):
        if isprime(p):
            prime_step += 1
            k *= p
            steps.append({
                'iteration': prime_step,
                'phase': 'prime',
                'description': f"Простое p = {p}, k = {k}, вычисляем {k}·P на E mod {N}"
            })
            d = try_scalar_mult(k, P, a, N, steps, prime_step)
            if d != 1:
                d = int(d)
                if 1 < d < N:
                    steps.append({
                        'iteration': 'Результат',
                        'description': f"Найден нетривиальный делитель! N = {d} × {N // d}"
                    })
                    return (d, N // d), steps
                else:
                    steps.append({
                        'iteration': 'Ошибка',
                        'description': f"GCD = {d} (тривиальный делитель), попробуйте другие параметры"
                    })
                    return None, steps

    steps.append({
        'iteration': 'Ошибка',
        'description': f'Делитель не найден при B = {B}. Увеличьте B или смените параметры кривой.'
    })
    return None, steps


@app.route('/')
def index():
    return render_template('index.html')

@app.route('/notes')
def notes():
    return render_template('notes.html')

@app.route('/gcd')
def gcd_page():
    return render_template('gcd.html')

@app.route('/dlog')
def dlog_page():
    return render_template('dlog.html')

@app.route('/extra')
def extra_page():
    return render_template('extra.html')

@app.route('/tools')
def tools_page():
    return render_template('tools.html')

@app.route('/advanced')
def advanced_page():
    return render_template('advanced.html')


@app.route('/api/poly', methods=['POST'])
def api_poly():
    """API для многочленов в F_2[X]. method ∈ {'euclid_binary', 'karatsuba'}"""
    try:
        data = request.get_json()
        method = data.get('method', 'euclid_binary')
        a_str = data.get('a', '')
        b_str = data.get('b', '')
        a = poly_from_string(a_str)
        b = poly_from_string(b_str)

        if method == 'euclid_binary':
            result, steps = poly_euclid_binary(a, b)
        elif method == 'karatsuba':
            threshold = int(data.get('threshold', 2))
            result, steps = poly_karatsuba(a, b, threshold=threshold)
        else:
            return jsonify({'error': 'Неизвестный метод'}), 400

        out = {
            'success': result is not None,
            'result': result,
            'result_bin': poly_to_string(result) if result is not None else None,
            'result_algebra': poly_to_algebra(result) if result is not None else None,
            'steps': steps
        }
        return jsonify(out)
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/montgomery', methods=['POST'])
def api_montgomery():
    """API для умножения по Монтгомери."""
    try:
        data = request.get_json()
        x = int(data.get('x', 0))
        y = int(data.get('y', 0))
        N = int(data.get('N', 0))
        R_raw = data.get('R')
        R = int(R_raw) if R_raw not in (None, '', 'null') else None

        if N < 2:
            return jsonify({'error': 'N должно быть ≥ 2'}), 400

        result, steps = montgomery_multiply(x, y, N, R=R)
        return jsonify({
            'success': result is not None,
            'result': result,
            'steps': steps
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


# ─────────── API: расширенный Евклид ───────────
@app.route('/api/euclid_ext', methods=['POST'])
def api_euclid_ext():
    try:
        data = request.get_json()
        a = int(data.get('a', 0))
        b = int(data.get('b', 0))
        result, steps = euclid_extended(a, b)
        if result is None:
            return jsonify({'success': False, 'steps': steps})
        d, x, y = result
        return jsonify({
            'success': True,
            'd': d, 'x': x, 'y': y,
            'steps': steps
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


# ─────────── API: цепные дроби ───────────
@app.route('/api/continued_fraction', methods=['POST'])
def api_cf():
    try:
        data = request.get_json()
        num = int(data.get('num', 0))
        den = int(data.get('den', 1))
        if den == 0:
            return jsonify({'error': 'Знаменатель не может быть 0'}), 400
        result, steps = continued_fraction_rational(num, den)
        if result is None:
            return jsonify({'success': False, 'steps': steps})
        return jsonify({
            'success': True,
            'quotients': result['quotients'],
            'convergents': result['convergents'],
            'steps': steps
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


# ─────────── API: быстрое возведение в степень ───────────
@app.route('/api/fastpow', methods=['POST'])
def api_fastpow():
    try:
        data = request.get_json()
        x = int(data.get('x', 0))
        e = int(data.get('e', 0))
        N = int(data.get('N', 1))
        result, steps = fast_pow_mod(x, e, N)
        return jsonify({
            'success': result is not None,
            'result': result,
            'steps': steps
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


# ─────────── API: index calculus ───────────
@app.route('/api/index_calculus', methods=['POST'])
def api_index_calculus():
    try:
        data = request.get_json()
        g = int(data.get('g', 2))
        p = int(data.get('p', 0))
        if not isprime(p):
            return jsonify({'error': f'p = {p} не простое'}), 400

        h = data.get('h')
        x_known = data.get('x_known')
        base_bound = data.get('base_bound')
        seed = data.get('seed')

        kwargs = {'g': g, 'p': p}
        if base_bound not in (None, '', 'null'):
            kwargs['base_bound'] = int(base_bound)
        if seed not in (None, '', 'null'):
            kwargs['seed'] = int(seed)
        if h not in (None, '', 'null'):
            kwargs['h'] = int(h)
        elif x_known not in (None, '', 'null'):
            kwargs['x'] = int(x_known)
        else:
            return jsonify({'error': 'Задайте h или x_known'}), 400

        result, steps = index_calculus(**kwargs)
        return jsonify({
            'success': result is not None,
            'result': result,
            'steps': steps
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


# ─────────── API: эллиптические кривые ───────────
@app.route('/api/elliptic', methods=['POST'])
def api_elliptic():
    """
    method ∈ {'add', 'double', 'mult', 'order', 'ecdlp'}
    Параметры кривой: a, b, p. Точки: P = (Px, Py), Q = (Qx, Qy).
    Для mult: k.
    Для ecdlp: k_known (тогда Q вычислим сами) или Q.
    """
    try:
        data = request.get_json()
        method = data.get('method', 'add')
        a = int(data.get('a', 0))
        b = int(data.get('b', 0))
        p = int(data.get('p', 0))
        if p < 2 or not isprime(p):
            return jsonify({'error': f'p = {p} должно быть простым ≥ 2'}), 400

        def parse_point(prefix):
            px = data.get(prefix + 'x')
            py = data.get(prefix + 'y')
            if px in (None, '', 'null') or py in (None, '', 'null'):
                return None
            return (int(px) % p, int(py) % p)

        P = parse_point('P')
        Q_input = parse_point('Q')

        steps = []
        if P is not None and not _is_on_curve(P, a, b, p):
            return jsonify({'error': f'P = {P} не лежит на кривой y² = x³ + {a}x + {b} (mod {p})'}), 400

        if method == 'add':
            if P is None or Q_input is None:
                return jsonify({'error': 'Для сложения нужны обе точки P и Q'}), 400
            if not _is_on_curve(Q_input, a, b, p):
                return jsonify({'error': f'Q = {Q_input} не лежит на кривой'}), 400
            R = ec_add(P, Q_input, a, b, p)
            steps.append({
                'iteration': 'init',
                'phase': 'init',
                'description': f"Кривая E: y² = x³ + {a}x + {b} (mod {p}).  "
                               f"P = {P}, Q = {Q_input}."
            })
            if P == Q_input:
                steps.append({'iteration': 'case',
                              'description': "P = Q ⇒ применяем удвоение (см. соотв. алгоритм)."})
            elif P and Q_input and P[0] == Q_input[0] and (P[1] + Q_input[1]) % p == 0:
                steps.append({'iteration': 'case',
                              'description': "P = −Q ⇒ P + Q = O (бесконечность)."})
            else:
                if P is not None and Q_input is not None:
                    lam_num = (Q_input[1] - P[1]) % p
                    lam_den = (Q_input[0] - P[0]) % p
                    lam = (lam_num * pow(lam_den, -1, p)) % p
                    x3 = (lam * lam - P[0] - Q_input[0]) % p
                    y3 = (lam * (P[0] - x3) - P[1]) % p
                    steps.append({
                        'iteration': 1,
                        'description': f"λ = (y₂ − y₁)/(x₂ − x₁) = "
                                       f"({Q_input[1]} − {P[1]}) · ({Q_input[0]} − {P[0]})⁻¹ mod {p} = {lam}"
                    })
                    steps.append({
                        'iteration': 2,
                        'description': f"x₃ = λ² − x₁ − x₂ = {lam}² − {P[0]} − {Q_input[0]} mod {p} = {x3}"
                    })
                    steps.append({
                        'iteration': 3,
                        'description': f"y₃ = λ·(x₁ − x₃) − y₁ = {lam}·({P[0]} − {x3}) − {P[1]} mod {p} = {y3}"
                    })
            steps.append({'iteration': 'Результат',
                          'description': f"P + Q = {R if R else 'O (бесконечность)'}"})
            return jsonify({'success': True, 'result': list(R) if R else None, 'steps': steps})

        if method == 'double':
            if P is None:
                return jsonify({'error': 'Нужна точка P'}), 400
            R = ec_double(P, a, b, p)
            steps.append({'iteration': 'init', 'phase': 'init',
                          'description': f"Кривая y² = x³ + {a}x + {b} (mod {p}). P = {P}."})
            if P[1] % p == 0:
                steps.append({'iteration': 'case',
                              'description': "y₁ = 0 ⇒ касательная вертикальна ⇒ 2P = O."})
            else:
                lam_num = (3 * P[0] * P[0] + a) % p
                lam_den = (2 * P[1]) % p
                lam = (lam_num * pow(lam_den, -1, p)) % p
                x3 = (lam * lam - 2 * P[0]) % p
                y3 = (lam * (P[0] - x3) - P[1]) % p
                steps.append({'iteration': 1,
                              'description': f"λ = (3x₁² + a)/(2y₁) = "
                                             f"(3·{P[0]}² + {a}) · (2·{P[1]})⁻¹ mod {p} = {lam}"})
                steps.append({'iteration': 2,
                              'description': f"x₃ = λ² − 2x₁ = {lam}² − 2·{P[0]} mod {p} = {x3}"})
                steps.append({'iteration': 3,
                              'description': f"y₃ = λ·(x₁ − x₃) − y₁ = "
                                             f"{lam}·({P[0]} − {x3}) − {P[1]} mod {p} = {y3}"})
            steps.append({'iteration': 'Результат',
                          'description': f"2P = {R if R else 'O'}"})
            return jsonify({'success': True, 'result': list(R) if R else None, 'steps': steps})

        if method == 'mult':
            if P is None:
                return jsonify({'error': 'Нужна точка P'}), 400
            k = int(data.get('k', 0))
            R, m_steps = ec_scalar_mult(k, P, a, b, p, verbose=True)
            steps = m_steps
            return jsonify({'success': True, 'result': list(R) if R else None, 'steps': steps})

        if method == 'order':
            if P is None:
                return jsonify({'error': 'Нужна точка P'}), 400
            order = ec_order(P, a, b, p)
            steps.append({'iteration': 'init', 'phase': 'init',
                          'description': f"Ищем порядок точки P = {P} на кривой y² = x³ + {a}x + {b} (mod {p}) "
                                         f"перебором: вычисляем k·P при k = 1, 2, … до первого k·P = O."})
            steps.append({'iteration': 'Результат',
                          'description': f"ord(P) = {order}"})
            return jsonify({'success': True, 'result': order, 'steps': steps})

        if method == 'ecdlp':
            if P is None:
                return jsonify({'error': 'Нужна точка P'}), 400
            k_known = data.get('k_known')
            if k_known not in (None, '', 'null'):
                k_known = int(k_known)
                Q = ec_scalar_mult(k_known, P, a, b, p)
                steps.append({'iteration': 0, 'phase': 'init',
                              'description': f"Вычислили Q = {k_known}·P = {Q}"})
            else:
                Q = Q_input
                if Q is None:
                    return jsonify({'error': 'Задайте Q или k_known'}), 400

            result, ec_steps = ec_bsgs(P, Q, a, b, p)
            steps.extend(ec_steps)
            return jsonify({'success': result is not None, 'result': result, 'steps': steps})

        return jsonify({'error': 'Неизвестный метод'}), 400

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/gcd', methods=['POST'])
def api_gcd():
    """API для НОД. method ∈ {'classic', 'binary', 'multi'}"""
    try:
        data = request.get_json()
        method = data.get('method', 'classic')

        if method == 'multi':
            raw = data.get('numbers', '')
            nums = [int(x) for x in str(raw).replace(',', ' ').split() if x.strip()]
            if len(nums) < 2:
                return jsonify({'error': 'Нужно минимум 2 числа'}), 400
            result, steps = euclid_multi(nums)
        else:
            a = int(data.get('a', 0))
            b = int(data.get('b', 0))
            if a == 0 and b == 0:
                return jsonify({'error': 'a и b не могут быть оба равны 0'}), 400
            if method == 'classic':
                result, steps = euclid_classic(a, b)
            elif method == 'binary':
                result, steps = euclid_binary(a, b)
            else:
                return jsonify({'error': 'Неизвестный метод'}), 400

        return jsonify({
            'success': result is not None,
            'result': result,
            'steps': steps
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/dlog', methods=['POST'])
def api_dlog():
    """API для дискретного логарифма. method ∈ {'bsgs', 'rho', 'pohlig_hellman'}"""
    try:
        data = request.get_json()
        method = data.get('method', 'bsgs')
        g = int(data.get('g', 0))
        p = int(data.get('p', 0))

        if g < 2 or p < 3:
            return jsonify({'error': 'Нужно: g ≥ 2, p ≥ 3'}), 400
        if not isprime(p):
            return jsonify({'error': f'p = {p} не простое'}), 400

        # h можно задать двумя способами:
        # 1) напрямую ('h')
        # 2) через степень ('x_known') — тогда h = g^x mod p (формат препода)
        h = data.get('h')
        x_known = data.get('x_known')

        kwargs = {'g': g, 'p': p}
        if h not in (None, '', 'null'):
            kwargs['h'] = int(h)
        elif x_known not in (None, '', 'null'):
            kwargs['x'] = int(x_known)
        else:
            return jsonify({'error': 'Задайте либо h, либо степень x'}), 400

        if method == 'bsgs':
            result, steps = bsgs(**kwargs)
        elif method == 'rho':
            result, steps = rho_pollard_dlog(**kwargs)
        elif method == 'pohlig_hellman':
            result, steps = pohlig_hellman(**kwargs)
        else:
            return jsonify({'error': 'Неизвестный метод'}), 400

        return jsonify({
            'success': result is not None,
            'result': result,
            'steps': steps
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/factorize', methods=['POST'])
def factorize():
    try:
        data = request.get_json()
        method = data.get('method')
        N = int(data.get('N', 0))

        if N < 2:
            return jsonify({'error': 'N должно быть ≥ 2'}), 400

        if isprime(N):
            return jsonify({'error': f'{N} — простое число, факторизация невозможна'}), 400

        if method == 'fermat':
            result, steps = fermat_factorization(N)

        elif method == 'p1_pollard':
            B = int(data.get('B', 10))
            if B < 2:
                return jsonify({'error': 'B должно быть ≥ 2'}), 400
            result, steps = p1_pollard(N, B)

        elif method == 'rho_pollard':
            poly_str = data.get('polynomial', 'x**2 + 1')
            x0 = int(data.get('x0', 2))
            result, steps = rho_pollard(N, poly_str, x0)

        elif method == 'lenstra':
            B = int(data.get('B', 10))
            x_coord = int(data.get('x_coord', 0))
            y_coord = int(data.get('y_coord', 1))
            a_param = int(data.get('a_param', 1))
            result, steps = lenstra(N, B, x_coord, y_coord, a_param)

        elif method == 'quadratic_sieve':
            B = int(data.get('B', 20))
            max_iterations = int(data.get('max_iterations', 20000))
            if B < 2:
                return jsonify({'error': 'B должно быть ≥ 2'}), 400
            if max_iterations < 1:
                return jsonify({'error': 'max_iterations должно быть ≥ 1'}), 400
            result, steps = quadratic_sieve_verbose(N, B, max_iterations)

        else:
            return jsonify({'error': 'Неизвестный метод'}), 400

        if result:
            return jsonify({
                'success': True,
                'factors': list(result),
                'steps': steps
            })
        else:
            return jsonify({
                'success': False,
                'steps': steps,
                'error': 'Факторизация не удалась с данными параметрами'
            })

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


if __name__ == '__main__':
    import socket
    import ipaddress

    def _is_real_lan(ip_str):
        """True, если IP принадлежит домашним/офисным подсетям:
        10.0.0.0/8, 172.16.0.0/12, 192.168.0.0/16.
        Отбрасываем: 127.*, 169.254.* (APIPA), CGNAT 100.64.0.0/10,
        зарезервированные (включая 198.18.0.0/15 — его любят раздавать
        VPN-клиенты типа Cloudflare Warp, AnyConnect, ClearVPN и т.п.)."""
        try:
            ip = ipaddress.ip_address(ip_str)
        except ValueError:
            return False
        if ip.is_loopback or ip.is_link_local or ip.is_multicast:
            return False
        # Явно типично-домашние диапазоны
        home_nets = [
            ipaddress.ip_network('10.0.0.0/8'),
            ipaddress.ip_network('172.16.0.0/12'),
            ipaddress.ip_network('192.168.0.0/16'),
        ]
        return any(ip in net for net in home_nets)

    def _get_all_ips():
        """Собирает все IPv4-адреса хоста. Сначала пробуем netifaces
        (если установлен), потом — getaddrinfo по hostname, потом —
        connect-трюк как запасной вариант."""
        ips = set()

        # Способ 1: netifaces — лучший вариант, но это опц. зависимость
        try:
            import netifaces
            for iface in netifaces.interfaces():
                addrs = netifaces.ifaddresses(iface).get(netifaces.AF_INET, [])
                for a in addrs:
                    if 'addr' in a:
                        ips.add(a['addr'])
        except ImportError:
            pass

        # Способ 2: getaddrinfo — работает всегда, но не на всех системах
        # видит все интерфейсы
        try:
            hostname = socket.gethostname()
            for info in socket.getaddrinfo(hostname, None, socket.AF_INET):
                ips.add(info[4][0])
        except Exception:
            pass

        # Способ 3: connect-трюк к Google DNS (как было), как последний шанс
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            s.settimeout(0.2)
            s.connect(('8.8.8.8', 80))
            ips.add(s.getsockname()[0])
            s.close()
        except Exception:
            pass

        return ips

    def pick_lan_ip():
        """Выбираем наиболее подходящий LAN IP: сначала из домашних подсетей,
        потом — любой не-loopback, потом — None."""
        ips = _get_all_ips()
        # 1) 192.168.*, 10.*, 172.16-31.*
        lan = [ip for ip in ips if _is_real_lan(ip)]
        if lan:
            # Предпочтём 192.168.* как типичный Wi-Fi
            lan.sort(key=lambda ip: (not ip.startswith('192.168.'), ip))
            return lan[0], ips
        # 2) хоть что-нибудь не-loopback
        fallback = [ip for ip in ips
                    if not ip.startswith('127.') and not ip.startswith('169.254.')]
        if fallback:
            return fallback[0], ips
        return None, ips

    port = 8080
    lan_ip, all_ips = pick_lan_ip()

    print()
    print("  ┌─ numtheory-tools ─────────────────────────────────────────────")
    print(f"  │  Локально:  http://127.0.0.1:{port}   http://localhost:{port}")
    if lan_ip:
        print(f"  │  По сети:   http://{lan_ip}:{port}")
        print(f"  │             (с телефона/другого ПК в той же Wi-Fi)")
    else:
        print(f"  │  LAN IP не найден — проверь подключение к Wi-Fi")
    print("  └──────────────────────────────────────────────────────────")

    # Покажем все найденные IP — если VPN раздаёт странный адрес,
    # пользователь сразу увидит истинный
    other = sorted(ip for ip in all_ips
                   if ip != lan_ip and not ip.startswith('127.'))
    if other:
        print(f"     Другие интерфейсы (на всякий случай): {', '.join(other)}")
        print(f"     Если выше показан VPN-адрес, отключи VPN или выбери IP из этого списка.")
        print()
    else:
        print()

    app.run(host='0.0.0.0', port=port, debug=False)