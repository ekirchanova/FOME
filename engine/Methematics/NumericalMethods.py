def newton_method(x0, f, f_d, Nmax=30):
    x = x0

    REPEAT_AMOUNT = 1
    while Nmax > REPEAT_AMOUNT:
        x = x - f(x) / f_d(x)
        REPEAT_AMOUNT += 1

    # print("repeats: " + str(REPEAT_AMOUNT))
    return (x, REPEAT_AMOUNT)


def dichotomy_method(x_left, x_right, f, Nmax=30):
    a = x_left
    b = x_right
    c = (a + b) / 2

    f_a = f(a)
    f_b = f(b)
    f_c = f(c)

    REPEAT_AMOUNT = 1
    while Nmax > REPEAT_AMOUNT:
        if f_c * f_b < 0:
            a = c
            f_a = f_c
        if f_c * f_a < 0:
            b = c
            f_b = f_c
        c = (a + b) / 2
        f_c = f(c)
        REPEAT_AMOUNT += 1

    return (c, REPEAT_AMOUNT)