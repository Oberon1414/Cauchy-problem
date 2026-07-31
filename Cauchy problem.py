import matplotlib.pyplot as plt
from math import pi, cos, sin

global ksi, A, B
ksi = 1/20
A = 1/30
B = 1/15


def check_value(x0, x, y1_0, y2_0):
    res = []
    q, p = (A * B) ** 0.5, (B / A) ** 0.5
    c = (p * y1_0 * sin(q * x0) + y2_0 * cos(q * x0))/p
    c1 = (y1_0 - c * sin(q * x0)) / cos(q * x0)
    res.append(c * sin(q * x) + c1 * cos(q * x))
    res.append(p * c * cos(q * x) - p * c1 * sin(q * x))
    return res


def constant_step(y1_0, y2_0, xn, h):
    # n = int(xn//h + 2)
    n = 1
    while xn > (n - 1) * h:
        n += 1
    y1 = [y1_0]
    y2 = [y2_0]

    b2 = 1 / (2 * ksi)
    b1 = 1 - b2

    for i in range(1, n):
        if i == n - 1:
            h = xn - (i - 1) * h
        else:
            h = h
        k1_1 = h * A * y2[i - 1]
        k1_2 = -h * B * y1[i - 1]
        k2_1 = h * (A * y2[i - 1] - ksi * h * A * B * y1[i - 1])
        k2_2 = h * (-B * y1[i - 1] - ksi * h * A * B * y2[i - 1])
        y1.append(y1[-1] + b1 * k1_1 + b2 * k2_1)
        y2.append(y2[-1] + b1 * k1_2 + b2 * k2_2)

    res = []
    res.append(y1[n - 1])
    res.append(y2[n - 1])
    return res, y1, y2


def with_Runge_abs(y1_0, y2_0, xn, s=2, eps=10**(-4)):
    h = (eps / ((1 / xn)**(s + 1) + ((B * B * y1_0 * y1_0 + A * A * y2_0 * y2_0)**0.5)**(s + 1)))**(1 / (s + 1))
    flag = True
    res = []
    while flag:
        res_prev, y1, y2 = constant_step(y1_0, y2_0, xn, h)
        res, y1, y2 = constant_step(y1_0, y2_0, xn, h/2)
        err = (((res[0] - res_prev[0]) / (2**s - 1))**2 + ((res[1] - res_prev[1]) / (2**s - 1))**2)**0.5
        if err < eps:
            flag = False
        else:
            h = h / 2
    return res


def y_with_a_step(x, y1, y2, h):
    # y(x+h) = ...
    b2 = 1 / (2 * ksi)
    b1 = 1 - b2

    k1_1 = h * A * y2
    k1_2 = -h * B * y1
    k2_1 = h * (A * y2 - ksi * h * A * B * y1)
    k2_2 = h * (-B * y1 - ksi * h * A * B * y2)

    res = []
    res.append(y1 + b1 * k1_1 + b2 * k2_1)
    res.append(y2 + b1 * k1_2 + b2 * k2_2)
    res.append(x + h) #???
    return res


def auto_step(x0, y1_0, y1_n, y2_0, y2_n, h0, s=2):
    res_first = y_with_a_step(x0, y1_0, y2_0, h0 / 2)
    res_second = y_with_a_step(res_first[2], res_first[0], res_first[1], h0 / 2)

    delta0 = res_second[0] - y1_n
    delta1 = res_second[1] - y2_n

    res = []
    res.append(res_second[0])
    res.append(res_second[1])
    res.append(res_second[2])
    res.append(((delta0 * delta0 + delta1 * delta1)**0.5) / (1 - 2**(-s)))
    return res


def with_Runge_loc(x0, xn, y1_0, y2_0, s=2, eps=10**(-7)):
    steps = []
    y1 = [y1_0]
    y2 = [y2_0]
    x = [x0]
    locerr = []
    count = 0

    h = (eps / ((1 / xn)**(s + 1) + ((B**2 * y1_0**2 + A**2 * y2_0**2)**0.5)**(s + 1))) ** (1 / (s + 1))
    i = 0
    point = 0
    while x[i] + h < xn:
        flag = True
        while flag:
            count += 6
            res_prev = y_with_a_step(point, y1[i], y2[i], h)
            res = auto_step(point, y1[i], res_prev[0], y2[i], res_prev[1], h)
            err = res[-1]
            if err > eps * 2**s:
                h /= 2
            elif (err <= eps * 2**s) and (err > eps):
                steps.append(h)
                x.append(x[i] + h)
                locerr.append(err)
                y1.append(res[0])
                y2.append(res[1])
                h /= 2
                i += 1
                flag = False
            elif (err <= eps) and (err >= eps / (2**(s+1))):
                steps.append(h)
                x.append(x[i] + h)
                locerr.append(err)
                y1.append(res_prev[0])
                y2.append(res_prev[1])
                i += 1
                flag = False
            elif err < eps / (2**(s+1)):
                steps.append(h)
                x.append(x[i] + h)
                locerr.append(err)
                y1.append(res_prev[0])
                y2.append(res_prev[1])
                h *= 2
                i += 1
                flag = False
    h = xn - x[i]
    res_prev = y_with_a_step(point, y1[i], y2[i], h)
    res = auto_step(point, y1[i], res_prev[0], y2[i], res_prev[1], h)
    err = res[-1]

    if (err <= eps * 2**s) and (err > eps):
        steps.append(h)
        x.append(x[i] + h)
        locerr.append(err)
        y1.append(res[0])
        y2.append(res[1])
    elif (err <= eps) and (err >= eps / (2**(s + 1))):
        steps.append(h)
        x.append(x[i] + h)
        locerr.append(err)
        y1.append(res_prev[0])
        y2.append(res_prev[1])
    elif err < eps / (2**(s+1)):
        steps.append(h)
        x.append(h + x[i])
        locerr.append(err)
        y1.append(res_prev[0])
        y2.append(res_prev[1])

    res = []
    res.append(y1[-1])
    res.append(y2[-1])
    return res, x, steps, y1, y2, locerr, count


def opponent_constant_step(y1_0, y2_0, xn, h0):
    # n = xn//h + 2
    n = 1
    while xn > (n - 1) * h0:
        n += 1
    y1 = [y1_0]
    y2 = [y2_0]
    for i in range(1, n):
        if i == n - 1:
            h = xn - (i - 1) * h0
        else:
            h = h0
        k1_1 = h * A * y2[i - 1]
        k1_2 = -h * B * y1[i - 1]
        k2_1 = h * (A * y2[i - 1] - 0.5 * h * A * B * y1[i - 1])
        k2_2 = h * (-B * y1[i - 1] - 0.5 * h * A * B * y2[i - 1])
        k3_1 = h * (A * y2[i - 1] - h * A * B * y1[i - 1] - h * h * A * A * B * y2[i - 1])
        k3_2 = h * (-B * y1[i - 1] - h * A * B * y2[i - 1] + h * h * A * B * B * y1[i - 1])

        y1.append(y1[-1] + (k1_1 + 4 * k2_1 + k3_1) / 6)
        y2.append(y2[-1] + (k1_2 + 4 * k2_2 + k3_2) / 6)

    res = []
    res.append(y1[-1])
    res.append(y2[-1])
    return res, y1, y2


def opponent_with_Runge_abs(y1_0, y2_0, xn, s=3, eps=10**(-4)):
    h = (eps / ((1 / xn)**(s + 1) + ((B * B * y1_0 * y1_0 + A * A * y2_0 * y2_0)**0.5)**(s + 1)))**(1 / (s + 1))
    flag = True
    res = []
    while flag:
        res_prev, y1, y2 = opponent_constant_step(y1_0, y2_0, xn, h)
        res, y1, y2 = opponent_constant_step(y1_0, y2_0, xn, h / 2)
        err = (((res[0] - res_prev[0]) / (2 ** s - 1)) ** 2 + ((res[1] - res_prev[1]) / (2 ** s - 1)) ** 2) ** 0.5
        if err < eps:
            flag = False
        else:
            h = h / 2
    return res


def opponent_y_with_a_step(x, y1, y2, h):
    k1_1 = h * A * y2
    k1_2 = -h * B * y1
    k2_1 = h * (A * y2 - 0.5 * h * A * B * y1)
    k2_2 = h * (-B * y1 - 0.5 * h * A * B * y2)
    k3_1 = h * (A * y2 - h * A * B * y1 - h * h * A * A * B * y2)
    k3_2 = h * (-B * y1 - h * A * B * y2 + h * h * A * B * B * y1)

    res = []
    res.append(y1 + (k1_1 + 4 * k2_1 + k3_1) / 6)
    res.append(y2 + (k1_2 + 4 * k2_2 + k3_2) / 6)
    res.append(x + h)
    return res


def opponent_auto_step(x0, y1_0, y1_n, y2_0, y2_n, h0, s=3):
    res_first = opponent_y_with_a_step(x0, y1_0, y2_0, h0/2)
    res_less_temp = opponent_y_with_a_step(res_first[2], res_first[0], res_first[1], h0/2)

    delta0 = res_less_temp[0] - y1_n
    delta1 = res_less_temp[1] - y2_n

    res = []
    res.append(res_less_temp[0])
    res.append(res_less_temp[1])
    res.append(res_less_temp[2])
    res.append(((delta0 * delta0 + delta1 * delta1)**0.5)/(1 - 2**(-s)))
    return res


def opponent_with_Runge_loc(x0, xn, y1_0, y2_0, s=3, eps=10**(-7)):
    steps = []
    y1 = [y1_0]
    y2 = [y2_0]
    x = [x0]
    locerr = []
    count = 0

    h = (eps / ((1 / xn) ** (s + 1) + ((B ** 2 * y1_0 ** 2 + A ** 2 * y2_0 ** 2) ** 0.5) ** (s + 1))) ** (1 / (s + 1))
    i = 0
    point = 0
    while x[i] + h < xn:
        flag = True
        while flag:
            res_prev = opponent_y_with_a_step(point, y1[i], y2[i], h)
            res = opponent_auto_step(point, y1[i], res_prev[0], y2[i], res_prev[1], h)
            err = res[-1]
            count += 12
            if err > eps * 2**s:
                h /= 2
            elif (err <= eps * 2**s) and (err > eps):
                steps.append(h)
                x.append(h + x[i])
                locerr.append(err)
                h /= 2
                y1.append(res[0])
                y2.append(res[1])
                i += 1
                flag = False
            elif (err <= eps) and (err >= eps / 2**(s + 1)):
                steps.append(h)
                x.append(h + x[i])
                locerr.append(err)
                y1.append(res_prev[0])
                y2.append(res_prev[1])
                i += 1
                flag = False
            elif err < eps / (2**(s + 1)):
                steps.append(h)
                x.append(h + x[i])
                locerr.append(err)
                y1.append(res_prev[0])
                y2.append(res_prev[1])
                i += 1
                flag = False
                h *= 2

    h = xn - x[i]
    res_prev = opponent_y_with_a_step(point, y1[i], y2[i], h)
    res = auto_step(point, y1[i], res_prev[0], y2[i], res_prev[1], h)
    err = res[-1]

    if (err <= eps*(2**s)) and (err > eps):
        steps.append(h)
        x.append(h + x[i])
        locerr.append(err)
        y1.append(res[0])
        y2.append(res[1])
    elif (err <= eps) and (err >= eps / (2**(s + 1))):
        steps.append(h)
        x.append(h + x[i])
        y1.append(res_prev[0])
        y2.append(res_prev[1])
        locerr.append(err)
    elif err < eps / (2**(s + 1)):
        steps.append(h)
        x.append(h + x[i])
        y1.append(res_prev[0])
        y2.append(res_prev[1])
        locerr.append(err)

    res = []
    res.append(y1[-1])
    res.append(y2[-1])
    return res, x, steps, y1, y2, locerr, count


y1_0 = B*pi
y2_0 = A*pi
x0 = 0
xn = pi
temp = 0
check = check_value(x0, xn, y1_0, y2_0)

#первая часть
# res = with_Runge_abs(y1_0, y2_0, xn)
# res, temp, temp, temp, temp, temp, temp = with_Runge_loc(x0, xn, y1_0, y2_0, eps=10**(-5))
# print(res)
# print(check)
# print(abs(res[0] - check[0]), abs(res[1] - check[1]))

#абсолютная ошибка от икса
# h = 0.07290766052127917
# # h = 0.5
# abs1, y1, y2 = constant_step(y1_0, y2_0, xn, h)
# abs2, y3, y4 = opponent_constant_step(y1_0, y2_0, xn, h)
# x = []
# err = []
# xopp = []
# erropp = []
# for i in range(len(y1) - 1):
#     x.append(x0 + i * h)
#     xopp.append(x[-1])
#     sheck = check_value(x0, x[-1], y1_0, y2_0)
#     err.append(((sheck[0] - y1[i])**2 + (sheck[1] - y2[i])**2)**0.5)
#     erropp.append(((sheck[0] - y3[i])**2 + (sheck[1] - y4[i])**2)**0.5)
# x.append(xn)
# xopp.append(xn)
# err.append(((check[0] - abs1[0])**2 + (check[1] - abs1[1])**2)**0.5)
# erropp.append(((check[0] - abs2[0])**2 + (check[1] - abs2[1])**2)**0.5)
#
# plt.plot(x, err, color='green')
# plt.plot(xopp, erropp, color='red')
# plt.show()

#шаг от икса
# loc1, x1, steps1, temp, temp, temp, temp = with_Runge_loc(x0, xn, y1_0, y2_0)
# # print(loc1)
# loc2, x2, steps2, temp, temp, temp, temp = opponent_with_Runge_loc(x0, xn, y1_0, y2_0)
# # print(loc2)
# h1 = []
# h2 = []
# print(steps1)
# for i in range(1, len(x1)):
#     h1.append(x1[i] - x1[i-1])
# for i in range(1, len(x2)):
#     h2.append(x2[i] - x2[i-1])
# plt.plot(x1[1:], h1, color='green')
# plt.plot(x2[1:], h2, color='red')
# plt.show()

#отношение оценки к погрешности
loc1, x1, steps1, y1, y2, locerr1, temp = with_Runge_loc(x0, xn, y1_0, y2_0, eps=10**(-7))
loc2, x2, steps2, y3, y4, locerr2, temp = opponent_with_Runge_loc(x0, xn, y1_0, y2_0, eps=10**(-10))
err1 = []
err2 = []
for i in range(1, len(y1)-1):
    temp_true1 = check_value(x1[i-1], x1[i], y1[i-1], y2[i-1])
    err1.append((((temp_true1[0] - y1[i]) * (temp_true1[0] - y1[i]) + (temp_true1[1] - y2[i]) * (temp_true1[1] - y2[i]))**0.5) / abs(locerr1[i]))
for i in range(1, len(y3)-1):
    temp_true2 = check_value(x2[i-1], x2[i], y3[i-1], y4[i-1])
    err2.append((((temp_true2[0] - y3[i]) * (temp_true2[0] - y3[i]) + (temp_true2[1] - y4[i]) * (temp_true2[1] - y4[i]))**0.5) / abs(locerr2[i]))
plt.plot(x1[2:], err1, color='green')
plt.plot(x2[2:], err2, color='red')
plt.show()

#количество пересчётов от погрешности
# count1 = []
# count2 = []
# q = []
# for i in range(1, 6):
#     temp, temp, temp, temp, temp, temp, c1 = with_Runge_loc(x0, xn, y1_0, y2_0, eps=10**(-i))
#     temp, temp, temp, temp, temp, temp, c2 = opponent_with_Runge_loc(x0, xn, y1_0, y2_0, eps=10**(-i))
#     count1.append(c1)
#     count2.append(c2)
#     q.append(i)
# plt.plot(q, count1, color='green')
# plt.plot(q, count2, color='red')
# plt.show()
