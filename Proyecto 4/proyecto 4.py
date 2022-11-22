import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np


def RK4(F, y0, t0, tN, h=0.001, Grafica='N', Tabla='N'):
    n = int((tN - t0) / h)
    t = np.zeros((n + 1, 1))
    y = np.zeros((n + 1, 1))
    t[0] = t0
    y[0] = y0

    for i in range(0, n):
        t[i + 1] = t[i] + h
        F1 = F(t[i], y[i])
        F2 = F(t[i] + 0.5 * h, y[i] + 0.5 * h * F1)
        F3 = F(t[i] + 0.5 * h, y[i] + 0.5 * h * F2)
        F4 = F(t[i] + h, y[i] + h * F3)
        y[i + 1] = y[i] + (h / 6) * (F1 + 2 * F2 + 2 * F3 + F4)  # Modifique sólo esta línea

    if Tabla == 'S':
        print('Resultados aproximados, Método del Trapecio')
        print('tiempo   valor y')
        print(t[0].item(), y[0].item())
        for i in range(0, n):
            print(t[i + 1].item(), y[i + 1].item())

    if Grafica == 'S':
        plt.plot(t, y)
        plt.title('Solución numérica Método Runge Kutta 4')
        plt.show()

    return t, y


def MetodoRK5(F, w0, t0, tn, Tol, h=0.01):
    w = []
    t = []
    w.append(w0)
    t.append(t0)
    H = []
    H.append(h)
    i = -1
    primero = True
    while t[i + 1] <= tn:
        repetir = False
        i = i + 1
        s1 = F(t[i], w[i])
        s2 = F(t[i] + 0.25 * h, w[i] + 0.25 * h * s1)
        s3 = F(t[i] + 3 * h / 8, w[i] + 3 * h * s1 / 32 + 9 * h * s2 / 32)
        s4 = F(t[i] + 12 * h / 13, w[i] + 1932 * h * s1 / 2197 - 7200 * h * s2 / 2197 + 845 * h * s3 / 2197)
        s5 = F(t[i] + h, w[i] + 439 * h * s1 / 216 - 8 * h * s2 + 3680 * h * s3 / 513 - 845 * h * s4 / 4104)
        s6 = F(t[i] + 0.5 * h,
               w[i] - 8 * h * s1 / 27 + 2 * h * s2 - 3544 * h * s3 / 2565 + 1859 * h * s4 / 4104 - 11 * h * s5 / 40)
        w1 = w[i] + h * (25 * s1 / 216 + 1408 * s3 / 2565 + 2197 * s4 / 4104 - 0.2 * s5)
        z1 = w[i] + h * (16 * s1 / 135 + 6656 * s3 / 12825 + 28561 * s4 / 56430 - 9 * s5 / 50 + 2 * s6 / 55)
        e1 = abs(w1 - z1)
        h1 = 0.8 * (Tol * abs(w[0]) / e1) ** (0.2) * h
        if (i == 0):
            if (e1 / abs(w1) < Tol):
                w1 = z1
            else:
                i = -1
                repetir = True
                if (primero == True):
                    primero = False
                    h = h1
                    H[0] = h1
                else:
                    h = h / 2
        if repetir == False:
            w.append(w1)
            h = h1
            H.append(h)
            t.append(t[i] + h1)
    return t, w, H



def f1(t, y):
    fun = (y ** 2 - 1) * (y ** 2 - 9)
    return fun


y0 = 1.5
t0 = 0
tn = 6
Tol = 0.000001
h = 0.1

[t2, y2, H] = MetodoRK5(f1, y0, t0, tn, Tol, h)
[t1, y1] = RK4(f1, y0, t0, tn, h)
plt.plot(t1, y1)
plt.plot(t2, y2)
plt.legend(["RK4_clase", "RK4_proyecto"])
plt.show()
print(H)


def f2(t, y):
    fun = -5 * t ** 4 * y
    return fun


[t2, y2, H] = MetodoRK5(f2, y0, t0, tn, Tol, h)
[t1, y1] = RK4(f1, y0, t0, tn, h)
plt.plot(t1, y1)
plt.plot(t2, y2)
plt.legend(["RK4_clase", "RK4_proyecto"])
plt.show()
print(H)


# def f3(t, y):
#     fun = 3 *y * t ** 2 + 6 * t ** 2
#     return fun
#
#
# [t2, y2, H] = MetodoRK5(f3, y0, t0, tn, Tol, h)
# [t1, y1] = RK4(f1, y0, t0, tn, h)
# plt.plot(t1, y1)
# plt.plot(t2, y2)
# plt.legend(["RK4_clase", "RK4_proyecto"])
# plt.show()
# print(H)


def f4(t, y):
    fun = (np.power(y, 3) - 2) * (np.power(y, 2) - 7)
    return fun


[t2, y2, H] = MetodoRK5(f4, y0, t0, tn, Tol, h)
[t1, y1] = RK4(f1, y0, t0, tn, h)
plt.plot(t1, y1)
plt.plot(t2, y2)
plt.legend(["RK4_clase", "RK4_proyecto"])
plt.show()
print(H)
