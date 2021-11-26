import numpy as num
import matplotlib.pyplot as pyplot
import sympy as simp
from sympy import *
import control.matlab as matlab
import colorama as color

def choice():

    need_new_choice = True

    while need_new_choice:

        print(color.Style.RESET_ALL)

        user_input = input('Введите номер команды: \n'
                        '1 - Переходная характеристика;\n' 
                        '2 - Критерий Найквиста;\n' 
                        '3 - Логарифмические характеристики;\n'
                        '4 - Критерий Михайлова;\n' 
                        '5 - Нахождение граничного Кос;\n'
                        '6 - Нахождение полюсов;\n' 
                        '7 - Всё сразу;\n'
                        '8 - Закончить выполнение;\n')
        # Проверяем, число ли было введено

        if not user_input.isdigit() or int(user_input)>8 or int(user_input) == 0:
            print(color.Fore.RED + '\nНедопустимое значение')
        else:
            user_input = int(user_input)
            need_new_choice = False

    return user_input

def step(W):
    x = []
    for i in range(0, 100000):
        x.append(i/1000)
    ampl = 'Амплитуда'
    time = 'Время (с)'
    pyplot.figure()
    [y, x] = matlab.step(W, x)
    pyplot.grid(True)
    pyplot.plot(x, y, 'purple')
    pyplot.axis([-10.0, 100.0, -500.0, 500.0])
    pyplot.ylabel(ampl)
    pyplot.xlabel(time)
    pyplot.title('Переходная характеристика')
    return

def log(W):
    x = []
    for i in range(0, 100000):
        x.append(i/1000)
    pyplot.figure()
    matlab.bode(W)
    pyplot.grid(True)
    pyplot.plot()
    return

def polus(W):
    flag = True
    ex = matlab.pole(W)
    print(ex)
    for i in range(len(ex)):
        if ex[i]>0:
            flag = False

    # построение корней на плоскости Re,Im
    pyplot.figure()
    pyplot.grid(True)
    root = matlab.pzmap(W)
    pyplot.axis([-1, 1, -1, 1])  # задаем границы графика по осям x и y
    pyplot.grid(True)
    if flag:
        print('Система устойчива')
    else:
        print('Система не устойчива')
    return

def gurvits(dzam):
    flag = True
    k=0.5
    while flag:
        k=k-0.000001
        denom1 = dzam[0]
        denom = denom1[0]

        denom[4]= 1 + 46*k

        i1 = [denom[1], denom[3], 0, 0]
        i2 = [denom[0], denom[2], denom[4], 0]
        i3 = [0, denom[1], denom[3], 0]
        i4 = [0, denom[0], denom[2], denom[4]]

        matrix_full = num.array([i1, i2, i3, i4])
        det_full = num.linalg.det(matrix_full)

        i1.pop()
        i2.pop()
        i3.pop()

        matrix3 = num.array([i1, i2, i3])
        det_3 = num.linalg.det(matrix3)

        i1.pop()
        i2.pop()
        matrix2 = num.array([i1, i2])
        det_2 = num.linalg.det(matrix2)


        if det_full>=0:
        # if det_full <= 0:
            flag = False
            print("Матрица Гурвица:\n", matrix_full)
            # print('Определители миноров: ' % det_full, % det_3, % det_2)
            print('Определители миноров: %s' %det_full,det_3, det_2)
            print('Граничное значение: %s' %k)
            return k
    return

def nyq(wraz):
    pyplot.figure()
    pyplot.title('Диаграмма Найквиста ')
    matlab.nyquist(wraz)
    pyplot.grid(True)
    pyplot.plot()
    return

def miha(dzam):
    denom1 = dzam[0]
    denom = denom1[0]
    w = simp.symbols('w', real=True)
    # p = simp.factor(denom[0]*(w*I)**6 + denom[1]*(w*I)**5 + denom[2]*(w*I)**4+denom[3]*(w*I)**3+denom[4]*(w*I)**2 + denom[5]*(w*I)**1  + denom[6]*(w*I) + denom[7])

    p = simp.factor(denom[0] * (w * I) ** 4 + denom[1] * (w * I) ** 3 + denom[2] * (w * I) ** 2 + denom[3] * (w * I) ** 1 + denom[4])
    print("Характеристический многочлен замкнутой системы -\n%s" % p)

    pr = re(p)
    pm = im(p)
    print("Действительная часть Re= %s" % pr)
    print("Мнимая часть Im= %s" % pm)

    pyplot.figure()
    pyplot.title('Годограф Михайлова ')
    x = [pr.subs({w: q}) for q in num.arange(0, 100, 0.1)]
    y = [pm.subs({w: q}) for q in num.arange(0, 100, 0.1)]
    pyplot.axis([-100.0, 100.0, -100.0, 100.0])
    pyplot.plot(x, y)
    pyplot.grid(True)
    return

# выбор звена
need_choice = True
while need_choice:
    # kos = 0.12182300000656068
    kos = 0.5
    Tos = 1
    Tg = 6.4
    kp = 2
    Tp = 4
    ku = 23
    Tu = 5

    unit_1 = matlab.tf([kos], [Tos, 1])
    unit_2 = matlab.tf([1], [Tg, 1])
    unit_3 = matlab.tf([kp], [Tp, 1])
    unit_4 = matlab.tf([ku], [Tu, 1])

    unit_5 = unit_2*unit_3*unit_4

    unit_raz =unit_1*unit_2*unit_3*unit_4

    unit_zam_pr = unit_5/(1+unit_1*unit_5)

    unit_zam = matlab.tf([46, 46], [128, 205.6, 93, 16.4, 24])

    unit_zam = matlab.tf([46, 46], [128, 205.6, 93, 16.4, 6.603858])

    dzam = unit_zam.den

    print (unit_raz)
    print (unit_zam)

    task = choice()

    if task == 1:
        step(unit_zam)
    elif task == 2:
        nyq(unit_raz)
    elif task == 3:
        log(unit_raz)
    elif task == 4:
        miha(dzam)
    elif task == 5:
        kos = gurvits(dzam)
    elif task == 6:
        polus(unit_zam)
    elif task == 7:
        step(unit_zam)
        nyq(unit_raz)
        log(unit_raz)
        miha(dzam)
        kos = gurvits(dzam)
        polus(unit_zam)
    else:
        need_choice = False


    pyplot.show()
