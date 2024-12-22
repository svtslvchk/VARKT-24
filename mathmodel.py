import numpy as np
import matplotlib.pyplot as plt
from math import  sqrt
from math import radians


stages = [

    {   #первая ступень
        'thrust': 2816500, #сила тяги (Н)
        'start_mass': 149186, #масса с топливом (кг)
        'end_mass': 48129, #сухая масса (кг)
        'time': 90, #время работы ступени (с)
    },

    {   #вторая ступень
        'thrust': 875820, #сила тяги (Н)
        'start_mass': 31656, #масса с топливом (кг)
        'end_mass': 16243, #сухая масса (кг)
        'time': 47, #время работы ступени (с)
    },

    {   #третья ступень
        'thrust': 302950, #сила тяги (Н)
        'start_mass': 23081, #масса с топливом (кг)
        'end_mass': 16000, #сухая масса (кг)
        'time': 142, #время работы ступени (с)
    },

    {   #четвертая ступень
        'thrust': 16560, #сила тяги (Н)
        'start_mass': 5623, #масса с топливом (кг)
        'end_mass': 3123, #сухая масса (кг)
        'time': 392, #время работы ступени (с)
    },

]

# константы
a0 = 90  # рандом

I = 280
t1 = 1
t_max = 105 # рандом
R = 6371000
C_d = 0.7  # леха сказал (хз)
e = 2.7
S = 5.56

M1 = 0.29
R1 = 8.314
L = 0.0065

g0 = 9.81

def simulator(stages):
    p0 = 1.25
    t = 0
    vx, vy = 0, 0
    h = 0
    g = g0
    x, y = 0, 0
    ax, ay = 0, 0
    v0 = 0
    h0 = 0
    T = 288
    M = sum([stage['start_mass'] for stage in stages])

    h_v, t_v, v_v, m_v = [], [], [], []
    k = 0

    for stage in stages:
        F_t = stage['thrust']
        M_fuel = stage['start_mass'] - stage['end_mass']
        time = stage['time']
        M_dry = stage['end_mass']
        m = M_fuel / time

        for _ in range(time):
            a = radians(a0 * e ** (-(t / t_max)))
            F_t_x = F_t * np.cos(a)
            F_t_y = F_t * np.sin(a)

            v = v0 + sqrt(ax ** 2 + ay ** 2) * t1

            h = h0 + v * t1


            x1 = x + vx * t1
            y1 = y + vy * t1


            p = p0 * e ** (-M1 * g * h / R1 / T)
            Fy = M * g
            F_s = C_d * p * v ** 2 * S / 2
            F_s_x = F_s * np.cos(a)
            F_s_y = F_s * np.sin(a)
            ax = (F_t_x - F_s_x) * np.sin(a) / M - g * np.cos(a)
            ay = (F_t_x - F_s_x) * np.cos(a) / M - g * np.sin(a)
            vx1 = vx + ax * t1
            vy1 = vy + ay * t1

            g = g0 * (R / (R + h)) ** 2
            M -= m


            h_v.append(h)
            v_v.append(v)
            t_v.append(t)
            m_v.append(M)


            x = x1
            y = y1
            h0 = h
            v0 = v
            p0 = p
            vx = vx1
            vy = vy1

            t += 1

        M -= M_dry

    return np.array(t_v), np.array(h_v), np.array(v_v), np.array(m_v)


t_v, h_v, v_v, m_v = simulator(stages)
k = 180
t_v, h_v, v_v, m_v = t_v[:k], h_v[:k], v_v[:k], m_v[:k]
s = len(t_v)

with open('KSP_time.txt', 'r', encoding='UTF-8') as f:
    t_values = np.array([float(x.rstrip('\n')) for x in f.readlines() if float(x) < t_v[-1]])

with open('KSP_heigh.txt', 'r', encoding='UTF-8') as f:
    h_values = np.array([float(x.rstrip('\n')) for x in f.readlines()])[:len(t_values)]

with open('KSP_speed.txt', 'r', encoding='UTF-8') as f:
    v_values = np.array([float(x.rstrip('\n')) for x in f.readlines()])[:len(t_values)]

with open('KSP_mass.txt', 'r', encoding='UTF-8') as f:
    m_values = np.array([float(x.rstrip('\n')) for x in f.readlines()])[:len(t_values)]



# График 1: Зависимость высоты от времени
plt.figure(figsize=(10, 5))
plt.plot(t_values, h_values, label="KSP")
plt.plot(t_v, h_v, label="math_model")
plt.xlabel("Время (с)")
plt.ylabel("Высота (м)")
plt.title('Сравнение зависимости высоты от времени')
plt.legend()
plt.grid()
plt.show()

# График 2: Зависимость скорости от времени
plt.figure(figsize=(10, 5))
plt.plot(t_values, v_values, label="KSP")
plt.plot(t_v, v_v, label="math_model")
plt.xlabel("Время (с)")
plt.ylabel("Скорость (м/с)")
plt.title('Сравнение зависимости скорости от времени')
plt.legend()
plt.grid()
plt.show()

# График 3: Зависимость массы от времени
plt.figure(figsize=(10, 5))
plt.plot(t_values, m_values, label="KSP")
plt.plot(t_v, m_v, label="math_model")
plt.xlabel("Время (с)")
plt.ylabel("Масса (кг)")
plt.title('Сравнение зависимости массы от времени')
plt.legend()
plt.grid()
plt.show()
