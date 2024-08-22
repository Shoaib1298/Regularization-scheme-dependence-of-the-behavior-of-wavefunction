import numpy as np
from scipy.optimize import bisect
import matplotlib.pyplot as plt
import csv

def find_convergance_point(
    v_n,
    ratio,
    u_inf,
    b,
    h,
    y,
    wavefunction_backward,
    u_0,
    ro_m,
    x,
    wavefunction_forward,
    lamda,
):

    
    v_0 = -1 * v_n[lamda - 8]
    gamma = v_0 / ratio

    E_min = 0 / v_0
    E_max = 5.0 / v_0
    E_trail = E_min

    E_bound = bisect(
        log_difference,
        E_min,
        E_max,
        args=(
            u_inf,
            b,
            h,
            y,
            wavefunction_backward,
            u_0,
            ro_m,
            x,
            wavefunction_forward,
            lamda,
            gamma,
        ),
        xtol=tolerance,
    )

    print("BE:", E_bound * v_0)

    BE = E_bound*v_0

    forward = claculate_wavefunction_backward_complete(
    u_inf, E_trail, ro_m, b, h, y, wavefunction_backward, lamda, gamma)
    wavefunction = list(forward[1])
    wavefunction.reverse()
    sum_wavefunction = sum(wavefunction)
    wavefunction = [i/sum_wavefunction for i in wavefunction]
    x = list(forward[0])
    x.reverse()
    exp =[np.exp(-gamma * E_bound * x) for x in x]

    # finding convergance point
    for i in range(0, len(wavefunction)):
        
        
        if np.abs(wavefunction[i] - exp[i]) <0.005 :
            
            return (x[i],E_bound * v_0)
    # error handling
    return ("no_converg",E_bound * v_0)


def calculate_reducedPotential(r, lamda):
    v = np.float128(np.exp(-lamda * (r**2)))
    return v


def log_difference(
    E_trail,
    u_inf,
    b,
    h,
    y,
    wavefunction_backward,
    u_0,
    ro_m,
    x,
    wavefunction_forward,
    lamda,
    gamma,
):
    # *  returns : x, wavefunction_forward, ro_m(matching point)
    forward = calculate_wavefunction_forward(
        u_0, E_trail, ro_m, h, x, wavefunction_forward, lamda, gamma
    )

    ro_mid = forward[2]

    # * returns : y, wavefunction_backward
    backward = claculate_wavefunction_backward(
        u_inf, E_trail, ro_mid, b, h, y, wavefunction_backward, lamda, gamma
    )

    # * finding logarithimic derivatives at matching point

    der_wavefunction_forward = (forward[1][-1] - forward[1][-2]) / h
    log_der_wavefunction_forward = der_wavefunction_forward / forward[1][-1]

    # there is still a small ambiguity about which points to choose for dervatives

    der_wavefunction_backward = (backward[1][-2] - backward[1][-1]) / h
    log_der_wavfunction_backward = der_wavefunction_backward / backward[1][-2]

    return log_der_wavfunction_backward - log_der_wavefunction_forward


def calculate_u_0(ro, l):
    return (ro) ** (l + 1)


def calculate_u_infinity(ro, gamma, E_trail):  # reduced energy = mod(E)/V_0
    return np.exp(-gamma * E_trail * ro)



def calculate_wavefunction_forward(
    u_0, E_trail, ro_m, h, x, wavefunction_forward, lamda, gamma
):

    # adding initial points in the wavefunction
    r_0 = 0
    wavefunction_current_1 = u_0
    wavefunction_current_2 = h * (l + 1) * (r_0**l) + wavefunction_current_1
    wavefunction_forward[0] = u_0
    wavefunction_forward[1] = wavefunction_current_2

    # important as x =2*h
    i = 2 * h
    # error handling due to appending
    if len(wavefunction_forward) != 2 and len(x) != 2:
        wavefunction_forward = [u_0, wavefunction_current_2]
        x = [0, h]

    while i <= ro_m:
        v = calculate_reducedPotential(i, lamda)
        k = np.float128(gamma * (v - E_trail))
        # finding next wafunction
        wavefunction_current_3 = np.float128(
            (-k * (h**2) + 2) * wavefunction_current_2 - wavefunction_current_1
        )

        # updating for finding the next point wavefunction

        wavefunction_current_1 = wavefunction_current_2
        wavefunction_current_2 = (
            wavefunction_current_3  # now current_2 is the most latest found value
        )
        x.append(i)
        wavefunction_forward.append(wavefunction_current_2)

        if wavefunction_current_2 - wavefunction_current_1 < 0:
            ro_m = i  # chaning ro_m to make sure i > ro_m if this if condition is satisfied

        i = i + h
    return (x, wavefunction_forward, ro_m)


def calculate_wavefunction_forward_complete(
    u_0, E_trail, b, h, x, wavefunction_forward, lamda, gamma
):

    # adding initial points in the wavefunction
    r_0 = 0
    wavefunction_current_1 = u_0
    wavefunction_current_2 = h * (l + 1) * (r_0**l) + wavefunction_current_1
    wavefunction_forward[0] = u_0
    wavefunction_forward[1] = wavefunction_current_2

    # important as x =2*h
    i = 2 * h
    # error handling due to appending
    if len(wavefunction_forward) != 2 and len(x) != 2:
        wavefunction_forward = [u_0, wavefunction_current_2]
        x = [0, h]

    while i <= b:
        v = calculate_reducedPotential(i, lamda)
        try:
            k = np.float128(gamma * (v - E_trail))
        except TypeError:
            k = np.float128(gamma * (-E_trail))

        # finding next wafunction
        wavefunction_current_3 = np.float128(
            (-k * (h**2) + 2) * wavefunction_current_2 - wavefunction_current_1
        )

        # updating for finding the next point wavefunction

        wavefunction_current_1 = wavefunction_current_2
        wavefunction_current_2 = (
            wavefunction_current_3  # now current_2 is the most latest found value
        )
        x.append(i)
        wavefunction_forward.append(wavefunction_current_2)
        i = i + h

    return (x, wavefunction_forward, ro_m)


def claculate_wavefunction_backward(
    u_inf, E_trail, ro_m, b, h, y, wavefunction_backward, lamda, gamma
):

    u_inf = calculate_u_infinity(b, gamma, E_trail)
    wavefunction_current_1 = u_inf
    wavefunction_current_2 = np.float128(
        gamma * E_trail * h * np.exp(-gamma * E_trail * b) + u_inf
    )
    wavefunction_backward = [wavefunction_current_1, wavefunction_current_2]
    y = [b, b - h]

    j = b - 2 * h
    while j >= ro_m:
        v = calculate_reducedPotential(j, lamda)
        k = np.float128(gamma * (v - E_trail))

        wavefunction_current_3 = (
            (2 - (h**2) * k) * wavefunction_current_2 -wavefunction_current_1 
        )

        # updating

        wavefunction_current_1 = wavefunction_current_2
        wavefunction_current_2 = wavefunction_current_3

        y.append(j)
        j = j - h
        wavefunction_backward.append(wavefunction_current_2)

    return (y, wavefunction_backward)

def claculate_wavefunction_backward_complete(
    u_inf, E_trail, ro_m, b, h, y, wavefunction_backward, lamda, gamma
):

    u_inf = calculate_u_infinity(b, gamma, E_trail)
    wavefunction_current_1 = u_inf
    wavefunction_current_2 = np.float128(
        gamma * E_trail * h * np.exp(-gamma * E_trail * b) + u_inf
    )
    wavefunction_backward = [wavefunction_current_1, wavefunction_current_2]
    y = [b, b - h]

    j = b - 2 * h
    while j >= 0:
        v = calculate_reducedPotential(j, lamda)
        k = np.float128(gamma * (v - E_trail))

        wavefunction_current_3 = (
            (2 - (h**2) * k) * wavefunction_current_2 -wavefunction_current_1 
        )

        # updating

        wavefunction_current_1 = wavefunction_current_2
        wavefunction_current_2 = wavefunction_current_3

        y.append(j)
        j = j - h
        wavefunction_backward.append(wavefunction_current_2)

    return (y, wavefunction_backward)
# main():

# b = 150
# # range[0,b]

# h = 10 ** (-5)
# # step size

ro_m = 0.8
# breaking point

tolerance = 1e-7

# constants

z = 0.5 * (938.211 + 939.505)
ratio = 197.31613 * 197.31613 / z

l = 0  # always

# intialisatons

u_0 = calculate_u_0(h, l)

wavefunction_forward = [0, 0]
x = [0, h]

u_inf = 1

y = [b, b - h]
wavefunction_backward = [0, 0]

# Energy range

v_n = [
    -974.44358871,
    -1090.78155295,
    -1206.87466956,
    -1322.66304077,
    -1438.39560818,
    -1553.96724232,
    -1668.97275703,
    -1784.88780437,
    -1900.64114357,
    -2014.10843202,
    -2129.00782100,
    -2243.69305358,
    -2358.72004282,
    -2472.99271987,
    -2586.82165122,
    -2701.28088046,
    -2815.62635980,
    -2931.31770607,
    -3044.80940351,
    -3158.37040404,
    -3275.43736806,
    -3387.01440991,
    -3500.29557543,
    -3615.82742806,
    -3728.48171651,
    -3842.84847331,
    -3956.24919463,
    -4070.52288747,
    -4185.96872588,
    -4298.41934049,
    -4412.59978550,
    -4525.86795233,
    -4639.58072838,
    -4753.17016559,
    -4866.74260853,
]


converging_ro = []


for lamda in range(8, 43):
    b = random.uniform(100, 120)
    h = random.uniform(10**-6, 10**-5)
    ro_c,BE = find_convergance_point(
        v_n,
        ratio,
        u_inf,
        b,
        h,
        y,
        wavefunction_backward,
        u_0,
        ro_m,
        x,
        wavefunction_forward,
        lamda,
    )

    print(ro_c)

    with open("ro_c.csv","a") as file :
        writer = csv.DictWriter(file,fieldnames=["BE","ro_c"])
        writer.writerow({"BE": BE ,"ro_c":ro_c})

    
    converging_ro.append(ro_c)


# graphs
# forward = calculate_wavefunction_forward(u_0,1.9/v_0,ro_m,h,x,wavefunction_forward)
# bacward =calulate_u_backward(u_inf,1.9/v_0,ro_m,b,h,y,wavefunction_backward)
# plt.plot(forward[0],forward[1])
# plt.plot(bacward[0],bacward[1])
# ]
"""
New output :

BE: 2.2325897216796875
50.05338000468811
BE: 2.236099243164063
45.750970003322465
BE: 2.240219116210937
no_converg
BE: 2.2408294677734375
47.679090003934476
BE: 2.2473907470703125
64.64120000931848
BE: 2.2541809082031254
45.03095000309392
BE: 2.243804931640625
no_converg
BE: 2.2749328613281246
48.09771000406735
BE: 2.3027038574218754
no_converg
BE: 2.2499084472656254
44.498130002924796
BE: 2.2541809082031254
49.630140004553766
BE: 2.2541809082031246
no_converg
BE: 2.2688293457031246
44.999930003084074
BE: 2.260589599609375
44.43585000290503
BE: 2.240753173828125
no_converg
BE: 2.2438049316406254
no_converg
BE: 2.245330810546875
no_converg
BE: 2.289276123046875
no_converg
BE: 2.2682189941406254
no_converg
BE: 2.25128173828125
no_converg
BE: 2.3373413085937504
no_converg
BE: 2.2671508789062504
42.51358000229487
BE: 2.24761962890625plt.plot(forward[0],forward[1])
43.55095000262415
BE: 2.2909545898437496
48.59268000422446
BE: 2.25738525390625
43.170580002503414
BE: 2.27081298828125
43.94847000275033
BE: 2.2598266601562504
43.865830002724095
BE: 2.27203369140625
no_converg
BE: 2.31536865234375
45.75272000332302
BE: 2.2830200195312496
44.261310002849626
BE: 2.2946166992187496
no_converg
BE: 2.2848510742187504
no_converg
BE: 2.2872924804687496
44.227470002838885
BE: 2.28668212890625
no_converg
BE: 2.28668212890625
50.30243000476716

"""

