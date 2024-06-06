import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from StarChart import StarChart as Chart

"""
基于matplot与矢量力学的太阳系天体运行程序
"""

sun = {
    "mass": 1.9885e30,
    "position": np.array([0.0, 0.0]),
    "v": np.array([0.0, 0.0])
}
earth_mass = 5.972e24  # 地球质量
earth_T = 365.242 * 86400  # 地球周期 = 1年
AU = 1.496e11  # 天文单位 = 地球轨道半径
earth_v = 2 * np.pi * AU / earth_T
planets_mass = [
    0.0553, 0.8150, 1,  # 水星, 金星, 地球
    0.1074, 317.94, 95.18,  # 火星, 木星, 土星
    14.63, 17.22  # 天王, 海王
]
planets_radial = [  # 轨道半长轴
    0.38, 0.72, 1,  # 水星, 金星, 地球
    1.52, 5.20, 9.54,  # 火星, 木星, 土星
    19.22, 30.06,  # 天王, 海王
]
planets_e = [   # 轨道偏心率
    0.2056, 0.0068, 0.0167,
    0.0934, 0.0483, 0.0560,
    0.0461, 0.0097
]
planets_period = []
planets_initial_vec = []
for i in range(len(planets_radial)):  # 逐一计算轨道周期和初始线速度
    r_p = planets_radial[i]
    T_p = 1 * np.sqrt(r_p * r_p * r_p)  # 开普勒第三定律
    planets_period.append(T_p)
    v_p = r_p / T_p  # 依旧以地球速度为归一化单位
    planets_initial_vec.append(v_p * np.sqrt(2 / (1 + planets_e[i]) - 1))

planets_names = [
    "Mercury", "Venus", "Earth",
    "Mars", "Jupiter", "Saturn",
    "Uranus", "Neptune"
]

label = ["sun",
         "Mercury", "Venus", "Earth",
         "Mars", "Jupiter", "Saturn",
         "Uranus", "Neptune"
         ]

color = [
    "red",  # 太阳
    "gold", "yellow", "blue",  # 水, 金, 地
    "firebrick", "goldenrod", "darkgoldenrod",  # 火, 木, 土
    "skyblue", "steelblue"
]

if __name__ == '__main__':
    solarChart = Chart(label=label, color=color)
    solarChart.addBodies(sun["mass"], sun["position"], sun["v"])  # 置入太阳
    for i in range(len(planets_names)):  # 置入行星
        solarChart.addBodies(planets_mass[i] * earth_mass,  # 行星质量
                             np.array([planets_radial[i] * AU * (1+planets_e[i]), 0.0]),  # 远日点距离
                             np.array([0.0, planets_initial_vec[i] * earth_v])   # 行星初速度
                             )
    solarChart.setRefLength(AU)  # 比例尺
    fig, a = plt.subplots()  # 用于绘图
    ax = a
    for body in solarChart.bodies:
        body.initPlotLines(ax)

    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.legend()
    ax.set_aspect('equal', adjustable='box')
    # 初始化动画
    ani = animation.FuncAnimation(fig, solarChart.animationUpdate, fargs=tuple([0.01]), frames=range(10),
                                  interval=1, blit=True)

    plt.ion()
    plt.pause(100)
