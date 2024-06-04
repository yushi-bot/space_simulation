import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from StarChart import StarChart as Chart

"""
基于matplot与矢量力学的拉格朗日点天体运行程序
将下方半径和质量，周期等参数换为其他天体(比如日地)
L4和L5也成立，L123三个点则需要重新算出坐标
"""

T = 27.32 * 86400  # 地月公转周期，单位秒
R = 3.84e8  # 月球轨道近似为圆形，平均轨道半径
v = 2 * np.pi / T * R  # 月球绕地线速度
earth = {  # 地球参数
    "mass": 5.972e24,
    "position": np.array([0, 0]),
    "v": np.array([0, 0]),
}
moon = {
    "mass": 7.342e22,
    "position": np.array([R, 0]),
    "v": np.array([0, v]),
}
Lagrangian_points = {  # 拉格朗日点坐标
    "L1": {
        "loc": np.array([(1-0.1596003)*R, 0]),
        "v":  np.array([0, 2 * np.pi / T * (1-0.1596003)*R])
    },
    "L2": {
        "loc": np.array([(1+0.1595926)*R, 0]),
        "v": np.array([0, 2 * np.pi / T * (1+0.1595926)*R])
    },
    "L3": {
        "loc": np.array([-0.992886*R, 0]),
        "v": np.array([0, 2 * np.pi / T * -3.816e8])
    },
    "L4": {
        "loc":np.array([R*0.5, R*1.732/2]),
        "v":np.array([-v*1.732/2, v*0.5])
    },
    "L5": {
        "loc": np.array([R*0.5, -R*1.732/2]),
        "v": np.array([v*1.732/2, v*0.5])
    }
}
LMass = 1e6  # 置入该点的天体质量，足够小则可认为对行星无干扰

label = ["earth", "moon", "L1", "L2", "L3", "L4", "L5"]

if __name__ == '__main__':
    myChart = Chart(label=label)
    myChart.addBodies(earth["mass"], earth["position"], earth["v"])  # 把地球放进星图
    myChart.addBodies(moon["mass"], moon["position"], moon["v"])  # 把月球放进星图

    myChart.setRefLength(R)  # 用轨道半径作参考长度

    fig, a = plt.subplots()     # 用于绘图
    ax = a
    for body in myChart.bodies:
        body.initPlotLines(ax)

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.legend()
    ax.set_aspect('equal', adjustable='box')
    # 初始化动画
    ani = animation.FuncAnimation(fig, myChart.animationUpdate, fargs=tuple([0.01]), frames=range(10),
                                  interval=1, blit=True)
    #plt.show()

        
        
    for i in range(1, 6):   # 逐个加入对应的拉格朗日空间站
        plt.ion()
        plt.pause(1)


        Lx = label[i+1]

        # 地月距离
        direction = myChart.bodies[1].location
        dis = np.sqrt(direction[0]*direction[0] + direction[1]*direction[1])
        # 地月方向
        direction = direction/dis

        # 坐标变换，算出该时刻空间站的速度和位置
        lp_loc = Lagrangian_points[Lx]["loc"]
        lp_v = Lagrangian_points[Lx]["v"]

        loc = np.array([direction[0] * lp_loc[0] + -direction[1] * lp_loc[1],
                       direction[1] * lp_loc[0] + direction[0] * lp_loc[1]])
        vec = np.array([direction[0] * lp_v[0] + -direction[1] * lp_v[1],
                       direction[1] * lp_v[0] + direction[0] * lp_v[1]])

        myChart.addBodies(LMass, loc, vec, checkRef=True)
        myChart.bodies[-1].initPlotLines(ax)
        ax.legend()

    plt.ion()
    plt.pause(100)



