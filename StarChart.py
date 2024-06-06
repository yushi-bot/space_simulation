import numpy as np

'''
星图类，用于存储星图上的天体
'''

tol = 1e-6
aTol = 10
colors = ['b', 'y', 'r', 'g', 'c', 'm', 'k']


class StarChart:
    """
    星图类，定义了天体，并用于计算天体位置、速度、加速度等
    """

    class Body:
        """
        天体类，定义了天体的坐标，质量，速度和加速度。
        """

        def __init__(self, mass, chart, index, initial_loc=np.array([0.0, 0.0]),
                     initial_vel=np.array([0.0, 0.0])):
            self.mass = mass  # 质量
            self.location = initial_loc  # 初位置
            self.velocity = initial_vel  # 初速度
            self.chart = chart  # 星图
            self.index = index  # 星图中的索引值
            self.color = chart.color[index]  # 绘图的颜色
            self.trail = [[], []]  # 轨迹线
            self.acceleration = np.array([0, 0])  # 加速度
            self.scatter = None  # 散点(当前位置)
            self.trace = None  # 轨迹曲线(plt.line2d类)
            self.broken = False  # 是否坠毁
            self.planet = self  # 坠毁后认为跟随行星一起运动

        def calcDistance(self, otherBody) -> list[float]:
            """
            计算距离和方向
            :param otherBody:  另一天体对象
            :return: [距离标量, 方向矢量]
            """
            direction = otherBody.location - self.location
            dis = np.sqrt(direction[0] * direction[0] + direction[1] * direction[1])
            direction = direction / dis
            return [dis, direction]

        def initPlotLines(self, ax):
            """
            初始化曲线
            :param ax: 动画参数
            :return:
            """
            self.scatter = ax.plot(self.location[0], self.location[1], self.color,
                                   marker='o',
                                   markersize=(20 + 0.8 * np.log10(self.mass)))[-1]
            self.trace = ax.plot([self.location[0]], [self.location[1]], self.color)[-1]
            if len(self.chart.label) > self.index:
                self.scatter.set_label(self.chart.label[self.index])

        def updatePlotLines(self):
            """
            更新轨迹和当前位置
            :return:
            """
            if not self.broken:
                self.scatter.set_data([self.location[0]], [self.location[1]])  # 更新点的位置
                self.trace.set_data(self.trail[0][-3000:], self.trail[1][-3000:])  # 更新轨迹
            if self.broken:
                self.scatter.set_data([], [])  # 撞毁后
                self.trace.set_data(self.trail[0][-3000:], self.trail[1][-3000:])  # 撞毁后跟着星星走

    def __init__(self, label=None, color=None, ref_length=1, refTime=1, refMass=1):
        if color is None:
            color = colors
        self.label = label
        self.color = color
        self.bodies = []  # 星图中的天体
        self.ref_length = ref_length  # 星图的长度比例尺 图中单位1 等于 ref_length 的长度
        self.refTime = refTime  # 星图参考时间 1现实时间 = refTime 模拟时间
        self.refMass = refMass  # 星图参考质量 星图内1的质量 = refMass kg 的质量
        self.G = 6.67e-11 / (ref_length * ref_length * ref_length) * (refTime * refTime) * refMass  # 引力常数
        self.time = 0  # 计时器
        self.reference = None  # 参考天体，假定该天体为惯性系(忽略惯性力)
        pass

    def setRefLength(self, ref_length):
        """
        长度比例尺，并根据该比例尺计算参考时间和参考质量
        :param ref_length: 比例尺长度 图中单位1 等于 ref_length 的长度
        :return:
        """
        self.ref_length = ref_length
        self.refTime = None
        self.refMass = None
        for body in self.bodies:

            body.location = body.location / ref_length
            dis = np.sqrt(body.location[0] * body.location[0] + body.location[1] * body.location[1])
            at_ref_point = -0.01 < dis - 1 < 0.01
            v = body.velocity
            if self.refTime is None or at_ref_point:  # 寻找一个参考时间系数，以令其中某一个天体的线速度为1,如果参考距离有天体，就用这个，找不到就随便找一个
                vv = np.sqrt(v[0] * v[0] + v[1] * v[1])  # 速度的模
                if vv < tol * 10:
                    continue
                self.refTime = ref_length / vv

        for body in self.bodies:
            body.velocity = body.velocity * (self.refTime / self.ref_length)

        # 寻找一个参考质量，以令其中某一个天体的质量为1(一般是最大的)
        if self.refMass is None:
            refMass = 0
            for body in self.bodies:
                if body.mass > refMass:
                    refMass = body.mass

            self.refMass = refMass

        for body in self.bodies:
            body.mass /= self.refMass
        self.G = 6.67e-11 / (self.ref_length * self.ref_length * self.ref_length) * \
                 self.refTime * self.refTime * self.refMass  # 引力常数

    def addBodies(self, mass, initial_loc=np.array([0.0, 0.0]), initial_vel=np.array([0.0, 0.0]), checkRef=False):
        """
        向星图中增加天体
        :param checkRef: 检查参考质量、参考坐标参考时间等因素
        :param mass: 质量
        :param initial_loc: 初始坐标
        :param initial_vel: 初始速度
        :return:
        """
        count = len(self.bodies)
        newBody = StarChart.Body(mass, self, count, initial_loc, initial_vel)
        if checkRef:
            newBody.location = newBody.location / self.ref_length
            newBody.velocity = newBody.velocity * (self.refTime / self.ref_length)
            newBody.mass = newBody.mass / self.refMass
        self.bodies.append(newBody)
        if self.reference is None:  # 如果尚未确定参考天体，以第一个加入的天体为参考系
            self.reference = self.bodies[0]

    def clone(self):
        """
        拷贝该星图
        :return: 该星图的拷贝(新内存)
        """
        res = StarChart(label=self.label,
                        color=self.color,
                        ref_length=self.ref_length,
                        refTime=self.refTime,
                        refMass=self.refMass)
        for body in self.bodies:
            res.addBodies(body.mass, body.location, body.velocity)
        return res

    def calcA(self):
        """
        计算每一个天体的加速度
        :return:
        """
        count = len(self.bodies)
        for i in range(count):

            body1 = self.bodies[i]
            if body1.broken:  # 已摧毁
                continue
            mass1 = body1.mass

            for j in range(i + 1, count):  # 计算每一个天体受每一个天体的引力加速度
                body2 = self.bodies[j]
                if body2.broken:  # 已摧毁
                    continue
                mass2 = body2.mass
                [dis, direction] = body1.calcDistance(body2)
                a1 = self.G * mass2 / dis / dis  # 计算body2 对 body1 产生的加速度
                a2 = self.G * mass1 / dis / dis  # 计算body1 对 body2 产生的加速度

                if a1 > tol:  # 忽略小于tol的加速度(即小天体的作用)
                    a1_vec = direction * a1
                    body1.acceleration = body1.acceleration + a1_vec

                if a2 > tol:  # 忽略小于tol的加速度(即小天体的作用)
                    a2_vec = -direction * a2
                    body2.acceleration = body2.acceleration + a2_vec

    def RKIntegrate(self, dt):
        """
        4阶龙格库塔积分，效果不好，不要用
        :param dt: 时间步长
        :return:
        """
        chart0 = self.clone()  # 初始星图 保险起见复制一份
        chart0.calcA()
        v0 = []
        loc0 = []
        a0 = []
        for body in chart0.bodies:
            v0.append(body.velocity)
            a0.append(body.acceleration)
            loc0.append(body.location)

        chart1 = self.clone()  # 根据初始星图的v, a预测得到的dt/2时间后的假想星图，星图1

        v1 = []
        a1 = []
        loc1 = []
        for bodyIndex in range(len(chart1.bodies)):  # 预测dt/2时间后的位置、速度、加速度
            body = chart1.bodies[bodyIndex]
            body.location = body.location + (dt * 0.5) * v0[bodyIndex] + 0.5 * (dt * 0.5) * (dt * 0.5) * a0[bodyIndex]
            body.velocity = body.velocity + a0[bodyIndex] * (dt * 0.5)
            loc1.append(body.location)
            v1.append(body.velocity)

        chart1.calcA()
        for body in chart1.bodies:
            a1.append(body.acceleration)

        chart2 = self.clone()  # 根据星图1的v, a预测得到的dt/2时间后的假想星图2

        v2 = []
        a2 = []
        loc2 = []
        for bodyIndex in range(len(chart2.bodies)):  # 预测dt/2时间后的位置、速度、加速度
            body = chart2.bodies[bodyIndex]
            body.location = body.location + (dt * 0.5) * v1[bodyIndex] + 0.5 * (dt * 0.5) * (dt * 0.5) * a1[bodyIndex]
            body.velocity = body.velocity + a1[bodyIndex] * (dt * 0.5)
            loc2.append(body.location)
            v2.append(body.velocity)

        chart2.calcA()
        for body in chart1.bodies:
            a2.append(body.acceleration)

        chart3 = chart1.clone()  # 根据星图2的v, a预测得到的dt时间后的假想星图3

        v3 = []
        a3 = []
        loc3 = []
        for bodyIndex in range(len(chart2.bodies)):  # 预测dt/2时间后的位置、速度、加速度
            body = chart2.bodies[bodyIndex]
            body.location = body.location + (dt * 0.5) * v2[bodyIndex] + 0.5 * (dt * 0.5) * (dt * 0.5) * a2[bodyIndex]
            body.velocity = body.velocity + a2[bodyIndex] * (dt * 0.5)
            v3.append(body.velocity)
            loc3.append(body.location)
        chart3.calcA()
        for body in chart1.bodies:
            a3.append(body.acceleration)

        v = (np.array(v0) + 2 * np.array(v1) + 2 * np.array(v2) + np.array(v3)) / 6
        loc = (np.array(loc0) + 2 * np.array(loc1) + 2 * np.array(loc2) + np.array(loc3)) / 6
        a = (np.array(a0) + 2 * np.array(a1) + 2 * np.array(a2) + np.array(a3)) / 6

        return v, loc

    def eulerIntegrate(self, dt):
        """
        改进欧拉积分
        :param dt: 时间步长
        :return:
        """
        chart0 = self.clone()  # 初始星图 保险起见复制一份
        chart0.calcA()
        v0 = []
        loc0 = []
        a0 = []
        for body in chart0.bodies:
            v0.append(body.velocity)
            a0.append(body.acceleration)
            loc0.append(body.location)

        chart1 = self.clone()  # 星图1 保险起见复制一份
        v1 = []
        loc1 = []
        for bodyIndex in range(len(chart1.bodies)):  # 预测dt时间后的位置、速度、加速度
            body = chart1.bodies[bodyIndex]
            body.location = body.location + dt * v0[bodyIndex] + 0.5 * dt * dt * a0[bodyIndex]
            body.velocity = body.velocity + a0[bodyIndex] * dt
            loc1.append(body.location)
            v1.append(body.velocity)

        # 计算更新后的速度和位置
        v = (np.array(v0) + np.array(v1)) / 2
        loc = (np.array(loc0) + np.array(loc1)) / 2
        return v, loc

    def updatePos(self, dt):
        """
        更新天体位置
        :param dt:时间增量步
        :return:
        """
        count = len(self.bodies)
        # v_vec, loc_vec = self.RKIntegrate(dt)
        v_vec, loc_vec = self.eulerIntegrate(dt)

        for k in range(count):  # 更新每一个天体的速度和位置
            body1 = self.bodies[k]
            if body1.broken:
                continue

            v = v_vec[k]
            loc = loc_vec[k]
            body1.location = loc
            body1.velocity = v

            if body1.index == self.reference.index:  # 参考天体的坐标和速度视为0
                continue
            body1.location -= self.reference.location
            body1.velocity -= self.reference.velocity
            # 更新轨迹
            body1.trail[0].append(body1.location[0])
            body1.trail[1].append(body1.location[1])

        self.reference.location = np.array([0, 0])  # 参考天体的坐标和速度视为0
        self.reference.velocity = np.array([0, 0])  # 参考天体的坐标和速度视为0
        for i in range(count):  # 检查空间站是否坠毁
            body1 = self.bodies[i]
            if body1.broken:  # 坠毁后跟着行星走
                body1.trail[0].append(body1.planet.location[0])
                body1.trail[1].append(body1.planet.location[1])
                continue

            for body2 in self.bodies:
                if body2.index == body1.index or body2.broken:
                    continue
                mass2 = body2.mass
                [dis, direction] = body1.calcDistance(body2)
                a1 = self.G * mass2 / dis / dis  # 计算body2 对 body1 产生的加速度
                if a1 > aTol:  # 受到加速度超过容许上限就认为body1坠毁了
                    body1.broken = True
                    print("L" + str(body1.index - 1) + " broken")
                    body1.planet = body2
                    break

    def animationUpdate(self, T, dt):
        """
        更新动画函数
        :param T: 当前帧数，没用到
        :param dt: 时间步长
        :return:
        """
        self.time += dt

        self.updatePos(dt)
        if self.time >= 8 * np.pi:
            self.time -= 8 * np.pi

        animateList = []
        for body in self.bodies:
            body.updatePlotLines()
            animateList.append(body.scatter)
            animateList.append(body.trace)

        return tuple(animateList)
