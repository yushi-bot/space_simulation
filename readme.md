# 天体模拟
# Space Simulation
通过python代码模拟天体运行，只使用了简单的牛顿力学模型和简单的欧拉积分
## 模拟结果

## simulation results
### 地月系模拟
<div align="center"> <img src="/results/地月拉格朗日点模拟结果.gif" width = 400 height = 400 /> </div>
<br>该模拟模拟了地月系的运行，图中月球绕地一圈相当于1个月的时间。
在绕地一定时间后，向地月系的L1至L5点附近添加小型天体，继续模拟该系统的行为(可认为是拉格朗日空间站)
其中L1点空间站运行数圈后脱离原有轨道
L2点空间站在约10天后坠入月球
L3点空间站稳定运行数月，但最终会仍会脱离轨道,并坠入月球
L4点和L5点空间站稳定运行
与理论预测结果相一致

### 内太阳系模拟
<div align="center"> <img src="/results/内太阳系模拟.gif" width = 400 height = 400 /> </div>
该模拟模拟了太阳系及其主要行星（未考虑偏心率的圆形轨道），图中地球绕太阳一圈即为1年时间。
图中长度单位为AU(天文单位)