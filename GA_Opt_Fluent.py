import matplotlib.pyplot as plt
import numpy as np
import geatpy as ea
import Problem_GA_Opt_Fluent
from Problem_GA_Opt_Fluent import MyProblem

NIND=Problem_GA_Opt_Fluent.NIND

heat_flux_array = np.zeros((NIND, 309))
heat_source=np.zeros((NIND, 7)) #7个分区，NIND种群规模
Min=np.zeros((NIND,309))#流入质量
Mout=np.zeros((NIND,309))#流入质量



if __name__=='__main__':
	problem=MyProblem()#实例化问题对象
	Encoding='RI'#编码方式
	NIND=100 #种群规模
	Field=ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders)#区域描述器
	population=ea.Population(Encoding, Field, NIND)#实例化种群对象
	myAlgorithm = ea.soea_DE_rand_1_bin_templet(problem, population)#实例化算法模板对象
	myAlgorithm.MAXGEN=5000#最大遗传代数
	myAlgorithm.mutOper.F = 0.8  #差分进化中的参数F，变异缩放因子
	myAlgorithm.recOper.XOVR = 0.6  #重组概率，交叉概率
	myAlgorithm.logTras = 1  #设置每隔多少代记录日志，若设置成0则表示不记录日志
	myAlgorithm.verbose = True  #设置是否打印输出日志信息
	myAlgorithm.drawing = 2  #设置绘图方式（0：不绘图；1：绘制结果图；2：绘制目标空间过程动画；3：绘制决策空间过程动画）
	[BestIndi, population, population2]=myAlgorithm.run()#执行算法模板，得到最优个体以及最后一代种群
	BestIndi.save()#把最优个体的信息保存到文件中，BestIndi是population类
	"""==========================输出结果==========================="""
	print('评价次数：%s' % myAlgorithm.evalsNum)
	print('时间已过 %s 秒' % myAlgorithm.passTime)
	if BestIndi.sizes != 0:
		print('最优的目标函数值为：%s' % BestIndi.ObjV[0][0])#BestIndi:最优个体
		# print('最优的控制变量值为：')
		# for i in range(BestIndi.Phen.shape[1]):#k.shape[0]输出矩阵k的行数,k.shape[1]输出矩阵k的列数。Phen种群表现型矩阵。
		# 	print(BestIndi.Phen[0, i])#输出最后一代种群的第一个个体，也就是最后一代种群表现型的第一行
		heat_source=BestIndi.Phen
		for i in range(0, 7):
			heat_source[:, [i]] = BestIndi.Phen[:, [i]]  # 每一行是一条染色体，也就是一个个体，一个具体问题的解


		for i in range(0,309):
			heat_flux_array[:,[i]]=three_cell_test.aim(x[:, [i]], Mout[:,[i]], Min[:,[i]], area[:,[i]],h_air[:,[i]], pressure[:,[i]], beta[:,[i]])
		np.savetxt("Result/heat_flux_array.csv", heat_flux_array, delimiter=',')
		np.savetxt("Result/Mout.csv", Mout, delimiter=',')
		print('中间网格热流密度: ', heat_flux_array[0, 150])
		print('后两个网格温度差: ', x[0, 307]-x[0,308])
		print('最后网格流出质量: ', Mout[0, 308])
	else:
		print('没找到可行解。')