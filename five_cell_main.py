import numpy as np
import geatpy as ea
from three_cell_test import MyProblem

if __name__=='__main__':
	problem=MyProblem()
	Encoding='RI'
	NIND=400
	Field=ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders)
	population=ea.Population(Encoding, Field, NIND)
	myAlgorithm = ea.soea_DE_rand_1_bin_templet(problem, population)
	myAlgorithm.MAXGEN=100
	myAlgorithm.mutOper.F = 0.5  # 差分进化中的参数F
	myAlgorithm.recOper.XOVR = 0.7  # 重组概率
	myAlgorithm.logTras = 1  # 设置每隔多少代记录日志，若设置成0则表示不记录日志
	myAlgorithm.verbose = True  # 设置是否打印输出日志信息
	myAlgorithm.drawing = 2  # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制目标空间过程动画；3：绘制决策空间过程动画）
	[BestIndi, population]=myAlgorithm.run()
	BestIndi.save()

	print('评价次数：%s' % myAlgorithm.evalsNum)
	print('时间已过 %s 秒' % myAlgorithm.passTime)
	if BestIndi.sizes != 0:
		print('最优的目标函数值为：%s' % BestIndi.ObjV[0][0])
		print('最优的控制变量值为：')
		for i in range(BestIndi.Phen.shape[1]):
			print(BestIndi.Phen[0, i])
	else:
		print('没找到可行解。')