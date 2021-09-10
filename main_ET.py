# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 20:56:16 2021

@author: guoxf
"""
import matplotlib.pyplot as plt
import numpy as np
import geatpy as ea
import problem_ET
from problem_ET import MyProblem


if __name__ == '__main__':
	problem=MyProblem()
	Encoding='RI'
	NIND=problem_ET.NIND
	Field=ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders)#区域描述器
	population=ea.Population(Encoding, Field, NIND)#实例化种群对象
	myAlgorithm = ea.soea_DE_rand_1_bin_templet(problem, population)#实例化算法模板对象
	myAlgorithm.MAXGEN=200#最大遗传代数
	myAlgorithm.mutOper.F = 0.6  #差分进化中的参数F，变异缩放因子
	myAlgorithm.recOper.XOVR = 0.4  #重组概率，交叉概率
	myAlgorithm.logTras = 1  #设置每隔多少代记录日志，若设置成0则表示不记录日志
	myAlgorithm.verbose = True  #设置是否打印输出日志信息
	myAlgorithm.drawing = 2  #设置绘图方式（0：不绘图；1：绘制结果图；2：绘制目标空间过程动画；3：绘制决策空间过程动画）

	[BestIndi, obj_trace, var_trace]=myAlgorithm.run()#执行算法模板，得到最优个体以及最后一代种群
	BestIndi.save()#把最优个体的信息保存到文件中，BestIndi是population类
	"""==========================输出结果==========================="""
	print('评价次数：%s' % myAlgorithm.evalsNum)
	print('时间已过 %s 秒' % myAlgorithm.passTime)
	if BestIndi.sizes != 0:
		print('最优的目标函数值为：%s' % BestIndi.ObjV[0][0])#BestIndi:最优个体
		print('最优的控制变量值为：')
		for i in range(BestIndi.Phen.shape[1]):#k.shape[0]输出矩阵k的行数,k.shape[1]输出矩阵k的列数。Phen种群表现型矩阵。
			print(BestIndi.Phen[0, i])#输出最后一代种群的第一个个体，也就是最后一代种群表现型的第一行
	else:
    		print('没找到可行解。')