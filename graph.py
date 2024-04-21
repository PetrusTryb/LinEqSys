from sys import argv
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
	with open("B_jacobi.txt", "r") as f:
		data_j = f.readlines()
		data_j = [float(x.strip()) for x in data_j]
	with open("B_gauss_seidel.txt", "r") as f:
		data_gs = f.readlines()
		data_gs = [float(x.strip()) for x in data_gs]
	print(data_j)
	print(data_gs)
	plt.yscale("log")
	plt.plot(data_j, label="Jacobi")
	plt.plot(data_gs, label="Gauss-Seidel")
	plt.legend()
	plt.xlabel("Iteration")
	plt.ylabel("Residual")
	plt.xticks(np.arange(0, max(len(data_j),len(data_gs)), 1))
	plt.savefig("B.png")
	plt.close()
 
	with open("C_jacobi.txt", "r") as f:
		data_j = f.readlines()
		data_j = [float(x.strip()) for x in data_j]
	with open("C_gauss_seidel.txt", "r") as f:
		data_gs = f.readlines()
		data_gs = [float(x.strip()) for x in data_gs]
	print(data_j)
	print(data_gs)
	plt.yscale("log")
	plt.plot(data_j, label="Jacobi")
	plt.plot(data_gs, label="Gauss-Seidel")
	plt.legend()
	plt.xlabel("Iteration")
	plt.ylabel("Residual")
	plt.xticks(np.arange(0, max(len(data_j),len(data_gs)), 4))
	plt.savefig("C.png")
	plt.close()
	
	with open("E_jacobi.txt", "r") as f:
		data_j = f.readlines()
		data_j = [float(x.strip()) for x in data_j]
	with open("E_gauss_seidel.txt", "r") as f:
		data_gs = f.readlines()
		data_gs = [float(x.strip()) for x in data_gs]
	with open("E_lu.txt", "r") as f:
		data_lu = f.readlines()
		data_lu = [float(x.strip()) for x in data_lu]
	print(data_j)
	print(data_gs)
	print(data_lu)
	x_axis = [100, 500, 1000, 2000, 3000, 5000]
	plt.plot(x_axis, data_j, label="Jacobi")
	plt.plot(x_axis, data_gs, label="Gauss-Seidel")
	plt.plot(x_axis, data_lu, label="LU")
	plt.legend()
	plt.xlabel("Matrix Size")
	plt.ylabel("Time (s)")
	plt.xticks(x_axis)
	plt.yticks(np.arange(0, max(max(data_j),max(data_gs),max(data_lu)), 100))
	plt.savefig("E.png")
	plt.close()
	