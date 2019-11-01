import math
import pylab
import numpy
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
from matplotlib.colors import LogNorm
import math

DPI = 200

if __name__ == '__main__':
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	x = []
	y = []
	z = []
	with open('out.txt') as f:
		for line in f:
			x1, y1, z1 = [float(x) for x in line.split()]
			x.append(x1)
			y.append(y1)
			z.append(math.log(z1))

	last = x[-1]
	with open('out3.txt') as f:
		for line in f:
			x1, y1, z1 = [float(x) for x in line.split()]
			x.append(x1 + last)
			y.append(y1)
			z.append(math.log(z1))

	last = x[-1]
	with open('out2.txt') as f:
		for line in f:
			x1, y1, z1 = [float(x) for x in line.split()]
			x.append(x1 + last)
			y.append(y1)
			z.append(math.log(z1))

	x_list = np.array(x)
	y_list = np.array(y)
	z_list = np.array(z)

	#z = z_list.reshape(11, 31)
	#plt.imshow(z, extent=(np.amin(x_list), np.amax(x_list), np.amin(y_list), np.amax(y_list)), norm=LogNorm(), aspect = 'auto')
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_trisurf(x_list, y_list, z_list, color='white', edgecolors='grey', alpha=0.5,  cmap=cm.coolwarm)
	#plt.colorbar()
	plt.show()