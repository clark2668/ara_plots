# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt

# very helpful: http://icecube.berkeley.edu/calibration/IceCubeCoordinateSystem.pdf
	
def main():

	data = np.genfromtxt("icecube_dom60.csv", delimiter=',',skip_header=1, names=['x','y','z','string','dom'])

	xs = data['x']
	ys = data['y']
	zs = data['z']
	strings = data['string']
	doms = data['dom']

	icecube_center=np.array([46500,52200])
	icecube_center=np.multiply(icecube_center,0.3048/1000.)
	# icecube_center*=.3048/1000.

	xs/=1000.
	ys/=1000.
	xs+=icecube_center[0]
	ys+=icecube_center[1]

	A1=np.array([38800.69, 51066.22])
	A2=np.array([35528.80, 45382.35])
	A3=np.array([32246.30, 51065.59])
	A4=np.array([35382.53, 56916.55])
	A5=np.array([32252.22, 39718.52])
	A1 = np.multiply(A1, 0.3048/1000.);
	A2 = np.multiply(A2, 0.3048/1000.);
	A3 = np.multiply(A3, 0.3048/1000.);
	A4 = np.multiply(A4, 0.3048/1000.);
	A5 = np.multiply(A5, 0.3048/1000.);

	# A1*=np.multiply(.3048/1000.)
	# A2*=.3048/1000.
	# A3*=.3048/1000.
	# A4*=.3048/1000.
	# A5*=.3048/1000.

	fig = plt.figure(figsize=(10*1.5,10*1.5))
	ax=fig.add_subplot(111)

	ax.plot(xs,ys,'o',color='gray',markersize=6,linewidth=0)
	ax.plot(A1[0],A1[1], "bo", label="A1", markersize=24)
	ax.plot(A2[0],A2[1], "go", label="A2", markersize=24)
	ax.plot(A3[0],A3[1], "ro", label="A3", markersize=24)
	ax.plot(A4[0],A4[1], "co", label="A4", markersize=24)
	ax.plot(A5[0],A5[1], "mo", label="A5", markersize=24)
	ax.plot(A5[0],A5[1], "mo", label="A5", markersize=35, fillstyle='none')
	# ax.patches.Circle(A5[0],A5[1],radius=2)
	# ax.scatter(A5[0],A5[1], s=40, facecolor='none', edgecolor='m')

	ax.text(13.7 ,16.6, 'IceCube', horizontalalignment='left', verticalalignment='center', size=24, color='gray', fontweight='bold')
	ax.text(12 ,15.7, 'A1', horizontalalignment='left', verticalalignment='center', size=24, color='b', fontweight='bold')
	ax.text(11 ,14.2, 'A2', horizontalalignment='left', verticalalignment='center', size=24, color='g', fontweight='bold')
	ax.text(10 ,15.7, 'A3', horizontalalignment='left', verticalalignment='center', size=24, color='r', fontweight='bold')
	ax.text(11 ,17.6, 'A4', horizontalalignment='left', verticalalignment='center', size=24, color='c', fontweight='bold')
	ax.text(10 ,12.3, 'A5+PA', horizontalalignment='left', verticalalignment='center', size=24, color='m', fontweight='bold')

	# ax.legend()

	ax.set_xlabel('Easting (km)',size=30)
	ax.set_ylabel('Northing (km)',size=30)
	ax.tick_params(labelsize=30,pad=10)

	ax.set_xlim(9,15)
	ax.set_ylim(11.5,18)

	fig.savefig("ara_map.png",edgecolor='none',bbox_inches="tight",dpi=300)

#actually execute the main function
main()


			
		
