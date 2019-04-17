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
	icecube_center*=.3048/1000.

	xs/=1000.
	ys/=1000.
	xs+=icecube_center[0]
	ys+=icecube_center[1]

	IC22=np.array([-492.43,-230.16])
	IC22/=1000.
	IC22[0]+=icecube_center[0]
	IC22[1]+=icecube_center[1]

	IC1=np.array([-256.14,-521.08])
	IC1/=1000.
	IC1[0]+=icecube_center[0]
	IC1[1]+=icecube_center[1]

	A1=np.array([38800.69, 51066.22])
	A2=np.array([35528.80, 45382.35])
	A3=np.array([32246.30, 51065.59])
	A4=np.array([35382.53, 56916.55])
	A5=np.array([32252.22, 39718.52])
	A1*=.3048/1000.
	A2*=.3048/1000.
	A3*=.3048/1000.
	A4*=.3048/1000.
	A5*=.3048/1000.

	fig = plt.figure(figsize=(10*1.5,10*1.5))
	ax=fig.add_subplot(111)

	ax.plot(xs,ys,'o',color='gray',markersize=6,linewidth=0)
	ax.plot(A1[0],A1[1], "bo", label="A1", markersize=24)
	ax.plot(A2[0],A2[1], "go", label="A2", markersize=24)
	ax.plot(A3[0],A3[1], "ro", label="A3", markersize=24)
	ax.plot(A4[0],A4[1], "co", label="A4", markersize=24)
	ax.plot(A5[0],A5[1], "mo", label="A5", markersize=24)
	ax.plot(IC22[0],IC22[1],"ko", label="IC22",markersize=18)
	ax.plot(IC1[0],IC1[1],"ko", label="IC1",markersize=18)

	print "IC22 to A4: ",np.sqrt((A4[0]-IC22[0])**2 + (A4[1]-IC22[1])**2)
	print "IC22 to A5: ",np.sqrt((A5[0]-IC22[0])**2 + (A5[1]-IC22[1])**2)

	# ax.legend()

	ax.set_xlabel('Easting (km)',size=30)
	ax.set_ylabel('Northing (km)',size=30)
	ax.tick_params(labelsize=30,pad=10)

	ax.set_xlim(9,15)
	ax.set_ylim(11.5,17.5)

	fig.savefig("x-y-projection.pdf",edgecolor='none',bbox_inches="tight")

#actually execute the main function
main()


			
		
