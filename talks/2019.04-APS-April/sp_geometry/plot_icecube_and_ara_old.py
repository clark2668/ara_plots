# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

def drawSphere(xCenter, yCenter, zCenter, r):
    #draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x=np.cos(u)*np.sin(v)
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    # shift and scale sphere
    x = r*x + xCenter
    y = r*y + yCenter
    z = r*z + zCenter
    return (x,y,z)
	
def main():

	fig = plt.figure(figsize=(11.5*1.5,8.5*1.5))
	ax = fig.add_subplot(111,projection='3d')

	data = np.genfromtxt("icecube_dom_locations.csv", delimiter=',',skip_header=1, names=['x','y','z','string','dom'])

	xs = data['x']
	ys = data['y']
	zs = data['z']
	strings = data['string']
	doms = data['dom']


	Xs=[0]
	Ys=[0]
	Zs=[0]

	iterator = xs.shape[0]
	i=0
	while i<iterator:
		if int(strings[i])<82:
			xs[i]+=3345
			ys[i]+=2079
			if i==0:
				Xs[0]=xs[0]
				Ys[0]=ys[0]
				Zs[0]=zs[0]
			else:
				Xs.append(xs[i])
				Ys.append(ys[i])
				Zs.append(zs[i])
		i+=1	

	#fig2 = plt.figure(figsize=(11.5*1.5,8.5*1.5))
	#ax2=fig2.add_subplot(111)
	#ax2.scatter(Xs,Ys)
	#fig2.savefig("x-y-projection.pdf",edgecolor='none',bbox_inches="tight")
	
	linex = np.linspace(2926,2926+(2000*.69),100)
	liney = np.linspace(1975,1975+(2000*.45),100)
	linez = np.linspace(-1432-20,-1432-20+(2000*-0.55),100)



	a2x=[-15,-15,15,15,-15,-15,15,15]
	a2y=[-15,15,-15,15,-15,-15,15,15]
	a2z=[-159,-159,-159,-159,-189,-189,-189,-189]	
	ax.scatter(a2x,a2y,a2z,c='r',marker='s')	
	ax.scatter(Xs,Ys,Zs,c='b',marker='o',depthshade=0)
	ax.scatter(linex,liney,linez,c='k',marker='o')
	
	(sx,sy,sz) = drawSphere(2926,1975,-1432-20,200)
	ax.plot_wireframe(sx,sy,sz,color="r")
	ax.view_init(0, 0)
	
	ax.set_ylim(-500,3000) #set the y limits of the plot
	ax.set_xlim(-500,4500) #set the y limits of the plot
	#ax.scatter(data['x'],data['y'],data['z'],c='b',marker='o')

	ax.set_xlabel('East')
	ax.set_ylabel('North')
	ax.set_zlabel('Elevation')
	fig.savefig("icecube_dom_locations.pdf",edgecolor='none',bbox_inches="tight") #save the figure

#actually execute the main function
main()


			
		
