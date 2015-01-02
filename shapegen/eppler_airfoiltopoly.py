from sys import argv
import numpy as np
import matplotlib.pyplot as plt

def readfile(fname):

	infile=open(fname,'r');

	infile.readline();
	line=infile.readline();

	nsucside  = int(float(line.split()[0]));
	npresside = int(float(line.split()[1]));

	print "nsucside:",nsucside
	print "npresside:",npresside

	infile.readline()

	suction_side_x = np.array([])
	suction_side_y = np.array([])
	pres_side_x = np.array([])
	pres_side_y = np.array([])

	for i in range(nsucside):

		line = infile.readline()
		suction_side_x=np.append(suction_side_x,float(line.split()[0]))
		suction_side_y=np.append(suction_side_y,float(line.split()[1]))

	infile.readline()
	
	for i in range(npresside):

		line = infile.readline()
		pres_side_x = np.append(pres_side_x,float(line.split()[0]))
		pres_side_y = np.append(pres_side_y,float(line.split()[1]))

	infile.close()

	return(suction_side_x,suction_side_y,pres_side_x,pres_side_y)




sx,sy,p1x,p1y=readfile(argv[1])

px = p1x[::-1];
py = p1y[::-1];

outfile=open("polypoints.inp",'w')

fact=1;
nshapes=1

outfile.write("%d\n"%(nshapes))
outfile.write("%d\n"%((sx.size+px.size-2)/fact))

for i in range(sx.size/fact):
	outfile.write("%f\t%f\n"%(sx[fact*i],sy[fact*i]))

for i in range((px.size-2)/fact):
	outfile.write("%f\t%f\n"%(px[fact*(i+1)],py[fact*(i+1)]))

outfile.close()

plt.plot(sx,sy,'r*-')
plt.plot(px,py,'bo-')

plt.show()

