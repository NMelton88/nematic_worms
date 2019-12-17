import numpy
import math

x, y = [], []
xnew, ynew = [], []

W = 1							# Number of Worms
N = 100		
r0 = 0.9
s = 1.0					
r_cut = 2**(1./6.)*s
dt = 0.001
g = 1.0
kT = 1.0
ks = 60.0
k = 1.0
e = 0.1
a = 20
box = 200

for i in range(0,N):
	x.append(15 + 0.9 * i)
	y.append(50.0)
	xnew.append(0.0)
	ynew.append(0.0)

myfile = open("myworm.xyz","w")
myfile.write("{}\n".format(N))
for i in range(0,N):
	myfile.write("H\t {}\t {}\t {}\t\n".format(x[i],y[i],0.0))

def LJ(r_test):
	f_lj = 24.0 * e * (-2.0*(s**12)/(r_test**13) + 1.0*(s**6)/(r_test**7))
	return f_lj

def Spring(r_test):
	f_sp = ks * (r_test - r0)
	return f_sp

def Bending(r_test):
	f_b = k/2.0 * r_test
	return f_b


for t in range(0,1000):
	for i in range(0,N):
		xold = x[i]
		yold = y[i]
		springx, springy = 0.0, 0.0
		ljx, ljy = 0.0, 0.0
		bendx, bendy = 0.0, 0.0
		for j in range(0,N):
			dx = x[j] - xold
			dy = y[j] - yold
			if math.fabs(dx) > box/2.0:  dx = dx - numpy.sign(dx) * box
			if math.fabs(dy) > box/2.0:  dy = dy - numpy.sign(dy) * box
			r1 = math.sqrt((dx)**2 + (dy)**2)
			if j == i-1:
				springx = springx + Spring(r1)*(dx)/r1
				springy = springy + Spring(r1)*(dy)/r1
			if j == i+1:
				springx = springx + Spring(r1)*(dx)/r1
				springy = springy + Spring(r1)*(dy)/r1

			if r1 < r_cut:
				if j != i:
					ljx = ljx + LJ(r1)*(dx)/r1
					ljy = ljy + LJ(r1)*(dy)/r1

		if (i < N-2):
			for j in range(i,i+2):
				dx = -x[j+1] - xold
				dy = -y[j+1] - yold
				if math.fabs(dx) > box/2.0:  dx = dx - numpy.sign(dx) * box
				if math.fabs(dy) > box/2.0:  dy = dy - numpy.sign(dy) * box
				r2 = math.sqrt((dx)**2 + (dy)**2)
				bendx = bendx + Bending(r2)*dx/r2
				bendy = bendy + Bending(r2)*dy/r2
		


		if (i < N-2):
			dx = x[i+1]-xold
			dy = y[i+1]-yold
			if math.fabs(dx) > box/2.0:  dx = dx - numpy.sign(dx) * box
			if math.fabs(dy) > box/2.0:  dy = dy - numpy.sign(dy) * box
			xnew[i] = xold + springx*dt/g + math.sqrt(2.0*kT*dt/g)*numpy.random.normal(0.0,0.1) + ljx*dt/g + a * (dx)/math.sqrt((dx)**2 + (dy)**2)*dt/g + bendx*dt/g
			ynew[i] = yold + springy*dt/g + math.sqrt(2.0*kT*dt/g)*numpy.random.normal(0.0,0.1) + ljy*dt/g + a * (dy)/math.sqrt((dx)**2 + (dy)**2)*dt/g + bendy*dt/g

		if (i < N-1):
			dx = x[i+1]-xold
			dy = y[i+1]-yold
			if math.fabs(dx) > box/2.0:  dx = dx - numpy.sign(dx) * box
			if math.fabs(dy) > box/2.0:  dy = dy - numpy.sign(dy) * box
			xnew[i] = xold + springx*dt/g + math.sqrt(2.0*kT*dt/g)*numpy.random.normal(0.0,0.1) + ljx*dt/g + a * (dx)/math.sqrt((dx)**2 + (dy)**2)*dt/g
			ynew[i] = yold + springy*dt/g + math.sqrt(2.0*kT*dt/g)*numpy.random.normal(0.0,0.1) + ljy*dt/g + a * (dy)/math.sqrt((dx)**2 + (dy)**2)*dt/g

		else:
			dx = x[i-1]-xold
			dy = y[i-1]-yold
			if math.fabs(dx) > box/2.0:  dx = dx - numpy.sign(dx) * box
			if math.fabs(dy) > box/2.0:  dy = dy - numpy.sign(dy) * box
			xnew[i] = xold + springx*dt/g + math.sqrt(2.0*kT*dt/g)*numpy.random.normal(0.0,0.1) + ljx*dt/g - a * (dx)/math.sqrt((dx)**2 + (dy)**2)*dt/g
			ynew[i] = yold + springy*dt/g + math.sqrt(2.0*kT*dt/g)*numpy.random.normal(0.0,0.1) + ljy*dt/g - a * (dy)/math.sqrt((dx)**2 + (dy)**2)*dt/g

		if xnew[i] < 0:  xnew[i] = xnew[i] + box
		if xnew[i] > box:  xnew[i] = xnew[i] - box
		if ynew[i] < 0:  ynew[i] = ynew[i] + box
		if ynew[i] > box:  ynew[i] = ynew[i] - box

	for i in range(0,N):
		x[i] = xnew[i]
		y[i] = ynew[i]
		#myfile.write("H\t {}\t {}\t {}\t\n".format(x[i],y[i],0.0))

	remainder = t%50
	if remainder == 0:
		for k in range(0,N):
			myfile.write("H\t {}\t {}\t {}\t\n".format(x[k],y[k],0.0))

myfile.close()
