
class Field():
	def __init__(self,time,Nx,Ny):
		NC = Ny/2+1
		self.Nx = Nx
		self.Ny = Ny
		self.NC = NC
		self.NR = 2*NC
		self.sig=fromfile('dens_'+str(time)+'.dat',dtype='complex').reshape((Nx,NC))
		self.u=fromfile('vx_'+str(time)+'.dat',dtype='complex').reshape((Nx,NC))
		self.v=fromfile('vy_'+str(time)+'.dat',dtype='complex').reshape((Nx,NC))
		coords=fromfile('coords.dat');
		self.x=coords[:Nx]
		self.k=coords[Nx:Nx+NC]
		self.y=coords[Nx+NC:]
		

def derivs(num,Nx,Ny):
	NC = Ny/2+1
	dxvx = fromfile('dxvx_'+str(num)+'.dat',dtype='complex').reshape((Nx,NC))
	dxvy = fromfile('dxvy_'+str(num)+'.dat',dtype='complex').reshape((Nx,NC))
	dxsig = fromfile('dxdens_'+str(num)+'.dat',dtype='complex').reshape((Nx,NC))
	return dxvx,dxvy,dxsig