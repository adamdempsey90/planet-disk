
class Field():
	def __init__(self,time,Nx,Ny,NG=0,dir=''):
		NC = Ny/2+1
		self.Nx = Nx
		self.Ny = Ny
		self.NC = NC
		self.NR = 2*NC
		self.NG = NG
		self.sig=fromfile(dir+'dens_'+str(time)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
		self.u=fromfile(dir+'vx_'+str(time)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
		self.v=fromfile(dir+'vy_'+str(time)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
		coords=fromfile(dir+'coords.dat');
		self.x=coords[:Nx+2*NG]
		self.k=coords[Nx+2*NG:Nx+2*NG+NC]
		self.y=coords[Nx+2*NG+NC:]
		

def derivs(num,Nx,Ny,NG=0):
	NC = Ny/2+1
	dxvx = fromfile('dxvx_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	dxvy = fromfile('dxvy_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	dxsig = fromfile('dxdens_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	return dxvx,dxvy,dxsig
	
def rhs(num,Nx,Ny):
	NC = Ny/2+1
	dtu = fromfile('dtvx_'+str(num)+'.dat',dtype='complex').reshape((Nx,NC))
	dtv = fromfile('dtvy_'+str(num)+'.dat',dtype='complex').reshape((Nx,NC))
	dtsig = fromfile('dtdens_'+str(num)+'.dat',dtype='complex').reshape((Nx,NC))
	return dtu,dtv,dtsig
	
def tens(num,Nx,Ny,NG=0):
	NC = Ny/2+1
	pixx = fromfile('pixx_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	pixy = fromfile('pixy_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	piyy = fromfile('piyy_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	dxpi = fromfile('divpix_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	dypi = fromfile('divpiy_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	return pixx,pixy,piyy,dxpi,dypi