
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
		self.dx = diff(self.x)[0]
		self.Lx = self.Nx*self.dx
		
def realspace(fld):
	x=hstack((-fld.x[::-1],fld.x))
	y=fld.y
	sig=vstack((conj(fld.sig[::-1,:]),fld.sig))
	u = vstack((-conj(fld.u[::-1,:]),fld.u))
	v = vstack((-conj(fld.v[::-1,:]),fld.v))
	vx = fft.irfft(u)*fld.Ny
	vy = fft.irfft(v)*fld.Ny
	dens = fft.irfft(sig)*fld.Ny
	return vx,vy,dens,x,y
	
def derivs(num,Nx,Ny,NG=0,dir=''):
	NC = Ny/2+1
	dxvx = fromfile(dir+'dxvx_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	dxvy = fromfile(dir+'dxvy_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	dxsig = fromfile(dir+'dxdens_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	return dxvx,dxvy,dxsig
	
def rhs(num,Nx,Ny,dir=''):
	NC = Ny/2+1
	dtu = fromfile(dir+'dtvx_'+str(num)+'.dat',dtype='complex').reshape((Nx,NC))
	dtv = fromfile(dir+'dtvy_'+str(num)+'.dat',dtype='complex').reshape((Nx,NC))
	dtsig = fromfile(dir+'dtdens_'+str(num)+'.dat',dtype='complex').reshape((Nx,NC))
	return dtu,dtv,dtsig
	
def tens(num,Nx,Ny,NG=0,dir=''):
	NC = Ny/2+1
	pixx = fromfile(dir+'pixx_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	pixy = fromfile(dir+'pixy_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	piyy = fromfile(dir+'piyy_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	dxpi = fromfile(dir+'divpix_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	dypi = fromfile(dir+'divpiy_'+str(num)+'.dat',dtype='complex').reshape((Nx+2*NG,NC))
	return pixx,pixy,piyy,dxpi,dypi
	
def plotfld(fld,q,i,xlims=(0,0)):
	
	tstr = 'mode #'+str(i)+', k = ' + str(fld.k[i])
	x = fld.x
	if q=='u':
		dat = fld.u
		ystr = 'u'
	elif q=='v':
		dat = fld.v
		ystr = 'v'
	elif q=='sig':
		dat = fld.sig
		ystr = '$\Sigma$'
	else:
		print "Not valid variable name"
		return

	
	figure()
	if xlims[0]!=0 and xlims[1]!=0:
		xlim(xlims)
	if i==0:
		plot(x,real(dat[:,i]))
		ylabel('<'+ystr+'>')
	else:
		plot(x,real(dat[:,i]),x,imag(dat[:,i]))
		ylabel(ystr)
		legend(('Re('+ystr+')','Im('+ystr+')'))

	xlabel('x')
	title(tstr)
	return
