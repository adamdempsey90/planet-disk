
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
		self.Ly = 2*self.Ny*self.y[0]/(1-self.Ny)
	
		
def realspace(fld):
#	x=hstack((-fld.x[::-1],fld.x))
#	y=fld.y
	y,x = meshgrid(fld.y,fld.x)
#	sig=vstack((conj(fld.sig[::-1,:]),fld.sig))
#	u = vstack((-conj(fld.u[::-1,:]),fld.u))
#	v = vstack((-conj(fld.v[::-1,:]),fld.v))
	vx = fft.irfft(fld.u)*fld.Ny
	vy = fft.irfft(fld.v)*fld.Ny
	dens = fft.irfft(fld.sig)*fld.Ny
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
	
def plotfld(fld,q,ilist,xlims=(0,0)):
	
	
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

	for i in ilist:
		tstr = 'mode #'+str(i)+', k = ' + str(fld.k[i])
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

def animate(q,k,t,Nx,Ny,dt=1,dir='',scale=True):
	
	dat= zeros((Nx,len(t)),dtype='complex')
	x = Field(t[0],Nx,Ny,dir=dir).x
	for i in t:
		if q=='u':
			dat[:,i] = Field(i,Nx,Ny,dir=dir).u[:,k]
		elif q=='v':
			dat[:,i] = Field(i,Nx,Ny,dir=dir).v[:,k]
		
		elif q=='sig':
			dat[:,i] = Field(i,Nx,Ny,dir=dir).sig[:,k]
		else:
			print 'Not a valid variable name'
			return
	
	fig=figure(figsize=(15,10))
	if scale:
		if k==0:
			qmin = real(dat.min())
			qmax = real(dat.max())
			
			print qmin,qmax
		else:
			qmin = dat.min()
			qmax = dat.max()
			if real(qmin) < imag(qmin):
				qmin = real(qmin)
			else:
				qmin = imag(qmin)
			
			if real(qmax) > imag(qmax):
				qmax = real(qmax)
			else:
				qmax = imag(qmax)	
			
	for i in t:
		clf()
		if k==0:
			plot(x,real(dat[:,i]))
			
		else:
			plot(x,real(dat[:,i]),x,imag(dat[:,i]))

		if scale:
			ylim((qmin,qmax))
		title('t='+str(i*dt)+'$\Omega$')
		
		fig.canvas.draw()
	return
	
	
def animate_real(q,t,Nx,Ny,dt=1,dir='',logscale=False):
	
	dat= zeros((Nx,Ny,len(t)))
	fld = Field(t[0],Nx,Ny,dir=dir)
	x,y = meshgrid(fld.y,fld.x)
	for i in t:
		fld = Field(i,Nx,Ny,dir=dir)
		vx,vy,dens,_,_ = realspace(fld)
		if q=='vx':
			dat[:,:,i] = vx
			tstr = 'vx at '
		elif q=='vy':
			dat[:,:,i] = vy
			tstr = 'vy at '
		elif q=='dens':
			dat[:,:,i] = dens
			tstr = '$\Sigma$ at '
		else:
			print 'Not a valid variable name'
			return
	
	
	d_min = dat.min(); d_max = dat.max();

	fig=figure(figsize=(15,10))
	for i in t:
		clf()
		if logscale:
			imshow(log10(dat[:,:,i]).transpose(),aspect='auto',vmin=log10(d_min),vmax=log10(d_max))
		else:
			imshow(dat[:,:,i].transpose(),aspect='auto',vmin=d_min,vmax=d_max)
		colorbar()
		title(tstr + 't='+str(i*dt)+'$\Omega$')
		fig.canvas.draw()
		
		
	return
def plotreal(fld,tstr=''):
	vx,vy,dens,x,y = realspace(fld)
	figure()
	pcolor(x,y,vx); colorbar()
	xlabel('x')
	title('v_x       ' + tstr)

	figure()
	pcolor(x,y,vy); colorbar()
	xlabel('x')
	title('\deltav_y     ' + tstr)
	
	figure()
	pcolor(x,y,log10(dens)); colorbar()
	xlabel('x')
	title('log10($\Sigma$)      ' + tstr)
	
	return

def calcTwd(fld,nu):
	
	dx = fld.dx
	dy = fld.Ly / fld.Ny
	
	vx,vy,dens,x,y = realspace(fld)
	
	dxu,dyu = gradient(vx,dx,dy)
	dxv,dyv = gradient(vy,dx,dy)
	dxd,dyd = gradient(dens,dx,dy)
	
	pixy = nu*(dxv+dyu)
	Ftot = dens*vx*vy
	Fb = Ftot.mean(axis=1)
	Fp = [Ftot[i,:]-f for i,f in enumerate(Fb)]
	
	
	
	return