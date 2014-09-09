
class Field():
	def __init__(self,time,np,dir=''):
		if dir[-1]!='/':
			dir += '/'
		vx,vy,dens,x,self.k,self.y,Nx,NC,splits = loadvars(np,time,dir);
		self.Nx = Nx.sum()
		self.NC = NC[0]
		self.Ny = 2*self.NC-2
		self.NR= 2*self.NC
		
		self.u = vx[0]
		self.v = vy[0]
		self.sig = dens[0]
		self.x = x[0]

		for i in range(1,np):
			self.u = vstack((self.u,vx[i]))
			self.v = vstack((self.v,vy[i]))
			self.sig = vstack((self.sig,dens[i]))
			self.x = hstack((self.x,x[i]))

		self.dx = diff(self.x)[0]
		self.Lx = self.Nx*self.dx
		self.Ly = 2*self.Ny*self.y[0]/(1-self.Ny)
		self.splits = self.x[splits]
			
 		dxv,_=gradient(self.v,self.dx,1)
 		self.kk=ones(self.u.shape)*2*pi/self.Ly
 		for j in range(self.NC):
 			self.kk[:,j] *= j
 		self.vort = dxv - 1j*self.kk*self.u	
	
	def plot(self,q,ilist,xlims=(0,0),scale=1,conj_flag=False):
		x = self.x
		if q=='u':
			dat = self.u
			ystr = 'u'
		elif q=='v':	
			dat = self.v
			ystr = 'v'
		elif q=='sig':
			dat = self.sig
			ystr = '$\Sigma$'
		else:
			print "Not valid variable name"
			return
		if conj_flag:
			dat = conj(dat)
		for i in ilist:
			tstr = 'mode #'+str(i)+', k = ' + str(self.k[i])
			figure()
			if xlims!=(0,0):
				xlim(xlims)
			if i==0:
				plot(x,scale*real(dat[:,i]))
				ylabel('<'+ystr+'>')
			else:
				plot(x,scale*real(dat[:,i]),x,scale*imag(dat[:,i]))
				ylabel(ystr)
				legend(('Re('+ystr+')','Im('+ystr+')'))

			xlabel('x')	
			title(tstr)
		return
	def plotreal(self,xlims=(0,0),ylims=(0,0)):
		vx,vy,dens,x,y = realspace(self)
		vortens = (fft.irfft(self.vort)+2)/dens
		
		exarg = (-self.Lx/2,self.Lx/2,-self.Ly/2,self.Ly/2)
		
		
		(endx,endy) = vx.shape
		(startx,starty) = (0,0)
		
		if xlims!=(0,0):
			temp = (self.x > xlims[0]) & (self.x < xlims[1])
			startx = next(i for i,j in enumerate(temp) if j)
			endx = next(len(temp)-i for i,j in enumerate(temp[::-1]) if j)
			exarg = xlims + exarg[-2:]
		
		if ylims!=(0,0):
			temp = (self.y > ylims[0]) & (self.y < ylims[1])
			starty = next(i for i,j in enumerate(temp) if j)
			endy = next(len(temp)-i for i,j in enumerate(temp[::-1]) if j)
			exarg = exarg[:2] + ylims


		print (startx,endx,starty,endy)
		
	
		figure()
		imshow(vx[startx:endx,starty:endy].transpose(), aspect='auto', origin='lower',extent=exarg,interpolation='bilinear')
		colorbar()
		title('vx'); xlabel('x'); ylabel('y')
	
		figure()
		imshow(vy[startx:endx,starty:endy].transpose(), aspect='auto', origin='lower', extent=exarg,interpolation='bilinear')
		colorbar()
		title('vy'); xlabel('x'); ylabel('y')
	
		figure()
		imshow(dens[startx:endx,starty:endy].transpose(), aspect='auto', origin='lower', extent=exarg,interpolation='bilinear')
		colorbar()
		title('$\Sigma$'); xlabel('x'); ylabel('y')
		
		figure()
		imshow(vortens[startx:endx,starty:endy].transpose(), aspect='auto', origin='lower', extent=exarg,interpolation='bilinear')
		colorbar()
		title('Vortensity'); xlabel('x'); ylabel('y')
		return
	
def loadvars(np,time,dir):
	vx = []
	vy = []
	dens = []
	x = []
	y = []
	k = []
	
	Nx = arange(np)
	NC = arange(np)

	for i in range(np):
		temp = fromfile(dir+'id'+str(i)+'/vx_'+str(time)+'.dat',dtype='complex')
				

		Nx[i] = int(real(temp[0]))
		NC[i] = int(real(temp[1]))/2+1
		vx.append(temp[2:].reshape((Nx[i],NC[i])))
		
		temp = fromfile(dir+'id'+str(i)+'/vy_'+str(time)+'.dat',dtype='complex')
		vy.append(temp[2:].reshape((Nx[i],NC[i])))
		temp = fromfile(dir+'id'+str(i)+'/dens_'+str(time)+'.dat',dtype='complex')
		dens.append(temp[2:].reshape((Nx[i],NC[i])))
		
		temp = fromfile(dir+'id'+str(i)+'/coords.dat')
		temp=temp[2:]
		x.append(temp[:Nx[i]])
		if i==0:
			k=temp[Nx[i]:Nx[i]+NC[i]]
			y=temp[Nx[i]+NC[i]:]
	
	splits=[]
	for i in range(1,np):
		temp = 0
		for j in range(i):
			temp += Nx[j]
		splits.append(temp)
	
	return vx,vy,dens,x,k,y,Nx,NC,splits
			
def realspace(fld):
#	x=hstack((-fld.x[::-1],fld.x))
#	y=fld.y
	y,x = meshgrid(fld.y,fld.x)
#	sig=vstack((conj(fld.sig[::-1,:]),fld.sig))
#	u = vstack((-conj(fld.u[::-1,:]),fld.u))
#	v = vstack((-conj(fld.v[::-1,:]),fld.v))
	vx = fft.irfft(fld.u)
	vy = fft.irfft(fld.v)
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

def plotvortens(fld,tstr=''):
	
	rvort=fft.irfft(fld.vort)
	vx,vy,dens,x,y=realspace(fld)
	
	vortens = (rvort+2)/dens
	vbar = vortens.mean(axis=1)
	vortensbar = ones(dens.shape)
	for j in range(fld.Ny):
		vortensbar[:,j] *= vbar
	
	figure()
	plot(fld.x,vbar)
	xlabel('x'); ylabel('$<\\xi>$')
	
	figure()
	pcolor(x,y,vortens); colorbar()
	xlabel('x'); ylabel('y')
	title('$\\xi$      ' + tstr)

	figure()
	pcolor(x,y,(vortens-vortensbar)/vortensbar,cmap='RdBu'); colorbar()
	xlabel('x'); ylabel('y')
	title('$\Delta \\xi / \\xi$       ' + tstr)

	
	return vortens,vortensbar
	
def plotfld(fld,q,ilist,xlims=(0,0),scale=1,conj_flag=False):
	print xlims[0],xlims[1]
	
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
	if conj_flag:
		dat = conj(dat)
	for i in ilist:
		tstr = 'mode #'+str(i)+', k = ' + str(fld.k[i])
		figure()
		if ~(xlims[0]==0 and xlims[1]==0):
			xlim(xlims)
		if i==0:
			plot(x,scale*real(dat[:,i]))
			ylabel('<'+ystr+'>')
		else:
			plot(x,scale*real(dat[:,i]),x,scale*imag(dat[:,i]))
			ylabel(ystr)
			legend(('Re('+ystr+')','Im('+ystr+')'))

		xlabel('x')
		title(tstr)
	return

def animate(q,k,t,np,dt=1,xlims=(0,0),dir='',scale=True,norm=1,conj_flag=False):
	
	
	temp = Field(t[0],np,dir=dir)
	x = temp.x
	Nx = temp.Nx; Ny = temp.Ny; NC = temp.NC;
	dat= zeros((Nx,len(t)),dtype='complex')
	for i in t:
		if q=='u':
			dat[:,i] = norm*Field(i,np,dir=dir).u[:,k]
		elif q=='v':
			dat[:,i] = norm*Field(i,np,dir=dir).v[:,k]
		
		elif q=='sig':
			dat[:,i] = norm*Field(i,np,dir=dir).sig[:,k]
		else:
			print 'Not a valid variable name'
			return
	
	if conj_flag:
		dat = conj(dat)
		
	fig=figure(figsize=(15,10))
	if scale:
		if k==0:
			qmin = real(dat.min())
			qmax = real(dat.max())
			
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
		if ~(xlims[0]==0 and xlims[1]==0):
			xlim(xlims)
		title('t$\Omega$='+str(i*dt)+'\t k='+str(k*2*pi/temp.Ly))
		
		fig.canvas.draw()
	return
	
	
def animate_real(q,t,np,dt=1,dir='',logscale=False):
	temp = Field(t[0],np,dir=dir)
	Nx=temp.Nx; Ny = temp.Ny; 
	NC=temp.NC
	x,y = meshgrid(temp.x,temp.y)
	dat= zeros((Nx,Ny,len(t)))
	for i in t:
		fld = Field(i,np,dir=dir)
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
			imshow(log10(dat[:,:,i]).transpose(),aspect='auto',vmin=log10(d_min),vmax=log10(d_max), \
			extent=(temp.x[0],temp.x[-1],temp.y[0],temp.y[-1]))
		else:
			imshow(dat[:,:,i].transpose(),aspect='auto',vmin=d_min,vmax=d_max, \
			extent=(temp.x[0],temp.x[-1],temp.y[0],temp.y[-1]))
		colorbar()
		title(tstr + 't$\Omega$='+str(i*dt))
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
	title('\delta v_y     ' + tstr)
	
	figure()
	pcolor(x,y,log10(dens)); colorbar()
	xlabel('x')
	title('log10($\Sigma$)      ' + tstr)
	
	return

def calcTwd(fld,nu):
	
# 	dx = fld.dx
# 	dy = fld.Ly / fld.Ny
	
# 	vx,vy,dens,x,y = realspace(fld)
# 	
# 	dxu,dyu = gradient(vx,dx,dy)
# 	dxv,dyv = gradient(vy,dx,dy)
# 	dxd,dyd = gradient(dens,dx,dy)
# 	
# 	pixy = nu*(dxv+dyu)
# 	Ftot = dens*vx*vy
# 	Fb = Ftot.mean(axis=1)
# 	Fp = [Ftot[i,:]-f for i,f in enumerate(Fb)]
	
	dxu,_ = gradient(fld.u,fld.dx,1)
	dxv,_ = gradient(fld.v,fld.dx,1)
	dxsig,_ = gradient(fld.sig,fld.dx,1)
	
	k=ones(fld.u.shape)*2*pi/fld.Ly
	for j in range(fld.NC):
		k[:,j] *= j
	
	Twd = 2*real( -fld.sig[:,0]*(conj(fld.u[:,1:])*dxv[:,1:]).sum(axis=1) \
					+(dxv[:,0]+2*fld.x)*(conj(fld.sig[:,1:])*fld.u[:,1:]).sum(axis=1))
					
	divv = dxu + 1j*k*fld.v				
	Txy = dxv+1j*k*fld.u
	Tyy = 2*1j*fld.v - (2./3)*divv
	
	dens = fft.irfft(fld.sig,fld.Ny)
	rTxy = fft.irfft(Txy,fld.Ny)
	rTyy = fft.irfft(Tyy,fld.Ny)
	Pixy = nu*fft.rfft(dens*rTxy)
	Piyy = -fld.sig + nu*fft.rfft(dens*rTyy)
	
	divPi,_ = gradient(Pixy,fld.dx,1)
	divPi += 1j*k*Piyy
	
	rdivPi = fft.irfft(divPi,fld.Ny)

	visc = fft.rfft(rdivPi/dens)
	
	Twd += 2*real((conj(fld.sig[:,1:])*visc[:,1:]).sum(axis=1) + fld.sig[:,0]*visc[:,0])
	
	return Twd
	
def loadphi(np,dir=''):
	x=[]
	phi=[]
	for i in range(np):
		temp = fromfile(dir+'id'+str(i)+'/dtvy_0.dat',dtype='complex')
		Nx = int(real(temp[0]))
		NC = int(real(temp[1]))/2 + 1
		print Nx, NC
		phi.append(temp[2:].reshape((Nx,NC)))
		
		temp = fromfile(dir+'id'+str(i)+'/coords.dat')
		
		x.append(temp[2:Nx+2])
	
	return phi,x
	
	
def writerun(fld,fname=''):

	with open(fname+'_rdens.dat',"w") as rd:
		with open(fname+'_rvx.dat',"w") as rg:
			with open(fname+'_rvy.dat',"w") as rh:
				with open(fname+'_idens.dat',"w") as id:
					with open(fname+'_ivx.dat',"w") as ig:
						with open(fname+'_ivy.dat',"w") as ih:	
							with open(fname+'_x.dat','w') as fx:		
								for i in range(fld.Nx):
									fx.write(str(fld.x[i])+'\n')
									for j in range(fld.NC):
										rd.write(str(real(fld.sig[i,j]))+'\t')
										rg.write(str(real(fld.u[i,j]))+'\t')
										rh.write(str(real(fld.v[i,j]))+'\t')
										id.write(str(imag(fld.sig[i,j]))+'\t')
										ig.write(str(imag(fld.u[i,j]))+'\t')
										ih.write(str(imag(fld.v[i,j]))+'\t')									
									rd.write('\n')
									rg.write('\n')
									rh.write('\n')
									id.write('\n')
									ig.write('\n')
									ih.write('\n')	

def pspec(fld):
	
	u2 = zeros((fld.u.shape[1]))
	v2 = zeros((fld.u.shape[1]))
	sig2 = zeros((fld.u.shape[1]))
	uv = zeros((fld.u.shape[1]))
	usig = zeros((fld.u.shape[1]))

	
	u2 = (2*real(conj(fld.u)*fld.u)).mean(axis=0)
	v2 = (2*real(conj(fld.v)*fld.v)).mean(axis=0)
	sig2 = (2*real(conj(fld.sig)*fld.sig)).mean(axis=0)
	uv = (2*real(conj(fld.u)*fld.v)).mean(axis=0)
	usig = (2*real(conj(fld.u)*fld.sig)).mean(axis=0)
	
	figure()
	semilogy(range(fld.NC),u2,'-x')
	xlabel('k (2 $\pi$ / Ly)')
	ylabel('|u|^2')
	
	figure()
	semilogy(range(fld.NC),v2,'-x')
	xlabel('k (2$\pi$ / Ly)')
	ylabel('|v|^2')

	figure()
	semilogy(range(fld.NC),abs(sig2),'-x')
	xlabel('k (2 $\pi$ / Ly)')
	ylabel('|sig|^2')

	figure()
	semilogy(range(fld.NC),abs(uv),'-x')
	xlabel('k (2 $\pi$ / Ly)')
	ylabel('$Re(u^* v)$')

	figure()
	semilogy(range(fld.NC),abs(usig),'-x')
	xlabel('k (2 $\pi$ / Ly)')
	ylabel('$Re(u^* \sigma)$')
	return u2,v2,sig2,uv,usig
		

def amf(fld,nu,mp,xs):
	vxbar = fld.u[:,0]
	vybar = fld.v[:,0]
	dbar = fld.sig[:,0]	
	
	x = fld.x
	up = zeros(fld.u.shape)
	vp = zeros(fld.v.shape)
	sigp = zeros(fld.sig.shape)
	
	up[:,1:] = fld.u[:,1:]
	vp[:,1:] = fld.v[:,1:]
	sigp[:,1:] = fld.sig[:,1:] 
	
	(Nx,NC) = fld.u.shape
	norm = fld.Ny
	
	

	dxu,_ = gradient(fld.u,fld.dx,1)
	dxv,_ = gradient(fld.v,fld.dx,1)
	dxsig,_ = gradient(fld.sig,fld.dx,1)
	
	dyu = fld.u * 1j*fld.kk
	dyv = fld.v * 1j*fld.kk
	dysig = fld.sig * 1j*fld.kk
	
	dxvxbar = dxu[:,0]
	dxvybar = dxv[:,0]
	dxdbar = dxsig[:,0]
	
	dxu[:,0] = 0
	dxv[:,0] = 0
	dxsig[:,0] = 0
	
	Txy = nu*(dxv + dyu) 
	Tyy = nu*(2*dyv - (2./3)*(dxu+dyv))
	
	rTxy = fft.irfft(Txy)
	rTyy = fft.irfft(Tyy)
	vxp = fft.irfft(up)
	dxvxp = fft.irfft(dxu)
	vyp = fft.irfft(vp)
	dxvyp = fft.irfft(dxv)
	dp  = fft.irfft(sigp)
	dxdp = fft.irfft(dxsig)
	
	rPixy = fft.irfft(fld.sig)*rTxy
	rPiyy = -fft.irfft(fld.sig) + fft.irfft(fld.sig)*rTyy
	
	Pixy = fft.rfft(rPixy)/norm
	Piyy = fft.rfft(rPiyy)/norm
	
	dxPixy,_ = gradient(Pixy,fld.dx,1)
	dyPiyy = 1j*fld.kk*Piyy
	
	Pixybar = Pixy[:,0]
	dxPixybar = dxPixy[:,0]
	
	divPip = dxPixy + dyPiyy
	divPip[:,0] = 0
	
	rdivPip = fft.irfft(divPip)
	
	Twd = (dp*vxp).mean(axis=1)*(dxvybar+2) - dbar*(vxp*dxvyp).mean(axis=1) \
			+ (dp*rdivPip).mean(axis=1) + dxPixybar*(dp*dp).mean(axis=1)/(dbar*dbar)

	print Nx,NC
# 	Twd = zeros(Nx)
# 	for i in range(Nx):
# 		Twd[i] = 0
# 		for j in range(NC):
# 			Twd[i] += 2*real(conj(sigp[i,j])*up[i,j]*(dxvybar[i]+2) - dbar[i]*conj(vxp[i,j])*dxvyp[i,j])

	
	Fp = dbar*(vxp*vyp).mean(axis=1) + vxbar*(dp*vyp).mean(axis=1) + (dp*vxp*vyp).mean(axis=1)
	Fb = (dbar*vxbar + (dp*vxp).mean(axis=1))*(vybar+2*x) - Pixybar
	dxFp = gradient(Fp,fld.x)
	dxFb = gradient(Fb,fld.x)
	
	yy,xx = meshgrid(fld.y,fld.x)
	Phi = -mp/sqrt(xs**2+xx**2+yy**2)
	phi = fft.rfft(Phi)/norm
	
	phi[:,0] = 0
	dyphi = 1j*fld.kk*phi
	
	dyPhip = fft.irfft(dyphi)
	
	Th = (dp*dyPhip).mean(axis=1)
	
	figure()
	plot(x,Twd,x,-(dxFp+Th))
	legend(('Twd','-LHS'))
	xlabel('x')
	title('Wave Steady State')
	
	figure()
	plot(x,Twd,x,dxFb)
	legend(('Twd','dxF'))
	xlabel('x')
	title('Viscous Steady State')
	
	return Twd,dxFp,Th, dxFb, Fp,Fb
	

			
