
class Field():
	def __init__(self,time,np,dir=''):
		if dir != '':
			if dir[-1]!='/':
				dir += '/'
		vx,vy,dens,x,self.k,self.y,Nx,NC,splits = loadvars(np,time,dir);
		self.Nx = Nx.sum()
		self.NC = NC
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
		self.Ny = int(self.Ly/self.dx)
		self.Ny += mod(self.Ny,2)
		self.y = -self.Ly/2+(.5+arange(self.Ny))*(self.Ly/self.Ny)
		
 		dxv,_=gradient(self.v,self.dx,1)
 		self.kk=ones(self.u.shape)*2*pi/self.Ly
 		for j in range(self.NC):
 			self.kk[:,j] *= j
 		self.vort = dxv - 1j*self.kk*self.u	
	
		self.xs,self.mp,self.nu,self.h,self.q,self.om,self.c = readparams(dir)	
		self.phi = calcPhi(self)
		
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
	def plotreal(self,xlims=(0,0),ylims=(0,0),filter='None'):
		vx,vy,dens,x,y = realspace(self,filter)
# 		vortens = (fft.irfft(self.vort)+2)/dens
		
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
		
# 		figure()
# 		imshow(vortens[startx:endx,starty:endy].transpose(), aspect='auto', origin='lower', extent=exarg,interpolation='bilinear')
# 		colorbar()
# 		title('Vortensity'); xlabel('x'); ylabel('y')
		return

def calcPhi(fld):
	yy,xx=meshgrid(fld.y,fld.x)
	Phi = -fld.mp/sqrt(xx**2 + yy**2 + fld.xs**2)
	phi = fft.rfft(Phi)/fld.Ny
	return phi	
def loadvars(np,time,dir):
	vx = []
	vy = []
	dens = []
	x = []
	y = []
	k = []
	
	Nx = arange(np)
	for i in range(np):
		temp = fromfile(dir+'id'+str(i)+'/coords.dat')
		Nx[i] = int(temp[0])
		Ny = int(temp[1])
		NC = Ny/2+1
		
		x.append(temp[2:Nx[i]+2])
		if i==0:
			k = temp[2+Nx[i]:2+Nx[i]+NC]
			y = temp[-Ny:]
		
		temp = fromfile(dir+'id'+str(i)+'/vx_'+str(time)+'.dat',dtype='complex')
				

		
		vx.append(temp[2:].reshape((Nx[i],NC)))
		
		temp = fromfile(dir+'id'+str(i)+'/vy_'+str(time)+'.dat',dtype='complex')
		vy.append(temp[2:].reshape((Nx[i],NC)))
		temp = fromfile(dir+'id'+str(i)+'/dens_'+str(time)+'.dat',dtype='complex')
		dens.append(temp[2:].reshape((Nx[i],NC)))
		
	
	splits=[]
	for i in range(1,np):
		temp = 0
		for j in range(i):
			temp += Nx[j]
		splits.append(temp)
	
	return vx,vy,dens,x,k,y,Nx,NC,splits

def readparams(dir):
	
	with open(dir+'id0/params.txt','r') as f:
		for line in f.readlines():
			if 'xsoft' in line:
				xs = float(line.split('=')[-1])	
			if 'Mp' in line:
				mp = float(line.split('=')[-1])
			if 'nu =' in line:
				nu = float(line.split('=')[-1])
			if 'h =' in line:
				h = float(line.split('=')[-1])
			if 'q =' in line:
				q = float(line.split('=')[-1])
			if 'omega' in line:
				om = float(line.split('=')[-1])	
	
	c = h*om
	nu *= om*h*h
	return xs,mp,nu,h,q,om,c
		
def realspace(fld,filter='None'):
#	x=hstack((-fld.x[::-1],fld.x))
#	y=fld.y

	
	
	Ny = int(fld.Ly/fld.dx)
	Ny += (1-mod(Ny,2))
	NC = fld.u.shape[1]
	y = -fld.Ly/2 + (.5 + arange(Ny))*fld.Ly/Ny
	
	y,x = meshgrid(fld.y,fld.x)

	u = zeros((fld.Nx,fld.Ny/2+1),dtype='complex')
	v = zeros(u.shape,dtype='complex')
	sig = zeros(u.shape,dtype='complex')

	u[:,:NC] = fld.u
	v[:,:NC] = fld.v
	sig[:,:NC] = fld.sig
	
	if filter != 'None':
		Nmax = int(NC*2./3)
		mask = zeros(u.shape)
		mask[:,:Nmax] = 1
		if filter == 'Lanczos' or filter == 'lanczos':
			mask *= sinc(fld.kk*fld.Ly/(2*Nmax))
		
		elif filter == 'Cosine' or filter == 'cosine':
			mask *= .5*(1+cos(fld.kk*fld.Ly/(2*Nmax)))
		else:
			print 'Not a valid filter.'
	else:
		mask = ones(u.shape)
		
#	sig=vstack((conj(fld.sig[::-1,:]),fld.sig))
#	u = vstack((-conj(fld.u[::-1,:]),fld.u))
#	v = vstack((-conj(fld.v[::-1,:]),fld.v))
	vx = fft.irfft(u*mask)*Ny
	vy = fft.irfft(v*mask)*Ny
	dens = fft.irfft(sig*mask)*Ny
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
	for i in range(len(t)):
		print 'loading ',i
		if q=='u':
			dat[:,i] = norm*Field(t[i],np,dir=dir).u[:,k]
		elif q=='v':
			dat[:,i] = norm*Field(t[i],np,dir=dir).v[:,k]
		
		elif q=='sig':
			dat[:,i] = norm*Field(t[i],np,dir=dir).sig[:,k]
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
			
	for i in range(len(t)):
		clf()
		if k==0:
			plot(x,real(dat[:,i]))			
		else:
			plot(x,real(dat[:,i]),x,imag(dat[:,i]))
	
		if scale:
			ylim((qmin,qmax))
		if xlims!=(0,0):
			xlim(xlims)
		title('t$\Omega$='+str(t[i]*dt)+'\t k='+str(k*2*pi/temp.Ly))
		
		fig.canvas.draw()
	return
	
def animate_real2(dat,exarg,dt=1,logscale=False):
	
	dmin=dat.min()
	dmax=dat.max()
	fig=figure(figsize=(15,10))
	xlabel('x'); ylabel('y')
	for i in range(dat.shape[-1]):
		title('t = '+str(i*dt))
		imshow(dat[:,:,i].transpose(), aspect='auto', origin='lower',extent=exarg,interpolation='bilinear',vmin=dmin,vmax=dmax)
		if i==0: colorbar()
		draw()
	

	return	
def animate_real(q,t,np,dt=1,dir='',logscale=False,xlims=(0,0),ylims=(0,0)):
	fld = Field(t[0],np,dir=dir)

	exarg = (-fld.Lx/2,fld.Lx/2,-fld.Ly/2,fld.Ly/2)
	vx,vy,dens,_,_ = realspace(fld,'None')
	
	(endx,endy) = vx.shape
	(startx,starty) = (0,0)
	
	if xlims!=(0,0):
		temp = (fld.x > xlims[0]) & (fld.x < xlims[1])
		startx = next(i for i,j in enumerate(temp) if j)
		endx = next(len(temp)-i for i,j in enumerate(temp[::-1]) if j)
		exarg = xlims + exarg[-2:]
	
	if ylims!=(0,0):
		temp = (fld.y > ylims[0]) & (fld.y < ylims[1])
		starty = next(i for i,j in enumerate(temp) if j)
		endy = next(len(temp)-i for i,j in enumerate(temp[::-1]) if j)
		exarg = exarg[:2] + ylims

	dat = zeros(vx[startx:endx,starty:endy].shape + (len(t),))

	for i in range(len(t)):
		vx,vy,dens,_,_ = realspace(Field(t[i],np,dir=dir),'None')
		if q=='u':
			dat[:,:,i] = vx[startx:endx,starty:endy]
		elif q=='v':
			dat[:,:,i] = vy[startx:endx,starty:endy]
		elif q=='sig':
			dat[:,:,i] = dens[startx:endx,starty:endy]
 
		else:
			print 'Not a valid variable name'
			return
	
# 	dmin=dat.min()
# 	dmax=dat.max()
# 	fig=figure(figsize=(15,10))
# 	xlabel('x'); ylabel('y')
# 	for i in t:
# 		title(tstr[i])
# 		imshow(dat[:,:,i].transpose(), aspect='auto', origin='lower',extent=exarg,interpolation='bilinear',vmin=dmin,vmax=dmax)
# 		if i==0: colorbar()
# 		draw()
	return dat,exarg

def plotreal(fld,tstr=''):
	vx,vy,dens,x,y = realspace(fld)
	figure()
	pcolormesh(x,y,vx); colorbar()
	xlabel('x')
	title('v_x       ' + tstr)

	figure()
	pcolormesh(x,y,vy); colorbar()
	xlabel('x')
	title('\delta v_y     ' + tstr)
	
	figure()
	pcolormesh(x,y,log10(dens)); colorbar()
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
	udxv = zeros((fld.u.shape[1]))
	dxv,_ = gradient(fld.v,fld.dx,1)
	
	u2 = (2*real(conj(fld.u)*fld.u)).mean(axis=0)
	v2 = (2*real(conj(fld.v)*fld.v)).mean(axis=0)
	sig2 = (2*real(conj(fld.sig)*fld.sig)).mean(axis=0)
	uv = (2*real(conj(fld.u)*fld.v)).mean(axis=0)
	usig = (2*real(conj(fld.u)*fld.sig)).mean(axis=0)
	udxv = (2*real(conj(fld.u)*dxv)).mean(axis=0)
	
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
	
	figure()
	semilogy(range(fld.NC),abs(udxv),'-x')
	xlabel('k (2 $\pi$ / Ly)')
	ylabel('$ |u \partial_x v| $')
	
	return u2,v2,sig2,uv,usig,udxv
		

def amf(fld):
	om = fld.om
	vx,vy,dens,xx,yy = realspace(fld)
	
#	vy -= fld.q*fld.om*xx
	
	dxvx,dyvx = gradient(vx,fld.dx,fld.Ly/fld.Ny)
	dxvy,dyvy = gradient(vy,fld.dx,fld.Ly/fld.Ny)
	
	Phi = -fld.mp/sqrt(xx**2 + yy**2 + fld.xs**2)
	_,dyPhi = gradient(Phi,fld.dx,fld.Ly/fld.Ny)
	
	vxbar = real(fld.u[:,0])
	vybar = real(fld.v[:,0] - fld.q*fld.om*fld.x)
	dbar = real(fld.sig[:,0])
	
	u = zeros(vx.shape)
	v = zeros(vy.shape)
	sig = zeros(dens.shape)

	Pixy = fld.nu*(dyvx+dxvy)*dens
	Piyy = fld.nu*dens*(2*dyvy - (2./3)*(dxvx+dyvy))
	pixybar = real(Pixy.mean(axis=1))
	piyybar = real(Piyy.mean(axis=1))
	pixyp = zeros(Pixy.shape)
	piyyp = zeros(Piyy.shape)
	
	for j in range(u.shape[1]):
		u[:,j] = vx[:,j] - vxbar
		v[:,j] = vy[:,j] - vybar
		sig[:,j] = dens[:,j] - dbar
		pixyp[:,j] = Pixy[:,j] - pixybar
		piyyp[:,j] = Piyy[:,j] - piyybar

	dxu,_ = gradient(u,fld.dx,1)
	dxv,_ = gradient(v,fld.dx,1)
	
	Twd = -dbar*(u*dxv).mean(axis=1) + (dxv[:,0]+2*fld.om)*(sig*u).mean(axis=1)
	
	divP,_ = gradient(pixyp,fld.dx,1)
	_,dyP = gradient(piyyp,fld.dx,fld.Ly/fld.Ny)
	
	divP += dyP
	
	dxPixybar= gradient(pixybar,fld.dx)
	
	
	
	visc = (sig**2).mean(axis=1) * dxPixybar/(dbar**2) - (sig*divP).mean(axis=1)/dbar
	Twd += visc
	
 		
	Ftot =(dens*vx*(vy+2*fld.om*xx) - Pixy).mean(axis=1)
	dxFtot = gradient(Ftot,fld.dx)
	
	Fb = (dbar*vxbar+(sig*u).mean(axis=1))*(vybar+2*om*fld.x) - Pixy.mean(axis=1)
	dxFb = gradient(Fb,fld.dx)
	Fp = Ftot-Fb
	dxFp = dxFtot - dxFb
#	Fp = vxbar*(sig*v).mean(axis=1) + dbar*(u*v).mean(axis=1)
	
	Th = (sig*dyPhi).mean(axis=1)
	  
	figure()
	plot(fld.x,Twd,label='$T_{wd}$')
#	plot(fld.x,-dxFp,label='$\\partial_x F_p$')
#	plot(fld.x,-Th,label='$\Sigma \\partial_y \Phi$')
	plot(fld.x,dxFp+Th,label='$-\\partial_x F_p - \Sigma \\partial_y \\Phi $')
	xlabel('x')
	legend()
	
	figure()
	plot(fld.x,Twd,label='$T_{wd}$')
	plot(fld.x,dxFb,label='$\\partial_x F_b $')
	xlabel('x')
	legend()
	return Twd,dxFp,Th,Fp



# def amf(fld,xlims=(0,0)):
# 	nu = fld.nu
# 	mp = fld.mp
# 	xs = fld.xs
# 	fld.v[:,0] -= fld.q*fld.om*fld.x
# 	sigu = real(conj(fld.sig)*fld.u)
# 	sigu = 2*sigu[:,1:].sum(axis=1)
# 	
# 	sigv = real(conj(fld.sig)*fld.v)
# 	sigv = 2*sigv[:,1:].sum(axis=1)
# 	
# 	vxvy = real(conj(fld.u)*fld.v)
# 	vxvy = 2*vxvy[:,1:].sum(axis=1)
# 	
# 	dxv,_ = gradient(fld.v,fld.dx,1)
# 	dxu,_ = gradient(fld.u,fld.dx,1)
# 	dyu = 1j*fld.kk*fld.u
# 	dyv = 1j*fld.kk*fld.v
# 	
# 	Txy = nu*(dxv+dyu)
# 	Tyy = nu*(2*dyv-(2./3)*(dxu+dyv))
# 	
# 	Pixy = zeros(fld.sig.shape)
# 	Piyy = zeros(fld.sig.shape)
# 	temp = real(conj(fld.sig)*Txy)
# 	Pixy[:,0] = 2*temp[:,1:].sum(axis=1)
# 	
# 	temp = real(conj(fld.sig)*Tyy)
# 	Piyy[:,0] = 2*temp[:,1:].sum(axis=1)
# 	
# 	udxv = real(conj(fld.u)*dxv)
# 	udxv = 2*udxv[:,1:].sum(axis=1)
# 	
# #	Pixy=fft.rfft(fft.irfft(Txy)*fft.irfft(fld.sig)*fld.Ny)/fld.Ny
# #	Piyy=fft.rfft(fft.irfft(Tyy)*fft.irfft(fld.sig)*fld.Ny)/fld.Ny
# 	
# 	for j in range(Pixy.shape[0]):
# 		Pixy[j,1:] = Txy[j,0]*fld.sig[j,1:] + Txy[j,1:]*fld.sig[j,0]
# 		Piyy[j,1:] = Tyy[j,0]*fld.sig[j,1:] + Tyy[j,1:]*fld.sig[j,0]
# 	dxPixy,_ = gradient(Pixy,fld.dx,1)
# 	dyPiyy = 1j*fld.kk*Piyy
# 	
# 	Tviscp = real(conj(fld.sig)*(dxPixy+Piyy))
# 	Tviscp = 2*Tviscp[:,1:].sum(axis=1)
# 	Tviscp /= real(fld.sig[:,0])
# 	
# 	Tviscb = real(conj(fld.sig)*fld.sig)
# 	Tviscb = 2*Tviscb[:,1:].sum(axis=1)
# 	Tviscb *= real(dxPixy[:,0])/(fld.sig[:,0]**2)
# 	
# 	
# 	Twd = sigu*(real(dxv[:,0])+2*fld.om) -fld.sig[:,0]*udxv + Tviscp + Tviscb
# 	
# 	
# 	if xlims!=(0,0):
# 		inds = (fld.x > xlims[0]) & (fld.x < xlims[1])
# 	else:
# 		inds = ones(fld.x.shape).astype('bool')
# 		
# 	figure()
# 	plot(fld.x[inds],(sigu*real(dxv[:,0]+2))[inds],label='$\langle \Sigma\' u\' \\rangle (\partial_x \langle v_y \\rangle + 2\Omega)$') 
# 	plot(fld.x[inds],(-fld.sig[:,0]*udxv)[inds],label='$-\langle u\' \partial_x v\' \\rangle \langle \Sigma \\rangle$')
# 	plot(fld.x[inds],(Tviscp+Tviscb)[inds],label='$\langle \\frac{\\nabla \cdot (\Sigma T)}{\Sigma} \\rangle$') 
# 	plot(fld.x[inds],Twd[inds],label='$T_{wd}$')
# 	legend()
# 	
# 		
# 	Fp = fld.sig[:,0]*vxvy + fld.u[:,0]*sigv
# 	
# 	
# 	Fb = (fld.v[:,0]+2*fld.x)*(fld.sig[:,0]*fld.u[:,0] + sigu) - Pixy[:,0]
# 			
# 	vx,vy,dens,x,y = realspace(fld)
# 
# 	Phi = -mp/sqrt(xs**2 + x**2 + y**2)
# 	phi = fft.rfft(Phi)/fld.Ny
# 	
# 	Th = real(conj(fld.sig)*phi*1j*fld.kk)
# 	Th = 2*Th[:,1:].sum(axis=1)
# 	
# 	sig = zeros(dens.shape)
# 	u = zeros(vx.shape)
# 	v = zeros(vy.shape)
# 	dbar = dens.mean(axis=1)
# 	vxbar = vx.mean(axis=1)
# 	vybar = vy.mean(axis=1)
# 	
# 	for i in range(dens.shape[0]):
# 		sig[i,:] = dens[i,:] - dbar[i]
# 		u[i,:] = vx[i,:] - vxbar[i]
# 		v[i,:] = vy[i,:] - vybar[i]
# 			
# 	trip = (sig*u*v).mean(axis=1)
# 	
# 	Fp += trip
# 	
# 	dxFp = gradient(Fp,fld.dx)
# 	dxFb = gradient(Fb,fld.dx)
# 	
# 	figure()
# 	plot(fld.x[inds],Twd[inds],label='Twd')
# 	plot(fld.x[inds],(-dxFp-Th)[inds],label='-(dxFp+Th)')
# 	title('Wave Steady State')
# 	xlabel('x')
# 	legend()
# 	
# 		
# 		
# 	figure()
# 	plot(fld.x[inds],Twd[inds],label='Twd')
# 	plot(fld.x[inds],dxFb[inds],label='dxFb')
# 	title('Viscous Steady State')
# 	xlabel('x')
# 	legend()
# 
# 		
# 	fld.v[:,0] += fld.q*fld.om*fld.x
# 
# 	return Twd,dxFp,Th,dxFb

def plotdenswake(fld,xvals,filter='None'):
	x = fld.x
	y = fld.y
	_,_,dens,_,_=realspace(fld,filter)
	ind = [where((x <i+fld.dx) & (x>i-fld.dx))[0][0] for i in xvals]

	figure()
	for i in ind:
   		ysh = .75*x[i]**2
		while abs(ysh > fld.Ly/2):
			ysh -= fld.Ly
		
		
		plot(y[::-1]-ysh,(dens[i,:]-1)/(fld.mp*sqrt(x[i])),label='$x=$'+str(round(x[i],2)))
	xlim(-6,6)
	legend()
	xlabel('$y-.75 x^2$',fontsize='large')
	ylabel('$\Delta\Sigma / \\sqrt{x}$',fontsize='large')
	return
	
def conv(input1,input2,linear=True):
	Nx, NC = input1.shape
	Ny = 2*NC-2
	if ~linear:
		mask = zeros((Nx,NC))
		Nmax = (2./3)*NC
		mask[:,:Nmax] = 1
	else:
		mask = ones((Nx,NC))
		Nmax = NC-1
	out = fft.rfft(fft.irfft(input1*mask)*fft.irfft(input2*mask))*Ny
	return out

	
def apply_filter(fld,filter='None',overwrite=False):
	vx,vy,dens,_,_=realspace(fld,filter)
	u = fft.rfft(vx)
	v = fft.rfft(vy)
	sig = fft.rfft(dens)/fld.Ny
	
	if overwrite:
		fld.u=u; fld.sig=sig; fld.v=v
		return
	else:
		return u,v,sig	
		
	