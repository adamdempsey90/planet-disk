#include "planetdisk.h"

double *lbuff, *rbuff;

void set_bc(Field *fld) {
/* Set the boundary conditions for the right boundary */
	int i,j;
//	double sig0 = (fld->Params->sig0);
	
	MPI_Status status;
	
/* Zero Gradient B.C */	

	if (np == 1) {
		for(i=0;i<NG;i++) {
#ifndef EXTRAP
			j = NC*(2*NG-(1+i));
			memcpy(&(fld->u[i*NC]),&(fld->u[j]),sizeof(double complex)*NC);
			memcpy(&(fld->v[i*NC]),&(fld->v[j]),sizeof(double complex)*NC);
			memcpy(&(fld->sig[i*NC]),&(fld->sig[j]),sizeof(double complex)*NC);

			j = NC*(2*(NG+Nx)-(i+Nx+NG+1));
			memcpy(&(fld->u[(i+Nx+NG)*NC]),&(fld->u[j]),sizeof(double complex)*NC);
			memcpy(&(fld->v[(i+Nx+NG)*NC]),&(fld->v[j]),sizeof(double complex)*NC);
			memcpy(&(fld->sig[(i+Nx+NG)*NC]),&(fld->sig[j]),sizeof(double complex)*NC);
#else 
			for(j=0;j<NC;j++) {
				fld->u[CINDX] = fld->u[j+NC*NG];
				fld->v[CINDX] = fld->v[j+NC*NG];
				fld->sig[CINDX] = fld->sig[j+NC*NG];

				fld->u[j+NC*(i+Nx+NG)] = fld->u[j+NC*(NG+Nx-1)];
				fld->v[j+NC*(i+Nx+NG)] = fld->v[j+NC*(NG+Nx-1)];
				fld->sig[j+NC*(i+Nx+NG)] = fld->sig[j+NC*(NG+Nx-1)];
			}
#endif
		}
	}
	else {
	if (rank==0) {
/* Proc 0 sends its right b.c to proc 1 */

		memcpy(&rbuff[0],(double *)&(fld->u[iend-NC*NG]),sizeof(double)*NR*NG);
		memcpy(&rbuff[NG*NR],(double *)&(fld->v[iend-NC*NG]),sizeof(double)*NR*NG);
		memcpy(&rbuff[2*NG*NR],(double *)&(fld->sig[iend-NC*NG]),sizeof(double)*NR*NG);
		
		MPI_Sendrecv_replace(rbuff,3*NR*NG,MPI_DOUBLE,rank+1,rank,
						rank+1,rank+1,MPI_COMM_WORLD,&status);

		memcpy((double *)&(fld->u[iend]), &rbuff[0],sizeof(double)*NR*NG);
		memcpy((double *)&(fld->v[iend]),&rbuff[NG*NR],sizeof(double)*NR*NG);
		memcpy((double *)&(fld->sig[iend]),&rbuff[2*NG*NR],sizeof(double)*NR*NG);	
		
/* Set left b.c for proc 0 */	
	
		for(i=0;i<NG;i++) {
#ifndef EXTRAP
			j = NC*(2*NG-(1+i));
			memcpy(&(fld->u[i*NC]),&(fld->u[j]),sizeof(double complex)*NC);
			memcpy(&(fld->v[i*NC]),&(fld->v[j]),sizeof(double complex)*NC);
			memcpy(&(fld->sig[i*NC]),&(fld->sig[j]),sizeof(double complex)*NC);
#else
			for(j=0;j<NC;j++) {
				fld->u[CINDX] = fld->u[j+NC*NG];
				fld->v[CINDX] = fld->v[j+NC*NG];
				fld->sig[CINDX] = fld->sig[j+NC*NG];
			}
#endif
		}
	}
	
	else {
	
	if (rank==np-1) {
/* Proc np-1 receives its left b.c from proc np-2 */

		memcpy(&lbuff[0],(double *)&(fld->u[istart]),sizeof(double complex)*NC*NG);
		memcpy(&lbuff[NG*NR],(double *)&(fld->v[istart]),sizeof(double complex)*NC*NG);
		memcpy(&lbuff[2*NG*NR],(double *)&(fld->sig[istart]),sizeof(double complex)*NC*NG);
		
		MPI_Sendrecv_replace(lbuff,3*NR*NG,MPI_DOUBLE,rank-1,rank,
					 rank-1,rank-1,MPI_COMM_WORLD,&status);

/* Copy over from buffer */

		memcpy((double *)&(fld->u[0]),&lbuff[0],sizeof(double complex)*NC*NG);
		memcpy((double *)&(fld->v[0]),&lbuff[NG*NR],sizeof(double complex)*NC*NG);
		memcpy((double *)&(fld->sig[0]),&lbuff[2*NG*NR],sizeof(double complex)*NC*NG);		

/* Set right b.c for proc np-1 */


		for(i=0;i<NG;i++) {
#ifndef EXTRAP
			j = NC*(2*(NG+Nx)-(i+Nx+NG+1));
			memcpy(&(fld->u[(i+Nx+NG)*NC]),&(fld->u[j]),sizeof(double complex)*NC);
			memcpy(&(fld->v[(i+Nx+NG)*NC]),&(fld->v[j]),sizeof(double complex)*NC);
			memcpy(&(fld->sig[(i+Nx+NG)*NC]),&(fld->sig[j]),sizeof(double complex)*NC);
#else
			for(j=0;j<NC;j++) {
				fld->u[j+NC*(i+Nx+NG)] = fld->u[j+NC*(NG+Nx-1)];
				fld->v[j+NC*(i+Nx+NG)] = fld->v[j+NC*(NG+Nx-1)];
				fld->sig[j+NC*(i+Nx+NG)] = fld->sig[j+NC*(NG+Nx-1)];
			}
#endif
		}	
	}
	else {
/* Everyone else sends right b.c to proc rank+1 
		and receives their right b.c from proc rank+1
   They also send their left b.c to proc rank-1 
   		and receive their left b.c from proc rank-1
*/

/* 		Left BC 		*/
		memcpy(&lbuff[0],(double *)&(fld->u[istart]),sizeof(double complex)*NC*NG);
		memcpy(&lbuff[NG*NR],(double *)&(fld->v[istart]),sizeof(double complex)*NC*NG);
		memcpy(&lbuff[2*NG*NR],(double *)&(fld->sig[istart]),sizeof(double complex)*NC*NG);
		
		MPI_Sendrecv_replace(lbuff,3*NR*NG,MPI_DOUBLE,rank-1,rank,
					 rank-1,rank-1,MPI_COMM_WORLD,&status);


		memcpy((double *)&(fld->u[0]), &lbuff[0],sizeof(double complex)*NC*NG);
		memcpy((double *)&(fld->v[0]),&lbuff[NG*NR],sizeof(double complex)*NC*NG);
		memcpy((double *)&(fld->sig[0]),&lbuff[2*NG*NR],sizeof(double complex)*NC*NG);		


/*		Right BC		*/
		memcpy(&rbuff[0],(double *)&(fld->u[iend-NC*NG]),sizeof(double complex)*NC*NG);
		memcpy(&rbuff[NG*NR],(double *)&(fld->v[iend-NC*NG]),sizeof(double complex)*NC*NG);
		memcpy(&rbuff[2*NG*NR],(double *)&(fld->sig[iend-NC*NG]),sizeof(double complex)*NC*NG);	
		
		MPI_Sendrecv_replace(rbuff,3*NR*NG,MPI_DOUBLE,rank+1,rank,
						rank+1,rank+1,MPI_COMM_WORLD,&status);
						
		memcpy((double *)&(fld->u[iend]), &rbuff[0],sizeof(double complex)*NC*NG);
		memcpy((double *)&(fld->v[iend]),&rbuff[NG*NR],sizeof(double complex)*NC*NG);
		memcpy((double *)&(fld->sig[iend]),&rbuff[2*NG*NR],sizeof(double complex)*NC*NG);					
						

			
	}
	}	
	}
//	 for(i=0;i<NG;i++) {
// 		for(j=0;j<NC;j++) {
// 			fld->u[CINDX] = 0;
// 			fld->v[CINDX] = 0;
// 			
// 			if(j==0) fld->sig[CINDX] = sig0;
// 			else	fld->sig[CINDX] = 0;
// 			
// 
// 		}
// 	}

	return;
}
void set_pi_bc(Field *fld) {
	int i;
	MPI_Status status;
	double qom = (fld->Params->q)*(fld->Params->omega);
	double c = (fld->Params->c)*(fld->Params->c);
	double nu = (fld->Params->nu);
	
	
	
	if (np == 1) {
	
		for(i=0;i<istart;i++) {
				fld->Tens->Pixx[i] = -c*(fld->sig[i]);
				fld->Tens->Pixy[i] = -nu*qom*(fld->sig[i]);
				fld->Tens->Piyy[i] = -c*(fld->sig[i]);	
		}
		for(i=iend;i<NTOTC;i++) {
				fld->Tens->Pixx[i] = -c*(fld->sig[i]);
				fld->Tens->Pixy[i] = -nu*qom*(fld->sig[i]);
				fld->Tens->Piyy[i] = -c*(fld->sig[i]);	
		}
	
	
	}
	else {
	if (rank==0) {
/* Proc 0 sends its right b.c to proc 1 */
		memcpy(&rbuff[0],(double *)&(fld->Tens->Pixx[iend-NC*NG]),sizeof(double complex)*NC*NG);
		memcpy(&rbuff[NG*NR],(double *)&(fld->Tens->Pixy[iend-NC*NG]),sizeof(double complex)*NC*NG);
		memcpy(&rbuff[2*NG*NR],(double *)&(fld->Tens->Piyy[iend-NC*NG]),sizeof(double complex)*NC*NG);
		
//		printf("%d sending rbc to %d\n",rank,rank+1);
		MPI_Sendrecv_replace(rbuff,3*NR*NG,MPI_DOUBLE,rank+1,rank,
						rank+1,rank+1,MPI_COMM_WORLD,&status);
//		printf("%d sent rbc to %d\n",rank,rank+1);
		memcpy(&(fld->Tens->Pixx[iend]), (double complex *)&rbuff[0],sizeof(double complex)*NC*NG);
		memcpy(&(fld->Tens->Pixy[iend]),(double complex *)&rbuff[NG*NR],sizeof(double complex)*NC*NG);
		memcpy(&(fld->Tens->Piyy[iend]),(double complex *)&rbuff[2*NG*NR],sizeof(double complex)*NC*NG);	
		
/* Set left b.c for proc 0 */		
		for(i=0;i<istart;i++) {
			fld->Tens->Pixx[i] = -c*(fld->sig[i]);
			fld->Tens->Pixy[i] = -nu*qom*(fld->sig[i]);
			fld->Tens->Piyy[i] = -c*(fld->sig[i]);	
		}
	}
	else {
	if (rank==np-1) {
/* Proc np-1 receives its left b.c from proc np-2 */
		memcpy(&lbuff[0],(double *)&(fld->Tens->Pixx[istart]),sizeof(double complex)*NC*NG);
		memcpy(&lbuff[NG*NR],(double *)&(fld->Tens->Pixy[istart]),sizeof(double complex)*NC*NG);
		memcpy(&lbuff[2*NG*NR],(double *)&(fld->Tens->Piyy[istart]),sizeof(double complex)*NC*NG);
		
//		printf("%d receiving lbc from %d\n",rank,rank-1);
		MPI_Sendrecv_replace(lbuff,3*NR*NG,MPI_DOUBLE,rank-1,rank,
					 rank-1,rank-1,MPI_COMM_WORLD,&status);
//		printf("%d received lbc from %d\n",rank,rank-1);
/* Copy over from buffer */
		memcpy(&(fld->Tens->Pixx[0]), (double complex *)&lbuff[0],sizeof(double complex)*NC*NG);
		memcpy(&(fld->Tens->Pixy[0]),(double complex *)&lbuff[NG*NR],sizeof(double complex)*NC*NG);
		memcpy(&(fld->Tens->Piyy[0]),(double complex *)&lbuff[2*NG*NR],sizeof(double complex)*NC*NG);		

/* Set right b.c for proc np-1 */
		for(i=iend;i<NTOTC;i++) {
			fld->Tens->Pixx[i] = -c*(fld->sig[i]);
			fld->Tens->Pixy[i] = -nu*qom*(fld->sig[i]);
			fld->Tens->Piyy[i] = -c*(fld->sig[i]);	
		}
	}
	else {
/* Everyone else sends right b.c to proc rank+1 
		and receives their right b.c from proc rank+1
   They also send their left b.c to proc rank-1 
   		and receive their left b.c from proc rank-1
*/


/* 		Left BC 		*/
		memcpy(&lbuff[0],(double *)&(fld->Tens->Pixx[istart]),sizeof(double complex)*NC*NG);
		memcpy(&lbuff[NG*NR],(double *)&(fld->Tens->Pixy[istart]),sizeof(double complex)*NC*NG);
		memcpy(&lbuff[2*NG*NR],(double *)&(fld->Tens->Piyy[istart]),sizeof(double complex)*NC*NG);
		
		MPI_Sendrecv_replace(lbuff,3*NR*NG,MPI_DOUBLE,rank-1,rank,
					 rank-1,rank-1,MPI_COMM_WORLD,&status);


		memcpy((double *)&(fld->Tens->Pixx[0]), &lbuff[0],sizeof(double complex)*NC*NG);
		memcpy((double *)&(fld->Tens->Pixy[0]),&lbuff[NG*NR],sizeof(double complex)*NC*NG);
		memcpy((double *)&(fld->Tens->Piyy[0]),&lbuff[2*NG*NR],sizeof(double complex)*NC*NG);		


/*		Right BC		*/
		memcpy(&rbuff[0],(double *)&(fld->Tens->Pixx[iend-NC*NG]),sizeof(double complex)*NC*NG);
		memcpy(&rbuff[NG*NR],(double *)&(fld->Tens->Pixy[iend-NC*NG]),sizeof(double complex)*NC*NG);
		memcpy(&rbuff[2*NG*NR],(double *)&(fld->Tens->Piyy[iend-NC*NG]),sizeof(double complex)*NC*NG);	
		
		MPI_Sendrecv_replace(rbuff,3*NR*NG,MPI_DOUBLE,rank+1,rank,
						rank+1,rank+1,MPI_COMM_WORLD,&status);
						
		memcpy((double *)&(fld->Tens->Pixx[iend]), &rbuff[0],sizeof(double complex)*NC*NG);
		memcpy((double *)&(fld->Tens->Pixy[iend]),&rbuff[NG*NR],sizeof(double complex)*NC*NG);
		memcpy((double *)&(fld->Tens->Piyy[iend]),&rbuff[2*NG*NR],sizeof(double complex)*NC*NG);					

	
	}
	}	
	}
//	 for(i=0;i<NG;i++) {
// 		for(j=0;j<NC;j++) {
// 			fld->u[CINDX] = 0;
// 			fld->v[CINDX] = 0;
// 			
// 			if(j==0) fld->sig[CINDX] = sig0;
// 			else	fld->sig[CINDX] = 0;
// 			
// 
// 		}
// 	}

	return;
}

void wavekillbc(Field *fld,double dt)
{
	int i;
	double R,tau,x,dtdtau;
	const double x_inf = -(fld->Params->Lx)*.5*0.8;
	const double x_sup = (fld->Params->Lx)*.5*0.8;
	const double tau0 = .1*(fld->k[1])/(fld->Params->omega);
	
	for(i=istart;i<iend;i++) {
		x = fld->xx[i-istart];
		R=0;
		if (x > x_sup) R = (x-x_sup)/(fld->x[Nx+NG-1] - x_sup);
		if (x < x_inf) R = (x_inf - x)/(x_inf - fld->x[NG]);

		R *= R;
		tau = tau0;

		if (R>0.0) {
			tau /= R; 
			dtdtau = dt/tau;
			if (fld->kk[i-istart]==0) {
#ifdef BACKEVOLVE	
				fld->u[i] = (fld->u[i])/(1+dtdtau );
				fld->v[i] = (fld->v[i])/(1+dtdtau );
				fld->sig[i] = ((fld->sig[i]) + (fld->Params->sig0)*dtdtau )/(1+dtdtau);
#endif
			}
			else {
#ifdef WAVEEVOLVE	
					fld->u[i] = (fld->u[i])/(1+dt/tau);
					fld->v[i] = (fld->v[i])/(1+dt/tau);
					fld->sig[i]=(fld->sig[i])/(1+dt/tau);
#endif
			}
		}
	
	}
	
	return;
}

void init_buff(void) {

	lbuff=(double *)malloc(sizeof(double complex)*NG*NC*3);
	if (lbuff == NULL) printf("ERROR with allocating lbuff for processor %d\n",rank);
	rbuff=(double *)malloc(sizeof(double complex)*NG*NC*3);
	if (lbuff == NULL) printf("ERROR with allocating rbuff for processor %d\n",rank);

	return;
}

void free_buff(void) {

	free(lbuff); free(rbuff);
	return;
}