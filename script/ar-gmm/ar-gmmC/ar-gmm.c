/*
 * =====================================================================================
 *
 *       Filename:  ar-gmm.c
 *
 *    Description:  -m mixNum 
 *    							-l dimension 
 *    							-o order 
 *    							-I initial model file 
 *    							-O initial model order
 *    									-1:0 order with zero Transformation matrix
 *    									k>=0: k order
 *    							-t 收敛条件 
 *    							-v vfloor 
 *    							-T dataType    0 prosody 1 spectrum
 *    							-N maxItrationNum
 *									datafile   training data
 *									outFile    output model 
 *
 *        Version:  1.0
 *        Created:  2013年08月01日 20时11分27秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Li Xian (), shysian@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <cblas.h>
//#include "ar-gmm_sub.h"

#define PI 3.1415926
#define ABS(x) ( (x)>0 ?(x):(-(x)) )

//enum DataType{PROSODY,SPECTRUM};
typedef struct _FVECTOR
{
	int len;
	float *data;
}FVECTOR;

typedef struct _FMATRIX
{
	int row;
	int col;
	float *data;
}FMATRIX;


typedef struct _JGAUSS{
	float weight;
	float gconst;
	FVECTOR miu;
	FVECTOR var;
	FMATRIX ar;
	FMATRIX bt;
}JGAUSS;

typedef struct _GMM{
	JGAUSS *gauss;
	int M;
	int L;
	int order;
	
} GMM;

typedef struct _ACCML
{
	double data;
	double *pdata;
}ACCML;

typedef struct _ACCM
{
  ACCML * dim;
}ACCM;

typedef struct _ACC
{
	int M;
	int L;
	int order;
	int type;
	ACCM *mix;	
}ACC;



void getmemACC(ACC *acc)
{
	int m,l,iO;
	int M=acc->M,L=acc->L,order=acc->order,type=acc->type;
	acc->mix=(ACCM *)calloc(M,sizeof(ACCM));
	for(m=0;m<M;m++)
	{
		acc->mix[m].dim=(ACCML *)calloc(L,sizeof(ACCML));
		for(l=0;l<L;l++)
		{
			if(type==0)
				;
			else if(type==1)
				acc->mix[m].dim[l].pdata=(double *)calloc(order,sizeof(double));
			else if(type==2)
				acc->mix[m].dim[l].pdata=(double *)calloc((order)*(order),sizeof(double));			
		}
	}
}

void freeACC(ACC *acc)
{
	
	int m,l,iO;
	int M=acc->M,L=acc->L,order=acc->L,type=acc->type;
	for(m=0;m<M;m++)
	{
		for(l=0;l<L;l++)
		{
			if(type==0)
				;
			else 
				free(acc->mix[m].dim[l].pdata);
		}
		free(acc->mix[m].dim);
	}
	free(acc->mix);
}

void zfillf(float *a,int L)
{
	memset(a,0,sizeof(float)*L);
}
void zfilld(double *a,int L)
{
	memset(a,0,sizeof(double)*L);
}
void zfillM(FMATRIX *matrix)
{
	int m=matrix->row,n=matrix->col;
	zfillf(matrix->data,m*n);
}

void zfillACC(ACC *acc)
{
	int m,l,iO;
	int M=acc->M,L=acc->L,order=acc->order,type=acc->type;
	for(m=0;m<M;m++)
	{
		for(l=0;l<L;l++)
		{
			if(type==0)
				acc->mix[m].dim[l].data=0.0;
			else if(type==1)
				zfilld(acc->mix[m].dim[l].pdata,order);
			else if(type==2)
				zfilld(acc->mix[m].dim[l].pdata,(order)*(order));
//			else if(type==3)
//				zfilld(acc->mix[m].dim[l].pdata,order+1+L);
//			else if(type==4)
//				zfilld(acc->mix[m].dim[l].pdata,(order+1+L)*(order+1+L));
		}
	}
}


void getmemV(FVECTOR *v)
{
	int len=v->len;
	v->data=(float *)calloc(len,sizeof(float));
}
void freeV(FVECTOR *v)
{
	free(v->data);
}
void getmemM(FMATRIX *matrix)
{
	int m=matrix->row;
	int n=matrix->col;
	matrix->data=(float *)calloc(m*n,sizeof(float ));
}

void freeM(FMATRIX *matrix)
{
	free(matrix->data);
}

void freadM(FMATRIX *matrix,FILE *fdata)
{
	int m=matrix->row,n=matrix->col;
	fread(matrix->data,sizeof(float),m*n,fdata);
}

void fwriteM(FMATRIX *matrix,FILE *fdata)
{
	int m=matrix->row,n=matrix->col;
	fwrite(matrix->data,sizeof(float),m*n,fdata);
}

void fwriteGMM(GMM *gmm,FILE *fgmm)
{
	int M=gmm->M,L=gmm->L,order=gmm->order;
	int m,l,iO;
	for(m=0;m<M;m++)
	{
		//fwrite(&(gmm->gauss[m].gconst),sizeof(float),1,fgmm);
		fwrite(&(gmm->gauss[m].weight),sizeof(float),1,fgmm);
	}
	for(m=0;m<M;m++)
	{
		fwrite(gmm->gauss[m].miu.data,sizeof(float),L,fgmm);
//		for(iO=0;iO<order;iO++)
//		{
//			fwrite(gmm->gauss[m].ar.data,sizeof(float),L,fgmm);
//		}
		fwriteM(&(gmm->gauss[m].ar),fgmm);
		fwrite(gmm->gauss[m].var.data,sizeof(float),L,fgmm);
		fwriteM(&(gmm->gauss[m].bt),fgmm);
	}
}

void getmemGMM(GMM *gmm,int M,int L,int order)
{
	int m;
	gmm->M=M;
	gmm->L=L;
	gmm->order=order;
	gmm->gauss=(JGAUSS *)malloc(sizeof(JGAUSS)*M);
	for (m=0;m<M;m++)
	{
		gmm->gauss[m].miu.len=L;
		gmm->gauss[m].var.len=L;
		gmm->gauss[m].ar.row=order;
		gmm->gauss[m].ar.col=L;
		gmm->gauss[m].bt.row=L/2;
		gmm->gauss[m].bt.col=L/2;
		getmemV(&(gmm->gauss[m].miu));
		getmemV(&(gmm->gauss[m].var));
		getmemM(&(gmm->gauss[m].ar));
		getmemM(&(gmm->gauss[m].bt));
	}

}

void freeGMM(GMM *gmm)
{
	int m;
	for(m=0;m<gmm->M;m++)
	{
		freeV(&(gmm->gauss[m].miu));
		freeV(&(gmm->gauss[m].var));
		freeM(&(gmm->gauss[m].ar));
		freeM(&(gmm->gauss[m].bt));
	}
	free(gmm->gauss);
}

void cal_gconst(JGAUSS *gauss)
{
	int L=gauss->miu.len;
	double det=1.0;
	int i;
	for (i=0;i<L;i++)
	{
	
		det*=(double)gauss->var.data[i];
	}
	//printf("det:%1.20f\n",det);
	gauss->gconst=log(gauss->weight)-0.5*log(det)-L/2*log(2.0*PI);
}


void readInitModel(FILE *fgmm,GMM *gmm,int initOrder)
{
	int m,M=gmm->M;
	int l,L=gmm->L;
	fseek(fgmm,0L,SEEK_END);
	int length=ftell(fgmm)/sizeof(float);
	float *data=(float *)malloc(length*sizeof(float));
	fseek(fgmm,0L,SEEK_SET);
	fread(data,sizeof(float),length,fgmm);

	for(m=0;m < gmm->M;m++)
	{
		gmm->gauss[m].weight=data[m];
	}
  if(initOrder==-1)
	{
		for(m=0;m< gmm->M;m++)
		{
			memcpy(gmm->gauss[m].miu.data,data+gmm->M+m*2*L,sizeof(float)*(gmm->L));
			memcpy(gmm->gauss[m].var.data,data+gmm->M+m*2*L+L,sizeof(float)*(gmm->L));
			cal_gconst(gmm->gauss+m);
			printf("gconst[m]:%f\n",gmm->gauss[m].gconst);
		}
	}
	else
	{
		float *tempp=data+gmm->M;
		for(m=0;m < gmm->M;m++)
		{
			memcpy(gmm->gauss[m].miu.data,tempp,sizeof(float)*(gmm->L));
			tempp+=L;
			memcpy(gmm->gauss[m].ar.data,tempp,sizeof(float)*(gmm->L)*initOrder);
			tempp+=L*initOrder;
			memcpy(gmm->gauss[m].var.data,tempp,sizeof(float)*(gmm->L));
			tempp+=L;
			memcpy(gmm->gauss[m].bt.data,tempp,sizeof(float)*(gmm->L/2)*(gmm->L/2));
			tempp+=L/2*L/2;
		}
		
	}
}
/* void fwriteM(FMATRIX *matrix,FILE *fdata)
 * {
 * 	int m=matrix->row,n=matrix->col;
 * 	fwrite(matrix->data,sizeof(float),m*n,fdata);
 * }
 * 
 * void fwriteGMM(GMM *gmm,FILE *fgmm)
 * {
 * 	int M=gmm->M,L=gmm->L,order=gmm->order;
 * 	int m,l,iO;
 * 	for(m=0;m<M;m++)
 * 	{
 * 		//fwrite(&(gmm->gauss[m].gconst),sizeof(float),1,fgmm);
 * 		fwrite(&(gmm->gauss[m].weight),sizeof(float),1,fgmm);
 * 	}
 * 	for(m=0;m<M;m++)
 * 	{
 * 		fwrite(gmm->gauss[m].miu.data,sizeof(float),L,fgmm);
 * //		for(iO=0;iO<order;iO++)
 * //		{
 * //			fwrite(gmm->gauss[m].ar.data,sizeof(float),L,fgmm);
 * //		}
 * 		fwriteM(&(gmm->gauss[m].ar),fgmm);
 * 		fwrite(gmm->gauss[m].var.data,sizeof(float),L,fgmm);
 * 		fwriteM(&(gmm->gauss[m].bt),fgmm);
 * 	}
 * }
 */

double cal_postProb(float *miu,float *x,JGAUSS *gauss)
{
	int L=gauss->miu.len;
	int i;
	double postprob=gauss->gconst;
	for(i=0;i<L;i++)
	{
		postprob+=-0.5*(x[i]-miu[i])*(x[i]-miu[i])/gauss->var.data[i];
	}
	return postprob;
}

//from SPTK _gmm.c
#define LZERO  (-1.0E10)        /* ~log(0) */
#define LSMALL (-0.5E10) 
double log_add(double logx,double logy)
{
   double swap, diff, minLogExp, z;

   if (logx < logy) {
      swap = logx;
      logx = logy;
      logy = swap;
   }

   diff = logy - logx;
   minLogExp = -log(-LZERO);

   if (diff < minLogExp)
      return ((logx < LSMALL) ? LZERO : logx);
   else {
      z = exp(diff);
      return (logx + log(1.0 + z));
   }
}

// blas
void matrixMultiply(float* A,int ra,int ca,float * B,int cb,float *C)
{ 
		const enum CBLAS_ORDER Order=CblasRowMajor;
    const enum CBLAS_TRANSPOSE TransA=CblasNoTrans;
    const enum CBLAS_TRANSPOSE TransB=CblasNoTrans;
    const double alpha=1;
    const double beta=1;
    const int lda=ca;//A的列
    const int ldb=cb;//B的列
    const int ldc=cb;//C的列
		cblas_sgemm(Order,TransA,TransB,ra,cb,ca,alpha,A,lda,B,ldb,beta,C,ldc);
}
//lapack 
//solve linear equation
int solve(double *A,double *B,int N)
{
	 int *ipiv=(int *)malloc(sizeof(int)*N);
	 const enum CBLAS_ORDER Order=CblasRowMajor;
	 int info=clapack_dgesv(Order,N,1,A,N,ipiv,B,N);
	 free(ipiv);
	 return info;	 
}


int main(int argc,char **argv)
{
	int M=12;
	int order=1;
	int L=16,ITE=20;
	float thr=1.0E-2,vfloor=0.001;
	FILE *fgmm;
	int dataType,initOrder=-1;
	GMM gmm;
//	gmm.gauss=(JGAUSS *)malloc(sizeof(JGAUSS)*M);
//	gmm.M=M;
//	gmm.L=L;
//	for (m=0;m<M;m++)
//	{
//		gmm.gauss[m].miu.len=L;
//		gmm.gauss[m].var.len=L;
//		gmm.gauss[m].ar.row=order;
//		gmm.gauss[m].ar.col=L;
//		gmm.gauss[m].bt.row=L/2;
//		gmm.gauss[m].bt.col=L/2;
//		getmemV(&(gmm.gauss[m].miu));
//		getmemV(&(gmm.gauss[m].var));
//		getmemM(&(gmm.gauss[m].ar));
//		getmemM(&(gmm.gauss[m].bt));
//	}
  while(--argc>2)
			if(**++argv=='-')
					switch(*(*argv+1))
					{
						case 'm':
										M=atoi(*++argv);
										--argc;
										break;
						case 'l':
										L=atoi(*++argv);
										--argc;
										break;
						case 'o':
										order=atoi(*++argv);
										--argc;
										break;
						case 't':
										thr=atof(*++argv);
										--argc;
										break;
						case 'v':
										vfloor=atof(*++argv);
										--argc;
										break;
						case 'I':
										fgmm=fopen(*++argv,"rb");
										--argc;
										break;
						case 'O':
										initOrder=atoi(*++argv);
										--argc;
										break;
						case 'T':
										dataType= atoi(*++argv);
										--argc;
										break;
						case 'N':
										ITE= atoi(*++argv);
										--argc;
										break;
										
					}

	getmemGMM(&gmm,M,L,order);
//	FILE *fgmm=fopen("neutral_sad.tone.gmm.f","rb");
	readInitModel(fgmm,&gmm,initOrder);
	fclose(fgmm);
  
/* 	double *testD=(double *)malloc(sizeof(double)*T*L);
 * 	fread(testD,sizeof(double),T*L,fdata);
 * 	for(m=0;m<T*L;m++)
 * 	{
 * 		printf("%f\n",testD[m]);
 * 	}
 * 	exit(0);
 */
	int ite=0;
	int iT,jT,iSen=0,iSeg=0,iO,jO,l,m,iL;

	FILE *fdata;
	int T;
	char p_order[10];
	char dataFile[400];
	FMATRIX trainData;
  FMATRIX *preData=(FMATRIX *)calloc(order,sizeof(FMATRIX));
	if(dataType == 0)
	{
		printf("Prosody\n");
		fdata=fopen(*++argv,"rb");

		fseek(fdata,0,SEEK_END);
		T=ftell(fdata)/sizeof(float)/L;
		fseek(fdata,0,SEEK_SET);
		printf("T=%d\n",T);
		trainData.row=T;
		trainData.col=L;
  	getmemM(&trainData);
		freadM(&trainData,fdata);
//		for(iO=0;iO<order,iO++)
//		{
//			preData[iO].row=T;
//			preData[iO].col=L;
//			for(iT=0;iT<T;iT++)
//			{
//				if(trainData[iT*L][0]==0)
//				{
//					for(jT=iT+1;iT<T;jT++)
//					{
//						if(trainData[jT*L][0]==0)
//					}
//				}
//			}
			// }
		}

	else if(dataType == 1)
	{
		printf("Spectrum\n");
		strcpy(dataFile,*++argv);
		sprintf(p_order,".%d-ar",0);
		strcat(dataFile,p_order);
		fdata=fopen(dataFile,"rb");
		fseek(fdata,0,SEEK_END);
		T=ftell(fdata)/sizeof(float)/L;
		fseek(fdata,0,SEEK_SET);
		printf("T=%d\n",T);
		trainData.row=T;
		trainData.col=L;
  	getmemM(&trainData);
		freadM(&trainData,fdata);
		fclose(fdata);
		for(iO=1;iO<=order;iO++)
		{
			
			strcpy(dataFile,*argv);
			sprintf(p_order,".%d-ar",iO);
			strcat(dataFile,p_order);
			fdata=fopen(dataFile,"rb");
			preData[iO-1].row=T;
			preData[iO-1].col=L;
			getmemM(preData+iO-1);
			freadM(preData+iO-1,fdata);
			fclose(fdata);
		}
	}


	
	

	// acc 

	ACC A;
	A.M=M;A.L=L/2;A.order=order+1;A.type=2;
	getmemACC(&A);
	ACC C;
	C.M=M;C.L=L/2;C.order=order+1;C.type=1;
  getmemACC(&C);
	ACC B;
	B.M=M;B.L=L/2;B.order=order+1+L/2;B.type=2;
	getmemACC(&B);
	ACC D;
	D.M=M;D.L=L/2;D.order=order+1+L/2;D.type=1;
  getmemACC(&D);
	
	ACC accZZ; // z*z
	accZZ.M=M;accZZ.L=L;accZZ.type=0;
	getmemACC(&accZZ);
	ACC accZ; // z
	accZ.M=M;accZ.L=L;accZ.type=0;
	getmemACC(&accZ);
	ACC accPC; //previous * current
	accPC.M=M;accPC.L=L;accPC.order=order;accPC.type=1;
	getmemACC(&accPC);
	ACC accXY; // x *y
	accXY.M=M;accXY.L=1;accXY.order=L/2;accXY.type=2;
	getmemACC(&accXY);
 	double *accW=(double *)calloc(M,sizeof(double));
//	ACC accX; // X*others
//	accX.M=M;accX.L=L/2;accX.order=order+1;accX.type=2;
//	getmemACC(&accX);
//	ACC accY; // X*others
//	accY.M=M;accY.L=L/2;accY.order=order+1+L/2;accY.type=2;
//	getmemACC(&accY);
//
//	double *tempx=(double *)calloc(order+1,sizeof(double));
//	double *tempy=(double *)calloc(order+1+L/2,sizeof(double));
//	double *tempC=(double *)calloc(L/2*L/2,sizeof(double));
		
  FMATRIX curMiu;
	curMiu.row=M;
	curMiu.col=L;
	getmemM(&curMiu);
	double *postProb=(double *)calloc(M,sizeof(double));	
	int maxMix=0;
	double maxProb;
	double sumProb;
	double likliProb;
	double preLikli;
	float *pd;

	float **pre,*pre1,*pre2;
	pre=(float **)malloc(sizeof(float*)*order);
	double *vFloor=(double *)malloc(sizeof(double)*L);
	double det;
	double *pA,*pB,*pC,*pD,*pXY;
	ACCML *pPC,*pZZ,*pZ;
	JGAUSS *pGauss;
	clock_t start,finish;
	double temp0,temp1,temp2;
 	int templ1,templ2;
	//start
	while(ite<ITE)
	{
		for(iSeg=0,iSen=0,iT=0,jT=0,likliProb=0.0,pd=trainData.data;iT<T;iT++,jT++,pd+=L)
		{
		  if(pd[0]==0)
			{
				//printf("%dth sentence\n",++iSen);
				fflush(stdout);
				jT=0;
				finish=clock();
			//	printf("%fsecond\n",(double)(finish-start)/CLOCKS_PER_SEC);
				start=clock();
				continue;
			}
			for(iO=0;iO<order;iO++)
			{
							if(dataType==0)
							{
								if(jT>iO)
						  		pre[iO]=trainData.data+(iT-iO-1)*L;
								else
									pre[iO]=trainData.data+(iT-jT)*L;
							}
							else if(dataType=1)
								pre[iO]=preData[iO].data+iT*L;
			}
		//	zfillM(&curMiu);	
			for(m=0,sumProb=-1.0E10,maxProb=-1.0E10;m<M;m++)
			{
				
				//matrixMultiply(gmm.gauss[m].bt.data,L/2,L/2,pd,1,curMiu.data+m*L+L/2);
				for(l=0;l<L;l++)
				{
						curMiu.data[m*L+l]=gmm.gauss[m].miu.data[l];
					  if(l>=L/2)
						{
							for(iL=0;iL<L/2;iL++)
							curMiu.data[m*L+l] += gmm.gauss[m].bt.data[(l-L/2)*L/2+iL]*pd[iL];
						}	
						for (iO=0;iO<order;iO++)
						{
							//printf("m=%d,i=%d,iO=%d\n",m,i,iO);
						  //	pre=trainData.data+(iT-iO-1)*L;
								//if(jT>
//							if(dataType==0)
//							{
//								if(jT>iO)
//						  		pre=trainData.data+(iT-iO-1)*L;
//								else
//									pre=trainData.data+(iT-jT)*L;
//							}
//							else if(dataType=1)
//								pre=preData[iO].data+iT*L;
							curMiu.data[m*L+l] +=gmm.gauss[m].ar.data[iO*L+l] * pre[iO][l]; //pre
						}
				}
				postProb[m]=cal_postProb(curMiu.data+m*L,pd,gmm.gauss+m);	
			
				if(ite==0)
				{
					if(postProb[m]>maxProb)
					{
						maxMix=m;
						maxProb=postProb[m];
					}
				}	

				sumProb=log_add(sumProb,postProb[m]);
			}
			likliProb +=sumProb;
			if(ite==0)
			{
				for(l=0;l<L;l++)
						//vFloor[l]+=pow((pd[l]-gmm.gauss[maxMix].miu.data[l]),2.0);
						vFloor[l]+=pow((pd[l]-curMiu.data[m*L+l]),2.0);
			}
			for(m=0;m<M;m++)
			{
				temp0=exp(postProb[m]-sumProb);
				accW[m]+=temp0;
				pPC=accPC.mix[m].dim;
				pZZ=accZZ.mix[m].dim;
				pZ=accZ.mix[m].dim;
				pXY=accXY.mix[m].dim[0].pdata;
					for(l=0;l<L/2;l++)
							{
									templ1=l+L/2;
									templ2=l*L/2;
									temp1=temp0*pd[l];
									temp2=temp0*pd[templ1];
									pA=A.mix[m].dim[l].pdata;
									pB=B.mix[m].dim[l].pdata;
									pC=C.mix[m].dim[l].pdata;
									pD=D.mix[m].dim[l].pdata;
									pZZ[l].data+=temp1*pd[l];
									//pZZ[l].data+=postProb[m]*(double)pd[l]*(double)pd[l];
									pZZ[templ1].data+=temp2*pd[templ1];
									//pZZ[L/2+l].data+=postProb[m]*(double)pd[templ1]*(double)pd[templ1];
									pZ[l].data+=temp1;
									pZ[templ1].data+=temp2;
									
									for(iO=0;iO<B.order;iO++)
									{
										for(jO=0;jO<=iO;jO++)
										{
											if(iO>order && jO>order)
											{
												pB[iO*(B.order)+jO]+=temp0* pd[iO-A.order]*pd[jO-A.order];
											}

											else if(iO>order && jO<order)
											{
												//pre=trainData.data+(iT-jO-1)*L;
//															if(dataType==0)
//															{
//																if(jT>jO)
//																	pre=trainData.data+(iT-jO-1)*L;
//																else
//																	pre=trainData.data+(iT-jT)*L;
//															}
//															else if(dataType=1)
//																pre=preData[jO].data+iT*L;
//
												pB[iO*(B.order)+jO]+=temp0* pd[iO-A.order]*pre[jO][templ1];
											}
											else if(iO>order && jO==order)
											{
												pB[iO*(B.order)+jO]+=temp0* pd[iO-A.order];
											}
											
											else if(iO<order)
											{
//												if(dataType==0)
//												{
//													if(jT-iO>0)
//														pre1=trainData.data+(iT-iO-1)*L;
//													else
//														pre1=trainData.data+(iT-jT)*L;
//													if(jT-jO>0)
//														pre2=trainData.data+(iT-jO-1)*L;
//													else 
//														pre2=trainData.data+(iT-jT)*L;
//												}
//												else if(dataType==1)
//												{				
//													pre1=preData[iO].data+iT*L;
//													pre2=preData[jO].data+iT*L;
//												}
												pA[iO*A.order+jO]+=temp0* pre[iO][l]*pre[jO][l];
										//		printf("mix:%d,dim:%d,iO:%d,jO:%d,iT:%d\n",m,l,iO,jO,iT);
												pB[iO*(B.order)+jO]+=temp0* pre[iO][templ1]*pre[jO][templ1];
											}
											else if(iO==order && jO<order)
											{

												//pre=trainData.data+(iT-jO-1)*L;
										//		pre=preData[jO].data+iT*L;
//															if(dataType==0)
//															{
//																if(jT>jO)
//																	pre=trainData.data+(iT-jO-1)*L;
//																else
//																	pre=trainData.data+(iT-jT)*L;
//															}
//															else if(dataType=1)
//																pre=preData[jO].data+iT*L;
												pA[iO*A.order+jO]+=temp0* pre[jO][l];
												pB[iO*(B.order)+jO]+=temp0* pre[jO][templ1];
											}
											else if(iO==order && jO==order)
											{
												pA[iO*(A.order)+jO]+=temp0;
												pB[iO*(B.order)+jO]+=temp0;
											}
									
										}// for jO
										if(iO>order)
										{
											pXY[templ2+(iO-A.order)]+=temp2*pd[(iO-A.order)];
											pD[iO]+=temp2*pd[iO-A.order];
										}

										else if(iO<order)
										{
											// pre=trainData.data+(iT-iO-1)*L;
//															if(dataType==0)
//															{
//																if(jT>iO)
//																	pre=trainData.data+(iT-iO-1)*L;
//																else
//																	pre=trainData.data+(iT-jT)*L;
//															}
//															else if(dataType=1)
//																pre=preData[iO].data+iT*L;
											pC[iO]+=temp1*pre[iO][l];
											pD[iO]+=temp2*pre[iO][templ1];
											pPC[l].pdata[iO]+=temp1*pre[iO][l];
											pPC[templ1].pdata[iO]+=temp2*pre[iO][templ1];
										}
										else if(iO==order)
										{

											pC[iO]+=temp1;
											pD[iO]+=temp2;
										}
																	} // for iO
							
							
							} // for l
			
			} // for m
					
	  
		iSeg++;
		} // end for iT

		if(ite==0)
		{
			for(l=0;l<L;l++)
				vFloor[l] *= vfloor/iSeg;
		}

		for(m=0;m<M;m++)
		{
			for(l=0;l<L/2;l++)
			{		
				for(iO=0;iO<order+L/2+1;iO++)
					for(jO=iO+1;jO<order+L/2+1;jO++)
					{
						if(iO>order)
							B.mix[m].dim[l].pdata[iO*(B.order)+jO]=B.mix[m].dim[l].pdata[jO*(B.order)+iO];
						else if(iO==order)
							B.mix[m].dim[l].pdata[iO*(B.order)+jO]=B.mix[m].dim[l].pdata[jO*(B.order)+iO];
						else if(iO<order && jO<order)
						{
							A.mix[m].dim[l].pdata[iO*(A.order)+jO]=A.mix[m].dim[l].pdata[jO*(A.order)+iO];
							B.mix[m].dim[l].pdata[iO*(B.order)+jO]=B.mix[m].dim[l].pdata[jO*(B.order)+iO];
						}
						else if(iO<order && jO==order)
						{
							A.mix[m].dim[l].pdata[iO*(A.order)+jO]=A.mix[m].dim[l].pdata[jO*(A.order)+iO];
							B.mix[m].dim[l].pdata[iO*(B.order)+jO]=B.mix[m].dim[l].pdata[jO*(B.order)+iO];
						}
						else if(iO<order && jO>order)
						{
							B.mix[m].dim[l].pdata[iO*(B.order)+jO]=B.mix[m].dim[l].pdata[jO*(B.order)+iO];
						}
					}

				}
		}
		printf("%dth iteration,liklihood:%f/%d=%f\n",ite+1,likliProb,iSeg,likliProb/iSeg);

		if(ite>=1 && ABS(likliProb/iSeg-preLikli)< thr)
		{
			
				printf("收敛\n");
				break;
		}
		
		preLikli=likliProb/iSeg;
		//printf("update parameters...\n");
		
		
		for(m=0;m<M;m++)
		{
			gmm.gauss[m].weight=(float)(accW[m]/(float)iSeg);
			for(l=0;l<L/2;l++)
			{

//				printf("A:\n");
//				for(iO=0;iO<A.order;iO++)
//				{
//					for(jO=0;jO<A.order;jO++)
//					{
//						printf("%f ",A.mix[m].dim[l].pdata[iO*A.order+jO]);
//					}
//					printf("\n");
//				}
//				printf("C:\n");
//				for(iO=0;iO<A.order;iO++)
//				{
//						printf("%f ",C.mix[m].dim[l].pdata[iO]);
//				}
//				printf("\n");

				int info=solve(A.mix[m].dim[l].pdata,C.mix[m].dim[l].pdata,A.order);
				if(info)
				  printf("can't not solve at %dmix %ddim\n",m,l);				
				else
				{
				
//								printf("C:\n");
//								for(iO=0;iO<A.order;iO++)
//								{
//										printf("%f ",C.mix[m].dim[l].pdata[iO]);
//								}
//								printf("\n");
					
						for(iO=0;iO<order;iO++)
						gmm.gauss[m].ar.data[iO*L+l]=C.mix[m].dim[l].pdata[iO];
					gmm.gauss[m].miu.data[l]=C.mix[m].dim[l].pdata[order];
				}


//				printf("B:\n");
//				for(iO=0;iO<B.order;iO++)
//				{
//					for(jO=0;jO<B.order;jO++)
//					{
//						printf("%f ",B.mix[m].dim[l].pdata[iO*B.order+jO]);
//					}
//					printf("\n");
//				}
//
//				printf("D:\n");
//				for(iO=0;iO<B.order;iO++)
//				{
//						printf("%f ",D.mix[m].dim[l].pdata[iO]);
//				}
//				printf("\n");
				info=solve(B.mix[m].dim[l].pdata,D.mix[m].dim[l].pdata,B.order);
				if(info)
				  printf("can't not solve at %dmix %ddim\n",m,l);				
				else
				{
				
//				printf("D:\n");
//				for(iO=0;iO<B.order;iO++)
//				{
//						printf("%f ",D.mix[m].dim[l].pdata[iO]);
//				}
//				printf("\n");
				
				}
				for(iO=0;iO<order;iO++)
					gmm.gauss[m].ar.data[iO*L+l+L/2]=D.mix[m].dim[l].pdata[iO];
				gmm.gauss[m].miu.data[l+L/2]=D.mix[m].dim[l].pdata[order];
				for(iL=0;iL<L/2;iL++)
					gmm.gauss[m].bt.data[L/2*l+iL]=D.mix[m].dim[l].pdata[order+1+iL];
				
			 gmm.gauss[m].var.data[l]=accZZ.mix[m].dim[l].data - accZ.mix[m].dim[l].data * gmm.gauss[m].miu.data[l];
			 gmm.gauss[m].var.data[l+L/2]=accZZ.mix[m].dim[l+L/2].data - accZ.mix[m].dim[l+L/2].data * gmm.gauss[m].miu.data[l+L/2];
			 for(iO=0;iO<order;iO++)
			 {
			 	gmm.gauss[m].var.data[l]-=accPC.mix[m].dim[l].pdata[iO]*gmm.gauss[m].ar.data[iO*L+l];
			 	gmm.gauss[m].var.data[l+L/2]-=accPC.mix[m].dim[l+L/2].pdata[iO]*gmm.gauss[m].ar.data[iO*L+l+L/2];
			 }

			//void matrixMultiply(float* A,int ra,int ca,float * B,int cb,float *C)
			//matrixMultiply(gmm.gauss[m].bt.data,L/2,L/2,accXY.mix[m].dim[0].pdata,L/2,tempC);
			for(iL=0;iL<L/2;iL++)
					gmm.gauss[m].var.data[l+L/2]-=gmm.gauss[m].bt.data[L/2*l+iL]*accXY.mix[m].dim[0].pdata[l*L/2+iL];
			gmm.gauss[m].var.data[l] /= accW[m];
			gmm.gauss[m].var.data[l+L/2] /= accW[m];
			
			}//end l
			det=1.0;

			for(l=0;l<L;l++)
			{
				det *= gmm.gauss[m].var.data[l];
			}
			if(det<=0)
			{
					printf("add vfloor at %dmix\n",m);
					for(l=0;l<L;l++)
						gmm.gauss[m].var.data[l]+=vFloor[l];
			}
			cal_gconst(gmm.gauss+m);
			
		}// end m
		//next iteration
			ite+=1;	
		zfillACC(&A);
		zfillACC(&B);
		zfillACC(&C);
		zfillACC(&D);
		zfillACC(&accZZ);
		zfillACC(&accZ);
		zfillACC(&accPC);
		zfillACC(&accXY);
		zfilld(accW,M);
	} // while
		FILE *fogmm=fopen(*++argv,"wb");
		fwriteGMM(&gmm,fogmm);
		fclose(fogmm);
		//fwrite(accXY.mix[0].dim[0].pdata,sizeof(double),L/2*L/2,fogmm);
  free(postProb);
	free(accW);
	freeACC(&A);
	freeACC(&B);
	freeACC(&C);
	freeACC(&D);
	freeGMM(&gmm);
	free(vFloor);
	freeM(&trainData);
	for(iO=0;iO<order;iO++)
	{
		freeM(preData+iO);
	}
	free(pre);
	//free(tempx);
	//free(tempy);
	return 0;
}












