/*
 * =====================================================================================
 *
 *       Filename:  ar-gmm_convert.c
 *
 *    Description:  
 *    						-m mixture number
 *    						-l joint dimension length
 *    						-o order
 *    						-M model file
 *    						-g global variance model file
 *    						-G 0:no global variance 1: with global variance
 *    						-s step size
 *    						-i max iteration number
 *    						input  (L/2+1)*T
 *    						output (L/2+1)*T
 *
 *        Version:  1.0
 *        Created:  2013年08月21日 10时03分17秒
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
#include <cblas.h>
#include <math.h>
#include "ar-gmm_sub.h"

#define ABS(x) ((x)>0 ? (x):(-(x)))
#define PI 3.1415926

typedef struct _FVECTOR
{
	int len;
	float *data;
}FVECTOR;

typedef struct _FMATRIX
{
	int row;
	int col;
	float **data;
}FMATRIX;


typedef struct _JGAUSS{
	float weight;
	FVECTOR miu;
	FVECTOR var;
	FMATRIX ar;
	FMATRIX bt;
//	FMATRIX preMiu;
	float gconst;
	float gconsty_;
	FVECTOR gconstx;
	FVECTOR gconsty;
}JGAUSS;

typedef struct _GMM{
	int M;
	int L;
	int order;
	JGAUSS *gauss;
	
} GMM;

typedef struct _GV{
	FVECTOR miu;
	FVECTOR var;
	double *gconst;
	double gconst_;
}GV;
void zfill(float *a,int L)
{
	memset(a,0,sizeof(float)*L);
}
void zfilld(double *a,int L)
{
	memset(a,0,sizeof(double)*L);
}
void zfillV(FVECTOR *v)
{
	int L=v->len;
	zfill(v->data,L);
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
	int i;
	matrix->data=(float **)calloc(m,sizeof(float *));
	for(i=0;i<m;i++)
	{
		matrix->data[i]=(float *)calloc(n,sizeof(float));
	}
}

void freeM(FMATRIX *matrix)
{
	int m=matrix->row;
	int n=matrix->col;
	int i;
	for(i=0;i<m;i++)
		free(matrix->data[i]);
	free(matrix->data);
	
}

void cal_gconst(JGAUSS *gauss)
{
	int L=gauss->miu.len;
	double det=1.0,dety=1.0;
	int i;
	for (i=0;i<L/2;i++)
	{
		det*=(double)gauss->var.data[i];
		dety*=(double)gauss->var.data[i+L/2];
	}
	//printf("det:%f\n",det);
	gauss->gconst=log(gauss->weight)-0.5*log(det)-L/4*log(2.0*PI);
	gauss->gconsty_=-0.5*log(dety)-L/4*log(2.0*PI);
}
void cal_gconstx(JGAUSS *gauss)
{
	int L=gauss->miu.len;
	double det=1.0;
	int i;
	for (i=0;i<L/2;i++)
	{
	
		det=(double)gauss->var.data[i];
		gauss->gconstx.data[i]=-0.5*log(det)-0.5*log(2*PI);
	}
	//printf("det:%f\n",det);
//	gauss->gconsty=-0.5*log(det)-L/4*log(2.0*PI);

}

void cal_gconsty(JGAUSS *gauss)
{
	int l=gauss->miu.len;
	double det=1.0;
	int i;
	for (i=l/2;i<l;i++)
	{
	
		det=(double)gauss->var.data[i];
		gauss->gconsty.data[i-l/2]=-0.5*log(det)-0.5*log(2*PI);
	}
	//printf("det:%f\n",det);
//	gauss->gconsty=-0.5*log(det)-l/4*log(2.0*pi);

}

void cal_gv_gconst(GV *gv)
{
	int L=gv->miu.len;
	double det=1.0,det_=1.0;
	int i;
	for (i=0;i<L;i++)
	{
		det=(double)gv->var.data[i];
	gv->gconst[i] = -0.5*log(det)- 0.5*log(2.0*PI);
	//printf("gv->gconst[%d]:%f\n",i,gv->gconst[i]);
		det_*=(double)gv->var.data[i];
	}
	gv->gconst_=-0.5*log(det_)-0.5*L*log(2.0*PI);
	printf("gv->gconst_:%f\n",gv->gconst_);
}


void readModel(FILE *fgmm,GMM *gmm)
{
	int m,l,iO,pos;
	int L=gmm->L,M=gmm->M,order=gmm->order;
	float *p,**pp;
	fseek(fgmm,0L,SEEK_END);
	int length=ftell(fgmm)/sizeof(float);
	float *data=(float *)malloc(length*sizeof(float));
	fseek(fgmm,0L,SEEK_SET);
	fread(data,sizeof(float),length,fgmm);
	fclose(fgmm);
	//printf("length:%d\n",length);
	for(m=0;m < M;m++)
	{
		gmm->gauss[m].weight=data[m];
	}
	
	pos=M;
	//	printf("order=%d\n",order);
	for(m=0;m< M;m++)
	{
		//pos=M+m*(2*L+order*L+L/2*L/2);
		memcpy(gmm->gauss[m].miu.data,data+pos,sizeof(float)*(L));
		pos += L;

		for(iO=0;iO<order;iO++)
		{
	//	printf("iO=%d,m=%d\n",iO,m);
		  pp=gmm->gauss[m].ar.data;
			memcpy(*(pp+iO),data+pos,sizeof(float)*L);
			pos+=L;
		}
//		for(iO=0;iO<order;iO++)
//		{
//		  
//		  pp=gmm->gauss[m].ar.data;
//			memcpy(*(pp+iO)+L/2,data+pos,sizeof(float)*L/2);
//			pos+=L/2;
//	//	printf("iO=%d,m=%d,pos=%d\n",iO,m,pos);
//		}
		
		memcpy(gmm->gauss[m].var.data,data+pos,sizeof(float)*L);
		pos+=L;

		for(l=0;l<L/2;l++)
		{
			pp=gmm->gauss[m].bt.data;				
			memcpy(*(pp+l),data+pos,sizeof(float)*L/2);
			pos+=L/2;
		}
		
	
		cal_gconst(gmm->gauss+m);
		cal_gconstx(gmm->gauss+m);
		cal_gconsty(gmm->gauss+m);
		//printf("gconst[%d]:%f\n",m,gmm->gauss[m].gconst);
	}

}

double cal_postprob(float *miu,float *x,JGAUSS *gauss)
{
	int L=gauss->miu.len;
	int i;
	double postprob=gauss->gconst;
	for(i=0;i<L/2;i++)
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



void convert(float *src,float *conv,FMATRIX *pre,GMM *gmm,int *pMix,double *postProb)
{
	int i,jL,m,iO,L=gmm->L,M=gmm->M,order=gmm->order,maxMix=0;
	double maxProb=-1.0e10,sumProb=-1.0E10;
	double *postprob=(double *)calloc(M,sizeof(double));
	FMATRIX curMiu;
	curMiu.row=M;
	curMiu.col=L/2;
	getmemM(&curMiu);

	for(m=0;m<M;m++)
	{
		for(i=0;i<L/2;i++)
		{
			curMiu.data[m][i]=gmm->gauss[m].miu.data[i];
			for (iO=0;iO<order;iO++)
			{
				//printf("m=%d,i=%d,iO=%d\n",m,i,iO);
				curMiu.data[m][i]+=gmm->gauss[m].ar.data[iO][i] * pre->data[iO][i];
			}
		}
		postprob[m]=cal_postprob(*(curMiu.data+m),src,gmm->gauss+m);
		if(postprob[m]>maxProb)
		{
			maxMix=m;
			maxProb=postprob[m];
		}
		sumProb=log_add(sumProb,postprob[m]);
	}
	
	for(i=0;i<L/2;i++)
	{
		conv[i]=gmm->gauss[maxMix].miu.data[L/2+i];
			for (iO=0;iO<order;iO++)
			{
				conv[i]+=gmm->gauss[maxMix].ar.data[iO][i+L/2]*pre->data[iO][i+L/2];
			}
		for(jL=0;jL<L/2;jL++)
			conv[i]+=gmm->gauss[maxMix].bt.data[i][jL]*src[jL];


		//*(postProb+i) = gmm->gauss[maxMix].gconstx.data[i]-0.5*(src[i]-gmm->gauss[maxMix].miu.data[i])*(src[i]-gmm->gauss[maxMix].miu.data[i])/gmm->gauss[maxMix].var.data[i];
	}
	*postProb=postprob[maxMix]-sumProb;
	*pMix=maxMix;
}

void readGV(FILE *fgv,GV *gv)
{
	int l=gv->miu.len;			
	fseek(fgv,l*sizeof(float),SEEK_SET);
	fread(gv->miu.data,sizeof(float),l,fgv);
	fseek(fgv,l*sizeof(float),SEEK_CUR);
	fread(gv->var.data,sizeof(float),l,fgv);
	int i;
	for(i=0;i<l;i++)
	{
		printf("gv.miu[%d]=%f\n",i,gv->miu.data[i]);
	}

	for(i=0;i<l;i++)
	{
		printf("gv.var[%d]=%f\n",i,gv->var.data[i]);
	}
}
/* float **getmemM(int m,int n)
 * {
 * 	int i;
 * 	float **p=NULL;
 * 	p=(float **)calloc(m,sizeof(float *));
 * 	for(i=0;i<m;i++)
 * 	{
 * 		p[m]=(float *)calloc(n,sizeof(float));
 * 	}
 * 	return p;
 * }
 */


//double cal_likli(int T,int l,int order,FMATRIX *omcep,FMATRIX *cmcep,JGAUSS *gauss,int *maxMix,FVECTOR *curMiu,GV *senSta,GV *gv,double *postProb,double w1,double w2)
//{
//	zfillV(&(senSta->miu));
//	zfillV(&(senSta->var));
//
//	int iT,iL,jL,L=gauss->miu.len,iO;
//	double likli=0.0,vLikli=0.0;
//	for(iT=0;iT<T;iT++)
//	{
//		likli+=gmm.gauss[maxMix[iT]].gconsty_;
//		for(iL=0;iL<L/2;iL++)
//		{			
//			//likli+=gmm.gauss[maxMix[iT]].gconsty.data[l];
//			senSta.miu.data[iL]+=cmcep.data[iT][l+1];
//			senSta.var.data[iL]+=cmcep.data[iT][l+1]*cmcep.data[iT][iL+1];
//			
//			curMiu.data[iL]=gmm.gauss[maxMix[iT]].miu.data[iL+L/2];
//			for (iO=0;iO<order;iO++)
//			{
//				//printf("m=%d,i=%d,iO=%d\n",m,i,iO);
//				if(iT-iO>=0)
//					curMiu.data[iL]+=gmm.gauss[maxMix[iT]].ar.data[iO][iL+L/2] * cmcep.data[iT-iO][iL+1];
//			}
//			for(jL=0;jL<L/2;jL++)
//				curMiu.data[iL]+=gmm.gauss[maxMix[iT]].bt.data[iL][jL] * omcep.data[iT][jL+1];
//			likli += -0.5*(cmcep.data[iT][l+1]-curMiu.data[iL])*(cmcep.data[iT][l+1]-curMiu.data[iL])/gmm.gauss[maxMix[iT]].var.data[l+L/2];
//			
//		}	
//		likli += postProb[iT];
//	}
//	for(iL=0;iL<L/2;iL++)
//	{
//			senSta.miu.data[l]/=(float)T;
//			senSta.var.data[l]/=(float)T;
//			senSta.var.data[l]-=senSta.miu.data[l]*senSta.miu.data[l];
//			
//			vLikli+= -0.5*( senSta.var.data[l] - gv.miu.data[l] )*( senSta.var.data[l] - gv.miu.data[l] ) / gv.var.data[l];
//		//	printf("vLikli:%f\n",vLikli);
//	}
//	//vLikli += gv.gconst[l];
//	vLikli += gv.gconst_;
//	likli=w1*likli/((double)T) + w2*vLikli;
//	return likli;
//}

int main(int argc,char **argv)
{
	int M=64;
	int order=3;
	int L=48;
	int l;
	int ngvite=20,withgv;
	float step=0.01;
	FILE *fgmm,*fgv;
	while(--argc>2)
	{
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
							case 'M':
											fgmm=fopen(*++argv,"rb");
											--argc;
											break;
							case 'g':
											fgv=fopen(*++argv,"rb");
											--argc;
											break;
						 case 'i':
											ngvite=atoi(*++argv);
											--argc;
											break	;
					 	 case 'G':
											withgv=atoi(*++argv) ;
											--argc;
											break;

						case 's':
											step=atof(*++argv);
											--argc;
											break;

						}
	}
	ngvite*=withgv;
	GMM gmm;
	gmm.gauss=(JGAUSS *)malloc(sizeof(JGAUSS)*M);
	int m;
	gmm.M=M;
	gmm.L=L;
	gmm.order=order;
	for (m=0;m<M;m++)
	{
		gmm.gauss[m].miu.len=L;
		gmm.gauss[m].var.len=L;
		gmm.gauss[m].ar.row=order;
		gmm.gauss[m].ar.col=L;
		gmm.gauss[m].bt.row=L/2;
		gmm.gauss[m].bt.col=L/2;
		gmm.gauss[m].gconsty.len=L/2;
		gmm.gauss[m].gconstx.len=L/2;
	//	gmm.gauss[m].preMiu.row=order;
		//gmm.gauss[m].preMiu.col=L;
		getmemV(&(gmm.gauss[m].miu));
		getmemV(&(gmm.gauss[m].var));
		getmemM(&(gmm.gauss[m].ar));
		getmemM(&(gmm.gauss[m].bt));
		getmemV(&(gmm.gauss[m].gconsty));
		getmemV(&(gmm.gauss[m].gconstx));
	//	getmemM(gmm.gauss[m].preMiu);
	}

//	FILE *fgmm=fopen("neutral-sad-order3.ar-gmm","rb");
	readModel(fgmm,&gmm);

	GV gv;
	gv.miu.len=L/2;
	gv.var.len=L/2;
	getmemV(&(gv.miu));
	getmemV(&(gv.var));
	gv.gconst=(double *)calloc(L/2,sizeof(double));
	readGV(fgv,&gv);
	fclose(fgv);
  
	cal_gv_gconst(&gv);


	FILE *fomcep=fopen(*++argv,"rb");
	fseek(fomcep,0,SEEK_END);
	int T=ftell(fomcep)/sizeof(float)/25;
	fseek(fomcep,0,SEEK_SET);
	FMATRIX omcep;
	omcep.row=T;
  omcep.col=25;
	getmemM(&omcep);
	int iT,iO;
	for(iT=0;iT<T;iT++)
	{				
		fread(*(omcep.data+iT),sizeof(float),25,fomcep);
	}
	FMATRIX cmcep;
	cmcep.row=T;
  cmcep.col=25;
	getmemM(&cmcep);	
	FMATRIX pre;
	pre.row=order;
  pre.col=48;
	getmemM(&pre);	
	int *maxMix=(int *)malloc(sizeof(int)*T);
	double *postProb=(double *)malloc(sizeof(double)*T);
	for (iT=0;iT<T;iT++)
	{
		//printf("iT=%d\n",iT);
	  cmcep.data[iT][0]=omcep.data[iT][0];
		convert(*(omcep.data+iT)+1,*(cmcep.data+iT)+1,&pre,&gmm,maxMix+iT,postProb+iT);
		//printf("pre%f pre%f\n",pre.data[0][0],pre.data[1][0]);
		for(iO=0;iO<order-1;iO++)
		{
			memcpy(*(pre.data+iO+1),*(pre.data+iO),48*sizeof(float));
			//printf("pre%f pre%f\n",pre.data[0][0],pre.data[1][0]);
		}
		if(order>0)
		{
		memcpy(*(pre.data),*(omcep.data+iT)+1,24*sizeof(float));
		memcpy(*(pre.data)+24,*(cmcep.data+iT)+1,24*sizeof(float));
		}
		//printf("omcep%f cmcep%f\n",omcep.data[iT][1],cmcep.data[iT][1]);
		//printf("pre%f pre%f %f %f \n",pre.data[0][0],pre.data[1][0]);
	}
	
	//calculate global variance
	GV senSta;
	senSta.miu.len=L/2;
	senSta.var.len=L/2;
	getmemV(&(senSta.miu));
	getmemV(&(senSta.var));
	FVECTOR curMiu;
	curMiu.len=L/2;
	getmemV(&curMiu);

	int gvite,iL,jL;
	float temp;
	//int maxMix=0;
	double maxProb=-1.0e10;
	double *postprob=(double *)calloc(M,sizeof(double));
	double *preDelta=(double *)calloc(T*L/2,sizeof(double));
	double *curDelta=(double *)calloc(T*L/2,sizeof(double));
  double likli,pLikli=0.0,vLikli;
	
if(ngvite>0)
{
//	double cal_likli(int T,JGAUSS *gauss,FVECTOR *curMiu,GV *senSta)
	for(iT=0;iT<T;iT++)
	{
		likli+=gmm.gauss[maxMix[iT]].gconsty_;
		for(l=0;l<L/2;l++)
		{
			senSta.miu.data[l]+=cmcep.data[iT][l+1];
			senSta.var.data[l]+=cmcep.data[iT][l+1]*cmcep.data[iT][l+1];
			
			curMiu.data[l]=gmm.gauss[maxMix[iT]].miu.data[l+L/2];
			for (iO=0;iO<order;iO++)
			{
				//printf("m=%d,i=%d,iO=%d\n",m,i,iO);
				if(iT-iO>=0)
					curMiu.data[l]+=gmm.gauss[maxMix[iT]].ar.data[iO][l+L/2] * cmcep.data[iT-iO][l+1];
			}
			for(jL=0;jL<L/2;jL++)
				curMiu.data[l]+=gmm.gauss[maxMix[iT]].bt.data[l][jL] * omcep.data[iT][jL+1];
			likli += -0.5*(cmcep.data[iT][l+1]-curMiu.data[l])*(cmcep.data[iT][l+1]-curMiu.data[l])/gmm.gauss[maxMix[iT]].var.data[l+L/2];
			
		}
	}
	for(l=0;l<L/2;l++)
	{
			senSta.miu.data[l]/=(float)T;
			senSta.var.data[l]/=(float)T;
			senSta.var.data[l]-=senSta.miu.data[l]*senSta.miu.data[l];
			
			vLikli+= -0.5*( senSta.var.data[l] - gv.miu.data[l] )*( senSta.var.data[l] - gv.miu.data[l] ) / gv.var.data[l];
		//	printf("vLikli:%f\n",vLikli);
	}
	vLikli += gv.gconst_;
	likli=likli/((double)T) + vLikli;
	

	for(iT=0;iT<T;iT++)
		for(l=0;l<L/2;l++)
		{
			cmcep.data[iT][l+1]=	sqrt(gv.miu.data[l]/senSta.var.data[l])*(cmcep.data[iT][l+1]-senSta.miu.data[l]) +senSta.miu.data[l];
		}
}

double w1=1.0;///(double)order;
double w2=1.0;
//for(l=0;l<L/2;l++)
//{
	//likli=cal_likli(T,l,order,&omcep,&cmcep,gmm.gauss,maxMix,&curMiu,&senSta,&gv,postProb,w1,w2);
  //pLikli=likli;
	for(gvite=0;gvite<ngvite;gvite++)
	{
					likli=0.0;
					vLikli=0.0;
					zfillV(&(senSta.miu));
					zfillV(&(senSta.var));
					zfilld(curDelta,T*L/2);
					for(iT=0;iT<T;iT++)
					{
						likli+=gmm.gauss[maxMix[iT]].gconsty_;
						for(iL=0;iL<L/2;iL++)
						{			
							//likli+=gmm.gauss[maxMix[iT]].gconsty.data[l];
							senSta.miu.data[iL]+=cmcep.data[iT][iL+1];
							senSta.var.data[iL]+=cmcep.data[iT][iL+1]*cmcep.data[iT][iL+1];
							
							curMiu.data[iL]=gmm.gauss[maxMix[iT]].miu.data[iL+L/2];
							for (iO=0;iO<order;iO++)
							{
								//printf("m=%d,i=%d,iO=%d\n",m,i,iO);
								if(iT-iO>=0)
									curMiu.data[iL]+=gmm.gauss[maxMix[iT]].ar.data[iO][iL+L/2] * cmcep.data[iT-iO][iL+1];
							}
							for(jL=0;jL<L/2;jL++)
								curMiu.data[iL]+=gmm.gauss[maxMix[iT]].bt.data[iL][jL] * omcep.data[iT][jL+1];
							likli += -0.5*(cmcep.data[iT][iL+1]-curMiu.data[iL])*(cmcep.data[iT][iL+1]-curMiu.data[iL])/gmm.gauss[maxMix[iT]].var.data[iL+L/2];
							
						}	
						likli += postProb[iT];
					}
					for(iL=0;iL<L/2;iL++)
					{
							senSta.miu.data[iL]/=(float)T;
							senSta.var.data[iL]/=(float)T;
							senSta.var.data[iL]-=senSta.miu.data[iL]*senSta.miu.data[iL];
							
							vLikli+= -0.5*( senSta.var.data[iL] - gv.miu.data[iL] )*( senSta.var.data[iL] - gv.miu.data[iL] ) / gv.var.data[iL];
						//	printf("vLikli:%f\n",vLikli);
					}
					//vLikli += gv.gconst[l];
					vLikli += gv.gconst_;
					likli=w1*likli/((double)T) + w2*vLikli;
					
					if(gvite>=1)
					{				
									if(ABS(likli-pLikli)<1E-4)
									{
													printf("%dth iteration 收敛,likli:%f,pLikli:%f,%f<%f\n",gvite,likli,pLikli,ABS(likli-pLikli),1E-4);
													break;
									}

									if(likli>pLikli)
													step*=1.2;
									else
									{
										for(iT=0;iT<T;iT++)
										{

											for(l=0;l<L/2;l++)
											{
												cmcep.data[iT][l+1] -= step*preDelta[iT*L/2+l];
											}
										}

										step*=0.5;

										for(iT=0;iT<T;iT++)
										{

											for(l=0;l<L/2;l++)
											{
												cmcep.data[iT][l+1] += step*preDelta[iT*L/2+l];
											}
										}

													printf("%dth iteration  go back likli:%f vLikli:%f likli+vLikli:%f \n",gvite,likli-vLikli,vLikli,likli);
													continue;
									}
					}
					pLikli=likli;
					printf("%dth iteration likli:%f vLikli:%f likli+vLikli:%f \n",gvite,likli-vLikli,vLikli,likli);

					for(iT=0;iT<T;iT++)
					{

						for(l=0;l<L/2;l++)
						{
							temp=cmcep.data[iT][l+1];
							curDelta[iT*L/2+l]=(-1*w1/((double)T)/gmm.gauss[maxMix[iT]].var.data[l+L/2]*(temp-curMiu.data[l]) + w2*(-2)/(double)T / gv.var.data[l] *( senSta.var.data[l] - gv.miu.data[l] ) * (temp-senSta.miu.data[l])  );
							cmcep.data[iT][l+1] += step* curDelta[iT*L/2+l];
							preDelta[iT*L/2+l] = curDelta[iT*L/2+l];
						}
//							postprob[m]=gmm.gauss[m].gconst;
//							for(l=0;l<L/2;l++)
//							{
//								postprob[m]+=-0.5*(cmcep.data[iT][l+1]-curMiu.data[m][l])*(cmcep.data[iT][l+1]-curMiu.data[m][l])/gmm.gauss[m].var.data[l+L/2];
//							}
//
//						if(postprob[m]>maxProb)
//						{
//							maxMix=m;
//							maxProb=postprob[m];
//						}
						

											
				}
				

	}

//} //end for l
//	printf("%s\n",*(argv+1));
	FILE *fcmcep=fopen(*++argv,"wb");
	for (iT=0;iT<T;iT++)
	{
		fwrite(*(cmcep.data+iT),sizeof(float),25,fcmcep);
	}
	free(postprob);
	//convert()
	return 0;
}

