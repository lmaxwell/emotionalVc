/*
 * =====================================================================================
 *
 *       Filename:  ar-gmm_convert.c
 *
 *    Description:  
 *    						-m mixture number
 *    						-l dimension
 *    						-o order
 *    						-M model file
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
}JGAUSS;

typedef struct _GMM{
	int M;
	int L;
	int order;
	JGAUSS *gauss;
	
} GMM;


void zfill(float *a,int L)
{
	memset(a,0,sizeof(float)*L);
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
	double det=1.0;
	int i;
	for (i=0;i<L/2;i++)
	{
	
		det*=(double)gauss->var.data[i];
	}
	//printf("det:%f\n",det);
	gauss->gconst=log(gauss->weight)-0.5*log(det)-L/4*log(2.0*PI);
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
		printf("gconst[%d]:%f\n",m,gmm->gauss[m].gconst);
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

void convert(float *src,float *conv,FMATRIX *pre,GMM *gmm)
{
	int i,jL,m,iO,L=gmm->L,M=gmm->M,order=gmm->order,maxMix=0;
	double maxProb=-1.0e10;
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


int main(int argc,char **argv)
{
	int M=64;
	int order=3;
	int L=48;
	FILE *fgmm;
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

						}
	}
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
	//	gmm.gauss[m].preMiu.row=order;
		//gmm.gauss[m].preMiu.col=L;
		getmemV(&(gmm.gauss[m].miu));
		getmemV(&(gmm.gauss[m].var));
		getmemM(&(gmm.gauss[m].ar));
		getmemM(&(gmm.gauss[m].bt));
	//	getmemM(gmm.gauss[m].preMiu);
	}

//	FILE *fgmm=fopen("neutral-sad-order3.ar-gmm","rb");
	readModel(fgmm,&gmm);
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
	
	for (iT=0;iT<T;iT++)
	{
		//printf("iT=%d\n",iT);
	  cmcep.data[iT][0]=omcep.data[iT][0];
		convert(*(omcep.data+iT)+1,*(cmcep.data+iT)+1,&pre,&gmm);
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
	FILE *fcmcep=fopen(*++argv,"wb");
	for (iT=0;iT<T;iT++)
	{
		fwrite(*(cmcep.data+iT),sizeof(float),25,fcmcep);
	}
	//convert()
	return 0;
}

