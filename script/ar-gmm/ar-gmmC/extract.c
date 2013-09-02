/*
 * =====================================================================================
 *
 *       Filename:  extract.c
 *
 *    Description:  extract -l dim -o order -t 阈值 findex fsrc ftrg fpowsrc fpowtrg fout
 *
 *        Version:  1.0
 *        Created:  2013年08月02日 16时38分01秒
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

int main(int argc,char **argv)
{
	FILE *findex,*fsrc,*ftrg,*fpowsrc,*fpowtrg,*fjoin;
	char outFile[10][400];	
	int l,order,length,i,j,io,lenPowsrc,lenPowtrg,srcl,trgl;
	int *index[2];
	float thr,*joinData,*powsrc,*powtrg,*src,*trg,*beginZero;
	while(--argc > 6)
		if(**++argv=='-')
			switch(*(*argv+1))
			{
				case 'l':
					l=atoi(*++argv);
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
			}

   findex=fopen(*++argv,"rb");
   fsrc=fopen(*++argv,"rb");
   ftrg=fopen(*++argv,"rb");
	 fpowsrc=fopen(*++argv,"rb");
	 fpowtrg=fopen(*++argv,"rb");
   ++argv;
   fseek(findex,0L,SEEK_END);
   length=ftell(findex)/sizeof(float)/2;
   fseek(findex,0L,SEEK_SET);

   fseek(fsrc,0L,SEEK_END);
	 srcl=ftell(fsrc)/sizeof(float);
	 fseek(fsrc,0L,SEEK_SET);
   
	 fseek(ftrg,0L,SEEK_END);
	 trgl=ftell(ftrg)/sizeof(float);
	 fseek(ftrg,0L,SEEK_SET);
   
	 src=(float *)malloc(sizeof(float)*srcl);
	 trg=(float *)malloc(sizeof(float)*trgl);
	 fread(src,sizeof(float),srcl,fsrc);
	 fread(trg,sizeof(float),trgl,ftrg);
	
	//printf("size:%d,length:%d\n",ftell(findex),length);
	for (i=0;i<2;i++)
	{
		index[i]=(int *)malloc(sizeof(int)*length);
	}
	for (j=0;j<length;j++)
	{
		fread(&index[0][j],sizeof(int),1,findex);
		fread(&index[1][j],sizeof(int),1,findex);
	//	printf("%d \n",index[0][j]);
	//	printf("%d \n",index[1][j]);
	}
	printf("here\n");
	fseek(fpowsrc,0L,SEEK_END);
	lenPowsrc=ftell(fpowsrc)/sizeof(float);
	fseek(fpowsrc,0L,SEEK_SET);
	fseek(fpowtrg,0L,SEEK_END);
	lenPowsrc=ftell(fpowtrg)/sizeof(float);
	fseek(fpowtrg,0L,SEEK_SET);
	powsrc=(float *)malloc(sizeof(float)*lenPowsrc);
	powtrg=(float *)malloc(sizeof(float)*lenPowtrg);
	for(i=0;i<lenPowsrc;i++)
	{
		fread(&powsrc[i],sizeof(float),1,fpowsrc);
	}
	
	for(i=0;i<lenPowtrg;i++)
	{
		fread(&powtrg[i],sizeof(float),1,fpowtrg);
	}

	fclose(fpowsrc);
	fclose(fpowtrg);
   char p_order[10];

		beginZero=(float *)calloc(l,sizeof(float));
   for (io=0;io<=order;io++)
   {
   		strcpy(&outFile[io][0],*argv);
		sprintf(p_order,"-%d",io);
		strcat(&outFile[io][0],p_order);
   

   fjoin=fopen(outFile[io],"wb");
   
   	//joinData=(float *)malloc(sizeof(float)*l*2*length);
	for (i=0,j=0;i<length;i++)
	{


		if(powsrc[index[0][i]]<thr || powtrg[index[1][i]]<thr) 
		{
			continue;
		}			
		if(index[0][i]-io<0)
			//memset(&joinData[j*l*2],0,sizeof(float)*l);
			fwrite(beginZero,sizeof(float),l,fjoin);
		else
		{
	   	//	fseek(fsrc,sizeof(float)*(index[0][i]-io)*l,SEEK_SET);
			//fread(&joinData[j*l*2],sizeof(float),l,fsrc);
			//memcpy(&joinData[j*l*2],src+(index[0][i]-io)*l,sizeof(float)*l);
			fwrite(src+(index[0][i]-io)*l,sizeof(float),l,fjoin);
		}
		if(index[1][i]-io<0)
			//memset(&joinData[j*l*2+l],0,sizeof(float)*l);
			fwrite(beginZero,sizeof(float),l,fjoin);
		else
		{
			//fseek(ftrg,sizeof(float)*(index[1][i]-io)*l,SEEK_SET);
			//fread(&joinData[j*l*2+l],sizeof(float),l,ftrg);
			//memcpy(&joinData[j*l*2+l],trg+(index[1][i]-io)*l,sizeof(float)*l);
			fwrite(trg+(index[1][i]-io)*l,sizeof(float),l,fjoin);
		}
		j++;
		//exit(0);
	 }                   
	//fwrite(joinData,sizeof(float),l*2*j,fjoin);
	//free(joinData);
	fclose(fjoin);
   }
	//for (j=0;j<l*2*length;j++)
	//	printf("%f\n",joinData[j]);
	//
// 	free(index);
  free(powsrc);
	free(powtrg);
	free(src);
	free(trg);
	free(beginZero);
	fclose(findex);
	fclose(fsrc);
	fclose(ftrg);
	for (i=0;i<2;i++)
	{
		free(index[i]);
	}

	return 0;
}

