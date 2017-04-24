#include<stdio.h>
#include<stdlib.h>

int main()
{
	FILE *fp,*rst;
	if((fp=fopen("S2LFFS15%.txt","r"))==NULL)
	{
		printf("cannot open file\n");
		exit(1);
	}
	if((rst=fopen("s2lffs15%.dat","w"))==NULL)
	{
		printf("cannot open file\n");
		exit(1);
	}
	int linenumber=40000,binnum=1000;
	int i=0,j=0;
	int num=linenumber/binnum;
	double *tmp1=(double*)malloc(sizeof(double)*linenumber);
	double *tmp2=(double*)malloc(sizeof(double)*linenumber);
	double *boxx=(double*)malloc(sizeof(double)*num);
	double *boxy=(double*)malloc(sizeof(double)*num);
	int *count=(int*)malloc(sizeof(int)*num);
	for(i=0;i<linenumber;i++)
	{
		fscanf(fp,"%lf %lf",&tmp1[i],&tmp2[i]);
	}
	double min=tmp1[0];
	double max=tmp1[linenumber-1];
	//double max=1.4;
	double delta=(max-min)/num;
	double lb=0,ub=0;
	//printf("%lf %lf\n",tmp1[0],tmp2[0]);
	for(i=0;i<num;i++)
	{
		boxx[i]=0;
		boxy[i]=0;
		count[i]=0;
	}

	for(i=0;i<num;i++)
	{
		lb=min+i*delta;
		ub=min+(i+1)*delta;
		boxx[i]=0;
		boxy[i]=0;
		for(j=0;j<linenumber;j++)
		{
			if(tmp1[j]>=lb && tmp1[j]<=ub)
			{
				boxx[i]+=tmp1[j];
				boxy[i]+=tmp2[j];
				count[i]++;
			}
		}
		if(count[i]!=0)
		{
			boxx[i]/=count[i];
			boxy[i]/=count[i];
		}
	}
	
	for(i=0;i<num;i++)
	{
		if(count[i]!=0)
			fprintf(rst,"%lf %lf\n",boxx[i],boxy[i]);
	}
	
	fclose(rst);
	fclose(fp);
	free(tmp1);
	free(tmp2);
	free(boxx);
	free(boxy);
	free(count);
	return 0;
}