#include<stdio.h>
#include<stdlib.h>
typedef struct {
	int id;
	int n3;
	int n4;
	int n5;
	int n6;
}Atom;

int main()
{
	FILE *fp;
	if((fp=fopen("cnindex15%.dat","r"))==NULL)
	{
		printf("cannot open file\n");
		exit(1);
	}
	int atomnum=0,i=0,conf=1,j=0;
	atomnum=40000;
	Atom **atom;
	atom=(Atom**)malloc(sizeof(Atom*)*conf);
	for(i=0;i<conf;i++)
		atom[i]=(Atom*)malloc(sizeof(Atom)*atomnum);

	for(i=0;i<conf;i++)
	{
		for(j=0;j<atomnum;j++)
			fscanf(fp,"%d %d %d %d %d",&atom[i][j].id,&atom[i][j].n3,&atom[i][j].n4,&atom[i][j].n5,&atom[i][j].n6);
	}
	fclose(fp);
	if((fp=fopen("LFFS15%.dat","w"))==NULL)
	{
		printf("cannot open file\n");
		exit(1);
	}
	
	double **lffs=(double**)malloc(sizeof(double*)*conf);
	for(i=0;i<conf;i++)
		lffs[i]=(double*)malloc(sizeof(double)*atomnum);

	for(i=0;i<conf;i++)
	{
		for(j=0;j<atomnum;j++)
		{
			lffs[i][j]=double(atom[i][j].n5)/double(atom[i][j].n3+atom[i][j].n4+atom[i][j].n5+atom[i][j].n6);
		}
	}
	
	double *lf=(double*)malloc(sizeof(double)*atomnum);

	for(i=0;i<atomnum;i++)
	{
		lf[i]=0.0;
		for(j=0;j<conf;j++)
		{
			lf[i]+=lffs[j][i];
		}
		lf[i]/=conf;
		fprintf(fp,"%lf\n",lf[i]);
	}

	free(atom);
	free(lffs);
	free(lf);
	
	fclose(fp);
	return 0;
}