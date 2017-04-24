#include<stdio.h>
#include<stdlib.h>
#include<algorithm>
#include<math.h>
# define PI       3.14159265358979323846  /* pi */
using namespace std;
//define Atom type including id,type,x,y,z
typedef struct {
	int id;
	int type;
	double x;
	double y;
	double z;
}Atom;

//compare function for sorting
int cmp(Atom a,Atom b)
{
	return a.id < b.id;
}

//calculate the distance of two atoms
double distance(double x1,double y1,double z1,double x2,double y2,double z2)
{
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}


int main()
{
	FILE *fp;
	char str[30];
	int i=0,atomnum=0,j=0;
	if((fp=fopen("296.atom","r"))==NULL)
	{
		printf("cannot open file\n");
		exit(1);
	}
	for(i=0;i<7;i++)
	{
		fscanf(fp,"%s",str);
	}
	fscanf(fp,"%d",&atomnum);
	//printf("atom num is %d\n",atomnum);
	for(i=0;i<6;i++)
	{
		fscanf(fp,"%s",str);
	}
	double xl,xh,yl,yh,zl,zh;
	fscanf(fp,"%lf %lf %lf %lf %lf %lf",&xl,&xh,&yl,&yh,&zl,&zh);
	//printf("xl,xh,yl,yh,zl,zh:%lf %lf %lf %lf %lf %lf",xl,xh,yl,yh,zl,zh);
	for(i=0;i<7;i++)
	{
		fscanf(fp,"%s",str);
	}
	Atom *atom=(Atom*)malloc(sizeof(Atom)*atomnum);
	//read the atom id,type,x,y,z
	for(i=0;i<atomnum;i++)
	{
		fscanf(fp,"%d %d %lf %lf %lf",&atom[i].id,&atom[i].type,&atom[i].x,&atom[i].y,&atom[i].z);
	}
	fclose(fp);
	//order the atoms with id
	sort(atom,atom+atomnum,cmp);
	//number density of type1 and type 2
	double numdensity1=0;
	double numdensity2=0;
	double volume=(xh-xl)*(yh-yl)*(zh-zl);
	for(i=0;i<atomnum;i++)
	{
		if(atom[i].type==1)
			numdensity1++;
		else 
			numdensity2++;
	}
	
	numdensity1=numdensity1/volume;
	numdensity2=numdensity2/volume;
	//printf("%lf %lf\n",numdensity1,numdensity2);

	int binnum=100;
	double radius=6.5;
	double deltr=radius/(double)binnum;
	double deltV=4.0/3.0*PI*deltr*deltr*deltr;
	double *nv=(double*)malloc(sizeof(double)*binnum);
	double **dis1;//different type of atoms
	double **dis2;//same type of atoms
	double r=0;//store two particles distance
	int index=0;//store which shell should it be stored 
	double *s2=(double*)malloc(sizeof(double)*atomnum);
	dis1=(double** )malloc(sizeof(double*)*atomnum);
	dis2=(double** )malloc(sizeof(double*)*atomnum);
	for(i=0;i<atomnum;i++)
	{
		dis1[i]=(double*)malloc(sizeof(double)*binnum);
		dis2[i]=(double*)malloc(sizeof(double)*binnum);
		s2[i]=0;
	}

	for(i=0;i<atomnum;i++)
	{
		for(j=0;j<binnum;j++)
		{
			dis1[i][j]=0;
			dis2[i][j]=0;
		}
	}
	//printf("%d %lf %lf %lf\n",atom[0].id,atom[0].x,atom[0].y,atom[0].z);
	//printf("%d %lf %lf %lf\n",atom[1].id,atom[1].x,atom[1].y,atom[1].z);
	//printf("%lf",distance(atom[0].x,atom[0].y,atom[0].z,atom[1].x,atom[1].y,atom[1].z));
	for(i=0;i<binnum;i++)
		nv[i]=deltV*((i+1)*(i+1)*(i+1)-i*i*i);

	for(i=0;i<atomnum;i++)
	{
		for(j=0;j<atomnum;j++)
		{
			if(i!=j)
			{
				r=distance(atom[i].x,atom[i].y,atom[i].z,atom[j].x,atom[j].y,atom[j].z);
				if(r<radius)
				{
					index=int(r/deltr);
					if(atom[i].type != atom[j].type)
					{
						dis1[i][index]=dis1[i][index]+1.0;
					}
					else
					{
						dis2[i][index]=dis2[i][index]+1.0;
					}
				}
			}
		}
	}
	/*
	fp=fopen("dis.dat","w");
	if(fp==NULL)
	{
		printf("file cannot open!\n");
		exit(1);
	}
	for(i=0;i<atomnum;i++)
	{
		fprintf(fp,"id:%4d\n",atom[i].id);
		for(j=0;j<binnum;j++)
			fprintf(fp,"%d %lf %lf\n",atom[i].id,dis1[i][j],dis2[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	*/
	for(i=0;i<atomnum;i++)
	{
		if(atom[i].type==1)
		{
			for(j=0;j<binnum;j++)
			{
				dis1[i][j]=dis1[i][j]/nv[j]/numdensity2;
				dis2[i][j]=dis2[i][j]/nv[j]/numdensity1;
			}
		}
		else
		{
			for(j=0;j<binnum;j++)
			{
				dis1[i][j]=dis1[i][j]/nv[j]/numdensity1;
				dis2[i][j]=dis2[i][j]/nv[j]/numdensity2;
			}
		}
	}
	/*
	fp=fopen("gr.data","w");
	if(fp==NULL)
	{
		printf("file cannot open!\n");
		exit(1);
	}

	for(i=0;i<atomnum;i++)
	{
		fprintf(fp,"id:%4d\n",atom[i].id);
		for(j=0;j<binnum;j++)
			fprintf(fp,"%lf\n",dis1[i][j]+dis2[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	*/
	//compute every atom's S2:
	for(i=0;i<atomnum;i++)
	{
		if(atom[i].type==1)
		{
			for(j=0;j<binnum;j++)
			{
				if(dis1[i][j]!=0)
					s2[i]+=pow((j+0.5)*deltr,2)*numdensity2*(deltr*(dis1[i][j]*log(dis1[i][j])-(dis1[i][j]-1)));//?
				if(dis2[i][j]!=0)
					s2[i]+=pow((j+0.5)*deltr,2)*numdensity1*(deltr*(dis2[i][j]*log(dis2[i][j])-(dis2[i][j]-1)));
			}
		}
		else
		{
			for(j=0;j<binnum;j++)
			{
				if(dis1[i][j]!=0)
					s2[i]+=pow((j+0.5)*deltr,2)*numdensity1*(deltr*(dis1[i][j]*log(dis1[i][j])-(dis1[i][j]-1)));
				if(dis2[i][j]!=0)
					s2[i]+=pow((j+0.5)*deltr,2)*numdensity2*(deltr*(dis2[i][j]*log(dis2[i][j])-(dis2[i][j]-1)));
			}
		}	
	}
	fp=fopen("296.data","w");
	if(fp==NULL)
	{
		printf("file cannot open!\n");
		exit(1);
	}
	for(i=0;i<atomnum;i++)
	{
		fprintf(fp,"%lf\n",(-0.5)*s2[i]);
	}
	fclose(fp);
	free(s2);
	free(nv);
	free(dis1);
	free(dis2);
	
	return 0;
}