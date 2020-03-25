/*    This file is part of KVFinder.

    KVFinder is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    KVFinder is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with KVFinder.  If not, see <http://www.gnu.org/licenses/>.
    
    
    The KVFinder software was developed by:
    Saulo Henrique Pires de Oliveira
    Felipe Augusto Nunes Ferraz
    Rodrigo Vargas Honorato
    Jose Xavier Neto
    Tiago Jose Paschoal Sobreira
    Paulo Sergio Lopes de Oliveira
    
    National Center of Energy and Material Research - CNPEM
    National Laboratory of Biosciences - LNBio
    Campinas, Brazil - P.O. Box 6192 - CEP 13083-970, Campinas - SP

    Contact: paulo.oliveira@lnbio.cnpem.br
    KVFinder website: http://lnbio.cnpem.br/bioinformatics/main/software/KVFinder

This file contains all the funtions used to process PDB files used as input for KVFinder*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"dictionary.h"
#include"PDBprocessing.h"
#include<math.h>

double max2(double a,double b)
{
	if(a<b)
		return b;
	else
		return a;
}

double min2(double a,double b)
{
	if(a>b)
		return b;
	else
		return a;
}

double max4(double a, double b,double c,double d)
{
	return max2(max2(a,b),max2(c,d));
}

double min4(double a,double b, double c, double d)
{
	return min2(min2(a,b),min2(c,d));
}
	
void prints(char S[50])
{
	int i;
	for(i=0;S[i] != '\0' && i <50;i++)
		printf("%c",S[i]);
	printf("\n");
}

void init(char S[50])
{
	int i;
	for(i=0; i <50; i++)
		S[i]='\0';
	
}

int get_line(FILE *arq, char LINE[100])
{
	int i;

	for(i=0 , LINE[i]=getc(arq); LINE[i]!=EOF && LINE[i]!='\n' && i < 100 ; i++, LINE[i]=getc(arq));
	
	if(i==100 && LINE[i]!='\n')
		for(;LINE[i]!='\n'&&LINE[i]!=EOF;LINE[i]=getc(arq));

	if(LINE[i]==EOF)
		return 0;
	else
		return 1;
}

void extract(char LINE[100],int a,int b,char S[50])
{
	int i;

	if(b-a>50)
		return;
	if(b>100)
		return;

	for(i=a;i<b;i++)
		S[i-a]=LINE[i];
	if(i<50)
		S[i]='\0';
}

char extractChain(char LINE[100])
{
	return LINE[21];
}


double convert(char S[50])
{
	double resp=0.0,aux=0.001;
	int i;

	for(i=0;S[i]!='\0' && i<50; i++);
	
	for(i=i-1;i>=0&&S[i]!=' ';i--)
	{
		if(S[i]!='-' && S[i]!='.')
		{
			resp+=aux*(double)(S[i]-48);
			aux*=10.0;
		}
		else if(S[i]=='-')
			resp=-resp;
	}
	
	return resp;
}

int PDB_load(dict *DIC[TABLE_SIZE],int tablesize,char TABLE[TABLE_SIZE][RES_SIZE],char PDB_NAME[NAME_MAX],double probe, double m, double n, double o,double h)
{
	int flag=1,i,j,k,aux,number;
	char AUX[50],LINE[100],c1,c2,RES[4],ATOM_TYPE[4],X[50],Y[50],Z[50],chain,ATOM_SYMBOL[2];	  /* Contains the name of the PDB file */
	FILE *arqPDB;
	double x,y,z, radius,M,N,O,x1,y1,z1;
	double xaux,yaux,zaux,xaux2,yaux2,zaux2,minX,minY,minZ,maxX,maxY,maxZ;
	v = NULL; 
	M=(double)m; 	N=(double)(n); 	O=(double)(o);
	arqPDB = fopen(PDB_NAME,"r");
		
	for(i=0;i<50;i++)
	{
		AUX[i]=0;
		LINE[i]=0;
		LINE[50+i]=0;
		X[i]=0;	
		Y[i]=0;
		Z[i]=0;
	}

	if(arqPDB == NULL)
	{
		flag = 0;	
		printf("Reading Error: PDB file not found!\n");
	}
	else
	{
		/* PDB Parsing */
		while(get_line(arqPDB,LINE))
		{
			
			extract(LINE,0,6,AUX);
			if(!strcmp(AUX,"ATOM  ")||!strcmp(AUX,"HETATM"))
			{
				for(i=12,j=0,k=0;i<17;i++)
				{
					if(LINE[i]!=' ')
					{
						ATOM_TYPE[j]=LINE[i];
						j++;
					}
				}

				ATOM_TYPE[j]='\0';
				RES[0]=LINE[17]; RES[1]=LINE[18]; RES[2]=LINE[19]; RES[3]='\0';				
				number=0;

				if(LINE[22]>47&&LINE[22]<58)
					number+=(LINE[22]-48)*1000;				
				if(LINE[23]>47&&LINE[23]<58)
					number+=(LINE[23]-48)*100;
				if(LINE[24]>47&&LINE[24]<58)
					number+=(LINE[24]-48)*10;
				if(LINE[25]>47&&LINE[25]<58)
					number+=(LINE[25]-48);


				init(AUX);
				extract(LINE,76,78,AUX);
				extract(LINE,30,38,X); 	
				extract(LINE,38,46,Y);
				extract(LINE,46,54,Z);
				chain = extractChain(LINE);
				strcpy(ATOM_SYMBOL,AUX);	
				trim(X,' ');  trim(Y,' ');  trim(Z,' ');
				x=convert(X); y=convert(Y); z=convert(Z); 

				radius=dictionary_consult_radius(RES,ATOM_TYPE,DIC,tablesize,TABLE,ATOM_SYMBOL);
				//charge=dictionary_consult_charge(RES,ATOM_TYPE,DIC,tablesize,TABLE,ATOM_SYMBOL);
	
	            x1=(x-X1)/h;
				y1=(y-Y1)/h;
				z1=(z-Z1)/h;
			
				xaux= x1*cosb + z1*sinb;
				yaux= y1;
				zaux= -x1*sinb + z1*cosb;

				x1= xaux;	
				y1= yaux*cosa - zaux*sina;
				z1= yaux*sina + zaux*cosa;
	
				if(x1 > 0.0-(probe+radius)/h  && x1 < (double)m + (probe + radius)/h && y1 > 0.0 -(probe +radius)/h && y1 < (double)n + (probe + radius)/h && z1 > 0.0 -(probe +radius)/h && z1 < (double)o + (probe + radius)/h)
			    insert_atom(x,y,z,radius,number,chain);
				
			}

		}
	}

	fclose(arqPDB);

	return flag;
} 

int PDB_load2(dict *DIC[TABLE_SIZE],int tablesize,char TABLE[TABLE_SIZE][RES_SIZE],char PDB_NAME[NAME_MAX])
{
	int flag=1,i;
	char AUX[50],LINE[100],X[50],Y[50],Z[50],chain;	  
	FILE *arqPDB;				/* Contains the name of the PDB file */
	double x,y,z;
	v = NULL; 

	arqPDB = fopen(PDB_NAME,"r");
		
	for(i=0;i<50;i++)
	{
		AUX[i]=0;
		X[i]=0;	
		Y[i]=0;
		Z[i]=0;
	}

	if(arqPDB == NULL)
	{
		flag = 0;	
		printf("Reading Error: PDB file not found!\n");
	}
	else
	{
		/* PDB Parsing */
		while(get_line(arqPDB,LINE))
		{
			
			extract(LINE,0,6,AUX);
			if(!strcmp(AUX,"ATOM  ")||!strcmp(AUX,"HETATM"))
			{
				extract(LINE,30,38,X); 	
				extract(LINE,38,46,Y);
				extract(LINE,46,54,Z);
				chain = extractChain(LINE);
				trim(X,' ');  trim(Y,' ');  trim(Z,' ');
				x=convert(X); y=convert(Y); z=convert(Z); 
				insert_atom(x,y,z,0.0,0,chain);
			}
		}
	}

	fclose(arqPDB);
	return flag;
} 


void insert_atom(double x,double y, double z, double radius,int resnumber,char chain)
{
	atom *p;
	
	if(v==NULL)
	{
 		v=malloc(sizeof(atom));
		v->x=x; v->y=y; v->z=z; v->radius=radius;// v->charge=charge;
		v->next=NULL;
	}
	else
	{
		p=malloc(sizeof(atom));
		p->x=x; p->y=y; p->z=z; p->radius=radius;/* p->charge=charge;*/ p->resnum=resnumber; p->chain=chain;
		p->next=v->next;
		v->next=p;
	}
}

void atom_list_print()
{
	int cont=0;
	atom *p;
	
	for(p=v;p!=NULL;p=p->next)
	{
		cont++;
		printf("%d:    %lf %lf %lf   %lf\n",cont,p->x,p->y,p->z,p->radius);//,p->charge);
	}
}

void free_atom()
{
	atom *p;
	
	while(v!=NULL)
	{
	   p=v;
	   v=v->next;
	   free(p);
	}
}	
