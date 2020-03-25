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

This file contains all the grid functions utilized by KVFinder. The KVFinder grid data is represented by a matrix.
The most complex functions are briefly explained in the source code.*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"dictionary.h"
#include"PDBprocessing.h"
#include"Matrixprocessing.h"
#include<string.h>
#include<unistd.h>

double max( double a, double b)
{
	if(a>b)
		return a;
	else
		return b;
}

double min( double a, double b)
{
	if(a<b)
		return a;
	else
		return b;
}

/*Checks if a given cavity point on the grid is next to a protein point*/
int check_pos(int ***A, int i, int j, int k, int m, int n, int o)
{
	int a,b,c,flag=0;
	
	
	for(a=i-1; a<=i+1 && !flag ; a++)
	   for(b=j-1; b<=j+1 && !flag; b++)
		for(c=k-1; c<=k+1 && !flag; c++)
		{  
	        if(a < 0 || b < 0 || c < 0 || a > m-1 || b > n-1 || c > o-1)
		            ;
		    else    
		        if(A[a][b][c]==0 || A[a][b][c]==-2)
		            flag=1;
	    }
	if(flag)
		return 1;
	else
		return 0;
}

/*Checks if a given cavity point on the grid is next to a protein point*/
int check_pos2(int ***A, int i, int j, int k, int m, int n, int o)
{
	int a,b,c,flag=0;
	
	
	
	for(a=i-1; a<=i+1 && !flag ; a++)
	   for(b=j-1; b<=j+1 && !flag; b++)
		for(c=k-1; c<=k+1 && !flag; c++){
		    if(a < 0 || b < 0 || c < 0 || a > m-1 || b > n-1 || c > o-1)
		        ;
		    else    
		        if(A[a][b][c]==0)
			        flag=1;
	    }
	if(flag)
		return A[i][j][k];
	else
		return 0;
}

/*When the cavities are separated, each cavity is identified by a tag. This checks, for a given identified cavity point, if there is a unidenfied cavity point around it*/
int check_pos3(int ***A, int i, int j, int k, int m, int n, int o)
{
	int a,b,c,flag=0;
	int a2,b2,c2;
	for(a=i-1; a<=i+1 && !flag ; a++)
	   for(b=j-1; b<=j+1 && !flag; b++)
		for(c=k-1; c<=k+1 && !flag; c++)
		{
		    if(a < 0 || b < 0 || c < 0 || a > m-1 || b > n-1 || c > o-1)
		            ;
		    else    
	    		if(A[a][b][c]>1){
	    		     flag=1;
	    		     a2 = a; b2 = b; c2 = c;
	    		     }
	    }
	if(flag)
		return A[a2][b2][c2];
	else
		return 0;
}

/*double getcharge(double ***M, int i, int j, int k)
{

    double Ke=8.9875517873681764;
    double charge=0.0;
    
    charge = Ke*M[i][j][k];
	return charge;	
}
*/
void Matrix_initialize(int ***A,int m,int n, int o)
{
	int i,j,k;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			for(k=0;k<o;k++)
				A[i][j][k]=1;
}

void Matrix_initialize2(double ***A,int m,int n, int o)
{
	int i,j,k;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			for(k=0;k<o;k++)
				A[i][j][k]=0.0; 
}

/*Marks the grid leaving a probe size around the protein*/


void Matrix_fill2(int ***A,double ***M, int m, int n, int o, double h,double probe)
{
	int i,j,k,flag=1,aux;
	double distance,x,y,z,xaux,yaux,zaux,H;
	double x1,y1,z1;
	atom *p;
	
	for(p=v;p!=NULL&&flag;p=p->next)
	{
		x1=(p->x-X1)/h;
		y1=(p->y-Y1)/h;
		z1=(p->z-Z1)/h;
			
		xaux= x1*cosb + z1*sinb;
		yaux= y1;
		zaux= -x1*sinb + z1*cosb;

		x1= xaux;	
		y1= yaux*cosa - zaux*sina;
		z1= yaux*sina + zaux*cosa;
		
		H=(probe+p->radius)/h;
		
		for(i=floor(x1-H);i<=ceil(x1+H);i++)
			for(j=floor(y1-H);j<=ceil(y1+H);j++)		
				for(k=floor(z1-H)+1;k<=ceil(z1+H)+1;k++)	
				{
					distance=sqrt((i-x1)*(i-x1) + (j-y1)*(j-y1) + (k-z1)*(k-z1));
					if(distance < H && i>=0 && i<m && j>=0 && j<n && k>=0 && k<o )
					{
						A[i][j][k]=0;
						//M[i][j][k]+=p->charge/(distance*h);
						
					}
				}
	}
}
/*Creates a list of residues contacting each cavity*/
void Matrix_search(int ***A, int m, int n, int o, double h,double probe)
{
	int i,j,k,flag=1,aux,oldnum, cont;
	double distance,x,y,z,xaux,yaux,zaux,H;
	double x1,y1,z1;
	char lastChain = '!';
	atom *p;
	FILE *arquivo, *arq_aux;
	char filename[31] = "KVFinderTemporaryCavityFileK";
	char LINE[20];
	char *LINE2;
	LINE2 = (char*  )malloc(20*sizeof(char));
    cont = 0;  
	oldnum=0;
	for(p=v;p!=NULL&&flag;p=p->next)
	{
	    x1=(p->x-X1)/h;
	    y1=(p->y-Y1)/h;
	    z1=(p->z-Z1)/h;
			
	    xaux= x1*cosb + z1*sinb;
	    yaux= y1;
	    zaux= -x1*sinb + z1*cosb;

	    x1= xaux;	
	    y1= yaux*cosa - zaux*sina;
	    z1= yaux*sina + zaux*cosa;
		
	    H=(probe+p->radius)/h;
		  
	    for(i=floor(x1-H);i<=ceil(x1+H);i++)
		    for(j=floor(y1-H);j<=ceil(y1+H);j++)		
			    for(k=floor(z1-H);k<=ceil(z1+H);k++)
				    if(i<m-1&&i>0&&j<n-1&&j>0&&k<o-1&&k>0)
					    //if(p->resnum!=oldnum && A[i][j][k]>=2)
					    if(A[i][j][k]>=2)
					    {
					        filename[28] = 65+(((A[i][j][k]-2)/26)%26);
					        filename[29] = 65+((A[i][j][k]-2)%26);
					        filename[30] = '\0';
					        arq_aux = fopen(filename , "a+");
					        fprintf(arq_aux,"Cavity K%c%c %d %c\n", 65+(((A[i][j][k]-2)/26)%26),65+((A[i][j][k]-2)%26),p->resnum, p->chain);
							oldnum=p->resnum;
							fclose(arq_aux);
							cont++;
							strcpy(LINE, "");
						}
	}
	
	//merge files!
	arquivo = fopen("cavres.txt","w");
	aux = 2;
	while(cont > 0){
        filename[28] = 65+(((aux-2)/26)%26);
		filename[29] = 65+((aux-2)%26);
		if (access(filename, F_OK) == 0){
		    arq_aux=fopen(filename, "r");
		    while(get_line(arq_aux,LINE)){   
		        if(strncmp(LINE2, LINE,15)!=0){
		            for(i=0; LINE[i]!='\n';i++)
		                fprintf(arquivo,"%c",LINE[i]);
		            fprintf(arquivo, "%c", '\n'); 
		            cont--;
		            free(LINE2);
		            LINE2 = (char*  )malloc(20*sizeof(char));   
		            strcpy(LINE2, LINE);
		        }   
		    }
		    fclose(arq_aux);
		    remove(filename);
		}
		cont--;
		aux++;
	}  
	fclose(arquivo);

}


/*Passes the probe around the surface of the protein*/
void Matrix_surf(int ***A, int m, int n, int o, double h, double radius)
{
	int i,j,k,i2,j2,k2,aux;
	double norm;
	aux = ((int)radius/h)+1;
	
	
	for(i=0;i<m;i++)
	  for(j=0;j<n;j++)
	    for(k=0;k<o;k++)
	    {
		if(A[i][j][k]==1 && check_pos(A,i,j,k,m,n,o))
		{
			for(i2=i-aux;i2<=i+aux;i2++) 
				for(j2=j-aux;j2<=j+aux;j2++) 
					for(k2=k-aux;k2<=k+aux;k2++)
						{
						    if(i2 > 0 && j2 > 0 && k2 > 0 && i2 < m-1 && j2 < n-1 && k2 < o-1)
						    {
							    norm=(i-i2)*(i-i2)+(j-j2)*(j-j2)+(k-k2)*(k-k2);
							    norm=sqrt(norm);
							    if(norm < (radius/h) && A[i2][j2][k2] == 0)
								    A[i2][j2][k2]=-2;
							}
						}
						
                
		}		
	    }
     
	    for(i=0;i<m;i++)
    	  	for(j=0;j<n;j++)
	    		for(k=0;k<o;k++)
			{
				if(A[i][j][k]==-2)
					A[i][j][k]=1;
			}
}

/*Marks the surface points on the grid*/
void Matrix_surface(int ***A, int ***S, int m, int n, int o, double h)
{
	int i,j,k;
	
	
	for(i=0;i<m;i++)
    	  for(j=0;j<n;j++)
	    for(k=0;k<o;k++)
	    {
		if(A[i][j][k]>=1)
			S[i][j][k]=check_pos2(A,i,j,k, m, n, o);
		else
			S[i][j][k]=0;	
	    }


}


void Matrix_adjust(int ***A, int m, int n, int o, double h,double limit)
{
	int i,j,k,flag,aux;
	double distance,x,y,z,xaux,yaux,zaux;
	atom *p;
	for(i=0;i<m;i++)
    	  for(j=0;j<n;j++)
	    for(k=0;k<o;k++)
	    {
		flag=0;
		for(p=v;p!=NULL&&!flag;p=p->next)
		{
			x=i*h;
			y=j*h;
			z=k*h;
			xaux= x*cosb + y*sina*sinb - z*cosa*sinb;
			yaux=          y*cosa      + z*sina;
			zaux= x*sinb - y*sina*cosb + z*cosa*cosb;	
			xaux+=X1; yaux+=Y1; zaux+=Z1;
			distance=sqrt((xaux-p->x)*(xaux-p->x) + (yaux-p->y)*(yaux-p->y) + (zaux-p->z)*(zaux-p->z));
			if(distance<limit)
				flag=1;
			if(distance<(-1.0*limit) && limit < 0.0)
				flag=1;
		}
		if(flag==0 && A[i][j][k] && limit > 0.0)
			A[i][j][k]=-1;
		if(flag==1 && A[i][j][k] && limit < 0.0)
			A[i][j][k]=-1;

	    }

}

/*Makes a spatial filter, used on the probe out adjutment, that requires a analysis of points outside the user defined search space. These points are excluded here*/
void Matrix_filter(int ***A, int m, int n, int o, double h, double bX1, double bY1, double bZ1, double bX2, double bY2, double bZ2, double norm1)
{

    int i, j, k;
    double aux, normB;
    
    
    normB=sqrt( (bX2-bX1)*(bX2-bX1) + (bY2-bY1)*(bY2-bY1) + (bZ2-bZ1)*(bZ2-bZ1));
    
    aux = (norm1 - normB)/2;
    
    aux = (int)(aux/h);
 
    
    for(i=0; i<=aux;i++)
        for(j=0;j<n;j++)
            for(k=0;k<o;k++)
                A[i][j][k] = -3;
    
    for(i=m-1;i>=m-aux-1;i--)
        for(j=0;j<n;j++)
            for(k=0;k<o;k++)
                A[i][j][k] = -3;
    
    
    for(j=0; j<=aux;j++)
         for(i=0;i<m;i++)
            for(k=0;k<o;k++)
                A[i][j][k] = -3;
    
    for(j=n-1;j>=n-aux-1;j--)
         for(i=0;i<m;i++)
            for(k=0;k<o;k++)
                A[i][j][k] = -3;
    
    for(k=0; k<=aux;k++)
         for(j=0;j<n;j++)
            for(i=0;i<m;i++)
                A[i][j][k] = -3;
    for(k=o-1;k>=o-aux-1;k--)
         for(j=0;j<n;j++)
            for(i=0;i<m;i++)
                A[i][j][k] = -3;                
              
}

/*Prints the grid points in a PDB file*/
void Matrix_export(int ***A, int m, int n, int o,double h,char Output_name[NAME_MAX])
{
	int i,j,k,count=1,flag=1, controle=1, tag=2;
	double x,y,z,xaux,yaux,zaux;
	FILE *output;
	
	
	for(i=0;Output_name[i]!='\0';i++);
	
	Output_name[i]='.'; Output_name[i+1]='p'; Output_name[i+2]='d'; Output_name[i+3]='b'; Output_name[i+4]='\0';
	output = fopen(Output_name,"w");

    while(controle==1){
        controle=0;
	    for(i=0;i < m;i++)
	        for(j=0;j < n;j++)
	            for(k=0;k < o ; k++)
	            {
		            if(A[i][j][k] == tag)
		            {   
		                controle=1;
			            x=i*h;
			            y=j*h;
			            z=k*h;
			            xaux= x*cosb + y*sina*sinb - z*cosa*sinb;
			            yaux=          y*cosa      + z*sina;
			            zaux= x*sinb - y*sina*cosb + z*cosa*cosb;	
			            xaux+=X1; yaux+=Y1; zaux+=Z1;

			            fprintf(output,"ATOM  %5.d  H   K%c%c   259    %8.3lf%8.3lf%8.3lf   1.00  0.00\n", count, 65+(((A[i][j][k]-2)/26)%26),65+((A[i][j][k]-2)%26), xaux,yaux,zaux);
  
			            count++;
			            if(count == 100000) count = 1;
		            }
	            }
    tag++;
    
    }
}



/*Prints the grid points in a PQR file*/
/*void Matrix_export2(double ***M, int ***A, int m, int n, int o,double h,char Output_name[NAME_MAX])
{
	int i,j,k,count=1,flag=1, controle = 1, tag=2;
	FILE *output;
	double charge,radius=1.4;	
	double x,y,z,xaux,yaux,zaux;

	for(i=0;Output_name[i]!='\0';i++);
	Output_name[i-2]='q'; Output_name[i-1]='r';
	output = fopen(Output_name,"w");
    while(controle == 1){
        controle = 0;
        for(i=0;i < m;i++)
	        for(j=0;j < n;j++)
	            for(k=0;k < o ; k++)
	            {
		            if(A[i][j][k]==tag)
		            {
		                controle = 1;
			            charge=getcharge(M,i,j,k);
			            x=i*h;
			            y=j*h;
			            z=k*h;
			            xaux= x*cosb + y*sina*sinb - z*cosa*sinb;
			            yaux=          y*cosa      + z*sina;
			            zaux= x*sinb - y*sina*cosb + z*cosa*cosb;	
			            xaux+=X1; yaux+=Y1; zaux+=Z1;
			            //fprintf(output,"ATOM  %5.d  H   ASP %c 259    %8.3lf%8.3lf%8.3lf  1.00%6.2lf\n", count, 65, xaux,yaux,zaux,(charge+1.0)*50);
			           
			                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                fprintf(output,"ATOM  %5.d  H   K%c%c   259    %8.3lf%8.3lf%8.3lf%7.2lf%6.2lf\n", count,65+(((A[i][j][k]-2)/26)%26),65+((A[i][j][k]-2)%26), xaux,yaux,zaux,charge,radius);
			            
			           
			            count++;
			            if(count == 100000) count = 1;
		            }
	            }
	    tag++;            
	}
    fclose(output);
	
}
*/		
void free_matrix(int ***A, int m, int n, int o)
{	
	int i, j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
			free(A[i][j]);
		
		free(A[i]);
	}

	free(A);
}

void free_matrix2(double ***A, int m, int n, int o)
{	
	int i, j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
			free(A[i][j]);
		
		free(A[i]);
	}

	free(A);
}

/*Mark as invalid all the cavities with volume smaller than a user defined threshold*/
void remove_cavity(int ***A, int m, int n, int o, int tag)
{
    int i, j, k;
    for(i=0; i < m;i++)
        for(j=0; j<n; j++)
            for(k=0; k<o; k++)
                if(A[i][j][k] == tag)
                    A[i][j][k] = 0;
}

/*Auxiliar function of the recursive search*/
void DFS(int ***A,int a,int b, int c, int m, int n, int o, int tag,double h)
{

    int i, j, k;
    //printf("%d %d %d %d               %d %d %d                  %f           %d\n", tag, a, b,c,m,n,o, volume, flagr);
    if(a==0 || a==m-1 || b==0 || b==n-1|| c==0 || c==o-1)
        return;
        
    if(A[a][b][c]==1 && !flagr){
        A[a][b][c] = tag;
        volume+=1.0;
        if(volume == 20000){
            flagr = 1;
        }
        if(!flagr){
            for(i=a-1 ; i<= a+1 ; i++)
			    for(j=b-1 ; j<= b+1 ; j++)
				    for(k=c-1; k<=c+1; k++)
					    DFS(A,i,j,k,m,n,o,tag,h);
					    }
    }        
    

}

/*makes a recursive search grouping the points in different cavities and calculating its volume
Due to memory restrictions, the recursion is divided for big cavities*/
int DFS_search(int ***A,int m, int n, int o, double h,double filter)
{
       int i,j,k,r,t,y,tag;
       flagr = 0;
       double vol_aux = 0.0;          
     
       
       node *p;
       V=NULL;
       
       tag=1;
       for(i=0;i<m;i++)
          for(j=0;j<n;j++)
	      for(k=0;k<o;k++)
	 	   if(A[i][j][k]==1)
		   {	
			
			volume=0.0;
			tag++;
			DFS(A,i,j,k,m,n,o,tag,h);
			vol_aux = volume;
			while(flagr){
			    vol_aux = 0;
                 for(r=1;r<m;r++)
                    for(t=0;t<n;t++)
	                    for(y=0;y<o;y++){
			            flagr = 0;
			            vol_aux+=volume;
			            volume=0.0;
			            if(A[r][t][y] == 1 && check_pos3(A,r,t,y,m,n,o) == tag)
			            {
			                DFS(A,r,t,y,m,n,o,tag,h);
			    
			            }
			    }
			}
			
            volume = vol_aux;
         
			if(volume*h*h*h < filter)
			{
				remove_cavity(A,m,n,o,tag);
				tag--;
			}
			else
			{

			    printf("Cavity K%c%c: Volume =  %lf Angstrons^3\n",65+(((A[i][j][k]-2)/26)%26),65+((A[i][j][k]-2)%26),(volume*h*h*h));
 
			   
			      			
			    p=malloc(sizeof(node));
				p->volume=volume*h*h*h;		
				p->pos=tag;
				if(V!=NULL)
				{
					p->next=V;
					V=p;
				}
				else
				{	
					p->next=NULL;
					V=p;
				}
						   
			}
		 }
		 

	return tag;
}

/*Non recursive version, counts tags in the surface matrix*/  
int Area_search(int ***A,int m, int n, int o, double h,double filter,int tag)
{

    int *vol, i, j, k; 
    
    vol = (int*)malloc(sizeof(int)*tag);
    for(i = 0; i < tag; i++)
    {
        vol[i] = 0;    
    }
    for(i=0;i<m;i++) 
       for(j=0;j<n;j++)
          for(k=0;k<o;k++)
          {

            if(A[i][j][k] >1)
            {
                vol[A[i][j][k]-2]++;
            }
          
          }
    for(i = 0; i < tag; i++)            
            printf("Cavity K%c%c: Area =  %lf Angstrons^2\n",65+(((i-2)/26)%26),65+((i)%26),(vol[i]*h*h));


	return tag;
}


/*Compares two matrices, selecting the points where the smaller probe passed and the bigger one doesn't*/
void Matrix_subtract(int ***S, int ***A, int m, int n, int o)
{
	int i,j,k,i2,j2,k2;
	double probe=4.0;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			for(k=0;k<o;k++)
			{	
				if(S[i][j][k])	
				{
					for(i2=floor(i-probe);i2<=i+probe;i2++)
						for(j2=floor(j-probe);j2<=j+probe;j2++)	
							for(k2=floor(k-probe);k2<=k+probe;k2++)
								if(i2>=0&&i2<m&&j2>=0&&j2<n&&k2>=0&&k2<o)
									A[i2][j2][k2]=-17;
				}
			}

}

void free_node()
{
	node *p;
	
	while(V!=NULL)
	{
	   p=V;
	   V=V->next;
	   free(p);
	}
}	

