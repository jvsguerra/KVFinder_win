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

This file contains the main functionality of KVFinder. It is responsible for controlling the execution flux of KVFinder software, according to the user definitions.
A brief explanation of the most important parts may be found in the source code*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "dictionary.h"
#include "PDBprocessing.h"
#include "Matrixprocessing.h"


int main()
{

	double h,***M;		/* A hash that will contain the radii dictionary, the step variable h and the grid */ 
	char TABLE[TABLE_SIZE][RES_SIZE];               /* A table that will hold the 3 letters name of the residues in the dictionary     */ 
	int ***A, ***S, m,n,o,i,j,k, ncav,tablesize;    /* Two grids representing empty spaces and surface points along with several counters     */
	dict *DIC[TABLE_SIZE];			/* A table that will represent the Dictionary                                      */
	FILE *parameters_file;		  	        /* A file pointer for the parameter input dataset                                         */ 		
	char AUX[50];					/* An auxiliar string variable                                                            */
	char dictionary_name[DIC_NAME_MAX],table_name[DIC_NAME_MAX]; /* File names                                 */
	char PDB_NAME[NAME_MAX],LIGAND_NAME[NAME_MAX],OUTPUT[NAME_MAX];
	int mode,cont,surface,whole_protein,redimension,boxmode,kvpmode;
	double cavfilter,probe,limit,norm1,norm2,norm3,probe_out;
	atom *p;
	double bX1, bY1, bZ1, bX2, bY2, bZ2, bX3, bY3, bZ3, bX4, bY4, bZ4;

	/************* Parameters input **************/

	parameters_file = fopen ("Parameters.txt","r");
		
	if(parameters_file == NULL)
			printf("ERROR! Parameters file not found!\n");
	else
	{
		fgets(dictionary_name, 500 , parameters_file);
		strtok(dictionary_name, "\n");
	
			/*      Dictionaries Loading     */
	    
		tablesize=define_table(TABLE,dictionary_name);
		if(dictionary_load(DIC,tablesize,dictionary_name))      	
		{
			fgets(PDB_NAME, 500 ,parameters_file);
			strtok(PDB_NAME, "\n");	
			for(cont=0;cont<1;cont++)
			{   /*  Reads some of the user defined parameters */
				fgets(OUTPUT, 500 ,parameters_file);
			    strtok(OUTPUT, "\n");
				fscanf(parameters_file, "%d",&whole_protein);
				fscanf(parameters_file, "%lf",&h);
				fscanf(parameters_file, "%d",&redimension);
				fscanf(parameters_file, "%lf %lf %lf",&bX1,&bY1,&bZ1);
				fscanf(parameters_file, "%lf %lf %lf %lf %lf %lf",&bX2,&bY2,&bZ2,&bX3,&bY3,&bZ3);
				fscanf(parameters_file, "%lf %lf %lf",&bX4,&bY4,&bZ4);
				fscanf(parameters_file, "%lf %lf %lf",&X1,&Y1,&Z1);
				fscanf(parameters_file, "%lf %lf %lf %lf %lf %lf",&X2,&Y2,&Z2,&X3,&Y3,&Z3);
				fscanf(parameters_file, "%lf %lf %lf",&X4,&Y4,&Z4);
				fscanf(parameters_file, "%lf",&probe);
				fscanf(parameters_file, "%lf",&probe_out);
				fscanf(parameters_file, "%d",&boxmode);
				if(whole_protein)
				{   /* If the whole protein mode is on, defines the search space around the protein */
					X1=999999;
					Y1=999999;
					Z1=999999;
					X2=0.0;
					Y3=0.0;
					Z4=0.0;
					PDB_load2(DIC,tablesize,TABLE,PDB_NAME);
					for(p=v;p!=NULL;p=p->next)
					{
						if(p->x < X1)
							X1=(p->x);
						if(p->x > X2)
							X2=(p->x);
						if(p->y < Y1)
							Y1=(p->y);
						if(p->y > Y3)
							Y3=(p->y);
						if(p->z < Z1)
							Z1=(p->z);
						if(p->z > Z4)
							Z4=(p->z);
					}
					
					X1-=(1.5*probe_out);	X2+=(1.5*probe_out);	
					Y1-=(1.5*probe_out);	Y3+=(1.5*probe_out);
					Z1-=(1.5*probe_out);	Z4+=(1.5*probe_out);
					X3=X1;		X4=X1;	
					Y2=Y1;		Y4=Y1;
					Z2=Z1;		Z3=Z1;		
					free_atom();	
				}
				/* If it will be used a user defined search space, withou the Probe Out Adjustment define its limits*/
				if(boxmode == 0 && whole_protein == 0){
				
				    X1 = bX1; X2 = bX2; X3 = bX3; X4 = bX4;
				    Y1 = bY1; Y2 = bY2; Y3 = bY3; Y4 = bY4;
				    Z1 = bZ1; Z2 = bZ2; Z3 = bZ3; Z4 = bZ4;
				    	    
				}
                /*Defines the grid inside the search space*/		
				norm1=sqrt( (X2-X1)*(X2-X1) + (Y2-Y1)*(Y2-Y1) + (Z2-Z1)*(Z2-Z1));
				norm2=sqrt( (X3-X1)*(X3-X1) + (Y3-Y1)*(Y3-Y1) + (Z3-Z1)*(Z3-Z1)); 			
				norm3=sqrt( (X4-X1)*(X4-X1) + (Y4-Y1)*(Y4-Y1) + (Z4-Z1)*(Z4-Z1));
				
				m=(int)(norm1/h)+2; 	n=(int)(norm2/h)+2; 		o=(int)(norm3/h)+2;
                /*If step redimensioning is on, increases the grid step size until the number of points on the grid be smaller than 80*80*80*/
				if(redimension)
				{
					while(m*n*o > 80*80*80)
					{
 						h+=0.1;	
						m=(int)(norm1/h); 	
						n=(int)(norm2/h); 		
						o=(int)(norm3/h);
						printf("Step size too small! Redimensioning... New step size = %lf\n",h);
					}
				}
				
				/*Calculates data used on the spatial manipulation of the protein*/				
				sina=(Y4-Y1)/norm3;	cosa=(Y3-Y1)/norm2;
				cosb=(X2-X1)/norm1;	sinb=(Z2-Z1)/norm1;

				printf("m: %d n: %d o: %d\n",m,n,o);
				printf("sina: %.4lf\tsinb: %.4lf\t\ncosa: %.4lf\tcosb: %.4lf\n",sina,sinb,cosa,cosb); 
				
					/*      Protein Coordinates Extraction         */

				if(PDB_load(DIC,tablesize,TABLE,PDB_NAME,probe,m,n,o,h))
		        	{
					
							/*        Matrix Allocation          */
					A = (int ***) malloc(sizeof(int**) * m);
					for(i=0;i<m;i++)
					{
						A[i] = (int **) malloc (sizeof (int *) * n);
						for(j=0;j<n;j++)
							A[i][j]=(int *) malloc (sizeof (int) * o);
					}
				
					S = (int ***) malloc(sizeof(int**) * m);
					for(i=0;i<m;i++)
					{
						S[i] = (int **) malloc (sizeof (int *) * n);
						for(j=0;j<n;j++)
							S[i][j]=(int *) malloc (sizeof (int) * o);
					}
			
					M = (double ***) malloc(sizeof(double**) * m);
					for(i=0;i<m;i++)
					{
						M[i] = (double **) malloc (sizeof (double *) * n);
						for(j=0;j<n;j++)
							M[i][j]=(double *) malloc (sizeof (double) * o);
					}
	
	
		
						/*          Matrix Initialization           */
	
	
					Matrix_initialize(A,m,n,o);
					Matrix_initialize(S,m,n,o);
					Matrix_initialize2(M,m,n,o);
					
					/*Reads more parameters defined by the user*/
					fscanf(parameters_file, "%lf",&cavfilter);
					fscanf(parameters_file, "%d",&surface);
					fscanf(parameters_file, "%d",&kvpmode);					
					fscanf(parameters_file, "%d",&mode);					
					
		
						/*************** Matrix Filling ***************/
				    int i, j, k;
					if(whole_protein)	
					{
					
					 Matrix_fill2(S,M,m,n,o,h,probe_out);
					 Matrix_initialize2(M,m,n,o);
					 Matrix_surf(S,m,n,o,h,probe_out);
                    				 
					}
					
					if(boxmode && !whole_protein)	
					{
					
					 Matrix_fill2(S,M,m,n,o,h,probe_out);
					 Matrix_initialize2(M,m,n,o);
					 Matrix_surf(S,m,n,o,h,probe_out);
					 			 
					}
					
					Matrix_fill2(A,M,m,n,o,h,probe);
									
					if(surface){
						Matrix_surf(A,m,n,o,h,probe);
				    }
				                     
					if(whole_protein)
					{
						Matrix_subtract(S,A,m,n,o);
						Matrix_initialize(S,m,n,o);	
					}
					
					if(boxmode && !whole_protein)
					{
						Matrix_subtract(S,A,m,n,o);
						Matrix_initialize(S,m,n,o);	
					}

					/*************** Matrix Adjusting **************/


					if(mode) 			
					{
						free_atom();
						fscanf(parameters_file, "%s",LIGAND_NAME);
						fscanf(parameters_file, "%lf",&limit);
						PDB_load(DIC,tablesize,TABLE,LIGAND_NAME,probe,m,n,o,h);
				        Matrix_adjust(A,m,n,o,h,limit);
				        PDB_load(DIC,tablesize,TABLE,PDB_NAME,probe,m,n,o,h);
				    }
                 

                    /*************Apply Mask to the Matrix*****************/
                    
                    if(boxmode)
                    {
                        Matrix_filter(A,m,n,o,h,bX1,bY1,bZ1,bX2,bY2,bZ2,norm1);
                        Matrix_filter(S,m,n,o,h,bX1,bY1,bZ1,bX2,bY2,bZ2,norm1);
                    }

					/******** Computing Volume and Grouping Cavities ********/

	
	
					ncav = DFS_search(A,m,n,o,h,cavfilter)-1;
   
					
								
					/*************** Obtain Matrix Surface  **************/
					
				    if(!kvpmode){
				        
					    Matrix_surface(A,S,m,n,o,h);
					    ncav = Area_search(S,m,n,o,h,cavfilter,ncav);
					    Matrix_search(A,m,n,o,h,probe);
					    
					    
					    Matrix_export(S,m,n,o,h,OUTPUT);
					    //Matrix_export2(M,S,m,n,o,h,OUTPUT);
					    
	                }
	                else{

					    Matrix_export(A,m,n,o,h,OUTPUT);

					      
	                }
					free_matrix(A,m,n,o);
					free_matrix(S,m,n,o);
					free_matrix2(M,m,n,o);
				}
			}
		}

		fclose(parameters_file);
		free_node();
		free_atom();
	}

	return 0;
}
