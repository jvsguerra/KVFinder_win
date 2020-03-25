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



This file contains the functions used by KVFinder to process the user defined dictionaries, fundamental data to the software*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "dictionary.h"

void trim(char S[50], char c)
{
	int i,j;
	
	for(i=0;S[i]!='\0';i++)
	{
		if(S[i]==c)
			for(j=i;S[j]!='\0';j++)
				S[j]=S[j+1];
	}
}

int define_table(char TABLE[TABLE_SIZE][RES_SIZE],char TABLE_NAME[DIC_NAME_MAX])
{
	int i=0;
	int j=0;
	FILE *arq;
	char AUX[50];

	arq = fopen (TABLE_NAME,"r");
	if(arq == NULL)
		printf("Reading Error: Residues dictionary file not found!\n");
	else
		while(fscanf(arq,"%s",AUX)!=EOF){
		    if(AUX[0] == '>'){
		        for(j=1;j<4;j++)
		            TABLE[i][j-1] = AUX[j];
		        TABLE[i][j-1]='\0'; 
		        i++;    
		    } 
		}

	fclose(arq);
	return i;
}



int dictionary_load(dict *DIC[TABLE_SIZE],int tablesize, char dictionary_name[500])
{
	int i;
	FILE *dictionary_file;
	char AUX[50];
	dict *p;
	int flag=1;

	dictionary_file = fopen (dictionary_name,"r");
		
	if(dictionary_file == NULL)
	{
		flag = 0;
		printf("Dictionary reading error! Please select a valid dictionary filename and try again.\n");
	}

	for(i=0;i<tablesize;i++)
		DIC[i]=NULL;

	for(i=-1;fscanf(dictionary_file,"%s",AUX)!=EOF && i<tablesize;)
	{
		if(AUX[0]=='>')
			i++;
		else if(AUX[0]=='#')
		    ;	
		else
		{
			p=malloc(sizeof(dict));
			trim(AUX,' ');
			//fscanf(dictionary_file,"%lf",&p->charge);
			fscanf(dictionary_file,"%lf",&p->radius);
			strcpy(p->symbol,AUX);
			p->next=NULL;
			if(DIC[i]==NULL)
			{
				DIC[i]=malloc(sizeof(dict));
			 	DIC[i]->next=p;
			}
			else
			{
				p->next=DIC[i]->next;
				DIC[i]->next=p;
			}
		}	
	}


	fclose(dictionary_file);
	return flag;
}



double dictionary_consult_radius(char RES[4],char ATOM_TYPE[4],dict *DIC[TABLE_SIZE],int tablesize,char TABLE[TABLE_SIZE][RES_SIZE],char ATOM_SYMBOL[2]){
    
    int i,aux;
	double value=0.0;
	dict *p;
	char ATOM_TYPE_AUX[4];
    char RES_AUX[4];
    
	for(i=0;i<tablesize && value == 0.0;i++)
	{  
		if(!strcmp(RES,TABLE[i]))
		{
			for(p=DIC[i]->next;p!=NULL&&value==0.0;p=p->next)
				if(!strcmp(ATOM_TYPE,p->symbol))
					value=p->radius;

		}
	}
		
	if(!value){
		printf("Warning: Atom radius not found in dictionary:%s %s!\n",ATOM_TYPE,RES);
		strcpy(RES_AUX,"GEN");
		ATOM_TYPE_AUX[0] = ATOM_SYMBOL[0];
		ATOM_TYPE_AUX[1] = ATOM_SYMBOL[1];
		ATOM_TYPE_AUX[2] ='\0';
		if(ATOM_TYPE_AUX[0] > 90 || ATOM_TYPE_AUX[0] < 65){
		    ATOM_TYPE_AUX[0] = ATOM_TYPE_AUX[1];
		    ATOM_TYPE_AUX[1] = ATOM_TYPE_AUX[2];
		}
	    for(i=0;i<tablesize && value == 0.0;i++)
	    {
		    if(!strcmp(RES_AUX,TABLE[i]))
		    {   
			    for(p=DIC[i]->next;p!=NULL&&value==0.0;p=p->next)
			    {
			       	if(!strcmp(ATOM_TYPE_AUX,p->symbol))
			    	{	
			    	    value=p->radius;
			    	    if(value)
			    		    printf("Warning: Using generic atom %s radius value %f\n", ATOM_TYPE_AUX, p->radius);
			    	}    
                }
		    }
	    }
	    if(value > 0.0)
	        return value;
	    else
	        printf("Warning: Radius data not found for atom %s. This atom will be excluded from analysis.\n", ATOM_TYPE_AUX); 
 	}
	return value;
}

/*double dictionary_consult_charge(char RES[4],char ATOM_TYPE[4],dict *DIC[TABLE_SIZE],int tablesize,char TABLE[TABLE_SIZE][RES_SIZE],char ATOM_SYMBOL[2]){
    
    int i,aux;
	double value=0.0;
	dict *p;
	char ATOM_TYPE_AUX[4];
    char RES_AUX[4];
	for(i=0;i<tablesize && value == 0.0;i++)
	{
		if(!strcmp(RES,TABLE[i]))
		{
			for(p=DIC[i]->next;p!=NULL&&value==0.0;p=p->next)
				if(!strcmp(ATOM_TYPE,p->symbol))
					value=p->charge;

		}
	}
		
	if(!value){
		printf("Warning: Atom charge not found in dictionary:%s %s!\n",ATOM_TYPE,RES);
		strcpy(RES_AUX,"GEN");
		ATOM_TYPE_AUX[0] = ATOM_SYMBOL[0];
		ATOM_TYPE_AUX[1] = ATOM_SYMBOL[1];
		ATOM_TYPE_AUX[2] ='\0';
		if(ATOM_TYPE_AUX[0] > 90 || ATOM_TYPE_AUX[0] < 65){
		    ATOM_TYPE_AUX[0] = ATOM_TYPE_AUX[1];
		    ATOM_TYPE_AUX[1] = ATOM_TYPE_AUX[2];
		}
	    for(i=0;i<tablesize && value == 0.0;i++)
	    {
		    if(!strcmp(RES_AUX,TABLE[i]))
		    {   
			    for(p=DIC[i]->next;p!=NULL&&value==0.0;p=p->next)
			    {
			       	if(!strcmp(ATOM_TYPE_AUX,p->symbol))
			    	{	
			    	    value=p->charge;
			    		printf("Warning: Using generic %s charge value %f\n", ATOM_TYPE_AUX, p->charge);
			    	}    
                }
		    }
	    } 
	}
	return value;
}
*/

