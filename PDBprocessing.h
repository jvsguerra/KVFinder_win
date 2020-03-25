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
    KVFinder website: http://lnbio.cnpem.br/bioinformatics/main/software/KVFinder*/

#ifndef PDBPROCESSING_H
#define PDBPROCESSING_H

#define NAME_MAX 500		      /* Maximum length of PDB filename             */

typedef struct ATOM atom;

struct ATOM
{
	double x;
	double y;
	double z;
	double radius;
	//double charge;
	int resnum;
	char chain;
	struct ATOM *next;
};

atom *v;
void prints(char S[50]);
double X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4,X5,Y5,Z5,sina,sinb,cosa,cosb;
int get_line(FILE *arq, char LINE[100]);
void extract(char LINE[100],int a,int b,char S[50]);
double convert(char S[50]);
int PDB_load(dict *DIC[TABLE_SIZE],int tablesize,char TABLE[TABLE_SIZE][RES_SIZE], char PDB_NAME[NAME_MAX],double probe, double m, double n, double o, double h);
int PDB_load2(dict *DIC[TABLE_SIZE],int tablesize,char TABLE[TABLE_SIZE][RES_SIZE], char PDB_NAME[NAME_MAX]);
void insert_atom(double x,double y, double z, double radius,int resnumber,char chain);
void atom_list_print();
void free_atom();
char extractChain(char LINE[100]);
#endif
