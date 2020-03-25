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

#ifndef MATRIXPROCESSING_H
#define MATRIXPROCESSING_H

typedef struct NODE node;

struct NODE
{
	double volume;
	int pos;
	struct NODE *next;
};

node *V;
double volume;
int ip, jp, kp, flagr;


double max( double a, double b);
double min( double a, double b);
int count(int ***A, int i,int j,int k);
int check_pos(int ***A, int i, int j, int k,int m, int n, int o);
int check_pos2(int ***A, int i, int j, int k, int m, int n, int o);
int check_pos3(int ***A, int i, int j, int k, int m, int n, int o);
//double getcharge(double ***M, int i, int j, int k);
void Matrix_initialize(int ***A,int m,int n, int o);
void Matrix_initialize2(double ***A,int m,int n, int o);
void Matrix_fill(int ***A,double ***M, int m, int n, int o, double h,double radius);
void Matrix_fill2(int ***A,double ***M, int m, int n, int o, double h,double probe);
void free_matrix(int ***A, int m, int n, int o);
void free_matrix2(double ***A, int m, int n, int o);
void Matrix_export(int ***A, int m, int n, int o, double h,char Output_name[NAME_MAX]);
void Matrix_export2(double ***M, int ***A, int m, int n, int o,double h,char Output_name[NAME_MAX]);
void Matrix_surf(int ***A, int m, int n, int o, double h, double radius);
void Matrix_surface(int ***A, int ***S, int m, int n, int o, double h);
void Matrix_search(int ***A, int m, int n, int o, double h,double probe);
void Matrix_adjust(int ***A, int m, int n, int o, double h, double limit);
int Area_search(int ***A,int m, int n, int o, double h,double filter, int tag);
void DFS(int ***A,int a,int b, int c, int m, int n, int o, int tag, double h);
int DFS_search(int ***A,int m, int n, int o, double h,double filter);
void Matrix_subtract(int ***S, int ***A, int m, int n, int o);
void remove_cavity(int ***A, int m, int n, int o, int tag);
void free_node();
void Matrix_filter(int ***A, int m, int n, int o, double h, double bX1, double bY1, double bZ1, double bX2, double bY2, double bZ2, double norm1);


#endif
