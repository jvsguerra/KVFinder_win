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


#ifndef DICTIONARY_H
#define DICTIONARY_H

#define DIC_NAME_MAX 500
#define TABLE_SIZE 500
#define RES_SIZE 4

typedef struct DICTIONARY dict;

struct DICTIONARY
{
	//double charge;
	double radius;
	char symbol[6];
	struct DICTIONARY *next;
};


void trim(char S[50], char c);
int dictionary_load(dict *DIC[TABLE_SIZE],int tablesize, char dictionary_name[500]);
double dictionary_consult_radius(char RES[4],char ATOM_TYPE[4],dict *DIC[TABLE_SIZE],int tablesize,char TABLE[TABLE_SIZE][RES_SIZE],char ATOM_SYMBOL[2]);
//double dictionary_consult_charge(char RES[4],char ATOM_TYPE[4],dict *DIC[TABLE_SIZE],int tablesize,char TABLE[TABLE_SIZE][RES_SIZE],char ATOM_SYMBOL[2]);
int define_table(char TABLE[TABLE_SIZE][RES_SIZE],char TABLE_NAME[DIC_NAME_MAX]);

#endif
