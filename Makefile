Main : dictionary.o Matrixprocessing.o PDBprocessing.o Analisa Main.c
	gcc -g -o KVFinder.exe dictionary.o Matrixprocessing.o PDBprocessing.o Main.c -static -lm
     
Matrixprocessing.o : Matrixprocessing.c Matrixprocessing.h
	gcc -c Matrixprocessing.c -lm

dictionary.o : dictionary.c dictionary.h 
	gcc -c dictionary.c -lm

PDBprocessing.o : PDBprocessing.c PDBprocessing.h 
	gcc -c PDBprocessing.c -lm

Analisa : 
	g++ KVFinder_win.cpp -o kvfinder_win.exe
    
clean :
	rm Matrixprocessing.o dictionary.o PDBprocessing.o
	
	

