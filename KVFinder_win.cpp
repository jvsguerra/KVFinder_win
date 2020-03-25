#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <map>
#include <vector>
#include <sstream>
#include <string.h>

using namespace std;




typedef std::map <string,string > MapType;
MapType hash;
MapType::iterator it;




MapType opts;
MapType::iterator it_opt;

int readNewParameters(const char *file){

   	string line, line_old;
	ifstream np;
	
	
	np.open(file);

    getline(np, line);
    while (line.length() == 0 or line[0] != '#')
		getline(np, line);
	while(!np.eof()){
		if (line[0] == '>'){
			line_old = line;
			getline(np, line);
			while (!np.eof() and (line.length() == 0 or line[0] == '#'))
				getline(np, line);
		    if (line[0] != '>')
		        hash[line_old] = line;
			if (line_old == ">BOX_COORDINATES"){
				getline(np, line);
				while (line.length() == 0 or line[0] == '#')
					getline(np, line);
			    if (line[0] != '>')
				    hash[">BOX_COORDINATES1"] = line;
				getline(np, line);
				while (line.length() == 0 or line[0] == '#')
					getline(np, line);
			    if (line[0] != '>')
				    hash[">BOX_COORDINATES2"] = line;
				getline(np, line);
				while (line.length() == 0 or line[0] == '#')
					getline(np, line);
				if (line[0] != '>')
				    hash[">BOX_COORDINATES3"] = line;
				getline(np, line);
				while (line.length() == 0 or line[0] == '#')
					getline(np, line);
				if (line[0] != '>')
				    hash[">BOX_COORDINATES4"] = line;
				getline(np, line);
				while (line.length() == 0 or line[0] == '#')
					getline(np, line);
				if (line[0] != '>')
				    hash[">BOX_COORDINATES5"] = line;
				getline(np, line);
				while (line.length() == 0 or line[0] == '#')
					getline(np, line);
				if (line[0] != '>')
				    hash[">BOX_COORDINATES6"] = line;
				getline(np, line);
				while (line.length() == 0 or line[0] == '#')
					getline(np, line);
				if (line[0] != '>')
				    hash[">BOX_COORDINATES7"] = line;
				}
			}
		if (line[0] != '>')	
		    getline(np, line);
	}
}


int writeParam(string file){


  	ofstream p;
	p.open(file.c_str());
	p << hash[">DIC_PATH"] << endl;
	p << hash[">INPUT_PATH"] << endl;
	p << hash[">OUTPUT_PATH"] << endl;
	p << hash[">WHOLE_MODE"] << endl;
	p << hash[">GRID"] << endl;
	p << hash[">STEP_RED"] << endl;
	p << hash[">BOX_COORDINATES"] << endl;
	p << hash[">BOX_COORDINATES1"] << endl;
	p << hash[">BOX_COORDINATES2"] << endl;
	p << hash[">BOX_COORDINATES3"] << endl;
	p << hash[">BOX_COORDINATES4"] << endl;
	p << hash[">BOX_COORDINATES5"] << endl;
	p << hash[">BOX_COORDINATES6"] << endl;
	p << hash[">BOX_COORDINATES7"] << endl;
	p << hash[">PROBE_IN"] << endl;
	p << hash[">PROBE_OUT"] << endl;
	p << hash[">PROBE_OUT_ADJUSTMENT"] << endl;
	p << hash[">VOLUME_FILTER"] << endl;
	p << hash[">SURFACE_MODE"] << endl;
	p << "0" << endl;
	p << hash[">LIGAND_ADJUSTMENT"] << endl;
	p << hash[">LIGAND_FILE"] << endl;
	p << hash[">DISTANCE_CUTOFF"] << endl;
	p.close();
	
}

int writeNewParam(string file){

	ofstream p;
	p.open(file.c_str());
    p << "#This is a configuration file used by KVFinder. The configurable parameters are identified by >. All parameters must have a value.\n\n"; 
	p << "#Path for the KVFinder dictionary, with atoms information. For a example, check the default dictionary file on the KVFinder folder.\n";
	p << ">DIC_PATH" << endl;
	p << hash[">DIC_PATH"] << endl;
	p << "\n#Path for the input PDB. When using KVFinder on the command line mode, it must be specified using a command line parameter.\n";
	p << ">INPUT_PATH" << endl;
	p << hash[">INPUT_PATH"] << endl;
	p << "\n#Path and base name for the output files (A PDB and results file). On the command line mode, it will be defined automatically.\n";
	p << ">OUTPUT_PATH" << endl;
	p << hash[">OUTPUT_PATH"] << endl;
	p << "\n#Whole Protein mode, work as boolean (1 on and 0 off). Defines the search space as the whole protein.\n";
	p << ">WHOLE_MODE" << endl;
	p << hash[">WHOLE_MODE"] << endl;
	p << "\n#Defines the size between grid points. Directly affects the precision. Also has a effect on the running time.\n";	
	p << ">GRID" << endl;
	p << hash[">GRID"] << endl;
    p << "\n#If set (1), it limits the number of grid points, avoiding extremely long analysis. But may be set off (0) for detailed analysis.\n";
	p << ">STEP_RED" << endl;
	p << hash[">STEP_RED"] << endl;
	p << "\n>BOX_COORDINATES\n#Coordinates of the box vertices, that defines the search space when the whole protein mode is off.\n#The first 4 points represents the visible box on the interface, while the rest are used for internal computations only.\n#When using the command line mode, this values can be easily defined using the pymol interface and the generate parameters option.\n";
	p << ">BOX_COORDINATES" << endl;	
	p << hash[">BOX_COORDINATES"] << endl;
	p << hash[">BOX_COORDINATES1"] << endl;
	p << hash[">BOX_COORDINATES2"] << endl;
	p << hash[">BOX_COORDINATES3"] << endl;
	p << hash[">BOX_COORDINATES4"] << endl;
	p << hash[">BOX_COORDINATES5"] << endl;
	p << hash[">BOX_COORDINATES6"] << endl;
	p << hash[">BOX_COORDINATES7"] << endl;
	p << "\n#KVFinder works with a two sized probe system. A smaller probe, the Probe In, and a bigger one, the Probe Out, rolls around the protein. The points reached by the Probe In but not by the Probe Out are considered cavity points.\n#The Probe In radius.\n";
	p << ">PROBE_IN" << endl;
	p << hash[">PROBE_IN"] << endl;
	p << "\n#The Probe Out radius"<< endl;
	p << ">PROBE_OUT" << endl;
	p << hash[">PROBE_OUT"] << endl;
	p << "\n#This boolean defines if the Probe Out will be used (1) or not (0). Using only one probe might be useful for fast evaluations of areas. Volume results will be directly affected by the search space defined.\n";
	p << ">PROBE_OUT_ADJUSTMENT" << endl;
	p << hash[">PROBE_OUT_ADJUSTMENT"] << endl;
	p << "\n#Sets a filter on the KVFinder output, excluding cavities with smaller volumes than this parameter.\n";
	p << ">VOLUME_FILTER" << endl;
	p << hash[">VOLUME_FILTER"] << endl;
	p << "\n#Selects surface type to be considered, (0) Solvent Accessible Surface (SAS) or (1) Molecular Surface (VdW).\n";
	p << ">SURFACE_MODE" << endl;
	p << hash[">SURFACE_MODE"] << endl;
	p << "\n#This option is used when the user wants to limit the search space around a ligand. (1) for on and (0) for off.\n";
	p << ">LIGAND_ADJUSTMENT" << endl;
	p << hash[">LIGAND_ADJUSTMENT"] << endl;
	p << "\n#If the ligand adjustment is set, this specifies the path for the ligand file.\n";
	p << ">LIGAND_FILE" << endl;
	p << hash[">LIGAND_FILE"] << endl;
    p << "\n#Defines the seach radius for the ligand adjustmen mode.\n";
    p << ">DISTANCE_CUTOFF" << endl;	
	p << hash[">DISTANCE_CUTOFF"] << endl;
	p.close();
}

/*converts char to string*/
string convertChar(char a){
   stringstream ss;
   ss << a;
   return ss.str();
}


string cut_path(string path){
    string aux;
    int i;
    i = 0;
    i = path.length();
    while ((path[i] != '\\' and path[i]!= '/') and i > 0){
        i--;
    }
    if (path[i] == '\\' or path[i] == '/')
        i++;
    if (i == 0){
        return path;
    }
    else{
        while (i < path.size()){
            aux.append(convertChar(path[i]));
            i++;
        }
        return aux;
    }
}

int main (int argc, char **argv)
{

    if (getenv("KVFinder_PATH") == NULL){
        printf("System Variable KVFinder_PATH is not defined\n");
    return 0;
    }
    
    
    hash[">DIC_PATH"] = getenv("KVFinder_PATH");
    hash[">DIC_PATH"].append("/dictionary");

    hash[">WHOLE_MODE"] = "1";
    hash[">GRID"] = "1.0";
    hash[">STEP_RED"] = "0";
    hash[">BOX_COORDINATES"] = "0 0 0";
    hash[">BOX_COORDINATES1"] = "0 0 0";
    hash[">BOX_COORDINATES2"] = "0 0 0";
    hash[">BOX_COORDINATES3"] = "0 0 0";
    hash[">BOX_COORDINATES4"] = "0 0 0";
    hash[">BOX_COORDINATES5"] = "0 0 0";
    hash[">BOX_COORDINATES6"] = "0 0 0";
    hash[">BOX_COORDINATES7"] = "0 0 0";
    hash[">PROBE_IN"] = "1.4";
    hash[">PROBE_OUT"] = "5.0";
    hash[">PROBE_OUT_ADJUSTMENT"] = "0";
    hash[">VOLUME_FILTER"] = "5.0";
    hash[">SURFACE_MODE"] = "1";
    hash[">LIGAND_ADJUSTMENT"] = "0";
    hash[">LIGAND_FILE"] = "-";
    hash[">DISTANCE_CUTOFF"] = "0.0";



    char *lfile = NULL;
    int index;
    int c;
    string command, line, com2, aux, newp;
    ofstream results;
    int p_flag = 0;
    int l_flag = 0;
    int c_flag = 0;
    int i_flag = 0;
    int g_flag = 0;
    int v_flag = 0;
    int r_flag = 0;
    int t_flag = 0;
    command.clear();
    command.append(getenv("KVFinder_PATH"));
    command.append("\\KVFinder > KVFinder.output.tmp");   
    opterr = 0;
    
    if (argc==1){
        printf ("use -h for help\n");
        return 0;
    }

    while ((c = getopt (argc, argv, "hp:r:g:v:i:t:l:")) != -1)
        switch (c){
        case 'h':
            printf("KVFinder Help\n-i PDB File for a single run\n-l File with a list of PDBs for several runs\n-p Parameters File\n-g Define grid size\n-r Probe In size\n-t Probe Out size\n-v Volume filter\n");
            return 1;
        case 'g':
            opts["g"] = optarg;
            g_flag = 1;
            break;
        case 'v':
            opts["v"] = optarg;
            v_flag = 1;
            break;
        case 'i':
            i_flag = 1;
            opts["i"] = optarg;
            break;
        case 't':
            t_flag = 1;
            opts["t"] = optarg;
            break;    
        case 'l':
            l_flag = 1;
            opts["l"] = optarg;
            break;     
        case 'p':
            opts["p"] = optarg;
            p_flag = 1;
            break;
        case 'r':
            opts["r"] = optarg;
            r_flag = 1;
            break;
        case '?':
            if (isprint (optopt))
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else
                fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
            return 1;
        default:
            abort ();    
        }
 

    if (!i_flag and !l_flag){
        cout << "Input file not specified. Use -h for help" << endl;
        return 0;
    } 

    if (p_flag)
        readNewParameters(opts["p"].c_str());

    if (g_flag)
        hash[">GRID"] = opts["g"];
    if (v_flag)
        hash[">VOLUME_FILTER"] = opts["v"];
    if (r_flag)
        hash[">PROBE_IN"] = opts["r"];       
    if (t_flag)
        hash[">PROBE_OUT"] = opts["t"];
    
    system("md KV_Files");

        


    if (i_flag){
        cout << "Running KVFinder for "<< opts["i"] << "..." << flush;     
        hash[">INPUT_PATH"] = opts["i"];
        com2.clear();
        
        opts["i"] = cut_path(opts["i"]);
        newp.clear();
        newp.append("KV_Files\\Parameters_");
        newp.append(opts["i"].erase(opts["i"].size()-4));
        
        com2.append(opts["i"]);
        hash[">OUTPUT_PATH"] = opts["i"].append(".KVFinder.output");
        writeNewParam(newp);
        writeParam("Parameters.txt");   
        system(command.c_str());
        com2.append(".KVFinder.results.txt");
        results.open(com2.c_str());
        results << "KVFinder Input = " << hash[">INPUT_PATH"] << endl;
        results << "KVFinder Output = " << hash[">OUTPUT_PATH"] << ".pdb" << endl << endl << "#KVFinder Results:" << endl;
        results.close();
        aux.clear();
        aux.append("findstr \"Cavity\" KVFinder.output.tmp >> ");
        aux.append(com2);
        system(aux.c_str());              
        system("type KVFinder.output.tmp >> KV_Files/KVFinder.log");
        system("del KVFinder.output.tmp");
        results.open(com2.c_str(),std::fstream::app);
        results << endl << "#Interface Residues for Each Cavity: " << endl;       
        results.close();
        aux.clear();
        aux.append("type cavres.txt >> ");
        aux.append(com2);
        system(aux.c_str());
        system("del cavres.txt");    
        system("del Parameters.txt");
        cout << "done!" << endl << flush;
   }
   
   if (l_flag){
        ifstream files;
        files.open(opts["l"].c_str());
        getline(files, line);
        while (line.length() == 0)
            getline(files, line);
        cout << "Running KVFinder for "<< line << "..." << flush; 
        hash[">INPUT_PATH"] = line;
        com2.clear();
        line = cut_path(line);
        newp.clear();
        newp.append("KV_Files\\Parameters_");
        newp.append(line.erase(line.size()-4));        
        com2.append(line);
        hash[">OUTPUT_PATH"] = line.append(".KVFinder.output");
        writeNewParam(newp);
        writeParam("Parameters.txt"); 
        system(command.c_str());
        com2.append(".KVFinder.results.txt");
        results.open(com2.c_str());
        results << "KVFinder Input = " << hash[">INPUT_PATH"] << endl;
        results << "KVFinder Output = " << hash[">OUTPUT_PATH"] << ".pdb" << endl << endl << "#KVFinder Results:" << endl;
        results.close();
        aux.clear();
        aux.append("findstr \"Cavity\" KVFinder.output.tmp >> ");
        aux.append(com2);
        system(aux.c_str());
        system("type KVFinder.output.tmp >> KV_Files/KVFinder.log");
        system("del KVFinder.output.tmp");      
        results.open(com2.c_str(),std::fstream::app);  
        results << endl << "#Interface Residues for Each Cavity: " << endl;
        results.close();
        aux.clear();
        aux.append("type cavres.txt >> ");
        aux.append(com2);
        system(aux.c_str());
        system("del cavres.txt");    
        system("del Parameters.txt");
        cout << "done!" << endl << flush;
        while(!files.eof()){
            getline(files, line);
            if (line.length() > 1){
                cout << "Running KVFinder for "<< line << "..." << flush; 
                hash[">INPUT_PATH"] = line;
                com2.clear();
                line = cut_path(line);
                newp.clear();
                newp.append("KV_Files/Parameters_");
                newp.append(line.erase(line.size()-4));
                com2.append(line);
                hash[">OUTPUT_PATH"] = line.append(".KVFinder.output");
                writeNewParam(newp);
                writeParam("Parameters.txt"); 
                system(command.c_str());
                com2.append(".KVFinder.results.txt");
                results.open(com2.c_str());
                results << "KVFinder Input = " << hash[">INPUT_PATH"] << endl;
                results << "KVFinder Output = " << hash[">OUTPUT_PATH"] << ".pdb" << endl << endl << "#KVFinder Results:" << endl;
                results.close();
                aux.clear();
                aux.append("findstr Cavity KVFinder.output.tmp >> ");
                aux.append(com2);
                system(aux.c_str());      
                results.open(com2.c_str(),std::fstream::app);
                system("type KVFinder.output.tmp >> KV_Files/KVFinder.log");
                system("del KVFinder.output.tmp");  
                results << endl << "#Interface Residues for Each Cavity: " << endl;
                results.close();
                aux.clear();
                aux.append("type cavres.txt >> ");
                aux.append(com2);
                system(aux.c_str());
                system("del cavres.txt");    
                system("del Parameters.txt");
                cout << "done!" << endl << flush;
            }
        }
        
   }
}
