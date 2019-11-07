#ifndef SIMULARITYCOMP_H
#define SIMULARITYCOMP_H

#include "../Utilities/types.h"
#include <string>
#include <cstdlib>
#include <map>
#include <vector>

#define MATCH 1 
#define MISMATCH -3 
const int GAP=-5;
#define EXACT 0
#define GAPPED 1
#define UNGAPPED 2
const int X_DROP=15;
const int X_DROP_GAP=5;
const int BREAK_SCORE=-1000;

namespace sarv{
    Alignment simularityCompUngapped(int a,int b,const std::string &query,const std::string &database, int kmerLength){
        printf("Base function simularityCompUngapped");
        Alignment tmp; 
        return tmp;
    }
    Alignment simularityCompGapped(Alignment vt,const std::string &query,const std::string &database){
        printf("Base function simularityCompGapped");
        Alignment tmp; 
        return tmp;
    }
    Alignment simularityCompUngappedCPU(int a,int b,const std::string &query,const std::string &database, int kmerLength);
    Alignment simularityCompGappedCPU(Alignment vt,const std::string & query,const std::string & database);
    

    // void match_maps(std::multimap<std::string,int>&, std::multimap<std::string,int>& , std::vector <std::pair<int,int> > &);
    // bool extend_left_bp(int &qlefto,int &dblefto,int &,int &,int &curr_Score,int &max_Score,const char *,const char *,int);
    // bool extend_right_bp(int &qrighto,int &dbrighto,int &,int &,int &curr_Score,int &max_Score,const char *, const char *,int,int,int);
    // Alignment simularityCompCPU(Index hit, const std::string &query, const std::string &database, int ungapped);
    // Alignment simularityCompCPU(Alignment it, const std::string &query, const std::string &database, int gapped);
    // Alignment simularityCompGPU(Index hit, const std::string &query, const std::string &database, int ungapped);
    // Alignment simularityCompGPU(Alignment it, const std::string &query, const std::string &database, int gapped);
    // Alignment simularityCompCLUSTER(Index hit, const std::string &query, const std::string &database, int ungapped);
    // Alignment simularityCompCLUSTER(Alignment it, const std::string &query, const std::string &database, int gapped);
}
#endif