#include "indexGen.h"

namespace sarv{
    /******CPU*******/
    // void indexGenCPU(std::multimap <std::string,int> &kmers, std::set<Index> &lookupTable){
    //  printf("In indexGenCPU\n");
    // }

    void indexGenCPU(std::multimap<std::string,int>& Q, std::multimap<std::string,int>& D,
         std::vector<int> &lookup_qoffset, std::vector<int> &lookup_doffset){
        for (std::multimap<std::string, int>::iterator ite = D.begin(); ite != D.end();++ite) {
            std::pair<std::multimap<std::string,int>::iterator, std::multimap<std::string,int>::iterator> range = Q.equal_range((*ite).first); //find matching query locations
            for(std::multimap<std::string,int>::iterator it=range.first; it!=range.second; ++it) {
                // std::cout << "\n" << (*it).second << " " << (*ite).second ;
                lookup_qoffset.push_back((*it).second);
                lookup_doffset.push_back((*ite).second);
            }
        }
    }

    // void indexGenCPU(std::map<unsigned int,std::set<int> >& Q, std::map<unsigned int,std::set<int> >& D,
    //         std::vector<int> &lookup_qoffset, std::vector<int> &lookup_doffset){
    //     //iterate over every entry in Q, match to entry in D
    //     for(auto q_it : Q){
    //         for(auto q_pos_it : q_it.second){
    //             for(auto d_pos_it : D[q_it.first]){
    //                 // std::cout<<"q_pos_it: "<<q_pos_it<<"\n";
    //                 // std::cout<<"d_pos_it: "<<d_pos_it<<"\n";
    //                 lookup_qoffset.push_back(q_pos_it);
    //                 lookup_doffset.push_back(d_pos_it);
    //             }
    //         }
    //     }
    // }


    /******GPU*******/
    void indexGenGPU(std::multimap <std::string,int> &kmers, std::set<Index> &lookupTable){
        printf("In indexGenGPU\n");
    }
    /*****CLUSTER****/
    void indexGenCLUSTER(std::multimap <std::string,int> &kmers, std::set<Index> &lookupTable){
        printf("In indexGenCLUSTER\n");
    }
}