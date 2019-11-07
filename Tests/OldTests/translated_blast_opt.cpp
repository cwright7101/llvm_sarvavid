#include"kernels.hpp"
#include <tuple>
#include <map>
#include <multimap>
#include <set>
#include <multiset>
#include <vector>
#include <string>
std::string get_file_contents(std::string path);
{
    return "ACGTACGT...";
}
long strncpy(std::string query_word,std::string s3,long W);
typedef struct{
    long qleftoffset;
    long qrightoffset;
    long dbleftoffset;
    long dbrightoffset;
    int score;
}alignment;
typedef struct{
    char *str;
    int len;
}dna, protein;

int main(){
    std::vector<long>  diag_last_hit1;
    diag_last_hit1.reserve(50000000);
    std::vector<long>  diag_last_hit2;
    diag_last_hit2.reserve(50000000);
    std::vector<std::tuple<long,long> >  lookups1s3;
    std::vector<std::tuple<long,long> >  lookups2s3;
    std::multimap<std::string,long>  kmer_s1;
    std::multimap<std::string,long>  kmer_s2;
    alignment vt_ug;
    alignment vt_ug1;
    alignment vt_ug2;
    long W;
    long qleftoffset;
    long dbleftoffset;
    long qrightoffset;
    long dbrightoffset;
    long q1l;
    long db1l;
    long q1r;
    long db1r;
    long q2l;
    long db2l;
    long q2r;
    long db2r;
    long score1;
    long score2;
    dna s2 = get_file_contents("/home/min/a/kmahadik/compilertrans/data/chr_1oran.fa", &);
    dna s3 = get_file_contents("/home/min/a/kmahadik/compilertrans/data/chr_1chimp.fa", &);
    dna s1 = get_file_contents("/home/min/a/kmahadik/compilertrans/data/chr_1human.fa", &);
    long act_S1 = s1.length();
    long act_S2 = s2.length();
    long act_S3 = s3.length();
    std::string query_word;
    long interval;
    long status;
    k_merize2(s1,act_S1,W,1, &kmer_s1);
    k_merize2(s2,act_S2,W,1, &kmer_s2);
    long extent = act_S3 - W;
    long interval = 1;
    k_merize2(s1,act_S1,W,1, &kmer_s1);
    k_merize2(s2,act_S2,W,1, &kmer_s2);
    for(int scan = 0;scan<s3.len;scan+=interval){
        s3.str = s3.str + scan;
        status = strncpy(query_word,s3.str,W);
        lookup(query_word,kmer_s1, &lookups1s3);
        lookup(query_word,kmer_s2, &lookups2s3);
    }
    for(std::vector<std::tuple<long,long> > ::iterator it1=lookups1s3.begin();it1!=lookups1s3.end();++it1){
        q1l = std::get<0>(*it1);
        db1l = std::get<1>(*it1);
        q1r = q1l + W - 1;
        db1r = db1l + W - 1;
        if(q1l - diag_last_hit1[act_S3 - db1r + q1r] < 0){
            continue;
        }
        ungapped(q1l,db1l,s1,s3,act_S1,act_S3, &vt_ug1);
        q1l = vt_ug1.qleftoffset;
        q1r = vt_ug1.qrightoffset;
        db1l = vt_ug1.dbleftoffset;
        db1r = vt_ug1.dbrightoffset;
        score1 = vt_ug1.score;
        diag_last_hit1[q1r - db1r + act_S3] = q1r;
    }
    for(std::vector<std::tuple<long,long> > ::iterator it2=lookups2s3.begin();it2!=lookups2s3.end();++it2){
        q2l = std::get<0>(*it2);
        db2l = std::get<1>(*it2);
        q2r = q2l + W - 1;
        db2r = db2l + W - 1;
        if(q2l - diag_last_hit2[act_S3 - db2r + q2r] < 0){
            continue;
        }
        ungapped(q2l,db2l,s2,s3,act_S2,act_S3, &vt_ug2);
        q2l = vt_ug2.qleftoffset;
        q2r = vt_ug2.qrightoffset;
        db2l = vt_ug2.dbleftoffset;
        db2r = vt_ug2.dbrightoffset;
        score2 = vt_ug2.score;
        diag_last_hit2[q2r - db2r + act_S3] = q2r;
    }
}