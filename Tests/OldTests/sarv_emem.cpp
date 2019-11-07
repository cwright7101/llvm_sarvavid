/*E-MEM*/

get_file_contents(string:path):string{
    return "ACGTACGT...";
};
strncpy(string:query_word,string:s3,long:W):long;
int main(){
    vector <[long,long]> lookups1s3;
    vector <[long,long]> lookups2s3;
    vector<[string,long]>kmer_s1;
    vector<[string,long]>kmer_s2;
    multimap <string,long> index_s1;//THIS MIGHT NEED CHANGED
    multimap <string,long> index_s2;//THIS MIGHT NEED CHANGED
    alignment vt_ug;
    alignment vt_ug1;
    alignment vt_ug2;
    long W;
    long q1l;
    long db1l;
    long q1r;
    long db1r;
    long q2l;
    long db2l;
    long q2r;
    long db2r;
    
    dna s2 = get_file_contents("/home/min/a/kmahadik/compilertrans/data/chr_1oran.fa");
    dna s3 = get_file_contents("/home/min/a/kmahadik/compilertrans/data/chr_1chimp.fa");
    dna s1 = get_file_contents("/home/min/a/kmahadik/compilertrans/data/chr_1human.fa");
    
    long act_S1=s1.length();
    long act_S2=s2.length();
    long act_S3=s3.length();
    
    string query_word;
    long interval = 10;
    long status;
    dna str1;
    dna str2;
    
    kmer_s1 = k_merize2(s1,act_S1,W,interval);
    index_s1 = index_generation(kmer_s1);
    kmer_s2 = k_merize2(s2,act_S2,W,interval);
    index_s2 = index_generation(kmer_s2);
    
    for_range(scan in [0, s3.len] 1){
        str1.str = s3.str + scan;
        status = strncpy(query_word,str1.str,W);
        lookups1s3 = lookup(query_word,index_s1);
    }
    
    for_range(scan in [0,s3.len] 1){
        str2.str = s3.str + scan;
        status = strncpy(query_word,str2.str,W);
        lookups2s3 = lookup(query_word,index_s2);
    }
    /*optimized loop fusion code    
     for_range(scan in [0, s3.len] 1){
     str1.str = s3.str + scan;
     status = strncpy(query_word,str1.str,W);
     lookups1s3 = lookup(query_word,index_s1);
     lookups2s3 = lookup(query_word,index_s2);
     }
     */
    
    for(it1 in lookups1s3){
        q1l=it1[0];
        db1l=it1[1];
        q1r=q1l+W-1;
        db1r=db1l+W-1;
        
        vt_ug1 = similarity_computation("exact",q1l,db1l,s1,s3,act_S1,act_S3);
    }
    for(it2 in lookups2s3) {
        q2l=it2[0];
        db2l=it2[1];
        q2r=q2l+W-1;
        db2r=db2l+W-1;

        vt_ug2 = similarity_computation("exact",q2l,db2l,s2,s3,act_S2,act_S3);
    }
}















