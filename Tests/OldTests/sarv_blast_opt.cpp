get_file_contents(string:path):string{
    return "ACGTACGT...";
};
strncpy(string:query_word,string:s3,long:W):long;
int main(){
    vector<long> diag_last_hit1;
    diag_last_hit1.reserve(50000000);
    vector<long> diag_last_hit2;
    diag_last_hit2.reserve(50000000);
    
    vector <[long,long]> lookups1s3;
    vector <[long,long]> lookups2s3;
    
    multimap <string,long> kmer_s1;//THIS MIGHT NEED CHANGED
    multimap <string,long> kmer_s2;//THIS MIGHT NEED CHANGED
    
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
    
    dna s2 = get_file_contents("/home/min/a/kmahadik/compilertrans/data/chr_1oran.fa");
    dna s3 = get_file_contents("/home/min/a/kmahadik/compilertrans/data/chr_1chimp.fa");
    dna s1 = get_file_contents("/home/min/a/kmahadik/compilertrans/data/chr_1human.fa");
    
    long act_S1=s1.length();
    long act_S2=s2.length();
    long act_S3=s3.length();
    
    string query_word;
    long interval;
    long status;
    
    kmer_s1 = k_merize2(s1,act_S1,W,1);

    kmer_s2 = k_merize2(s2,act_S2,W,1);
    
    long extent = act_S3-W;
    long interval=1;
    //Loop fusion

    kmer_s1 = k_merize2(s1,act_S1,W,1);
    kmer_s2 = k_merize2(s2,act_S2,W,1);
    dna str1;
    for_range(scan in [0, s3.len] interval){
        str1.str = s3.str + scan;
        status = strncpy(query_word,s3.str,W);
        lookups1s3 = lookup(query_word,kmer_s1);
//    }
//    
//    for_range(scan in [0,s3.len] interval){
//        s3.str = s3.str + scan;
//        status = strncpy(query_word,s3.str,W);
        lookups2s3 = lookup(query_word,kmer_s2);
    }
    
    for(it1 in lookups1s3){
        q1l=it1[0];
        db1l=it1[1];
        q1r=q1l+W-1;
        db1r=db1l+W-1;
        
        if(q1l - diag_last_hit1[act_S3 - db1r + q1r]<0){  // overlap
            continue;
        }
        vt_ug1 = ungapped(q1l,db1l,s1,s3,act_S1,act_S3);
        
        q1l = vt_ug1.qleftoffset;
        q1r = vt_ug1.qrightoffset;
        db1l = vt_ug1.dbleftoffset;
        db1r = vt_ug1.dbrightoffset;
        score1 = vt_ug1.score;
        
        diag_last_hit1[q1r - db1r + act_S3] = q1r;
    }
    for(it2 in lookups2s3) {
        q2l=it2[0];
        db2l=it2[1];
        q2r=q2l+W-1;
        db2r=db2l+W-1;
        
        if(q2l-diag_last_hit2[act_S3-db2r+q2r]<0){
            continue;
        }
        vt_ug2 = ungapped(q2l,db2l,s2,s3,act_S2,act_S3);
        
        q2l = vt_ug2.qleftoffset;
        q2r = vt_ug2.qrightoffset;
        db2l = vt_ug2.dbleftoffset;
        db2r = vt_ug2.dbrightoffset;
        score2 = vt_ug2.score;
        
        diag_last_hit2[q2r - db2r + act_S3] = q2r;
    }
}
