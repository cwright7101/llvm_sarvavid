/*MUMmer*/
configure(string:config_type):config;
int main(){
    vector<[long,long,long]> match_list;
    vector<[long,long,long,long]> cluster_list;
    config setup = configure("MUMmer");
    dna subjectmultiseq=get_file_contents(setup);
    dna querymultiseq=get_file_contents(setup);
    suffixtree stree = index_generation(setup, subjectmultiseq);
    match_list = index_lookup(setup,subjectmultiseq,stree, querymultiseq);
    cluster_list = clustering(setup, match_list);
    status = similarity_computation(subjectmultiseq, querymultiseq, cluster_list,"outfilename");
}













