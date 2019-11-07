// tuple example
#include <iostream>     // std::cout
#include <tuple>        // std::tuple, std::get, std::tie, std::ignore
#include <vector>
typedef struct{
    long qleftoffset;
    long qrightoffset;
    long dbleftoffset;
    long dbrightoffset;
    int score;
}alignment;
typedef struct{
    char *s1;
    int len;
}dna, protein;
int main ()
{
    std::vector<int> diag_last_hit1;
    diag_last_hit1.reserve(50000000);
    std::vector<int> diag_last_hit2;
    diag_last_hit2.reserve(50000000);
    std::vector<std::tuple<int,int> > tup_vec;
    std::tuple<int,int>  var(1,2);
    std::tuple<int,int>  var1(10,10);
    tup_vec.push_back(var);
    tup_vec.push_back(var1);
    std::cout<<"First we have: ";
    std::cout << std::get<0>(tup_vec[0]) << ' ';
    std::cout << std::get<1>(tup_vec[0]) << '\n';
    std::cout<<"Resetting both elements to 0\n";
    std::get<0>(tup_vec[0]) = 0;
    std::get<1>(tup_vec[0]) = 0;
    std::cout<<"Now we have: ";
    std::cout << std::get<0>(tup_vec[0]) << ' ';
    std::cout << std::get<1>(tup_vec[0]) << '\n';
    std::cout<<"Setting both elements to 10\n Now we have: ";
    tup_vec[0]=var1;
    std::cout << std::get<0>(tup_vec[0]) << ' ';
    std::cout << std::get<1>(tup_vec[0]) << '\n';
    
    
    
    std::tuple<int,char> foo (10,'x');
    auto bar = std::make_tuple ("test", 3.1, 14, 'y');
    std::get<2>(bar) = 100;                                    // access element
    
    long myint; char mychar;
    
    std::tie (myint, mychar) = foo;                            // unpack elements
    std::tie (std::ignore, std::ignore, myint, mychar) = bar;  // unpack (with ignore)
    
    mychar = std::get<3>(bar);
    
    std::get<0>(foo) = std::get<2>(bar);
    std::get<1>(foo) = mychar;
    
    std::cout << "foo contains: ";
    std::cout << std::get<0>(foo) << ' ';
    std::cout << std::get<1>(foo) << '\n';
    
    
    for(std::vector<std::tuple<int,int> >::iterator it1=tup_vec.begin();it1!=tup_vec.end();++it1){
        myint = std::get<0>(*it1);
        std::cout << "myint contains: "<<myint;
//        myint = std::get<1>*it1;
    }
    
    
    return 0;
}
