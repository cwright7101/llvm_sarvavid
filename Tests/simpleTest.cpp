#include "Utilities/IO.h"

#include <iostream>


int main(int argc, char* argv[]){
	std::cout<<"Hello\n";
	std::string query = sarv::fileToString("helloPassed");
	std::cout<<query<<"\n";
}