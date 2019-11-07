all: 
	nvcc -w -std=c++11 noError.cu -o noErr -Xcompiler -fopenmp
	nvcc -w -std=c++11 kmer_less14_opt.cu -o kmer_less14 -Xcompiler -fopenmp

#IO Running
runNoErrIO:
	./noErr 14 /home/wrigh338/fasta_files/ecoliFixed.fasta noError.fa

runSmallErrIO:
	./smallErr 14 /home/wrigh338/fasta_files/ecoliFixed.fasta smallError.fa

#Desktop Running
runNoErrDesktop:
	./noErr 14 ~/Development/fasta_files/ecoliFixed.fasta noError.fa

runSmallErrDesktop:
	./smallErr 14 ~/Development/fasta_files/ecoliFixed.fasta smallError.fa

clean:
	rm -f *.fa smallErr noErr 