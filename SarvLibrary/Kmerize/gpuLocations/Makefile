all: 
	nvcc -w -std=c++11 kmerLocations.cu -o kmerLocs14 -Xcompiler -fopenmp

#IO Running
runIO:
	./kmerLocs14 14 /home/wrigh338/fasta_files/ecoliFixed.fasta noError.fa

#Desktop Running
runNoErrDesktop:
	./kmerLocs14 14 ~/Development/fasta_files/ecoliFixed.fasta noError.fa

clean:
	rm -f *.fa smallErr noErr 