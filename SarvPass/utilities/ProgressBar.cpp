#include "ProgressBar.h"
#include "llvm/Support/raw_ostream.h"
// using namespace SarvOpts;
// using namespace llvm;

void SarvOpts::ProgressBar::printProgress(){
    ith++;
    double progress = 0;
    if (ith <= num)
        progress = (double)ith / (double)num;
    else
        progress = 1.0;
    
    int barWidth = 70;
    
    llvm::errs() << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            llvm::errs() << "=";
        else if (i == pos)
            llvm::errs() << ">";
        else
            llvm::errs() << " ";
    }
    int per = int(progress * 100.0);
    llvm::errs() << "] " << per << " %\r";
    if (per == 100)
        llvm::errs() << "\n";
    llvm::errs().flush();
}

void SarvOpts::ProgressBar::printDone(){
    ith = num;
    printProgress();
}
