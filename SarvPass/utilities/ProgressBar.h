#ifndef CODE_SRC_PROGRESSBAR_H_
#define CODE_SRC_PROGRESSBAR_H_

namespace SarvOpts {
    class ProgressBar{
    private:
        unsigned num; // number of items or iterations
        unsigned ith; //
    public:
        ProgressBar(unsigned n) : num(n), ith(0) {};
        void printProgress();
        void printDone();
    };
}

#endif /* CODE_SRC_PROGRESSBAR_H_ */
