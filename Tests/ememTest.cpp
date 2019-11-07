


int main(int argc, char * argv[])
{
	int32_t i=0, n=1;
    uint32_t options=0, revComplement=0;
    seqFileReadInfo RefFile, QueryFile;

    checkCommandLineOptions(options);

    tmpFilesInfo arrayTmpFile(IS_MATCH_BOTH_DEF(options)?(2*NUM_TMP_FILES+2):NUM_TMP_FILES+2);
    arrayTmpFile.openFiles(ios::out|ios::binary, IS_MATCH_BOTH_DEF(options)?(2*NUM_TMP_FILES+2):NUM_TMP_FILES+2);

    RefFile.generateRevComplement(0); // This routine also computers size and num sequences
    QueryFile.generateRevComplement((IS_MATCH_REV_DEF(options) || IS_MATCH_BOTH_DEF(options))); // Reverse complement only for query

    if (IS_MATCH_REV_DEF(options)){
        QueryFile.setReverseFile();
        SET_MATCH_REV(revComplement);
    }
    arrayTmpFile.setNumMemsInFile(QueryFile.allocBinArray(), QueryFile.getNumSequences());
    RefFile.allocBinArray();
    RefFile.clearFileFlag();

    while (true){
        for (i=0; i<commonData::d; i++) {
            if(RefFile.readChunks()){
                processReference(RefFile, QueryFile, arrayTmpFile, revComplement);
                RefFile.setCurrPos();
                RefFile.clearMapForNs();
            }
            else
                break;
        }

        /*
         * Process MemExt list 
         */ 

        arrayTmpFile.mergeMemExtVector(revComplement);

        if (revComplement)
            break;
        if (IS_MATCH_BOTH_DEF(options)){
            SET_MATCH_BOTH(revComplement);
            //revComplement=1;
            RefFile.clearFileFlag();
            RefFile.resetCurrPos();
            RefFile.totalBases=0;
            QueryFile.setReverseFile();
            QueryFile.totalBases=0;
        }
        else
            break;
    }

    /*
     * Free up the allocated arrays
     */
    arrayTmpFile.closeFiles(IS_MATCH_BOTH_DEF(options)?(2*NUM_TMP_FILES):NUM_TMP_FILES);
    RefFile.destroy();
    QueryFile.destroy();

    /* 
     * Populate sequence information in vectors. Use this to get MEM
     * positions relative to the original sequences.
     */
    vector<seqData> refSeqInfo;
    vector<seqData> querySeqInfo;
    refSeqInfo.reserve(RefFile.getNumSequences());
    querySeqInfo.reserve(QueryFile.getNumSequences());
    RefFile.generateSeqPos(refSeqInfo);
    QueryFile.generateSeqPos(querySeqInfo);
    RefFile.closeFile();
    QueryFile.closeFile();

    arrayTmpFile.removeDuplicates(refSeqInfo, querySeqInfo, revComplement);
    fflush(0);
    return 0;












	E_MEM e_mem;

	system_init("e-mem", argc, argv, (void *) (&e_mem));

	k_merize("e-mem", (void *) (&e_mem));
	index_generation("e-mem", (void *) (&e_mem));
	lookup("e-mem", (void *) (&e_mem));


	return EXIT_SUCCESS;
}