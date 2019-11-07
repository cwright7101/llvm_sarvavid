#include "simularityComp.h"
#include <cstdio>

namespace sarv{
    /******HELPER FUNCTIONS FIRST!!!********/
    inline int givemax(int a,int b,int c){
        return a >= b ? (a >= c ? a : c) : (b >= c ? b : c); 
    }
    inline int evaluate(char a,char b){
        if(a==b)
            return MATCH;
        return MISMATCH;
    }
    bool extend_left_bp(int &qlefto, int &dblefto, int &maxqleft, int &maxdbleft, int &curr_Score, int &max_Score,
                        const std::string &query, const std::string &db, int type){
        bool leftext=true;
        if(type==UNGAPPED) {
            if(qlefto==0||dblefto==0)
                leftext=false;
            else {
                if (query[--qlefto] == db[--dblefto])
                    curr_Score = curr_Score + MATCH;
                else
                    curr_Score = curr_Score + MISMATCH;
                if(curr_Score > max_Score) {
                    max_Score=curr_Score;
                    maxqleft=qlefto;
                    maxdbleft=dblefto;
                }
                else if(max_Score-curr_Score >= X_DROP||curr_Score<=0) {
                    leftext=false;
                }
                else
                    leftext=true;
            }
        }
        else if (type==EXACT)
        {
            if(qlefto==0||dblefto==0)
                leftext=false;
            else {
                if (query[--qlefto] == db[--dblefto]) {
                    curr_Score = curr_Score + MATCH;
                    max_Score=curr_Score;
                    maxqleft=qlefto;
                    maxdbleft=dblefto;
                }
                else
                    leftext=false;
            }
        }
        return leftext;
    }
    bool extend_right_bp(int &qrighto, int &dbrighto, int &maxqright, int &maxdbright, int &curr_Score,
                    int &max_Score, const std::string &query, const std::string &db, int type)
    {
        int Qlength = query.length();
        int Dlength = db.length();
        bool rightext=true;
        if(type==UNGAPPED){
            if(qrighto+1==Qlength||dbrighto+1==Dlength)
                rightext=false;
            else{
                if (query[++qrighto]==db[++dbrighto])
                    curr_Score = curr_Score + MATCH;
                else
                    curr_Score = curr_Score + MISMATCH;
                if(curr_Score > max_Score){
                    max_Score=curr_Score;
                    maxqright=qrighto;
                    maxdbright=dbrighto;
                }
                else if(max_Score-curr_Score >= X_DROP||curr_Score<=0)
                    rightext=false;
                else
                    rightext=true;
            }
        }
        else if(type==EXACT){
            if (query[++qrighto]==db[++dbrighto]) {
                curr_Score = curr_Score + MATCH;
                max_Score=curr_Score;
                maxqright=qrighto;
                maxdbright=dbrighto;
            }
            else
                rightext=false;
        }
        return rightext;
    }

    /******END OF HELPER FUNCTIONS*******************************************************************/
    /************************************************************************************************/
    /**NOW THE DIFFERENT KERNEL VERSIONS*************************************************************/

    Alignment simularityCompUngappedCPU(int a, int b, const std::string &query, const std::string &database, int kmerLength){
    // Alignment simularityCompUngapped(int a,int b,const std::string &query,const std::string database,int Qlength,int Dlength, int kmerLength){
        bool leftext=true, rightext=true;
        int qleftoffset = a,dbleftoffset = b,qrightoffset = a+kmerLength-1,dbrightoffset = b+kmerLength-1,max_Score = kmerLength*MATCH,curr_Score = max_Score;
        int maxqleft=qleftoffset;
        int maxdbleft=dbleftoffset;
        int maxqright=qrightoffset;
        int maxdbright=dbrightoffset;
        while(leftext)
            leftext=extend_left_bp(qleftoffset,dbleftoffset,maxqleft,maxdbleft,curr_Score,max_Score,query,database,UNGAPPED);
        curr_Score=max_Score;
        while(rightext)
            rightext=extend_right_bp(qrightoffset,dbrightoffset,maxqright,maxdbright,curr_Score,max_Score,query,database,UNGAPPED);
        //,query.length(),database.length());
        Alignment vt;
        vt.qleftoffset=maxqleft;
        vt.qrightoffset=maxqright;
        vt.dbleftoffset=maxdbleft;
        vt.dbrightoffset=maxdbright;
        vt.score=max_Score;
        return vt;
    }

    /*Inputs: An ungapped alignment, the query and database strings */
    /*Outputs: A gapped alignment */
    /*None of the inputs are changed, they are passed by references to save space*/
    // Alignment simularityCompGappedCPU(const Alignment &vt,const std::string &query,const std::string &database){// int **score2dArray){
    //     //get the string from the query and Alignment vt
    //     //find the highest-scoring group in the current Alignment vt (highest scoring 11 contiguous letters)
    //     //the middle of this region is the start point for gapped alignment. (For protein it needs to have at least 4 exact matches)
    //     //convert the Alignment vt into a gapped Alignment
    //          //this is done using xdrop dynamic programming
    //     //return the gapped Alignment
    // }

    Alignment simularityCompGappedCPU(Alignment vt,const std::string &query,const std::string &database){
        // printf("File %s, Line: %d\n", __FILE__, __LINE__);
        //Taken from Kanak's code
        int gscore[200][200];
        const char *charquery = query.c_str(), *chardb = database.c_str();
        int Qlength = query.length(); 
        int Dlength = database.length();

        int max_Score,curr_Score;

        int i,j,k,l;
        int gql,gdl,maxdbleft,maxqleft;
        int gqr,gdr,maxdbright,maxqright;
        int break_db,break_query;
        int diag1_x1,diag1_y1;
        int diag2_x1,diag2_y1;


        int flag=1;

        break_db=0;


        int score1,score2,score3,score4;
        int score_left,score_right;
        int eval;


        gql=vt.qleftoffset;
        gqr=vt.qrightoffset;
        gdl=vt.dbleftoffset;
        gdr=vt.dbrightoffset;
        max_Score=vt.score;


        k=0,l=0;


        curr_Score=max_Score;

        maxqright=gqr;
        maxdbright=gdr;

        maxdbleft=gdl;
        maxqleft=gql;
        //initialize
        int count_iter;
        //initialize
        int startq_right,startd_right;
        startq_right=Qlength;
        startd_right=Dlength;
        break_query=Qlength;
        break_db=Dlength;
        for(i=gqr;i<startq_right;i++)
        {
            gscore[i-gqr][0]=max_Score + (i-gqr)*GAP;
            curr_Score=gscore[i-gqr][0];
            if((max_Score-curr_Score)>=X_DROP_GAP||curr_Score<0)
            {
                break_query=i;
                break;
            }
        }
        diag1_x1=break_query;
        diag1_y1=gqr;
        break_query=diag1_y1-diag1_x1;

        for(j=gdr;j<startd_right;j++)
        {
            gscore[0][j-gdr]=max_Score+(j-gdr)*GAP ;
            curr_Score=gscore[0][j-gdr];
            if((max_Score-curr_Score)>=X_DROP_GAP||curr_Score<0)
            {
                break_db=j;
                break;
            }
        }
        diag2_x1=gdr;
        diag2_y1=break_db;
        break_db=break_db-gdr;
        count_iter=0;
        int iter=0;
        char q,d,t1,t2;
        int rdiff,ldiff,common_exp;
        rdiff=gqr-gdr;
        //cout<<"\nCase1 db "<<break_db<<" Query "<<break_query;

        for(i=gqr+1;i<startq_right;i++)
        {
            for(j=gdr+1;j<startd_right;j++)
            {
                common_exp=j-i+rdiff;
                if(break_db-common_exp>1)
                {
                    if (break_query-common_exp<-1)
                    {
                        count_iter++;
                        score2=gscore[i-gqr-1][j-gdr]+GAP;
                        score3=gscore[i-gqr][j-gdr-1]+GAP;
                        if(score2>=score3)
                            score4=score2;
                        else
                            score4=score3;

                        if (charquery[i]==chardb[j])
                            score1=gscore[i-gqr-1][j-gdr-1]+MATCH;
                        else
                            score1=gscore[i-gqr-1][j-gdr-1]+MISMATCH;

                        if(score1>=score4)
                            curr_Score=score1;
                        else
                            curr_Score=score4;
                        gscore[i-gqr][j-gdr]=curr_Score;

                        if(max_Score<curr_Score)
                        {
                            max_Score=curr_Score;
                            maxqright=i;
                            maxdbright=j;
                        }
                        else if(max_Score-curr_Score>=X_DROP_GAP||curr_Score<0)
                        {
                            if(break_db-break_query==2)
                            {

                                flag=0;
                                break;
                            }
                            //gscore[i-gqr][j-gdr]=BREAK_SCORE;

                            if(break_db+break_query>2*common_exp)
                                break_query=common_exp;
                            else if(break_db+break_query<2*common_exp)
                                break_db=common_exp;
                        }

                    }
                    else if (break_query-common_exp==-1)
                    {
                        count_iter++;
                        score2=gscore[i-gqr-1][j-gdr]+GAP;
                        score3=BREAK_SCORE+GAP;
                        if(score2>=score3)
                            score4=score2;
                        else
                            score4=score3;

                        if (charquery[i]==chardb[j])
                            score1=gscore[i-gqr-1][j-gdr-1]+MATCH;
                        else
                            score1=gscore[i-gqr-1][j-gdr-1]+MISMATCH;
                        if(score1>=score4)
                            curr_Score=score1;
                        else
                            curr_Score=score4;
                        gscore[i-gqr][j-gdr]=curr_Score;
                        if(max_Score<curr_Score)
                        {
                            max_Score=curr_Score;
                            maxqright=i;
                            maxdbright=j;
                        }
                        else if(max_Score-curr_Score>=X_DROP_GAP||curr_Score<0)
                        {
                            if(break_db-break_query==2)
                            {

                                flag=0;
                                break;
                            }
                            //gscore[i-gqr][j-gdr]=BREAK_SCORE;
                            if(break_db+break_query>2*common_exp)
                                break_query=common_exp;
                            else if(break_db+break_query<2*common_exp)
                                break_db=common_exp;


                        }

                    }
                }
                else if(break_db-common_exp==1)
                {
                    if (break_query-common_exp<-1)
                    {

                        count_iter++;
                        score2=BREAK_SCORE+GAP;
                        score3=gscore[i-gqr][j-gdr-1]+GAP;
                        if(score2>=score3)
                            score4=score2;
                        else
                            score4=score3;

                        if (charquery[i]==chardb[j])
                            score1=gscore[i-gqr-1][j-gdr-1]+MATCH;
                        else
                            score1=gscore[i-gqr-1][j-gdr-1]+MISMATCH;

                        if(score1>=score4)
                            curr_Score=score1;
                        else
                            curr_Score=score4;
                        gscore[i-gqr][j-gdr]=curr_Score;
                        if(max_Score<curr_Score)
                        {
                            max_Score=curr_Score;
                            maxqright=i;
                            maxdbright=j;
                        }
                        else if(max_Score-curr_Score>=X_DROP_GAP||curr_Score<0)
                        {
                            if(break_db-break_query==2)
                            {

                                flag=0;
                                break;
                            }
                            if(break_db+break_query>2*common_exp)
                                break_query=common_exp;
                            else if(break_db+break_query<2*common_exp)
                                break_db=common_exp;
                        }

                    }
                    if (break_query-common_exp==-1)
                    {
                        count_iter++;
                        score2=BREAK_SCORE+GAP;
                        score3=BREAK_SCORE+GAP;
                        if(score2>=score3)
                            score4=score2;
                        else
                            score4=score3;
                        if (charquery[i]==chardb[j])
                            score1=gscore[i-gqr-1][j-gdr-1]+MATCH;
                        else
                            score1=gscore[i-gqr-1][j-gdr-1]+MISMATCH;
                        if(score1>=score4)
                            curr_Score=score1;
                        else
                            curr_Score=score4;
                        gscore[i-gqr][j-gdr]=curr_Score;

                        if(max_Score<curr_Score)
                        {
                            max_Score=curr_Score;
                            maxqright=i;
                            maxdbright=j;
                        }
                        else if(max_Score-curr_Score>=X_DROP_GAP||curr_Score<0)

                        {
                            if(break_db-break_query==2)
                            {

                                flag=0;
                                break;
                            }//gscore[i-gqr][j-gdr]=BREAK_SCORE;
                            if(break_db+break_query>2*common_exp)
                                break_query=common_exp;
                            else if(break_db+break_query<2*common_exp)
                                break_db=common_exp;
                            //break

                        }

                    }
                }
                else
                    break;
            }
            if(flag==0)
                break;
        }
        //cout<<"\n Count iter 1: "<<count_iter<<" max i: "<<i-gqr<<" max j: "<<j-gdr;
        score_right=max_Score;

        //max_Score=0;
        //left extension
        flag=1;
        for(i=0;i<gql;i++)
        {
            gscore[i][0]=max_Score + (i)*GAP;
            curr_Score=gscore[i][0];
            if(max_Score-curr_Score >= X_DROP_GAP||curr_Score<0 )
            {
                break_query=-i;
                break;
            }
        }


        for(j=0;j<gdl;j++)
        {
            gscore[0][j]=max_Score+(j)*GAP ;
            curr_Score=gscore[0][j];
            if(max_Score-curr_Score >= X_DROP_GAP||curr_Score<0)
            {
                break_db=j;
                break;
            }
        }


        //cout<<"\nCase2 db "<<break_db<<"Query "<<break_query;
        count_iter=0;
        for(i=1;i<=gql;i++)
        {
            for(j=1;j<=gdl;j++)
            {
                if((break_db-(j-i))>1)
                {
                    if ((break_query)-(j-i)<-1)
                    {

                        score2=gscore[i-1][j]+GAP;
                        score3=gscore[i][j-1]+GAP;
                        if(score2>=score3)
                            score4=score2;
                        else
                            score4=score3;
                        count_iter++;

                        if (charquery[gql-i]==chardb[gdl-j])
                            score1=gscore[i-1][j-1]+MATCH;
                        else
                            score1=gscore[i-1][j-1]+MISMATCH;

                        if(score1>=score4)
                            curr_Score=score1;
                        else
                            curr_Score=score4;
                        gscore[i][j]=curr_Score;

                        if(max_Score<curr_Score)
                        {
                            max_Score=curr_Score;
                            maxqleft=gql-i;
                            maxdbleft=gdl-j;
                        }
                        else if(max_Score-curr_Score>=X_DROP_GAP||curr_Score<0)
                        {
                            if(break_db-break_query==2)
                            {

                                flag=0;
                                break;
                            }
                            //gscore[i-gqr][j-gdr]=BREAK_SCORE;
                            if((break_db-(j-i))>-(break_query-(j-i)))
                                break_query=(j-i);
                            else if((break_db-(j-(i)))<=-(break_query-(-(i)+(j))))
                                break_db=(j-i);

                        }
                    }
                    else if ((break_query)-(j-i)==-1)
                    {
                        count_iter++;
                        score2=gscore[i-1][j]+GAP;
                        score3=BREAK_SCORE+GAP;
                        if(score2>=score3)
                            score4=score2;
                        else
                            score4=score3;

                        if (charquery[gql-i]==chardb[gdl-j])
                            score1=gscore[i-1][j-1]+MATCH;
                        else
                            score1=gscore[i-1][j-1]+MISMATCH;

                        if(score1>=score4)
                            curr_Score=score1;
                        else
                            curr_Score=score4;
                        gscore[i][j]=curr_Score;

                        if(max_Score<curr_Score)
                        {
                            max_Score=curr_Score;
                            maxqleft=gql-i;
                            maxdbleft=gdl-j;

                        }
                        else if(max_Score-curr_Score>=X_DROP_GAP||curr_Score<0)
                        {
                            if(break_db-break_query==2)
                            {

                                flag=0;
                                break;
                            }
                            //gscore[i-gqr][j-gdr]=BREAK_SCORE;
                            if((break_db-(j-i))>-(break_query-(j-i)))
                                break_query=(j-i);
                            else if((break_db-(j-(i)))<=-(break_query-(-(i)+(j))))
                                break_db=(j-i);

                        }

                    }
                }
                else if((break_db-(j-i))==1)
                {
                    if ((break_query)-(j-i)<-1)
                    {


                        count_iter++;
                        score2=BREAK_SCORE+GAP;
                        score3=gscore[i][j-1]+GAP;
                        if(score2>=score3)
                            score4=score2;
                        else
                            score4=score3;

                        if (charquery[gql-i]==chardb[gdl-j])
                            score1=gscore[i-1][j-1]+MATCH;
                        else
                            score1=gscore[i-1][j-1]+MISMATCH;

                        if(score1>=score4)
                            curr_Score=score1;
                        else
                            curr_Score=score4;


                        gscore[i][j]=curr_Score;

                        if(max_Score<curr_Score)
                        {
                            max_Score=curr_Score;
                            maxqleft=gql-i;
                            maxdbleft=gdl-j;
                        }
                        else if(max_Score-curr_Score>=X_DROP_GAP||curr_Score<0)
                        {
                            if(break_db-break_query==2)
                            {

                                flag=0;
                                break;
                            }
                            //gscore[i-gqr][j-gdr]=BREAK_SCORE;
                            if((break_db-(j-i))>-(break_query-(j-i)))
                                break_query=(j-i);
                            else if((break_db-(j-(i)))<=-(break_query-(-(i)+(j))))
                                break_db=(j-i);

                        }

                    }
                    else if ((break_query)-(j-i)==-1)
                    {
                        count_iter++;
                        score2=BREAK_SCORE+GAP;
                        score3=BREAK_SCORE+GAP;
                        if(score2>=score3)
                            score4=score2;
                        else
                            score4=score3;


                        if (charquery[gql-i]==chardb[gdl-j])
                            score1=gscore[i-1][j-1]+MATCH;
                        else
                            score1=gscore[i-1][j-1]+MISMATCH;
                        if(score1>=score4)
                            curr_Score=score1;
                        else
                            curr_Score=score4;
                        gscore[i][j]=curr_Score;

                        if(max_Score<curr_Score)
                        {
                            max_Score=curr_Score;
                            maxqleft=gql-i;
                            maxdbleft=gdl-j;
                        }
                        else if(max_Score-curr_Score>=X_DROP_GAP||curr_Score<0)
                        {
                            if(break_db-break_query==2)
                            {

                                flag=0;
                                break;
                            }
                            //gscore[i-gqr][j-gdr]=BREAK_SCORE;
                            if((break_db-(j-i))>-(break_query-(j-i)))
                                break_query=(j-i);
                            else if((break_db-(j-(i)))<=-(break_query-(-(i)+(j))))
                                break_db=(j-i);

                        }

                    }
                }
                else
                    break;


            }
            if(flag==0)
                break;
        }
        Alignment temp;
        temp.qleftoffset=maxqleft;
        temp.qrightoffset=maxqright;
        temp.dbleftoffset=maxdbleft;
        temp.dbrightoffset=maxdbright;

        temp.score=max_Score;
        return temp;
    }
}
