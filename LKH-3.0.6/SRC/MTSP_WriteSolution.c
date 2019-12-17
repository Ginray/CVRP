#include "LKH.h"

void MTSP_WriteSolution(char *FileName, GainType Penalty, GainType Cost)
{
    FILE *SolutionFile;
    Node *N, *NextN;
    int Size, Forward;
    char *FullFileName;
    GainType Sum;

    if (FileName == 0)
        return;
    FullFileName = FullName(FileName, Cost);
    if (TraceLevel >= 1)
        printff("Writing MTSP_SOLUTION_FILE: \"%s\" ... ", FullFileName);
    assert(SolutionFile = fopen(FullFileName, "w"));
    fprintf(SolutionFile, "%s, Cost: " GainFormat "_" GainFormat "\n",
            Name, Penalty, Cost);
    fprintf(SolutionFile, "The tours traveled by the %d salesmen are:\n",
            Salesmen);
    N = Depot;
    Forward = N->Suc->Id != N->Id + DimensionSaved;
    int temp_c = 0;
    do {
        Sum = 0;
        Size = -1;
        long long int demand_sum = 0;
        temp_c = 0;
        do {
            if(N->Capacity == Capacity)
                temp_c ++;
            else
                temp_c --;
            temp_c = temp_c<0?0:temp_c;
            if(temp_c == 0)
                Capacity = N->Capacity;    //code changed

            fprintf(SolutionFile, "%d ", N->Id <= Dim ? N->Id : Depot->Id);
            NextN = Forward ? N->Suc : N->Pred;
            Sum += C(N, NextN) - N->Pi - NextN->Pi;
            demand_sum += N->Demand ; // add code here
            Size++;
            if (NextN->Id > DimensionSaved)
                NextN = Forward ? NextN->Suc : NextN->Pred;
            N = NextN;
        } while (N->DepotId == 0);
        // code change here
        N->Capacity = Capacity;

        fprintf(SolutionFile, "%d (#%d)  Cost: %lld, Demand sum = %lld , Deport Capacity = %d\n ",
                Depot->Id, Size, Sum / Precision,demand_sum,N->Capacity);
    } while (N != Depot);
    fclose(SolutionFile);
    if (TraceLevel >= 1)
        printff("done\n");
}
