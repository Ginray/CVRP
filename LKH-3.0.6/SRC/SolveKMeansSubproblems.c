#include "LKH.h"
#include "GeoConversion.h"

/*
 * The SolveKMeansSubproblems function attempts to improve a given tour
 * by means of a partitioning scheme based on K-means clustering.
 *
 * The overall region containing the nodes is subdivided into K clusters,
 * where K = ceil(Dimension/SubproblemSize). Each cluster together with
 * the given tour induces a subproblem consisting of all nodes in the
 * cluster and with edges fixed between nodes that are connected by tour
 * segments whose interior points do not belong to the cluster.
 *  
 * If an improvement is found, the new tour is written to TourFile.
 * The original tour is given by the SubproblemSuc references of the nodes.
 */

static void KMeansClustering(int K);
void Read_clustering_file();

double lambda = 5e-7;   //capacity 的惩罚值，值越大 惩罚越小    //2e-7
double lambda2 = 1e-8;   //demand 的惩罚值，1e-8
int max_capacity_each_car = 1700;  // 每个聚类中心最大的容量
int max_demand_each_car = 300000;  // 每个聚类中心最大的需求量
double eps = 1e-20;   //防止除零错误
int lambda3 = 1;      //当出现超出最大限制的时候 给予的更大惩罚
char Cluster_filename[100];
char Cluster_filepath[100];

void SolveKMeansSubproblems()
{
    Node *N;
    GainType GlobalBestCost, OldGlobalBestCost;
    double EntryTime = GetTime();
    int CurrentSubproblem, Subproblems;
    
    AllocateStructures();
    ReadPenalties();

    /* Compute upper bound for the original problem */
    // GlobalBestCost = 0;
    // N = FirstNode;
    // do {
    //     if (!Fixed(N, N->SubproblemSuc))
    //         GlobalBestCost += Distance(N, N->SubproblemSuc);
    //     N->Subproblem = 0;
    // }
    // while ((N = N->SubproblemSuc) != FirstNode);
    // if (TraceLevel >= 1) {
    //     if (TraceLevel >= 2)
    //         printff("\n");
    //     printff("*** K-means partitioning *** [Cost = " GainFormat "]\n",
    //             GlobalBestCost);
    // }

    KMeansClustering(CLUSTERING_K);
    // for (CurrentSubproblem = 1;
    //      CurrentSubproblem <= Subproblems; CurrentSubproblem++) {
    //     OldGlobalBestCost = GlobalBestCost;
    //     SolveSubproblem(CurrentSubproblem, Subproblems, &GlobalBestCost);
    //     if (SubproblemsCompressed && GlobalBestCost == OldGlobalBestCost)
    //         SolveCompressedSubproblem(CurrentSubproblem, Subproblems,
    //                                   &GlobalBestCost);
    // }
    // printff("\nCost = " GainFormat, GlobalBestCost);
    // if (Optimum != MINUS_INFINITY && Optimum != 0)
    //     printff(", Gap = %0.4f%%",
    //             100.0 * (GlobalBestCost - Optimum) / Optimum);
    // printff(", Time = %0.2f sec. %s\n", fabs(GetTime() - EntryTime),
    //         GlobalBestCost < Optimum ? "<" : GlobalBestCost ==
    //         Optimum ? "=" : "");
    // if (SubproblemBorders && Subproblems > 1)
    //     SolveSubproblemBorderProblems(Subproblems, &GlobalBestCost);
}

#define M Beta

/*
 * The KMeansClustering function performs K-means clustering. Each node
 * is given a unique cluster number in its Subproblem field.
 *
 * The algorithm is accellerated using ideas from
 *
 *     Dan Judd, Philip K. McKinley, Anil K. Jain:
 *     Large-Scale Parallel Data Clustering.
 *     IEEE Transactions on Pattern Analysis and
 *     Machine Intelligence 20(8):871-876 (1998)
 */


int Distance_kmeans(Node * Na, Node * Nb)
{
    double xd = Na->X - Nb->X, yd = Na->Y - Nb->Y;
    int Scale_temp = 1000; // 防止距离过小
    return (int) (Scale_temp * sqrt(xd * xd + yd * yd) + 0.5);
}

int Distance_capacity_demand(Node * N,Node *Center,int *Count,int *demand_count,int num_of_center){

    int distance1 = (int)Distance(N, &Center[num_of_center]) * Precision ;
    int distance2 = (int)1.0/(eps+lambda*(max_capacity_each_car-Count[num_of_center])) ;
    int distance3 = (int)1.0/(eps+lambda2*(max_demand_each_car-demand_count[num_of_center]));
    
    if (distance2 <0){ //当小于零的时候说明都超出最大限制的容量了，给更高的惩罚
        distance2 = (int)1.0/(eps+lambda*(lambda3)) ;
    }
    if (distance3 <0){
        distance3 = (int)1.0/(eps+lambda2*(lambda3)) ;
    }

    printf("distance = %d\n",distance1 );
    printf("capacity cost = %d\n",distance2);
    printf("demand cost = %d\n",distance3);

    return  distance1+ distance2+ distance3;
}

void Read_clustering_file()
{
    if (TraceLevel >= 1)
        printf("Clustering_file_name = %s\n",Clustering_file_name);
    FILE *Cluster_problemfile;
    if (!(Cluster_problemfile = fopen(Clustering_file_name, "r")))
        eprintf("Cannot open PROBLEM_FILE: \"%s\"", Clustering_file_name);
    char *Line,*Keyword;
    char Delimiters[] = " :=\n\t\r\f\v\xef\xbb\xbf";
    while ((Line = ReadLine(Cluster_problemfile))) {
        if (!(Keyword = strtok(Line, Delimiters)))
            continue;
        for (int i = 0; i < (int) strlen(Keyword); i++)
            Keyword[i] = (char) toupper(Keyword[i]);
        if (!strcmp(Keyword, "LAMBDA1")){
            char *Token = strtok(0, Delimiters);
            if (!Token || !sscanf(Token, "%lf", &lambda))
                eprintf("LAMBDA1: Double expected");
        }
        else if (!strcmp(Keyword, "LAMBDA2")){
            char *Token = strtok(0, Delimiters);
            if (!Token || !sscanf(Token, "%lf", &lambda2))
                eprintf("LAMBDA2: Double expected");
        }
        else if (!strcmp(Keyword, "LAMBDA3")){
            char *Token = strtok(0, Delimiters);
            if (!Token || !sscanf(Token, "%d", &lambda3))
                eprintf("LAMBDA3: Integer expected");
        }
        else if (!strcmp(Keyword, "MAX_POINT_NUM")){
            char *Token = strtok(0, Delimiters);
            if (!Token || !sscanf(Token, "%d", &max_capacity_each_car))
                eprintf("MAX_POINT_NUM: Integer expected");
        }
        else if (!strcmp(Keyword, "MAX_DEMAND")){
            char *Token = strtok(0, Delimiters);
            if (!Token || !sscanf(Token, "%d", &max_demand_each_car))
                eprintf("MAX_DEMAND: Integer expected");
        }
        else if (!strcmp(Keyword, "CLUSTER_FILENAME")){
            char *Token = strtok(0, Delimiters);
            if (!Token || !sscanf(Token, "%s", Cluster_filename))
                eprintf("CLUSTER_FILENAME: String expected");
        }
        else if (!strcmp(Keyword, "CLUSTER_FILEPATH")){
            char *Token = strtok(0, Delimiters);
            if (!Token || !sscanf(Token, "%s", Cluster_filepath))
                eprintf("CLUSTER_FILEPATH: String expected");
        }
    }
    fclose(Cluster_problemfile);
}

static void KMeansClustering(int K)
{   
    Read_clustering_file();
    Node *Center, **CenterRef, **Perm, *N, Old;
    int *Count, i, j, d, OldSubproblem;
    double *SumXc, *SumYc, *SumZc, Xc, Yc, Zc;
    int *Movement, *MMax, Max;
    int Moving = 0;
    CostFunction OldDistance = Distance;

    int *demand_count ;
    assert(Center = (Node *) calloc((K + 1), sizeof(Node)));
    assert(CenterRef = (Node **) calloc((K + 1), sizeof(Node *)));
    assert(SumXc = (double *) calloc(K + 1, sizeof(double)));
    assert(SumYc = (double *) calloc(K + 1, sizeof(double)));
    assert(SumZc = (double *) calloc(K + 1, sizeof(double)));
    assert(Count = (int *) calloc(K + 1, sizeof(int)));
    assert(Movement = (int *) calloc(K + 1, sizeof(int)));
    assert(MMax = (int *) calloc(K + 1, sizeof(int)));
    assert(Perm = (Node **) malloc(Dimension * sizeof(Node *)));
    assert(demand_count = (int *) calloc(K + 1, sizeof(int)));

    /* Pick random initial centers */

    srand((int)time(0));

    for (i = 0; i < Dimension; i++)
        Perm[i] = &NodeSet[i + 1];
    for (j = 1; j < Dimension; j++) {
        i = Random() % j;
        N = Perm[j];
        Perm[j] = Perm[i];
        Perm[i] = N;
    }
    for (i = 1; i <= K; i++) {
        Center[i].X = Perm[i - 1]->X;
        Center[i].Y = Perm[i - 1]->Y;
        Center[i].Z = Perm[i - 1]->Z;
    }

    /* Assign each node to center 0 (a ghost cluster) */
    N = FirstNode;
    do {
        N->BestPi = N->Pi;
        N->Pi = 0;
        if (WeightType == GEO || WeightType == GEO_MEEUS)
            GEO2XYZ(N->X, N->Y, &N->Xc, &N->Yc, &N->Zc);
        else if (WeightType == GEOM || WeightType == GEOM_MEEUS)
            GEOM2XYZ(N->X, N->Y, &N->Xc, &N->Yc, &N->Zc);
        else {
            N->Xc = N->X;
            N->Yc = N->Y;
            N->Zc = N->Z;
        }
        N->Cost = INT_MAX / 2;
        N->M = INT_MIN;
        N->Subproblem = N->LastV = 0;
        SumXc[0] += N->Xc;
        SumYc[0] += N->Yc;
        SumZc[0] += N->Zc;
        Count[0]++;
        demand_count[0] += N->Demand;
    } while ((N = N->Suc) != FirstNode);
    Xc = Center[0].Xc = SumXc[0] / Count[0];
    Yc = Center[0].Yc = SumYc[0] / Count[0];
    Zc = Center[0].Zc = SumZc[0] / Count[0];

    if (WeightType == GEO || WeightType == GEO_MEEUS)
        XYZ2GEO(Xc, Yc, Zc, &Center[0].X, &Center[0].Y);
    if (WeightType == GEOM || WeightType == GEOM_MEEUS)
        XYZ2GEOM(Xc, Yc, Zc, &Center[0].X, &Center[0].Y);
    else {
        Center[0].X = Xc;
        Center[0].Y = Yc;
        Center[0].Z = Zc;
    }
    if (Distance == Distance_TOR_2D)
        Distance = Distance_EUC_2D;
    else if (Distance == Distance_TOR_3D)
        Distance = Distance_EUC_3D;
    
    // change for cluster
    Distance = Distance_kmeans;

    for (int step = 0;step<20;step++) {
        /* Assign each node to its closest center */
        Moving = 0;
        do {
            N->M -= Movement[N->Subproblem] + MMax[N->Subproblem];
            if (N->M < 0) {
                OldSubproblem = N->Subproblem;
                if (OldSubproblem == 0)
                    N->Cost = N->NextCost = INT_MAX;
                else {
                    N->Cost = Distance_capacity_demand(N,Center,Count,demand_count,N->Subproblem);
                    // Distance(N, &Center[N->Subproblem]) * Precision  + 1.0/(eps+lambda*(max_capacity_each_car-Count[N->Subproblem])) + 1.0/(eps+lambda2*(max_demand_each_car-demand_count[N->Subproblem]));
                    N->NextCost = Distance_capacity_demand(N,Center,Count,demand_count,N->LastV);
                    // Distance(N, &Center[N->LastV]) * Precision + 1.0/(eps+lambda*(max_capacity_each_car-Count[N->LastV])) + 1.0/(eps+lambda2*(max_demand_each_car-demand_count[N->LastV]));
                    if (N->Cost > N->NextCost) {
                        i = N->LastV;
                        N->LastV = N->Subproblem;
                        N->Subproblem = i;
                        d = N->NextCost;
                        N->NextCost = N->Cost;
                        N->Cost = d;
                    }
                }
                for (i = 1; i <= K; i++) {
                    if (i == N->Subproblem || i == N->LastV)
                        continue;
                    d = INT_MIN;
                    if ((!c || c(N, &Center[i]) <= N->Cost) &&
                        (d = Distance_capacity_demand(N,Center,Count,demand_count,i)
                        // Distance(N, &Center[i]) * Precision+ 1.0/(eps+lambda*(max_capacity_each_car-Count[i])) + 1.0/(eps+lambda2*(max_demand_each_car-demand_count[i]))
                        ) <= N->Cost) {
                        N->NextCost = N->Cost;
                        N->Cost = d;
                        N->LastV = N->Subproblem;
                        if (d < N->NextCost)
                            N->Subproblem = i;
                    } else if (d < N->NextCost &&
                               (!c || c(N, &Center[i]) < N->NextCost) &&
                               (d != INT_MIN || (d =Distance_capacity_demand(N,Center,Count,demand_count,i)
                            //    Distance(N, &Center[i]) *Precision + 1.0/(eps+lambda*(max_capacity_each_car-Count[i])) + 1.0/(eps+lambda2*(max_demand_each_car-demand_count[i]))
                               ) <
                                N->NextCost)) {
                        N->NextCost = d;
                        N->LastV = i;
                    }
                }
                N->M = N->NextCost - N->Cost;
                if (N->Subproblem != OldSubproblem) {
                    Moving++;
                    SumXc[OldSubproblem] -= N->Xc;
                    SumYc[OldSubproblem] -= N->Yc;
                    SumZc[OldSubproblem] -= N->Zc;
                    Count[OldSubproblem]--;
                    demand_count[OldSubproblem] -= N->Demand;
                    SumXc[N->Subproblem] += N->Xc;
                    SumYc[N->Subproblem] += N->Yc;
                    SumZc[N->Subproblem] += N->Zc;
                    Count[N->Subproblem]++;
                    demand_count[N->Subproblem] +=  N->Demand;
                }
            }
        } while ((N = N->Suc) != FirstNode);
        if (!Moving)
            break;
        if (TraceLevel >= 1)
            printff("Moving %d %s\n", Moving,
                    Moving > 1 ? "points" : "point");
        /* Move centers */
        Max = INT_MIN;
        for (i = 1; i <= K; i++) {
            if (Count[i] > 0) {
                Old.X = Center[i].X;
                Old.Y = Center[i].Y;
                Old.Z = Center[i].Z;
                Xc = Center[i].Xc = SumXc[i] / Count[i];
                Yc = Center[i].Yc = SumYc[i] / Count[i];
                Zc = Center[i].Zc = SumZc[i] / Count[i];
                if (WeightType == GEO || WeightType == GEO_MEEUS)
                    XYZ2GEO(Xc, Yc, Zc, &Center[i].X, &Center[i].Y);
                else if (WeightType == GEOM || WeightType == GEOM_MEEUS)
                    XYZ2GEOM(Xc, Yc, Zc, &Center[i].X, &Center[i].Y);
                else {
                    Center[i].X = Xc;
                    Center[i].Y = Yc;
                    Center[i].Z = Zc;
                }
                Movement[i] = Distance(&Old, &Center[i]) * Precision;
                if (Movement[i] > Max)
                    Max = Movement[i];
            } else
                Movement[i] = 0;
        }
        for (i = 1; i <= K; i++) {
            if (Movement[i] != Max)
                MMax[i] = Max;
            else {
                MMax[i] = INT_MIN;
                for (j = 1; j <= K; j++)
                    if (j != i && Movement[j] > MMax[i])
                        MMax[i] = Movement[j];
            }
        }
    }
    Distance = OldDistance;
    N = FirstNode;
    do
        N->Pi = N->BestPi;
    while ((N = N->Suc) != FirstNode);
    for (i = 1; i <= K; i++)
        CenterRef[i] = &Center[i];
    
    int K_num[10];
    int demain_sum[10];
    int sum_of_point = 0;
    memset(K_num,0,sizeof(K_num));
    memset(demain_sum,0,sizeof(demain_sum));

    N = FirstNode;
    do
    {
        // printff("%d ",N->Subproblem);
        if(N->DepotId>0)   // 去除DEPORT点的影响
            continue;
        K_num[N->Subproblem] ++;
        demain_sum[N->Subproblem] += N->Demand;
        sum_of_point++;
    }while ((N = N->Suc) != FirstNode);
    printff("\n");
    for(int  temp = 1;temp<6;temp++){
        printff("capacity %d = %d\n",temp,K_num[temp]);
    }
    printf("--------- each deport demand = \n");
    for(int  temp = 1;temp<6;temp++){
        printff("demain %d = %d\n",temp,demain_sum[temp]);
        assert(demand_count[temp] == demain_sum[temp]);
    }
    for(int  temp = 1;temp<6;temp++){
        printf("------------------\n");
        printff("center X= %lf\n",Center[temp].X);
        printff("center Y= %lf\n",Center[temp].Y);
    }


    printff("sum_of_point = %d\n",sum_of_point);

    char full_file_path[200];
    strcat(full_file_path,Cluster_filepath);
    strcat(full_file_path,Cluster_filename);
    char FileName[10][200];
    FILE *SolutionFile[10];

    for(int i = 0;i<K;i++){
        char temp_path[200] ;
        sprintf(temp_path,"%s_%d.vrp",full_file_path,i);
        printf("write to %s\n",temp_path);
        strcpy(FileName[i],temp_path);
    }

    for(int i = 0;i<K;i++){
        assert(SolutionFile[i] = fopen(FileName[i], "w"));
    }
    for(int i = 0;i<K;i++){
        fprintf(SolutionFile[i], "NAME : %s_%d\n",Cluster_filename,i+1);
        fprintf(SolutionFile[i], "TYPE : CVRP\n");
        fprintf(SolutionFile[i], "COMMENT : 555.43\n");
        fprintf(SolutionFile[i], "DIMENSION : %d\n",K_num[i+1]+1);  //加上DEPORT点
        fprintf(SolutionFile[i], "CAPACITY : 7500\n");
        fprintf(SolutionFile[i], "CAPACITY_LIST : 7000,7500,8500\n");
        fprintf(SolutionFile[i], "CAPACITY_LIST_NUM : 20,20,21\n"); 
        fprintf(SolutionFile[i], "SERVICE_TIME : 0\n");
        fprintf(SolutionFile[i], "EDGE_WEIGHT_TYPE : GEO\n");
        fprintf(SolutionFile[i], "NODE_COORD_SECTION\n");
    }
    N = FirstNode;
    int temp_ID[10] = {1,1,1,1,1,1,1,1,1,1} ;
    do {
        if(N->DepotId>0)   // 去除DEPORT点的影响
            continue;
        fprintf(SolutionFile[N->Subproblem-1], "%d\t%lf\t%lf\n", temp_ID[N->Subproblem-1]++,N->X,N->Y);
    }while ((N = N->Suc) != FirstNode);

    for(int i = 0;i<K;i++){
        fprintf(SolutionFile[i], "%d\t%lf\t%lf\n",temp_ID[i],Depot->X,Depot->Y);
        fprintf(SolutionFile[i], "DEMAND_SECTION\n");
    }

    int temp_ID_K[10] = {1,1,1,1,1,1,1,1,1,1} ;
    N = FirstNode;
    do {
        if(N->DepotId>0)   // 去除DEPORT点的影响
            continue;
        fprintf(SolutionFile[N->Subproblem-1], "%d %d\n", temp_ID_K[N->Subproblem-1]++ ,N->Demand);
    }while ((N = N->Suc) != FirstNode);

    for(int i = 0;i<K;i++){
        fprintf(SolutionFile[i], "%d %d\n",temp_ID_K[i],0);
        fprintf(SolutionFile[i], "DEPOT_SECTION\n");
        fprintf(SolutionFile[i], "%d\n",temp_ID_K[i]);
        fprintf(SolutionFile[i], "-1\n");
        fprintf(SolutionFile[i], "EOF\n");
    }


    for(int i = 0;i<K;i++){
        fclose(SolutionFile[i]);
    }


    printff("write to file done\n");






    // AdjustClusters(K, CenterRef);  // 调整聚类中心的值，避免有的中心的点个数不平衡
    free(Center);
    free(CenterRef);
    free(SumXc);
    free(SumYc);
    free(SumZc);
    free(Count);
    free(Perm);
    free(Movement);
    free(MMax);
}
