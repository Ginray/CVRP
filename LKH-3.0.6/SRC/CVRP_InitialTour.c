#include "LKH.h"

/* The CVRP_InitialTour function computes an initial tour using the 
 * Clarke and Wright savings algorithm.
 */

#define ITERATIONS 1

#define Degree V
#define Size LastV
#define Load Loc

static Node *Find(Node * v);
static void Union(Node * x, Node * y);
static void MakeSets();
static void Distribute(int Constrained, double R);
static void StartDistribute(int Constrained, double R );
static void MultiDistribute(int Constrained, double R,int now_Capacity,int now_Capacity_ID);
static void SequenceDistribute();
static void MultiChange();
static int compareValue(const void *s1, const void *s2);
static int compareValue2(int s1, int s2);
static void CreateS();

typedef struct Saving {
    int Value, i, j;
} Saving;

static Saving *S;
static int Sets, SSize;

int remainCapacity[1000];  // 剩余的容量，用于判断是否要到下一辆车

int num_of_capacity = 0;
int nowCapacityNum = 0;
int meanCapacity = 0;

GainType CVRP_InitialTour()
{
    /*
        New feature here
    */
    num_of_capacity = 0;
    nowCapacityNum = 0;
    meanCapacity = 0;

    for(int i = 0 ;i<CapacityNumberSum ; i++){
        for(int j = 0 ;j<CapacityNum[i];j++){
            meanCapacity += CapacityList[i];
            OldCapacityList[num_of_capacity] = CapacityList[i];
            remainCapacity[num_of_capacity++] = CapacityList[i];
        }
    }
    meanCapacity = meanCapacity/num_of_capacity;
    qsort(OldCapacityList,num_of_capacity,sizeof(int),compareValue2);
    qsort(remainCapacity,num_of_capacity,sizeof(int),compareValue2);
    /*
        END
    */
    Node *N, *Last, *Next, *Tour;
    int s, Dim = Dimension - Salesmen + 1, it;
    GainType Cost, BestCost = PLUS_INFINITY, BestPenalty = PLUS_INFINITY;
    double EntryTime = GetTime();

    assert(!Asymmetric);
    if (!S)
        CreateS();  // 按照两个点原先与DEPORT的总距离-两个点两两之间的距离，create S，从大往小排序
    for (it = 0; it < ITERATIONS; it++) {
        printf("==========it = %d\n",it);
        for (s = 0; s < Salesmen; s++) {
            Tour = s == 0 ? Depot : &NodeSet[Dim + s];
            if (Tour == FirstNode)
                FirstNode = Tour->Suc;
            Follow(Tour, Depot);
            Tour->Dad = Tour; //类似并查集，root的dad为自己
            Tour->Prev = Tour->Next = 0;
            Tour->Degree = 0;
            Tour->Size = 1;
        }
        MakeSets();  //初始化
        // if (it > 0)   // 逐渐减少惩罚力度，保证所有点都在tour上了
            // StartDistribute(1, 0.01);
        if (Sets > Salesmen)
            StartDistribute(1, 0);
        assert(Sets <= Salesmen);
        // 强制分配
        if (Sets > Salesmen) {  // 一般 sets = salesmen，只剩下deport 没有分配   
            if (BestPenalty == 0)
                continue;
            StartDistribute(0, 0);
        }
        Sets += Salesmen;
        N = FirstNode;
        do {   //这是把N个点随机分给Salesmen个deport了?
            if (N->Degree <= 1 && (s = N->Special) > 0) {   // degree <=1 ，那么一般是端点了
                Tour = s == 1 ? Depot : &NodeSet[Dim + s - 1];
                // printf("N->capacity = %d \n ",N->Capacity);
                Union(N, Tour);
                Tour->Capacity = N->Capacity;
                // printf("Tour->capacity = %d \n ",Tour->Capacity);
                s = s == Salesmen ? 1 : s + 1;
                if (N->Degree <= 1) {                       // 如果还是degree == 1 ，说明由于demand 太大，需要单独配送
                    Tour = s == 1 ? Depot : &NodeSet[Dim + s - 1];
                    Union(N, Tour);
                    Tour->Capacity = N->Capacity;
                    // printf("Tour->capacity = %d \n ",Tour->Capacity);
                }
            }
        } while ((N = N->Suc) != FirstNode);
        while (Sets > 0) {   //Sets  = 2* salesmen - deport 和端点的联合数 ,跑两遍应该就是两段进行拼接
            if (N->Degree <= 1) {
                // printf("N->capacity = %d \n ",N->Capacity);
                for (s = 1; s <= Salesmen; s++) {
                    Tour = s == 1 ? Depot : &NodeSet[Dim + s - 1];
                    if (Tour->Degree <= 1 &&
                        (Sets == 1 || Find(N) != Find(Tour))) {
                        Union(N, Tour);
                        Tour->Capacity = N->Capacity;
                        // printf("Tour->capacity = %d \n ",Tour->Capacity);
                        if (N->Degree <= 1)
                            N = N->Pred;
                        break;
                    }
                }
            }
            N = N->Suc;
        }
        Last = N = FirstNode = Depot;
        do {      //  Suc更改成next？ 
            Next = N->Next != Last ? N->Next : N->Prev;
            Follow(N, Last);
            Last = N;
        } while ((N = Next) != FirstNode);
        N = FirstNode = Depot;
        Cost = 0;
        do
            Cost += C(N, N->Suc) - N->Pi - N->Suc->Pi;
        while ((N = N->Suc) != FirstNode);
        Cost /= Precision;
        Cost += ServiceTime * (Dim - 1);
        CurrentPenalty = PLUS_INFINITY;
        CurrentPenalty = Penalty();
        if (CurrentPenalty < BestPenalty ||
            (CurrentPenalty == BestPenalty && Cost < BestCost)) {
            N = FirstNode;
            while ((N = N->OldSuc = N->Suc) != FirstNode);   // 更新为新的tour
            BestCost = Cost;
            BestPenalty = CurrentPenalty;
        }
    }
    N = FirstNode; 
    do        // suc = old suc : 保存迭代的结果
        (N->Suc = N->OldSuc)->Pred = N;
    while ((N = N->Suc) != FirstNode);
    Cost = BestCost;
    CurrentPenalty = BestPenalty;
    if (TraceLevel >= 1) {
        printff("CVRP = ");
        if (Salesmen > 1 || ProblemType == SOP)
            printff(GainFormat"_"GainFormat, CurrentPenalty, Cost);
        else
            printff(GainFormat, Cost);
        if (Optimum != MINUS_INFINITY && Optimum != 0)
            printff(", Gap = %0.2f%%", 100.0 * (Cost - Optimum) / Optimum);
        printff(", Time = %0.2f sec.\n", fabs(GetTime() - EntryTime));
    }
    if (Run == Runs) {
        free(S);
        S = 0;
    }
    return Cost;
}

void CreateS()
{
    int Dim = Dimension - Salesmen + 1, i, j;
    Node *Ni, *Nj;
    SSize = 0;
    assert(S =
           (Saving *) malloc((Dim - 2) * (Dim - 1) / 2 * sizeof(Saving)));
    /* Compute savings */
    for (i = 1; i < Dim; i++) {
        Ni = &NodeSet[i];
        if (Ni == Depot)
            continue;
        for (j = i + 1; j <= Dim; j++) {
            Nj = &NodeSet[j];
            if (Nj == Depot || Forbidden(Ni, Nj))
                continue;
            S[SSize].Value = FixedOrCommon(Ni, Nj) ? INT_MAX :
                OldDistance(Ni, Depot) +
                OldDistance(Depot, Nj) - OldDistance(Ni, Nj);
            S[SSize].i = i;
            S[SSize].j = j;
            SSize++;
        }
    }
    /* Rank the savings in descending order */
    qsort(S, SSize, sizeof(Saving), compareValue);
}

static Node *Find(Node * v)
{
    if (v != v->Dad)
        v->Dad = Find(v->Dad);
    return v->Dad;
}

static void Union(Node * x, Node * y)
{
    Node *u = Find(x), *v = Find(y);
    if (u->Size < v->Size) {
        u->Dad = v;
        v->Size += u->Size;
        v->Load += u->Load;
        v->Cost += u->Cost + OldDistance(x, y) -
            OldDistance(x, Depot) - OldDistance(y, Depot);
    } else {
        v->Dad = u;
        u->Size += v->Size;
        u->Load += v->Load;
        u->Cost += v->Cost + OldDistance(x, y) -
            OldDistance(x, Depot) - OldDistance(y, Depot);   // 变成从u到v的路线，而不是经过DEPORT了
    }
    if (x->Degree++ == 0)
        x->Prev = y;
    else
        x->Next = y;
    if (y->Degree++ == 0)
        y->Prev = x;
    else
        y->Next = x;
    Sets--;
}

static void MakeSets()
{
    Node *N = FirstNode;
    Sets = 0;
    do {
        N->Capacity_ID = 0;
        N->Capacity = 0;
        N->Dad = N;
        N->Prev = N->Next = 0;
        N->Degree = 0;
        N->Size = 1;
        N->Load = N->Demand;
        N->Cost = 2 * OldDistance(N, Depot);
        Sets++;
    } while ((N = N->Suc) != FirstNode);
}

static void Distribute(int Constrained, double R)
{
    Node *Ni, *Nj, *u, *v;
    int i;

    for (i = 0; i < SSize && Sets > Salesmen; i++) {
        if (R > 0 && Random() % 1000 <= 1000 * R)
            continue;
        Ni = &NodeSet[S[i].i];  //S的两边节点
        Nj = &NodeSet[S[i].j];
        if (Ni->Degree < 2 && Nj->Degree < 2) {
            u = Find(Ni);
            v = Find(Nj);
            if (u == v)
                continue;
            if (!Constrained ||
                (u->Load + v->Load <= Capacity &&
                 (DistanceLimit == DBL_MAX ||
                  u->Cost + v->Cost + OldDistance(Ni, Nj) -
                  OldDistance(Ni, Depot) - OldDistance(Nj, Depot) +
                  (u->Size + v->Size) * ServiceTime <= DistanceLimit)))
                Union(Ni, Nj);
        }
    }
}

/*
    New feature here 
*/


static void StartDistribute(int Constrained, double R )
{

    printf("before set number = %d ,SSize = %d\n",Sets,SSize);
    MultiDistribute( Constrained, R ,OldCapacityList[num_of_capacity-1],0);  // 用最小的进行测试
    // SequenceDistribute();
    printf("after set number = %d \n",Sets);

}


static void SequenceDistribute()
{
    Node* Ni = FirstNode;
    int Constrained = 0;
    int st_capacity_num = 0;
    int now_capaciy = OldCapacityList[st_capacity_num];
    do {
        Node* Nj = Ni->Suc;
        if (Ni->Degree < 2 && Nj->Degree < 2) {
            Node *u = Find(Ni);
            Node *v = Find(Nj);
            if (u == v){
                printf("something seems wrong\n");
                continue;
            }

            if (Constrained ||(u->Load + v->Load <= now_capaciy &&
                 (DistanceLimit == DBL_MAX ||
                  u->Cost + v->Cost + OldDistance(Ni, Nj) -
                  OldDistance(Ni, Depot) - OldDistance(Nj, Depot) +
                  (u->Size + v->Size) * ServiceTime <= DistanceLimit))){
                        Union(Ni, Nj);
            }else{   // 超出容量限制了，需要新的车
                st_capacity_num++;
                // printf("st_capacity_num = %d\n",st_capacity_num);
                if(st_capacity_num >= num_of_capacity){   // 需要加车了
                    st_capacity_num  = num_of_capacity-1;  
                    Constrained = 1 ; //强制了
                }
                now_capaciy = OldCapacityList[st_capacity_num];
            }
        }
        else{
            printf("Ni->Degree = %d ,Nj->Degree = %d\n",Ni->Degree,Nj->Degree);
        }
    } while ((Ni = Ni->Suc) != FirstNode);
}


static void MultiDistribute(int Constrained, double R ,int now_Capacity ,int now_Capacity_ID)
{
    Node *Ni, *Nj, *u, *v;
    Capacity = now_Capacity;
    for (int ii = 0; ii < SSize && Sets > Salesmen; ii++) {
        if (R > 0 && Random() % 1000 <= 1000 * R)
            continue;
        Ni = &NodeSet[S[ii].i];  //S的两边节点
        Nj = &NodeSet[S[ii].j];
        if (Ni->Degree < 2 && Nj->Degree < 2) {
            u = Find(Ni);
            v = Find(Nj);
            if (u == v)
                continue;

            if (!Constrained ||
                (u->Load + v->Load <= Capacity &&
                 (DistanceLimit == DBL_MAX ||
                  u->Cost + v->Cost + OldDistance(Ni, Nj) -
                  OldDistance(Ni, Depot) - OldDistance(Nj, Depot) +
                  (u->Size + v->Size) * ServiceTime <= DistanceLimit))){
                        Union(Ni, Nj);
                  }
        }
    }
    MultiChange();
}

static void MultiChange()
{
    Node *Ni, *Nj, *v;
    int i;
    int st_capacity_num = 1;
    Node * t = FirstNode;
    int unset_node_num = 0;
    Node *FirsetUnset = NULL;
    Node *NextUnset = NULL;
    printf("sets = %d\n",Sets);
    Sets = Salesmen;
    do {
        Node* u = Find(t);
        if(u->Capacity_ID == 0){
            // printf("u->load  =  %d\n",u->Load);
            u->Capacity_ID = st_capacity_num ;
            t->Capacity_ID = st_capacity_num ;
            if(st_capacity_num<=num_of_capacity){
                u->Capacity = OldCapacityList[st_capacity_num-1] ;
                t->Capacity = OldCapacityList[st_capacity_num-1] ;
            }else{   // 重置这个点，为未分配的点
                t->Dad = t;
                t->Prev = t->Next = 0;
                t->Degree = 0;
                t->Size = 1;
                t->Load = t->Demand;
                t->Cost = 2 * OldDistance(t, Depot);
                unset_node_num++;
                if(FirsetUnset == NULL){
                    FirsetUnset = t;
                    NextUnset = t;
                }else{
                    NextUnset->next_unset = t;
                    NextUnset = t;
                }
                Sets++;
            }
            st_capacity_num ++;
        }else{
            t->Capacity_ID = u->Capacity_ID;
            t->Capacity = u->Capacity;

            if(u->Capacity_ID > num_of_capacity){   // 重置这个点，为未分配的点
                t->Dad = t;
                t->Prev = t->Next = 0;
                t->Degree = 0;
                t->Size = 1;
                t->Load = t->Demand;
                t->Cost = 2 * OldDistance(t, Depot);
                unset_node_num++;
                if(FirsetUnset == NULL){
                    FirsetUnset = t;
                    NextUnset = t;
                }else{
                    NextUnset->next_unset = t;
                    NextUnset = t;
                }
                Sets++;
            }
        }
    } while ((t = t->Suc) != FirstNode);
    printf("sets = %d\n",Sets);
    printf("unset_node_num++ = %d\n",unset_node_num);
    
    NextUnset = FirsetUnset;
    int Constrained = 0;
    for(int ii = 0 ;ii<2;ii++){   // 第二遍开始就强制相加
        printf("the number of MultiChange = %d\n",ii);
        if(ii == 1)
            Constrained =1;
        t = FirstNode;
        do{
            if(t->Degree < 2 && t->Capacity_ID<= num_of_capacity){  // 是一个deport的端点
                Ni = t;
                

                while(1){

                    int min_dis_with_ni  = 99999999;
                    Node* temp = FirsetUnset;
                    Node* min_node = NULL;

                    while(temp!=NULL){
                        if(temp->Capacity_ID > num_of_capacity && temp->Degree < 2){
                            int now_dis = OldDistance(Ni, temp);
                            if(now_dis<min_dis_with_ni){
                                min_dis_with_ni = now_dis;
                                min_node = temp;
                            }
                        }
                        temp = temp->next_unset;
                    }

                    if(min_node == NULL|| min_dis_with_ni  == 99999999)  // 已经全部分配
                        break;

                    Nj = min_node;
                    if (Ni->Degree < 2 && Nj->Degree < 2) {
                        int now_capaciy = Ni->Capacity;
                        Node *u = Find(Ni);
                        Node *v = Find(Nj);
                        // printf("u->load = %d , v->load = %d ,u->capacity = %d\n",u->Load,v->Load,u->Capacity);
                        if(Constrained|| (u->Load + v->Load <= now_capaciy &&
                            (DistanceLimit == DBL_MAX ||
                            u->Cost + v->Cost + OldDistance(Ni, Nj) -
                            OldDistance(Ni, Depot) - OldDistance(Nj, Depot) +
                            (u->Size + v->Size) * ServiceTime <= DistanceLimit))){
                            Union(Ni, Nj);
                            Nj ->Capacity = Ni->Capacity;
                            Nj ->Capacity_ID = Ni->Capacity_ID;
                            Ni = Nj;
                        }
                        else{
                            break;
                        }
                    }else{
                        // printf("Ni->Degree = %d ,Nj->Degree = %d\n",Ni->Degree,Nj->Degree);
                        break;
                    }
                }

                // while(NextUnset!=NULL){
                //     Nj = NextUnset;
                //     if (Ni->Degree < 2 && Nj->Degree < 2) {
                //         int now_capaciy = Ni->Capacity;
                //         Node *u = Find(Ni);
                //         Node *v = Find(Nj);
                //         // printf("u->load = %d , v->load = %d ,u->capacity = %d\n",u->Load,v->Load,u->Capacity);
                //         if(Constrained|| (u->Load + v->Load <= now_capaciy &&
                //             (DistanceLimit == DBL_MAX ||
                //             u->Cost + v->Cost + OldDistance(Ni, Nj) -
                //             OldDistance(Ni, Depot) - OldDistance(Nj, Depot) +
                //             (u->Size + v->Size) * ServiceTime <= DistanceLimit))){
                //             Union(Ni, Nj);
                //             Nj ->Capacity = Ni->Capacity;
                //             Nj ->Capacity_ID = Ni->Capacity_ID;
                //             NextUnset = NextUnset->next_unset;
                //             Ni = Nj;
                //         }
                //         else{
                //             break;
                //         }
                //     }else{
                //         // printf("Ni->Degree = %d ,Nj->Degree = %d\n",Ni->Degree,Nj->Degree);
                //         break;
                //     }
                        
                // }
                
            }
        }while ((t = t->Suc) != FirstNode);
    }
    // assert(NextUnset==NULL);

    t = FirstNode;
    int degree_1_num = 0;
    do{
        if(t->Degree <=1 )
            degree_1_num++;
    }while ((t = t->Suc) != FirstNode);
    printf("degree_1_num = %d\n",degree_1_num);
}

/*
    END
*/




static int compareValue(const void *s1, const void *s2)
{
    int v1 = ((Saving *) s1)->Value;
    int v2 = ((Saving *) s2)->Value;
    return v1 > v2 ? -1 : v1 == v2 ? 0 : 1;
}

static int compareValue2(int s1, int s2)
{
    return s2-s1;
}