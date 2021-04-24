/*
 This version is the same as version 1 with this difference
 * it reads lp file as input instead of A, b, C, d.
 */
#include <fstream>
#include <iostream>
#include <string>
#include <ilcplex/ilocplex.h>
#include <vector>        
#include <cstdlib>
#include <valarray>
#include <time.h>
#include <sys/time.h>
#include <algorithm>
#define EPS .000001
#define INF 1e16
#define EqualityEPS .000001


using namespace std;

// Declarations-------------------------------------------------------
#include "Declaration.h"

// Subroutines -------------------------------------------------------
#include "Subroutines.h"

int main(int argc, char** argv) {
    // Main specifications of the problem
    Solver.setOut(Env.getNullStream());
    Solver.setParam(IloCplex::Param::Threads,nThread);
    Solver.setParam(IloCplex::Param::MIP::Tolerances::MIPGap,0);
    Solver.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap,0);
    P = 2;//atoi(argv[1]); //
    Sense = "Max";//argv[2];  //
    SenseSign = +1;  // This must change to -1 if Sense is Min 
    lpFileName = argv[1]; //
    
    HYPO_CUTS = true;
    WSM_Enhance = false;
    McCormick_Enhance1 = true;
    McCormick_Enhance2 = true;
    HYPO_CUTS_In_McCormick = true;
    if (McCormick_Enhance1 == false){HYPO_CUTS_In_McCormick=false;}
    int McBoundPieceNum = atoi(argv[2]);// Mc1
    int McCutPieceNum = atoi(argv[3]); // Mc2
    if (McBoundPieceNum == 0){McCormick_Enhance1=false;}
    if (McCutPieceNum == 0){McCormick_Enhance2=false;}
    int McCounter1 = 0;int McCounter2 = 0;
    
    //if(Sense[1] == 'i' || Sense[1] == 'I'){SenseSign = -1;}
    
    ReadData();
    Time = get_wall_time();
    FixlpModel();
    
    //// Main Algorithm
    // Finding the extreme points -------------------------------------   
    
    BestSol = new double[N];
    double Sol[N+2];
    double YI[P];
    double *ObjVal;
    double MultiObj;
    for(int i = 0;i<P;i++)
    {
        ObjVal = new double[P];
        
        ObjVal = Lex_Optimize(i, Sol);
        // Time commands ******************************************
        if (TimeLim-(get_wall_time()-Time)<0){break;}
        // ********************************************************
        YI[i] = ObjVal[i];
        MultiObj = CalcObj(ObjVal);
        
        if(MultiObj > GLB)
        {
            GLB = MultiObj;
            copy(Sol,Sol+N, BestSol);
        }
    }
    delete ObjVal;
    
    GUB = CalcObj(YI);
    
    double Refp[P] ;
    for(int i = 0;i<P;i++)
    {
        Refp[i] = Cons[P+i].getUB();
//        cout << Refp[i] << endl;
    }
    InsertNode(YI,Refp,GLB,GUB,YI, false);
    NodeCounter++;
    
    
    if(WSM_Enhance)
    {
        double Yp[P];
        double stat = WSO(Powers, Yp, YI);
        double GeoMeanUB = GeoMeanDual(Yp);
        if (GeoMeanUB<GUB)
        {
            GUB = GeoMeanUB;
            Tree.at(0)->UB = GUB;
            Tree.at(0)->HasYp = true;
            copy(Yp,Yp+P, Tree.at(0)->Yp);
            cout<< "Geometric mean is better in the first iter!"<< endl;
        }
    }
    
    if(McCormick_Enhance1)
    {
        double McUB, McLB;
        bool stat = McCormick(YI,Refp,&McUB, &McLB,McBoundPieceNum);
        McCormickCalls++;
        if (stat == true)
        {
            if (McUB<GUB)
            {
                GUB = McUB;
                Tree.at(0)->UB = GUB;
//                cout<< "McCormick gave a better Dual"<< endl;
                McCounter1++;
            }
            if (McLB>GLB)
            {
                GLB = McLB;
                // I need to update BestSol here. Do not FORGET
                Tree.at(0)->LB = GLB;
//                cout<< "McCormick gave a better Primal"<< endl;
                McCounter2++;
            }
        }
    }
    // Main Loop --------------------------------------------------------
    int ChosenNode;
    double Y_tilde[P];
    double Y_bar[P];
    bool IsExplored;
    NODE *node;
    double ub;
    
    FindUB(&ChosenNode);
    
    while(TimeLim-(get_wall_time()-Time)>.0001 && Tree.size()>0 && GUB-GLB>aEPS && ((GUB-GLB)/GUB > rEPS))
    { 
        // Pop out the chosen node;
        node = new NODE();
        node = PopOut(ChosenNode);
                
        // Explore the node
        if(McCormick_Enhance2)
        {
            IsExplored = Explore(node->Ref, node->Refp, node->UB, Y_bar, Sol, Y_tilde,McCutPieceNum);    
        }else{
            IsExplored = Explore(node->Ref, Y_bar, Sol, Y_tilde);
        }
        Iteration++;
        // Time commands ******************************************
        if (TimeLim-(get_wall_time()-Time)<0){break;}
        // ********************************************************
        
        if (!IsExplored ){
            if (Tree.size()>0){
                GUB = FindUB(&ChosenNode);
            }
            continue;
        }
        

        // Update the primal bound
        if (CalcObj(Y_bar) > GLB) ///////////////////////
        {
            GLB = CalcObj(Y_bar);
            copy(Sol,Sol+N,BestSol);
        }
        
        
        // Branch
        for(int p=0;p<P;p++)
        {
            double NewRef[P];
            double NewRefp[P];
            double Yp[P];
            copy(node->Ref,node->Ref+P, NewRef);
            NewRef[p] = Y_tilde[p];
            
            copy(Y_tilde,Y_tilde+P,NewRefp);
            NewRefp[p] = node->Refp[p];
            
            ub = min(node->UB, CalcObj(NewRef));
            
            if (WSM_Enhance)
            {                
                double stat=0;
                if(node->HasYp && Dominates(NewRef, node->Yp))
                {                    
                    copy(node->Yp,node->Yp+P,Yp);                    
                }else{
                    stat= WSO(Powers, Yp, NewRef);
                    // Time commands ******************************************
                    if (TimeLim-(get_wall_time()-Time)<0){break;}
                    // ********************************************************
                    if (stat<0) continue;                    
                    GLB = max(CalcObj(Yp), GLB);
                }
                ub = min(ub,GeoMeanDual(Yp));
                
            }
          
            
            if(McCormick_Enhance1)
            {
                double McUB, McLB;
                bool stat = McCormick(NewRef,NewRefp, &McUB, &McLB, McBoundPieceNum);
                McCormickCalls++;
                if (stat == true){
                    if (McUB<ub)
                    {
                        ub = McUB;
//                        cout<< "McCormick gave a better Dual"<< endl;
                        McCounter1++;
                    }
                    if (McLB>GLB)
                    {
                        GLB = McLB;
                        // I need to update BestSol here. Do not FORGET
//                        cout<< "McCormick gave a better Primal"<< endl;
                        McCounter2++;
                    }
                }
            }
            
            
            if ((abs(ub-GLB)>aEPS) && ((ub-GLB)/ub>rEPS) && !DominatedByAll(NewRef) )//////////////////////
            {
                InsertNode(NewRef,NewRefp,GLB,ub, Yp,true);///////////////////////
                NodeCounter++;
            }
            
        }
//        cout << GLB << " " << GUB  << " GAP = " << ((GUB-GLB)/GUB ) << endl;
        if (!(TimeLim-(get_wall_time()-Time)>.0001 &&Tree.size()>0 && GUB-GLB>aEPS && ((GUB-GLB)/GUB > rEPS)))
        {
            break;
        }
        // update Global Dual Bound
        if (Tree.size()==0)
        {
            break;
        }
        GUB = FindUB(&ChosenNode);
        
                
    }
//    
    // Result ----------------------------------------------------------
    Time = get_wall_time()-Time;
    cout << GLB << " " << GUB <<"    Time = " << Time <<"   y1= "<< BestSol[0]<<"   y2= "<< BestSol[1] << endl;

    
    Output.open("Output.txt", ios_base::app | ios_base::out);

    if(IsFileEmpty())
    {
        Output << "Instance nCons nVars GLB GUB NodeNumber Time McDualImp McCutImp "
       "McDualPieceNum McCutPieceNum Iteration McDualCallNum CplexDualNodeNum CplexPrimalNodeNum" <<endl;
    }
    
    string Instance = lpFileName;
    Output << Instance.substr(0,Instance.size()-3)<< " " << M << " " << N << " "  << GLB  << " " << GUB << " " << NodeCounter << " " << Time<< " " << McCounter1<< " " << McCounter2 << 
            " "<<McBoundPieceNum<<" "<< McCutPieceNum << " " << Iteration<< " "<< McCormickCalls <<" "<< CplexNodeDual<< " " << CplexNodePrimal<<endl;
    Output.close();
    
    return 0;
}

void ReadData()
{
    Solver.importModel(Model, lpFileName, InitObj, X, Cons); 
    M = Cons.getSize()- nExtraFormatCons;
    N = X.getSize() - P;
    cout << M<< "  "<< N << "  " << lpFileName<< endl;
    Powers = new double[P];
    Data.open("Powers.txt");
    for(int i = 0; i<P;i++)
    {
        Data >> Powers[i];
    }
    Data.close();
}

void FixlpModel()
{
    for(int p = 0;p<P;p++)
    {
        Model.remove(Cons[p]);
    }
    Model.remove(InitObj);
 }

double *Lex_Optimize(int p, double *Sol)
{
    double *out = new double[P]; 
    IloExtractableArray Lex_Cons(Env);
    if(Sense[1] == 'a' || Sense[1] == 'A')
    {
        for(int i=0;i<P;i++)
        {  
            IloObjective TempObj(Env);
            TempObj = IloMaximize(Env,X[p]); //-------------------------*

            Model.add(TempObj);
            Solver.extract(Model);
//            Solver.exportModel("Test.lp");
            // Time commands ******************************************
            if (TimeLim-(get_wall_time()-Time)<0){return out;}
            // ********************************************************
            Solver.setParam(IloCplex::TiLim,TimeLim-(get_wall_time()-Time));
            int Stat = Solver.solve();
            // Time commands ******************************************
            if (TimeLim-(get_wall_time()-Time)<0){return out;}
            // ********************************************************
            if (Stat != 1)
            {
                //throw runtime_error(strcat("The model is not solved to optimality in " , __func__));
            }            
            
            out[p] = Solver.getObjValue();
            for(int j=0;j<N+P;j++)
            {
                Sol[j] = Solver.getValue(X[j]);
            }                    
                    
            Model.remove(Lex_Cons);
            Lex_Cons.add(X[p] >= out[p]-EPS); //-------------------------*
            Model.add(Lex_Cons);
            Model.remove(TempObj);
            TempObj.end();
            Solver.clear();
            p++;
            if(p>P-1){
                p -= P;
            }
        }
    }else if (Sense[1] == 'i' || Sense[1] == 'I'){
        for(int i=0;i<P;i++)
        {            

            IloObjective TempObj(Env);
            TempObj = IloMinimize(Env,X[p]); //-------------------------*
            
            Model.add(TempObj);
            Solver.extract(Model);
//            Solver.exportModel("Test.lp");
            // Time commands ******************************************
            if (TimeLim-(get_wall_time()-Time)<0){return out;}
            // ********************************************************
            Solver.setParam(IloCplex::TiLim,TimeLim-(get_wall_time()-Time));
            // Time commands ******************************************
            if (TimeLim-(get_wall_time()-Time)<0){return out;}
            // ********************************************************
            int Stat = Solver.solve();
            if (Stat != 1)
            {
                //throw runtime_error(strcat("The model is not solved to optimality in " , __func__));
            }
            out[p] = Solver.getObjValue();
            for(int j=0;j<N+P;j++)
            {
                Sol[j] = Solver.getValue(X[j]);
            } 
            Model.remove(Lex_Cons);
            Lex_Cons.add(X[p] <= out[p]+EPS); //-------------------------*
            Model.add(Lex_Cons);
            Model.remove(TempObj);
            TempObj.end();
            Solver.clear();
            p++;
            if(p>P-1){
                p -= P;
            }
        }
    }else{
        printf("Sense must be either max or min --> %s", Sense.c_str());
        //throw runtime_error("The sense of the objective is not recognized");
    }
    Model.remove(Lex_Cons);
    return out;
}


double FindUB()
{
    double UB = Tree.at(0)->UB;
    for(int i = 1;i<Tree.size();i++)
    {
        if (Tree.at(i)->UB > UB)
        {
            UB = Tree.at(i)->UB ; 
        }
    }
    return UB;
}

double FindUB(int *ind)
{
    double UB = Tree.at(0)->UB;
    *ind = 0;
    for(int i = 1;i<Tree.size();i++)
    {
        if (Tree.at(i)->UB > UB)
        {
            UB = Tree.at(i)->UB ; 
            *ind = i;
        }
    }
    return UB;
}

double FindLB()
{
    double LB = Tree.at(0)->LB;
    for(int i = 1;i<Tree.size();i++)
    {
        if (Tree.at(i)->LB < LB)
        {
            LB = Tree.at(i) -> LB ; 
        }
    }
    return LB;
}

double FindLB(int *Ind)
{
    *Ind = 0;
    double LB = Tree.at(0)->LB;
    for(int i = 1;i<Tree.size();i++)
    {
        if (Tree.at(i)->LB < LB)
        {
            LB = Tree.at(i) -> LB ; 
            *Ind = i;
        }
    }
    return LB;
}

void InsertNode(double *Ref,double *Refp, double LB,double UB, double *Yp, bool HasYp)
{
    NODE *node = new NODE();
    copy(Ref, Ref+P, node->Ref);
    copy(Refp, Refp+P, node->Refp);
    copy(Yp, Yp+P, node->Yp);
    node->LB = LB;
    node->UB = UB;
    node->HasYp = HasYp;
    Tree.insert(Tree.begin(), node); 
}

NODE *PopOut(int i)
{
    NODE *node = new NODE();
    copy(Tree.at(i)->Ref,Tree.at(i)->Ref+P, node->Ref);
    copy(Tree.at(i)->Refp,Tree.at(i)->Refp+P, node->Refp);
    copy(Tree.at(i)->Yp,Tree.at(i)->Yp+P, node->Yp);

    node->LB = Tree.at(i)->LB;
    node->UB = Tree.at(i)->UB;
    node->HasYp = Tree.at(i)->HasYp;
    
    Tree.at(i)->~NODE();
    Tree.erase(Tree.begin()+i);
    return node;
}

#include "MiniMax.h"

#include "Explore.h"

#include "WSO.h"

#include "McCormick.cpp"



double CalcObj(double *y)
{
    double obj=1;
    for(int p=0;p<P;p++)
    {
        obj = obj*pow(y[p],Powers[p]);
    }
    return obj;
}

bool DominatedByAll(double *Ref)
{
    bool key = false;
    for(int i=0;i<Tree.size();i++)
    {
        key = Dominates(Tree.at(i)->Ref, Ref);
        
        if (key== true){
            return true;
        }
    }
    return false;
}
bool Dominates(double *x, double *y)
{
 
    bool key = true;
    for(int p= 0;p<P;p++)
    {
        if (SenseSign * x[p] < SenseSign*y[p]){
            key = false;
            break;
        }
    }
    return key;
    
}

double GeoMeanDual(double Y_bar[])
{
    double dual = 0;
    for(int p = 0;p<P;p++)
    {
        dual += Powers[p]*Y_bar[p];
    }
    dual = dual / SUM(Powers,P);
    dual = pow(dual,SUM(Powers,P));
    
    return dual;
//    if (Tree.size()==0  )
//    {
//        GUB = dual;
//    }else if(dual < FindUB()){
//        GUB = dual;
//    }
    
}

void ShowResult()
{
//    Output.open();
//    Output << ;
//    double Ys[P];
//    for(int p = 0;p<P;p++)
//    {
//        
//    }
}



NODE:: NODE()
{
    Ref = new double[P];
    Refp = new double[P];
    LB = 0.;
    UB = 0.;
    Yp = new double[P];
    HasYp = false;
}

NODE:: ~NODE()
{
    delete [] Ref;
    delete [] Yp;
}

string Str(int x)
{
	stringstream S;
	S << x;
	return S.str();
}

double SUM(double *Array, int size)
{
    double S = 0;
    for(int i = 0;i< size;i++)
    {
        S += Array[i];
    }
   
    return S;
}

double get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

bool IsFileEmpty()
{
    return Output.tellp() == 0;
}

double MinMax(double *Ref, double *Refp, double DB, double *Y_bar, double *Sol, int PieceNum)
{
    IloNumVar Z(Env, 0,IloInfinity, ILOFLOAT);
    
    IloNumVarArray z(Env,PieceNum, 0, 1, ILOBOOL);
    IloNumVarArray Y1(Env,PieceNum,0, IloInfinity, ILOFLOAT);
    IloNumVarArray Y2(Env,PieceNum,0, IloInfinity, ILOFLOAT);
    IloNumVar w(Env, GLB,GUB, ILOFLOAT);

    // Add the objective function
    IloObjective TempObj(Env);
    TempObj = IloMinimize(Env, Z);
    Model.add(TempObj); //
    
    // Add the related constraints
    IloExtractableArray MinMaxCons(Env); 
    IloExtractableArray McCormickCons(Env);
    for(int p=0;p<P;p++)
    {
        MinMaxCons.add(SenseSign*(Ref[p] - X[p])  <= Z);
        MinMaxCons.add(SenseSign*(Ref[p] - X[p])  >= 0);
    }
    //General McCormick Cons
    IloExpr tempExp(Env);
    double Y1_l[PieceNum],Y2_l[PieceNum],Y1_u[PieceNum],Y2_u[PieceNum];
    for(int i = 0;i<PieceNum;i++)
    {
        Y1_l[i] = Refp[0] + i * (Ref[0]-Refp[0]) / PieceNum;
        Y2_l[i] = Refp[1];// - (i+1) * (Ref[1]-Refp[1]) / PieceNum;
        Y1_u[i] = Refp[0] + (i+1) * (Ref[0]-Refp[0]) / PieceNum;
        Y2_u[i] = Ref[1] ;//- i * (Ref[1]-Refp[1]) / PieceNum;       
        McCormickCons.add(Y1[i] <= Y1_u[i]*z[i]);
        McCormickCons.add(Y1[i] >= Y1_l[i]*z[i]);
        McCormickCons.add(Y2[i] <= Y2_u[i]*z[i]);
        McCormickCons.add(Y2[i] >= Y2_l[i]*z[i]);
        tempExp += z[i];
    }
    McCormickCons.add(tempExp == 1);
    tempExp.end();
    
    IloExpr tempExp1(Env);
    IloExpr tempExp2(Env);
    IloExpr tempExp3(Env);
    IloExpr tempExp4(Env);
    IloExpr tempExpY1(Env);
    IloExpr tempExpY2(Env);
    for(int i = 0;i<PieceNum;i++)
    {
        tempExp1 += Y2_l[i]*Y1[i]+Y1_l[i]*Y2[i]-Y1_l[i]*Y2_l[i]*z[i];
        tempExp2 += Y2_u[i]*Y1[i]+Y1_u[i]*Y2[i]-Y1_u[i]*Y2_u[i]*z[i];
        tempExp3 += Y2_u[i]*Y1[i]+Y1_l[i]*Y2[i]-Y2_u[i]*Y1_l[i]*z[i];
        tempExp4 += Y2_l[i]*Y1[i]+Y1_u[i]*Y2[i]-Y2_l[i]*Y1_u[i]*z[i];
        tempExpY1 += Y1[i];
        tempExpY2 += Y2[i];
    }
    
    McCormickCons.add(tempExp1<= w);
    McCormickCons.add(tempExp2<= w);
    McCormickCons.add(tempExp3>= w);
    McCormickCons.add(tempExp4>= w);
    McCormickCons.add(tempExpY1 == X[0]);
    McCormickCons.add(tempExpY2 == X[1]);
    tempExp1.end();
    tempExp2.end();
    tempExp3.end();
    tempExp4.end();
    tempExpY1.end();
    tempExpY2.end();
    
    Model.add(McCormickCons);
    Model.add(MinMaxCons);
    Solver.extract(Model);
    // Time commands ******************************************
    if (TimeLim-(get_wall_time()-Time)<0){return -1;}
    // ********************************************************
//    Solver.exportModel("test.lp");
    Solver.setParam(IloCplex::TiLim,TimeLim-(get_wall_time()-Time));
    int Stat = Solver.solve();
    // Time commands ******************************************
    if (TimeLim-(get_wall_time()-Time)<0){return -1;}
    // ********************************************************
    if (Stat != 1)
    {
        Solver.clear();
        Model.remove(MinMaxCons);
        Model.remove(TempObj);
        Model.remove(McCormickCons);
        McCormickCons.end();
        TempObj.end();
        MinMaxCons.end();
        return -1;
//        throw runtime_error(strcat("The model is not solved to optimality in " , __func__));
    }
    
    // Get the Solution and Point and Max 
    for(int j=0;j<N+P;j++)
    {        
        Sol[j] = Solver.getValue(X[j]);
    }
    for(int p = 0; p<P;p++)
    {
        Y_bar[p] = Solver.getValue(X[p]);
    }
    
    double Max = Solver.getObjValue();
    CplexNodePrimal += Solver.getNnodes();
    cout << Solver.getNnodes64() << endl;
    // Clear up the model
    Solver.clear();
    Model.remove(MinMaxCons);
    Model.remove(TempObj);
    Model.remove(McCormickCons);
    McCormickCons.end();
    TempObj.end();
    MinMaxCons.end();

    return Max;
}
