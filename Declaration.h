const int nThread = 1; // Number of threads in the used solver
const double TimeLim = 3600.; // The time limit for the algorithm.
const double rEPS = 1e-6;
const double aEPS = 1e-6;
const int nExtraFormatCons = 4;
double GUB = INF;
double GLB = 0;
int SenseSign = 1;

int NodeCounter = 0;
int Iteration = 0;
int McCormickCalls = 0;
int CplexNodePrimal = 0;
int CplexNodeDual = 0;


// Enhancements
bool HYPO_CUTS = true;
bool WSM_Enhance = true;
bool McCormick_Enhance1;
bool McCormick_Enhance2;
bool HYPO_CUTS_In_McCormick;


double *BestSol;
class NODE{
public:
    double *Ref;
    double *Refp;
    double LB;
    double UB;
    double *Yp;
    bool HasYp;
    
    NODE();
    virtual ~NODE();
};
vector <NODE*> Tree;


const char *lpFileName;
int P, M, N; //, VariableDef[3];
double *Powers;
//double *A, *b, *C, *d;
double Time;
string Sense;



ifstream Data;
fstream Output;

IloEnv Env;
IloModel Model(Env);
IloCplex Solver(Model);
IloRangeArray Cons(Env);
IloObjective InitObj(Env);
IloExtractableArray Cuts(Env);
IloNumVarArray X(Env);
//IloNumVarArray Y(Env);
typedef IloArray<IloNumVarArray> D2NumVarMat;