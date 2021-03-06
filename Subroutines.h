void ReadData();
void FixlpModel();
double *Lex_Optimize(int p, double *Sol);
double FindUB();
double FindUB(int *Ind);
double FindLB();
double FindLB(int *Ind);
bool McCormick(double *Top, double *Down, double *Dual, double *Primal, int PieceNum);
void InsertNode(double *Ref, double *Refp, double LB, double UB, double *Yp, bool HasYp);
NODE *PopOut(int i);
double CalcObj(double *y);
double MinMax(double *Ref, double *Y_bar, double *Sol);
double MinMax(double *Ref, double *Refp, double DB, double *Y_bar, double *Sol, int PieceNum);
bool Explore(double *Ref, double *Y_bar, double *Sol, double *Y_tilde);
bool Explore(double *Ref, double *Refp, double DB, double *Y_bar, double *Sol, double *Y_tilde, int PieceNum);
bool DominatedByAll(double *Ref);
bool Dominates(double *x, double *y);
double GeoMeanDual(double Y_bar[]);
void ShowResult();
double WSO(double Weights[], double *Y_bar, double *Ref);
string Str(int x);
double SUM(double *Array, int size); 
double get_wall_time();
bool IsFileEmpty();