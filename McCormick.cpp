bool McCormick(double *Ref, double *Refp, double *Dual, double *Primal, int PieceNum)
{
    //Declaration
    IloExtractableArray McCormickCons(Env);
    IloNumVar W(Env,0,GUB, ILOFLOAT);
    IloNumVarArray z(Env,PieceNum, 0, 1, ILOBOOL);
    IloNumVarArray Y1(Env,PieceNum,0, IloInfinity, ILOFLOAT);
    IloNumVarArray Y2(Env,PieceNum,0, IloInfinity, ILOFLOAT);
    

    IloObjective TempObj(Env);
    TempObj = IloMaximize(Env,W);
    
    
//    Refp[0] = Refp[0] - EPS;
//    Refp[1] = Refp[1] - EPS;
//    Ref[0] = Ref[0] + EPS;
//    Ref[1] = Ref[1] + EPS;
    
    //General McCormick Cons
    IloExpr tempExp(Env);
    double Y1_l[PieceNum],Y2_l[PieceNum],Y1_u[PieceNum],Y2_u[PieceNum];
    for(int i = 0;i<PieceNum;i++)
    {
        Y1_l[i] = Refp[0] + i * (Ref[0]-Refp[0]) / PieceNum;
        Y2_l[i] = Refp[1];// - (i+1) * (Ref[1]-Refp[1]) / PieceNum;
        Y1_u[i] = Refp[0] + (i+1) * (Ref[0]-Refp[0]) / PieceNum;
        Y2_u[i] = Ref[1] ;//- i * (Ref[1]-Refp[1]) / PieceNum; 
        
//        cout << Y1_l[i]<< "---"<< Y2_l[i]<< "---"<< Y1_u[i]<< "---"<< Y2_u[i]<< endl;
        McCormickCons.add(Y1[i] <= Y1_u[i]*z[i]);
        McCormickCons.add(Y1[i] >= Y1_l[i]*z[i]);
        McCormickCons.add(Y2[i] <= Y2_u[i]*z[i]);
        McCormickCons.add(Y2[i] >= Y2_l[i]*z[i]);
        tempExp += z[i];
    }
    McCormickCons.add(tempExp == 1);
    tempExp.end();
    
//    IloExpr tempExp1(Env);
//    IloExpr tempExp2(Env);
    IloExpr tempExp3(Env);
    IloExpr tempExp4(Env);
    IloExpr tempExpY1(Env);
    IloExpr tempExpY2(Env);
    for(int i = 0;i<PieceNum;i++)
    {
//        tempExp1 += Y2_l[i]*Y1[i]+Y1_l[i]*Y2[i]-Y1_l[i]*Y2_l[i]*z[i];
//        tempExp2 += Y2_u[i]*Y1[i]+Y1_u[i]*Y2[i]-Y1_u[i]*Y2_u[i]*z[i];
        tempExp3 += Y2_u[i]*Y1[i]+Y1_l[i]*Y2[i]-Y2_u[i]*Y1_l[i]*z[i];
        tempExp4 += Y2_l[i]*Y1[i]+Y1_u[i]*Y2[i]-Y2_l[i]*Y1_u[i]*z[i];
        tempExpY1 += Y1[i];
        tempExpY2 += Y2[i];
    }
//    McCormickCons.add(tempExp1<= W);
//    McCormickCons.add(tempExp2<= W);
    McCormickCons.add(tempExp3>= W);
    McCormickCons.add(tempExp4>= W);
    McCormickCons.add(tempExpY1 == X[0]);
    McCormickCons.add(tempExpY2 == X[1]);
//    tempExp1.end();
//    tempExp2.end();
    tempExp3.end();
    tempExp4.end();
    tempExpY1.end();
    tempExpY2.end();
    
    
    
    
    
    //--------------------------------------------------------------------
    // Old simple McCprmick
//    McCormickCons.add(W <= Top[0]*X[1]+Down[1]*X[0]-Top[0]*Down[1]);
//    McCormickCons.add(W <= Top[1]*X[0]+Down[0]*X[1]-Top[1]*Down[0]);
//    McCormickCons.add(X[0] <= Top[0]);
//    McCormickCons.add(X[1] <= Top[1]);

    Model.add(McCormickCons);      
    Model.add(TempObj);
            
    // Add the setting of the solver here
    Solver.extract(Model);
//    Solver.exportModel("test.lp");
    
    int Stat=0;
    try{
        Stat = Solver.solve();   
        CplexNodeDual += Solver.getNnodes();
    }catch(IloException &e){
        cout << e << endl;
    }

    
    //---------------------

    if (Stat == 1){
        *Dual = Solver.getObjValue();
        *Primal = Solver.getValue(X[0])*Solver.getValue(X[1]);
        
        // Add Hypotenuse Cut
        if (HYPO_CUTS_In_McCormick){
            IloExpr term(Env);
            for(int p = 0;p<P;p++)
            {
                term += Powers[p]/Solver.getValue(X[p])*X[p]-Powers[p];           
            }
            Cuts.add(term >=0);

            Model.add(Cuts);
        }
        
//        cout << *Primal << "----" << *Dual << endl;
    }else{
//        cout << "INFEAS"<< endl;
        Model.remove(McCormickCons);
        McCormickCons.end();

        Model.remove(TempObj);
        TempObj.end();
        return false;
    }
    //Removing added objects from Model
    Model.remove(McCormickCons);
    McCormickCons.end();
    
    Model.remove(TempObj);
    TempObj.end();
    
    
    
    return true;
}