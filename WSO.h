double WSO(double Weights[], double *Y_bar, double *Ref)
{
        
    // Add the objective function
    IloObjective TempObj(Env);
    IloExpr term(Env);
    for(int p = 0; p<P; p++)
    {
        term += Weights[p]*X[p];
    }
    TempObj = IloMaximize(Env,term);
    Model.add(TempObj); //
    
    // Add the related constraints
    IloExtractableArray MinMaxCons(Env);    
    for(int p=0;p<P;p++)
    {
        MinMaxCons.add(SenseSign*(Ref[p] - X[p])  >= EPS); // Or 0???
    }
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
        TempObj.end();
        MinMaxCons.end();
//        cout << Stat << endl;
        return -1;
//        throw runtime_error(strcat("The model is not solved to optimality in " , __func__));
    }
    
    // Get the Solution and Point and Max 

    for(int p = 0; p<P;p++)
    {
        Y_bar[p] = Solver.getValue(X[p]);
    }
    
    double Z = Solver.getObjValue();
    
    // Clear up the model
    Solver.clear();
    Model.remove(MinMaxCons);
    Model.remove(TempObj);
    TempObj.end();
    MinMaxCons.end();

    return Z;
}