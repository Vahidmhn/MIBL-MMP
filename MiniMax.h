double MinMax(double *Ref, double *Y_bar, double *Sol)
{
    IloNumVar Z(Env, 0,IloInfinity, ILOFLOAT);
    
    // Add the objective function
    IloObjective TempObj(Env);
    TempObj = IloMinimize(Env, Z);
    Model.add(TempObj); //
    
    // Add the related constraints
    IloExtractableArray MinMaxCons(Env);    
    for(int p=0;p<P;p++)
    {
        MinMaxCons.add(SenseSign*(Ref[p] - X[p])  <= Z);
        MinMaxCons.add(SenseSign*(Ref[p] - X[p])  >= 0);
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
            
    // Clear up the model
    Solver.clear();
    Model.remove(MinMaxCons);
    Model.remove(TempObj);
    TempObj.end();
    MinMaxCons.end();

    return Max;
}

