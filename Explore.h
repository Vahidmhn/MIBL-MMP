bool Explore(double *Ref, double *Y_bar, double *Sol, double *Y_tilde)
{
    double Max = MinMax(Ref, Y_bar, Sol);
    if (Max == -1){
        return false;
    }

    // Time commands ******************************************
    if (TimeLim-(get_wall_time()-Time)<0){return false;}
    // ********************************************************
    
    for(int p = 0;p<P; p++)
    {
        Y_tilde[p] = Ref[p] - SenseSign * Max;
    }
    
    
    // Add Hypotenuse Cut
    if (HYPO_CUTS){
        // Version 2 from Dr. Payman's paper
        IloExpr term(Env);
        for(int p = 0;p<P;p++)
        {
            term += Powers[p]/Y_bar[p]*X[p]-Powers[p];           
        }
        Cuts.add(term >=0);

        Model.add(Cuts);
    }
    
    
    return true;
}

bool Explore(double *Ref, double *Refp, double DB, double *Y_bar, double *Sol, double *Y_tilde, int PieceNum)
{
    double Max = MinMax(Ref,Refp,DB, Y_bar, Sol, PieceNum);
    if (Max == -1){
        return false;
    }

    // Time commands ******************************************
    if (TimeLim-(get_wall_time()-Time)<0){return false;}
    // ********************************************************
    
    for(int p = 0;p<P; p++)
    {
        Y_tilde[p] = Ref[p] - SenseSign * Max;
    }
    
    
    // Add Hypotenuse Cut
    if (HYPO_CUTS){
        // Version 2 from Dr. Payman's paper
        IloExpr term(Env);
        for(int p = 0;p<P;p++)
        {
            term += Powers[p]/Y_bar[p]*X[p]-Powers[p];           
        }
        Cuts.add(term >=0);

        Model.add(Cuts);
    }
    
    return true;
}

