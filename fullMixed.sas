proc template;
    edit Stat.Mixed.CovParms;
        edit ProbZ;
             format=E12.;
        end;
    end;
    
    edit Stat.Mixed.Ftests;
        edit ProbF;
             format=E12.;
        end;
    end;

    edit Stat.Mixed.Ttests;
        edit Probt;
             format=E12.;
        end;
    end;
run;

data pheno;
    infile 'life.all.csv' dlm = ',' lrecl = 1000000;
    input line sex $ rep age block temp $ vial $;
run;

proc mixed data = pheno method = reml covtest;
    class line sex rep block temp;
    model age = sex temp sex*temp / solution;
    random line temp*line sex*line sex*temp*line rep(temp*line) sex*rep(temp*line);
run;
