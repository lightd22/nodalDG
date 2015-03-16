% Data Extraction and plotting function for 2d unsplit modal DG
% By Devin Light 5/1/14
% ---

function out = plot_2dadv(methname,which_test,ncfilename,res)
                      
    Qname = strcat('Q',res{1});
    xname = strcat('x',res{1});
    yname = strcat('y',res{1});
    muname = strcat('mu',res{1});
    maxPolyName = 'maxPoly';
    tname = 'time';
    
    out.data = nc_varget(ncfilename, Qname);
    out.x = nc_varget(ncfilename, xname);
    out.y = nc_varget(ncfilename, yname);
    out.t = nc_varget(ncfilename, tname);
    out.mu = nc_varget(ncfilename,muname);
    out.N = nc_varget(ncfilename,maxPolyName);
    out.method = methname;
    out.test = which_test;

end