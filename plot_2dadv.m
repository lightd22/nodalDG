% Data Extraction and plotting function for 2d unsplit modal DG
% By Devin Light 5/1/14
% ---

function out = plot_2dadv(methname,which_test,ncfilename,res,meqn)
                      
    Qname = strcat('Q',res{1});
    xname = strcat('x',res{1});
    yname = strcat('y',res{1});
    muname = strcat('mu',res{1});
    maxPolyName = 'maxPoly';
    tname = 'time';
    
    data = nc_varget(ncfilename, Qname);
    if(ndims(data) > 3)
        for m=1:meqn
            qname = ['q',num2str(m)];
            out.(qname) = squeeze(data(m,:,:,:));
        end
    else
        qname = ['q',num2str(1)];
        out.(qname) = data;
    end
    out.x = nc_varget(ncfilename, xname);
    out.y = nc_varget(ncfilename, yname);
    out.t = nc_varget(ncfilename, tname);
    out.mu = nc_varget(ncfilename,muname);
    out.N = nc_varget(ncfilename,maxPolyName);
    out.method = methname;
    out.test = which_test;

end