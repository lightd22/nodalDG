% Data Extraction and plotting function for 2d unsplit modal DG
% By Devin Light 5/1/14
% ---

function out = plot_2dadv(methname,which_test,ncfilename,res,file_out,stat,...
                          contourlvl,axis_lim,print_extrema,saveOutput)
    FS = 'FontSize';
                      
    Qname = strcat('Q',res{1});
    xname = strcat('x',res{1});
    yname = strcat('y',res{1});
    muname = strcat('mu',res{1});
    tname = 'time';
    
    out.data = nc_varget(ncfilename, Qname);
    out.x = nc_varget(ncfilename, xname);
    out.y = nc_varget(ncfilename, yname);
    out.t = nc_varget(ncfilename, tname);
    out.mu = nc_varget(ncfilename,muname);
    out.N = length(out.nodes) - 1;
    out.method = methname;
    out.test = which_test;
    
    nt = size(out.t,1);
    
    if(stat == 1 || stat == 2)
        subPlotWidth = 380; subPlotHeight = 380;
        xwidth = 2*subPlotWidth; ywidth = 1.2*subPlotHeight;
        fig = figure();
        set(gcf, 'PaperUnits', 'points');
        set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
        set(fig, 'Position', [0 0 xwidth ywidth])
        r = xwidth/ywidth;
       
    end
    if(stat==0) % Just output data
    elseif(stat==1) % Make plot at half time and final time

        nlvls = [round(nt/2) nt]; % What time levels to plot at
        nCol = 2; nRow = 1; numPlot = length(nlvls);
        pSpc = 0.08; 
        colW = (1-(nCol+1)*pSpc)/nCol; rowH = colW*r; %rowH = (1-(nRow+1)*titleH)/nRow;         
        titleH = 1-rowH*nRow-pSpc*(nRow+1);
        colX = pSpc + linspace(0,nCol*(colW+pSpc),nCol+1); colX = colX(1:end-1);
        rowY = pSpc + linspace(0,nRow*(rowH+pSpc),nRow+1); rowY = rowY(1:end-1);
      
        x = out.x; y = out.y;
        for dId = 1 : numPlot;
            rowId = ceil( dId / nCol ) ;
            colId = dId - (rowId - 1) * nCol ;
            axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
            
            tmp = squeeze(out.data(nlvls(dId),:,:));
            contourf(x,y,tmp,contourlvl); caxis(axis_lim);axis([0 1 0 1]);
            xlabel('x',FS,18); ylabel('y',FS,18);
            set(gca,'XTick',0:.2:1,'YTick',0:.2:1,'XTickLabel',[0:.2:1],'YTickLabel',[0:.2:1]);
            if(print_extrema == 1)
                ftitle = ['Max:',num2str(max(tmp(:))),' ; Min:',num2str(min(tmp(:)))];
                title(ftitle,'FontSize',12);
            end
            if( colId == nCol )
                pos = get(gca,'Position');
                opos = get(gca,'OuterPosition');
                colorbar; caxis(axis_lim);
                set(gca,'Position',pos);
            end
            
        end
        axes( 'Position', [0, 1-0.5*titleH, 1, titleH] ) ;

        titlemeth = strrep(out.test,'_',' ');
        title1 = strcat(out.method,'--',titlemeth);
        text(0.5,0.2,title1,FS,14,'FontWeight','Bold','HorizontalAlignment','Center') ;
        title2 = ['N=',num2str(out.N),', ', num2str(length(out.x)) , '\times', num2str(length(out.y)), ' DOF'];       
        text(0.5,0,title2,FS,14,'FontWeight','Bold','HorizontalAlignment','Center') ;
        set(gca,'Visible','off');

        if(saveOutput == 1)
            name = strcat(file_out, res{1},'.pdf');
            saveas(fig, name, 'pdf');
        end
        
    elseif(stat==2) % Make error plot
        contourf(out.x,out.y,squeeze(out.data(end,:,:))-squeeze(out.data(1,:,:)));
        colorbar('location','EastOutside');
    elseif(stat==3) % Make animation
        scrsz = get(0,'ScreenSize');
        fig = figure('Position',[1 scrsz(4)/2 3*scrsz(4)/4 3*scrsz(4)/4]);
        for i=1:length(out.t)
            tmp = squeeze(out.data(i,:,:)); 
            contourf(out.x,out.y,tmp,contourlvl);
            colorbar('location','EastOutside')
            axis image; caxis(axis_lim);
            ftitle = ['t=',num2str(out.t(i))];title(ftitle);
            pause(0.05);
        end
    end

end