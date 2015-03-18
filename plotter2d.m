% Plot Advection Tests using plot_2dadv.m
% By: Devin Light
% ------

clear all;
close all;
clc;
%%
cd('/Users/Devin/Desktop/R/ModalDG/2d_adv/terminatorTest');

tests = {
         'def_cosinebell', ... % 1, LeVeque deformation test cosinebell                          
         'def_cyl',... % 2, Deformation flow applied to slotted cylinder
         'consistency',... %3 uniform field deformation flow
         'reactive',... % 4 Reactive half plane flow
         };
res = {'1','2','3','4'};
methods = { 'modal',...
            'modalPDt',...
            'modalPDr',...
            'modalPDs',...
          };

whichTest = tests(1);
whichRes = res(2);

ncfilename = strcat('spltMod2d_' ,whichTest{1}, '.nc');
%% Read in data
ntest = 4;
meqn = 2;
whichRes = res(2);
whichTest = tests{ntest};

subDir = '';
whichMethods = [2 3];
ncfilename = strcat('spltMod2d_' ,whichTest, '.nc');

for imethod=1:length(whichMethods)
    nmethod = whichMethods(imethod);
    methName = methods{nmethod};
    if(nmethod == 1)
        methname = 'Modal Unlimited';
        nc = ['_modal/' subDir ncfilename];
        out = plot_2dadv(methname,whichTest,nc,whichRes,meqn);
        out.figLabel = 'a';
    elseif(nmethod == 2)
        methname = 'Modal PD (trunc)';
        nc = ['_pdModal/trunc/' subDir ncfilename];
        out = plot_2dadv(methname,whichTest,nc,whichRes,meqn);
        out.figLabel = 'b';
    elseif(nmethod == 3)
        methname = 'Modal PD (rescale)';
        nc = ['_pdModal/rescale/' subDir ncfilename];
        out = plot_2dadv(methname,whichTest,nc,whichRes,meqn);
        out.figLabel = 'c';
    elseif(nmethod == 4)
        methname = 'Modal PD (strictest)';
        nc = ['_pdModal/equivscale/' subDir ncfilename];
        out = plot_2dadv(methname,whichTest,nc,whichRes,meqn);
        out.figLabel = 'd';
    end
    meth.(methName) = out;
    
    for m=1:meqn
        error = [];
        einf = [];

        qname = ['q' num2str(m)];
        ic = squeeze(out.(qname)(1,:,:));
        final = squeeze(out.(qname)(end,:,:));

        % Compute error            
        nError = sqrt(mean( (ic(:)-final(:)).^2 ));
        error = [error nError];
        errName = [qname 'l2'];
        meth.(methName).(errName) = error;

        nError = max(abs(ic(:)-final(:)));
        errName = [qname 'inf'];
        einf = [einf nError];
        meth.(methName).(errName) = einf;
    end
%}
end
%% Make Comparison Figure for this resolution
FS = 'FontSize';
cd('/Users/Devin/Desktop/R/ModalDG/2d_adv/terminatorTest/');
close all

print_errors    = 1;
print_label     = 1;
plotICs         = 1;
numRows         = 1;
figsPerRow      = 1;
saveFigure      = 0;

xwidth = 400; 
ywidth = 400;
outDir = ['_figs/_' whichTest];

for imethod=1:length(whichMethods)
    nmethod = whichMethods(imethod);
    if(ntest == 1)
        xloc1 = 0.1; xloc2 = 0.50;
        yloc1 = 0.65; yloc2 = 0.75;
        excontours = [.05 .75]; clvls = 0.1:0.1:1.0;
    elseif(ntest == 4)
        xloc1 = 0.1; xloc2 = 0.6;
        yloc1 = 0.1; yloc2 = 0.2;
        excontours = [.05 1.0]; clvls = 0.1:0.1:1.0;
    end

    methName = methods{nmethod};
    currMeth = meth.(methName);
    disp(['Plotting: ' methName]);

    x = currMeth.x;
    y = currMeth.y;
    
    for m=1:meqn
        fig = figure();
        set(gcf, 'PaperUnits', 'points');
        set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
        set(fig, 'Position', [0 0 xwidth ywidth])

        qname = ['q' num2str(m)];
        e2name = [qname 'l2'];
        einfname = [qname 'inf'];

        ics = squeeze(currMeth.(qname)(1,:,:));
        final = squeeze(currMeth.(qname)(end,:,:));
        
        final(final < 0) = final(final<0)+10^(-15);
    
        err2 = sprintf('%#.3g',currMeth.(e2name));
        einf = sprintf('%#.3g',currMeth.(einfname));
        
        hold on
        [C,h] = contour(x,y,final,clvls);
        set (h, 'LineWidth', 1,'LineColor','k');

        negclvls = -1*[1 1]*10^(-14);
        [C,h] = contour(x,y,final,negclvls,'LineWidth',0.5,'LineColor',[0.75 0.75 0.75]);
        %set(h,'LineWidth', 0.2,'LineColor',[0.95 0.95 0.95]);

        if(plotICs)
            [C,h] = contour(x,y,ics,excontours);
            set (h, 'LineWidth', 2,'LineColor','k');        
            axis square
        end
        hold off

        xlabel('x',FS,18); ylabel('y',FS,18);
        set(gca,'XTick',[0:.2:1],'YTick',[0:0.2:1]);

        if(print_errors == 1)
            text(xloc1,yloc1,['E_2= ' err2],FS,18 );
            text(xloc1,yloc2,['E_{\infty}= ' einf],FS,18 );
            text(xloc2,yloc1,sprintf('Max_ = %4.3f', max(final(:))),FS,18 );
            text(xloc2,yloc2,sprintf('Min_ = %4.3f', min(final(:))),FS,18 );
        end

        if(print_label)
%            if(ntest == 1 || nTe )
            hLu = text(0.05,0.95,[currMeth.figLabel ') ' currMeth.method '--' qname],FS,18); 
            axis([-0.005 1.005 -0.005 1.005]);
%            end
        end

        opos = get(gca,'OuterPosition');
        pos = get(gca,'Position');

        currLabel = currMeth.figLabel;
        if( (figsPerRow*numRows-(currLabel-'a')) > 2)
            xtl = '';
            set(gca,'XTickLabel',xtl,FS,8);
            xlabel('');
        end
        if( mod(currLabel-'a',figsPerRow) ~= 0)
            ytl = ''; yaxlab = '';
            set(gca,'YTickLabel',ytl,FS,8);
            ylabel('');
        end

        set(gca,FS,16,'Position',pos,'OuterPosition',opos);
        box on;

        pow = str2double(whichRes);
        nelem = length(ics)/(currMeth.N+1);
        name = [methName '/' subDir methName '_2d',qname, num2str(nelem),'e','.pdf'];
        name = [outDir name];

        if(saveFigure == 1)
            print(fig,'-dpdf',name);
        end
        pause(0.5);
    end % m
end
%% For reactive flows, plot change in 2*q1+q2 (should be constant)
xwidth = 400; 
ywidth = 400;

for imethod=1:length(whichMethods)
    nmethod = whichMethods(imethod);
    methName = methods{nmethod};
    currMeth = meth.(methName);
    
    disp(['Reading: ' methName]);
    
    q1init=squeeze(currMeth.q1(1,:,:));
    q2init=squeeze(currMeth.q2(1,:,:));
    qTinit=2*q1init+q2init;

    q1final=squeeze(currMeth.q1(end,:,:));
    q2final=squeeze(currMeth.q2(end,:,:));
    qT=2*q1final+q2final;
    
    fig = figure();
    set(gcf, 'PaperUnits', 'points');
    set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
    set(fig, 'Position', [0 0 xwidth ywidth])

    x = currMeth.x;
    y = currMeth.y;

    clvls = 0.1:0.1:2.0;
    contourf(x,y,qT,clvls);
    colorbar;
    hLu = text(0.05,0.95,[currMeth.figLabel ') ' currMeth.method '-- qT'],FS,18); 
    axis([-0.005 1.005 -0.005 1.005]);

    
    normFac = sqrt(mean(qTinit(:).^2));
    el2 = sqrt(mean( (qT(:)-qTinit(:)).^2 ))/normFac;
    normFac = max(abs(qTinit(:)));
    einf = max(abs(qT(:)-qTinit(:)))/normFac;
    
    disp(['el2=' num2str(el2)]);
    disp(['einf=' num2str(einf)]);
    disp('');
end
