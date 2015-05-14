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
whichRes = res(3);
whichTest = tests{ntest};

subDir = '';
whichMethods = [1 2 3 4];
ncfilename = strcat('spltMod2d_' ,whichTest, '.nc');

for imethod=1:length(whichMethods)
    nmethod = whichMethods(imethod);
    methName = methods{nmethod};
    if(nmethod == 1)
        methname = 'Modal Unlimited';
        nc = ['_modal/' subDir ncfilename];
        out = plot_2dadv(methname,whichTest,nc,whichRes,meqn);
        out.figLabel = 'a';
        out.pltStyle = 'k-';
    elseif(nmethod == 2)
        methname = 'Modal PD (trunc)';
        nc = ['_pdModal/trunc/' subDir ncfilename];
        out = plot_2dadv(methname,whichTest,nc,whichRes,meqn);
        out.figLabel = 'b';
        out.pltStyle = 'r--';
    elseif(nmethod == 3)
        methname = 'Modal PD (rescale)';
        nc = ['_pdModal/rescale/' subDir ncfilename];
        out = plot_2dadv(methname,whichTest,nc,whichRes,meqn);
        out.figLabel = 'c';
        out.pltStyle = 'b-.';
    elseif(nmethod == 4)
        methname = 'Modal PD (strictest)';
        nc = ['_pdModal/eqscale/' subDir ncfilename];
        out = plot_2dadv(methname,whichTest,nc,whichRes,meqn);
        out.figLabel = 'd';
        out.pltStyle = 'g:';
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
FS = 'FontSize'; LW = 'LineWidth';
figl2 = figure(1);
figinf = figure(2);
for imethod=1:length(whichMethods)
    nmethod = whichMethods(imethod);
    methName = methods{nmethod};
    currMeth = meth.(methName);
    
    disp(['Reading: ' methName]);
    
    q1init=squeeze(currMeth.q1(1,:,:));
    q2init=squeeze(currMeth.q2(1,:,:));
    qTinit=2*q1init+q2init;
    
    fprintf('qT min=%6.4e \n',min(qTinit(:)));
    disp('');
    fprintf('qT max=%6.4e \n' ,max(qTinit(:)));
    disp('');

    nt = length(currMeth.t);
    el2errs = [];
    einferrs = [];
    
    normFac1 = sqrt(mean(qTinit(:).^2));
    normFac2 = max(abs(qTinit(:)));
    for n=1:nt
        q1=squeeze(currMeth.q1(n,:,:));
        q2=squeeze(currMeth.q2(n,:,:));
        qT=2*q1+q2;

        el2 = sqrt(mean( (qT(:)-qTinit(:)).^2 ))/normFac1;
        einf = max(abs(qT(:)-qTinit(:)))/normFac2;
        
        el2errs = [el2errs el2];
        einferrs = [einferrs einf];
    end

    fig = figure(figl2);
    set(gcf, 'PaperUnits', 'points');
    set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
    set(fig, 'Position', [0 0 xwidth ywidth])
    hold on,plot(currMeth.t,el2errs,currMeth.pltStyle,LW,2);
    disp('');
    
    fig = figure(figinf);
    set(gcf, 'PaperUnits', 'points');
    set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
    set(fig, 'Position', [xwidth 0 xwidth ywidth])
    hold on,plot(currMeth.t,einferrs,currMeth.pltStyle,LW,2);
end
fig = figure(figl2);
set(gca,'YScale','log');
xlabel('time',FS,14); ylabel('L_2 Error',FS,14);
hLu = text(0.25,3.5,['a) L_2 Normalized Errors'],FS,18); 
set(gca,FS,14);
axis([0 5 10^(-6) 10]); box on;

name = '_figs/term_N4E48_l2.pdf';
print(fig,name,'-dpdf');

fig = figure(figinf);
set(gca,'YScale','log');
xlabel('time',FS,14); ylabel('L_{\infty} Error',FS,14);
hLu = text(0.25,3.5,['b) L_{\infty} Normalized Errors'],FS,18); 
set(gca,FS,14);
axis([0 5 10^(-6) 10]); box on;

name = '_figs/term_N4E48_einf.pdf';
print(fig,name,'-dpdf');

%% Make IC plots
FS = 'FontSize';
qT = 4*10^(-6);

x = meth.modal.x; y = meth.modal.y;
fig = figure();
set(gcf, 'PaperUnits', 'points');
set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
set(fig, 'Position', [0 0 xwidth ywidth])

q1 = squeeze(meth.modal.q1(1,:,:));
contourf(x,y,q1);
colorbar;
caxis([0 qT]);
hLu = text(0.05,1.03,['a) Cl_2(x,y,0)'],FS,18); 
xlabel('x',FS,14); ylabel('y',FS,14);
set(gca,FS,14);

name = '_figs/cl2_ics.pdf';
print(fig,name,'-dpdf');

fig = figure();
set(gcf, 'PaperUnits', 'points');
set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
set(fig, 'Position', [0 0 xwidth ywidth])

q2 = squeeze(meth.modal.q2(1,:,:));
contourf(x,y,q2);
hLu = text(0.05,1.05,['b) Cl(x,y,0)'],FS,18); 
colorbar;
caxis([0 qT]);
set(gca,FS,14);
set(gca,'YTickLabel',''); xlabel('x',FS,14); 

name = '_figs/cl_ics.pdf';
print(fig,name,'-dpdf');

%% Coeff k1
fig = figure();
set(gcf, 'PaperUnits', 'points');
set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
set(fig, 'Position', [0 0 xwidth ywidth])

k1 = zeros(length(x),length(y));
xloc = x<0.5;
for j=1:length(y)
    k1(j,xloc) = cos(2.0*pi*(x(xloc)-0.25));
end

contourf(x,y,k1);
colorbar;
xlabel('x',FS,14);ylabel('y',FS,14);set(gca,FS,14);
title('k1(x,y)',FS,14);
name = '_figs/k1.pdf';
print(fig,name,'-dpdf');

%% 
qT = 4*10^(-6);
for imethod = 1:length(whichMethods)
    nmethod = whichMethods(imethod);
    methName = methods{nmethod};
    currMeth = meth.(methName);
    
    fig = figure();
    set(gcf, 'PaperUnits', 'points');
    set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
    set(fig, 'Position', [0 0 xwidth ywidth])

    
    q1 = squeeze(currMeth.q1(end,:,:));
    q2 = squeeze(currMeth.q2(end,:,:));
    qTf = 2.*q1+q2;
    max(qTf(:))
    min(qTf(:))
    contourf(x,y,qTf);
    colorbar; caxis([0 5*10^(-6)]);
    hLu = text(0.05,1.03,[currMeth.figLabel ') Cl_T(x,y,5)'],FS,18); 
    xlabel('x',FS,14); ylabel('y',FS,14);
    set(gca,FS,14);
    name = ['_figs/' methName '_N4E48_final.pdf'];
    print(fig,name,'-dpdf');

end
