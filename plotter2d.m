% Plot Advection Tests using plot_2dadv.m
% By: Devin Light
% ------

clear all;
close all;
clc;
%%
cd('/Users/Devin/Desktop/R/ModalDG/2d_adv/terminatorTest');

tests = {'adv_sine', ... % 1, Uniform adv of sine^4
         'def_cosinebell', ... % 2, LeVeque deformation test cosinebell
         'def_smth_cosbell', ... % 3, Smoother version of LeVeque test
         'fdef_sqwave', ... % 4, LeVeque deformation test square wave
         'hadv_cosinebell', ... % 5, Horizontal advection of steep cosinebell
         'uniform',... % 6, Uniform field in swirling flow
         'def_cyl',... % 7, Deformation flow applied to slotted cylinder
         'rot_cylinder',... % 8, Solid body rotation of a cylinder
         'rot_cylinder_modified',... %9 solid body rotation for comparison to frank's code
         'consistency',... %10uniform field deformation flow
         };
res = {'1','2','3','4'};

which_test = tests(2);
which_res = res(2);

ncfilename = strcat('spltMod2d_' ,which_test{1}, '.nc');
%% Read in data
nc = ['_pdModal/', ncfilename];
methname = 'Modal PD';
out = plot_2dadv(methname,which_test,nc,which_res);

%%
close all;
x = out.x; y = out.y; t = out.t;
%{
ics = squeeze(out.data(1,:,:));
fig = figure();
contourf(x,y,ics);
title(['t=' num2str(t(1))]);

fig = figure();
final = squeeze(out.data(2,:,:));
contourf(x,y,final);
title(['t=' num2str(t(2))]);

fig = figure();
final = squeeze(out.data(end,:,:));
contourf(x,y,final);
title(['t=' num2str(t(end))]);
%}

fig = figure();
for n=1:length(t)
    plt = squeeze(out.data(n,:,:));
    contourf(x,y,plt);
    title(['t=' num2str(t(n))]);
    pause(1.0);
end
pause(0.5);close all;
