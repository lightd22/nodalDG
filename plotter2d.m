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

which_test = tests(6);
which_res = res(1);

ncfilename = strcat('spltMod2d_' ,which_test{1}, '.nc');
%%
nc = ['_pdModal/', ncfilename];
methname = 'Modal PD';
out = plot_2dadv(methname,which_test,nc,which_res);
