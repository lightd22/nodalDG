%% Reactive chemistry exact solution at t=0,T,2T,...
function [q1,q2] = reactiveExact(reactionCoeff,q1_ic,qT,t)
    q1 = (q1_ic.*qT)./(q1_ic+(qT-q1_ic).*exp(reactionCoeff.*qT*t));
    q2 = qT-q1;
end