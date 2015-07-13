%% Reactive chemistry exact solution at t=0,T,2T,...
function qOut = reactiveExact(r,q_ic,t)
    qOut = 0.*q_ic;
    stat = size(q_ic,3);
    disp(['Reading ' num2str(stat) ' equations...']);
    if(stat==2)
        q1_ic = q_ic(:,:,1);
        q2_ic = q_ic(:,:,2);
        
        a = r.*(q1_ic+q2_ic);
%        qOut(:,:,1) = q1_ic.*(q1_ic+q2_ic)./(q1_ic+q2_ic.*exp(-a.*t));
%        qOut(:,:,2) = -a.*q2_ic./(r.*(q2_ic+q1_ic.*exp(a.*t)));

%        qOut(:,:,1) = q1_ic.*(q1_ic+q2_ic)./(q1_ic+q2_ic.*exp(a.*t));
%        qOut(:,:,2) = q1_ic+q2_ic-qOut(:,:,1);

      % ====
      % dq1/dt = r*q2^2
      % dq2/dt = -dq1/dt
      % ====
        c = 1./q2_ic;
        qOut(:,:,2) = 1./(r.*t+c);
        qOut(:,:,1) = (q1_ic+q2_ic)-qOut(:,:,2);
        
        %q1 = (q1_ic.*qT)./(q1_ic+(qT-q1_ic).*exp(r.*qT*t));
        %q2 = qT-q1;
    elseif(stat==3)
        q1_ic = q_ic(:,:,1);
        q2_ic = q_ic(:,:,2);
        qT = q1_ic+q2_ic+q_ic(:,:,3);
        
        q1Out = 0.*q1_ic;
        q2Out = 0.*q2_ic;
        q3Out = 0.*q2_ic;
        
        qBar = q1_ic-q2_ic;
        loc = qBar == 0;
        q1Out(loc) = q1_ic(loc)./(q1_ic(loc).*r.*t+1);
        q2Out(loc) = q2_ic(loc)./(q2_ic(loc).*r.*t+1);
        
        a = q2_ic(~loc)./q1_ic(~loc);
        beta = r.*qBar(~loc);
        q1Out(~loc) = (qBar(~loc))./(1-a.*exp(-beta.*t));
        q2Out(~loc) = -a.*qBar(~loc)./(a-exp(beta.*t));
        
        q3Out = qT-(q1Out+q2Out);
        
        qOut(:,:,1) = q1Out;
        qOut(:,:,2) = q2Out;
        qOut(:,:,3) = q3Out;
        
    end
end