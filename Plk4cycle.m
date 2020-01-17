function [Aout,tout,PPTout,Plk4out] = Plk4cycle(xPlk4, t, varargin)
    
    % Input dimensional time, but the function outputs dimensionless time,
    % by non-dimensionalising with respect to tperiod. The plot option will plot
    % using dimensional time
    
    defaultOffset         = 'on';
    defaultPlot           = 'off';
    defaultNormalize      = 'off';
    defaultPlk4stage      = 'on';
    defaultPPTstage       = 'M';
    defaultPlk4shape      = 'switch';
    defaultPPTshape       = 'switch';
    defaultCdk1           = 'on';
    defaultmultCycle      = 'off';
    defaultnumReceptors   = 30;
    defaultnumPlk4        = 10^5;
    defaulttperiod        = 900;
    defaultSfrac          = 0.7;
    defaultN              = 10;
    defaultblockNum       = 0; % Proportion [0,1] of Asl receptors with restricted phosphorylation/unbinding capability
    defaultunbindBlockVal = 1; % Fraction [0,1] of blocked unbind rate to normal unbind rate
    defaultkinBlockVal    = 9; % Maximum number of phosphorylations on blocked sites
    
    expectedOffset    = {'on','off'};
    expectedPlot      = {'on','off'};
    expectedNormalize = {'on','off'};
    expectedPlk4stage = {'S','M1','M2','M','on','off'};
    expectedPPTstage  = {'S','M1','M2','M','on','off'};
    expectedPlk4shape = {'switch','spike'};
    expectedPPTshape  = {'switch','spike'};
    expectedCdk1      = {'on','off'};
    expectedmultCycle = {'on','off'};
    
    p = inputParser;
    addRequired(p,'xPlk4');
    addRequired(p,'t');
    addParameter(p,'tperiod',defaulttperiod);
    addParameter(p,'Sfrac',defaultSfrac);
    addParameter(p,'N',defaultN);
    addParameter(p,'numReceptors',defaultnumReceptors);
    addParameter(p,'numPlk4',defaultnumPlk4);
    addParameter(p,'blockNum',defaultblockNum);
    addParameter(p,'unbindBlockVal',defaultunbindBlockVal);
    addParameter(p,'kinBlockVal',defaultkinBlockVal);
    addParameter(p,'Plk4stage',defaultPlk4stage,@(x) any(validatestring(x,expectedPlk4stage)));
    addParameter(p,'PPTstage',defaultPPTstage,@(x) any(validatestring(x,expectedPPTstage)));
    addParameter(p,'Plk4shape',defaultPlk4shape,@(x) any(validatestring(x,expectedPlk4shape)));
    addParameter(p,'PPTshape',defaultPPTshape,@(x) any(validatestring(x,expectedPPTshape)));
    addParameter(p,'Cdk1',defaultCdk1,@(x) any(validatestring(x,expectedCdk1)));
    addParameter(p,'offset',defaultOffset,@(x) any(validatestring(x,expectedOffset)));
    addParameter(p,'plot',defaultPlot,@(x) any(validatestring(x,expectedPlot)));
    addParameter(p,'normalize',defaultNormalize,@(x) any(validatestring(x,expectedNormalize)));
    addParameter(p,'multCycle',defaultmultCycle,@(x) any(validatestring(x,expectedmultCycle)));
    parse(p,xPlk4,t,varargin{:});
    
    Plk4stage      = p.Results.Plk4stage;
    PPTstage       = p.Results.PPTstage;
    Plk4shape      = p.Results.Plk4shape;
    PPTshape       = p.Results.PPTshape;
    Cdk1           = p.Results.Cdk1;
    offset         = p.Results.offset;
    plotquery      = p.Results.plot;
    normalize      = p.Results.normalize;
    tperiod        = p.Results.tperiod;
    numReceptors   = p.Results.numReceptors;
    numPlk4        = p.Results.numPlk4;
    multCycle      = p.Results.multCycle;
    Sfrac          = p.Results.Sfrac;
    N              = p.Results.N;
    blockNum       = p.Results.blockNum;
    unbindBlockVal = p.Results.unbindBlockVal;
    kinBlockVal    = p.Results.kinBlockVal;
    
    
    t = t/tperiod; % Convert time (seconds) to dimensionless timesteps. One unit equals one period of the cycle (default 900 seconds).
    xPlk4 = xPlk4*tperiod;
    
%     xPlk4 = xPlk4*(9/(kinBlockVal))^(0.3);
    
    if strcmp(Cdk1,'on') == 1
        if strcmp(Plk4stage,'on') == 1
            tPlk4start = 0;
            tPlk4end   = 1;
        elseif strcmp(Plk4stage,'off') == 1
            tPlk4start = 0;
            tPlk4end   = -1;
        elseif strcmp(Plk4stage,'S') == 1
            tPlk4start = 0;
            tPlk4end   = Sfrac;
        elseif strcmp(Plk4stage,'M') == 1
            tPlk4start = Sfrac;
            tPlk4end   = 1;
        elseif strcmp(Plk4stage,'M1') == 1
            tPlk4start = Sfrac;
            tPlk4end   = 0.5*(1 + Sfrac);
        elseif strcmp(Plk4stage,'M2') == 1
            tPlk4start = 0.5*(1 + Sfrac);
            tPlk4end   = 1;
        end
        
        if strcmp(PPTstage,'on') == 1
            tPPTstart = 0;
            tPPTend   = 1;
        elseif strcmp(PPTstage,'off') == 1
            tPPTstart = 0;
            tPPTend   = -1;
        elseif strcmp(PPTstage,'S') == 1
            tPPTstart = 0;
            tPPTend   = Sfrac;
        elseif strcmp(PPTstage,'M') == 1
            tPPTstart = Sfrac;
            tPPTend   = 1;
        elseif strcmp(PPTstage,'M1') == 1
            tPPTstart = Sfrac;
            tPPTend   = 0.5*(1 + Sfrac);
        elseif strcmp(PPTstage,'M2') == 1
            tPPTstart = 0.5*(1 + Sfrac) - 0.125*(1-Sfrac) ;
            tPPTend   = 1 - 0.125*(1-Sfrac) ;
        end
        
        if strcmp(Plk4shape,'switch') == 1
            Plk4 = @(t) (mod(t,1) >= tPlk4start).*(mod(t,1) < tPlk4end);
        elseif strcmp(Plk4shape,'spike') == 1
            Plk4 = @(t) 0.5*(mod(t,1) >= tPlk4start).*(mod(t,1) < tPlk4end).*(1 - cos(2*pi*(mod(t,1) - tPlk4start)/(tPlk4end-tPlk4start)));
        end
        if strcmp(PPTshape,'switch') == 1
            PPT = @(t) (mod(t,1) >= tPPTstart).*(mod(t,1) < tPPTend);
        elseif strcmp(PPTshape,'spike') == 1
            PPT  = @(t) 0.5*(mod(t,1) >= tPPTstart).*(mod(t,1) < tPPTend).*(1 - cos(2*pi*(mod(t,1) - tPPTstart)/(tPPTend-tPPTstart)));
        end
        
    else
        Plk4 = @(t) t.^0;
        PPT  = @(t) 0.1*t.^0;
    end
    
    kKin        = xPlk4(3)*ones(N-1,1);
    kPPT        = xPlk4(4)*ones(N-1,1);
    kBind       = xPlk4(1)*0.5.^(0:N-1)';
    kUnbind     = xPlk4(2)*0.5.^(N-1:-1:0)';
    
    kKinBlock    = [xPlk4(3)*ones(kinBlockVal,1); zeros(N - kinBlockVal - 1,1)];
    kUnbindBlock = unbindBlockVal*xPlk4(2)*0.5.^(N-1:-1:0)';
    kBind        = zeros(N,1);
    kUnbind      = zeros(N,1);
    kUnbind(end) = xPlk4(2);
    kBind(1)     = xPlk4(1);
    
    kKin      = [kKin; 0];
    kPPT      = [0; kPPT];
    kKinBlock = [kKinBlock; 0];
    
    waveTime = 0.1;
    
    v1 = 1:N-1;
    v2 = 2:N;
    
    A0 = zeros(4*N,1);
    
    A0(1)     = 1 - blockNum;
    A0(2*N+1) = blockNum;
    
    if strcmp(offset,'on') == 1
        tspan = sort([(0:1:11) , (1:1:11)-0.01 ,11 + t(t>0)]);
    else
        tspan = t;
    end
    
    [tout,Avec] = ode45(@DifEq, tspan, A0);
    
    function out = numCentrioles(t)
        
        out = (2.^(floor(t))).*(1 + max(((t - floor(t)) - 1 + waveTime)/waveTime,0));
        
    end

    function out = numCentriolesdt(t)
        
        if strcmp(multCycle,'on') == 1
            out = (2.^(floor(t)))*heaviside(t - floor(t) - 1 + waveTime)/waveTime;
        else
            out = 0;
        end
        
    end

    function dAvdt = DifEq(t,Av)
        
        if strcmp(multCycle,'on') == 1
            Plk4     = max(1 - numCentrioles(t)*numReceptors*sum([Av(N+1:2*N);Av(3*N+1:4*N)])/numPlk4,0);
            kBindnew = kBind*Plk4;
        else
            kBindnew = kBind*Plk4(t(end));
        end
        
        kPPTnew  = kPPT*PPT(t(end));
        MA       = -diag(kPPTnew + kBindnew) + sparse(v1,v2,kPPTnew(2:end),N,N);
        vA       = kUnbind;
        MAblock  = -diag(kPPTnew + kBindnew) + sparse(v1,v2,kPPTnew(2:end),N,N);
        vAblock  = kUnbindBlock;
        MAP      = -diag(kPPTnew + kKin + kUnbind) + sparse(v1,v2,kPPTnew(2:end),N,N) + sparse(v2,v1,kKin(1:end-1),N,N);
        vAP      = kBindnew;
        MAPblock = -diag(kPPTnew + kKinBlock + kUnbindBlock) + sparse(v1,v2,kPPTnew(2:end),N,N) + sparse(v2,v1,kKinBlock(1:end-1),N,N);
        vAPblock = kBindnew;
        dA       = MA*Av(1:N) + vA.*Av(N+1:2*N) + numCentriolesdt(t)*Av(N+1:2*N)/numCentrioles(t);
        dAP      = MAP*Av(N+1:2*N) + vAP.*Av(1:N) - numCentriolesdt(t)*Av(N+1:2*N)/numCentrioles(t);
        dAblock  = MAblock*Av(2*N+1:3*N) + vAblock.*Av(3*N+1:4*N) + numCentriolesdt(t)*Av(3*N+1:4*N)/numCentrioles(t);
        dAPblock = MAPblock*Av(3*N+1:4*N) + vAPblock.*Av(2*N+1:3*N) - numCentriolesdt(t)*Av(3*N+1:4*N)/numCentrioles(t);
        
        dAvdt = [dA ; dAP; dAblock; dAPblock];
        
    end

    if strcmp(offset,'on') == 1
        indx = tout >= 11;
        Avec = Avec(indx,:)';
        Aout = sum([Avec(N+1:2*N,:);Avec(3*N+1:4*N,:)],1);
        tout = tout(indx)';
    else
        Avec = Avec';
        Aout = sum([Avec(N+1:2*N,:);Avec(3*N+1:4*N,:)],1);
        tout = tout' - tout(1);
    end
    
    if strcmp(multCycle,'on') == 1
        Plk4out = 1 - numCentrioles(t).*numReceptors.*Aout/numPlk4;
    else
        Plk4out = Plk4(tout);
    end
    
    PPTout  = PPT(tout);
    
    Aout = 100*Aout;
    
    if strcmp(plotquery,'on') == 1
        figure
        clf
        hold on
        plot(tperiod*tout,Aout,'k','linewidth',2)
        %         plot(tperiod*tout,100*PPTout,'r--','linewidth',2)
        if strcmp(multCycle,'on') == 1
            plot(tperiod*tout,100*Plk4out,'k:','linewidth',2)
        end
        for i = 1:ceil(t(end))
            hold on
            tn = [i - (1-Sfrac), i, i, i - (1-Sfrac)];
            yn = [0 0 100 100];
            hatchfill(patch(tn*tperiod,yn,'r'));
        end
        ylim([0 100])
        xlabel('Time (s)')
        ylabel('% Plk4-bound Asl')
        set(gca,'linewidth',2,'Fontsize',20)
    end
    
    if strcmp(normalize,'on') == 1
        Aout = Aout/max(Aout);
    end
    
end




