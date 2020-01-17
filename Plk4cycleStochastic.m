function [Aout,tout,ARout,PPTout,Plk4out] = Plk4cycleStochastic(xPlk4, t, varargin)

% Input dimensional time, but the function outputs dimensionless time,
% by non-dimensionalising with respect to tperiod. The plot option will plot
% using dimensional time

defaultPlk4stage      = 'on';
defaultPPTstage       = 'M';
defaultPlk4shape      = 'switch';
defaultPPTshape       = 'switch';
defaultCdk1           = 'off';
defaultPlot           = 'off';
defaultmultCycle      = 'off';
defaultnumPlk4        = 10^5;
defaulttperiod        = 900;
defaultSfrac          = 0.7;
defaultN              = 10;
defaultnumTrials      = 10;
defaultnumReceptors   = 30;
defaultblockNum       = 0; % Proportion [0,1] of Asl receptors with restricted phosphorylation/unbinding capability
defaultunbindBlockVal = 1; % Fraction [0,1] of blocked unbind rate to normal unbind rate
defaultkinBlockVal    = 9; % Maximum number of phosphorylations on blocked sites

expectedPlk4stage = {'S','M1','M2','M','on','off'};
expectedPPTstage  = {'S','M1','M2','M','on','off'};
expectedPlk4shape = {'switch','spike'};
expectedPPTshape  = {'switch','spike'};
expectedCdk1      = {'on','off'};
expectedPlot      = {'on','off'};
expectedmultCycle = {'on','off'};

p = inputParser;
addRequired(p,'xPlk4');
addRequired(p,'t');
addParameter(p,'tperiod',defaulttperiod);
addParameter(p,'Sfrac',defaultSfrac);
addParameter(p,'N',defaultN);
addParameter(p,'numTrials',defaultnumTrials);
addParameter(p,'numReceptors',defaultnumReceptors);
addParameter(p,'blockNum',defaultblockNum);
addParameter(p,'unbindBlockVal',defaultunbindBlockVal);
addParameter(p,'kinBlockVal',defaultkinBlockVal);
addParameter(p,'numPlk4',defaultnumPlk4);
addParameter(p,'Plk4stage',defaultPlk4stage,@(x) any(validatestring(x,expectedPlk4stage)));
addParameter(p,'PPTstage',defaultPPTstage,@(x) any(validatestring(x,expectedPPTstage)));
addParameter(p,'Plk4shape',defaultPlk4shape,@(x) any(validatestring(x,expectedPlk4shape)));
addParameter(p,'PPTshape',defaultPPTshape,@(x) any(validatestring(x,expectedPPTshape)));
addParameter(p,'Cdk1',defaultCdk1,@(x) any(validatestring(x,expectedCdk1)));
addParameter(p,'plot',defaultPlot,@(x) any(validatestring(x,expectedPlot)));
addParameter(p,'multCycle',defaultmultCycle,@(x) any(validatestring(x,expectedmultCycle)));
parse(p,xPlk4,t,varargin{:});

Plk4stage      = p.Results.Plk4stage;
PPTstage       = p.Results.PPTstage;
Plk4shape      = p.Results.Plk4shape;
PPTshape       = p.Results.PPTshape;
Cdk1           = p.Results.Cdk1;
plotquery      = p.Results.plot;
tperiod        = p.Results.tperiod;
Sfrac          = p.Results.Sfrac;
N              = p.Results.N;
numPlk4        = p.Results.numPlk4;
multCycle      = p.Results.multCycle;
numTrials      = p.Results.numTrials;
numReceptors   = p.Results.numReceptors;
blockNum       = p.Results.blockNum;
unbindBlockVal = p.Results.unbindBlockVal;
kinBlockVal    = p.Results.kinBlockVal;

t = t/tperiod; % Convert time (seconds) to dimensionless timesteps. One unit equals one period of the cycle (default 900 seconds).
xPlk4 = xPlk4*tperiod;

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
        tPPTstart = 0.5*(1 + Sfrac);
        tPPTend   = 1;
    end
    
    if strcmp(Plk4shape,'switch') == 1
        Plk4 = @(t) (mod(t,1) > tPlk4start).*(mod(t,1) < tPlk4end);
    elseif strcmp(Plk4shape,'spike') == 1
        Plk4 = @(t) 0.5*(mod(t,1) > tPlk4start).*(mod(t,1) < tPlk4end).*(1 - cos(2*pi*(mod(t,1) - tPlk4start)/(tPlk4end-tPlk4start)));
    end
    if strcmp(PPTshape,'switch') == 1
        PPT = @(t) (mod(t,1) > tPPTstart).*(mod(t,1) < tPPTend);
    elseif strcmp(PPTshape,'spike') == 1
        PPT  = @(t) 0.5*(mod(t,1) > tPPTstart).*(mod(t,1) < tPPTend).*(1 - cos(2*pi*(mod(t,1) - tPPTstart)/(tPlk4end-tPPTstart)));
    end
    
else
    Plk4 = @(t) t.^0;
    PPT  = @(t) 0.1*t.^0;
end

kKin         = xPlk4(3)*ones(N-1,1);
kPPT         = xPlk4(4)*ones(N-1,1);
kBind        = xPlk4(1)*0.5.^(0:N-1)';
kUnbind      = xPlk4(2)*0.5.^(N-1:-1:0)';
kBind        = zeros(N,1);
kUnbind      = zeros(N,1);
kUnbind(end) = xPlk4(2);
kBind(1)     = xPlk4(1);
 
kKin = [kKin; 0];
kPPT = [0; kPPT];

v1  = 1:N-1;
v2  = 2:N;
vN1 = N+1:2*N-1;
vN2 = N+2:2*N;

w1 = 1:N;
wN = N+1:2*N;

ARout = ones(length(t),numReceptors,numTrials);

block = blockNum < rand(numReceptors,numTrials);

x0        = zeros(numReceptors,numTrials,2*N);
x0(:,:,1) = 1;

Plk4out = ones(length(t),1);
PPTout  = zeros(length(t),1);

timesplit = (rand(numTrials,1) + 9)/10;

for j = 2:length(t)
    
    dt = t(j) - t(j-1);
    if strcmp(multCycle,'on') == 1
        numCentrioles = 2.^(floor(t(j) - timesplit + 1));
        Plk4          = max(0,1 - sum(numCentrioles.*reshape(sum(ARout(j-1,:,:)>N,2),[numTrials,1])/(numPlk4*numTrials)));
        Plk4out(j)    = Plk4;
        kBindnew      = kBind*Plk4;
    else
        kBindnew     = kBind*Plk4(t(j));
        Plk4out(j) = Plk4(t(j)); 
    end
    kPPTnew     = kPPT*PPT(t(j));
    PPTout(j) = PPT(t(j));
    
    for k = 1:numTrials
        if strcmp(multCycle,'on') == 1
            if mod(t(j) - timesplit(k),1) < mod(t(j-1) - timesplit(k),1)
                splitVec          = repmat(rand(numReceptors,1,1) > 0.5,1,1,N);
                x0(:,k,1:N)       = x0(:,k,1:N) + splitVec.*x0(:,k,N+1:2*N);
                x0(:,k,N+1:2*N)   = (1 -  splitVec).*x0(:,k,N+1:2*N);
            end
        end
        for i = 1:numReceptors
            kUnbindnew = kUnbind*(block(i,k) + (1 - block(i,k))*unbindBlockVal);
            kKinnew    = kKin*block(i,k) + (1-block(i,k))*[kKin(1:kinBlockVal+1); zeros(N - kinBlockVal - 1,1)];
            
            P = sparse(w1,w1,exp(-(kBindnew + kPPTnew)*dt),2*N,2*N)...
                + sparse(w1,wN,kBindnew.*(1 - exp(-(kBindnew + kPPTnew)*dt))./(kBindnew + kPPTnew),2*N,2*N)...
                + sparse(v2,v1,kPPTnew(2:end).*(1 - exp(-(kBindnew(2:end) + kPPTnew(2:end))*dt))./(kBindnew(2:end) + kPPTnew(2:end)), 2*N,2*N)...
                + sparse(wN,wN,exp(-(kUnbindnew + kKinnew + kPPTnew)*dt),2*N,2*N)...
                + sparse(wN,w1,kUnbindnew.*(1 - exp(-(kUnbindnew + kKinnew + kPPTnew)*dt))./(kUnbindnew + kKinnew + kPPTnew),2*N,2*N)...
                + sparse(vN1,vN2,kKinnew(1:end-1).*(1 - exp(-(kUnbindnew(1:end-1) + kKinnew(1:end-1) + kPPTnew(1:end-1))*dt))./(kUnbindnew(1:end-1) + kKinnew(1:end-1) + kPPTnew(1:end-1)),2*N,2*N)...
                + sparse(vN2,vN1,kPPTnew(2:end).*(1 - exp(-(kUnbindnew(2:end) + kKinnew(2:end) + kPPTnew(2:end))*dt))./(kUnbindnew(2:end) + kKinnew(2:end) + kPPTnew(2:end)),2*N,2*N);
            
            P(isnan(P)) = 0;
            
            M               = dtmc(P);
            Anew            = simulate(M,1,'X0',reshape(x0(i,k,:),[1,2*N]));
            ARout(j,i,k)    = Anew(2);
            x0(i,k,:)       = zeros(1,2*N);
            x0(i,k,Anew(2)) = 1;
        end
    end
end

Aout = 100*sum(ARout>N,2)/numReceptors;
Aout = reshape(Aout,length(t),numTrials);
tout = t;

if strcmp(plotquery,'on') == 1
    figure
%     subplot(2,1,1)
    hold on
    for k = 1:numTrials
        plot(tperiod*t,Aout(:,k),'linewidth',2)
    end
    if strcmp(multCycle,'on') == 1
        for i = 1:ceil(t(end))
            hold on
            tn = [i - (1-Sfrac), i, i, i - (1-Sfrac)];
            yn = [0 0 100 100];
            hatchfill(patch(tn*tperiod,yn,'r'));
        end
    end
    ylim([0 100])
    xlabel('Time (s)')
    ylabel('% Plk4-bound Asl')
    set(gca,'linewidth',2,'Fontsize',20)
    
    %     subplot(2,1,2)
    figure
    plot(tperiod*t,mean(Aout,2),'k','linewidth',2)
    if strcmp(multCycle,'on') == 1
        for i = 1:ceil(t(end))
            hold on
            tn = [i - (1-Sfrac), i, i, i - (1-Sfrac)];
            yn = [0 0 100 100];
            hatchfill(patch(tn*tperiod,yn,'r'));
        end
        plot(tperiod*tout,100*Plk4out,'k:','linewidth',2)
        yyaxis right
        ax = gca;
        ax.YAxis(2).Color = 'k';
        ylim([0 100])
        ylabel('% Plk4 in cytoplasm')
    end
    ylim([0 100])
    xlabel('Time (s)')
    ylabel('% Plk4-bound Asl')
    set(gca,'linewidth',2,'Fontsize',20)
    
    %     subplot(3,1,3)
%     plot(tperiod*t(1:50:end),mean(Aout(1:50:end,:),2),'k')
%     ylim([0 100])
%     xlabel('Time (s)')
%     ylabel('% Plk4-bound Asl')
end
end