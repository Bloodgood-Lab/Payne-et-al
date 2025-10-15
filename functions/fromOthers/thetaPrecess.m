function [Circ,lin,circFit]=thetaPrecess(SpkPhase,SpkPosition,SlopeRange)
% Code originally from Li Yuan in S. Leutgeb Lab

    %if length(SpkPhase)>= 5
    if length(SpkPhase)>= 0 % moved this setting to external loops
        [Circ.Alpha,Circ.Phi0] = CircLinFit(SpkPhase,SpkPosition,SlopeRange); 
        [Circ.Coeff,Circ.pValue] = CircularCoeff(SpkPhase,SpkPosition,Circ.Alpha,Circ.Phi0);
        [r,p] =corrcoef(SpkPosition,SpkPhase);
        lin.r = r(1,2); lin.p = p(1,2);
        [P,S] = polyfit(SpkPosition,SpkPhase',1);  %matlab default function
        lin.Alpha = P(1); circFit.Alpha = P(1); 
        lin.Phi0 = P(2); circFit.Phi0 = P(2); 
        [circFit.rho, circFit.p] = circ_corrcl(SpkPhase, SpkPosition);
    else
        Circ.Alpha = NaN;
        Circ.Phi0 = NaN;
        Circ.Coeff = NaN;
        Circ.pValue = NaN;
        lin.r = NaN;
        lin.Alpha = NaN;
        lin.Phi0 = NaN;
    end
    
    % alpha: slope of fitting. range:[-2*2pi 2*2pi] (Robert Schmidt,2009,Single-Trial
    % Place Precession in the Hippocampus)
    % phi0: phase offset
    function [alpha,phi0] = CircLinFit(phi,x,a)
        %phi=2*pi.*phi/360;
        %a=-2:0.01:2; %feel free to change
        %a = -5:0.001:5;
        C=zeros(length(a),1);
        S=zeros(length(a),1);
        C_slope=zeros(length(a),1);
        S_slope=zeros(length(a),1);
        for j=1:length(a)
            for i = 1:length(phi)
                C(j) = C(j) + cos(phi(i)-2*pi*a(j)*x(i));
                S(j) = S(j) + sin(phi(i)-2*pi*a(j)*x(i));
            end
            R(j) = sqrt((C(j)/length(phi)).^2+(S(j)/length(phi)).^2); %fit of the line
        end
        MaxR = max(R); %get the best fit
        % In case there is over 1 max value
        % Did not solve the problem that there might be several maximum
        idx = find(R==MaxR);
        alpha=a(idx(1));
        phi0 = atan2(S(idx(1)),C(idx(1)));
        % Shift phi0 so that it is between 0 and 2pi
        phi0 = mod(phi0, 2*pi);

    end

    function [coeff,pValue] = CircularCoeff(phi,x,alpha,phi0)
        %phi=2*pi.*phi/360;
        S=0;
        C=0;
        SX=0;
        CX=0;
        YY=0;
        XX1=0;
        XX2=0;
        X = mod(2*pi*abs(alpha).*x,2*pi);
        for j = 1:length(phi)
            S=S+sin(phi(j));
            C=C+cos(phi(j));
            SX = SX+sin(X(j));
            CX=CX+cos(X(j));
        end

        phiEst = atan2(S,C);
        XEst = atan2(SX,CX);

        for i=1:length(phi)
            YY = YY+sin(phi(i)-phiEst)*sin(X(i)-XEst);
            XX1=XX1+(sin(phi(i)-phiEst))^2;
            XX2=XX2+(sin(X(i)-XEst))^2;
        end
        XX = sqrt(XX1*XX2);

        coeff = YY/XX;

        lamda20 = NormDist(2,0,phi,X,phiEst,XEst);
        lamda02 = NormDist(0,2,phi,X,phiEst,XEst);
        lamda22 = NormDist(2,2,phi,X,phiEst,XEst);
        z=coeff*sqrt(length(phi)*lamda20*lamda02/lamda22);
        pValue = 1-erf(abs(z)/sqrt(2));
    end

    function lamda = NormDist(i,j,phi,th,phiEst,thEst)
        lamda = 0;
        for k=1:length(phi)
            lamda = lamda+((sin(phi(k)-phiEst)^i)*(sin(th(k)-thEst)^j));
        end
        lamda=lamda/length(phi);
    end

end
