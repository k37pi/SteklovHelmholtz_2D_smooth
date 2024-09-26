% This function takes the curve name and returns
% the parameterization x, its derivatie dx and second derivative
% d2x
% Inputs : curve name, point set = t, normalize length = len1, curve params if any
% Outputs : x, dx, d2x at given point set
function [x,dx,d2x,nx,len] = curve(name,t,len1,varargin)%curve(name,t,len1,a,varargin)
    % w = @(s) pi*(s.^2)./(s.^2+-2*pi*s+2*pi^2); % p = 2;
    % t = w(t);   
    p = 8; %t1 = 0:2*pi/N:(2*pi);
    v1 = @(s) (1/p-0.5)/pi^3*(pi-s).^3+(s-pi)/(p*pi)+0.5;
    dv1 = @(s) -3*(1/p-0.5)/pi^3*(pi-s).^2+1/(p*pi);
    w = @(s) 2*pi*v1(s).^p./(v1(s).^p+v1(2*pi-s).^p);
    dw = @(s) 2*pi*(p*v1(s).^(p-1).*dv1(s).*(v1(s).^p+v1(2*pi-s).^p)-v1(s).^p.*(p*v1(s).^(p-1).*dv1(s)-p*v1(2*pi-s).^(p-1).*dv1(2*pi-s)))./(v1(s).^p+v1(2*pi-s).^p).^2;
    nn = (length(t)-1)/2;
    switch lower(name)
        case 'circle'
            r = varargin{1};
            x = r*[cos(t);sin(t)];
            dx = r*[-sin(t);cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            d2x = -x;
            len = 2*pi*r;
        case 'circle2'
            r = varargin{1};
            x = [r*cos(t)+2;r*sin(t)];
            dx = r*[-sin(t);cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            d2x = -r*[cos(t);sin(t)];
            len = 2*pi*r;   
            % x = [cos(t);sin(t)];
            % dx = [-sin(t);cos(t)]; %nx = [dx(2);-dx(1)];
            % d2x = -x;
            % len = 2*pi;
        % case 'annulus'
        %     r1 = varargin{1};
        %     r2 =   varargin{2}; 
        %     R = max(r1,r2); r = min(r1,r2);
        %     x = r*[cos(t);sin(t)];
        %     dx = r*[-sin(t);cos(t)]; %nx = [dx(2,:);-dx(1,:)];
        %     d2x = -x;
        %     len = 2*pi*r;
            % x = [cos(t);sin(t)];
            % dx = [-sin(t);cos(t)]; %nx = [dx(2);-dx(1)];
            % d2x = -x;
            % len = 2*pi;
        case 'ellipse'
            a = varargin{1}; b = varargin{2};
            x = [a*cos(t);b*sin(t)];
            dx = [-a*sin(t);b*cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            d2x = -x;
            len = integral(@(t)sqrt(a^2*sin(t).^2+b^2*cos(t).^2),0,2*pi);
        case 'kite'
            x = [cos(t)+0.65*cos(2*t)-0.65;1.5*sin(t)];
            dx = [-sin(t)-1.3*sin(2*t);1.5*cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            d2x = [-cos(t)-2.6*cos(2*t);-1.5*sin(t)];
            len = integral(@(t) sqrt((sin(t)+1.3*sin(2*t)).^2+(1.5*cos(t)).^2),0,2*pi);
        % case 'coolcurve'
        %     % ab = varargin{1}; 
        %     % x = [ab(1)*cos(ab(2)*t)+ab(3)*sin(ab(4)*t);ab(5)*cos(ab(6)*t)+ab(7)*sin(ab(8)*t)];
        %     % dx = [-ab(2)*ab(1)*sin(ab(2)*t)+ab(4)*ab(3)*cos(ab(4)*t);-ab(6)*ab(5)*sin(ab(6)*t)+ab(8)*ab(7)*cos(ab(8)*t)]; %nx = [dx(2);-dx(1)];
        %     % d2x = [-ab(2)^2*ab(1)*cos(ab(2)*t)-ab(4)^2*ab(3)*sin(ab(4)*t);-ab(6)^2*ab(5)*cos(ab(6)*t)-ab(8)^2*ab(7)*sin(ab(8)*t)];
        %     % len = integral(@(t) sqrt((-ab(2)*ab(1)*sin(ab(2)*t)+ab(4)*ab(3)*cos(ab(4)*t)).^2+(-ab(6)*ab(5)*sin(ab(6)*t)+ab(8)*ab(7)*cos(ab(8)*t)).^2),0,2*pi);
        %     ab = varargin{1};
        %     % A = ab(1)*cos(ab(2)*t)+ab(3)*sin(ab(4)*t); At = @(t) ab(1)*cos(ab(2)*t)+ab(3)*sin(ab(4)*t);
        %     % dA = -ab(2)*ab(1)*sin(ab(2)*t)+ab(4)*ab(3)*cos(ab(4)*t); dAt = @(t) -ab(2)*ab(1)*sin(ab(2)*t)+ab(4)*ab(3)*cos(ab(4)*t);
        %     % d2A = -ab(2)^2*ab(1)*cos(ab(2)*t)-ab(4)^2*ab(3)*sin(ab(4)*t);
        %     A = ab(1)*cos(ab(2)*t)+ab(3); At = @(t) ab(1)*cos(ab(2)*t)+ab(3);
        %     dA = -ab(2)*ab(1)*sin(ab(2)*t); dAt = @(t) -ab(2)*ab(1)*sin(ab(2)*t);
        %     d2A = -ab(2)^2*ab(1)*cos(ab(2)*t);
        %     x = A.*[cos(t);sin(t)];
        %     dx = [-A.*sin(t)+dA.*cos(t); A.*cos(t)+dA.*sin(t)]; 
        %     d2x = [-A.*cos(t)-2*dA.*sin(t)+d2A.*cos(t); -A.*sin(t)+2*dA.*cos(t)+d2A.*sin(t)]; 
        %     %nx = [dx(2,:);-dx(1,:)];
        %     len = integral(@(t) sqrt((-At(t).*sin(t)+dAt(t).*cos(t)).^2+(At(t).*cos(t)+dAt(t).*sin(t)).^2),0,2*pi); 
        case 'coolcurve'
            ab = varargin{1};
            nn = size(ab,2);
            x = [1+ab(1,1)*cos(t);1+ab(2,1)*sin(t)];
            dx = [-ab(1,1)*sin(t);ab(2,1)*cos(t)];
            d2x = [ab(1,1)*cos(t);ab(2,1)*sin(t)];         
            %len = integral(@(t) sqrt((ab(1)*sin(t)+2*ab(2)*sin(2*t)).^2+(ab(3)*cos(t)+2*ab(4)*cos(2*t)).^2),0,2*pi); 
            for ii = 2:nn
            x = [x(1,:)+ab(1,ii)*cos(ii*t);x(2,:)+ab(2,ii)*sin(ii*t)];
            dx = [dx(1,:)-ii*ab(1,ii)*sin(ii*t);dx(2,ii)+ii*ab(2,ii)*cos(ii*t)];
            d2x = [d2x(1,ii)+ii^2*ab(1,ii)*cos(ii*t);d2x(2,ii)+ii^2*ab(2,ii)*sin(ii*t)];         
            %len = integral(@(t) sqrt((ab(1)*sin(t)+2*ab(2)*sin(2*t)).^2+(ab(3)*cos(t)+2*ab(4)*cos(2*t)).^2),0,2*pi); 
            end
            d2x = -d2x;
            len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2));
        case 'smoothsquare'
            x = [sqrt(cos(t));sqrt(sin(t))];
            dx = (1/2)*[-sin(t);cos(t)]./x; %nx = [dx(2,:);-dx(1,:)];
            d2x = [-2*cos(2*t);2*cos(2*t)];    
        case 'disk_def'
            k = varargin{1};
            x = [cos(t)+k*cos(2*t)-k;1.5*sin(t)];
            dx = [-sin(t)-2*k*sin(2*t);1.5*cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            d2x = [-cos(t)-4*k*cos(2*t);-1.5*sin(t)];
            len = integral(@(t) sqrt((sin(t)+2*k*sin(2*t)).^2+(1.5*cos(t)).^2),0,2*pi);
            v = [cos(2*t)-1;zeros(1,length(t))]; vn = v.*[dx(2,:);-dx(1,:)]./sqrt(dx(1,:).^2+dx(2,:).^2);
            kappa = (dx(1,:).*d2x(2,:)-d2x(1,:).*dx(2,:))./(dx(1,:).^2+dx(2,:).^2).^1.5; 
        case 'acorn'
            A = sqrt(4.25+2*cos(3*t)); At = @(t) sqrt(4.25+2*cos(3*t));
            dA = -3*sin(3*t)./A; dAt = @(t) -3*sin(3*t)./At(t);
            d2A = -9*cos(3*t)./A+3*sin(3*t).*dA./A.^2;
            x = 0.6*A.*[cos(t);sin(t)];
            dx = 0.6*[-A.*sin(t)+dA.*cos(t); A.*cos(t)+dA.*sin(t)]; 
            d2x = 0.6*[-A.*cos(t)-2*dA.*sin(t)+d2A.*cos(t); -A.*sin(t)+2*dA.*cos(t)+d2A.*sin(t)]; 
            %nx = [dx(2,:);-dx(1,:)];
            len = integral(@(t) 0.6*sqrt((-At(t).*sin(t)+dAt(t).*cos(t)).^2+(At(t).*cos(t)+dAt(t).*sin(t)).^2),0,2*pi);
        case 'polar2'
            a = varargin{1};
            A = 1+0.3*cos(a*(t+0.2*sin(t))); 
            dA = -0.3*sin(a*(t+0.2*sin(t))).*(a*(1+0.2*cos(t))); 
            d2A = -0.3*cos(a*(t+0.2*sin(t))).*(a*(1+0.2*cos(t))).^2-0.3*sin(a*(t+0.2*sin(t))).*(-a*(0.2*sin(t)));
            x = A.*[cos(t);sin(t)];
            dx = [-A.*sin(t)+dA.*cos(t); A.*cos(t)+dA.*sin(t)]; 
            d2x = [-A.*cos(t)-2*dA.*sin(t)+d2A.*cos(t); -A.*sin(t)+2*dA.*cos(t)+d2A.*sin(t)]; 
            %nx = [dx(2,:);-dx(1,:)];
            len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2));
        
        case 'asypolar'
            a = varargin{1};
            A = 1+0.3*cos(a*(t+0.2*cos(t))); 
            dA = -0.3*sin(a*(t+0.2*cos(t))).*(a*(1-0.2*sin(t))); 
            d2A = -0.3*cos(a*(t+0.2*cos(t))).*(a*(1-0.2*sin(t))).^2-0.3*sin(a*(t+0.2*cos(t))).*(-a*(0.2*cos(t)));
            x = A.*[cos(t);sin(t)];
            dx = [-A.*sin(t)+dA.*cos(t); A.*cos(t)+dA.*sin(t)]; 
            d2x = [-A.*cos(t)-2*dA.*sin(t)+d2A.*cos(t); -A.*sin(t)+2*dA.*cos(t)+d2A.*sin(t)]; 
            %nx = [dx(2,:);-dx(1,:)];
            len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2));   

        case 'star'
            aa = 5;
            A = 1+0.3*cos(aa*t); At = @(t) 1+0.3*cos(aa*t);
            dA = -aa*0.3*sin(aa*t); dAt = @(t) -aa*0.3*sin(aa*t);
            d2A = -aa^2*0.3*cos(aa*t);
            x = A.*[cos(t);sin(t)];
            dx = [-A.*sin(t)+dA.*cos(t); A.*cos(t)+dA.*sin(t)]; 
            d2x = [-A.*cos(t)-2*dA.*sin(t)+d2A.*cos(t); -A.*sin(t)+2*dA.*cos(t)+d2A.*sin(t)]; 
            %nx = [dx(2,:);-dx(1,:)];
            len = integral(@(t) sqrt((-At(t).*sin(t)+dAt(t).*cos(t)).^2+(At(t).*cos(t)+dAt(t).*sin(t)).^2),0,2*pi);
        case 'flower'
            aa = varargin{1};
            A = 1+0.3*cos(aa*t); At = @(t) 1+0.3*cos(aa*t);
            dA = -aa*0.3*sin(aa*t); dAt = @(t) -aa*0.3*sin(aa*t);
            d2A = -aa^2*0.3*cos(aa*t);
            x = A.*[cos(t);sin(t)];
            dx = [-A.*sin(t)+dA.*cos(t); A.*cos(t)+dA.*sin(t)]; 
            d2x = [-A.*cos(t)-2*dA.*sin(t)+d2A.*cos(t); -A.*sin(t)+2*dA.*cos(t)+d2A.*sin(t)]; 
            %nx = [dx(2,:);-dx(1,:)];
            len = integral(@(t) sqrt((-At(t).*sin(t)+dAt(t).*cos(t)).^2+(At(t).*cos(t)+dAt(t).*sin(t)).^2),0,2*pi);            
        case 'butterfly'
            A = exp(cos(t)).*cos(2*t).^2+exp(sin(t)).*sin(2*t).^2; 
            At = @(t) exp(cos(t)).*cos(2*t).^2+exp(sin(t)).*sin(2*t).^2;
            dA = exp(cos(t)).*(-sin(t).*cos(2*t).^2-2*sin(4*t))+exp(sin(t)).*(cos(t).*sin(2*t).^2+2*sin(4*t)); 
            dAt = @(t) exp(cos(t)).*(-sin(t).*cos(2*t).^2-2*sin(4*t))+exp(sin(t)).*(cos(t).*sin(2*t).^2+2*sin(4*t));
            d2A = exp(cos(t)).*(-sin(t)).*(-sin(t).*cos(2*t).^2-2*sin(4*t))+exp(cos(t)).*(-cos(t).*cos(2*t).^2+2*sin(t).*sin(4*t)-8*cos(4*t))+...
                  exp(sin(t)).*(cos(t)).*(cos(t).*sin(2*t).^2+2*sin(4*t))+exp(sin(t)).*(-sin(t).*sin(2*t).^2+2*cos(t).*sin(4*t)+8*cos(4*t));
            x = A.*[cos(t);sin(t)];
            dx = [-A.*sin(t)+dA.*cos(t); A.*cos(t)+dA.*sin(t)]; 
            d2x = [-A.*cos(t)-2*dA.*sin(t)+d2A.*cos(t); -A.*sin(t)+2*dA.*cos(t)+d2A.*sin(t)]; 
            %nx = [dx(2,:);-dx(1,:)];
            len = integral(@(t) sqrt((-At(t).*sin(t)+dAt(t).*cos(t)).^2+(At(t).*cos(t)+dAt(t).*sin(t)).^2),0,2*pi);

        case 'insect'
            AA = varargin{1};
            aa = AA(1); bb = AA(2); cc = AA(3); dd = AA(4);
            A = exp(cos(aa*t)).*cos(bb*t).^2+exp(sin(cc*t)).*sin(dd*t).^2; 
            At = @(t) exp(cos(aa*t)).*cos(bb*t).^2+exp(sin(cc*t)).*sin(dd*t).^2;

            dA = exp(cos(aa*t)).*(-aa*sin(aa*t).*cos(bb*t).^2-bb*sin(2*bb*t))+exp(sin(cc*t)).*(cc*cos(cc*t).*sin(dd*t).^2+dd*sin(2*dd*t)); 
            dAt = @(t) exp(cos(aa*t)).*(-aa*sin(aa*t).*cos(bb*t).^2-bb*sin(2*bb*t))+exp(sin(cc*t)).*(cc*cos(cc*t).*sin(dd*t).^2+dd*sin(2*dd*t));

            d2A = exp(cos(aa*t)).*(-aa*sin(aa*t)).*(-aa*sin(aa*t).*cos(bb*t).^2-bb*sin(2*bb*t))+...
                  exp(cos(aa*t)).*(-aa^2*cos(aa*t).*cos(bb*t).^2+aa*bb*sin(aa*t).*sin(2*bb*t)-2*bb^2*cos(2*bb*t))+...
                  exp(sin(cc*t)).*(cc*cos(cc*t)).*(cc*cos(cc*t).*sin(dd*t).^2+dd*sin(2*dd*t))+...
                  exp(sin(cc*t)).*(-cc^2*sin(cc*t).*sin(dd*t).^2+cc*dd*cos(cc*t).*sin(2*dd*t)+2*dd^2*cos(2*dd*t));
            x = A.*[cos(t);sin(t)];
            dx = [-A.*sin(t)+dA.*cos(t); A.*cos(t)+dA.*sin(t)]; 
            d2x = [-A.*cos(t)-2*dA.*sin(t)+d2A.*cos(t); -A.*sin(t)+2*dA.*cos(t)+d2A.*sin(t)]; 
            %nx = [dx(2,:);-dx(1,:)];
            len = integral(@(t) sqrt((-At(t).*sin(t)+dAt(t).*cos(t)).^2+(At(t).*cos(t)+dAt(t).*sin(t)).^2),0,2*pi);
    

        case 'cavity'
            A1 = 0.25*(cos(t)+2*cos(2*t)); %A1t = @(t) 0.25*cos(t)+2*cos(2*t);
            A2 = 0.5*(sin(t)+sin(2*t)+0.5*sin(3*t)); %A2t = @(t)0.5*(sin(t)+sin(2*t)+0.5*sin(3*t));
            A3 = 1/48*(4*sin(t)-7*sin(2*t)+6*sin(3*t)-2*sin(4*t));
            %A3t = @(t) 1/48*(4*sin(t)-7*sin(2*t)+6*sin(3*t)-2*sin(4*t));

            dA1 = -0.25*(sin(t)+4*sin(2*t)); dA1t = @(t) -0.25*(sin(t)+4*sin(2*t));
            dA2 = 0.5*(cos(t)+2*cos(2*t)+1.5*cos(3*t)); dA2t = @(t)0.5*(cos(t)+2*cos(2*t)+1.5*cos(3*t));
            dA3 = 1/48*(4*cos(t)-14*cos(2*t)+18*cos(3*t)-8*cos(4*t));
            dA3t = @(t) 1/48*(4*cos(t)-14*cos(2*t)+18*cos(3*t)-8*cos(4*t));

            d2A1 = -0.25*(cos(t)+8*cos(2*t)); %d2A1t = @(t) -0.25*cos(t)-8*cos(2*t);
            d2A2 = -0.5*(sin(t)+4*sin(2*t)+4.5*sin(3*t)); %d2A2t = @(t)-0.5*(sin(t)+4*sin(2*t)+4.5*sin(3*t));
            d2A3 = 1/48*(-4*sin(t)+28*sin(2*t)-54*sin(3*t)+32*sin(4*t));
            %d2A3t = @(t) 1/48*(-4*sin(t)+28*sin(2*t)-54*sin(3*t)+32*sin(4*t));

            x = [A1;A2+A3]; dx = [dA1;dA2+dA3]; d2x = [d2A1;d2A2+d2A3];
            len = integral(@(t) sqrt((dA1t(t)).^2+(dA2t(t)+dA3t(t)).^2),0,2*pi);
        case 'rectangle'
            a = varargin{1}; b = varargin{2};
            x = zeros(2,length(t)); dx = zeros(2,length(t)); d2x = zeros(2,length(t)); 
            % x(1,t <= pi/4) = a; x(2,t <= pi/4) = b*t(t <= pi/4)/(pi/4);
            % dx(2,t <= pi/4) = b/(pi/4);
            %x(2,t <= pi/4) = b*polycov(t(t <= pi/4),p)/(pi/4);
                        
            % [w,dw] = polycov(t(t >= pi/4 & t <= 3*pi/4),p);
            % x(1,t >= pi/4 & t <= 3*pi/4) = -2*a/(pi/2)*(w-pi/4)+a;
            % x(2,t >= pi/4 & t <= 3*pi/4) = b;
            % dx(1,t >= pi/4 & t <= 3*pi/4) = -2*a/(pi/2)*dw;
            % W = w; DW = dw;
            % 
            % [w,dw] = polycov(t(t >= 3*pi/4 & t <= 5*pi/4),p);
            % x(1,t >= 3*pi/4 & t <= 5*pi/4) = -a; 
            % x(2,t >= 3*pi/4 & t <= 5*pi/4) = -2*b/(pi/2)*(w-3*pi/4)+b;
            % dx(2,t >= 3*pi/4 & t <= 5*pi/4) = -2*b/(pi/2)*dw;
            % W = [W,w]; DW = [DW,dw];
            % 
            % [w,dw] = polycov(t(t >= 5*pi/4 & t <= 7*pi/4),p);
            % x(1,t >= 5*pi/4 & t <= 7*pi/4) = 2*a/(pi/2)*(w-5*pi/4)-a;
            % x(2,t >= 5*pi/4 & t <= 7*pi/4) = -b;
            % dx(1,t >= 5*pi/4 & t <= 7*pi/4) = 2*a/(pi/2)*dw;  
            % W = [W,w]; DW = [DW,dw];
            % 
            % % x(1,t >= 7*pi/4) = a; x(2,t > 7*pi/4) = b*(t(t > 7*pi/4)-2*pi)/(pi/4);
            % % 
            % % x(2,t >= 7*pi/4) = b*(polycov(t(t >= 7*pi/4),p)-2*pi)/(pi/4);
            % 
            % [w,dw] = polycov([t(t >= 7*pi/4),t(t <= pi/4)+2*pi],p); 
            % cross_x = 2*b/(pi/2)*(w-7*pi/4)-b;
            % x(1,t >= 7*pi/4) = a;
            % x(1,t <= pi/4) = a;
            % x(2,t <= pi/4) = cross_x((end/2+1):end);
            % x(2,t > 7*pi/4) = cross_x(1:end/2);
            % dx(2,t <= pi/4) = 2*b/(pi/2)*dw((end/2+1):end);
            % dx(2,t > 7*pi/4) = 2*b/(pi/2)*dw(1:end/2);
            % W = [w((end/2+1):end)-2*pi,W,w(1:end/2)]; DW = [dw((end/2+1):end),DW,dw((1:end/2))];
            % [W,finind] = unique(W); DW = DW(finind);

            ttemp = t(t<=pi/2);
            x(1,t<=pi/2) = a; x(2,t<=pi/2) = b*ttemp/(pi/2);
            ttemp = t(t>=pi/2 & t<=pi);
            x(1,t>=pi/2 & t<=pi) = a*(pi-ttemp)/(pi/2); x(2,t>=pi/2 & t<=pi) = b;
            ttemp = t(t>=pi & t<=3*pi/2);
            x(1,t>=pi & t<=3*pi/2) = 0; x(2,t>=pi & t<=3*pi/2) = b*(3*pi/2-ttemp)/(pi/2);
            ttemp = t(t>=3*pi/2);
            x(1,t>=3*pi/2) = a*(ttemp-3*pi/2)/(pi/2); x(2,t>=3*pi/2) = 0;

            %nx = [dx(2,:);-dx(1,:)];
            len = 2*(a+b);
        case 'triangle'
            a = varargin{1}; 
            angs = t >= 0 & t <= pi/2;
            [w,dw] = polycov(t(angs),p);
            x(1,angs) = -a*w/(pi/2)+a;
            x(2,angs) = sqrt(3)*a/(pi/2)*w;
            dx(1,angs) = -a/(pi/2)*dw; dx(2,angs) = sqrt(3)*a/(pi/2)*dw;
            d2x(1,angs) = -a/(pi/2)*d2w; d2x(2,angs) = sqrt(3)*a/(pi/2)*d2w;
        
        case 'egg'    
          eps = varargin{1};  
          x = [0.75*cos(t)+eps*cos(2*t);sin(t)];
          dx = [-0.75*sin(t)-2*eps*sin(2*t);cos(t)]; %nx = [dx(2,:);-dx(1,:)];
          d2x = [-0.75*cos(t)-4*eps*cos(2*t);-sin(t)];
          len = integral(@(t) sqrt((-0.75*sin(t)-2*eps*sin(2*t)).^2+(cos(t)).^2),0,2*pi);

        case 'coolpolar'
        
            ab = varargin{1};
            ab = [ab(1:length(ab)/2);ab(length(ab)/2+1:end)];
            a0 = sum(sum(abs(ab)))*2+.1;
            %ab1 = [0.1 0 0 0 0.1;0 0.1 0 0 -0.1];
            %ab1(1,1) = ab;
            %ab = ab1;
            nn = size(ab,2);
            rho = a0;
            drho = 0;
            d2rho = 0;

            for ii = 1:nn
            rho = rho + ab(1,ii)*cos(ii*t) + ab(2,ii)*sin(ii*t); 
            drho = drho - ii*ab(1,ii)*sin(ii*t) + ii*ab(2,ii)*cos(ii*t); 
            d2rho = d2rho + ii^2*ab(1,ii)*cos(ii*t) + ii^2*ab(2,ii)*sin(ii*t); 
            end
            d2rho = -d2rho;
            x1 = rho.*cos(t); x2 = rho.*sin(t); 
            x = [x1;x2]; clear x1 x2

            x1 = drho.*cos(t)-rho.*sin(t);
            x2 = drho.*sin(t)+rho.*cos(t);
            dx = [x1;x2]; clear x1 x2
            
            x1 = d2rho.*cos(t)-2*drho.*sin(t)-rho.*cos(t);
            x2 = d2rho.*sin(t)+2*drho.*cos(t)-rho.*sin(t);
            d2x = [x1;x2]; clear x1 x2  

            len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2));
        
        case 'fourier'
            asx = varargin{1}; bsx = varargin{2};
            asy = varargin{3}; bsy = varargin{4}; fN = length(bsy);
         
            fsx = @(t) asx(1)+sum(asx(2:end).*cos((1:fN)*t)+bsx(1:end).*sin((1:fN)*t));
            fsy = @(t) asy(1)+sum(asy(2:end).*cos((1:fN)*t)+bsy(1:end).*sin((1:fN)*t));
            fxr = arrayfun(@(x) fsx(x),t); fyr = arrayfun(@(x) fsy(x),t);
            x = [fxr;fyr];
            
            dfsx = @(t) sum( (1:fN).*(-asx(2:end).*sin((1:fN)*t)+bsx(1:end).*cos((1:fN)*t)) );
            dfsy = @(t) sum( (1:fN).*(-asy(2:end).*sin((1:fN)*t)+bsy(1:end).*cos((1:fN)*t)) );
            dfxr = arrayfun(@(x) dfsx(x),t); dfyr = arrayfun(@(x) dfsy(x),t);
            dx = [dfxr;dfyr];
      
            d2fsx = @(t) -sum( (1:fN).^2.*(asx(2:end).*cos((1:fN)*t)+bsx(1:end).*sin((1:fN)*t)) );
            d2fsy = @(t) -sum( (1:fN).^2.*(asy(2:end).*cos((1:fN)*t)+bsy(1:end).*sin((1:fN)*t)));
            d2fxr = arrayfun(@(x) d2fsx(x),t); d2fyr = arrayfun(@(x) d2fsy(x),t);
            d2x = [d2fxr;d2fyr];

            len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2));
        
        case 'opt_disk'

            rho = 2.5 + 0.057475351612645*cos(50*t) + ...
            0.002675998736772*cos(100*t) - 0.002569287572637*cos(150*t);
            Rho = @(t) 2.5 + 0.057475351612645*cos(50*t) + ...
            0.002675998736772*cos(100*t) - 0.002569287572637*cos(150*t);

            drho = -50*0.057475351612645*sin(50*t) + ...
            -100*0.002675998736772*sin(100*t) + 150*0.002569287572637*sin(150*t);
            Drho = @(t) -50*0.057475351612645*sin(50*t) + ...
            -100*0.002675998736772*sin(100*t) + 150*0.002569287572637*sin(150*t);


            d2rho = -50^2*0.057475351612645*cos(50*t) + ...
            -100^2*0.002675998736772*cos(100*t) + 150^2*0.002569287572637*cos(150*t);
            
            
            [x1,x2] = pol2cart(t,rho); x = [x1;x2]; clear x1 x2
            
            x1 = drho.*cos(t)-rho.*sin(t);
            x2 = drho.*sin(t)+rho.*cos(t);
            Dx1 = @(t) Drho(t).*cos(t)-Rho(t).*sin(t);
            Dx2 = @(t) Drho(t).*sin(t)+Rho(t).*cos(t);
            dx = [x1;x2]; clear x1 x2
            
            x1 = d2rho.*cos(t)-2*drho.*sin(t)-rho.*cos(t);
            x2 = d2rho.*sin(t)+2*drho.*cos(t)-rho.*sin(t);
            d2x = [x1;x2]; clear x1 x2  

            len = integral(@(t) sqrt(Dx1(t).^2+Dx2(t).^2),0,2*pi);

        case 'bunimovich'
            bb = 18; deltab = 1/4; sb = 2;
            rb = sb*deltab; ab = pi*(rb-1)+bb; cb = ab-bb;
            t1 = pi*rb; t2 = t1+bb-ab; t3 = t2+t1;
            tone = t(t<=t1); ttwo = t(t1<t&t<=t2); tthree = t(t2<t&t<=t3);
            tfour = t(t3<t);
            if ~isempty(tone)
                x = [bb+rb*sin(tone/rb);sb-rb*cos(tone/rb)]; 
                dx = [cos(tone/rb);sin(tone/rb)];
                d2x = [-sin(tone/rb);cos(tone/rb)]/rb;
                len = t1;
            end

            if ~isempty(ttwo)                
                x = [x,[bb+t1-ttwo;(sb+rb)*ones(1,length(ttwo))]];
                dx = [dx,[-ones(1,length(ttwo));zeros(1,length(ttwo))]];
                d2x = [d2x,zeros(2,length(ttwo))];
                len = len + t2-t1;
            end
            
            if ~isempty(tthree)               
                x = [x,[ab+rb*sin((tthree+cb)/rb);sb-rb*cos((tthree+cb)/rb)]]; 
                dx = [dx,[cos((tthree+cb)/rb);sin((tthree+cb)/rb)]];
                d2x = [d2x,[-sin((tthree+cb)/rb);cos((tthree+cb)/rb)]/rb];
                len = len + t3 - t2;
            end
 
            if ~isempty(tfour)                          
                x = [x,[ab-t3+tfour;(sb-rb)*ones(1,length(tfour))]];
                dx = [dx,[ones(1,length(tfour));zeros(1,length(tfour))]];
                d2x = [d2x,zeros(2,length(tfour))];
                len = len + 2*pi-t3;
            end

        case 'sponge'    
          sa = varargin{1};  
          x = [cos(t);sin(t)+sa*sin(3*t)];
          dx = [-sin(t);cos(t)+3*sa*cos(3*t)]; %nx = [dx(2,:);-dx(1,:)];
          d2x = [-cos(t);-sin(t)-9*sa*sin(3*t)];
          len = integral(@(t) sqrt(sin(t).^2+(cos(t)+3*sa*cos(3*t)).^2),0,2*pi);
        case 'sponge2'    
          sa = varargin{1};  
          x = [cos(t)+2.5;sin(t)+sa*sin(3*t)];
          dx = [-sin(t);cos(t)+3*sa*cos(3*t)]; %nx = [dx(2,:);-dx(1,:)];
          d2x = [-cos(t);-sin(t)-9*sa*sin(3*t)];
          len = integral(@(t) sqrt(sin(t).^2+(cos(t)+3*sa*cos(3*t)).^2),0,2*pi);
        case 'pert_circle'
         %t = w(t);   
         t = [t((nn+1):end)-pi t(2:nn+1)+pi];
         f = abs(t-pi).^3; g = exp(-2*(t-pi).^2); h = cos(t);
         df = 3*abs(t-pi).*(t-pi); dg = -4*g.*(t-pi); dh = -sin(t);
         d2f = 6*sign(t-pi).*(t-pi); d2g = 4*g.*(4*(t-pi).^2-1); d2h = -cos(t);
         x = [(f.*g+1).*h; sin(t)];
         dx = [(df.*g+f.*dg).*h+(f.*g+1).*dh; cos(t)];  
         d2x1 = (d2f.*g+2*df.*dg+f.*d2g).*h+2*(df.*g+f.*dg).*dh+(f.*g+1).*d2h;
         d2x = [d2x1;-sin(t)];
         
         %ttemp = [t(((length(t)-1)/2+1):end)-pi t(2:(length(t)-1)/2+1)+pi];
         % t = t(1:end-1);
         % x = x(:,1:end-1);
         % x = [x(:,(length(t)/2+1):end) x(:,1:(length(t)/2))];
         %xtemp = [x(:,((length(t)-1)/2+1):end) x(:,2:(length(t)-1)/2+1)];
         %dxtemp = [dx(:,((length(t)-1)/2+1):end) dx(:,2:(length(t)-1)/2+1)];
         %d2xtemp = [d2x(:,((length(t)-1)/2+1):end) d2x(:,2:(length(t)-1)/2+1)];

         % f1 = abs(t).^3; g1 = exp(-2*(t).^2); h1 = cos(t+pi);
         % df1 = 3*abs(t).*(t); dg1 = -4*g.*(t); dh1 = -sin(t+pi);
         % d2f1 = 6*sign(t).*(t); d2g1 = 4*g.*(4*(t).^2-1); d2h1 = -cos(t+pi);
         % 
         % f2 = abs(t-2*pi).^3; g2 = exp(-2*(t-2*pi).^2); h2 = cos(t-3*pi);
         % df2 = 3*abs(t-2*pi).*(t-2*pi); dg2 = -4*g.*(2*pi-t); dh2 = -sin(t-3*pi);
         % d2f2 = 6*sign(t-2*pi).*(t-2*pi); d2g2 = 4*g.*(4*(t-2*pi).^2-1); d2h2 = -cos(t-3*pi);
         % 
         % x1 = [(f1.*g1+1).*h1; sin(t+pi)];
         % dx1 = [(df1.*g1+f1.*dg1).*h1+(f1.*g1+1).*dh1; cos(t+pi)];  
         % d2x01 = (d2f1.*g1+2*df1.*dg1+f1.*d2g1).*h1+2*(df1.*g1+f1.*dg1).*dh1+(f1.*g1+1).*d2h1;
         % d2x1 = [d2x01;-sin(t+pi)];
         % 
         % x2 = [(f2.*g2+1).*h2; sin(t-3*pi)];
         % dx2 = [(df2.*g2+f2.*dg2).*h2+(f2.*g2+1).*dh2; cos(t-3*pi)];  
         % d2x02 = (d2f2.*g2+2*df2.*dg2+f2.*d2g2).*h2+2*(df2.*g2+f2.*dg2).*dh2+(f2.*g2+1).*d2h2;
         % d2x2 = [d2x02;-sin(t-3*pi)];
         % 
         % x = [x1(:,1:nn) x2(:,(nn+1):end)];
         % dx = [dx1(:,1:nn) dx2(:,(nn+1):end)];
         % d2x = [d2x1(:,1:nn) d2x2(:,(nn+1):end)];

         len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2));  
        case 'ellipse_opt'
            b = 1;
            %a = varargin{1}; 
            %c = 0.8;
            c = varargin{1};
            a = b/sqrt(1-c^2);
            %b = a*sqrt(1-c^2);

            x = [a*cos(t);b*sin(t)];%[b/sqrt(1-c^2)*cos(t);a*sqrt(1-c^2)*sin(t)];
            dx = [-a*sin(t);b*cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            d2x = -x;
            len = integral(@(t)sqrt(a^2*sin(t).^2+b^2*cos(t).^2),0,2*pi);
        case 'tear'
            %t1 = t; t = w(t);
            x = [2*sin(t/2);-sin(t)];
            dx = [cos(t/2);-cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            d2x = [-1/2*sin(t/2);sin(t)];
            len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2));
            %len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2).*dw(t1(1:end-1)));
        % case 'leaf'
        %     x = [;];
         case 'coolpolar_1d'
        
            a = varargin{1};

            a0 = 3;
            n = 5; ab = zeros(2,n)+0.1;
            ab(2,4) = a;
            nn = n;
            rho = a0;
            drho = 0;
            d2rho = 0;

            for ii = 1:nn
            rho = rho + ab(1,ii)*cos(ii*t) + ab(2,ii)*sin(ii*t); 
            drho = drho - ii*ab(1,ii)*sin(ii*t) + ii*ab(2,ii)*cos(ii*t); 
            d2rho = d2rho + ii^2*ab(1,ii)*cos(ii*t) + ii^2*ab(2,ii)*sin(ii*t); 
            end
            d2rho = -d2rho;
            x1 = rho.*cos(t); x2 = rho.*sin(t); 
            x = [x1;x2]; clear x1 x2

            x1 = drho.*cos(t)-rho.*sin(t);
            x2 = drho.*sin(t)+rho.*cos(t);
            dx = [x1;x2]; clear x1 x2
            
            x1 = d2rho.*cos(t)-2*drho.*sin(t)-rho.*cos(t);
            x2 = d2rho.*sin(t)+2*drho.*cos(t)-rho.*sin(t);
            d2x = [x1;x2]; clear x1 x2  

            len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2));
         case 'coolpolar_2d'
        
            a = varargin{1};            
            a0 = 3;
            n = 5; ab = zeros(2,n)+0.1;
            ab(1,5) = a(1);
            ab(2,5) = a(2);

            nn = n;
            rho = a0;
            drho = 0;
            d2rho = 0;

            for ii = 1:nn
            rho = rho + ab(1,ii)*cos(ii*t) + ab(2,ii)*sin(ii*t); 
            drho = drho - ii*ab(1,ii)*sin(ii*t) + ii*ab(2,ii)*cos(ii*t); 
            d2rho = d2rho + ii^2*ab(1,ii)*cos(ii*t) + ii^2*ab(2,ii)*sin(ii*t); 
            end
            d2rho = -d2rho;
            x1 = rho.*cos(t); x2 = rho.*sin(t); 
            x = [x1;x2]; clear x1 x2

            x1 = drho.*cos(t)-rho.*sin(t);
            x2 = drho.*sin(t)+rho.*cos(t);
            dx = [x1;x2]; clear x1 x2
            
            x1 = d2rho.*cos(t)-2*drho.*sin(t)-rho.*cos(t);
            x2 = d2rho.*sin(t)+2*drho.*cos(t)-rho.*sin(t);
            d2x = [x1;x2]; clear x1 x2  

            len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2));
            
        otherwise
             disp('curve not found, try another')
     end    
        if len1 ~=0
            x = x/len*len1; dx = dx/len*len1; d2x = d2x/len*len1; len = len1;
        end    
        nx = [dx(2,:);-dx(1,:)];
        nx = nx./sqrt(dx(1,:).^2+dx(2,:).^2);
        kappa = (dx(1,:).*d2x(2,:)-d2x(1,:).*dx(2,:))./(dx(1,:).^2+dx(2,:).^2).^1.5;
        % if a~=0 && len1~=1
        %    x = a*x; dx = a*dx; d2x = a*d2x; %nx = a*%nx; len = a*len;
        % end    
        %disp("print this again")
end

% plot(t,x(1,:),'--*','linewidth',2)
% hold on
% plot(t,dx(1,:))
% plot(t,d2x(1,:))
% histogram(t,'Normalization','probability')

% t1 = t;
% t = w(t);
% 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2).*dw(t1(1:end-1)))