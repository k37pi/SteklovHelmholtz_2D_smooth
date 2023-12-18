% This function takes the curve name and returns
% the parameterization x, its derivatie dx and second derivative
% d2x
% Inputs : curve name, point set = t, normalize length = len1, curve params if any
% Outputs : x, dx, d2x at given point set
function [x,dx,d2x,nx,len] = curve(name,t,len1,varargin)%curve(name,t,len1,a,varargin)
    p = 2;
    switch lower(name)
        case 'circle'
            r = varargin{1};
            x = r*[cos(t);sin(t)];
            dx = r*[-sin(t);cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            d2x = -x;
            len = 2*pi*r;
            % x = [cos(t);sin(t)];
            % dx = [-sin(t);cos(t)]; %nx = [dx(2);-dx(1)];
            % d2x = -x;
            % len = 2*pi;
        % case 'annulus'
        %     r1 = varargin{1};
        %     r2 = varargin{2}; 
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
            x = [1+ab(1)*cos(t)+ab(2)*cos(2*t);1+ab(3)*sin(t)+ab(4)*sin(2*t)];
            dx = [-ab(1)*sin(t)-2*ab(2)*sin(2*t);ab(3)*cos(t)+2*ab(4)*cos(2*t)];
            d2x = -[ab(1)*cos(t)+4*ab(2)*cos(2*t);ab(3)*sin(t)+4*ab(4)*sin(2*t)];         
            len = integral(@(t) sqrt((ab(1)*sin(t)+2*ab(2)*sin(2*t)).^2+(ab(3)*cos(t)+2*ab(4)*cos(2*t)).^2),0,2*pi); 
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
        case 'star'
            A = 1+0.3*cos(5*t); At = @(t) 1+0.3*cos(5*t);
            dA = -1.5*sin(5*t); dAt = @(t) -1.5*sin(5*t);
            d2A = -7.5*cos(5*t);
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
                        
            [w,dw] = polycov(t(t >= pi/4 & t <= 3*pi/4),p);
            x(1,t >= pi/4 & t <= 3*pi/4) = -2*a/(pi/2)*(w-pi/4)+a;
            x(2,t >= pi/4 & t <= 3*pi/4) = b;
            dx(1,t >= pi/4 & t <= 3*pi/4) = -2*a/(pi/2)*dw;
            W = w; DW = dw;
            
            [w,dw] = polycov(t(t >= 3*pi/4 & t <= 5*pi/4),p);
            x(1,t >= 3*pi/4 & t <= 5*pi/4) = -a; 
            x(2,t >= 3*pi/4 & t <= 5*pi/4) = -2*b/(pi/2)*(w-3*pi/4)+b;
            dx(2,t >= 3*pi/4 & t <= 5*pi/4) = -2*b/(pi/2)*dw;
            W = [W,w]; DW = [DW,dw];

            [w,dw] = polycov(t(t >= 5*pi/4 & t <= 7*pi/4),p);
            x(1,t >= 5*pi/4 & t <= 7*pi/4) = 2*a/(pi/2)*(w-5*pi/4)-a;
            x(2,t >= 5*pi/4 & t <= 7*pi/4) = -b;
            dx(1,t >= 5*pi/4 & t <= 7*pi/4) = 2*a/(pi/2)*dw;  
            W = [W,w]; DW = [DW,dw];
                        
            % x(1,t >= 7*pi/4) = a; x(2,t > 7*pi/4) = b*(t(t > 7*pi/4)-2*pi)/(pi/4);
            % 
            % x(2,t >= 7*pi/4) = b*(polycov(t(t >= 7*pi/4),p)-2*pi)/(pi/4);

            [w,dw] = polycov([t(t >= 7*pi/4),t(t <= pi/4)+2*pi],p); 
            cross_x = 2*b/(pi/2)*(w-7*pi/4)-b;
            x(1,t >= 7*pi/4) = a;
            x(1,t <= pi/4) = a;
            x(2,t <= pi/4) = cross_x((end/2+1):end);
            x(2,t >= 7*pi/4) = cross_x(1:end/2);
            dx(2,t <= pi/4) = 2*b/(pi/2)*dw((end/2+1):end);
            dx(2,t >= 7*pi/4) = 2*b/(pi/2)*dw(1:end/2);
            W = [w((end/2+1):end)-2*pi,W,w(1:end/2)]; DW = [dw((end/2+1):end),DW,dw((1:end/2))];
            [W,finind] = unique(W); DW = DW(finind);

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

        otherwise 
            disp('curve not found, try another')
     end    
        if len1 ~=0
            x = x/len; dx = dx/len; d2x = d2x/len; len = 1;
        end    
        nx = [dx(2,:);-dx(1,:)];
        nx = nx./sqrt(dx(1,:).^2+dx(2,:).^2);
        kappa = (dx(1,:).*d2x(2,:)-d2x(1,:).*dx(2,:))./(dx(1,:).^2+dx(2,:).^2).^1.5;
        % if a~=0 && len1~=1
        %    x = a*x; dx = a*dx; d2x = a*d2x; %nx = a*%nx; len = a*len;
        % end    
        %disp("print this again")
end
