function [model] = get_model()
    model = model_info;
    model.f = @f;
    model.g = @g;
    model.a = @a;
end

function [out1,out2,out3,out4,out5] = f(s,x,snext,xnext,p)
    n = size(s,1);

    % f
      out1 = zeros(n,5);
      out1(:,1) = p(2).*(-x(:,1) + 1).^(-p(1)) - s(:,1).^p(5).*x(:,2).*x(:,1).^(-p(5) + 1).*(-p(5) + 1).*exp(s(:,3))./x(:,1);
      out1(:,2) = -p(4).*x(:,4).*xnext(:,2) + x(:,2);
      out1(:,3) = -p(4).*((xnext(:,3) + 1)./snext(:,2)).^p(3)./(xnext(:,3) + 1) + x(:,2);
      out1(:,4) = -p(5).*snext(:,1).^p(5).*xnext(:,1).^(-p(5) + 1).*exp(snext(:,3))./snext(:,1) + p(6) + x(:,4) - 1;
      out1(:,5) = x(:,5) - x(:,4).*(xnext(:,3) + 1);

if nargout >= 2

    % df/ds
      out2 = zeros(n,5,4);
      out2(:,1,1) = -p(5).*s(:,1).^p(5).*x(:,2).*x(:,1).^(-p(5) + 1).*(-p(5) + 1).*exp(s(:,3))./(s(:,1).*x(:,1)); % d eq_1 w.r.t. k
      out2(:,1,2) = 0; % d eq_1 w.r.t. mm
      out2(:,1,3) = -s(:,1).^p(5).*x(:,2).*x(:,1).^(-p(5) + 1).*(-p(5) + 1).*exp(s(:,3))./x(:,1); % d eq_1 w.r.t. z
      out2(:,1,4) = 0; % d eq_1 w.r.t. u
      out2(:,2,1) = 0; % d eq_2 w.r.t. k
      out2(:,2,2) = 0; % d eq_2 w.r.t. mm
      out2(:,2,3) = 0; % d eq_2 w.r.t. z
      out2(:,2,4) = 0; % d eq_2 w.r.t. u
      out2(:,3,1) = 0; % d eq_3 w.r.t. k
      out2(:,3,2) = 0; % d eq_3 w.r.t. mm
      out2(:,3,3) = 0; % d eq_3 w.r.t. z
      out2(:,3,4) = 0; % d eq_3 w.r.t. u
      out2(:,4,1) = 0; % d eq_4 w.r.t. k
      out2(:,4,2) = 0; % d eq_4 w.r.t. mm
      out2(:,4,3) = 0; % d eq_4 w.r.t. z
      out2(:,4,4) = 0; % d eq_4 w.r.t. u
      out2(:,5,1) = 0; % d eq_5 w.r.t. k
      out2(:,5,2) = 0; % d eq_5 w.r.t. mm
      out2(:,5,3) = 0; % d eq_5 w.r.t. z
      out2(:,5,4) = 0; % d eq_5 w.r.t. u

    % df/dx
      out3 = zeros(n,5,5);
      out3(:,1,1) = p(2).*p(1).*(-x(:,1) + 1).^(-p(1))./(-x(:,1) + 1) - s(:,1).^p(5).*x(:,2).*x(:,1).^(-p(5) + 1).*(-p(5) + 1).^2.*exp(s(:,3))./x(:,1).^2 + s(:,1).^p(5).*x(:,2).*x(:,1).^(-p(5) + 1).*(-p(5) + 1).*exp(s(:,3))./x(:,1).^2; % d eq_1 w.r.t. n
      out3(:,1,2) = -s(:,1).^p(5).*x(:,1).^(-p(5) + 1).*(-p(5) + 1).*exp(s(:,3))./x(:,1); % d eq_1 w.r.t. lamb
      out3(:,1,3) = 0; % d eq_1 w.r.t. PI
      out3(:,1,4) = 0; % d eq_1 w.r.t. R
      out3(:,1,5) = 0; % d eq_1 w.r.t. I
      out3(:,2,1) = 0; % d eq_2 w.r.t. n
      out3(:,2,2) = 1; % d eq_2 w.r.t. lamb
      out3(:,2,3) = 0; % d eq_2 w.r.t. PI
      out3(:,2,4) = -p(4).*xnext(:,2); % d eq_2 w.r.t. R
      out3(:,2,5) = 0; % d eq_2 w.r.t. I
      out3(:,3,1) = 0; % d eq_3 w.r.t. n
      out3(:,3,2) = 1; % d eq_3 w.r.t. lamb
      out3(:,3,3) = 0; % d eq_3 w.r.t. PI
      out3(:,3,4) = 0; % d eq_3 w.r.t. R
      out3(:,3,5) = 0; % d eq_3 w.r.t. I
      out3(:,4,1) = 0; % d eq_4 w.r.t. n
      out3(:,4,2) = 0; % d eq_4 w.r.t. lamb
      out3(:,4,3) = 0; % d eq_4 w.r.t. PI
      out3(:,4,4) = 1; % d eq_4 w.r.t. R
      out3(:,4,5) = 0; % d eq_4 w.r.t. I
      out3(:,5,1) = 0; % d eq_5 w.r.t. n
      out3(:,5,2) = 0; % d eq_5 w.r.t. lamb
      out3(:,5,3) = 0; % d eq_5 w.r.t. PI
      out3(:,5,4) = -xnext(:,3) - 1; % d eq_5 w.r.t. R
      out3(:,5,5) = 1; % d eq_5 w.r.t. I

    % df/dsnext
      out4 = zeros(n,5,4);
      out4(:,1,1) = 0; % d eq_1 w.r.t. k(1)
      out4(:,1,2) = 0; % d eq_1 w.r.t. mm(1)
      out4(:,1,3) = 0; % d eq_1 w.r.t. z(1)
      out4(:,1,4) = 0; % d eq_1 w.r.t. u(1)
      out4(:,2,1) = 0; % d eq_2 w.r.t. k(1)
      out4(:,2,2) = 0; % d eq_2 w.r.t. mm(1)
      out4(:,2,3) = 0; % d eq_2 w.r.t. z(1)
      out4(:,2,4) = 0; % d eq_2 w.r.t. u(1)
      out4(:,3,1) = 0; % d eq_3 w.r.t. k(1)
      out4(:,3,2) = p(3).*p(4).*((xnext(:,3) + 1)./snext(:,2)).^p(3)./(snext(:,2).*(xnext(:,3) + 1)); % d eq_3 w.r.t. mm(1)
      out4(:,3,3) = 0; % d eq_3 w.r.t. z(1)
      out4(:,3,4) = 0; % d eq_3 w.r.t. u(1)
      out4(:,4,1) = -p(5).^2.*snext(:,1).^p(5).*xnext(:,1).^(-p(5) + 1).*exp(snext(:,3))./snext(:,1).^2 + p(5).*snext(:,1).^p(5).*xnext(:,1).^(-p(5) + 1).*exp(snext(:,3))./snext(:,1).^2; % d eq_4 w.r.t. k(1)
      out4(:,4,2) = 0; % d eq_4 w.r.t. mm(1)
      out4(:,4,3) = -p(5).*snext(:,1).^p(5).*xnext(:,1).^(-p(5) + 1).*exp(snext(:,3))./snext(:,1); % d eq_4 w.r.t. z(1)
      out4(:,4,4) = 0; % d eq_4 w.r.t. u(1)
      out4(:,5,1) = 0; % d eq_5 w.r.t. k(1)
      out4(:,5,2) = 0; % d eq_5 w.r.t. mm(1)
      out4(:,5,3) = 0; % d eq_5 w.r.t. z(1)
      out4(:,5,4) = 0; % d eq_5 w.r.t. u(1)

    % df/dxnext
      out5 = zeros(n,5,5);
      out5(:,1,1) = 0; % d eq_1 w.r.t. n(1)
      out5(:,1,2) = 0; % d eq_1 w.r.t. lamb(1)
      out5(:,1,3) = 0; % d eq_1 w.r.t. PI(1)
      out5(:,1,4) = 0; % d eq_1 w.r.t. R(1)
      out5(:,1,5) = 0; % d eq_1 w.r.t. I(1)
      out5(:,2,1) = 0; % d eq_2 w.r.t. n(1)
      out5(:,2,2) = -p(4).*x(:,4); % d eq_2 w.r.t. lamb(1)
      out5(:,2,3) = 0; % d eq_2 w.r.t. PI(1)
      out5(:,2,4) = 0; % d eq_2 w.r.t. R(1)
      out5(:,2,5) = 0; % d eq_2 w.r.t. I(1)
      out5(:,3,1) = 0; % d eq_3 w.r.t. n(1)
      out5(:,3,2) = 0; % d eq_3 w.r.t. lamb(1)
      out5(:,3,3) = -p(3).*p(4).*((xnext(:,3) + 1)./snext(:,2)).^p(3)./(xnext(:,3) + 1).^2 + p(4).*((xnext(:,3) + 1)./snext(:,2)).^p(3)./(xnext(:,3) + 1).^2; % d eq_3 w.r.t. PI(1)
      out5(:,3,4) = 0; % d eq_3 w.r.t. R(1)
      out5(:,3,5) = 0; % d eq_3 w.r.t. I(1)
      out5(:,4,1) = -p(5).*snext(:,1).^p(5).*xnext(:,1).^(-p(5) + 1).*(-p(5) + 1).*exp(snext(:,3))./(snext(:,1).*xnext(:,1)); % d eq_4 w.r.t. n(1)
      out5(:,4,2) = 0; % d eq_4 w.r.t. lamb(1)
      out5(:,4,3) = 0; % d eq_4 w.r.t. PI(1)
      out5(:,4,4) = 0; % d eq_4 w.r.t. R(1)
      out5(:,4,5) = 0; % d eq_4 w.r.t. I(1)
      out5(:,5,1) = 0; % d eq_5 w.r.t. n(1)
      out5(:,5,2) = 0; % d eq_5 w.r.t. lamb(1)
      out5(:,5,3) = -x(:,4); % d eq_5 w.r.t. PI(1)
      out5(:,5,4) = 0; % d eq_5 w.r.t. R(1)
      out5(:,5,5) = 0; % d eq_5 w.r.t. I(1)

end

        
end

function [out1,out2,out3] = g(s,x,e,p)
    n = size(s,1);

    % g

      out1 = zeros(n,4);
      out1(:,1) = s(:,1).*(-p(6) + 1) + s(:,1).^p(5).*x(:,1).^(-p(5) + 1).*exp(s(:,3)) - s(:,2)./(x(:,3) + 1);
      out1(:,2) = s(:,2).*(p(8).*s(:,4) + p(7) + p(9).*s(:,3) + e(:,2) + 1)./(x(:,3) + 1);
      out1(:,3) = p(10).*s(:,3) + e(:,1);
      out1(:,4) = p(8).*s(:,4) + p(9).*s(:,3) + e(:,2);

if nargout >=2
    % dg/ds
          out2 = zeros(n,4,4);
      out2(:,1,1) = p(5).*s(:,1).^p(5).*x(:,1).^(-p(5) + 1).*exp(s(:,3))./s(:,1) - p(6) + 1; % d eq_1 w.r.t. k
      out2(:,1,2) = -1./(x(:,3) + 1); % d eq_1 w.r.t. mm
      out2(:,1,3) = s(:,1).^p(5).*x(:,1).^(-p(5) + 1).*exp(s(:,3)); % d eq_1 w.r.t. z
      out2(:,1,4) = 0; % d eq_1 w.r.t. u
      out2(:,2,1) = 0; % d eq_2 w.r.t. k
      out2(:,2,2) = (p(8).*s(:,4) + p(7) + p(9).*s(:,3) + e(:,2) + 1)./(x(:,3) + 1); % d eq_2 w.r.t. mm
      out2(:,2,3) = p(9).*s(:,2)./(x(:,3) + 1); % d eq_2 w.r.t. z
      out2(:,2,4) = p(8).*s(:,2)./(x(:,3) + 1); % d eq_2 w.r.t. u
      out2(:,3,1) = 0; % d eq_3 w.r.t. k
      out2(:,3,2) = 0; % d eq_3 w.r.t. mm
      out2(:,3,3) = p(10); % d eq_3 w.r.t. z
      out2(:,3,4) = 0; % d eq_3 w.r.t. u
      out2(:,4,1) = 0; % d eq_4 w.r.t. k
      out2(:,4,2) = 0; % d eq_4 w.r.t. mm
      out2(:,4,3) = p(9); % d eq_4 w.r.t. z
      out2(:,4,4) = p(8); % d eq_4 w.r.t. u
    % dg/dx
          out3 = zeros(n,4,5);
      out3(:,1,1) = s(:,1).^p(5).*x(:,1).^(-p(5) + 1).*(-p(5) + 1).*exp(s(:,3))./x(:,1); % d eq_1 w.r.t. n
      out3(:,1,2) = 0; % d eq_1 w.r.t. lamb
      out3(:,1,3) = s(:,2)./(x(:,3) + 1).^2; % d eq_1 w.r.t. PI
      out3(:,1,4) = 0; % d eq_1 w.r.t. R
      out3(:,1,5) = 0; % d eq_1 w.r.t. I
      out3(:,2,1) = 0; % d eq_2 w.r.t. n
      out3(:,2,2) = 0; % d eq_2 w.r.t. lamb
      out3(:,2,3) = -s(:,2).*(p(8).*s(:,4) + p(7) + p(9).*s(:,3) + e(:,2) + 1)./(x(:,3) + 1).^2; % d eq_2 w.r.t. PI
      out3(:,2,4) = 0; % d eq_2 w.r.t. R
      out3(:,2,5) = 0; % d eq_2 w.r.t. I
      out3(:,3,1) = 0; % d eq_3 w.r.t. n
      out3(:,3,2) = 0; % d eq_3 w.r.t. lamb
      out3(:,3,3) = 0; % d eq_3 w.r.t. PI
      out3(:,3,4) = 0; % d eq_3 w.r.t. R
      out3(:,3,5) = 0; % d eq_3 w.r.t. I
      out3(:,4,1) = 0; % d eq_4 w.r.t. n
      out3(:,4,2) = 0; % d eq_4 w.r.t. lamb
      out3(:,4,3) = 0; % d eq_4 w.r.t. PI
      out3(:,4,4) = 0; % d eq_4 w.r.t. R
      out3(:,4,5) = 0; % d eq_4 w.r.t. I
end
        
end

function [out1,out2,out3] = a(s,x,p)
    n = size(s,1);

    % a

      out1 = zeros(n,6);
      out1(:,1) = s(:,1).^p(5).*x(:,1).^(-p(5) + 1).*exp(s(:,3));
      out1(:,2) = p(7) + s(:,4);
      out1(:,3) = s(:,2)./(x(:,3) + 1);
      out1(:,4) = ((x(:,3) + 1)./s(:,2)).^p(3);
      out1(:,5) = p(2).*(-x(:,1) + 1).^(-p(1));
      out1(:,6) = s(:,2)./(x(:,3) + 1);

if nargout >=2
    % da/ds
          out2 = zeros(n,6,4);
      out2(:,1,1) = p(5).*s(:,1).^p(5).*x(:,1).^(-p(5) + 1).*exp(s(:,3))./s(:,1); % d eq_1 w.r.t. k
      out2(:,1,2) = 0; % d eq_1 w.r.t. mm
      out2(:,1,3) = s(:,1).^p(5).*x(:,1).^(-p(5) + 1).*exp(s(:,3)); % d eq_1 w.r.t. z
      out2(:,1,4) = 0; % d eq_1 w.r.t. u
      out2(:,2,1) = 0; % d eq_2 w.r.t. k
      out2(:,2,2) = 0; % d eq_2 w.r.t. mm
      out2(:,2,3) = 0; % d eq_2 w.r.t. z
      out2(:,2,4) = 1; % d eq_2 w.r.t. u
      out2(:,3,1) = 0; % d eq_3 w.r.t. k
      out2(:,3,2) = 1./(x(:,3) + 1); % d eq_3 w.r.t. mm
      out2(:,3,3) = 0; % d eq_3 w.r.t. z
      out2(:,3,4) = 0; % d eq_3 w.r.t. u
      out2(:,4,1) = 0; % d eq_4 w.r.t. k
      out2(:,4,2) = -p(3).*((x(:,3) + 1)./s(:,2)).^p(3)./s(:,2); % d eq_4 w.r.t. mm
      out2(:,4,3) = 0; % d eq_4 w.r.t. z
      out2(:,4,4) = 0; % d eq_4 w.r.t. u
      out2(:,5,1) = 0; % d eq_5 w.r.t. k
      out2(:,5,2) = 0; % d eq_5 w.r.t. mm
      out2(:,5,3) = 0; % d eq_5 w.r.t. z
      out2(:,5,4) = 0; % d eq_5 w.r.t. u
      out2(:,6,1) = 0; % d eq_6 w.r.t. k
      out2(:,6,2) = 1./(x(:,3) + 1); % d eq_6 w.r.t. mm
      out2(:,6,3) = 0; % d eq_6 w.r.t. z
      out2(:,6,4) = 0; % d eq_6 w.r.t. u
    % da/dx
          out3 = zeros(n,6,5);
      out3(:,1,1) = s(:,1).^p(5).*x(:,1).^(-p(5) + 1).*(-p(5) + 1).*exp(s(:,3))./x(:,1); % d eq_1 w.r.t. n
      out3(:,1,2) = 0; % d eq_1 w.r.t. lamb
      out3(:,1,3) = 0; % d eq_1 w.r.t. PI
      out3(:,1,4) = 0; % d eq_1 w.r.t. R
      out3(:,1,5) = 0; % d eq_1 w.r.t. I
      out3(:,2,1) = 0; % d eq_2 w.r.t. n
      out3(:,2,2) = 0; % d eq_2 w.r.t. lamb
      out3(:,2,3) = 0; % d eq_2 w.r.t. PI
      out3(:,2,4) = 0; % d eq_2 w.r.t. R
      out3(:,2,5) = 0; % d eq_2 w.r.t. I
      out3(:,3,1) = 0; % d eq_3 w.r.t. n
      out3(:,3,2) = 0; % d eq_3 w.r.t. lamb
      out3(:,3,3) = -s(:,2)./(x(:,3) + 1).^2; % d eq_3 w.r.t. PI
      out3(:,3,4) = 0; % d eq_3 w.r.t. R
      out3(:,3,5) = 0; % d eq_3 w.r.t. I
      out3(:,4,1) = 0; % d eq_4 w.r.t. n
      out3(:,4,2) = 0; % d eq_4 w.r.t. lamb
      out3(:,4,3) = p(3).*((x(:,3) + 1)./s(:,2)).^p(3)./(x(:,3) + 1); % d eq_4 w.r.t. PI
      out3(:,4,4) = 0; % d eq_4 w.r.t. R
      out3(:,4,5) = 0; % d eq_4 w.r.t. I
      out3(:,5,1) = p(2).*p(1).*(-x(:,1) + 1).^(-p(1))./(-x(:,1) + 1); % d eq_5 w.r.t. n
      out3(:,5,2) = 0; % d eq_5 w.r.t. lamb
      out3(:,5,3) = 0; % d eq_5 w.r.t. PI
      out3(:,5,4) = 0; % d eq_5 w.r.t. R
      out3(:,5,5) = 0; % d eq_5 w.r.t. I
      out3(:,6,1) = 0; % d eq_6 w.r.t. n
      out3(:,6,2) = 0; % d eq_6 w.r.t. lamb
      out3(:,6,3) = -s(:,2)./(x(:,3) + 1).^2; % d eq_6 w.r.t. PI
      out3(:,6,4) = 0; % d eq_6 w.r.t. R
      out3(:,6,5) = 0; % d eq_6 w.r.t. I
end
                    
end

function [out1] = model_info() % informations about the model

    mod = struct;
    mod.states = { 'k','mm','z','u' };
    mod.controls = { 'n','lamb','PI','R','I' };
    mod.auxiliaries = { 'y','theta','m','uc','un','c' };
    mod.parameters = { 'eta','Psi','Phi','beta','alpha','delta','theta_ss','gamma','varphi','rho','sigma_z','sigma_u' };
    mod.s_ss = [14.9561777181 ; 0.979354410572 ; 0.0 ; 0.0];
    mod.x_ss = [0.31 ; 1.04402658206 ; 0.0125 ; 1.0111223458 ; 1.02376137513];
    mod.params = [ 1.         1.8611691  2.         0.989      0.36       0.019      0.0125  0.5        0.         0.95       0.007      0.0089   ];
    mod.X = cell(2,1);
    mod.X{1} = [ 0.31        1.04402658  0.0125      1.01112235  1.02376138];
    mod.X{2} = [[ -9.31433760e-03  -1.11482973e-17   1.48083605e-01  -1.08921490e-02 ;  -5.05162819e-02  -5.47280051e-17  -6.40424720e-01  -2.96865517e-02 ;  -2.49543488e-02   1.03384433e+00  -3.01278025e-01   3.18737453e-01 ;  -1.79951012e-03  -7.99091242e-19   3.48497997e-02  -8.54976422e-04 ;  -8.93934046e-04  -3.19149246e-18   1.65496314e-02   3.36581242e-01]];
    out1 = mod;

end
