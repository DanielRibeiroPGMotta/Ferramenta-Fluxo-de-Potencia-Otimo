%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%     Programa para a solu��o do problema de Fluxo de Pot�ncia via       %
%                via m�todo de Newton-Raphson                            %
%          (formula��o convencional em coordenadas cartesianas)          %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

disp(' ')
disp('     ***************************************************')
disp('     *   Programa de Fluxo de Pot�ncia via  m�todo     *')
disp('     *   de Newton-Raphson - coordenadas retangulares  *')
disp('     *           (formula��o convencional)             *')
disp('     ***************************************************')
disp('     LABSPOT/UFSC                          Fev de 2008')
disp('     R. S. Salgado')
disp(' ')

caso = input(' nome do arquivo de dados (sem extens�o):   ', 's');

if ~exist(caso)
  disp(' ')
  disp('   ********************************')
  disp('   *     ARQUIVO INEXISTENTE!     *')
  disp('   ********************************')
  disp(' ')
  return;
end

   format short e;
   format compact;

disp(' ')
disp('   Itera��o   Max_Desbalan�o       miv          ||dx||')
disp(' ========================================================')

%-----  especifica��o das constantes  -----

   Sbase = 100;             % pot�ncia aparente base
   radeg = 180/pi;
   tol_res = 1.0e-04;
   tol_mi = 1.0e-04;
   maxit = 50;
   j = sqrt(-1);
   fc = 1.0;                % fator de carga

   eval(caso);

%-----  leitura dos dados  -----

   %-----  LEITURA DOS DADOS DO SISTEMA  -----


radeg = 180. / pi;
nlt = size(LT, 1);
nb = size(barras, 1);

no = LT(:, 1);
nd = LT(:, 2);
circ = LT(:, 3);
r = LT(:, 4) / 100.0;
x = LT(:, 5) / 100.0;
b = LT(:, 6) / 200.0;
a = LT(:, 7);
amin = LT(:, 8);
amax = LT(:, 9);

% faca a renumeracao das barras e dos circuitos

AuxBarr = barras(:, 1);

for ij = 1:nb
   barras(ij, 1) = ij;
end 

for ik = 1 : nlt
    for jk = 1:nb
        if (no(ik) == AuxBarr(jk))
            nsb(ik) = jk;
        end
       if (nd(ik) == AuxBarr(jk))
           neb(ik) = jk; 
       end
   end 
end
%
LT(:, 1) = nsb';
LT(:, 2) = neb';

lt_tr = find(amin < amax);
ntraf = length(lt_tr);

barra = barras(:, 1);
tipo = barras(:, 2);
V = barras(:, 3);
ang = barras(:, 4)/radeg;
Pg = barras(:, 5) / Sbase;
Qg = barras(:, 6) / Sbase;
Qgmin = barras(:, 7) / Sbase;
Qgmax = barras(:, 8) / Sbase;
Pd = fc * barras(:, 9) / Sbase;
Qd = fc * barras(:, 10) / Sbase;
shunt = barras(:, 11) / Sbase;
iAr  = barras(:, 12);

%-----  LEITURA DOS LIMITES de tens�o e pot�ncia ativa gerada e custo de gera��o -----

Vmin = Lim(:, 2);
Vmax = Lim(:, 3);
Pgmin = Lim(:, 4) / Sbase;
Pgmax = Lim(:, 5) / Sbase;
bg = Lim(:, 6);
cg = Lim(:, 7);

for ibb=1:nb
     Lim(ij,1) = ibb;
end 
   
%-----  especifica��o dos conjuntos de barras  -----

   bpv = find(tipo == 1 | tipo == 4);
   npv = length(bpv);
   bpq = find(tipo == 0 | tipo == 3);
   npq = length(bpq);
   bf  = find(tipo == 2);
   bsf = find(tipo ~= 2);
   nsf = nb - 1;
   bge = find(tipo == 2 | tipo == 1 | tipo == 4);
   nge = length(bge);

%-----  forma��o da lista de adjac�ncia e da matriz Y_barra -----

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%     Subrotina para a determina��o da matriz admit�ncia de barra        %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----  forma��o da lista de adjac�ncia  -----

l = 1;
for i = 1 : nb
  no_temp = 0;
  nd_temp = 0;
  adj(i) = l;
  for k = 1 : nlt
    if (nsb(k) ~= no_temp) | (neb(k) ~= nd_temp)
      if nsb(k) == i
        b_adj(l) = neb(k);
        l = l + 1;
      end
      if neb(k) == i
        b_adj(l) = nsb(k);
        l = l + 1;
      end
    end
    no_temp = nsb(k);
    nd_temp = neb(k);
  end
end
adj(nb + 1) = l;


%-----  determina��o da matriz Y_barra  -----

Y_barra = sparse(nb, nb);

for i = 1 : nb
  Y_barra(i, i) = j*shunt(i);
end

y = 1.0 ./ (r + j*x);

for l = 1 : nlt
  i = nsb(l);
  k = neb(l);
  if a(l) == 0.0
    Y_barra(i, i) = Y_barra(i, i) + y(l) + j*b(l);
    Y_barra(k, k) = Y_barra(k, k) + y(l) + j*b(l);
    Y_barra(i, k) = Y_barra(i, k) - y(l);
    Y_barra(k, i) = Y_barra(i, k);
  else
    tap = 1.0 / a(l);
    Y_barra(i, i) = Y_barra(i, i) + tap^2 * y(l);
    Y_barra(k, k) = Y_barra(k, k) + y(l);
    Y_barra(i, k) = Y_barra(i, k) - tap * y(l);
    Y_barra(k, i) = Y_barra(i, k);
  end
end

G = real(Y_barra);
B = imag(Y_barra);


%-----  inicializa��o das vari�veis  -----
   
   istart = 1;
   Vref = V;
   e = V;
   if (istart == 0)
      e(bpq) = ones(npq, 1); 
      f = sparse(nb, 1);
   else
      e = V.*cos(ang);
      f = V.*sin(ang);
   end
   mi = 1.0;
   delta = sparse(2*nb-2,1);

%-----  processo iterativo  -----

for iter = 1 : maxit

   %-----  c�lculo das inje��es de pot�ncia -----

   E = e + j * f;
   I = Y_barra * E;
   Ia = real(I);
   Ib = imag(I);
   S = E .* conj(I);
   P = real(S);
   Q = imag(S);

   %-----  verifica��o dos limites de pot�ncia reativa gerada -----

   if (npv > 0)
      Qg(bpv) = Q(bpv) + Qd(bpv);
      if (iter > 2)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%     Subrotina para verifica��o dos limites de pot�ncia reativa gerada  %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if (iter >= 2) & (npv > 0)
      flag1 = 0;   
      for i = 1 : npv
         if (flag1 == 0)   
            if Qg(bpv(i)) > 0.999*Qgmax(bpv(i))
               Qg(bpv(i)) = Qgmax(bpv(i));
               flag1 = bpv(i);
               tipo(flag1) = 0;
               disp(' **********************************************************************')
               fprintf(' * Limite m�ximo de pot�ncia reativa gerada atingido na barra * %d\n',flag1)
               disp(' **********************************************************************')
            end
            if Qg(bpv(i)) < 1.001*Qgmin(bpv(i))
               Qg(bpv(i)) = Qgmin(bpv(i));
               flag1 = bpv(i);
               tipo(flag1) = 0;
               disp(' **********************************************************************')
               fprintf(' * Limite m�nimo de pot�ncia reativa gerada atingido na barra * %d\n',flag1)
               disp(' **********************************************************************')
            end
         end
      end
   
      if flag1 > 0
         bpv = find(tipo == 1);
         npv = length(bpv);
         bpq = find(tipo == 0 | tipo == 3);
         npq = length(bpq);
      end
   end


      end
   end
   
%-----  c�lculo dos desbalan�os de pot�ncia  -----

  dP = P(bsf) - Pg(bsf) + Pd(bsf);
  dQ = Q(bpq) - Qg(bpq) + Qd(bpq);
  if (npv > 0)
     dV = (e(bpv).^2 + f(bpv).^2) - V(bpv).^2;  
     dsa = [ dP; dQ; dV ];
  else
     dsa = [ dP; dQ ];
  end
  
%-----  teste de converg�ncia -----  
  
   tst = max(abs(dsa));
   mdx = norm(mi*delta);
   fprintf(' %8.4d     %8.4d     %8.4d     %8.4d\n',iter,tst,mi,mdx)
   if (tst < tol_res | mi < tol_mi)
      break
   end

%-----  determina��o da matriz jacobiana -----

   J1 =  diag(e)*G + diag(f)*B + diag(Ia);
   J2 =  diag(f)*G - diag(e)*B + diag(Ib);
   J3 =  diag(f)*G - diag(e)*B - diag(Ib);
   J4 = -diag(e)*G - diag(f)*B + diag(Ia);
   J5 = diag(2*e);
   J6 = diag(2*f);

   if (npv > 0) 
      J = [  J1(bsf, bsf)   J2(bsf, bsf)
             J3(bpq, bsf)   J4(bpq, bsf)
             J5(bpv, bsf)   J6(bpv, bsf)  ];
   else
      J = [  J1(bsf, bsf)   J2(bsf, bsf)
             J3(bpq, bsf)   J4(bpq, bsf)  ];
   end
   

%-----  solu��o do sistema linear  -----

   delta = - J \ dsa;

%-----  c�lculo do fator de passo para evitar a diverg�ncia do processo iterativo -----

    aa = dsa;
    bb = - dsa;
  
    xx = zeros(nb, 1);
    xx(bsf) = delta(1 : nb-1) + j * delta(nb : 2*nb-2);
    yy = xx .* conj(Y_barra * xx);
  
    cc = [ real(yy(bsf))
    	   imag(yy(bpq))
    	   real(xx(bpv)).^2 + imag(xx(bpv)).^2 ];
    g0 = sum(aa .* bb);
    g1 = sum(bb .^ 2 + 2 * aa .* cc);
    g2 = 3 * sum(bb .* cc);
    g3 = 2 * sum(cc .^ 2);
  
    equacao = [ g3 g2 g1 g0 ];
    raizes = roots(equacao);
    valor1 = find(imag(raizes) == 0.0 & real(raizes) > 0.0);
    valor2 = find(imag(raizes) == 0.0 & real(raizes) < 0.0);
    
    if length(valor1) == 0
      mi = max(raizes(valor2));
    else
      mi = min(raizes(valor1));
    end
    
%-----  atualiza��o das vari�veis  -----

   if (mi > 1.0)
      mi = 1.0;
   end
  
  e(bsf) = e(bsf) + mi * delta(1 : nb-1);
  f(bsf) = f(bsf) + mi * delta(nb : 2*nb-2);
    
end


format loose;
disp(' ')
if iter == maxit
  disp(' ********************************************')
  disp(' *        Convergencia N�o Alcan�ada        *')
  disp(' ********************************************')
  disp(' ')
  disp(' ')
  disp(' ')
  return
else
disp('     *************************************************')
disp('     *        Solu��o do Fluxo De Pot�ncia           *')
disp('     *  Newton-Raphson em Coordenadas Retangulares   *')
disp('     *************************************************')
disp(' ')
end

%-----  c�lculos complementares  -----

V = abs(E);
ang = angle(E);

Pg(bf) = Pd(bf) + P(bf);
Qg(bf) = Qd(bf) + Q(bf);
Qg(bpv) = Qd(bpv) + Q(bpv);
for i = 1 : nb
  Qsh(i) = shunt(i) * V(i)^2;
end

Pgtot = sum(Pg);
Qgtot = sum(Qg);
Pdtot = sum(Pd);
Qdtot = sum(Qd);
Qshtot = sum(Qsh);

tipo(bge) = 1;
tipo(bf) = 2;
bpv = find(tipo == 1);
npv = length(bpv);
bpq = find(tipo == 0 | tipo == 3);
npq = length(bpq);

% ----- impress�o dos resultados 

ang = radeg * ang;
 
%----- impress�o dos resultados (barras) 

disp(' ')
fprintf( '\n==============================================================================');
fprintf( '\n|                             Resultados das barras                           |');
fprintf( '\n==============================================================================');
fprintf( '\n    barra   tens�o          gera��o              carga           shunt        ');
fprintf('\n #  tipo  V(pu)  ang(0)   Pg(MW)   Qg(Mvar)   Pd(MW)   Qd(Mvar)   (MVAr)         ');
fprintf('\n--- ----  -----  ------  --------  --------  --------  --------  --------      ');
for i = 1:nb
   fprintf('\n%3.0f %3.0f %7.3f  %6.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f\n',AuxBarr(i),tipo(i),V(i),ang(i), ...
      Pg(i)*Sbase,Qg(i)*Sbase,Pd(i)*Sbase,Qd(i)*Sbase,shunt(i)*Sbase);
end
fprintf('\n                          --------  --------  --------  --------    ');
fprintf( '\nTotal:                    %8.2f  %8.2f  %8.2f  %8.2f   ', Pgtot*Sbase, ... 
                                         Qgtot*Sbase, Pdtot*Sbase, Qdtot*Sbase);
fprintf('\n');
fprintf('\nPerda de pot�ncia ativa nas LTs = %8.2f (MW)',(Pgtot-Pdtot)*Sbase);
fprintf('\n');


%888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

%-------------------------------------------------------------------
% Metodo da continuacao para maximo carregamento

% CASO BASE
% e0=e;     %solucao base
% f0=f;     %nao borrar
% Pd0=Pd;
% Qd0=Qd;
% 
%                        
% e=e0;
% f=f0;
% Pd=Pd0;
% Qd=Qd0;
% 
% k=0

Vesp=V
k=0
rho=0

alpha=1
    
    % ------------------ * * * SOLUCAO PREDITA * * * ---------------------- 
    E=e+j*f;   
    I = Y_barra*E;
    Ia = real(I);
    Ib = imag(I);

    %----------------------------------------------------------------------
%-------------------------------------------------------------------
       %tipo 0= PQ
       %tipo 1= PV
       %tipo 2= folga
   
    v1=tipo;                      
    v2=tipo+10*ones(nb,1);        
    v3=tipo+20*ones(nb,1);       
    v=[v1; v2; v3];
    
    w1=find(v==2);          % Ubica�ao das Barra Slack que preciso quitar no Jacobiano
    w2=find(v==11 | v==12); % Ubica�ao das Barras PV e Slack que preciso quitar no Jacobiano
    w3=find(v==20 | v==22); % Ubica�ao das Barras PQ e Slack que preciso quitar no Jacobiano
    w=[w1; w2; w3];         % Vetor com as filas que debem ser eliminadas
    
    z1=find(v==2)  ;        % Barra Slack
    z2=find(v==12) ;        % Barra Slack
    z=[z1;z2];              % Vetor com as colunas que debem ser eliminadas
  %------------------------------------------------------------------- 
    %----------------------------------------------------------------------
    
    % Caso base com rho=0 -------------------------------------------------
    %rho=i-1;
    EFrho=[e;f;rho];
    EFrho(z,:)=[];
    %----------------------------------------------------------------------
   
    % Construcao da matriz Jacobiana Predita ------------------------------
    J1 =  diag(e)*G + diag(f)*B + diag(Ia);
    J2 =  diag(f)*G - diag(e)*B + diag(Ib);
    J3 =  diag(f)*G - diag(e)*B - diag(Ib);
    J4 = -diag(e)*G - diag(f)*B + diag(Ia);
    J5 =  diag(2*e);
    J6 =  diag(2*f);
    
    Jp=-[J1 J2;J3 J4;J5 J6];
    Jp(w,:)=[];              
    Jp(:,z)=[];  
   
    %----------------------
    
    dP_dp=-0.01*Pd; 
    dQ_dp=-0.01*Qd;
    dV_dp=zeros(1,nb)';
    dPQV_p=[dP_dp;dQ_dp; dV_dp];
    dPQV_p(w,:)=[];  %  
 
    Jp=[Jp dPQV_p];
    dimJp=size(Jp);
    u0=zeros(1,dimJp(1));
    u0=[u0 1];
    Jp=[Jp;u0]; % Jacobiano predita
    %----------------------------------------------------------------------
    
    % Vetor tangente da predicao ------------------------------------------
    betha=1; 
    b=zeros(1,dimJp(1));
    b=[b betha]';
   
    delta_EFrho_p=inv(Jp)*b;  
    %----------------------------------------------------------------------
    
    % Solucao Predita -----------------------------------------------------
    n=length(delta_EFrho_p);
    

    EFrho=EFrho+alpha*delta_EFrho_p;
    
    e=[e(1); EFrho(1:(0.5*(n-1)))];
    f=[f(1); EFrho((0.5*(n-1)+1):(n-1))];
    rho=EFrho(n);
    
    %----------------------------------------------------------------------
    
    % Conver��o de tensao rectangular - polar --------------------------
    for i =1:nb
        V_rect(i) = complex(e(i),f(i));
        V_mod_p(i) = sqrt((e(i))^2+(f(i))^2);
        V_ang(i) = rad2deg(angle(V_rect(i)));
    end
    V_mod_p
    %----------------------------------------------------------------------
    
    % Calculos de Potencia injectadas -------------------------------------
    Pef=[diag(e) diag(f)]*[G -B;B G]*[e;f];
    Qef=[-diag(e) diag(f)]*[B G;G -B]*[e;f]; %igual que arriba 
    Vef=[diag(e) diag(f)]*[e;f];
    %----------------------------------------------------------------------
    
        
    
%##########################################################################    
                                                      
        % ---------------  SOLUCAO CORRIGIDA  ------------------- 
 
        tol=1;
        while tol>0.001
             tol=0.01;
 
        E=e+j*f;
        I = Y_barra*E;
        Ia = real(I);
        Ib = imag(I);   
        
        nPV = find(tipo == 1);
        %------------------------------------------------------------------
        %-------------------------------------------------------------------
       %tipo 0= PQ
       %tipo 1= PV
       %tipo 2= folga
    v1=tipo;                      
    v2=tipo+10*ones(nb,1);        
    v3=tipo+20*ones(nb,1); 
    v=[v1; v2; v3];
    
    w1=find(v==2);          
    w2=find(v==11 | v==12); 
    w3=find(v==20 | v==22); 
    w=[w1; w2; w3];         %  filas que devem ser eliminadas
    
    z1=find(v==2)  ;        % Barra Slack
    z2=find(v==12) ;        % Barra Slack
    z=[z1;z2];              % colunas que devem ser eliminadas
  %------------------------------------------------------------------- 

        %------------------------------------------------------------------
        
        
        % Matriz Jacobiana Corregida --------------------------------------
        J1 =  diag(e)*G + diag(f)*B + diag(Ia);
        J2 =  diag(f)*G - diag(e)*B + diag(Ib);
        J3 =  diag(f)*G - diag(e)*B - diag(Ib);
        J4 = -diag(e)*G - diag(f)*B + diag(Ia);
        J5 =  diag(2*e);
        J6 =  diag(2*f);

        Jc=-[J1 J2;J3 J4;J5 J6];
        Jc(w,:)=[];              
        Jc(:,z)=[];  
      
        
        Jc=[Jc dPQV_p];   
        Jc=[Jc;delta_EFrho_p'];   
        %------------------------------------------------------------------

        % Matriz de desbalanco de Potencia ---------------------------------
        deltaP=Pg-(1+0.01*rho)*Pd-Pef;  
        deltaQ=Qg-(1+0.01*rho)*Qd-Qef;
        deltaV=Vesp.^2-Vef;
        deltaPQV=[deltaP; deltaQ; deltaV];
        deltaPQV=[deltaPQV;0];
        deltaPQV(w,:)=[];
        %------------------------------------------------------------------

        % Criterio de parada do metodo de Newton / Correcao ---------------
         tol=max(abs(deltaPQV))
        %------------------------------------------------------------------


        % Vetor tangente da correcao --------------------------------------
        delta_EFrho_c=-inv(Jc)*deltaPQV;
        %------------------------------------------------------------------

        % Solucao corregida -----------------------------------------------
        n=length(delta_EFrho_c);
        alpha=1;
        EFrho=EFrho+alpha*delta_EFrho_c;

        e=[e(1); EFrho(1:(0.5*(n-1)))];
        f=[f(1); EFrho((0.5*(n-1)+1):(n-1))];
        rho=EFrho(n);
        %------------------------------------------------------------------

        % Conver��o de tensiones rectangular - polar ----------------------
        for i =1:nb
            V_rect(i) = complex(e(i),f(i));
            V_mod_c(i) = sqrt((e(i))^2+(f(i))^2);
            V_ang(i) = rad2deg(angle(V_rect(i)));
        end
        
        %------------------------------------------------------------------

        
        % Calculos de Potencia injectadas ---------------------------------
        Pef=[diag(e) diag(f)]*[G -B;B G]*[e;f];
        Qef=[-diag(e) diag(f)]*[B G;G -B]*[e;f]; 
        Vef=[diag(e) diag(f)]*[e;f];
        %------------------------------------------------------------------

        % criterio para mudanca de uma barra PV a uma PQ ------------------
        Q1=Qef+Qd-Qgmax;
        Q2=Qef+Qd-Qgmin;  
        for i=1:length(nPV)             
               if Q1(nPV(i))>0        
               Qg(nPV(i))=Qgmax(nPV(i));
               tipo(nPV(i))=0
               end
               if Q2(nPV(i))<0 
               Qg(nPV(i))=Qgmin(nPV(i));
               tipo(nPV(i))=0
               end   
        end
        %------------------------------------------------------------------

        end   %----while final

for iter = 1 : maxit
    
    V = abs(E);
    ang = angle(E);


   %-----  c�lculo das inje��es de pot�ncia -----

   E = e + j * f;
   I = Y_barra * E;
   Ia = real(I);
   Ib = imag(I);
   S = E .* conj(I);
   P = real(S);
   Q = imag(S);
   A=0.01;
   %-----  verifica��o dos limites de pot�ncia reativa gerada -----

   if (npv > 0)
      Qg(bpv) = Q(bpv) + Qd(bpv);
      if (iter > 1)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%     Subrotina para verifica��o dos limites de pot�ncia reativa gerada  %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if (iter >= 2) & (npv > 0)
      flag1 = 0;   
      for i = 1 : npv
         if (flag1 == 0)   
            if Qg(bpv(i)) > 0.999*Qgmax(bpv(i))
               Qg(bpv(i)) = Qgmax(bpv(i));
               flag1 = bpv(i);
               tipo(flag1) = 0;
               disp(' **********************************************************************')
               fprintf(' * Limite m�ximo de pot�ncia reativa gerada atingido na barra * %d\n',flag1)
               disp(' **********************************************************************')
            end
            if Qg(bpv(i)) < 1.001*Qgmin(bpv(i))
               Qg(bpv(i)) = Qgmin(bpv(i));
               flag1 = bpv(i);
               tipo(flag1) = 0;
               disp(' **********************************************************************')
               fprintf(' * Limite m�nimo de pot�ncia reativa gerada atingido na barra * %d\n',flag1)
               disp(' **********************************************************************')
            end
         end
      end
   
      if flag1 > 0
         bpv = find(tipo == 1);
         npv = length(bpv);
         bpq = find(tipo == 0 | tipo == 3);
         npq = length(bpq);
      end
   end


      end
   end
   
%-----  c�lculo dos desbalan�os de pot�ncia  -----

  dP = P(bsf) - Pg(bsf) + (Pd(bsf) + rho*A*Pd(bsf));   %altera��o 
  dQ = Q(bpq) - Qg(bpq) + (Qd(bpq) + rho*A*Qd(bpq));   %altera��o
  
  if (npv > 0)
     dV = (e(bpv).^2 + f(bpv).^2) - V(bpv).^2;  
     dsa = [ dP; dQ; dV ];
  else
     dsa = [ dP; dQ ];
  end
   
  
%-----  teste de converg�ncia -----  
  
   tst = max(abs(dsa));
   mdx = norm(mi*delta);
   %fprintf(' %8.4d     %8.4d     %8.4d     %8.4d\n',iter,tst,mi,mdx);
   if (tst < tol_res | mi < tol_mi)
      break
   end

%-----  determina��o da matriz jacobiana -----

   J1 =  diag(e)*G + diag(f)*B + diag(Ia);
   J2 =  diag(f)*G - diag(e)*B + diag(Ib);
   J3 =  diag(f)*G - diag(e)*B - diag(Ib);
   J4 = -diag(e)*G - diag(f)*B + diag(Ia);
   J5 = diag(2*e);
   J6 = diag(2*f);

   if (npv > 0) 
      J = [  J1(bsf, bsf)   J2(bsf, bsf)
             J3(bpq, bsf)   J4(bpq, bsf)
             J5(bpv, bsf)   J6(bpv, bsf)  ];
   else
      J = [  J1(bsf, bsf)   J2(bsf, bsf)
             J3(bpq, bsf)   J4(bpq, bsf)  ];
   end
 %--------Defini��o da matriz expandida da corre��o --------    
   Jrho; %mesmo da predi��o 
  
   Mcorrigido = [J Jrho;  Vtang'];
   
%-----  solu��o do sistema linear  -----

   delta_corrigido = Mcorrigido \ -[dsa; 0] ; 


%-----  c�lculo do fator de passo para evitar a diverg�ncia do processo iterativo -----

    aa = dsa;
    bb = - dsa;
  
    xx = zeros(nb, 1);
    xx(bsf) = delta_corrigido(1 : nb-1) + j * delta_corrigido(nb : 2*nb-2);
    yy = xx .* conj(Y_barra * xx);
  
    cc = [ real(yy(bsf))
    	   imag(yy(bpq))
    	   real(xx(bpv)).^2 + imag(xx(bpv)).^2 ];
    g0 = sum(aa .* bb);
    g1 = sum(bb .^ 2 + 2 * aa .* cc);
    g2 = 3 * sum(bb .* cc);
    g3 = 2 * sum(cc .^ 2);
  
    equacao = [ g3 g2 g1 g0 ];
    raizes = roots(equacao);
    valor1 = find(imag(raizes) == 0.0 & real(raizes) > 0.0);
    valor2 = find(imag(raizes) == 0.0 & real(raizes) < 0.0);
    
    if length(valor1) == 0
      mi = max(raizes(valor2));
    else
      mi = min(raizes(valor1));
    end
    
%-----  atualiza��o das vari�veis  -----

   if (mi > 1.0)
      mi = 1.0;
   end
  
  e(bsf) = e(bsf) + mi * delta_corrigido(1 : nb-1);
  f(bsf) = f(bsf) + mi * delta_corrigido(nb : 2*nb-2);
    
%    e(bsf) = e(bsf) + delta_corrigido(1 : nb-1);
%    f(bsf) = f(bsf) +  delta_corrigido(nb : 2*nb-2);   
end

% Pegando valores do caso base 
nbpv=find(tipo==1);
rho = 0; %valor inicial de rho 
alpha = 1; %valor do passo na predi��o
A = 0.01; %percentual da carga 
tol_j = 0.01; 

e0folga = e(find(tipo==2),:);  %e da folga (n�o altera) 
f0folga = f(find(tipo==2),:);  % f da folga  (n�o altera)

%% Predi��o 
% Retirando das submatrizes cada componentes de PV ,PQ e Folga 
  for inter1=1:1000
      
 % barra PQ tipo 0 e 3  
 % barra PV tipo 1 e 4 
     
 J;  
  
   %derivar em fun��o de rho 
   deltaP= -A * Pd;     % A=0.01%  
   deltaQ = -A *Qd;     % Q=0.01%
   deltaV=0*V;          % V=0
   
   
     deltaP = deltaP(bsf);           %elimina linhas das barras Folga
     deltaQ = deltaQ(bpq);           %elimina linhas das barras Folga e PV
     deltaV = deltaV(bpv);          %elimina linhas das barras Folga e PQ
    
   % calculo do vetor tangente predi��o
   Jrho=[deltaP; deltaQ; deltaV]; 
   Vzero=zeros(1,length(J)); 
   
   Jtang=[J Jrho; Vzero 1];    
   Vtang = inv(Jtang)*[Vzero 1]';
         
   rho                     % valor de rho 
   rhoplot(inter1) =rho;   % rho para plotar 
   Normal_Vtang(inter1) = norm(Vtang); 
   
   %----criterio de parada--------------
   minJ=svd(J);               % valores singulares de J 
   min_J(inter1)=min(minJ);    % minimo valor singular 
   if min(minJ) < tol_j
               break
   end
     
%  Novo Ponto predito: 
e(find(tipo==2),:)=[]; %retira ponto barra de folga
f(find(tipo==2),:)=[]; %retira ponto barra de folga 
   
[Vpredito] = [e;f;rho] + alpha*[Vtang];         %pontos preditos (SOLU��O)
   
e = Vpredito([1:length(e)]);                     %retirando e de Vpredito
f = Vpredito([length(e)+1:length(Vpredito)-1]);  %retirando f de Vpredito
rho = Vpredito(length(Vpredito));                %retirando rho de Vpredito 

%acrescentando e0folga e f0folga para o fluxo 
e=[e0folga ; e]; 
f=[f0folga ; f]; 


%% Novo Fluxo para Corre��o 


for iter = 1 : maxit
    
    V = abs(E);
    ang = angle(E);


   %-----  c�lculo das inje��es de pot�ncia -----

   E = e + j * f;
   I = Y_barra * E;
   Ia = real(I);
   Ib = imag(I);
   S = E .* conj(I);
   P = real(S);
   Q = imag(S);

   %-----  verifica��o dos limites de pot�ncia reativa gerada -----

   if (npv > 0)
      Qg(bpv) = Q(bpv) + Qd(bpv);
      if (iter > 1)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%     Subrotina para verifica��o dos limites de pot�ncia reativa gerada  %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if (iter >= 2) & (npv > 0)
      flag1 = 0;   
      for i = 1 : npv
         if (flag1 == 0)   
            if Qg(bpv(i)) > 0.999*Qgmax(bpv(i))
               Qg(bpv(i)) = Qgmax(bpv(i));
               flag1 = bpv(i);
               tipo(flag1) = 0;
               disp(' **********************************************************************')
               fprintf(' * Limite m�ximo de pot�ncia reativa gerada atingido na barra * %d\n',flag1)
               disp(' **********************************************************************')
            end
            if Qg(bpv(i)) < 1.001*Qgmin(bpv(i))
               Qg(bpv(i)) = Qgmin(bpv(i));
               flag1 = bpv(i);
               tipo(flag1) = 0;
               disp(' **********************************************************************')
               fprintf(' * Limite m�nimo de pot�ncia reativa gerada atingido na barra * %d\n',flag1)
               disp(' **********************************************************************')
            end
         end
      end
   
      if flag1 > 0
         bpv = find(tipo == 1);
         npv = length(bpv);
         bpq = find(tipo == 0 | tipo == 3);
         npq = length(bpq);
      end
   end


      end
   end
   
%-----  c�lculo dos desbalan�os de pot�ncia  -----

  dP = P(bsf) - Pg(bsf) + (Pd(bsf) + rho*A*Pd(bsf));   %altera��o 
  dQ = Q(bpq) - Qg(bpq) + (Qd(bpq) + rho*A*Qd(bpq));   %altera��o
  
  if (npv > 0)
     dV = (e(bpv).^2 + f(bpv).^2) - V(bpv).^2;  
     dsa = [ dP; dQ; dV ];
  else
     dsa = [ dP; dQ ];
  end
   
  
%-----  teste de converg�ncia -----  
  
   tst = max(abs(dsa));
   mdx = norm(mi*delta);
   %fprintf(' %8.4d     %8.4d     %8.4d     %8.4d\n',iter,tst,mi,mdx);
   if (tst < tol_res | mi < tol_mi)
      break
   end

%-----  determina��o da matriz jacobiana -----

   J1 =  diag(e)*G + diag(f)*B + diag(Ia);
   J2 =  diag(f)*G - diag(e)*B + diag(Ib);
   J3 =  diag(f)*G - diag(e)*B - diag(Ib);
   J4 = -diag(e)*G - diag(f)*B + diag(Ia);
   J5 = diag(2*e);
   J6 = diag(2*f);

   if (npv > 0) 
      J = [  J1(bsf, bsf)   J2(bsf, bsf)
             J3(bpq, bsf)   J4(bpq, bsf)
             J5(bpv, bsf)   J6(bpv, bsf)  ];
   else
      J = [  J1(bsf, bsf)   J2(bsf, bsf)
             J3(bpq, bsf)   J4(bpq, bsf)  ];
   end
 %--------Defini��o da matriz expandida da corre��o --------    
   Jrho; %mesmo da predi��o 
  
   Mcorrigido = [J Jrho;  Vtang'];
   
%-----  solu��o do sistema linear  -----

   delta_corrigido = Mcorrigido \ -[dsa; 0] ; 

%-----  c�lculo do fator de passo para evitar a diverg�ncia do processo iterativo -----

    aa = dsa;
    bb = - dsa;
  
    xx = zeros(nb, 1);
    xx(bsf) = delta_corrigido(1 : nb-1) + j * delta_corrigido(nb : 2*nb-2);
    yy = xx .* conj(Y_barra * xx);
  
    cc = [ real(yy(bsf))
    	   imag(yy(bpq))
    	   real(xx(bpv)).^2 + imag(xx(bpv)).^2 ];
    g0 = sum(aa .* bb);
    g1 = sum(bb .^ 2 + 2 * aa .* cc);
    g2 = 3 * sum(bb .* cc);
    g3 = 2 * sum(cc .^ 2);
  
    equacao = [ g3 g2 g1 g0 ];
    raizes = roots(equacao);
    valor1 = find(imag(raizes) == 0.0 & real(raizes) > 0.0);
    valor2 = find(imag(raizes) == 0.0 & real(raizes) < 0.0);
    
    if length(valor1) == 0
      mi = max(raizes(valor2));
    else
      mi = min(raizes(valor1));
    end
    
%-----  atualiza��o das vari�veis  -----

   if (mi > 1.0)
      mi = 1.0;
   end
  
  e(bsf) = e(bsf) + mi * delta_corrigido(1 : nb-1);
  f(bsf) = f(bsf) + mi * delta_corrigido(nb : 2*nb-2);
    
%    e(bsf) = e(bsf) + delta_corrigido(1 : nb-1);
%    f(bsf) = f(bsf) +  delta_corrigido(nb : 2*nb-2);

     
end
 
% Resulta em 'delta_corrigido'

%% Altera��o do passo alpha 
E = e + j * f;       
Vmod = abs(E);           % Modulo das tens�es 
Angulo = angle(E);        % Angulo das tens�es

% Valores para curva PxV
Pcpv1(inter1) = -P(107); 
Vcpv1(inter1) = V(106); 

Pcpv2(inter1) = -P(107); 
Vcpv2(inter1) = V(106); 

Pcpv3(inter1) = -P(107); 
Vcpv3(inter1) = V(106); 
 


%Para verificar quando Arho muda de sinal
arho(inter1)=delta_corrigido(length(delta_corrigido));
testedelta(inter1) = arho(inter1) ;

  end
  %----- impress�o dos resultados (barras) 

disp(' ')
fprintf( '\n==============================================================================');
fprintf( '\n|                   Resultados das barras Ponto Critico                      |');
fprintf( '\n==============================================================================');
fprintf( '\n    barra   tens�o          gera��o              carga           shunt        ');
fprintf('\n #  tipo  V(pu)  ang(0)   Pg(MW)   Qg(Mvar)   Pd(MW)   Qd(Mvar)   (MVAr)         ');
fprintf('\n--- ----  -----  ------  --------  --------  --------  --------  --------      ');
for i = 1:nb
   fprintf('\n%3.0f %3.0f %7.3f  %6.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f\n',AuxBarr(i),tipo(i),V(i),ang(i), ...
      Pg(i)*Sbase,Qg(i)*Sbase,Pd(i)*Sbase,Qd(i)*Sbase,shunt(i)*Sbase);
end
fprintf('\n                          --------  --------  --------  --------    ');
fprintf( '\nTotal:                    %8.2f  %8.2f  %8.2f  %8.2f   ', Pgtot*Sbase, ... 
                                         Qgtot*Sbase, Pdtot*Sbase, Qdtot*Sbase);
fprintf('\n');
fprintf('\nPerda de pot�ncia ativa nas LTs = %8.2f (MW)',(Pgtot-Pdtot)*Sbase);
fprintf('\n');

for m=1:maxit
for m=1:maxit
for m=1:maxit
Pgm(1:nb)=Pg + rho*m*A*Pg;
Pdm(1:nb)=Pd + rho*m*A*Pd;
Qdm(1:nb)=Qd + rho*m*A*Qd;

Qgm=Qg;
Pdm(Pdm==0)=1;
Pgm(Pgm==0)=1;
angDM=atan(Qdm'./Pdm');
Rd=sqrt(Qdm'.^2+Pdm'.^2);
angGM=atan(Qgm./Pgm');
Rm=sqrt(Pgm'.^2+Qgm.^2);
PgM=Pgm'.*cos(angGM);
PdM=Pdm'.*cos(angDM);
QdM=Qdm'.*sin(angDM);
QgM=Qgm.*sin(angGM);
PgM(3:6)=0;
PgM(1)<= Pgmax(1);
if PgM(1) >= Pgmax(1) -0.05
    break
end
PgM(2)<= Pgmax(2);
if PgM(2) >= Pgmax(2) -0.05
    break
end
 PdM <=  PgM;
 if PdM >=  PgM;
     break
 end
end
end


Qgtot=sum(QgM);
Pgtot=sum(PgM);
Qdtot=sum(QdM);
Pdtot=sum(PdM);
end

rho=rho;
fprintf('Carregamento m�ximo: %4.2f\n',rho);
fprintf('%d itera��es\nmi de %1.3f\n',iter,tol_mi)  
disp('===================================================================')

disp('===================================================================')


fprintf( '\n==============================================================================');
fprintf( '\n|               Resultados das barras com Carregamento M�ximo                |');
fprintf( '\n==============================================================================');
fprintf( '\n    barra   tens�o          gera��o              carga           shunt        ');
fprintf('\n #  tipo  V(pu)  ang(0)   PgM(MW)   QgM(Mvar)   PdM(MW)   QdM(Mvar)   (MVAr)       ');
fprintf('\n--- ----  -----  ------  --------  --------  --------  --------  --------      ');
for i = 1:nb
   fprintf('\n%3.0f %3.0f %7.3f  %6.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f\n',AuxBarr(i),tipo(i),V(i),ang(i), ...
      PgM(i)*Sbase,QgM(i)*Sbase,PdM(i)*Sbase,QdM(i)*Sbase,shunt(i)*Sbase);
end

fprintf('\n                          --------  --------  --------  --------    ');
fprintf( '\nTotal:                  %8.2f   %8.2f  %8.2f  %8.2f',Pgtot*Sbase, Qgtot*Sbase, Pdtot*Sbase, Qdtot*Sbase);