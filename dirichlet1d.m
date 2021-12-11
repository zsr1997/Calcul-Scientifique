function [err, soldir] = dirichlet1d(a, b)

  disp('Choisir un nombre de points de discretisation N = 9 ou 99 ou 999 ou ...')

  N = input('Effectuez votre choix : N = 9 ou 19 ou ... ou 999   ')
  disp('Vous avez choisi N = ')
  disp(N)

  % Determination du pas de discretisation 

  h = 1/ (N + 1) ;

  % Saisie de la matrice du systeme 

  A   = 2 * speye(N) ;
  Aup = sparse(1:(N-1), 2:N, -ones(1,N-1), N, N, N-1) ;
  A   = Aup' + A + Aup ;

  % ou encore
  % A = toeplitz([2, -1, zeros(1,N-2)]) ;

  clear Aup ;

  A = full(A) ;

  % Repartition des coefficients non nuls de la matrice du systeme

  disp('Visualisez la repartition des coefficients non nuls de la matrice')
  disp('Et saisir return pour la suite des calculs')
  spy(A)

  keyboard

  % Saisie du second membre du systeme

  % Executer avec la saisie suivante 
  % (en utilisant ainsi la formule du Trapeze)
  for i = 1:N 
      %B(i)  = h^2 * f((i-1/3)*h) ;
      %B(i)= h^2 * (feval('f',(i-(1/2))*h) + feval('f', i*h) +feval('f',(i+(1/2))*h))/3;
      B(i)= (h^2)*f(i*h);
      %B(i)=(h^2*(feval('f',(i-1/3)*h)+4*feval('f',(i-1/3)*h)+feval('f',(i+1/3)*h)))/6
  end

  % Reprendre en effectuant la saisie basee sur la formule de Simpson

  % Calcul de la solution de A X = B :

  X = A\B' ;

  % Consideration de la solution sur tout le domaine 

  soldir = [0 ; X ; 0] ;

  % Consideration de la solution exacte dans la base discrete 

  for i = 1:(N+2)
      u(i) = feval('uex', (i-1)*h) ;
  end

  erreur = norm(u' - soldir, inf)/norm(u', inf) ;
  err    = 100 * erreur ;

  disp('Erreur commise en % avec la methode de discretisation est de :')
  disp(err)

  % Representation graphique de la solution

  for i = 1:(N+2)
      x(i) = (i-1)*h ;
  end

  plot(x, soldir, 'bo')
  hold on
  fplot('uex', [0 1],'r-')
  legend('Solution discrete', 'Solution exacte', -1)
  xlabel('Abscisses')
  ylabel('Valeurs prises par la solution discrete et la solution exacte')
  title('Representation de la solution discrete et de la solution exacte')

  %Impression en postscript
  print dirichlet1d_9.ps

  hold off

end 

