#
# alspgbox implements the augmented-lagrangian method on top of spgbox, to allow 
# the use of constraints
# 
# Implemented by J. M. Martínez, 2018 (In Fortran)
# To Julia: L. Martínez, 2021.
#


      implicit none
      integer nn, n, nrest, maxinner, maxouter, maxnef, ier, i, ipri
      integer npun, ciclo, exito
      parameter(nn=20000)
      double precision L(nn), u(nn), x(nn), lambda(nn), aux(4*nn)
      double precision epsfeas, epsopt, frac, factor, seed
      double precision hmax, hmed, save, f, fmin

      nrest = 1 


      epsfeas = 1.d-4
      epsopt = 1.d-4
      maxinner = 10000
      maxouter = 50 
      maxnef = 10000000
      ipri = 1
      frac = 0.5
      factor = 5


      call alspgbox(n, nrest, L, u, x, ier, lambda, epsfeas,
     * epsopt, maxinner, maxouter, maxnef, aux, frac, factor, ipri)

      write(*, *)' Solucao:' 
      do i = 1, n
      write(*, *) i, x(i)
      write(11,*) i, x(i)
      end do

 
      stop
      



      write(*, *)' puntos:'
      read(*, *) npun
      n = 2*npun
      do i = 1, npun
      L(2*i-1) = -100.d0
      u(2*i-1) = 100.d0
      L(2*i) = 0.d0
      u(2*i) = 300 .d0
      end do
      write(*, *)' Points:', npun
      write(*, *)' horizontal bounds:', L(1), u(1)

      seed = 111298277

      nrest = 0

      do i = 1, n
      call rando(seed, x(i))
      x(i) = L(i) + x(i)*(u(i)- L(i))
      end do


      epsfeas = 1.d-4
      epsopt = 0.01
      maxinner = 1000
      maxouter = 1 
      maxnef = 10000000
      ipri = 2

      frac = 0.5
      factor = 5

      write(*, *)' Tratando de bajar inicialmente:'
      call fun(n, x, fmin)
      ciclo = 1
2     exito = 0
      write(*, *)' Ciclo :', ciclo, ' fmin =', fmin
      do i = 1, npun
      save = x(2*i)
      x(2*i) = x(2*i)-1.
      call fun(n, x, f)
      if(f.ge.fmin.or.x(2*i).lt.0.d0) then
      x(2*i) = save
      else
c      write(*, *)' i = ', i,' f anterior =', fmin, ' f new = ', f
      exito = 1
      fmin = f
      endif
      end do
      if(exito.eq.0) go to 3
      ciclo = ciclo+1
      go to 2

3     write(*, *)' ciclos = ' , ciclo
      call fun(n, x, f)
      write(*, *)' f despues del primer ciclo:', f
 
 


      call alspgbox(n, nrest, L, u, x, ier, lambda, epsfeas,
     * epsopt, maxinner, maxouter, maxnef, aux, frac, factor, ipri)

      call fun(n, x, f)
      write(*, *)' f antes del ciclo:', f
 

      write(*, *)' Tratando de bajar mas:'
      call fun(n, x, fmin)
      ciclo =  1
4     exito = 0 
      write(*, *)' Ciclo :', ciclo
      do i = 1, npun
      save = x(2*i)
      x(2*i) = x(2*i)-1.
      call fun(n, x, f)
      if(f.ge.fmin.or.x(2*i).lt.0.d0) then
      x(2*i) = save
      else
      write(*, *)' i = ', i,' f anterior =', fmin, ' f new = ', f
      fmin = f
      exito = 1
      endif
      end do
      if(exito.eq.0) go to 5
      ciclo=ciclo+1
      go to 4


5      call fun(n, x, f)
      write(*, *)' f despues del ciclo:', f
 

 
      write(*, *)' Puntos salida de alspgbox :'
      hmax = 0.d0
      hmed = 0.d0
      do i = 1, npun
      write(*, *) i,  x(2*i-1), x(2*i)
      write(7, *) i, x(2*i-1), x(2*i)
      hmax = dmax1(hmax, x(2*i))
      hmed = hmed + x(2*i)
      end do
      hmed = hmed/dfloat(npun)
      write(*, *)' Altura maxima:', hmax,' Altura media: ', hmed
      write(*, *)' f menor:', fmin 



      stop           


      n = 3

      nrest = 1 

      do i = 1, n - 1  
      L(i) = -10.d0
      u(i) = 10.d0
      end do

      L(n) = 0.d0
      u(n) = 1.d3

      x(1) = -1.2d0
      x(2) = 1.d0
      x(3) = 0.d0




      epsfeas = 1.d-4
      epsopt = 1.d-4
      maxinner = 100000
      maxouter = 50 
      maxnef = 10000000
      ipri = 0

      frac = 0.5
      factor = 5


      call alspgbox(n, nrest, L, u, x, ier, lambda, epsfeas,
     * epsopt, maxinner, maxouter, maxnef, aux, frac, factor, ipri)

      stop
      end

      subroutine fun(n, x, f)
      implicit none
      integer n, i, j, npun
      double precision x(n), f, xi(2), xj(2), ener, sigma

      sigma = 10.

      f = 0.d0
      do i = 1, n
      f = f + x(i)
      end do
      f = sigma*f

      do i = 1, n-1
      f = f + (x(i+1)-x(i))**2
      end do
      
      return



      sigma = 10000

      f = 0.d0
      do i =  1, n
      f = f + x(i)**3
      end do

      do i = 2, n-1
      f = f + (x(i-1) - 2*x(i) + x(i+1))**2
      end do


      do i = 1, 449
      f = f + sigma*x(i)
      end do
      do i = 450, 550
      f = f + sigma*(200.d0 - x(i)) 
      end do
      do i = 551, n
      f = f + sigma*x(i)
      end do
 
      return
 


      sigma = 0.01
      f = 0.d0
      do i =  2, n-1
      f = f + (x(i-1) - x(i+1))**2
      end do
      do i = 1, 449
      f = f + sigma*x(i)
      end do
      do i = 450, 550
      f = f + sigma*(200.d0 - x(i)) 
      end do
      do i = 551, 1000
      f = f + sigma*x(i)
      end do
 
      return
      

      npun = n/2
      f = 0.d0
      do i = 1, npun
      f = f + x(2*i)
      end do
       
      do j = 1, npun-1
      do i = j+1, npun
      xi(1) = x(2*i-1)
      xi(2) = x(2*i)
      xj(1) = x(2*j-1)
      xj(2) = x(2*j)
      call energia(xi, xj, ener)
      f = f + ener
      end do
      end do

      
      return


      f = 0.d0
      do i = 1, n/2
      f = f + x(2*i)
      end do

      return
 



      f = 1.d2 * (x(2)-x(1)**2)**2 + (x(1) - 1.d0)**2

      f = 1.d0 * (x(2)-x(1)**2)**2 + (x(1) - 1.d0)**2
 

      return
      end

      subroutine energia(x, y, ener)
      implicit none
      double precision x(2), y(2), ener, d2, d
      d2 = (x(1) - y(1))**2 + (x(2) - y(2))**2
      d = dsqrt(d2) 
      call fdede(d, ener)
      return
      end

      subroutine fdede(d, ener)
      implicit none
      double precision d, ener
      ener = 1.d0
      if(d.le.1.d0) then
      ener = 1.d4 - 1.d4*d
      else
      if(d.le.2.d0) then
      ener = d - 1.d0
      endif
      endif
      return
      end

      subroutine dedede(d, deri)
      implicit none
      double precision d, deri
      deri = 0.d0
      if(d.le.1.d0) then
      deri = -1.d4
      else
      if(d.le.2.d0) then
      deri = 1.d0
      endif
      endif
      return
      end 


      
      

      subroutine graenergia(x, y, g)  
      implicit none
      integer i
      double precision x(2), y(2), g(4), d2, d, deri



      d2 = (x(1)-y(1))**2 + (x(2)-y(2))**2
      d = dsqrt(d2)
      if(d.eq.0.d0)d=1.d-8
      g(1) = (x(1)-y(1))/d
      g(2) = (x(2)-y(2))/d
      g(3) = -g(1)
      g(4) = -g(2)
      call dedede(d, deri)
      do i = 1, 4
      g(i) = deri*g(i)
      end do 


      return
      end 
      



      subroutine grad(n, x, g)
      implicit none
      integer n, i, j, npun
      double precision x(n), g(n), h, save, fmas, fmen
      double precision f, xi(2), xj(2), gij(4), sigma

      sigma = 10.


      
      do i = 1, n
c      g(i) = 0.5*sigma/dsqrt(x(i))
      g(i) = sigma
      end do

      do i = 1, n-1
      g(i) = g(i) + 2.d0*(x(i)-x(i+1))
      g(i+1) = g(i+1) + 2.d0*(x(i+1)-x(i))
      end do

      return

      f = sigma*f

      do i = 1, n-1
      f = f + (x(i+1)-x(i))**2
      end do
      
      return

 


      sigma = 10000

      do i = 1, n
      g(i) = 0.d0
      end do


      do i = 1, n
      g(i) = 3.d0*x(i)**2
      end do
      do i = 2, n-1
c      f = f + (x(i-1) - 2*x(i) + x(i+1))**2
      g(i-1) = g(i-1) + 2.d0*(x(i-1) - 2.d0*x(i) + x(i+1))
      g(i) = g(i) - 4.d0*(x(i-1) - 2*x(i) + x(i+1))
      g(i+1) = g(i+1) + 2.d0*(x(i-1) - 2.d0*x(i) + x(i+1))     
      end do
 



      do i = 1, 449
      g(i) = g(i)+sigma
      end do
      do i = 450, 550
      g(i) = g(i)-sigma
      end do
      do i = 551, n
      g(i) = g(i)+sigma 
      end do
 
      return
 
 


c      f = 0.d0
      do i = 1, n
      g(i) = 0.d0
      end do

      do i =  2, n-1
c      f = f + (x(i-1) - x(i+1))**2
      g(i-1) = g(i-1) + 2.d0*(x(i-1)-x(i+1))
      g(i+1) = g(i+1) + 2.d0*(x(i+1)-x(i-1))
      end do


      do i = 1, 449
      g(i) = g(i)+sigma
      end do
      do i = 450, 550
      f = f + sigma*dabs(x(i)-200.d0)
      g(i) = g(i)-sigma
      end do
      do i = 551, 1000
      g(i) = g(i)+sigma 
      end do
 
      return
 





c      go to 1

      npun = n/2
      do i = 1, npun
      g(2*i) = 1.d0
      end do
      do j = 1, npun-1
      do i = j+1, npun
      xi(1) = x(2*i-1)
      xi(2) = x(2*i)
      xj(1) = x(2*j-1)
      xj(2) = x(2*j)
      call graenergia(xi, xj, gij)
      g(2*i-1) = g(2*i-1) + gij(1)
      g(2*i) = g(2*i) + gij(2)
      g(2*j-1) = g(2*j-1) + gij(3)
      g(2*j) = g(2*j) + gij(4) 
      end do
      end do   

c      write(*, *)' gradiente continuo:' 
c      write(*, *)(g(i),i=1,n)
c      stop        

      return

      do i = 1, n/2
      g(2*i) = 1.d0
      g(2*i-1) = 0.d0
      end do

      return
 


      g(1) = -4.d2*x(1)*(x(2)-x(1)**2) + 2.d0*(x(1)-1.d0)
      g(2) = 200.d0*(x(2) - x(1)**2)


1      do i = 1, n
      save = x(i)
      h = 1.d-9 * dabs( save )
      if(h.lt.1.d-12) h = 1.d-12
      x(i) = save + h
      call fun(n, x, fmas)
      x(i) = save - h
      call fun(n, x, fmen)
      x(i) = save
      g(i) = (fmas - fmen)/(2.d0*h)
      end do

      write(*, *)' gradiente discreto:'
      write(*, *)(g(i),i=1,n)

      stop
      

 
      return
      end

      subroutine haga(i, n, x, h)
      implicit none
      integer i, n, npun, ii, j, k
      double precision x(n), h, z

      h = 0.d0
      do j = 1, n
      h = h + x(j)
      end do
      h = h - 101.*200.
      return

      npun = n/2
      k = 1
      do j = 1, npun-1
      do ii = j+1, npun
      if(k.eq.i) then
      z = (x(2*ii-1) - x(2*j-1))**2 +   (x(2*ii) - x(2*j))**2 
c        if (z.ge.1.d0) then
c        h = 0.d0
c        else
c        h = (z - 1)**2
      h = 100.d0/dsqrt(z)
c        endif
      return   
      endif
      k = k+1
      end do
      end do

c
c     
 



      h = x(1)**2 + x(2)**2 + x(3) -  1.d0
      return
      end

      subroutine jacob (i, n, x, g)
      implicit none
      integer i, j, n
      double precision x(n), g(n), h, save, hmas, hmen


      do j = 1, n
      g(j) = 1.d0
      end do
      return

c      g(1) = 2.d0*x(1) 
c      g(2) = 2.d0*x(2) 
c      g(3) = 1.d0


      do j  = 1, n
      save = x(j)
      h = 1.d-6 * dabs( save )
      if(h.lt.1.d-12) h = 1.d-12
      x(j) = save + h
      call haga(i, n, x, hmas)
      x(j) = save - h
      call haga(i, n, x, hmen)
      x(j) = save
      g(j) = (hmas - hmen)/(2.d0*h)
      end do
 


      return
      end
 

 


 

c************************************************************************
c************************************************************************
c************************************************************************      


c  Lagrangiano Aumentado para minimizar f(x) sujeita a restricoes de 
c  igualdade e caixas.
c  Para resolver os subproblemas usa spgbox (SPG en caixas)
c  Escrito com objetivos didaticos (MS 628) por:
c  J. M. Martínez, Setembro 2018.



      subroutine alspgbox(n, nrest, L, u, x, ier, lambda, epsfeas, 
     *  epsopt, maxinner, maxouter, maxnef, aux, frac, factor, ipri)
  
      implicit none
      integer n, outer, nef, inner, nrest, i, ier, maxinner, maxouter,
     *  maxnef, konmax, nafmax, kon, naf, ipri
      double precision x(n), L(n), u(n), f, hnor, h, epsfeas, rho, 
     *  lambda(nrest), faux, gnoraux, numer, deno, epsopt, aux(4*n),
     *  frac, factor, hnoran

      konmax = dfloat(maxinner)/2.d0  
      nafmax = dfloat(maxnef)/2.d0

      if(nrest.eq.0) then
      konmax = maxinner
      nafmax = maxnef
      endif
 

      outer = 0
      nef = 0
      inner = 0
      ier = 0

      do i = 1, nrest
      lambda(i) = 0.d0
      end do

      rho = 1.

10    write(*, *)
      write(*, *)' Outer Iteration:', outer
      write(*, *)' X = ', x(1), x(2), ', ... ,', x(n)
      call fun(n, x, f)
      write(*, *)' f(X) = ', f
      hnor = 0.d0
      do i = 1, nrest
      call haga(i, n, x, h)
      hnor = dmax1(hnor, dabs(h))
      end do
      write(*, *)' ||h(X)||_sup = ', hnor
      write(*, *)' Inner iterations:', inner
      write(*, *)' Function evaluations:', nef
      if(inner.gt.0) then
      write(*, *)' Outer obtained with rho = ', rho
      endif
 
      if(hnor.le.epsfeas) then
      write(*, *)' Constraints are satisfied'
      if(outer.gt.0) return
      endif

      if(outer.eq.0) then
      if(dabs(f).lt.1.d0) then
      deno = 1.d0
      else
      deno = dabs(f)
      endif
      if(hnor*hnor.gt.1.d2) then
      numer = 1.d2
      else
      numer = hnor*hnor
      endif
      rho = dmax1(rho,  numer/deno)
      else
      if(hnor.gt.frac*hnoran) rho = factor*rho
      endif

      hnoran = hnor
      


      if(outer.ge.maxouter) then
      write(*, *)' Outer iterations exceeded'
      ier = 10 + ier
      return
      endif

      if(inner.ge.maxinner) then
      write(*, *)' Inner iterations exceeded'
      ier = 10 + ier
      return
      endif

      if(nef.ge.maxnef) then
      write(*, *)' Function evaluations exceeded'
      ier = 10 + ier
      return
      endif


      call spgboxal(n, x, aux, L, u, 10, konmax, nafmax, epsopt, 
     * faux, gnoraux, ier, kon, naf, aux(n+1), aux(2*n+1), nrest,
     * lambda, rho, aux(3*n+1), ipri) 

      inner = inner + kon
      nef = naf + nef
     
      do i = 1, nrest
      call haga(i, n, x, h)
      lambda(i) = lambda(i) + rho*h
      end do   
      outer = outer + 1
      go to 10

      end




      subroutine spgboxal(n, x, g, L, u,  m, konmax, nafmax, eps, f,
     *   gnor, ier, kon, naf, xn, gn, nrest, lambda, rho, aux, ipri) 
c  Programa principal para testar alspgbox

      implicit none
      double precision x(n), g(n), xn(n), gn(n), f, gnor, tspg, L(n),
     *   u(n), t, fn, p, num, den, fant(100), fref, eps, lambda(nrest),
     *   rho, aux(n), z
      integer n, m, kon, naf, i, ier, konmax, nafmax, nrest, ipri

      kon = 0
      call funal(n, x, f, nrest, lambda, rho)
      naf = 1
      call gradal(n, x, g, nrest, lambda, rho, aux)
10    gnor = 0.d0 
      do i =  1, n
      z = dmax1(L(i), dmin1(u(i), x(i) - g(i))) - x(i)
      gnor = dmax1(gnor, dabs(z))
      end do         


      if(ipri.ne.0) then
      write(*, *)
      write(*, *)' Inner iteration ', kon
      write(*, *)' X = ', x(1), ' ... ', x(n)
      write(*, *)' f(X) = ', f

      write(*, *)' Norma do gradiente projetado = ', gnor
      write(*, *)' Avaliacoes de AL(x) = ', naf
      endif


      if(gnor.le.eps) then
      ier = 0
      write(*, *)' spg stopped by small projected gradient' 
      return
      endif
      if(kon.ge.konmax) then
      ier =1
      write(*, *)' spg iterations exhausted, kon =', kon
      return
      endif
      if(naf.ge.nafmax) then
      ier =  2
      write(*, *)' spg evaluations exhausted, kon =', kon
 
      return
      endif
      if(kon.eq.0) then
      tspg = 1.d0
      do i = 1, m
      fant(i) = f
      end do
      endif

      write(*, *)' tspg = ', tspg

      t = tspg
      fref = 0.d0
      do i = 1, m
      fref = dmax1(fref, fant(i))
      end do

c      write(*, *)' fref = ', fref
c      write(*, *)' fant = ', (fant(i),i=1,m)

c      write(*, *)' t = ', t

20    do i = 1, n
      xn(i) = x(i) - t*g(i)
      xn(i) = dmax1(L(i), dmin1(xn(i),u(i)))
      end do

c      write(*, *)' xn =', xn(1), xn(2)

      call funal(n, xn, fn, nrest, lambda, rho)

c      write(*, *)' f(xn) = ', fn, ' fref = ', fref

      naf = naf + 1
      if(naf.gt.nafmax) then
      ier = 2
      write(*, *)' spg evaluations exhausted'
      return
      endif
      if(fn.le.fref) then
      call gradal(n, xn, gn, nrest, lambda, rho, aux)
      den = 0.d0
      num = 0.d0
      do i = 1, n
      num = num + (xn(i)-x(i))**2
      den = den + (xn(i)-x(i))*(gn(i)-g(i))

c      write(*, *)'i, xn(i), x(i), gn(i), g(i):', i, xn(i), x(i), 
c     *            gn(i), g(i), den
 
      end do

c      write(*, *)' den = ', den
      if(den.le.0.d0) then
      tspg = 1.d2
      else
          tspg = num/den
          if(tspg.ge.1.d3) then
          tspg = 1.d3
          endif
      endif

c      write(*, *)' num, den, tspg:', num, den, tspg
      f = fn
      do i = 1, n
      x(i) = xn(i)
      g(i) = gn(i)
      end do
      do i = 1, m-1
      fant(i) = fant(i+1)
      end do
      fant(m) = f
      kon = kon +1
      go to 10
      endif

      t = t/2.d0
      go to 20
      end

     
      subroutine funal(n, x, f, nrest, lambda, rho)
      implicit none
      integer n, nrest, j 
      double precision x(n), h, lambda(nrest), f, rho
      call fun(n, x, f)
      do j = 1, nrest
      call haga(j, n, x, h)
      f = f + lambda(j)*h + (rho/2.d0)*h**2
      end do
      return
      end

      subroutine gradal(n, x, g, nrest, lambda, rho, gh)
      implicit none
      integer n, nrest, i, j
      double preci*  sion x(n), g(n), lambda(nrest), rho, gh(n), h
      call grad(n, x, g)
      do j = 1, nrest
      call haga(j, n, x, h)
      call jacob(j, n, x, gh)
      do i = 1, n
      g(i) = g(i) + (lambda(j) + rho*h)*gh(i)
      end do
      end do
      return
      end

 

  

 



 
                       
      subroutine rando(seed, x)

C     This is the random number generator of Schrage:
C
C     L. Schrage, A more portable Fortran random number generator, ACM
C     Transactions on Mathematical Software 5 (1979), 132-138.

      double precision seed, x 

      double precision a,p,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/

      xhi= seed/b16
      xhi= xhi - dmod(xhi,1.d0)
      xalo= (seed-xhi*b16)*a
      leftlo= xalo/b16
      leftlo= leftlo - dmod(leftlo,1.d0)
      fhi= xhi*a + leftlo
      k= fhi/b15
      k= k - dmod(k,1.d0)
      seed= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (seed.lt.0) seed = seed + p
      x = seed*4.656612875d-10

      return

      end 
             
 
  

 




