c Computes air shower timing in a Direct manner,
c meaing not using the zplus formula.
c Figure 2 data is output all Figures supported.
       program airshowerD
c  
       integer*4, parameter :: n=100000
       integer*4 i
       real*8 h, L, z_C, dz
       real*8 c, c_a
       real*8 t(n), t_total(n), t_shower(n), t_light(n)
       real*8 z(n), vim(n), vim_r(n), vim_t(n)
       real*8 th(n), dth(n), thdot(n)
       real*8 dt_total(n), d(n), b(n), bground
c
       open (1, file="airshowerD.csv", status="unknown")
c                             Distances in km.
       h=25.0d0
       L=0.100
       c=2.99792458d5
       c_a=c/1.00029d0
       v=c
c                             Set up for grand loop: First iter.
       dz=h/dble(n-1)
       z(1)=0.0d0
       th(1)=0.0d0
       d(1)=sqrt(L**2 + z(1)**2)
       t_shower(1)=h/v
       t_light(1)=L/c_a
       t_total(1)=t_shower(1) + t_light(1)
c                             Second iter.
       z(2)=dz
       th(2)=datan(z(2)/L)
       d(2)=sqrt(L**2 + z(2)**2)
       t_shower(2)=dble(h - z(2))/v
       t_light(2)=d(2)/c_a
       t_total(2)=t_shower(2) + t_light(2)
       dt_total(2)=t_total(2) - t_total(1)
       vim(2)=abs(dz/dt_total(2))
       vim_t(2)=vim(2)*cos(th(2))
       thdot(2)=vim_t(2)/d(2)
       b(2)=thdot(2)/d(2)**2
       bground=b(2)
       t_min=1.0d0
       write (*,*) ' bground = ', bground
c                             Grand loop: Compute from bottom up.
       do 100 i=3, n
c                             Compute base parameters.
       z(i)=dble(i-1)*dz
       d(i)=sqrt(L**2 + z(i)**2)
       th(i)=datan(z(i)/L)
       t_shower(i)=dble(h - z(i))/v
       t_light(i)=d(i)/c_a
       t_total(i)=t_shower(i) + t_light(i)
       if (t_total(i).lt.t_min) t_min=t_total(i)
c                             Compute differentials.
       dth(i)=th(i) - th(i-1)
       dt_total(i)=t_total(i) - t_total(i-1)
       vim(i)=abs(dz/dt_total(i))
       vim_r(i)=vim(i)*sin(th(i))
       vim_t(i)=vim(i)*cos(th(i))
       thdot(i)=vim_t(i)/d(i)
c                             Compute relative brightness.
       b(i)=thdot(i)/d(i)**2
       b(i)=b(i)/bground
c                             Tell the world
c      write (*,*) z(i), abs(t_total(i) - t_min)*1.0d9
       write (1,*) z(i), b(i)
 100   continue
c                             Tell the world (part 2).
c      write (*,*) ' t_min = ', t_min*1.0d9
       write (*,*) ' z(i) km , b(i) / b_ground '
c
       do 200 i=1, n
c      if (i.eq.1) write (*,*) z(i), b(i)*bground
       if ((i/1000)*1000.eq.i) write (*,*) z(i), b(i)
 200   continue
c                             Stop the madness.
       stop
       end
