c                                                                      c
c                                                                      c
c             this program calculates steady 3d flow through           c
c             multiple turbomachine blade rows                         c
c             by solving the unsteady continuity,momentum,             c
c             and energy equations until a steady state is reached.    c
c                                                                      c
c           this version allows for multiple blade rows. there         c
c           is no intrinsic limit on the number of rows which          c
c           is determined by the parameter  "nrs" .                    c
c           increase this if required.                                 c
c                                                                      c
c                      development history                             c
c                                                                      c
c   the original develoment started about 1973 at cegb marchwood       c
c   this was for inviscid flow through a single blade row with         c
c   extremely coarse grids, typically 10x25x10 points.                 c
c                                                                      c
c   the early versions used the "opposed difference" numerical scheme  c
c   with cell centre storage on overlapping cells.                     c
c                                                                      c
c   it was extended to deal with two or more blade rows using a simple c 
c   plane model in 1979, oon after the author moved to cambridge.      c
c                                                                      c
c   around 1980 the scheme was changed to use cell corner storage of   c
c   the variables with non-overlapping cells and multigrid was         c
c   introduced. both pf these gave major improvements in performance.  c
c                                                                      c
c   the first nmoves to include viscous effects (1985) was by a separate
c   boundary layer calculation with transpiration through the blade    c
c   surfaces to allow for the boundary layer blockage. this was        c
c   surprisingly useful for compressors.                               c
c                                                                      c
c   around 1986 the first approach to including viscous forces         c
c   in the calculation used a skin friction coefficient and an         c
c   empircal  distribution of the viscous viscous stresses through the c
c   flow. this was soon extended to obtain the skin friction from wall c
c   functions and a simple mixing length model for the turbulent       c
c   viscosity. the viscous forces are included via a body force which  c
c   only needs to be upadted about every 5 time steps, giving          c
c   significant savings in run time.                                   c
c                                                                      c
c   numerous minor improvements were included throughout the 1990's    c
c   especially to the mixing plane treatment, and the program was      c
c   widely used for multistage compressors and turbines.               c
c   turbines. cooling flow and bleed flows were added together with    c
c   the pinched tip model for tip leakage and a model for shroud       c
c   leakage flows.                                                     c
c                                                                      c
c   a major development around 1998 was the change from the "opposed   c
c   difference" scheme to the "scree" scheme. this was simpler and     c
c   more accurate but the mixing plane treatment needed changing.      c
c                                                                      c
c   a major "tidying up" was performed in 2005 and some of the older   c
c   loss routines were removed or combined into a single subroutine    c
c   "loss" . an option to use real gas properties was added in 2006    c
c   and an allowance for surface roughness was added in 2008.          c
c                                                                      c
c   an option to use the spalart-allmaras turbulence model was added   c
c  in 2010 and at the same time an updated mixing length model "newlos"c
c  was added. both of these are full navier-stokes models whilst       c
c  the original routine "loss" is a thin shear layer approximation.    c
c                                                                      c
c   in 2014 the option to perform q3d calculations on a blade-to blade c
c   stream surface was added together with a new solution algorith,    c
c   the "sss" scheme, which permits larger cfl numbers.                c
c                                                                      c
c   the mixing plane model was further improved in 2015 to allow betterc
c   interation with shock waves and to pernmit reverse flows across it.c
c                                                                      c
c   also in 2015 a major tidying up was performed to obtain the currentc
c   version. this included considerable rearranging of the input data  c
c   so that data from previous versions is no longer compatiable.      c
c                                                                      c
c======================================================================c
c                                                                      c
c   to change the dimensions of the arrays perform a global change     c
c   of the  parameter statements in "commall-open-19.2"                c
c                                                                      c
c   i.e of 'id, jd, kd, maxki, nrs, ig1, jg1, kg1, ig2, jg2, kg2, jg3. c
c                                                                      c
c************************************************************************
c
c
c        this is the main program which is just used to call the
c                            subroutines.
c   
c
      character*1 ansin
c
      open(unit=1,file='/dev/tty')
      open(unit=4,file ='stage.log')
      open(unit=7,file ='flow_out',   form = 'unformatted')
      open(unit=11,file='global.plt', form = 'unformatted')
      open(unit=3,file ='results.out')
      open(unit=12,file='stopit')
c
      ifstop = 0
      write(12,*) ifstop
      close(12)
c
c******************************************************************************
c     decide which input style to use
c
      open(unit=13,file='intype')
            read(13,*,end=10,err=10)   ansin
      go to 20
   10 write(6,*) 'stopping because file  "intype" does not exist '
      stop 
c
   20 continue
c
      close(13)
c
      if(ansin.eq.'n'.or.ansin.eq.'n') then
             write(6,*) ' new_readin data format specified.'
             call new_readin
      else
c
      if(ansin.eq.'o'.or.ansin.eq.'o')then
             write(6,*)' old_readin data format specified.'
             call old_readin
      else
c
      write(6,*) ' stopping because file "intype" does not contain  "o" 
     & or  "n" '
      stop
c
      end if
c
      end if
c******************************************************************************
c 
      call setup(ansin)
c
c   loop calls  many other subroutines, especially tstep .
      call loop
c
c*********************************************************************************
c
      stop
      end
c
c**************************************************************************************
c**************************************************************************************
c**************************************************************************************
c
      subroutine new_readin
c
c       this subroutine reads in the data in
c       ====================

c
      include 'commall-open-19.2'
c
      common/bkodds/
     &           xint(jd),yint(jd),rint(jd),sdist(jd),
     &           icusp(nrs),lcusp(nrs),lcuspup(nrs),
     &           fracnew(jd),betanew(jd),slope(jd),
     &           thickup(jd),thicklow(jd),
     &           xint1(kd),xint2(kd),xint3(kd),xint4(kd)
c
      dimension xhub(jd),rhub(jd),xtip(jd),rtip(jd),xqo(maxki),
     &          rqo(maxki),xnewhub(jd),rnewhub(jd),xnewtip(jd),
     &          rnewtip(jd)
c
c       start of input section
c       throughout the input section the variables are as follows
c       =====================
c
c       xsurf(j,k)    is the store for the input blade axial coordinates
c       rt_upp(j,k)   is the store for the input blade suction surface co-ord's
c       rt_thick(j,k) is the store for the input the blade tangential thickness
c       rsurf(j,k)    is the store for the input blade radial coordinates
c
c       =====================
 1200 format(a72)
 1700 format(8f12.6)
 1730 format(8f10.3)
 1800 format(8f10.1)
 1610 format(40i2)
 1000 format(10i5)
c
      pi     = 3.14159265
      degrad = pi/180.
      raddeg = 180./pi
c
c**************************************************************************************
c**************************************************************************************
c
c    			 start to input data.
c     first for quantities which are not dependent on the blade row
c
c**************************************************************************************
c**************************************************************************************
c
c     input a title for the run. any characters in rows 1 to 72.
c
      write(6,*)
      write(6,*)
      read(5,1200)  title
      write(6,1200) title
      write(6,*)
      write(6,*)
c*******************************************************************************
c     input the gas properties, defaults to cp =1005, gamma = 1.4
c
      ifgas = 0
      cp = 1005.0
      ga = 1.4
      read(5,*) dummy_input
      read(5,*,err=7000) cp,ga
 7000 continue
      write(6,*) ' gas properties:    cp = ', cp, ' gamma = ',ga
      write(6,*)
c
c    use real gas properties if cp is input as negative.
c    typical values for combustion products  are:cp1 = 1272.5, cp2 = 0.2125,
c    cp3 = 0.000015625, rgas = 287.15 at  tref = 1400 k.
c
      if(cp.lt.0.0) then
      cp1  = 1272.5
      cp2  = 0.2125
      cp3  = 0.000015625
      tref = 1400.0
      rgas = 287.15
      read(5,*,err=7001) cp1, cp2, cp3, tref, rgas
 7001 continue
c
      write(6,*) ' ideal gas properties read in '
      write(6,*) ' cp1, cp2, cp3, tref, rgas =',cp1,cp2,cp3,tref,rgas
      write(6,*)
c
      cpgas  = cp1
      gagas  = cp1/(cp1 - rgas)
      cp     = cp1
      ga     = gagas
      cv     = cp/ga
      ifgas  = 1
      call set_coeffs
      end if
c
c******************************************************************************
c    input the time stepping option, 3 , 4 , -4,  5 or 6 . default = 3 .
c
      itimst = 3
      read(5,*) dummy_input
      read(5,*,err= 7002) itimst
 7002 continue
      write(6,*) ' time step type , itimst = ', itimst
      write(6,*)
c
c    set the coefficients for the sss scheme
c
      if(itimst.eq.3.or.itimst.eq.5.or.itimst.eq.6) then
              f1          =  2.0000
              f2          = -1.000
              f3          =  0.00
              f2eff       = -1.0
              nrsmth      =  0
              rsmth       =  0.40
      end if
c
      if(itimst.eq.-3.or.itimst.eq.-5.or.itimst.eq.-6) then
              f1          =  2.0000
              f2          = -1.65
              f3          = -0.65
              f2eff       = -1.0
              nrsmth      =  1
              rsmth       =  0.40
              itimst      =  abs(itimst)
      end if
c
      if(itimst.eq.4.or.itimst.eq.-4) then
              read(5,*) f1, f2eff, f3 , rsmth, nrsmth
              write(6,*) ' f1, f2eff, f3 , rsmth, nrsmth ',
     &                     f1, f2eff, f3 , rsmth, nrsmth
              write(6,*)
              if(f2eff.gt.0.0) then 
              write(6,*)
              write(6,*) ' error,  f2eff  must be negative.'
              write(6,*) ' the sign of the input value will be changed.'
              write(6,*)
              f2eff = -f2eff
              end if
              if(f3.gt.0.0) then 
              write(6,*)
              write(6,*) ' error,  f3  must be negative.'
              write(6,*) ' the sign of the input value will be changed.'
              write(6,*)
              f3  = -f3
              end if
              f2     = f2eff*(1.0 - f3)
              itimst = 3
      end if
c
              write(6,*) ' f1, f2eff, f3 , rsmth, nrsmth ',
     &                     f1, f2eff, f3 , rsmth, nrsmth
              write(6,*)
c
c**********************************************************************************
c
c    input the artificial speed of sound if itimst = 5.
c    this should be about half the maximum relative velocity in the flow.
c
      if(itimst.ge.5) then
      vsound    = 150.
      rf_ptru   = 0.01
      rf_vsound = 0.002
      densty    = 1.20
      vs_vmax   = 2.0
      if(itimst.eq.5)  read(5,*,end= 2350)
     &  vsound, rf_ptru, rf_vsound, vs_vmax
      if(itimst.eq.6)  read(5,*,end= 2350)
     &  vsound, rf_ptru, rf_vsound, vs_vmax, densty
 2350 continue
         rf_ptru1    = 1.0 - rf_ptru
         rf_vsound1  = 1.0 - rf_vsound
         write(6,*)
         write(6,*) ' calculation using artificial compressibility '
         write(6,*) ' artificial speed of sound = ', vsound
         write(6,*) ' density relaxation factor = ', rf_ptru
         write(6,*) ' sound speed relaxation factor = ', rf_vsound
         write(6,*) ' ratio of sound speed to maximum speed = ',vs_vmax
         if(itimst.eq.6) write(6,*) 
     &   ' incompressible flow with density = ',densty
         write(6,*)
      end if
c
c******************************************************************************
c   input the cfl number, default value = 0.4, but increase to 0.7 if using the sss scheme.
c   also the damping factor, mach number limiter and fraction of pressure downwinding.
c
      cfl     = 0.4
      dampin  = 10.0
      machlim = 2.0
      f_pdown = 0.0
      read(5,*) dummy_input
      read(5,*,end=7003) cfl, dampin, machlim, f_pdown
 7003 continue
      write(6,*) ' cfl number = ', cfl,' damping factor = ',dampin,
     & ' machlim = ',machlim, ' f_pdown = ', f_pdown
      write(6,*)
c
c******************************************************************************
c    read in the restart and output file options.
c    set  "if_restart"  = 1 to start from a restart file.  
c    the combined restart and plotting file,  "flow_out"   is always written.
c
      if_restart  = 0
      read(5,*)            dummy_input
      read(5,*,err = 7017) if_restart
 7017 continue
      write(6,*) 'restart file options, if_restart= ',if_restart
      write(6,*)
c
c**************************************************************************
c   input the maximum number of time steps and convergence limit.
c
      nsteps_max = 5000
      conlim     = 0.005
      read(5,*) dummy_input
      read(5,*,err = 7004) nsteps_max, conlim
 7004 continue
      write(6,*) ' maximum time steps= ',nsteps_max,
     &           ' convergence limit = ',conlim
      write(6,*)
c
c****************************************************************************
c      input the smoothing factors, fraction of fourth order smoothing and steps over
c      which they are gradually decreased
c      default values are:   0.005,  0.005 , 0.8 and nsteps/4  .
c
      sfxin   = 0.005
      sftin   = 0.005
      fac_4th = 0.8
      nchange = nsteps_max/4
      read(5,*) dummy_input
      read(5,*,err=7005) sfxin,sftin,fac_4th,nchange
 7005 continue
c    8/4/2017.  set nchange = 100 if starting from a restart file.
      if(if_restart.ne.0) nchange = 100
      write(6,*) ' sfxin,sftin,fac_4th,nchange ',
     &             sfxin,sftin,fac_4th,nchange
      write(6,*)
c
c******************************************************************************
c******************************************************************************
c   input the number of blade rows to be calculated
c
      read(5,*)    dummy_input
      read(5,*)    nrows
      write(6,*) ' number of blade rows to be calculated = ',nrows
      write(6,*)
c
c**************************************************************************
c   input the number of grid points in the pitchwise (i) and spanwise (k) directions.
c
      read(5,*)    dummy_input
      read(5,*)    im,km
      write(6,*) ' number of pitchwise grid points = ',im
      write(6,*) ' number of spanwise grid points  = ',km
      imm1 = im-1
      imm2 = im-2
      if(im.eq.2) imm2=1
      kmm1 = km-1
      kmm2 = km-2
      write(6,*)
c
c
c   input the relative spacing of the grid points in the pitchwise direction. 
      read(5,*)     dummy_input
      write(6,*)    dummy_input
      read(5,*)    (fp(i),i=1,imm1)
      write(6,*) ' relative spacing of the grid points in the pitchwise
     & direction .'
      write(6,*) ' this is the same for all blade rows.'
      write(6,1700) (fp(i),i=1,imm1)
      write(6,*)
c
c   input the relative spacing of the grid points in the spanwise direction.'
      read(5,*)      dummy_input
      write(6,*)     dummy_input
      read(5,*)     (fr(k),k=1,kmm1)
      write(6,*) ' relative spacing of the grid points in the spanwise
     & direction .'
      write(6,*) ' this is the same for all blade rows.'
      write(6,1700) (fr(k),k=1,kmm1)
      write(6,*)
c**************************************************************************      
c    input the multigrid block sizes. defaults = 3  and 9
c
      ir   = 3
      jr   = 3
      kr   = 3
      irbb = 9
      jrbb = 9
      krbb = 9
      read(5,*) dummy_input
      read (5,*,err=7007) ir,jr,kr,irbb,jrbb,krbb
 7007 continue
      if(km.eq.2) then
           kr   = 1
           krbb = 1
      end if
      if(im.eq.2) then
           ir   = 1
           irbb = 1
      end if
c
      write(6,*) ' multigrid block sizes = ', ir,jr,kr,irbb,jrbb,krbb
      write(6,*)
c
c************************************************************************
c   input the multigrid time step factors, fblk1,  fblk2, fblk3
c
      fblk1 = 0.4
      fblk2 = 0.2
      fblk3 = 0.1
      read(5,*) dummy_input
      read(5,*,err=8007) fblk1,fblk2,fblk3
      write(6,*) ' multigrid time step factors= ',fblk1,fblk2,fblk3
 8007 continue
      write(6,*)
c
c******************************************************************************
c    read in the mixing plane parameters
c
      ifmix = 1
      read(5,*)    dummy_input
      read(5,*,err=7009) ifmix
 7009 continue
      write(6,*) ' mixing plane marker, ifmix = ', ifmix
      write(6,*)
c
      if(ifmix.ne.0) then
      rfmix    = 0.025
      fsmthb   = 1.0
      fextrap  = 0.80
      fangle   = 0.80
      read(5,*) dummy_input
      read(5,*,end=7010,err=7010)  rfmix,fextrap,fsmthb,
     &                             fangle
 7010 continue
      write(6,*)' the mixing plane parameters are '
      write(6,*)' rfmix, fextrap, fsmthb ,fangle',
     &            rfmix, fextrap, fsmthb, fangle
      write(6,*)
      end if
c
c******************************************************************************
c   read in markers for bleed flows, cooling flows and surface roughness
c
      ifcool   = 0
      ifbleed  = 0
      if_rough = 0
      read(5,*) dummy_input
      read(5,*,err = 7011) ifcool,ifbleed,if_rough
 7011 continue
      write(6,*) ' ifcool,ifbleed,if_rough ', ifcool,ifbleed,if_rough
      write(6,*)
c
c******************************************************************************
c   read in the number of blade sections on which the geometry is to be generated.
c 
      read(5,*)   dummy_input 
      write(6,*)  dummy_input 
      read(5,*)   nsecs_in
      write(6,*)
     &'number of sections used for blade geometry generation= ',nsecs_in
      write(6,*)
c
c******************************************************************************
c      read in the inlet boundary condition options
c
      in_press   = 0
      in_vtan    = 0
      in_vr      = 1
      in_flow     = 0
      if_repeat  = 0
      rfin       = 0.1
      read(5,*) dummy_input
      read(5,*,err=7012) in_press,in_vtan,in_vr,in_flow,if_repeat,rfin
 7012 continue
      write(6,*) ' in_press,in_vtan,in_vr,in_flow,if_repeat ',
     &             in_press,in_vtan,in_vr,in_flow,if_repeat
      write(6,*) ' rfin = ', rfin
      write(6,*)
c
c******************************************************************************
c
c      read in the exit boundary condition options
c      new  if_wave = 0  or 1  option added november 2019.
c      defaults to the original model, if_wave = 0 , if fracwave = 0.
      ipout          = 1
      sfexit         = 0.0 
      nsfexit        = 0
      fp_xtrap       = 0.9
      fracwave       = 0.0 
      if_wave        = 0
      read(5,*)          dummy_input
      read(5,*,end=8012,err=8012) ipout,sfexit,nsfexit,fp_xtrap,fracwave
      read(5,*)   dummy_input
 8012 continue      
      write(6,*)' sfexit= ',sfexit, ' nsfexit= ', nsfexit,' fp_xtrap =',
     &            fp_xtrap,' fracwave = ',fracwave
      if(fracwave.gt.0.001) if_wave = 1
      write(6,*) ' if_wave = ', if_wave
c
c****************************************************************************** 
c
      plate_loss    = 0.0
      throttle_exit = 0.0

      read(5,*,err=8013)   plate_loss, throttle_exit
 8013 continue
      write(6,*)' loss coefficient of a perforated plate at exit = ',
     &            plate_loss
      write(6,*)' marker for setting a throttle at the exit boundary= ',
     &            throttle_exit 
c
      if(throttle_exit.gt.0.001) then
            read(5,*)  throttle_pres,throttle_mas,rfthrotl
            write(6,*)'a throttle exit boundary condition has been set.'
            write(6,*)'throttle_pres, throttle_mas, rfthrotl=',
     &                 throttle_pres,throttle_mas,rfthrotl 
           rfthrotl1 = 1.0 - rfthrotl
      end if
      write(6,*)
c
c******************************************************************************
c     read in in the specified inlet flow and relaxation factor if in_flow
c     is not = zero
c
      if(in_flow.ne.0) then
           read(5,*)    dummy_input
           read(5,*)    flowin, rflow
           write(6,*) ' inlet flow forcing, flowin, rflow ',flowin,rflow
           write(6,*)
      end if
c
c******************************************************************************
c     read in the factors used for repeating flow conditions
c
      if(if_repeat.ne.0) then
           ninmod = 10
           rfinbc = 0.025
           read(5,*) dummy_input
           read(5,*,err=7013)  ninmod, rfinbc
 7013      continue
           write(6,*) ' repeating stage specified, ninmod = ', ninmod,
     &                ' rfinbc = ', rfinbc
           write(6,*)
      end if
c
c******************************************************************************
c    read in the choice of viscous model
c
      ilos   = 100
      nlos   = 5
      ibound = 0
      read(5,*) dummy_input
      read(5,*,err=7014) ilos, nlos, ibound
 7014 continue
      if(km.eq.2) ibound = 2
      write(6,*) ' viscous model set by ilos= ',ilos,' nlos= ',nlos,
     &           ' ibound= ', ibound
      write(6,*)
c
c******************************************************************************
c    read in the parameters for the viscous model.

      reyno       = 500000.0
      rf_vis      = 0.5
      ftrans      = 0.0001
      prandtl     = 1.0
      ypluswall   = 0.0
      turbvis_lim = 3000.0
c
      if(ilos.ne.0) then
      read(5,*)          dummy_input
      read(5,*,err=7015) reyno,rf_vis,ftrans,turbvis_lim,
     &                   prandtl,ypluswall
 7015 continue
c
      write(6,*)
      write(6,*) ' reynolds number                       = ',reyno
      write(6,*) ' viscous term relaxation factor        = ',rf_vis
      write(6,*) ' transition factor, ftrans             = ',ftrans
      write(6,*) ' limit on turbulent/laminar viscosity  = ',turbvis_lim
      write(6,*) ' prandtl number                        = ',prandtl
      write(6,*) ' ypluswall - if used (not often used)  = ',ypluswall
      write(6,*)
c
      end if
c
c   if using the sa turbulence model
c
      if(ilos.ge.200)then
           fac_stmix = 0.0
           fac_st0   = 1.0
           fac_st1   = 1.0
           fac_st2   = 1.0
           fac_st3   = 1.0
           fac_sfvis = 2.0
           fac_vort  = 0.0
           fac_pgrad = 0.0
           read(5,*) dummy_input
           read(5,*,err= 7016) fac_stmix, fac_st0, fac_st1,
     &                         fac_st2  , fac_st3, fac_sfvis,
     &                         fac_vort, fac_pgrad
 7016 continue
c
           write(6,*)
           write(6,*) ' spalart-allmaras turbulence model is being used'
           write(6,*) ' the s_a source term multipliers are ',
     &                  fac_stmix, fac_st0, fac_st1, fac_st2, fac_st3,
     &                  fac_vort, fac_pgrad 
           write(6,*) ' the turbulent viscosity smoothing factor is ', 
     &                  fac_sfvis
           write(6,*)
      end if
c
c**********************************************************************************
c    input the range of yplus values over which the turbulent viscosity will be reduced.
c
      if(ilos.ne.0) then
           yplam    = 5.0
           ypturb   = 25.0
           read(5,*) dummy_input
           read(5,*, err= 590) yplam,ypturb
  590 continue
           write(6,*)' the turbulent viscosity is reduced over the range 
     &  yplus = ',yplam,' to ',ypturb 
      end if
c
c******************************************************************************
c   input the forcing factor and smoothing factor if doing a throughflow  calculation
c
      if(im.eq.2) then
      q3dforce = 1.0
      sfpbld   = 0.1
      nsfpbld  = 2
      read(5,*) dummy_input
      read(5,*,err=591) q3dforce, sfpbld, nsfpbld
  591 continue
      sfpbld1 = 1.0 - sfpbld
      write(6,*)
      write(6,*) 'throughflow calculation requested, q3dforce= ',
     & q3dforce, 'smoothing factor=',sfpbld,' no of smoothings=',nsfpbld
      write(6,*)
      end if       
c
c******************************************************************************
c    read in the blade row and grid alignment options
c
      ishift      = 2
      nextrap_le  = 10
      nextrap_te  = 10
      read(5,*)          dummy_input
      read(5,*,err=7008) ishift,nextrap_le,nextrap_te
 7008 continue
      write(6,*) ' ishift ,nextrap_le, nextrap_te = ',
     &             ishift, nextrap_le, nextrap_te
      write(6,*)
c
c**************************************************************************************
c
c     read in the stage numbers and sort out the start and end of each stage.
c
c     first set defaults.
      do  n = 1,nrows
      nstage(n)    = 1 + (n-1)/2
      end do
c
c     now read in the actual stage number for each blade row.
c
      read(5,*) dummy_input
      read(5,*,err=588)(nstage(n),n=1,nrows)
  588 continue
c
c**************************************************************************
c**************************************************************************
c   input 5 time steps where output is requested.
c   set to a number above nsteps_max if none required.
c
      do n=1,5
           nout(n) = nsteps_max+10 
      end do
      read(5,*) dummy_input
      read(5,*,err=7006) (nout(n),n=1,5)
 7006 continue
      write(6,*) ' output requested at time steps = ', (nout(n),n=1,5)
      write(6,*)
c******************************************************************************
c     input a list of variables to be sent to the output file ' results.out'.
      do i=1,20
           iout(i) = 0
      end do
      read(5,*)      dummy_input
      read(5,*,err=7027) (iout(i),i=1,13)
 7027 continue
      write(6,*)'the variables to be output to the file results.out are'
      write(6,1610) (iout(i),i=1,13)
c***************************************************************************************
c     choose which k values (stream surfaces) are to be output to the file 'results.out'.
      do k=1,km
           kout(k) = 0
      end do
      read(5,*)      dummy_input
      read(5,*,err=2028) (kout(k),k=1,km)
 2028 continue
      write(6,*)  'the stream surfaces on which results are sent to the
     &output file are:'
      write(6,1610) (kout(k),k=1,km)
c
c******************************************************************************
c**************************************************************************************
c    end of data input which is not dependent on the blade row.
c**************************************************************************************
c**************************************************************************************
c
c     check that the dimensions are not too large.
c
      if(maxki.lt.kd)  write(6,*)
     &            ' stopping because maxki is less than kd',
     &            ' maxki = ',maxki, ' kd = ',kd
      if(maxki.lt.id)  write(6,*)
     &            ' stopping because maxki is less than id',
     &            ' maxki = ',maxki, ' id = ',id
      if(im.gt.id)  write(6,*) ' stopping because im too large.',
     &            ' im= ',im,  ' dimension limit = ',id
      if(km.gt.kd)  write(6,*) ' stopping because km too large.',
     &            ' km= ',km,  ' dimension limit = ',kd
c
      if(maxki.lt.kd.or.maxki.lt.id.or.im.gt.id.or.km.gt.kd) stop
c
      if(nrows.gt.nrs)  write(6,*) ' stopping because nrows too large.',
     &            ' nrows= ',nrows,' dimension limit = ',nrs
      if(nrows.gt.nrs) stop
c
c**************************************************************************************
c**************************************************************************************
c**************************************************************************************
c**************************************************************************************
c
c            now input the control parameters for each blade row
c            and set the values of some blade row variables
c
      j1        = 1
      ifshroud  = 0
c
      write(6,*)'******************************************************'
      write(6,*)'******************************************************'
c
      write(6,*)
      write(6,*) ' starting input for each blade row '
      write(6,*)
c
      do 1550 nr = 1,nrows
c
c    read in two blank lines to help space out the data.
      read(5,*) dummy_input
      read(5,*) dummy_input
      write(6,*)'******************************************************'
      write(6,*)'******************************************************'
c
c******************************************************************************
c     read in the row title. this is never used but helps to lay out the data
      read(5,1200) rowtyp
c
      write(6,*) ' inputting data for blade row number ',nr
      write(6,*) ' blade row title = ', rowtyp
      write(6,*)'******************************************************'
      write(6,*)'******************************************************'
c
c******************************************************************************
c     read in the number of blades in this row  
      read(5,*)    dummy_input    
      read(5,*)    nblades_in_row
      write(6,*)
      write(6,*) ' number of blades in row no.',nr,' =',nblades_in_row
      write(6,*)
c
c******************************************************************************
c    jmrow = number of j grid points on this row.
c    jlerow and jterow are measured from the start of this row. 
c    they are not the final overall values.
      read(5,*)    dummy_input
      read(5,*)    jmrow,jlerow,jterow
      write(6,*) ' jm = ',jmrow,'jle =',jlerow,'jte = ',jterow
      write(6,*)
c     set the overall j values of the leading and trailing edge points on this row.   
      jlee = j1+jlerow-1
      jtee = j1+jterow-1
c
c******************************************************************************
c
c      ktips is the k value of the point where the tip gap starts
c      ktipe is the k value of the point where the tip gap ends.
c      set ktips(nr) = 0  for no tip gap on this row.    
c
c      set ktips(nr) negative to use the shroud leakage model on this
c      row. extra input data is then needed at the end of the data file.
c
      read(5,*)    dummy_input
      read(5,*)    ktips(nr),ktipe(nr)
      write(6,*) ' k tip start= ',ktips(nr),' k tip end= ',ktipe(nr)
      write(6,*)
c     for q3d
      if(km.eq.2) then
           ktips(nr) = 0
           ktipe(nr) = 0
      end if
c     end q3d
c     set marker  "ifshroud" for a tip shrouud if ktips is negative on any blade row.
c
      if(ktips(nr).le.0) ktipe(nr) = 0
c     set   "ifshroud" = 1   if any blade row has a shroud.
c     if no blade row has a shroud   "ifshroud"   remains at 0 .
      if(ktips(nr).lt.0) ifshroud  = 1
c
c*******************************************************************************
c      input the tip clearance as a fraction of the blade span at the leading edge
c      and trailing edge.
       if(ktips(nr).gt.0) then
            read(5,*)   dummy_input
            read(5,*)   fractip1(nr),fractip2(nr)
            write(6,*)' fractip1, fractip2 = ',fractip1(nr),fractip2(nr)
            write(6,*)
      end if
c
c*******************************************************************************
c      fthick(nr,k) is the multiplying factor on the blade thickness so
c      that it can be reduced at and beyond the tip.
c
      do 2507 k=1,km
 2507 fthick(nr,k)=1.0
c
c      fthick  is assumed to be 1.0 unless input here.
c      set fthick = 0.0 for points in the tip gap, including the point at the blade tip.
      if(ktips(nr).gt.0.and.ktipe(nr).gt.0) then
           write(6,*)   ' blade tip thinning factors '
           read(5,*)      dummy_input
           read(5,*)     (fthick(nr,k),k=1,km)
           write(6,1700) (fthick(nr,k),k=1,km)
           write(6,*) 
      end if    
c
c******************************************************************************
c   boundary layer transition is specified separately on each surface.
c   but is also determined by ftrans which has already been input . 
c   jtran_i1, etc,  are measred from the start of this row.
c
      read(5,*)    dummy_input
      read(5,*)    jtran_i1(nr),jtran_im(nr),jtran_k1(nr),jtran_km(nr)
      write(6,*) ' boundary layer transition is specified as follows '
      write(6,*) ' transition on i= 1 surface fixed at j= ',jtran_i1(nr)
      write(6,*) ' transition on i=im surface fixed at j= ',jtran_im(nr)
      write(6,*) ' transition on k= 1 surface fixed at j= ',jtran_k1(nr)
      write(6,*) ' transition on k=km surface fixed at j= ',jtran_km(nr)
      write(6,*)
c
c******************************************************************************
c     set new_grid = 1  to generate a new grid in the "j" (streamwise) direction
c     data for this is read in later. set new_grid = 0 to use the grid input here.
c
      read(5,*)    dummy_input
      read(5,*)    new_grid
      write(6,*) ' newgrid = ', new_grid
      if(new_grid.ne.0) write(6,*)
     & ' a new grid will be generated using data to be input later. '
      write(6,*)
c
c
c******************************************************************************
c     read in the rpm and hub rotation of this row. 
c     the casing rotation is taken to be rpmrow between jrotts and jrotte and zero elsewhere.
      read(5,*)    dummy_input
      read(5,*)    rpmrow, rpmhub  
      write(6,*) ' blade row rpm = ',rpmrow,' hub rpm = ',rpmhub
      write(6,*)
c
c******************************************************************************
c     the hub is taken to be rotating at rpmhub between jroths and jrothe .
c     the casing is taken to be rotating at rpmrow between jrotts and jrotte.
c     outside these limits the hub and casing are not rotating.
c     jroths, jrothe, etc  are measured from the j value at start of this blade row
c
      read(5,*)    dummy_input
      read(5,*)    jroths,jrothe,jrotts,jrotte
c
c   set all rotation points to jmrow if they are outside the limits for this row.
      if(jroths.eq.0.or.jroths.gt.jmrow) jroths = jmrow
      if(jrothe.eq.0.or.jrothe.gt.jmrow) jrothe = jmrow
      if(jrotts.eq.0.or.jrotts.gt.jmrow) jrotts = jmrow
      if(jrotte.eq.0.or.jrotte.gt.jmrow) jrotte = jmrow
c
      write(6,*) ' hub rotation starts at j = ',jroths
      write(6,*) ' hub rotation ends   at j = ',jrothe
      write(6,*) ' tip rotation starts at j = ',jrotts
      write(6,*) ' tip rotation ends   at j = ',jrotte
      write(6,*)
c
c****************************************************************************** 
c     input an initial guess of the upstream, leading edge, trailing edge and 
c     downstream static pressures for this row. 
      read(5,*)    dummy_input
      read(5,*)    puprow,plerow,pterow,pdnrow
      write(6,*) ' initial guess of pressures = ',
     &             puprow,plerow,pterow,pdnrow 
      write(6,*)
c
c******************************************************************************
c     input the number of input sections for this blade row.
      read(5,*)   dummy_input
      nsecs_row  = nsecs_in
      insurf     = 0
      read(5,*,err=7019) nsecs_row,insurf
 7019 continue
      write(6,*)'number of input sections for this blade row=',nsecs_row
      write(6,*)'     new endwall generation option, insurf =', insurf
      if(nsecs_in.eq.1)  nsecs_row = 1
      if(nsecs_in.eq.1)  nsecs_in  = 2
      if(nsecs_row.eq.1) nsecs_row = 2
c
c******************************************************************************
c******************************************************************************
c     set the cusp generation parameters if if_cusp is not 0.
c     maintain the original grid and no cusp is generated if if_cusp = 0
c     the cusp is centred on the blade centre line if icusp = 0.
c     the cusp makes the i=1 surface continuous on the cusp if  icusp =  1.
c     the cusp makes the i=im surface continuous on the cusp if icusp = -1.
c     the cusp is of length  lcusp and starts lcuspup points before the
c     trailing edge.
c
      if_cusp(nr)    = 0
      if_angles(nr)  = 0
      read(5,*)  dummy_input
      read(5,*,err=7020) if_cusp(nr), if_angles(nr)
 7020 continue
      write(6,*)'           cusp generation option, if_cusp =',
     &           if_cusp(nr)
      write(6,*)'upstream and downstream angle input option =',
     &           if_angles(nr)
      write(6,*)
c
      if(if_cusp(nr).eq.1) then
	   read(5,*)    icusp(nr),lcusp(nr),lcuspup(nr)
           write(6,*) ' icusp, lcusp , lcuspup = ',
     &                  icusp(nr),lcusp(nr),lcuspup(nr)
           write(6,*)
      end if
c
c    if  ifcusp = 2 a body force is used to force separation at the trailing edge.
c    the force starts nup_i1 grid points upstream of the te  on the i=1 blade surface, 
c    and at nup_im points upstream on the i=im blade surface. it extends n_wake points
c    downstream of the te. the thickness of the affected layer is determined by sep_thik,
c    typical value  0.01, and the strength of the body force by sep_drag, typical value 0.99.
c
      if(if_cusp(nr).eq.2) then
           icusp(nr)   = 0
           lcusp(nr)   = 0
           lcuspup(nr) = 0
           read(5,*)  nup_i1(nr),nup_im(nr),n_wake(nr),
     &                  sep_thik(nr),sep_drag(nr)
           write(6,*)' body force used to force trailing edge separtion'
           write(6,*)' cusp body force points, nup_i1, nup_im, n_wake, 
     & sep_thick, sep_drag'
             write(6,*) nup_i1(nr),nup_im(nr),n_wake(nr),
     &                  sep_thik(nr),sep_drag(nr)
             write(6,*)
      end if
c******************************************************************************
c******************************************************************************
c
c   set the j value relative to the overall first grid point
      j2 = j1 + jmrow - 1
c
c*******************************************************************************
c*******************************************************************************
c*******************************************************************************
c
c      now  read in main geometrical data for the current blade row 
c      on  "nsecs_row"  blade sections.
c
c*******************************************************************************
c*******************************************************************************
c*******************************************************************************
c
      write(6,*)'*******************************************************
     &********************'
      write(6,*)
     &      ' starting to input the blade coordinates for row number',nr
      write(6,*)'*******************************************************
     &********************'
      write(6,*)
c
      do 1600 k=1,nsecs_row
c
      write(6,*)'*******************************************************
     &********************'
      write(6,*)'*******************************************************
     &********************'
c
      if_design  = 0
      if_stagger = 0
      if_lean    = 0
      read(5,*)  dummy_input
      read(5,*)  dummy_input 
      read(5,*,  err=1601) if_design,if_stagger,if_lean
 1601 continue 
c
      write(6,*) ' if_design = ', if_design
      if(if_design.ne.0) write(6,*)
     & ' this blade section is to be designed'
c
      write(6,*) ' if_stagger = ', if_stagger
      if(if_stagger.ne.0) write(6,*)
     & ' this blade section is to be re-staggered'
c
      write(6,*) ' if_lean = ', if_lean
      if(if_lean.ne.0) write(6,*)
     & ' this blade section is to be leaned'
c
c******************************************************************************
c     do not read in any geometry if "if_design" is non zero. 
c     instead call subroutine re_design to generate a new section.
c     then jump to 1595 .
c
      if(if_design.ne.0) then
          write(6,*)'calling re_design to generate a new blade geometry'
          call re_design(nr,k,j1,j2,jlerow,jterow) 
          go to 1595
      end if
c
c******************************************************************************
c******************************************************************************
c
c   now input the blade geometry on the quasi stream surface .
c******************************************************************************
c******************************************************************************
c
c******************************************************************************
c***********xsurf(j,k) is axial coordinate of points on the streamwise surface.
c
      write(6,*)  ' row number ',nr,' section number ',k,
     &            ' x  coordinates '      
      read(5,*)     fac1, xshift
      write(6,*)  ' fac1, xshift = ', fac1,xshift
      write(6,*)
      read(5,*)     (xsurf(j,k),j=j1,j2)
      write(6,1700) (xsurf(j,k),j=j1,j2)
      write(6,*)
      xrange = xsurf(j2,k) - xsurf(j1,k)      
c
c***********rt_upp(j,k) is the r-theta coordinate of points on the streamwise
c           surface on the blade surface with largest value of theta.
c           i.e. the upper surface of the blade and the lower surface of the
c           blade to blade passage.
c
      write(6,*)  ' row number ',nr,' section number ', k,
     &            ' r_theta of upper blade surface'
      read(5,*)     fac2, tshift
      write(6,*)  ' fac2, tshift = ', fac2,tshift
      write(6,*) 
      read(5,*)     (rt_upp(j,k),j=j1,j2)
      write(6,1700) (rt_upp(j,k),j=j1,j2)
      write(6,*)
c
c***********rt_thick(j,k) is blade thickness delta r-theta measured
c          in the tangential direction.
c 
      write(6,*)  ' row number ',nr,' section number ',k,
     &            ' blade tangential thickness. '
      read(5,*)     fac3
      write(6,*)  ' fac3= ', fac3
      write(6,*)
      read(5,*)     (rt_thick(j,k),j=j1,j2)
      write(6,1700) (rt_thick(j,k),j=j1,j2)
      write(6,*)
c
c      rsurf(j,k) is the radii of points on the streamwise surfaces on which data is input.
c
      write(6,*)   ' row number ',nr,' section number ',k,
     &             ' stream surface radius. '
      read(5,*)      fac4, rshift
      write(6,*)   ' fac4, rshift  = ', fac4,rshift
      write(6,*)
      read(5,*)     (rsurf(j,k),j=j1,j2)
      write(6,1700) (rsurf(j,k),j=j1,j2)
      write(6,*)
c    end of all geometrical input for this blade section
c
c******************************************************************************
c******************************************************************************
c     check that the input meridional grid spacings are not too small
c     and adjust them if they are.
c
      sdist(j1) = 0.0
      do 1450 j= j1+1,j2
           xd = xsurf(j,k)  - xsurf(j-1,k)
           rd = rsurf(j,k)  - rsurf(j-1,k)
           sdist(j) = sdist(j-1) + sqrt(xd*xd+rd*rd)
 1450 continue
      srange = sdist(j2) - sdist(j1)
      slim   = 0.0001*srange
      do 1460 j = j1+2,j2
        if((sdist(j)- sdist(j-1)).lt.slim) then
           write(6,*) 
     &   ' warning, the input grid spacings are too small at j= ',
     &     j,'k= ',k
           write(6,*) ' adjusting the spacings and continuing. '
           xd = slim*(xsurf(j-2,k)-xsurf(j-1,k))/(sdist(j-2)-sdist(j-1))
           rd = slim*(rsurf(j-2,k)-rsurf(j-1,k))/(sdist(j-2)-sdist(j-1))
           xsurf(j,k) = xsurf(j-1,k) + xd
           rsurf(j,k) = rsurf(j-1,k) + rd
        end if
 1460 continue
c
c******************************************************************************
c******************************************************************************
c     scale and shift the input blade sections by xshift, tshift or rshift.
c
      do 1500 j=j1,j2
      xsurf(j,k)    = fac1*(xshift + xsurf(j,k))
      rt_thick(j,k) = fac3*rt_thick(j,k)
      rsurf(j,k)    = fac4*(rshift + rsurf(j,k))
      rt_upp(j,k)   = fac2*(tshift + rt_upp(j,k))
 1500 continue
c
c******************************************************************************
c****************************************************************************** 
c   re enter here if a new section was designed,  i.e  if  if_design was non zero.
c
 1595 continue
c
c******************************************************************************
c******************************************************************************
c    now start to rotate or restagger this blade section .
c
c******************************************************************************
c******************************************************************************
c     restagger the blade section if if_stagger is non zero.
c     positive  "rotate"  means clockwise rotation.
c
      if(if_stagger.ne.0 ) then
           write(6,*) ' restaggering section number ',k
           call restagger(k,j1,j2,jlerow,jterow,rotate,fracx_rot)
      end if
c
c******************************************************************************
c******************************************************************************
c     lean the blade by anglean if "if_lean" is greater than zero.
c     if "anglean" is positive the hub is held fixed and the other sections are 
c     moved in the positive theta direction .
c
      if(if_lean.ne.0) then
      write(6,*) '  leaning section number ', k
      call lean(k,j1,j2,anglean)
      end if
c
c************************** ****************************************************
c******************************************************************************
c   q3d
c     call set ssthick to set the stream surface thickness if doing a q3d
c     blade to blade calculation. then jump to 1605 to end the geometry input
c     for this blade row as only a single stream surface is used for q3d.
c
      if(km.eq.2) then
           call set_ssthick(j1,j2)
c    jump out of the do 1600 loop as only one stream surface is needed.
           go to 1605
      end if
c   end q3d
c
c******************************************************************************
c******************************************************************************
c     end of all input for this blade section
 1600 continue
c
c  re-enter here if  km = 2
 1605 continue
c
c******************************************************************************
c******************************************************************************
c     the blade geometry has now been input for all sections.
c     all the following input applies to the whole row.
c*******************************************************************************
c*******************************************************************************
c
c     set up new hub and casing coordinates and blade sections if insurf > 0 .
c     use the input coordinates and sections if insurf = 0 .
c
c*******************************************************************************
c*******************************************************************************
c
      if(insurf.eq.1.or.insurf.gt.2) then
           read(5,*) nhub
           read(5,*) (xhub(n),n=1,nhub)
           write(6,*) ' xhub = ', (xhub(n),n=1,nhub)
           read(5,*) (rhub(n),n=1,nhub)
           write(6,*) ' rhub = ', (rhub(n),n=1,nhub)
           do n = 1,nhub
               xhub(n) = (xshift + xhub(n) ) * fac1
               rhub(n) = (rshift + rhub(n) ) * fac4
           end do
      end if
c
      if(insurf.ge.2) then
           read(5,*) ntip
           read(5,*) (xtip(n),n=1,ntip)
           write(6,*) ' xtip =', (xtip(n),n=1,ntip)
           read(5,*) (rtip(n),n=1,ntip)
           write(6,*) ' rtip =', (rtip(n),n=1,ntip)
           do n = 1,ntip
               xtip(n) = (xshift + xtip(n) ) * fac1
               rtip(n) = (rshift + rtip(n) ) * fac4
           end do
      end if
c
c******************************************************************************
c******************************************************************************
      do 1650 j = j1,j2
c
      do 1660 k = 1,nsecs_row
           xqo(k)  = xsurf(j,k)
           rqo(k)  = rsurf(j,k)
 1660 continue
c
c      calling  insect  to find the intersection of the hub curve with
c      the blade qo curve.
c
      if(insurf.eq.1.or.insurf.gt.2) then
            call insect(jd,maxki,nhub,xhub,rhub,nsecs_row,xqo,rqo,
     &           xnewhub(j),rnewhub(j))
      else
            xnewhub(j) = xsurf(j,1)
            rnewhub(j) = rsurf(j,1)
      end if
c
c      calling  insect  to find the intersection of the  casing curve with
c      the blade qo curve.
c
      if(insurf.ge.2) then
            call insect(jd,maxki,ntip,xtip,rtip,nsecs_row,xqo,rqo,
     &           xnewtip(j),rnewtip(j))
      else
            xnewtip(j) = xsurf(j,nsecs_row)
            rnewtip(j) = rsurf(j,nsecs_row)
      end if
c
 1650 continue
c
c   new hub and casing cordinates have been found.
c
c******************************************************************************
c******************************************************************************
c  interpolate to obtain  "nsecs_in"  new blade sections equally spaced between 
c  the hub and casing.
c
      do 1670 j = j1,j2
c
      xspan  = xnewtip(j) - xnewhub(j)
      rspan  = rnewtip(j) - rnewhub(j)
c
      do 1671 k = 1,nsecs_row
           xint1(k)  = xsurf(j,k)
           xint2(k)  = rsurf(j,k)
           xint3(k)  = rt_upp(j,k)
           xint4(k)  = rt_thick(j,k)
 1671 continue
c
c  first if the spanwise direction is predominantly axial.
c
      if(abs(xspan).gt.abs(rspan) ) then
c
      do 1672 k = 1, nsecs_in
      f_span     = float(k-1)/float(nsecs_in -1)
      xarg       = xnewhub(j) + f_span*(xnewtip(j) - xnewhub(j) )
      xsurf(j,k) = xarg 
      call intp(nsecs_row,xint1,xint2,xarg,rsurf(j,k) )
      call intp(nsecs_row,xint1,xint3,xarg,rt_upp(j,k) )
      call intp(nsecs_row,xint1,xint4,xarg,rt_thick(j,k) )
 1672 continue
c
      else    
c
c     next if the spanwise direction is predominantly radial
c
      do 1673 k = 1, nsecs_in
      f_span     = float(k-1)/float(nsecs_in -1)
      rarg       = rnewhub(j) + f_span*(rnewtip(j) - rnewhub(j) )
      rsurf(j,k) = rarg 
      call intp(nsecs_row,xint2,xint1,rarg,xsurf(j,k) )
      call intp(nsecs_row,xint2,xint3,rarg,rt_upp(j,k) )
      call intp(nsecs_row,xint2,xint4,rarg,rt_thick(j,k) )
 1673 continue
c
      end if 
c 
 1670 continue
c
c   finished interpolating in input sections to form nsecs_in sections.
c
c************************************************** ****************************
c******************************************************************************
c     interpolate in the upstream and downstream grid angles if if_angles > 0
c
      if(if_angles(nr).gt.0) then
c
      write(6,*)  ' inputting the upstream and downstream grid angles'
c
      read(5,*)    dummy_input
      read(5,*)    n_angles
      read(5,*)   (fracn_span(nk),nk=1,n_angles)
      read(5,*)   (angl_up(nk),   nk=1,n_angles)
      read(5,*)   (angl_dwn1(nk), nk=1,n_angles)
      read(5,*)   (angl_dwn2(nk), nk=1,n_angles)
c
      qspan(1) = 0.0
      j_mid     = j1 + (jlerow + jterow)/2
      do 1680 k = 2,nsecs_in
           xdif = xsurf(j_mid,k) - xsurf(j_mid,k-1)
           rdif = rsurf(j_mid,k) - rsurf(j_mid,k-1)
           qspan(k) = qspan(k-1) + sqrt(xdif*xdif + rdif*rdif)
 1680 continue
      do 1681 k = 1,nsecs_in
      qspan(k) = qspan(k)/qspan(nsecs_in)
 1681 continue 
c
      do 1685 k = 1,nsecs_in
      arg = qspan(k)
      call intp(n_angles,fracn_span,angl_up,arg,betaup(nr,k))
      call intp(n_angles,fracn_span,angl_dwn1,arg,betadwn1(nr,k))
      call intp(n_angles,fracn_span,angl_dwn2,arg,betadwn2(nr,k))
c
      write(6,*)'interpolated grid angles on the new stream surfaces',
     &           betaup(nr,k),betadwn1(nr,k),betadwn2(nr,k)
 1685 continue
      write(6,*)
c
      end if
c
c******************************************************************************
c******************************************************************************
c  tflow
c******************************************************************************
c******************************************************************************
c  if a throughflow calculation read in the required exit flow angle or deviation 
c  angle for this blade row. angl_typ = 'a' for exit angle, = 'd' for deviation angle.
      if(im.eq.2) then
      read(5,*)    dummy_input
      write(6,*) ' inputting the throughflow data '
      read(5,*)    angl_typ(nr), nangles(nr)
      read(5,*)  ( frac_span(nr,nk),nk=1,nangles(nr))
      read(5,*)  ( exit_angl(nr,nk),nk=1,nangles(nr))
      write(6,*) ' nangles = ',nangles(nr),' angl_typ = ',angl_typ(nr)
      write(6,9048) (frac_span(nr,nk),nk=1,nangles(nr))
      if(angl_typ(nr).eq.'a')
     &     write(6,9045) (exit_angl(nr,nk),nk=1,nangles(nr))
      if(angl_typ(nr).eq.'d')
     &     write(6,9046) (exit_angl(nr,nk),nk=1,nangles(nr))
 9048 format(' fraction of span',/,(10f10.3))
 9045 format(' exit flow angle in degrees ',/,(10f10.3))
 9046 format(' exit flow angle deviation from grid angle, in degrees.',
     & /,(10f10.3))
      write(6,*) ' end throughflow data '
      end if
c ******************************************************************************
c  end tflow
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c
c      call newgrid for the current blade row if requested
c      this reads in more data and generates new streamwise (j) grid points.
c
      if(new_grid.ne.0) then
            call newgrid(jmrow,j1,j2,
     &                   jlerow,jterow,jroths,jrothe,jrotts,jrotte)
      end if
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c
c      set the overall  "j"  markers to be relative to the first point on the first row.
c
      jle(nr) = j1+jlerow-1
      jte(nr) = j1+jterow-1
      j2      = j1 -1 + jmrow
      jrow    = 0
c
      do 100 j = j1,j2
      jm1      = j-1
      if(j.eq.1) jm1=1
      jrow     = jrow+1
c
c********************************************************************************
c      set the rotation of the hub and casing.
c      the default if jroths and jrotts are zero
c      is that the whole hub is rotating at the same speed as
c      the blade row and the casing is stationary.
c
      wrad(j) = rpmrow*pi/30.
      whub(j) = 0.0
      if(jrow.gt.jroths.and.jrow.lt.jrothe) whub(j) = rpmhub*pi/30.
      wtip(j) = 0.0
      if(jrow.gt.jrotts.and.jrow.lt.jrotte) wtip(j) = wrad(j)
c
c     set limits of solid part of the blade elements, ie where no flow.
c     these are from 1 to ktips-1  or from ktipe to km-1
c
      if(ktips(nr).le.0) then
           ks1(nr) = 1
           ks2(nr) = kmm1
      endif
      if(ktips(nr).eq.1) then
           ks1(nr) = ktipe(nr)
           ks2(nr) = kmm1
      endif
      if(ktips(nr).gt.1) then
           ks1(nr)   = 1
           ks2(nr)   = ktips(nr)-1
           ktipe(nr) = km
      endif
c
c     set some markers
c
      ind(j)    = 0
      indle(j)  = 0
      indte(j)  = 0
      indmix(j) = 0
      indmid(j) = 0
      nblade(j) = nblades_in_row
      nrow(j)   = nr
      jmidrw    = 0.5*(jlerow+jterow)
      if(jrow.gt.jlerow.and.jrow.le.jterow) ind(j) = 1
      if(jrow.eq.jlerow) indle(j)    = 1
      if(jrow.eq.jlerow) jled(nr)    = j
      if(jrow.eq.jterow) indte(j)    = 1
      if(jrow.eq.jlerow) indle(j-1) = -1
      if(jrow.eq.jlerow) indle(j-2) = -2
      if(jrow.eq.jmidrw) indmid(j)   = 1
c
c     set the initial guess of static pressure
c
      if(jrow.le.jlerow) pguess(j) = puprow+(plerow-puprow)*
     & (jrow-1)/(jlerow-1)
      if(jrow.gt.jlerow.and.jrow.le.jterow) pguess(j) = 
     &          plerow + (pterow-plerow)*(jrow-jlerow)/(jterow-jlerow)
      if(jrow.gt.jterow) pguess(j) =
     &          pterow + (pdnrow-pterow)*(jrow-jterow)/(jmrow-jterow)
  100 continue
c
       indmix(j2) = 1
       jmix(nr)   = j2
       jstart(nr) = j1
c
c     reset the j index to continue through the next blade row, if any.
c
      j1 = j2 + 1
c
c**********************************************************************************
c        end of data input on this blade row. return to input data on the
c        next row unless this is the last, ie unless nr = nrows.
c**********************************************************************************
c**********************************************************************************
 1550 continue
c**********************************************************************************
c**********************************************************************************
c
c   set the total number of "j" points = jm.
c
      jm         = j2
      jmm1       = jm-1
      jmm2       = jm-2
      indmix(jm) = 0
c
c     check that jm is not too large.
c
      if(jm.gt.jd)  write(6,*) 'stopping because jm too large.',
     &            ' jm= ',jm, ' dimension limit = ',jd
      if(jm.gt.jd) stop
c**********************************************************************************
c     set the start and end points of each stage
c
      jstg_start(1)   = 1
      nstages         = nstage(nrows)
      nstage(nrows+1) = nstages + 1
      jstart(nrows+1) = jm +1
c
      do n = 1,nrows
      nstg   = nstage(n)
      nstgp1 = nstage(n+1)
      if(nstg.ne.nstgp1) then
             jstg_end(nstg)     = jmix(n)
             jstg_start(nstgp1) = jstart(n+1)
      end if
      write(6,*) 'row number ',n,'stage number ',nstg,'jstart ',
     &             jstart(n),'jend ',jmix(n)
      end do
c
      do n = 1,nstages
      write(6,*) ' stage number ',n,' jstart= ',jstg_start(n),
     &           ' jend = ',jstg_end(n)
      end do
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c      automatically shift the blades to make the  x  spacing continuous
c      between the last grid point on one blade and the first grid
c      point on the next one on the hub streamline. if ishift = 1 or > 2.
c**********************************************************************************
c**********************************************************************************
c    start to set the grid between blade rows and make the mixing planes contiguous
c**********************************************************************************
c
c     skip the next part and do not move the blades or grid if ishift = 0 .
      if(ishift.eq.0) go to 140
c     
      if(ishift.ge.2) go to 133
c
c      next section only for ishift = 1
c      if  ishift = 1  automatically shift the blades to make the x spacing continuous
c      between the last grid point on one blade and the first grid point on the next
c      one on the hub streamline but do not change the grid spacings or the radii .
c
      xshift = 0.0
      do 130 j = 2,jm
      if(indmix(j-1).ne.1) go to 129
      diff1  = xsurf(j-1,1) - xsurf(j-2,1)
      diff2  = xsurf(j+1,1) - xsurf(j,1)
      xshift = xsurf(j-1,1) + 0.01*(diff1+diff2) - xsurf(j,1)
  129 continue
      do 131 k=1,km
      xsurf(j,k) = xsurf(j,k) + xshift
      x(j,k)     = x(j,k)     + xshift
  131 continue
  130 continue
c
c     make no more changes to the grid and jump to 140 if ishift = 1.
      if(ishift.eq.1) go to 140
c
c**********************************************************************************
c     re enter here if ishift = 2
  133 continue
c
c******************************************************************************c
c******************************************************************************c
c******************************************************************************c
c     start the loop to grid the gaps between blade rows and also the grids
c     upstream of the first row and downstream of the last row.
c     nrow_gap  is the number of regions to be gridded, including upstream of the
c     first row, all inter row gaps and downstream of the last row. 
c     so there are nrows + 1  gaps .
c
      do 6666 nrow_gap = 1,nrows+1
c
      if(nrow_gap.eq.1) then
           nrw  = 1
           jst  = 1
           jend = jle(1)
      end if
      if(nrow_gap.gt.1.and.nrow_gap.ne.nrows+1) then
          nrw  = nrow_gap
          jst  = jte(nrw-1)
          jend = jle(nrw)
          jmid = jmix(nrw-1)
      end if
      if(nrow_gap.eq.nrows+1) then
           nrw  = nrows + 1
           jst  = jte(nrows)
           jend = jm
      end if
c
      nrwm1 = nrw-1
c******************************************************************************c
c     note that this section sets the grid in the gap between blade rows also
c     upstream of the first row and downstream of the last row.
c     "nrw"  is the row number of the row downstream of the gap being gridded.
c******************************************************************************c
      do 8000 k=1,nsecs_in
c
c    set the meridional distance smerid
c
      smerid(1,k) = 0.0
      do 213 j=2,jm
      xd = xsurf(j,k)  - xsurf(j-1,k)
      rd = rsurf(j,k)  - rsurf(j-1,k)
      smerid(j,k) = smerid(j-1,k) + sqrt(xd*xd+rd*rd)
  213 continue
c
c    if ishift = 2 maintain the meridional curve in the gap between blade rows.
c    if ishift = 3 make the meridional view of the grid linear in the gap between blade rows.
c    if ishift = 4 same as 3 but do not change the hub and casing profiles
c
      dlast= -1.0
      npoints   = 1
      xgap = xsurf(jend,k)  - xsurf(jst,k)
      rgap = rsurf(jend,k)  - rsurf(jst,k)
      gap  = sqrt(xgap*xgap + rgap*rgap)
c
      do 110 j = jst,jend
      xd   = xsurf(j,k)  - xsurf(jst,k)
      rd   = rsurf(j,k)  - rsurf(jst,k)
      proj = (xd*xgap + rd*rgap)
      if(ishift.ge.3) proj = proj/gap
      dist = sqrt(xd*xd+rd*rd)
      if(proj.lt.0.0) go to 110
c
      if(dist.gt.(1.00001*gap)) then
           xint(npoints) = xsurf(jend,k)
           rint(npoints) = rsurf(jend,k)
           npoints = npoints+1
      go to 111
      endif
c
      if(dist.lt.dlast) go to 110
c
      dlast     = dist
c
      if(ishift.eq.3) then
           xint(npoints)  = xsurf(jst,k)   + proj*xgap/gap
           rint(npoints)  = rsurf(jst,k)   + proj*rgap/gap
      endif
c
      if(ishift.eq.2) then
           xint(npoints)  = xsurf(j,k)
           rint(npoints)  = rsurf(j,k)
      endif
c
      if(ishift.eq.4) then
           if(k.eq.1.or.k.eq.nsecs_in) then
                xint(npoints)  = xsurf(j,k)
                rint(npoints)  = rsurf(j,k)
           else
                xint(npoints)  = xsurf(jst,k)   + proj*xgap/gap
                rint(npoints)  = rsurf(jst,k)   + proj*rgap/gap
           end if
      end if
c
      npoints = npoints+1
c
  110 continue
c
  111 continue
c******************************************************************************c
c    smooth the x and r coordinates in the blade to blade gap.
c
      ngaps   = npoints - 1
      nsmooth = 10
      sfgap   = 0.2
      do 113 ns = 1,nsmooth
      do 114 nn = 2,ngaps-1
      xint(nn) = (1.-sfgap)*xint(nn) + sfgap*0.5*(xint(nn-1)+xint(nn+1))
      rint(nn) = (1.-sfgap)*rint(nn) + sfgap*0.5*(rint(nn-1)+rint(nn+1))
  114 continue
  113 continue
c******************************************************************************c
c     re calculate the meridional distance  "sdist "  as it will have been changed
c     by the smoothing.
c
      sdist(1)  = smerid(jst,k)
      do 112 nn = 2,ngaps
      xd = xint(nn) - xint(nn-1)
      rd = rint(nn) - rint(nn-1)
      sdist(nn) = sdist(nn-1) + sqrt(xd*xd+rd*rd)
  112 continue
c
c
      smid     = 0.5*(sdist(1) + sdist(ngaps))
      ste      = sdist(1)
      sle      = sdist(ngaps)
      anglup   = betaup(nrw,k)
      if(nrw.eq.1) then
           angldwn1 = betadwn1(1,k)
           angldwn2 = betadwn2(1,k)
      else
           angldwn1 = betadwn1(nrwm1,k)
           angldwn2 = betadwn2(nrwm1,k)
      end if
      if(nrw.eq.nrows+1) then
           angldwn1 = betadwn1(nrows,k)
           angldwn2 = betadwn2(nrows,k)
      end if
c
c******************************************************************************
c    start to call grid_up and grid_down to set the grid upstream and downstream of
c    all blade rows. also sets the trailing edge cusps.
c******************************************************************************
c    call grid_up to form the grid upstream of the first blade row.
c
            write(6,*)
      if(nrw.eq.1.and.nrw.ne.nrows+1.and.k.eq.1) then
            write(6,*)
     &'****************************************************************'
      write(6,*) ' calling grid_up and grid_down to set the grid and any
     & cusps in the gaps between blade rows.'
            write(6,*)
     &'****************************************************************'
            write(6,*)
            write(6,*)  ' starting to call grid_up and grid_down for row
     & number', 1
            write(6,*)
      end if
c
      if(nrw.eq.1) then
           write(6,*)
     &    'calling grid_up for row 1, k = ',k, ' jstart,jle1= ',jst,jend
           call grid_up(k,jst,jend,ste,sle,ngaps,sdist,xint,
     &                  rint,nextrap_le,anglup,if_angles(1) )
      end if
c
      if(nrw.eq.nrows+1) go to 7900
c
c     call grid_up and grid_down to form the grids upstream and downstream of the
c     mixing plane for interior blade rows.
c
      if( (nrw.ne.1.and.nrw.ne.nrows+1).and.k.eq.1) then
            write(6,*)
            write(6,*)
     &'****************************************************************'
            write(6,*)  ' starting to call grid_up and grid_down for row
     & number', nrw
            write(6,*)
      end if
c
      if (nrw.ne.1)  then
c
       write(6,*)
     &' calling grid_down for row no.',nrw,'k = ',k,' jte,jmid= ',
     &  jst,jmid, 'if_cusp=', if_cusp(nrwm1)
      call grid_down(k,jst,jmid,ste,smid,ngaps,sdist,xint,rint,
     &        nextrap_te,angldwn1,angldwn2,if_cusp(nrwm1),
     &        icusp(nrwm1),lcusp(nrwm1),lcuspup(nrwm1),if_angles(nrwm1))
c
      write(6,*)
     &' calling grid_up for row no.  ',nrw,'k = ',k,'jmid, jle= ',
     &  jmid+1,jend
      call grid_up(k,jmid+1,jend,smid,sle,ngaps,sdist,xint,rint,
     &             nextrap_le,anglup,if_angles(nrw))
c
      endif
c
      go to 8000
c
 7900 continue
c
c     call grid_down to form the grid downstream of the last blade row.
c
      if(nrw.eq.nrows+1.and.k.eq.1) then
            write(6,*)
     &'****************************************************************'
            write(6,*)'calling grid_down the last blade row, row number'
     &                 ,nrows
            write(6,*)
      end if
c
      if(nrw.eq.nrows+1)  then
c
      write(6,*)
     & ' calling grid_down for the last blade row, k = ',k,'jte,jend= ',
     &   jst,jend, 'if_cusp=', if_cusp(nrows)
      call grid_down(k,jst,jend,ste,sle,ngaps,sdist,
     &        xint,rint,nextrap_te,angldwn1,angldwn2,if_cusp(nrows),
     &        icusp(nrows),lcusp(nrows),lcuspup(nrows),if_angles(nrows))
      end if
c
 8000 continue
c
      write(6,*)
c
c*******************************************************************************
c*******************************************************************************
c
c     if  jend = jm gap grid generation completed
c     
      if(jend.eq.jm) then
             write(6,*) 'end of setting the grid and cusps upstream and 
     & downstream of all blade rows'
             write(6,*)
     &'****************************************************************'
             write(6,*)
     &'****************************************************************'
             write(6,*)
      end if
c
c******************************************************************************
c     jump back to 6666 to start on the next row unless jend = jm.
c
 6666 continue
c
c******************************************************************************c
c******************************************************************************c
c      reset smerid  which will have been changed by the new grid.
c
      do 216 k=1,nsecs_in
      do 216 j=2,jm
      xd = xsurf(j,k) - xsurf(j-1,k)
      rd = rsurf(j,k) - rsurf(j-1,k)
  216 smerid(j,k) = smerid(j-1,k) + sqrt(xd*xd+rd*rd)
c
c******************************************************************************
c******************************************************************************
c   end of all geometry input and manipulation.
c
  140 continue
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c
c    read in the spanwise variation in inlet and exit boundary conditions.
c
c******************************************************************************
c******************************************************************************     
c******************************************************************************
c   input the number of data points used for the inlet boundary conditions.
c
      read(5,*)    dummy_input
      read(5,*)    dummy_input
      read(5,*)    kin
      write(6,*) ' number of data points used for the inlet boundary'
      write(6,*) ' conditions = ', kin
      write(6,*)
c
      if(kin.gt.maxki) write(6,*) ' stopping because kin too large.',
     &            ' kin= ',kin,' dimension limit = ',maxki
      if(kin.gt.maxki) stop
c
      kin_mid=ifix(0.5*kin)
c  q3d
      if(km.eq.2) kin_mid = 1
c  end q3d
c
c
c*******************************************************************************
c   input the relative spanwise spacing of the points where the inlet boundary
c   conditions are given'
c
      read(5,*)      dummy_input
      read(5,*)     (fr_in(k),k=1,kin-1)
      write(6,*)  ' relative spacing of the points where the inlet bound
     &ary conditions are given'
      write(6,1700) (fr_in(k),k=1,kin-1)
      write(6,*)
c
c    input the inlet absolute stagnation pressure variation with span
      read(5,*)      dummy_input
      read(5,*)     (po1(k),k=1,kin)
      write(6,*)  ' inlet stagnation pressure variation with span.'
      write(6,1800) (po1(k),k=1,kin)
      po_in_mid = po1(kin_mid)
c
c    input the exit static pressure variation with span. only if  "ipout" = 3 .
      if(ipout.eq.3) then
           read(5,*)     dummy_input
           read(5,*)    (pd(k),k=1,kin)
           write(6,*)  ' exit static pressure variation with span.'
           write(6,1800)(pd(k),k=1,kin)
      end if
c
c     input the inlet absolute stagnation temperature variation with span.
      read(5,*)       dummy_input
      read(5,*)      (to1(k),k=1,kin)
      write(6,*)    ' inlet stagnation temperature variation with span.'
      write(6,1730)  (to1(k),k=1,kin)
      to_in_mid = to1(kin_mid)
c
c   input the inlet tangential velocity variation with span.
      read(5,*)      dummy_input
      read(5,*)     (vtin(k),k=1,kin)
      write(6,*)   ' inlet tangential velocity variation with span.'
      write(6,1730) (vtin(k),k=1,kin)
      vt_in_mid  = vtin(kin_mid)
c
c    input the inlet meridional velocity variation with span.
      read(5,*)      dummy_input
      read(5,*)     (vm1(k),k=1,kin)
      write(6,*)   ' inlet meridional velocity variation with span.'
      write(6,1730) (vm1(k),k=1,kin)
c
c    input the inlet yaw angle variation with span.
      read(5,*)      dummy_input
      read(5,*)     (bs(k),k=1,kin)
      write(6,*)  ' inlet yaw angle variation with span.'
      write(6,1730) (bs(k),k=1,kin)
      yaw_in_mid = bs(kin_mid)
c
c    input the inlet meridional pitch angle variation with span.
      read(5,*)      dummy_input
      read(5,*)     (br(k),k=1,kin)
      write(6,*)   ' inlet meridional pitch angle variation with span.'
      write(6,1730) (br(k),k=1,kin)
      pitch_in_mid  = br(kin_mid)
c
c    end of inputting the spanwise variation of boundary conditions.
c****************************************************************************** 
c****************************************************************************** 
c     input the exit static pressures on the hub and casing.
c     this is the main exit boundary condition.
c
      read(5,*)    dummy_input
      read(5,*)    pdown_hub,pdown_tip
      write(6,*) ' specified downstream pressures on hub and casing=' ,
     &             pdown_hub,pdown_tip
      write(6,*)
c
c*****************************************************************************
c******************************************************************************
c******************************************************************************
c     form the fp(i) and fr(k) and fr_in(k) into a geometric series if requested.
c     by setting fp(3) or fr(3) or fr_in(3) to zero.
c
      num3 = 3
c     
c   q3d
      if(km.eq.2) then
           fr(1)    = 1.0
           fr(2)    = 0.0
      else
c   end q3d
c
c   first for the spanwise grid spacings,  rf(k) .
      if(fr(num3).gt.0.00001) go to 555
      frrat = fr(1)
      frmax = fr(2)
      fr(1) = 1.0
      fr(km)  = 0.0
      do 556 k = 2,kmm1
  556 fr(k) = fr(k-1)*frrat
      do 557 k = 1,kmm1
      frev  = fr(km-k)
      if(frev.lt.fr(k))  fr(k) = frev
      if(fr(k).gt.frmax) fr(k) = frmax
  557 continue
c
  555 continue
c
      end if
c
c    next for the boundary condition spacings, fr_in(k) .
      if(kin.eq.2) then
           fr_in(1) = 1.0
           fr_in(2) = 0.0
      else
c
      if(fr_in(num3).gt.0.00001) go to 560
      frrat = fr_in(1)
      frmax = fr_in(2)
      fr_in(1)   = 1.0
      fr_in(kin) = 0.0
      do 558 k = 2,kin-1
      fr_in(k) = fr_in(k-1)*frrat
  558 continue
      do 559 k = 1,kin-1
      frev  = fr_in(kin-k)
      if(frev.lt.fr_in(k))  fr_in(k) = frev
      if(fr_in(k).gt.frmax) fr_in(k) = frmax
  559 continue
c
  560 continue
c
      end if
c
c   next for the pitchwise grid spacings, fp(i) .
c   throughflow
      if(im.eq.2) then
           fp(1) = 1.0
           fp(2) = 0.0
      else
c    end throughflow
c
      if(fp(num3).gt.0.00001) go to 563
      fprat = fp(1)
      fpmax = fp(2)
      fp(1) = 1.0
      fp(im)= 0.0
      do 564 i=2,imm1
      fp(i) = fp(i-1)*fprat
  564 continue
      do 565 i=1,imm1
      frev  = fp(im-i)
      if(frev.lt.fp(i))  fp(i) = frev
      if(fp(i).gt.fpmax) fp(i) = fpmax
  565 continue
c
  563 continue
c
      end if
c
c***************************************************************************************
c***************************************************************************************
c   call "inpint" to interpolate in the spanwise variation of inflow properties.
c   and set up values on the grid points.
c
      write(6,*) ' calling inpint to set up the inlet flow. '
           call inpint
      write(6,*) ' leaving inpint.'
c
c*************************************************************************************** 
c******************************************************************************
c      read in the mixing length limits if ilos is not zero.
c
      do 571 n = 1, nrows
           xllim_i1(n)  = 0.03
           xllim_im(n)  = 0.03
           xllim_k1(n)  = 0.03
           xllim_km(n)  = 0.03
           xllim_dwn(n) = 0.03
           xllim_up(n)  = 0.02
           xllim_in(n)  = 0.02
           xllim_le(n)  = 0.03
           xllim_te(n)  = 0.04
           xllim_dn(n)  = 0.05
           fsturb(n)    = 1.0
           turbvis_damp(n) = 0.5
  571 continue
c
c      read in the mixing length limits if ilos is not zero.
c
      if(ilos.ne.0) then
c
      read(5,*) dummy_input
c
      do 570 n = 1,nrows

      if(ilos.eq.10) then
c
           read(5,*,err=580)  xllim_i1(n),xllim_im(n),
     &                 xllim_k1(n),xllim_km(n),xllim_dwn(n),xllim_up(n)
      end if
c
      if(ilos.ge.100) then
c
           read(5,*,err=580)  xllim_in(n),xllim_le(n),
     &                xllim_te(n),xllim_dn(n),fsturb(n),turbvis_damp(n)
      end if
c           
  580 continue
c
           write(6,*) ' row number ', n
           if(ilos.eq.10) write(6,*)
     &     'ilos = 10 mixing length limits= ',
     &     xllim_i1(n),xllim_im(n),xllim_k1(n),xllim_km(n),xllim_dwn(n),
     &     xllim_up(n)
           if(ilos.ge.100) write(6,*) 'ilos > 100 mixing length limits',
     &     xllim_in(n),xllim_le(n),xllim_te(n),xllim_dn(n)
           if(ilos.ge.100) write(6,*) 
     &     ' free stream turbulent viscosity ratio',fsturb(n),
     &     ' mixing plane turbulence decay',turbvis_damp(n)
c
  570 continue
c
c     read in a factor to increase the turbulent viscosity for the first nmixup steps.
c
           facmixup = 2.0
           nmixup   = 1000
           read(5,*) dummy_input
           read(5,*, err= 585)     facmixup, nmixup
           if(facmixup.lt.1.0)     facmixup = 1.0
           if(if_restart.ne.0)     facmixup = 1.0
  585 continue
           write(6,*) ' facmixup = ', facmixup,' nmixup = ',nmixup
c
      endif
c
c******************************************************************************
c     write out the value of yplus at the wall if ypluswall > 5.0
c
      if(ypluswall.gt.5.0)  then
           write(6,*)
           write(6,*)
     &    ' wall shear stresses calculated using the ypluswall model. '
           write(6,*) ' yplus at the wall is taken as, ', ypluswall
           cfwall = 1/(ypluswall*ypluswall)
      endif
c******************************************************************************
c  read in the surface roughnesses in microns if  if_rough  >  0  .
c
      if(if_rough.gt.0) then
c
	write(6,*) ' non-hydraulic smooth surfaces specified'
	write(6,*) ' input surface roughness in microns for all 4'
	write(6,*) ' surfaces of each row in turn.'
           read(5,*) dummy_input
      do 586 n = 1,nrows
           read(5,*)  rough_h(n),rough_t(n),rough_l(n),rough_u(n)
c
           write(6,*) ' row number ', n
           write(6,*) ' surface roughnesses in microns ',
     &     rough_h(n),rough_t(n),rough_l(n),rough_u(n)
c
c   change to physical roughness in metres . 
           rough_h(n) = rough_h(n)*1.0e-06
           rough_t(n) = rough_t(n)*1.0e-06
           rough_l(n) = rough_l(n)*1.0e-06
           rough_u(n) = rough_u(n)*1.0e-06
c
  586 continue
c
      else
c
      do 587 n = 1,nrows
           rough_h(n)  = 0.0
           rough_t(n)  = 0.0
           rough_l(n)  = 0.0
           rough_u(n)  = 0.0
  587 continue
c
      endif
c
c******************************************************************************
c******************************************************************************
c
      write(6,*)
      write(6,*)  ' subroutine new_readin completed, input data ok.'
      write(6,*)
c
c
      return
      end
c
c*************************************************************************
c
      subroutine old_readin
c
c       this subroutine reads in the data in
c       ====================

c
      include 'commall-open-19.2'
c
      common/bkodds/
     &           xint(jd),yint(jd),rint(jd),sdist(jd),
     &           icusp(nrs),lcusp(nrs),lcuspup(nrs),
     &           fracnew(jd),betanew(jd),slope(jd),
     &           thickup(jd),thicklow(jd),
     &           xint1(kd),xint2(kd),xint3(kd),xint4(kd)
c
c       start of input section
c       throughout the input section the variables are as follows
c       =====================
c
c       rt_upp(j,k)  is temporary store for blade suction surface co-ord's
c       rt_thick(j,k) is the blade tangential thickness
c       rcyl(k)    is radius of k'th cylindrical surface if input is on
c                  cylindrical surfaces.
c       r(j,1)     is radii of hub
c       r(j,km)    is radii of casing
c
c       =====================
 1000 format(16i5)
 1100 format(8f10.6)
 1200 format(a72)
 1600 format(16i5)
 1700 format(8f10.6)
 1701 format(8f12.6)
 1720 format(8f10.1)
 1730 format(f10.3,4f10.1,2f10.3)
 1800 format(18a4)
c
      pi     = 3.14159265
      degrad = pi/180.
      raddeg = 180./pi
c
      read(5,1200)  title
      write(6,1800) title
c
c************first read in main integer control variables ************
c
      read(5,1000) im,jdum,km,if_rough,nsteps_max,ifcool,ifbleed,nosect,
     &             nrows,ifmix,ishift,kin,nextrap_le,nextrap_te,nchange
      write(6,1600)im,jdum,km,if_rough,nsteps_max,ifcool,ifbleed,nosect,
     &             nrows,ifmix,ishift,kin,nextrap_le,nextrap_te,nchange
c
      read(5,1000) in_vtan,insurf,in_vr,itimst,idummy,ipout,in_flow,
     &             ilos,nlos,if_restart,idumy,ibound,if_repeat
      write(6,1600)in_vtan,insurf,in_vr,itimst,idummy,ipout,in_flow,
     &             ilos,nlos,if_restart,idumy,ibound,if_repeat
c
      in_press = in_vr
      
c
c******************************************************************************
c    set the coefficients   f1, f2,  f3  for the time stepping scheme
c
      write(6,*)
      write(6,*) ' time step type , itimst = ', itimst
      write(6,*)
c
c    set the coefficients for the sss scheme
c
      if(itimst.eq.3.or.itimst.eq.5.or.itimst.eq.6) then
              f1          =  2.0000
              f2          = -1.000
              f3          =  0.00
              f2eff       = -1.0
              nrsmth      =  0
              rsmth       =  0.40
      end if
c
      if(itimst.eq.-3.or.itimst.eq.-5.or.itimst.eq.-6) then
              f1          =  2.0000
              f2          = -1.65
              f3          = -0.65
              f2eff       = -1.0
              nrsmth      =  1
              rsmth       =  0.40
              itimst      =  abs(itimst)
      end if
c
      if(itimst.eq.4.or.itimst.eq.-4) then
              read(5,*) f1, f2eff, f3 , rsmth, nrsmth
              write(6,*) ' f1, f2eff, f3 , rsmth, nrsmth ',
     &                     f1, f2eff, f3 , rsmth, nrsmth
              write(6,*)
              if(f2eff.gt.0.0) then 
              write(6,*)
              write(6,*) ' error,  f2eff  must be negative.'
              write(6,*) ' the sign of the input value will be changed.'
              write(6,*)
              f2eff = -f2eff
              end if
              if(f3.gt.0.0) then 
              write(6,*)
              write(6,*) ' error,  f3  must be negative.'
              write(6,*) ' the sign of the input value will be changed.'
              write(6,*)
              f3  = -f3
              end if
              f2     = f2eff*(1.0 - f3)
              itimst = 3
      end if
c
              write(6,*) ' f1, f2eff, f3 , rsmth, nrsmth ',
     &                     f1, f2eff, f3 , rsmth, nrsmth
              write(6,*)
c
c**********************************************************************************
c      
c  q3d
      if(km.eq.2) ibound = 2
c  end q3d
c
      read (5,1000) ir,jr,kr,irbb,jrbb,krbb
      write(6,1600) ir,jr,kr,irbb,jrbb,krbb
c
c  q3d
      if(km.eq.2) then
           kr   = 1
           krbb = 1
      end if
      if(im.eq.2) then
           ir   = 1
           irbb = 1
      end if
c  end q3d
c
      if_kint = 0
      if(kin.lt.0)     if_kint = 1
      kin = abs(kin)      
      if(kin.eq.0)     kin  = km
      if(kin.ne.km)    if_kint = 1
c
      if(nextrap_le.eq.0)     nextrap_le  = 5
      if(nextrap_te.eq.0)     nextrap_te  = 5
      if(nchange.eq.0) nchange = nsteps_max/4
      if(if_restart.ne.0) nchange = 100
      if(nlos.eq.0)    nlos = 5
      imm1 = im-1
      imm2 = im-2
      if(im.eq.2) imm2 =1
      kmm1 = km-1
      kmm2 = km-2
c
c
c     check that the dimensions are not too large.
c
      if(maxki.lt.kd)  write(6,*)
     &      ' stopping because maxki is less than kd',
     &      ' maxki = ',maxki, ' kd = ',kd
      if(maxki.lt.id)  write(6,*)
     &      ' stopping because maxki is less than id',
     &      ' maxki = ',maxki, ' id = ',id
      if(im.gt.id)  write(6,*) ' stopping because im too large.',
     &            ' im= ',im,  ' dimension limit = ',id
      if(km.gt.kd)  write(6,*) ' stopping because km too large.',
     &            ' km= ',km,  ' dimension limit = ',kd
      if(kin.gt.maxki) write(6,*) ' stopping because kin too large.',
     &            ' kin= ',kin,' dimension limit = ',maxki
c
      if(im.gt.id.or.km.gt.kd)  stop
      if(kin.gt.maxki.or.maxki.lt.kd.or.maxki.lt.id) stop
c
      if(nrows.gt.nrs)  write(6,*) 'stopping because nrows too large.',
     &            ' nrows= ',nrows,' dimension limit = ',nrs
      if(nrows.gt.nrs) stop
c
c     check if there are variable numbers of input sections, i.e.  if nosect is negative
c
      if_sects = 0
      if(nosect.lt.0) then
           if_sects = 1
           nosect   = abs(nosect)
      end if
c
c
c
c***********input the blade geometry on nosect streamwise surfaces******
c*********** the geometry is input for each blade row separately******
c
c            first input the control parameters for the blade row
c            and set the values of some blade row variables
c
      j1        = 1
      ifshroud  = 0
c
      do 1550 nr = 1,nrows
c
c      read in the row title. this is never used but helps to lay out the data
c
      read(5,1200) rowtyp
c
c
c      ktips is the k value of the point where the tip gap starts
c      ktipe is the k value of the point where the tip gap ends.
c      set ktips(nr) = 0  for no tip gap on this row.    
c
c      set ktips(nr) negative to use the shroud leakage model on this
c      row. extra input data is then needed at the end of the data file.
c
      read(5,1000) jmrow,jlerow,jterow,nblades_in_row,
     &   ktips(nr),ktipe(nr),jroths,jrothe,jrotts,jrotte,
     &   new_grid,jtran_i1(nr),jtran_im(nr),jtran_k1(nr),
     &   jtran_km(nr),if_cusp(nr)
c
c     for q3d
      if(km.eq.2) then
           ktips(nr) = 0
           ktipe(nr) = 0
      end if
c     end q3d
c
c     input the number of input sections if this is not constant.
c
      if(if_sects.eq.1) then
            read(5,*) nsecs_row
      else
            nsecs_row = nosect
      end if
c    nsecs_in  must be set as it is used in intpol
            nsecs_in = nsecs_row
c
c     end variable input sections
c   
      jlee = j1+jlerow-1
      jtee = j1+jterow-1
c
      if(ktips(nr).le.0) ktipe(nr) = 0
      if(ktips(nr).lt.0) ifshroud  = 1
c
      if(jroths.eq.0.or.jroths.gt.jmrow) jroths = jmrow
      if(jrothe.eq.0.or.jrothe.gt.jmrow) jrothe = jmrow
      if(jrotts.eq.0.or.jrotts.gt.jmrow) jrotts = jmrow
      if(jrotte.eq.0.or.jrotte.gt.jmrow) jrotte = jmrow
c
c
      write(6,1600) jmrow,jlerow,jterow,nblades_in_row,
     &   ktips(nr),ktipe(nr),jroths,jrothe,jrotts,jrotte,new_grid,
     &   jtran_i1(nr),jtran_im(nr),jtran_k1(nr),jtran_km(nr),if_cusp(nr)
c
c
c     set the cusp generation parameters if if_cusp is not 0.
c     maintain the original grid and no cusp is generated if if_cusp = 0
c     the cusp is centred on the blade centre line if icusp = 0.
c     the cusp makes the i=1 surface continuous on the cusp if  icusp =  1.
c     the cusp makes the i=im surface continuous on the cusp if icusp = -1.
c     the cusp is of length  lcusp and starts lcuspup points before the
c     trailing edge.
c
c
      if(if_cusp(nr).eq.0) write(6,*) '  no cusp specified '
c
      if(if_cusp(nr).eq.1) then
	   read(5,*)    icusp(nr),lcusp(nr),lcuspup(nr)
           write(6,*) ' icusp, lcusp , lcuspup = ',
     &                  icusp(nr),lcusp(nr),lcuspup(nr)
           write(6,*)
      end if
c
c    if  if_cusp = 2 a body force is used to force separation at the trailing edge.
c    the force starts nup_i1 grid points upstream of the te  om the i=1 blade surface, 
c    and at nup_im points upstream on the i=im blade surface. it extends n_wake points
c    downstream of the te. the thickness of the affected layer is determined bt sep_thik
c    , typical value  0.01, and the strength of the body force by sep_drag, typical value 0.99.
c
      if(if_cusp(nr).eq.2) then
           icusp(nr)   = 0
           lcusp(nr)   = 0
           lcuspup(nr) = 0
           read(5,*)  nup_i1(nr),nup_im(nr),n_wake(nr),
     &                  sep_thik(nr),sep_drag(nr)
           write(6,*)' body force used to force trailing edge separtion'
           write(6,*)' cusp body force points, nup_i1, nup_im, n_wake, 
     & sep_thick, sep_drag'
             write(6,*) nup_i1(nr),nup_im(nr),n_wake(nr),
     &                  sep_thik(nr),sep_drag(nr)
             write(6,*)
      end if
c
c
c     read in the rpm, tip clearance and initial guess of inlet and exit 
c     pressures for this blade row.
c
      read(5,1100) rpmrow,puprow,plerow,pterow,pdnrow,fractip(nr),rpmhub
      write(6,1730)rpmrow,puprow,plerow,pterow,pdnrow,fractip(nr),rpmhub
c
      if(fractip(nr).lt.0.0) then 
           read(5,*) fractip1(nr),fractip2(nr)
           write(6,*)' fractip1,  fractip2 = ',fractip1(nr),fractip2(nr)
      else
           fractip1(nr) = fractip(nr)
           fractip2(nr) = fractip(nr)
      end if
c
c
c      fthick(nr,k) is the multiplying factor on the blade thickness so
c      that it can be reduced at the tip. it is assumed to be 1.0 unless
c      input here.
c
      do 2507 k=1,km
 2507 fthick(nr,k)=1.0
      if(ktips(nr).gt.0)  read(5,1100) (fthick(nr,k),k=1,km)
      if(ktips(nr).gt.0) write(6,1700) (fthick(nr,k),k=1,km)
c
c*********************************************************************************
c*********************************************************************************
c*********************************************************************************
c
      j2 = j1 + jmrow - 1
      anglean = 0.0
      if_angles(nr) = 0
c
c*********************************************************************************
c*********************************************************************************
c      now  read in main geometrical data for the blade row 
c      on  nosect  blade sections of the current blade row
c*********************************************************************************
c**********************************************************************************
c
      do 1555 k=1,nsecs_row
c
c
      read(5,1111)  fac1,xshift,if_design,if_restagger,if_lean
      write(6,1111) fac1,xshift,if_design,if_restagger,if_lean
 1111 format(2f10.5,3i10)
c 
c******************************************************************************
c******************************************************************************
c    do not read in an existing blade section if  "if_design"  is non-zero. 
c    jump to 1510 to design a new section
c 
      if(if_design.ne.0) go to 1510
c
c******************************************************************************
c******************************************************************************
c***********xsurf(j,k) is axial coordinate of points on the streamwise surface.
c
      read(5,1100)  (xsurf(j,k),j=j1,j2)
      write(6,1700) (xsurf(j,k),j=j1,j2)
c
c
c     lean the whole blade by an angle anglean if anglean is greater than zero)
c
           read(5,1100)  fac2,tshift
           write(6,1700) fac2,tshift   
c
c***********rt_upp(j,k) is the r-theta coordinate of points on the streamwise
c           surface on the blade surface with largest value of theta.
c           i.e. the upper surface of the blade and the lower surface of the
c                blade to blade passage.
c
      read(5,1100)  (rt_upp(j,k),j=j1,j2)
      write(6,1700) (rt_upp(j,k),j=j1,j2)
c
c
      read(5,1100)  fac3,betaup(nr,k),betadwn1(nr,k),betadwn2(nr,k)
      write(6,1700) fac3,betaup(nr,k),betadwn1(nr,k),betadwn2(nr,k)
c
c
      if(abs(betadwn2(nr,k)).lt.0.0001) betadwn2(nr,k) = betadwn1(nr,k)
      if(abs(betadwn2(nr,k)).gt.0.0001) if_angles(nr) = 1
c
c
c***********rt_thick(j,k) is blade thickness delta r-theta measured
c           in the tangential direction.
c
      read(5,1100)  (rt_thick(j,k),j=j1,j2)
      write(6,1700) (rt_thick(j,k),j=j1,j2)
c
c
c      if insurf =1 or 2 rsurf(j,k) is radii of points on the streamwise
c      surfaces on which data is input.
c      if insurf =0 rcyl(j) is the radius of the cylindrical surfaces on
c       data is input
c
      read(5,1100)   fac4,rshift
      write(6,1700)  fac4,rshift
c
c
      if(insurf.ne.0) read(5,1100)  (rsurf(j,k),j=j1,j2)
      if(insurf.ne.0) write(6,1700) (rsurf(j,k),j=j1,j2)
      if(insurf.eq.0) read(5,1100)  rcyl(k)
      if(insurf.eq.0) write(6,1700) rcyl(k)
c
c******************************************************************************
c******************************************************************************
c     shift the blade section by xshift or tshift
c
      do 1500 j=j1,j2
      if(insurf.eq.0)  rsurf(j,k) = rcyl(k)
      xsurf(j,k)    = fac1*(xshift + xsurf(j,k))
      rt_thick(j,k) = fac3*rt_thick(j,k)
      rsurf(j,k)    = fac4*(rshift + rsurf(j,k))
      rt_upp(j,k)   = fac2*(tshift + rt_upp(j,k))
 1500 continue
c 
c******************************************************************************
c******************************************************************************
 1510 continue
c******************************************************************************
c******************************************************************************
c     jdd addition to allow redesign of this blade section.
c     changed to use subroutine re_design. august 2016.
c
      if(if_design.ne.0) then
          write(6,*)'calling re_design to generate a new blade geometry'
          call re_design(nr,k,j1,j2,jlerow,jterow) 
          betaup(nr,k)   = 0.0
          betadwn1(nr,k) = 0.0
          betadwn2(nr,k) = 0.0
          if_angles(nr)  = 0
      end if
c   
c    end of option to redesign the blade section
c******************************************************************************
c******************************************************************************
c
c     restagger the blade section if "if_restagger" is non-zero.
c
      if(if_restagger.ne.0) then
      write(6,*) ' calling subroutine restagger to rotate the blade '
      call restagger(k,j1,j2,jlerow,jterow,rotate,fracx_rot)
      end if
c
c******************************************************************************
c******************************************************************************
c 
c    lean the blade section if "if_lean"  is not zero.
c
      if(if_lean.ne.0) then
      write(6,*) ' calling subroutine lean to lean the blade section '
      call lean(k,j1,j2,anglean)
      end if
c
c******************************************************************************
c******************************************************************************
c   q3d
c     call set ssthick to set the stream surface thickness if doing a q3d
c     blade to blade calculation. then end the geometry input for this blade row
c     as only a single stream surface is used .
c
      if(km.eq.2) then
           call set_ssthick(j1,j2)
c    jump out of the do 1555 loop as only one stream surface is needed.
c    skip the separate hub and casing geometry input if doing a q3d calculation.
           go to 1504
      end if
c   end q3d
c
c******************************************************************************
c******************************************************************************
c     end if input for this blade section
c
 1555 continue
c
c
c******************************************************************************
c******************************************************************************
c  tflow
c******************************************************************************
c******************************************************************************
c  if a throughflow calculation read in the required exit flow angle or deviation 
c  angle for this blade row. angl_typ = 'a' for exit angle, = 'd' for deviation angle.
c
      if(im.eq.2) then  
      read(5,*)    dummy_input
      write(6,*) ' inputting the throughflow data '
      read(5,*)    angl_typ(nr), nangles(nr)
      read(5,*)  ( frac_span(nr,nk),nk=1,nangles(nr))
      read(5,*)  ( exit_angl(nr,nk),nk=1,nangles(nr))
      write(6,*) ' nangles = ',nangles(nr),' angl_typ = ',angl_typ(nr)
      write(6,9048) (frac_span(nr,nk),nk=1,nangles(nr))
      if(angl_typ(nr).eq.'a')
     &     write(6,9045) (exit_angl(nr,nk),nk=1,nangles(nr))
      if(angl_typ(nr).eq.'d')
     &     write(6,9046) (exit_angl(nr,nk),nk=1,nangles(nr))
 9048 format(' fraction of span',/,(10f10.3))
 9045 format(' exit flow angle in degrees ',/,(10f10.3))
 9046 format(' exit flow angle deviation from grid angle, in degrees.',
     & /,(10f10.3))
      write(6,*) ' end throughflow data '
      end if
c ******************************************************************************
c  end tflow
c******************************************************************************
c******************************************************************************
c
c          input the hub and casing geometry.
c          if insurf = 2 the hub and casing  are taken as the
c          first and last streamwise surfaces. otherwise read in
c          hub and casing coordinates at the ends of the quasi orthogonals.
c          this cannot be used if km = 2 for a q3d calculation.
c
      if(insurf.ne.2) then
c
      read(5,1100)  (x(j,1),j=j1,j2)
      read(5,1100)  (r(j,1),j=j1,j2)
      write(6,1700) (x(j,1),j=j1,j2)
      write(6,1700) (r(j,1),j=j1,j2)
      read(5,1100)  (x(j,km),j=j1,j2)
      read(5,1100)  (r(j,km),j=j1,j2)
      write(6,1700) (x(j,km),j=j1,j2)
      write(6,1700) (r(j,km),j=j1,j2)
      do 1505 j=j1,j2
      x(j,1)  = (xshift + x(j,1))*fac1
      r(j,1)  = (rshift + r(j,1))*fac4
      x(j,km) = (xshift + x(j,km))*fac1
 1505 r(j,km) = (rshift + r(j,km))*fac4
c
      else
c
      do 1503 j = j1,j2
      r(j,1)    = rsurf(j,1)
      r(j,km)   = rsurf(j,nsecs_row)
      x(j,1)    = xsurf(j,1)
      x(j,km)   = xsurf(j,nsecs_row)
 1503 continue
c
      end if
c
c*********************************************************************************
c*********************************************************************************
c     interpolate extra input sections if  nsecs_row  not equal to   nosect
c     this routine provided by s gallimore.
c
      if(nsecs_row.ne.nosect)then
        do j = j1,j2
	  fspan(1) = 0.0
          do k = 1,nsecs_row
            xint1(k) = xsurf(j,k)
            xint2(k) = rsurf(j,k)
            xint3(k) = rt_upp(j,k)
            xint4(k) = rt_thick(j,k)
            if (k.gt.1) then
              rd = rsurf(j,k)-rsurf(j,k-1)
              xd = xsurf(j,k)-xsurf(j,k-1)
	    fspan(k) = fspan(k-1) + sqrt(xd*xd+rd*rd)
            end if
          end do
c
          xarg = 0.0
          do k = 1,nosect
            if(k.gt.1) xarg = xarg+(fspan(nsecs_row)/(nosect-1))
            call intp(nsecs_row,fspan,xint1,xarg,xsurf(j,k))
            call intp(nsecs_row,fspan,xint2,xarg,rsurf(j,k))
            call intp(nsecs_row,fspan,xint3,xarg,rt_upp(j,k))
            call intp(nsecs_row,fspan,xint4,xarg,rt_thick(j,k))
          end do
        end do
      end if
c
c****************************************************************************
c  re-enter here if km = 2
 1504 continue
c
c****************************************************************************
c      call newgrid for the current blade row if requested
c      this reads in more data and generates new streamwise (j) grid points.
c
      if(new_grid.ne.0) call newgrid(jmrow,j1,j2,
     &                jlerow,jterow,jroths,jrothe,jrotts,jrotte)
c
c*****************************************************************************
c
c      set the j markers, initial guess of pressure, etc
c
      jle(nr) = j1+jlerow-1
      jte(nr) = j1+jterow-1
      j2      = j1 -1 + jmrow
      jrow    = 0
c
      do 100 j = j1,j2
      jm1      = j-1
      if(j.eq.1) jm1=1
      jrow     = jrow+1
c
c      set the rotation of the hub and casing.
c      the default if jroths and jrotts are zero
c      is that the whole hub is rotating at the same speed as
c      the blade row and the casing is stationary.
c
      wrad(j) = rpmrow*3.14159/30.
      whub(j) = 0.0
      if(jrow.gt.jroths.and.jrow.lt.jrothe) whub(j) = rpmhub*3.14159/30.
      wtip(j) = 0.0
      if(jrow.gt.jrotts.and.jrow.lt.jrotte) wtip(j) = wrad(j)
c
c     set limits of solid part of the blade elements, ie where no flow.
c     these are from 1 to ktips-1  or from ktipe to km-1
c
      if(ktips(nr).le.0) then
      ks1(nr) = 1
      ks2(nr) = kmm1
      endif
      if(ktips(nr).eq.1) then
      ks1(nr) = ktipe(nr)
      ks2(nr) = kmm1
      endif
      if(ktips(nr).gt.1) then
      ks1(nr)   = 1
      ks2(nr)   = ktips(nr)-1
      ktipe(nr) = km
      endif
c
      ind(j)    = 0
      indle(j)  = 0
      indte(j)  = 0
      indmix(j) = 0
      indmid(j) = 0
      nblade(j) = nblades_in_row
      nrow(j)   = nr
      jmidrw    = 0.5*(jlerow+jterow)
      if(jrow.gt.jlerow.and.jrow.le.jterow) ind(j)=1
      if(jrow.eq.jlerow) indle(j)    = 1
      if(jrow.eq.jlerow) jled(nr)    = j
      if(jrow.eq.jterow) indte(j)    = 1
      if(jrow.eq.jlerow) indle(j-1) = -1
      if(jrow.eq.jlerow) indle(j-2) = -2
      if(jrow.eq.jmidrw) indmid(j)   = 1
      if(jrow.le.jlerow) pguess(j) = puprow+(plerow-puprow)*
     & (jrow-1)/(jlerow-1)
      if(jrow.gt.jlerow.and.jrow.le.jterow) pguess(j) = 
     &          plerow + (pterow-plerow)*(jrow-jlerow)/(jterow-jlerow)
      if(jrow.gt.jterow) pguess(j) =
     &          pterow + (pdnrow-pterow)*(jrow-jterow)/(jmrow-jterow)
  100 continue
c
       indmix(j2) = 1
       jmix(nr)   = j2
       jstart(nr) = j1
c
c     reset the j index to continue through the next blade row, if any.
c
      j1 = j2 + 1
c
c        end of data input on this blade row. return to input data on the
c        next row unless this is the last, ie unless nr=nrows.
c
 1550 continue
c
c       end of input of blade geometry data
c
      jm         = j2
      jmm1       = jm-1
      jmm2       = jm-2
      indmix(jm) = 0
c
c
c     check that jm is not too large.
c
      if(jm.gt.jd)  write(6,*) 'stopping because jm too large.',
     &            ' jm= ',jm,' dimension limit = ',jd
      if(jm.gt.jd) stop
c
c**********************************************************************************
c**********************************************************************************
c    start to set the grid between blade rows and make the mixing planes contiguous
c**********************************************************************************
c
c     skip the next part and do not move the blades or grid if ishift = 0 .
      if(ishift.eq.0) go to 140
c     
      if(ishift.ge.2) go to 133
c
c      next section only for ishift = 1
c      if  ishift = 1  automatically shift the blades to make the x spacing continuous
c      between the last grid point on one blade and the first grid point on the next
c      one on the hub streamline but do not change the grid spacings.
c
      xshift = 0.0
      do 130 j = 2,jm
      if(indmix(j-1).ne.1) go to 129
      diff1  = xsurf(j-1,1) - xsurf(j-2,1)
      diff2  = xsurf(j+1,1) - xsurf(j,1)
      xshift = xsurf(j-1,1) + 0.01*(diff1+diff2) - xsurf(j,1)
  129 continue
      do 131 k=1,km
      xsurf(j,k) = xsurf(j,k) + xshift
      x(j,k)     = x(j,k)     + xshift
  131 continue
  130 continue
c
c     make no more changes to the grid if ishift = 1.
      if(ishift.eq.1) go to 140
c
c**********************************************************************************
c     re enter here if ishift = 2
  133 continue
c
c********************************************************************************
c      start to make the grid spacing vary geometrically between blade rows 
c      if ishift >= 2.
c*******************************************************************************
c
      nr   = 1
      jst  = 1
      jend = jle(1)
c
c******************************************************************************
c
 6666 continue
c
c******************************************************************************
      do 210 k=1,nosect
c
      smerid(1,k) = 0.0
      do 213 j=2,jm
      xd = xsurf(j,k)  - xsurf(j-1,k)
      rd = rsurf(j,k)  - rsurf(j-1,k)
  213 smerid(j,k) = smerid(j-1,k) + sqrt(xd*xd+rd*rd)
c
c    if ishift = 2 maintain the meridional curve in the gap between blade rows.
c    if ishift = 3 make the meridional view of the grid linear in the gap between blade rows.
c    if ishift = 4 same as 3 but do not change the hub and casing profiles
c
      dlast = -1.0
      np = 1
      xgap = xsurf(jend,k)  - xsurf(jst,k)
      rgap = rsurf(jend,k)  - rsurf(jst,k)
      gap  = sqrt(xgap*xgap + rgap*rgap)
c
      do 110 j = jst,jend
      xd   = xsurf(j,k)  - xsurf(jst,k)
      rd   = rsurf(j,k)  - rsurf(jst,k)
      proj = (xd*xgap + rd*rgap)

      if(ishift.ge.3) proj = proj/gap

      dist = sqrt(xd*xd+rd*rd)
      if(proj.lt.0.0) go to 110
      if(dist.gt.(1.00001*gap)) then
      xint(np) = xsurf(jend,k)
      rint(np) = rsurf(jend,k)
      np = np+1
      go to 111
      endif
      if(dist.lt.dlast) go to 110
      dlast     = dist

      if(ishift.eq.3) then
           xint(np)  = xsurf(jst,k)   + proj*xgap/gap
           rint(np)  = rsurf(jst,k)   + proj*rgap/gap
      endif

      if(ishift.eq.2) then
           xint(np)  = xsurf(j,k)
           rint(np)  = rsurf(j,k)
      endif

      if(ishift.eq.4) then
           if(k.eq.1.or.k.eq.nosect) then
                xint(np)  = xsurf(j,k)
                rint(np)  = rsurf(j,k)
           else
                xint(np)  = xsurf(jst,k)   + proj*xgap/gap
                rint(np)  = rsurf(jst,k)   + proj*rgap/gap
           end if
      end if
c
      np = np+1
c
  110 continue
c
  111 continue
c
c    smooth the x and r coordinates in the blade to blade gap
c
      ngap    = np - 1
      nsmooth = 10
      sfgap   = 0.2
      do 113 ns = 1,nsmooth
      do 114 nn = 2,ngap-1
      xint(nn) = (1.-sfgap)*xint(nn) + sfgap*0.5*(xint(nn-1)+xint(nn+1))
      rint(nn) = (1.-sfgap)*rint(nn) + sfgap*0.5*(rint(nn-1)+rint(nn+1))
  114 continue
  113 continue
c
c
      sdist(1)  = smerid(jst,k)
      do 112 nn = 2,ngap
      xd = xint(nn) - xint(nn-1)
      rd = rint(nn) - rint(nn-1)
  112 sdist(nn) = sdist(nn-1) + sqrt(xd*xd+rd*rd)
c
c
      smid   = 0.5*(sdist(1) + sdist(ngap))
      ste    = sdist(1)
      sle    = sdist(ngap)
      if(nr.le.nrows) anglup = betaup(nr,k)
c
      if(nr.ne.1)  then
                   angldwn1 = betadwn1(nr-1,k)
                   angldwn2 = betadwn2(nr-1,k)
      end if
c
c******************************************************************************
c    start to call grid_up and grid_down to set the grid upstream and downstream of
c    all blade rows.
c******************************************************************************
c    call grid_up to form the grid upstream of the first blade row.
c
      if(nr.eq.1.and.k.eq.1) then
            write(6,*)
            write(6,*) ' starting new blade row, row number', 1
            write(6,*)
      end if
      if(nr.eq.1) write(6,*)
     & ' calling grid_up for row 1, k = ',k, ' jstart,jle1= ',jst,jend
c
      if(nr.eq.1)  call grid_up(k,jst,jend,ste,sle,ngap,sdist,xint,
     &             rint,nextrap_le,anglup,if_angles(nr))
c
c     call grid_up and grid_down to form the grids upstream and downstream of the
c     mixing plane for interior blade rows.
c
      if( (nr.ne.1.and.nr.ne.nrows+1).and.k.eq.1) then
            write(6,*)
            write(6,*) ' starting new blade row, row number',nr
            write(6,*)
      end if
      if(nr.ne.1.and.nr.ne.nrows+1) then
c
      write(6,*)
     &' calling grid_down for row no.',nr,'k = ',k,' jte,jmid= ',
     &  jst,jmid
c
      call grid_down(k,jst,jmid,ste,smid,ngap,sdist,xint,rint,
     &         nextrap_te,angldwn1,angldwn2,if_cusp(nr-1),icusp(nr-1),
     &         lcusp(nr-1),lcuspup(nr-1),if_angles(nr-1) )
      write(6,*)
     &' calling grid_up for row no.  ',nr,'k = ',k,'jmid, jle= ',
     &  jmid+1,jend
      call grid_up(k,jmid+1,jend,smid,sle,ngap,sdist,xint,rint,
     &             nextrap_le,anglup,if_angles(nr) )
c
      endif
c
c     call grid_down to form the grid downsream of the last blade row.
c
      if(nr.eq.nrows+1.and.k.eq.1) then
            write(6,*)
            write(6,*) ' starting new blade row, row number',nrows
            write(6,*)
      end if
c
      if(nr.eq.nrows+1)  then
           angldwn1 = betadwn1(nrows,k)
           angldwn2 = betadwn2(nrows,k)
      write(6,*)
     & ' calling grid_down for the last blade row, k = ',k,'jte,jend= ',
     &   jst,jend
c
      call grid_down(k,jst,jend,ste,sle,ngap,sdist,
     &        xint,rint,nextrap_te,angldwn1,angldwn2,if_cusp(nrows),
     &        icusp(nrows),lcusp(nrows),lcuspup(nrows),if_angles(nrows))
c
      end if
c
      if(abs(anglup).gt.0.0.or.abs(angldwn1).gt.0.0) then
           write(6,*)
           write(6,*) 'row no ',nr,'section no ',k,'anglup= ',anglup,
     &                'angldwn1 & 2= ',angldwn1,angldwn2,' degrees.'
           write(6,*)
      end if
c
c     end of setting the grid upstream and downsream of all blade rows.
c
  210 continue
c
      write(6,*)
c
c********************************************************************************
c
      if(jend.eq.jm) go to 200
c
      nr = nr + 1
      if(nr.ne.nrows+1) then
      jst    = jte(nr-1)
      jend   = jle(nr)
      jmid   = jmix(nr-1)
      go to 6666
      endif
      if(nr.eq.nrows+1) then
      jst  = jte(nrows)
      jend = jm
      go to 6666
      endif
c
  200 continue
c
c      reset smerid
c
      do 216 k=1,nosect
      do 216 j=2,jm
      xd = xsurf(j,k) - xsurf(j-1,k)
      rd = rsurf(j,k) - rsurf(j-1,k)
  216 smerid(j,k) = smerid(j-1,k) + sqrt(xd*xd+rd*rd)

c      reset the hub and casing coordinates to new stream surface values if insurf = 2.
c
      if(insurf.eq.2) then
c
       do 214 j=1,jm
       x(j,1)  = xsurf(j,1)
       r(j,1)  = rsurf(j,1)
       x(j,km) = xsurf(j,nosect)
       r(j,km) = rsurf(j,nosect)
  214 continue
c
      else
c
c    reset the hub and casing coordinates if insurf is not = 2.
c
      do 218       nr = 1,nrows+1
      if(nr.ne.1)  jt = jte(nr-1)
      if(nr.eq.1)  jt = 1
      if(nr.ne.nrows+1) jl = jle(nr)
      if(nr.eq.nrows+1) jl = jm
      do 215 k = 1,km,kmm1
      m = 1
      if(k.eq.km) m = nosect
      do 333      j = jt,jl
      nj            = j - jt + 1
      if(nj.eq.1) then
      sdist(nj) = 0.0
      else
      xd = x(j,k) - x(j-1,k)
      rd = r(j,k) - r(j-1,k)
      sdist(nj) = sdist(nj-1) + sqrt(xd*xd+rd*rd)
      endif
      xint(nj) = x(j,k)
      rint(nj) = r(j,k)
  333 continue
c
      gap = sdist(nj) - sdist(1)
      do 334 nn = 2,nj
  334 sdist(nn) = (sdist(nn)-sdist(1))/gap
c
      sgap = smerid(jl,m) - smerid(jt,m)
      do 217 j = jt,jl
      arg = (smerid(j,m) - smerid(jt,m))/sgap
      call linint(nj,sdist,xint,arg,x(j,k))
      call linint(nj,sdist,rint,arg,r(j,k))
  217 continue
  215 continue
  218 continue
c
      endif
c
  140 continue
c
c   end of geometry input and manipulation.
c******************************************************************************
c******************************************************************************
c
c
c******************************************************************************
c******************************************************************************
c**********read in gas constants ,time step length,smoothing factor,etc.
c
      write(6,*)'reading in the gas constants,time step length,smoothing 
     &factor,etc. '
c
      ifgas = 0
      read(5,1100)  cp,ga,cfl,sftin,sfxin,machlim
      write(6,1701) cp,ga,cfl,sftin,sfxin,machlim
      if(machlim.lt.1.0) machlim = 2.0
c
c    use real gas properties if cp is input as negative.
c    typical values for combustion products  are:cp1 = 1272.5, cp2 = 0.2125, cp3 = 0.000015625, rgas = 287.15
c    at  tref = 1400 k.
c
      if(cp.lt.0.0) then
      read(5,*) cp1, cp2, cp3, tref, rgas
c
      write(6,*) ' ideal gas properties read in '
      write(6,*) ' cp1, cp2, cp3, tref, rgas =',cp1,cp2,cp3,tref,rgas
c
      cpgas  = cp1
      gagas  = cp1/(cp1 - rgas)
      cp     = cp1
      ga     = gagas
      cv     = cp/ga
      ifgas  = 1
      call set_coeffs
      end if
c
c******************************************************************************
c******************************************************************************
c
      read(5,1100)  dampin,dumm,fblk1,fblk2,fblk3,sfexit,conlim,rfin
      if(sfexit.gt.0.0001) read(5,*) nsfexit
      write(6,1700) dampin,dumm,fblk1,fblk2,fblk3,sfexit,conlim,rfin
      if(sfexit.gt.0.0001) write(6,*) nsfexit
c
      if(dampin.lt.0.001)    dampin = 10.
      if(conlim.lt.0.000001) conlim = 0.005
      if(rfin.lt.0.000001)   rfin = 0.10
      rfthrotl = 0.0
c
c******************************************************************************
c******************************************************************************
c      read in initial guess of upstream and downstream pressures
c      ,stagnation temperature , velocities,flow directions,etc.
c
      read(5,1100)  puphub,puptip,pdown_hub,pdown_tip,plate_loss,
     &              throttle_exit,dumm,f_pdown
      if(throttle_exit.gt.0.001) read(5,*) throttle_pres, throttle_mas,
     &                                      rfthrotl
      write(6,1720) puphub,puptip,pdown_hub,pdown_tip,plate_loss,
     &              throttle_exit,dumm,f_pdown
      if(throttle_exit.gt.0.001) write(6,*) throttle_pres,throttle_mas,
     &                                       rfthrotl
c
      rfthrotl1 = 1.0 - rfthrotl
c
      kmid=ifix(0.5*kin)
c  q3d
      if(km.eq.2) kmid = 1
c  end q3d
c
c******************************************************************************
c******************************************************************************
      write(6,*)
      write(6,*) ' reading in the inlet boundary conditions '
c
      read(5,1100)  (po1(k),k=1,kin)
      write(6,1720) (po1(k),k=1,kin)
      po_in_mid = po1(kmid)
c
      if(ipout.eq.3) read(5,1100)(pd(k),k=1,kin)
      if(ipout.eq.3) write(6,1720)(pd(k),k=1,kin)
c
      read(5,1100)  (to1(k),k=1,kin)
      write(6,1700) (to1(k),k=1,kin)
      to_in_mid = to1(kmid)
c
      read(5,1100)  (vtin(k),k=1,kin)
      write(6,1700) (vtin(k),k=1,kin)
      vt_in_mid  = vtin(kmid)
c
      read(5,1100)  (vm1(k),k=1,kin)
      write(6,1700) (vm1(k),k=1,kin)
c
      read(5,1100)  (bs(k),k=1,kin)
      write(6,1700) (bs(k),k=1,kin)
      yaw_in_mid = bs(kmid)
c
      read(5,1100)  (br(k),k=1,kin)
      write(6,1700) (br(k),k=1,kin)
      pitch_in_mid  = br(kmid)
c
      read(5,1100)  (fr_in(k),k=1,kin-1)
      write(6,1700) (fr_in(k),k=1,kin-1)
c
      read(5,1100)  (fp(i),i=1,imm1)
      write(6,1700) (fp(i),i=1,imm1)
c
      read(5,1000)  (nout(l),l=1,5)
      write(6,1600) (nout(l),l=1,5)
c
      read(5,1610)  (iout(i),i=1,13)
      write(6,1610) (iout(i),i=1,13)
c
      read(5,1610)  (kout(k),k=1,km)
      write(6,1610) (kout(k),k=1,km)
 1610 format(40i2)
c
c***************************************************************************************
c***************************************************************************************
c     form the fp(i) and fr_in(k) into a geometric series if requested.
c     by setting fp(3) or fr_in(3) to zero.
c
      num3 = 3
c
c   q3d
      if(km.eq.2) then
           fr_in(1) = 1.0
           fr_in(2) = 0.0
      else
c   end q3d
c
      if(fr_in(num3).gt.0.00001) go to 555
      frrat = fr_in(1)
      frmax = fr_in(2)
      fr_in(1) = 1.0
      fr_in(km)  = 0.0
      fr_in(kin) = 0.0
      do 556 k = 2,kin-1
           fr_in(k) = fr_in(k-1)*frrat
  556 continue
      do 557 k = 1,kin-1
           frev  = fr_in(kin-k)
           if(frev.lt.fr_in(k))  fr_in(k) = frev
           if(fr_in(k).gt.frmax) fr_in(k) = frmax
  557 continue
c
  555 continue
c
      end if
c
c   next for the pitchwise grid spacings, fp(i) .
c   throughflow
      if(im.eq.2) then
            fp(1) = 1.0
            fp(2) = 0.0
      else
c  end throughflow
c
      if(fp(num3).gt.0.00001) go to 563
      fprat = fp(1)
      fpmax = fp(2)
      fp(1) = 1.0
      fp(im)= 0.0
      do 564 i=2,imm1
      fp(i) = fp(i-1)*fprat
  564 continue
      do 565 i=1,imm1
      frev  = fp(im-i)
      if(frev.lt.fp(i))  fp(i) = frev
      if(fp(i).gt.fpmax) fp(i) = fpmax
  565 continue
c
  563 continue
c
      end if
c
c***************************************************************************************
c***************************************************************************************
c
c     read in in the specified inlet flow and relaxation factor
c     if in_flow not = zero
c
      if(in_flow.ne.0) then
           read(5,1100)  flowin, rflow
           write(6,*)  ' flowin, rflow ', flowin,rflow
           flowin = flowin/nblade(1)
      end if
c
c******************************************************************************
c******************************************************************************
c      read in data for viscous flow modelling.
c
      reyno       = 500000.0
      rf_vis      = 0.5
      ftrans      = 0.0001
      fac_4th     = 0.8
      prandtl     = 1.0
      ypluswall   = 0.0
      turbvis_lim = 3000.0
c
      if(ilos.ne.0) then
      read(5,*,end=448,err=448)reyno,rf_vis,ftrans,fac_4th,turbvis_lim,
     &                         prandtl,ypluswall
  448 continue
c
      if(turbvis_lim.lt.0.1) turbvis_lim = 3000.0
c
      write(6,*)
      write(6,*) ' reynolds number                       = ',reyno
      write(6,*) ' viscous term relaxation factor        = ',rf_vis
      write(6,*) ' transition factor, ftrans             = ',ftrans
      write(6,*) ' proportion of fourth order smoothing  = ', fac_4th
      write(6,*) ' limit on turbulent/laminar viscosity  = ',turbvis_lim
      write(6,*) ' prandtl number                        = ',prandtl
      write(6,*) ' ypluswall - if used (not often used)  = ',ypluswall
      write(6,*) ' mach number  limiter, usually 2.0     = ',machlim
      write(6,*)
c
      if(ilos.ge.200)then
           fac_stmix = 0.0
           fac_st0   = 1.0
           fac_st1   = 1.0
           fac_st2   = 1.0
           fac_st3   = 1.0
           fac_sfvis = 2.0
           read(5,*,end = 450,err= 450) fac_stmix, fac_st0, fac_st1,
     &                                  fac_st2, fac_st3, fac_sfvis,
     &                                  fac_vort,fac_pgrad 
  450 continue
c
           write(6,*)
           write(6,*) ' spalart-allmaras turbulence model is being used'
           write(6,*) ' the s_a source term multipliers are ',
     &                  fac_stmix, fac_st0, fac_st1, fac_st2, fac_st3,
     &                  fac_vort,fac_pgrad 
           write(6,*) ' the turbulent viscosity smoothing factor is ', 
     &                  fac_sfvis
           write(6,*)
      end if
c
      end if 
c
      write(6,*)
c
c***************************************************************************************
c**************************************************************************************
c     read in the mixing plane parameters
c
      rfmix    = 0.025
      fsmthb   = 1.0
      fextrap  = 0.95
      fangle   = 0.95
c
      if(ifmix.ne.0) then
      read(5,*,end=449,err=449)  rfmix,fextrap,fsmthb,
     &                           fangle
  449 continue
      write(6,*)' the mixing plane parameters are '
      write(6,*)' rfmix, fextrap, fsmthb ,fangle',
     &            rfmix, fextrap, fsmthb, fangle
      write(6,*)
      end if
c
c***************************************************************************************
c***************************************************************************************
c    set the spanwise grid spacings,  "fr(k)" .
c   these are the same as the inlet bc spacings, "fr_in(k)" ,  if "if_kint" = 0.
c
      if(if_kint.eq.1) then
c
c     read in the relative spanwise spacing of the grid points.
c     note that this is formatted input.
      read(5,10)    (fr(k),k=1,kmm1)
   10 format(8f10.5)
c
      write(6,*)
      write(6,*)       ' the relative spacings of the grid points in the
     &  spanwise direction are: '
      write(6,1700)(fr(k),k=1,kmm1)
      write(6,*)
c
c   make the spanwise spacings into a geometrical progression if fr(3) = 0 .
      num3 = 3
c
      if(km.eq.2) then
           fr(1) = 1.0
           fr(2) = 0.0
      else
c
      if(fr(num3).gt.0.00001) go to 652
      frrat = fr(1)
      frmax = fr(2)
      fr(1)=1.0
      do 650 k=2,km-1
      fr(k) = fr(k-1)*frrat
  650 continue
      do 651 k=1,kmm1
      frev = fr(km-k)
      if(frev.lt.fr(k))  fr(k) = frev
      if(fr(k).gt.frmax) fr(k) = frmax
  651 continue
c
  652 continue
c
      end if
c
c
      write(6,*) ' calling inpint to set up a new inlet flow. '
c
            call inpint
c
      else
c
      do 655 k=1,kmm1
           fr(k) = fr_in(k)
  655 continue
c  end of  if if_kint = 1 loop
      endif
c
c***************************************************************************************
c***************************************************************************************
c      read in the mixing length limits if ilos is not zero.
c
      do 2344 n = 1,nrows
           xllim_i1(n)  = 0.03
           xllim_im(n)  = 0.03
           xllim_k1(n)  = 0.03
           xllim_km(n)  = 0.03
           xllim_dwn(n) = 0.03
           xllim_up(n)  = 0.02
           xllim_in(n)  = 0.02
           xllim_le(n)  = 0.03
           xllim_te(n)  = 0.04
           xllim_dn(n)  = 0.05
           fsturb(n)    = 1.0
           turbvis_damp(n) = 0.5
 2344 continue
c
c
      if(ilos.ne.0) then
c
      do 2345 n = 1,nrows

           if(ilos.lt.100) read(5,*,err=2346)  xllim_i1(n),xllim_im(n),
     &                 xllim_k1(n),xllim_km(n),xllim_dwn(n),xllim_up(n)
           if(ilos.ge.100) read(5,*,err=2346)  xllim_in(n),xllim_le(n),
     &                xllim_te(n),xllim_dn(n),fsturb(n),turbvis_damp(n)
c           
 2346 continue
c
           write(6,*) ' row number ', n
           if(ilos.lt.100) write(6,*)'ilos = 9/10 mixing length limits',
     &     xllim_i1(n),xllim_im(n),xllim_k1(n),xllim_km(n),xllim_dwn(n),
     &     xllim_up(n)
           if(ilos.ge.100) write(6,*) 'ilos > 100 mixing length limits',
     &     xllim_in(n),xllim_le(n),xllim_te(n),xllim_dn(n)
           if(ilos.ge.100) write(6,*) 
     &     ' free stream turbulent viscosity ratio',fsturb(n),
     &     ' mixing plane turbulence decay', turbvis_damp(n)
c
 2345 continue
c
c******************************************************************************
c******************************************************************************
c     read in a factor to increase the turbulent viscosity for the first nmixup steps.
c
           facmixup = 2.0
           nmixup   = 1000
           read(5,*, err= 2341)    facmixup, nmixup
           if(facmixup.lt.1.0) facmixup = 1.0
           if(if_restart.ne.0)     facmixup = 1.0
 2341 continue
           write(6,*) ' facmixup = ', facmixup,' nmixup = ',nmixup
c
      endif
c
c******************************************************************************
c******************************************************************************
c     write out the value of yplus at the wall if ypluswall > 5.0.
c
      if(ypluswall.gt.5.0)  then
         write(6,*)'ypluswall is being used to obtain the skin friction'
         write(6,*)'yplus at the wall taken as, ', ypluswall
           cfwall = 1/(ypluswall*ypluswall)
      endif
c
c******************************************************************************
c******************************************************************************
c  read in the surface roughnesses in microns if  if_rough  >  0  .
c
      if(if_rough.gt.0) then
c
	write(6,*) ' non-hydraulic smooth surfaces specified'
	write(6,*) ' input surface roughness in microns for all'
	write(6,*) ' surfaces.'
      do 2347 n = 1,nrows
           read(5,*)  rough_h(n),rough_t(n),rough_l(n),rough_u(n)
c
           write(6,*) ' row number ', n
           write(6,*) ' surface roughnesses in microns ',
     &     rough_h(n),rough_t(n),rough_l(n),rough_u(n)
c
c   change to physical roughness in metres . 
           rough_h(n) = rough_h(n)*1.0e-06
           rough_t(n) = rough_t(n)*1.0e-06
           rough_l(n) = rough_l(n)*1.0e-06
           rough_u(n) = rough_u(n)*1.0e-06
c
 2347 continue
c
      else
c
      do 2348 n = 1,nrows
           rough_h(n)  = 0.0
           rough_t(n)  = 0.0
           rough_l(n)  = 0.0
           rough_u(n)  = 0.0
 2348 continue
c
      end if
c
c**********************************************************************************
c**********************************************************************************
c    input the artificial speed of sound if itimst = 5.
c    this should be about half the maximum relative velocity in the flow.
c
      if(itimst.ge.5) then
      vsound    = 150.
      rf_ptru   = 0.01
      rf_vsound = 0.002
      densty    = 1.20
      vs_vmax   = 2.0
      if(itimst.eq.5)  read(5,*,end= 2350)
     &  vsound, rf_ptru, rf_vsound, vs_vmax
      if(itimst.eq.6)  read(5,*,end= 2350)
     &  vsound, rf_ptru, rf_vsound, vs_vmax, densty
 2350 continue
         rf_ptru1    = 1.0 - rf_ptru
         rf_vsound1  = 1.0 - rf_vsound
         write(6,*)
         write(6,*) ' calculation using artificial compressibility '
         write(6,*) ' artificial speed of sound = ', vsound
         write(6,*) ' density relaxation factor = ', rf_ptru
         write(6,*) ' sound speed relaxation factor = ', rf_vsound
         write(6,*) ' ratio of sound speed to maximum speed = ',vs_vmax
         if(itimst.eq.6) write(6,*) 
     &   ' incompressible flow with density = ',densty
         write(6,*)
      end if
c
c******************************************************************************
c******************************************************************************
c     read in the option to use repeating flow conditions
      if(if_repeat.ne.0) then
           read(5,*)    ninmod, rfinbc
           write(6,*) ' repeating stage specified, ninmod, rfinbc = ',
     &                  ninmod, rfinbc 
      end if
c
c******************************************************************************
c******************************************************************************
c
c     read in the stage numbers and sort out the start and end of each stage.
c     first set defaults.
      do  n = 1,nrows
      nstage(n)    = 1 + (n-1)/2
      end do
c
c     now read in the actual stage number for each blade row.
c
      read(5,*,end=3456,err=3456)(nstage(n),n=1,nrows)
 3456 continue
c
      write(6,*) ' read in the stage numbers for each blade row'
      do n= 1,nrows
           write(6,*) ' row number ',n,'is in stage number ',nstage(n)
      end do
c
c  now set the start and end points of each stage
c
      jstg_start(1)   = 1
      nstages         = nstage(nrows)
      nstage(nrows+1) = nstages + 1 
      jstart(nrows+1) = jm +1
c
      do n = 1,nrows
      nstg   = nstage(n)
      nstgp1 = nstage(n+1)
      if(nstg.ne.nstgp1) then
             jstg_end(nstg)     = jmix(n)
             jstg_start(nstgp1) = jstart(n+1)
      end if
      write(6,*) 'row number ',n,'stage number ',nstg,'jstart ',
     &             jstart(n),'jend ',jmix(n)
      end do
c
      do n = 1,nstages
      write(6,*) ' stage number ',n,' jstart= ',jstg_start(n),
     &           ' jend = ',jstg_end(n)
      end do 
c
c**********************************************************************************
c**********************************************************************************
c   input the forcing factor and smoothing factor if doing a throghflow  calculation
c
      if(im.eq.2) then
      q3dforce = 1.0
      sfpbld   = 0.1
      nsfpbld  = 2
      read(5,*) dummy_input
      read(5,*,err=591) q3dforce, sfpbld, nsfpbld
  591 continue
      sfpbld1 = 1.0 - sfpbld
      write(6,*) 'throughflow calculation requested, q3dforce= ',
     & q3dforce, 'smoothing factor=',sfpbld,' no of smoothings=',nsfpbld
      write(6,*)
      end if       
c
c******************************************************************************
c******************************************************************************
c    input the range of yplus values over which the turbulent viscosity will be reduced.
c
      if(ilos.gt.0) then
           yplam    = 5.0
           ypturb   = 25.0
           read(5,*, err= 3458,end = 3458) yplam,ypturb
 3458 continue
           write(6,*) ' turbulent viscosity reduced over the range yplus
     & = ',yplam,' to',ypturb
      end if
c
c******************************************************************************
c******************************************************************************
      write(6,*)  ' subroutine old_readin completed, input data ok'
c
c
      return
      end
c
c*************************************************************************
c*************************************************************************
c*************************************************************************
c*************************************************************************
c
      subroutine loop
c
c      this is the main time stepping loop. it is executed many hundreds of
c      times and uses most of the cpu time. its main subroutine is 'tstep'
c      which also uses much of the time.
c
      include 'commall-open-19.2'
c
      dimension check_flow(jd),blade_flow(jd),rovmsq(kd),ang_inc(kd)
      save start,fini,check_flow,blade_flow
c
c****************************************************************************** 
c    set some reference values
c
      ro_ref = po1(kmid)/(rgas*to1(kmid))
      if(itimst.eq.6) ro_ref = densty
      romin  = amin1(ro(imid,1,kmid),ro(imid,jm,kmid))
      rolim  = 0.25*romin
      rijkm  = 2./(im*jm*km)
c
      do 4001 k=1,km
      plim(k) = 0.9*po1(k) + 0.1*p(1,1,k)
 4001 continue
c
c     set the reference pressure. this is subtracted fron the true pressure to
c     minimise rounding errors.
      preff = 0.5*(p(imid,1,kmid) + p(imid,jm,kmid))
      write(6,*) ' reference pressure in bar = ', preff/1.0e05
c
c****************************************************************************** 
c****************************************************************************** 
c    set more constants
c
c     make the smoothing proportional to the cfl number
      sfxin  = sfxin*cfl/0.5
      sftin  = sftin*cfl/0.5
      rfin1  = 1.0 - rfin
      rfin2  = 2.*rfin
c
      econt    = 1.
      eavg     = 1.0e06
      emax     = 1.0e06
      hbledtot = 0.0
      tcool_min = 1.0e6
c
      ifend    = 0
      nstep    = 0
c
      ibleed = 0
      kbleed = 0
      ncoolb = 0
      ncoolw = 0
c
c******************************************************************************
c     call coolin to set the coolant source terms
c
      if(ifcool.ne.0) then
           call cool_input
           write(6,*) ' called cool_input '
           if(ifcool.eq.1) call coolin_1
           if(ifcool.eq.2) call coolin_2
           write(6,*) ' called coolin '
      endif
c
c******************************************************************************
c     call bleedout to set the bleed flow terms
c
      if(ifbleed.ne.0) then
           call bleedout(ibleed,kbleed)
           write(6,*) ' called bleedout '
      endif
c
c
c******************************************************************************
c      call the timing routines. these are machine specific and may need changing.
c      call clock@(start)
       start = float(mclock())
c
c*************************************************************************
c**********************start of the main time stepping loop ********************
c*************************************************************************
c
      write(6,*)
      write(6,*)'*******************************************************
     &**********'
      write(6,*) ' starting the main time stepping loop '
c
c       return here after every time step
c
 5000 continue
c
c******************************************************************************
c     check if converged or if the maximum number of iterations has been reached.
c     set  "ifend" = 1 to stop the calculation .
c
      ifprint = 0
      nstep  = nstep + 1
      n      = nstep
      if(nstep.eq.nsteps_max) ifend = 1
      if(eavg.lt.conlim.and.econt.lt.0.01)          ifend = 1
      if(abs(emax).lt.(4*conlim).and.econt.lt.0.01) ifend = 1
c
c     check if printout required on this time step
c
      do 5010 i=1,5
 5010 if(nstep.eq.nout(i)) ifprint = 1
c
c    check if a request to stop has been received by "stopit" being set
c    greater than zero.
c
      if(mod(n,10).eq.0) then
           open(unit=12,file='stopit')
           read(12,*) ifstop
           if(ifstop.ge.1) ifend = 1
           close(12)
      endif 
c
c******************************************************************************
c      vary the smoothing and damping over the first  "nchange"  steps.
c      if starting from an initial guess.
c
c     8/4/2017  modify the next section so that damping and smoothing are increased
c     over the first nchange steps even when starting from a restart file.
      if(if_restart.ne.0.and.nstep.gt.100) go to 5555
c
      fchange     = 1.0 -  float(nstep-1)/nchange
c
      if(fchange.gt.0.0)   then
c     delete next line 8/4/2017
c           if(if_restart.ne.0) fchange = 0.0
           sft  = sftin   + 0.02*fchange
           sfx  = sfxin   + 0.02*fchange
           damp = dampin*(1. - 0.75*fchange)
      endif
c
      sft1    = 1.0-sft
      sfth    = 0.5*sft
      sfx1    = 1.-sfx
      sfxh    = 0.5*sfx
      sfx15   = 1.0-1.5*sfx
      sfx14   = 0.25*sfx
      sfxb    = fsmthb*sfx
      sfxb1   = 1.0 - sfxb
      sfxbh   = 0.5*sfxb
      sfxhm1  = 1.-sfxh
c
c    set fmixup to increase the turbulent viscosity over the first nmixup steps
c
      fmixup    = 1.0 + (facmixup-1.)*float(nmixup-nstep)/float(nmixup)
      if(fmixup.lt.1.0) fmixup = 1.0
      fmixlen   = fmixup
      fmixup    = fmixup*fmixup
c
c    jump to here if starting from a restart file
c
 5555 continue
c
c*******************************************************************************
c*******************************************************************************
c
c       calculate vx,vr,vt from rovx,rovr,rovt
c
      do 7450 k=1,km
      do 7450 j=2,jm
           rrnow  = 1/r(j,k)
      do 7450 i=1,im
           recip       = 1.0/ro(i,j,k)
           vx(i,j,k)   = rovx(i,j,k)*recip
           vr(i,j,k)   = rovr(i,j,k)*recip
           rovt(i,j,k) = rorvt(i,j,k)*rrnow
           vt(i,j,k)   = rovt(i,j,k)*recip
 7450 continue
c
c*******************************************************************************
c*******************************************************************************
c
      if(mod(nstep,5).eq.0.or.ifend.eq.1) then
c
      if(ifcool.eq.2) call coolin_2
c
c      every 5 steps sum the relative velocities to find a reference velocity
c      and the maximum relative velocity.
c
      vrms=0.0
      vmax=0.0
      do 5910 k=1,km
      do 5910 j=2,jm
      do 5910 i=1,im
      eke = vr(i,j,k)*vr(i,j,k) + wt(i,j,k)*wt(i,j,k) 
     &       + vx(i,j,k)*vx(i,j,k)
      if(eke.gt.vmax) then
	 vmax  = eke
	 ivmax = i
	 jvmax = j
	 kvmax = k
      endif
      vrms = vrms + eke
 5910 continue
      vmax = sqrt(vmax)
c
      endif
c
c*******************************************************************************
c*******************************************************************************
c
c   calculate new pressures and temperatures.the method depends on the flow model, itimst.
c
      if(itimst.lt.5) then
c
      do 5900 k=1,km
      do 5900 j=2,jmm1
      do 5900 i=1,im
           ronow = ro(i,j,k)
           eke   =.5*(vr(i,j,k)*vr(i,j,k) + vt(i,j,k)*vt(i,j,k) 
     &           + vx(i,j,k)*vx(i,j,k))
c
      if(ifgas.eq.0) then
           tstatic      = (roe(i,j,k)/ronow - eke)/cv
           p(i,j,k)     = ronow*rgas*tstatic
           ho(i,j,k)    = cp*tstatic + eke 
      else
           estat        = roe(i,j,k)/ronow - eke
           tstatic      = tfrome(estat,tref,eref,et1,et2,et3,et4)
c
           if(tstatic.lt.0.1*tref) tstatic = 0.1*tref
           p(i,j,k)     = ronow*rgas*tstatic
           ho(i,j,k)    = estat + rgas*tstatic + eke
      end if
c
           t_static(i,j,k) = tstatic
c
 5900 continue
c
      end if
c
c
      if(itimst.eq.5.or.itimst.eq.6) then
c
c     vary the artificial speed of sound so it very gradually tends to a value
c      = 2 x the maximum relative velocity.
c
      vs_new   = vs_vmax*vmax
c      vsound   = rf_vsound1*vsound + rf_vsound*vs_new
      fracdif  =  vs_new/vsound - 1.0
      dvsound  =  vsound*fracdif*(1.0 + 3*abs(fracdif) )
      if(abs(fracdif).lt.0.1) dvsound = 0.0
      vsound   =  vsound + rf_vsound*dvsound
      dp_dro   =  vsound*vsound
c
      do 5950 k=1,km
      do 5950 j=2,jm
      do 5950 i=1,im
      eke          =.5*(vr(i,j,k)*vr(i,j,k) + vt(i,j,k)*vt(i,j,k) 
     &               + vx(i,j,k)*vx(i,j,k))
      p(i,j,k)     = po_ref  - dp_dro*(ro_ref - rosub(i,j,k))
      tstatic      = (roe(i,j,k)/rosub(i,j,k) - eke)/cv
      ro(i,j,k)    = rf_ptru*p(i,j,k)/rgas/tstatic 
     &             + rf_ptru1*ro(i,j,k)
      if(itimst.eq.6) ro(i,j,k) = densty 
      ho(i,j,k)    = cp*tstatic  + eke 
      t_static(i,j,k) = tstatic
 5950 continue
c     end of itimst = 5  or 6 loop
      end if
c
c*******************************************************************************
c*******************************************************************************
c
c         call loss routines to update body forces every nlos iterations
c
c     if ilos = 10  use the original mixing length model.
      if((ilos.eq.10).and.mod(nstep,nlos).eq.0)  call loss
c
c     if ilos > 100 and < 200  use the new mixing length model.
      if((ilos.ge.100.and.ilos.lt.200).and.mod(nstep,nlos).eq.0)
     &    call new_loss
c
c     if(ilos > 200  use the spalart allmaras model.  
      if(ilos.ge.200.and.mod(nstep,nlos).eq.0)   call spal_loss
c
c*******************************************************************************
c*******************************************************************************
c      call shroudflow to set shroud leakage fluxes if 'ifshroud' ne zero.
c
      if((nstep.eq.1.or.mod(nstep,nlos).eq.0).and.(ifshroud.eq.1))
     &  call shroudflow
c
c*******************************************************************************
c*******************************************************************************
c
c          work out the mass fluxes and store as flowx, flowt, flowr.
c
c       flowr .... mass flow rate through the stream-surface (top) face
c       flowt .... mass flow rate through the bladewise (side) face
c       flowx .... mass flow rate through the quasi-orthogonal (upstream) face
c
      do 5020 k=1,kmm1
      do 5020 j=1,jm
      do 5020 i=1,imm1
      avgrvx = rovx(i,j,k)+rovx(i,j,k+1)+rovx(i+1,j,k)+rovx(i+1,j,k+1)
      avgrvr = rovr(i,j,k)+rovr(i,j,k+1)+rovr(i+1,j,k)+rovr(i+1,j,k+1)
      flowx(i,j,k)   = 0.25*(avgrvx*aqx(i,j,k) + avgrvr*aqr(i,j,k))
      source(i,j,k ) = 0.0
      xflux(i,j,k)   = flowx(i,j,k)
 5020 continue
c
c     calculate and store  rowt(i,j,k) and wt(i,j,k) 
c
      do 5025 k=1,km
      do 5025 j=1,jm
      do 5025 i=1,im
      rowt(i,j,k) = rovt(i,j,k) - ublade(j,k)*ro(i,j,k)
      wt(i,j,k)   = rowt(i,j,k)/ro(i,j,k)
 5025 continue
c
      do 5030 k=1,kmm1
      do 5030 j=2,jm
      do 5030 i=1,im
      avgrvx = rovx(i,j,k)+rovx(i,j-1,k)+rovx(i,j-1,k+1)+rovx(i,j,k+1)
      avgrvt = rowt(i,j,k)+rowt(i,j-1,k)+rowt(i,j-1,k+1)+
     &         rowt(i,j,k+1)
      avgrvr = rovr(i,j,k)+rovr(i,j-1,k)+rovr(i,j-1,k+1)+rovr(i,j,k+1)
      avgro  = ro(i,j,k)+ro(i,j-1,k)+ro(i,j-1,k+1)+ro(i,j,k+1)
      flowt(i,j,k) = (avgrvx*abx(i,j,k) + avgrvt*abt(j,k)
     &             +  avgrvr*abr(i,j,k))*0.25 
      tflux(i,j,k) = flowt(i,j,k) 
 5030 continue
c
c      balance mass fluxes into and out of cusps
c
      do 5050 j=2,jm
      if(ind(j).eq.1) go to 5050
      do 5051 k=1,kmm1
      avgflux        = 0.5*(flowt(1,j,k) + flowt(im,j,k))
      flowt(1,j,k)   = avgflux 
      flowt(im,j,k)  = avgflux
      tflux(1,j,k)   = avgflux
      tflux(im,j,k)  = avgflux
 5051 continue
 5050 continue
c
c**********************************************************************
c    tflow addition
c**********************************************************************
      if(im.ne.2) go to 5655
c
      csq = ga*rgas*to1(kmid) * q3dforce
c
c  set the flow angle within the blade passage using the deviation from the grid angle.
c  the deviation varies linearly with the j index.  
c  gradually change the incidence over "frac_inc" of the j values from the leading edge.
c  add the deviation over "frac_dev" of the j values from the trailing edge.
c
      frac_inc = 0.3333
      frac_dev = 0.5
c
      do 5552 k=1,kmm1
c
      do 5553 j= 2,jm
c
      nr    = nrow(j) 
      jlee  = jle(nr)
      jlem1 = jlee-1
      jlep1 = jlee+1
      jtee  = jte(nr) 
c
      fdevn = float(j-jlee)/float(jtee - jlee)
      fdevn = 1.0 - (1.0 - fdevn)/frac_dev
      if(fdevn.gt.1.00001)    fdevn = 1.0 
      if(fdevn.lt.1.0e-4)     fdevn = 0.0
c
      finc = float(j - jlee)/float(jtee - jlee)/frac_inc
      finc =  1.0 - finc
      if(finc.gt.1.00001)    finc = 1.0
      if(finc.lt.1.0e-4)     finc = 0.0
c 
      fload = 1.0
      if(j.le.jlem1)  fload = 0.0
      if(j.gt.jtee+1) fload = 0.0
c 
c   calculate the flow angle at the leading edge      
c
           rowtavg    = 0.25*(rowt(1,j,k)   + rowt(1,j,k+1)
     &                + rowt(1,j-1,k) + rowt(1,j-1,k+1) )
           xd = dx(j,k)
           rd = dr(j,k)
           sd = ds(j,k)
           rovmk = (xd*(rovx(1,j,k) + rovx(1,j-1,k))
     &           +  rd*(rovr(1,j,k) + rovr(1,j-1,k)))/sd
           xd = dx(j,k+1)
           rd = dr(j,k+1)
           sd = ds(j,k+1)
           rovmkp1  = (xd*(rovx(1,j,k+1) + rovx(1,j-1,k+1))
     &              +  rd*(rovr(1,j,k+1) + rovr(1,j-1,k+1)))/sd
           rovmavg  = 0.25*(rovmk + rovmkp1) 
c
           flo_angl    = atan(rowtavg/rovmavg)
c 
c   set the incidence angle based on the centre line angle at the leading edge
c    
      if(j.eq.jlee)   ang_inc(k)  = flo_angl  - alpha_cent(j,k)
c
      if(j.lt.jlee.or.j.gt.jtee) ang_inc(k)  = 0.0 

c    set the angle to which the flow is forced as a combination of the centre
c    line angle, the incidence angle and the deviation angle.
c
      grid_angl    = alpha_cent(j,k)
      dev_angl     = devn_angl(nr,k)*fdevn  
      angl_inc     = ang_inc(k)*finc
      alpha_forced = alpha_cent(j,k) - dev_angl + angl_inc
c  
c     calculate the change in blade surface pressure due to any flow across
c     the forced angle surface.
c 
      avgtflux1  = 0.5*(tflux(1,j,k) + tflux(2,j,k) ) 
      avgtflux2  = avgtflux1 + rovmavg*abt(j,k)*(tan(alpha_cent(j,k))
     &           - tan(alpha_forced))
      rovarn        = avgtflux2*step(1,j,k)
      pvarn         = csq*(2*rovarn - rovar_m1(j,k) )
      rovar_m1(j,k) = rovarn
      pblade(j,k)   = pvarn + pblade(j,k)
      pblade(j,k)   = pblade(j,k)*fload
      if(k.ge.ktips(nr).and.k.lt.ktipe(nr) ) pblade(j,k) = 0.0
 5553 continue
c
 5552 continue
c
c    smooth the blade surface pressure distribution.
c
      do 5651 n = 1,nsfpbld
      do 5652 k=1,kmm1
      do 5653 j=3,jm-3
      sleft  = smerid(j,k)   - smerid(j-1,k)
      sright = smerid(j+1,k) - smerid(j,k)
      avgpblade(j) = sfpbld1*pblade(j,k)
     &             + sfpbld*(sright*pblade(j-1,k) + sleft*pblade(j+1,k)) 
     &             /(sleft + sright)
 5653 continue
      do 5654 j=3,jm-3
      pblade(j,k) = avgpblade(j)    
 5654 continue
 5652 continue
 5651 continue
c
c
 5655 continue
c
c**********************************************************************
c   end tflow addition
c**********************************************************************
c
c      set zero mass flow through the blade surfaces - where no bleed.
c
      do 5060 j=2,jm
      if(ind(j).ne.1) go to 5060
      nr= nrow(j)      
      do 5061 k=ks1(nr),ks2(nr)
      tflux(1,j,k)  = 0.0
      tflux(im,j,k) = 0.0
      flowt(1,j,k)  = 0.0
      flowt(im,j,k) = 0.0
 5061 continue
 5060 continue
c
c    remove any bleed flows through the blade surfaces
c
      if(ibleed.ne.0) then
      do 5065 j = 2,jm
      if(ind(j).ne.1) go to 5065
      do 5064 k = 1,kmm1
      tflux(1,j,k)  = -blflowi1(j,k)
      tflux(im,j,k) =  blflowim(j,k)
      flowt(1,j,k)  = -blflowi1(j,k)
      flowt(im,j,k) =  blflowim(j,k)
 5064 continue
 5065 continue
      endif
c          
c    now set rflux to the mass flow through the streamwise faces.
c  q3d  
      if(km.ne.2) then
           k1 = 2
           k2 = kmm1
      else
c   end q3d
           k1 = 1
           k2 = 2
      end if
c
      do 5040 k=k1,k2
      do 5040 j=2,jm
      do 5040 i=1,imm1
      avgrvx = rovx(i,j,k)+rovx(i,j-1,k)+rovx(i+1,j-1,k)+rovx(i+1,j,k)
      avgrvr = rovr(i,j,k)+rovr(i,j-1,k)+rovr(i+1,j-1,k)+rovr(i+1,j,k)
      flowr(i,j,k)    = 0.25*(avgrvx*asx(i,j,k)+avgrvr*asr(i,j,k))
      rflux(i,j,k)    = flowr(i,j,k)
 5040 continue
c
c***************************************************************************************
c     set the body force to keep the flow on the stream surface if doing a q3d calculation.
c   q3d
      if(km.ne.2) go to 5042
c
      csq = ga*rgas*to1(1)*q3dforce
      do 5041 j=2,jm
      do 5041 i=1,imm1
      avrflux     = rflux(i,j,1) + rflux(i,j,2)
      rovarn      = avrflux*step(i,j,1)
      pvarn       = csq*rovarn
      xforce(i,j,1) = xforce(i,j,1) + pvarn*(asx(i,j,1)+asx(i,j,2))
      rforce(i,j,1) = rforce(i,j,1) + pvarn*(asr(i,j,1)+asr(i,j,2))
 5041 continue
c
 5042 continue
c      
c  end q3d
c***************************************************************************************
c      set zero mass flow through the hub and casing - where no bleed.
c
      do 5045 j=2,jm
      do 5045 i=1,imm1
      flowr(i,j,1)  = 0.0
      flowr(i,j,km) = 0.0
      rflux(i,j,1)  = 0.0
      rflux(i,j,km) = 0.0
 5045 continue
c
c   q3d
      if(km.eq.2) go to 5081
c   end q3d
c
c     remove any bleed flows through the hub or casing
c
      if(kbleed.ne.0) then
      do 5080 j=1,jm
      do 5080 i=1,imm1
      flowr(i,j,1)     = -blflowk1(i,j)
      flowr(i,j,km)    =  blflowkm(i,j)
      rflux(i,j,1)     =  flowr(i,j,1)
      rflux(i,j,km)    =  flowr(i,j,km)
 5080 continue
      endif
c
 5081 continue
c
c     add coolant mass fluxes through the blade surfaces.
c
      if(ncoolb.ne.0) then
      do 5550 nc=1,ncoolb
      do 5556 j=jcbs(nc)+1,jcbe(nc)
      do 5556 k=kcbs(nc),kcbe(nc)-1
      if(ic(nc).eq.1)  source(1,j,k)     =   - cflowi1(j,k)
      if(ic(nc).eq.im) source(imm1,j,k)  =   - cflowim(j,k)
 5556 continue
 5550 continue
      endif
c
c
c  q3d
      if(km.eq.2) go to 5557
c  end q3d
c
c     add coolant mass flows through the hub and casing
c
      if(ncoolw.ne.0) then
      do 5558 nc=1,ncoolw
      do 5559 j=jcws(nc)+1,jcwe(nc)
      do 5559 i=icws(nc),icwe(nc)-1     
      if(kc(nc).eq.1) source(i,j,1)    = source(i,j,1)    - cflowk1(i,j)
      if(kc(nc).eq.km)source(i,j,kmm1) = source(i,j,kmm1) - cflowkm(i,j)
 5559 continue
 5558 continue
      endif
c
 5557 continue
c
c***********************************************************************
c      set the shroud leakage mass fluxes
c
      if(ifshroud.eq.1) call shroudflux(shroudgr)
c
c*******************************************************************
c*******************************************************************
c       call subroutine tstep to update the density at all points
c
      if(itimst.ge.5) then
           call tstep(rosub,dro,1)
      else
           call tstep(ro,dro,1)
      end if
c
c***********************************************************************
c***********************************************************************
c   form the cell average density if using newloss or spal_loss.
c
      if(ilos.ge.100) then
c
      do 5600 k=1,kmm1
      do 5600 j=2,jm
      do 5600 i=1,imm1
      roavg = 0.125*(ro(i,j,k)+ro(i+1,j,k)+ro(i+1,j,k+1)+ro(i,j,k+1)
     &      + ro(i,j-1,k)+ro(i+1,j-1,k)+ro(i+1,j-1,k+1)+ro(i,j-1,k+1))
      roavg_cell(i,j,k) = roavg
 5600 continue
c
      end if
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c     set the fluxes of turbulent viscosity for the sa model only
c     note that these and the source term are doubled to give an effective
c     larger time step for the turbulent viscosity equation.
c
      if(ilos.ge.200) then
c
      do 6001 k = 1,kmm1
      do 6001 j = 1,jm
      jp1 = j+1
      if(j.eq.jm) jp1 = jm
      if(j.eq.1)  jp1  = 1
      if(indmix(j).eq.1) jp1 = j
      do 6001 i = 1,imm1
      xflux(i,j,k) = flowx(i,j,k)*(trans_kvis(i,j,k)
     &             + trans_kvis(i,jp1,k))
 6001 continue
c
      do 6008 k=1,kmm1
      do 6008 j=2,jm
      do 6008 i=1,imm1
      source(i,j,k) = -2.0*t_source(i,j,k)
 6008 continue
c
c
      do 6002 k=1,kmm1
      do 6002 j=2,jm
      do 6002 i=2,imm1
      tflux(i,j,k) = flowt(i,j,k)*
     &               (trans_kvis(i-1,j,k)+trans_kvis(i,j,k))
 6002 continue
      do 6003 k=1,kmm1
      do 6003 j= 2,jm
      tflux(1,j,k) = flowt(1,j,k)*
     &               (trans_kvis(imm1,j,k)+trans_kvis(1,j,k))
      tflux(im,j,k) = tflux(1,j,k)
 6003 continue
c
c  q3d
      if(km.eq.2) go to 6005
c  end q3d
c
      do 6004 k=2,kmm1
      do 6004 j=2,jm
      do 6004 i=1,imm1
      rflux(i,j,k) = flowr(i,j,k)*
     &               (trans_kvis(i,j,k-1) + trans_kvis(i,j,k))
 6004 continue
c
 6005 continue
c
      do 6006 j=2,jm
      do 6006 i=1,imm1
      rflux(i,j,1) = 0.0
      rflux(i,j,km)= 0.0
 6006 continue
c
c********************************************************************************
c********************************************************************************
c     call tstep to update the turbulent viscosity for the sa model only.
c
      call tstep(trans_dyn_vis, del_dynvis,2)
c
c********************************************************************************
c********************************************************************************
c  
      transvismin = 0.1*visc_lam(imid,2,kmid)
      do 6007 k=1,kmm1
      do 6007 j=2,jm
      do 6007 i=1,imm1
      if(trans_dyn_vis(i,j,k).lt.transvismin) 
     &  trans_dyn_vis(i,j,k) = transvismin
      trans_kvis(i,j,k) = trans_dyn_vis(i,j,k)/roavg_cell(i,j,k)
 6007 continue
c
c     end of setting the turbulent viscosity fluxes for the sa model.
c
      end if
c
c******************************************************************************
c********************************************************************************
c   now work out the fluxes of energy and store as xflux, tflux and rflux.
c
      do 5501 k=1,kmm1
      do 5501 j=1,jm
      do 5501 i=1,imm1
      avgho = ho(i,j,k)+ho(i,j,k+1)+ho(i+1,j,k)+ho(i+1,j,k+1)
      xflux(i,j,k) = flowx(i,j,k)*avgho*0.25
 5501 continue
c
      do 5505 j=2,jm
      do 5505 k=1,kmm1
      do 5505 i=1,imm1
      source(i,j,k) = qsource(i,j,k)
 5505 continue
c
c  tflow addition
      if(im.eq.2) then
      do 5506 k=1,kmm1
      do 5506 j=2,jm
      source(1,j,k) = source(1,j,k) + bforce_q(j,k)
 5506 continue
      end if
c   end tflow addition
c
      do 5502 k=1,km
      do 5502 j=2,jm
      do 5502 i=1,imm1
      avgho = 0.25*(ho(i,j,k)+ho(i,j-1,k)+ho(i+1,j-1,k)+ho(i+1,j,k))
      rflux(i,j,k) = flowr(i,j,k)*avgho
 5502 continue
c
c         
      do 5503 k=1,kmm1
      do 5503 j=2,jm
      do 5503 i=1,im
      avgho = ho(i,j,k)+ho(i,j-1,k)+ho(i,j-1,k+1)+ho(i,j,k+1)
      avgpb = peff(i,j,k)+peff(i,j-1,k)+peff(i,j-1,k+1)+peff(i,j,k+1) 
      tflux(i,j,k) = (flowt(i,j,k)*avgho + wrabt(j,k)*avgpb)*0.25
 5503 continue
c
c
c   calculate the enthalpy bled from the main flow - every 200 steps.
c
      if((mod(nstep,200).eq.0.or.ifprint.eq.1.or.ifend.eq.1)
     &    .and.ifbleed.ne.0)  then
c
      hbledtot = 0.0
      if(kbleed.ne.0) then
      do 5507 j=2,jm
      do 5507 i=1,imm1
      hbledtot  = hbledtot + nblade(j)*(rflux(i,j,km) - rflux(i,j,1))
 5507 continue
      endif
c
      if(ibleed.ne.0) then
      do 5508 j=2,jm
      do 5508 k=1,kmm1
      hbledtot = hbledtot + nblade(j)*(tflux(im,j,k)  - tflux(1,j,k))
 5508 continue
      endif
c
      endif
c
c     add any coolant energy fluxes through the blade surfaces.
c
      if(ncoolb.ne.0) then
      do 5650 nc=1,ncoolb
      do 5656 j = jcbs(nc)+1,jcbe(nc)
      do 5656 k = kcbs(nc),kcbe(nc)-1
      if(ic(nc).eq.1) source(1,j,k)    = source(1,j,k) -  hocwli1(j,k)
      if(ic(nc).eq.im)source(imm1,j,k) = source(imm1,j,k)-hocwlim(j,k)
 5656 continue
 5650 continue
      endif
c
c     add any coolant energy flows through the hub and casing.
c
c   q3d
      if(km.eq.2) go to 5657
c   end q3d
c
      if(ncoolw.ne.0) then
      do 5658 nc=1,ncoolw
      do 5659 j=jcws(nc)+1,jcwe(nc)
      do 5659 i=icws(nc),icwe(nc)-1     
      if(kc(nc).eq.1)  source(i,j,1)    = source(i,j,1) -  hocwlk1(i,j)
      if(kc(nc).eq.km) source(i,j,kmm1) = source(i,j,kmm1)-hocwlkm(i,j)
 5659 continue
 5658 continue
      endif
c
 5657 continue
c
c      balance the fluxes of energy across periodic boundaries
c
      do 900 j=2,jm
      if(ind(j).eq.1) go to 900
      do 901 k=1,kmm1
      tflux(1,j,k)  = 0.5*(tflux(1,j,k)+tflux(im,j,k))
      tflux(im,j,k) = tflux(1,j,k)
  901 continue
  900 continue
c
c      set shroud leakage energy fluxes
c
      if(ifshroud.eq.1) call shroudflux(shroudho)
c
c*******************************************************************************
c*******************************************************************************   
c         call tstep to update the internal energy at all points
c
            call tstep(roe,droe,3)   
c
c*******************************************************************************
c*******************************************************************************
c
c   now start to update the inflow and outflow boundary conditions.
c
c   first the exit boundary conditions
c
c*******************************************************************************
c*******************************************************************************
c
      if(plate_loss.gt.0.001) then
           do 5999 k=1,km
           rovmsq(k) = 0.0
           do 5999 i=1,im
           rovmsq(k) = rovmsq(k) +
     &            rovx(i,jm,k)*vx(i,jm,k) + rovr(i,jm,k)*vr(i,jm,k)
 5999      continue
           rovmsq(k) = 0.5*rovmsq(k)/im
      endif
c
c     vary the fixed hub or casing pressures if using the throttle exit condition.
c
      pthrottle = 0.0
      pdmid     = 0.0
      if(nstep.gt.100.and.throttle_exit.gt.0.001) then
           fflow      = check_flow(jm)/throttle_mas 
           pthrottle  = throttle_pres*fflow*fflow
           pdown_tip  = rfthrotl1*pdown_tip + rfthrotl*pthrottle
           pdown_hub  = rfthrotl1*pdown_hub + rfthrotl*pthrottle
           pdmid      = pthrottle - pd(kmid)
c   8/4/2017.  gradually relax throttle_press to the current value of pthrottle.
           rf_p       = 0.02*rfthrotl
           throttle_pres  = (1.-rf_p)*throttle_pres + rf_p*pthrottle
           if(mod(nstep,50).eq.0) then
           write(6,*) 'throttle set ',throttle_pres,
     &                'pthrottle= ',pthrottle,' flow ratio=', fflow
           end if
      end if
c
      pdiff = 0.0
      if(ipout.eq.-1) pdiff = pdown_tip - pd(km)
c
c************************************************************************
c************************************************************************
c   start a loop over all spanwise points
      do 6000 k=1,km
c
c  form some pitchwise averaged pressures
      pavg2   = 0.0
      pavg3   = 0.0
      pavgm1  = 0.0
      tavgout = 0.0
      sum     = 0.0
      do 6010 i=1,imm1
      pavg2   = pavg2   + 0.5*(p(i,2,k)+p(i+1,2,k))*fp(i)
      pavg3   = pavg3   + 0.5*(p(i,3,k)+p(i+1,3,k))*fp(i)  
      pavgm1  = pavgm1  + 0.5*(p(i,jmm1,k)+p(i+1,jmm1,k))*fp(i)
      tavgout=tavgout   + 
     &       0.5*(t_static(i,jmm1,k) + t_static(i+1,jmm1,k))*fp(i)
      sum     = sum     + 0.5*fp(i)*
     &(rovt(i,jm,k)*vt(i,jm,k)+rovt(i+1,jm,k)*vt(i+1,jm,k))/r(jm,k)
 6010 continue

c
c***********************************************************************
c************************************************************************
c       reset the pressure  "pd(k)" at the downstream (exit) boundary
c       using radial equilibrium and assuming no other radial acceleration
c       if  ipout = 0   or  = -1 .
c       or maintaining constant input values if ipout = 1, or = 3.
c
      if(ipout.le.0) then
           pd(1) = pdown_hub
           if(k.eq.1)     go to 6017
           dp    = (sum + sump)*0.5*(r(jm,k)-r(jm,k-1))
           pd(k) = 0.9*pd(k)+ 0.1*(pd(k-1) + dp)
 6017 continue
           sump  = sum
           if(ipout.eq.-1) pd(k) = pd(k) + pdiff
      end if  
c
c     vary the fixed exit pressure with the exit mass flow rate if throttle_exit gt 0.
c
      if(ipout.ge.1.and.throttle_exit.gt.0.001) then
             pd(k) = pd(k) + rfthrotl*pdmid
      end if
c
c   use a spanwise pressure loss proportional to rovm**2 to make the exit flow
c   more uniform if plate_loss > zero.
c
      if(plate_loss.gt.0.001) pd(k) = pd(k)
     &  + 0.01*plate_loss*(rovmsq(k) - rovmsq(kmid))
c
c    end of setting the exit pressure pd(k)
c********************************************************************************
c********************************************************************************
c
c     new exit boundary condition  option added november 2019.
c     apply characteristics based non_reflecting exit boundary condition.
c     if  if_wave = 1 .
c
      if(if_wave.eq.1) then
c
      do i=1,im
      p(i,jm,k) = p(i,jmm1,k) + fp_xtrap*(p(i,jmm1,k) - p(i,jmm2,k))
      end do
      pavgjm  = 0.0
      do i=1,imm1
      pavgjm  = pavgjm + 0.5*(p(i,jm,k)+p(i+1,jm,k))*fp(i)
      end do
c
      dpwave    = fracwave*(pd(k) - pavgjm)
      csound    = sqrt(ga*rgas*tavgout)
c
      xdiff  = x(jm,k) - x(jmm2,k)
      rdiff  = r(jm,k) - r(jmm2,k)
      dmer   = sqrt(xdiff*xdiff + rdiff*rdiff)
      cosout = xdiff/dmer
      sinout = rdiff/dmer
c
      do 6018 i = 1,im
      roold   =   ro(i,jm,k)
      vmersq  =   vx(i,jm,k)*vx(i,jm,k) + vr(i,jm,k)*vr(i,jm,k)
      eke     =   0.5*(vmersq + vt(i,jm,k)*vt(i,jm,k))
c    
      drowave =   dpwave/(csound*csound)
      dvwave  =  -dpwave/(roold*csound)
      dvxwave =   dvwave * cosout 
      dvtwave =   0.0
      dvrwave =   dvwave * sinout 
      p(i,jm,k)    = p(i,jm,k) + dpwave
      ronew        = roold + drowave
      ro(i,jm,k)   = ronew
      vx(i,jm,k)   = vx(i,jm,k) + dvxwave
      vt(i,jm,k)   = vt(i,jm,k) + dvtwave
      vr(i,jm,k)   = vr(i,jm,k) + dvrwave
      rovx(i,jm,k) = ronew*vx(i,jm,k)
      rovr(i,jm,k) = ronew*vr(i,jm,k)
      rovt(i,jm,k) = ronew*vt(i,jm,k)
      rorvt(i,jm,k) = rovt(i,jm,k)*r(jm,k)
      t_exit       = p(i,jm,k)/(rgas*ronew)
      if(ifgas.eq.0) then
           roe(i,jm,k)  = ronew*(cv*t_exit + eke)
           ho(i,jm,k)   = cp*t_exit + eke
      else
           hstat        = hfromt(t_exit,tref,href,cp1,cp2,cp3)
           ho(i,jm,k)   = hstat + eke
           roe(i,jm,k)  = ronew*(hstat - rgas*t_exit + eke) 
      end if
c
 6018 continue
c
c     end of if_wave = 1 loop.
      end if
c
c*********************************************************************
c*********************************************************************
c 
c     the next section allows the exit pressure to vary along the pitch by
c     extrapolating the variation from upstream. the stagnation enthalpy 
c     is also extrapolated.
c
      if(if_wave.eq.0)  then
c
      do 6029 i=1,im
c
       eke  =.5*(vr(i,jm,k)*vr(i,jm,k) + vt(i,jm,k)*vt(i,jm,k) 
     &      + vx(i,jm,k)*vx(i,jm,k))  
c
      if(itimst.ge.5) then
            p(i,jm,k)    = pd(k)  + fp_xtrap*(p(i,jmm1,k) - pavgm1)     
            rosub(i,jm,k)= ro_ref - (po_ref - p(i,jm,k))/dp_dro
            ho(i,jm,k)   = 2.0*ho(i,jmm1,k) - ho(i,jmm2,k)
            texit        = (ho(i,jm,k) - eke)/cp
            roe(i,jm,k)  = rosub(i,jm,k)*(cv*texit + eke)
	    ro(i,jm,k)   = p(i,jm,k)/rgas/texit 
            if(itimst.eq.6) ro(i,jm,k) = densty   
      else
c
            p(i,jm,k)    = pd(k)  + fp_xtrap*(p(i,jmm1,k) - pavgm1)
            ho(i,jm,k)   = 2.0*ho(i,jmm1,k) - ho(i,jmm2,k)
c
            if(ifgas.eq.0) then
                  texit        = (ho(i,jm,k) - eke)/cp
	          ro(i,jm,k)   = p(i,jm,k)/rgas/texit
                  roe(i,jm,k)  = ro(i,jm,k)*(cv*texit + eke)
c
            else
                  hexit        = ho(i,jm,k) - eke
                  texit        = tfromh(hexit,tref,href,ht1,ht2,ht3,ht4)
	          ro(i,jm,k)   = p(i,jm,k)/rgas/texit
                  roe(i,jm,k)  = ro(i,jm,k)*(hexit - rgas*texit + eke)               
            end if
c
c     end of itimst.ge.5 loop
      end if
c
 6029 continue
c
c     end of if_wave = 0  option
      end if
c
c   end of setting the exit pressure.
c********************************************************************************
c*********************************************************************************           
c     now update the inlet pressure.
c
c     start a loop over all pitchwise grid points
c
      do 6030 i=1,im
c
c         update flow conditions at inlet by extrapolating the pressure to the
c         inlet boundary and taking isentropic flow from the inlet
c         stagnation conditions
c
      if(itimst.ge.5) then
           pnew     = po_ref - (ro_ref - rosub(i,1,k))*dp_dro
           p(i,1,k) = rfin1*p(i,1,k) + rfin*pnew
      else
           ro_stag = po1(k)/rgas/to1(k)
           ro_in   = ro(i,1,k)
           if(ro_in.gt.0.9999*ro_stag) ro_in = 0.9999*ro_stag
c
      if(ifgas.eq.0) then
           pnew    = po1(k) * (ro_in/ro_stag)**ga
      else
           rorat  = ro_in/ro_stag
           trat   = trat_from_rorat(to1(k),rorat,ga1gas,eps0,
     &              eps1,eps2)
           tnew   = to1(k)*trat
           pnew   = ro_in*rgas*tnew
      end if
c
           if(in_press.eq.0) p(i,1,k) = rfin1*p(i,1,k) + rfin*pnew
c
           if(iabs(in_press).eq.1) p(i,1,k) =
     &     rfin1*p(i,1,k) + rfin2*p(i,2,k)- rfin*p(i,3,k)
c
           if(iabs(in_press).eq.3) p(i,1,k) = 
     &     rfin1*p(i,1,k) + rfin2*pavg2 - rfin*pavg3
c
           if(iabs(in_press).eq.4) p(i,1,k) =
     &     rfin1*p(i,1,k) + rfin2*p(imid,2,kmid) - rfin*p(imid,3,kmid)
c
      end if
c
      if(p(i,1,k).gt.plim(k))    p(i,1,k)  = 0.9*p(i,1,k) + 0.1*plim(k)
      if(p(i,1,k).gt.0.99999*po1(k)) p(i,1,k) = po1(k)*0.99999
c
c  next set the velocities and temperatures at inlet
c
      if(itimst.ne.6) then
             if(ifgas.eq.0) then
                  tstatic = to1(k)*(p(i,1,k)*rpo1(k))**fga
                  vabsq   = 2*cp*(to1(k) - tstatic)
             else
                  prat    = p(i,1,k)*rpo1(k)
                  trat    = trat_from_prat(to1(k),prat,fgagas,r_alpha,
     &                      balpha1,balpha2)
                  if(trat.gt.1.0) trat = 0.999999
                  tstatic = to1(k)*trat
                  hstag   = hfromt(to1(k),tref,href,cp1,cp2,cp3) 
                  hstat   = hfromt(tstatic,tref,href,cp1,cp2,cp3)
                  if(hstat.gt.hstag) hstat = 0.99999*hstag 
                  vabsq   = 2.0*(hstag - hstat)
             end if  
      else
             vabsq   = 2.0*(po1(k) - p(i,1,k))/densty
             tstatic = to1(k) - 0.5*vabsq/cp
      end if
c
      t_static(i,1,k) = tstatic
      vabs    = sqrt(vabsq)         
c
      if(in_vtan.eq.0) then
           vmer      = vabs*bscos(k)
           vt(i,1,k) = vabs*bssin(k)
      end if
c
      if(in_vtan.eq.1) then
           vmersq = vabsq - vtin(k)*vtin(k)
           if(vmersq.lt.0.01) vmersq = 0.01
           vmer      = sqrt(vmersq)
           vt(i,1,k) = vtin(k)
      end if
c
      if(in_vtan.eq.2) then
           wr        = ublade(1,k)
           vrelsq    = vabsq - wr*wr*bscos(k)*bscos(k)
           if(vrelsq.lt.0.01) vrelsq = 0.01
           vrel      = sqrt(vrelsq) - wr*bssin(k)
           vmer      = vrel*bscos(k)
           vt(i,1,k) = vrel*bssin(k) + wr
      end if
c
c
      if(in_vr.gt.0) vx(i,1,k) = vmer*brcos(k)
      vxsq = vmer*vmer - vr(i,2,k)*vr(i,2,k)
      if(vxsq.lt.0.0000001*to1(k)*cp) vxsq = 0.0000001*to1(k)*cp
      if(in_vr.le.0) vx(i,1,k) = sqrt(vxsq)
      if(in_vr.gt.0) vr(i,1,k) = vmer*brsin(k)
      if(in_vr.le.0) vr(i,1,k) = vr(i,2,k)
c
      if(itimst.ge.5) then
c      
           if(itimst.eq.5) ro(i,1,k) = rf_ptru*p(i,1,k)/rgas/tstatic 
     &                               + rf_ptru1*ro(i,1,k) 

           if(itimst.eq.6) ro(i,1,k) = densty
c	   
           rovt(i,1,k)  = ro(i,1,k)*vt(i,1,k)
           rovr(i,1,k)  = ro(i,1,k)*vr(i,1,k)
           rorvt(i,1,k) = rovt(i,1,k)*r(1,k)
           rovx(i,1,k)  = ro(i,1,k)*vx(i,1,k)
	   roe(i,1,k)   = rosub(i,1,k)*(cv*tstatic + 0.5*vabsq)
	   ho(i,1,k)    = cp*to1(k)
      else
           ro(i,1,k)    = p(i,1,k)/(rgas*tstatic)
           rovt(i,1,k)  = ro(i,1,k)*vt(i,1,k)
           rovr(i,1,k)  = ro(i,1,k)*vr(i,1,k)
           rorvt(i,1,k) = rovt(i,1,k)*r(1,k)
           rovx(i,1,k)  = ro(i,1,k)*vx(i,1,k)
c
           if(ifgas.eq.0) then
                roe(i,1,k)   = p(i,1,k)/ga1 + ro(i,1,k)*0.5*vabsq
   	        ho(i,1,k)    = cp*to1(k)
           else
                roe(i,1,k)   = ro(i,1,k)*(hstat-rgas*tstatic+0.5*vabsq)
                ho(i,1,k)    = hstag
           end if
c
      end if
c
c******************************************************************************
c******************************************************************************
c  update the turbulent viscosity at inlet if using the sa model.
c
      if(ilos.ge.200) then
           if(fsturb(1).lt.1.0) then 
           transvisin  = visc_lam(i,2,k)*4.35*fsturb(1)**0.25
           else
           transvisin  = visc_lam(i,2,k)*(3.5 + fsturb(1)*0.85)
           end if
           trans_dyn_vis(i,1,k) = transvisin
           trans_dyn_vis(i,2,k) = trans_dyn_vis(i,1,k)
           trans_kvis(i,1,k)    = trans_dyn_vis(i,1,k)/ro(imid,1,kmid)
           trans_kvis(i,2,k)    = trans_kvis(i,1,k)
      end if
c
c     end of setting the inlet boundary conditions.
c
c   end of the " do i "  loop
 6030 continue
c
c     end of the  " do  k " loop
 6000 continue
c
c******************************************************************************
c******************************************************************************
c    all inlet and exit boundary conditions are now updated
c
c******************************************************************************
c******************************************************************************
c       call setflo to adjust the rov's to try to achieve a specified mass
c       flow rate if in_flow=3. or to relax to an average flow if in_flow=2
c       in_flow=2  often gives improved convergence.
c
      if(in_flow.ne.0) call setflo
c
c    form peff to minimise rounding errors.
c    also allow downwinding of pressure to smear shocks using  fp_down .
c
      do 4321 k=1,km
      do 4321 j=1,jm
      jp1 = j+1
      if(j.eq.jm) jp1 = jm
      do 4321 i=1,im
      peff(i,j,k) = p(i,j,k) - preff + f_pdown*(p(i,jp1,k) - p(i,j,k))
 4321 continue
c
c*******************************************************************************
c   tflow addition
c*******************************************************************************
      if(im.eq.2) then
c
c     resolve the blade force perpendicular to the streamline so it does not
c     generate any loss.
c
      do 6050 k=1,kmm1
      kp1 = k+1
      do 6050 j=2,jm
      jm1 = j-1
      p_diff = 0.5*pblade(j,k)
      peff(1,j,k)    = peff(1,j,k)     - p_diff
      peff(1,j,kp1)  = peff(1,j,kp1)   - p_diff
      peff(1,jm1,k)  = peff(1,jm1,k)   - p_diff
      peff(1,jm1,kp1)= peff(1,jm1,kp1) - p_diff
      peff(2,j,k)    = peff(2,j,k)     + p_diff
      peff(2,j,kp1)  = peff(2,j,kp1)   + p_diff
      peff(2,jm1,k)  = peff(2,jm1,k)   + p_diff
      peff(2,jm1,kp1)= peff(2,jm1,kp1) + p_diff
 6050 continue
c
      do 6051 j=2,jm

      do 6052 k=1,kmm1
      avgp1 =(peff(1,j,k)+peff(1,j,k+1)+peff(1,j-1,k)+peff(1,j-1,k+1))/4
      avgp2 =(peff(2,j,k)+peff(2,j,k+1)+peff(2,j-1,k)+peff(2,j-1,k+1))/4
      abxavg = (abx(1,j,k) + abx(2,j,k))/2
      abravg = (abr(1,j,k) + abr(2,j,k))/2
      xfor   = (avgp1 - avgp2)*abxavg
      rfor   = (avgp1 - avgp2)*abravg
      tfor   = (avgp1 - avgp2)*abt(j,k)
c      
      vxavg = vx(1,j,k)+vx(1,j,k+1)+vx(1,j-1,k)+vx(1,j-1,k+1)
      vravg = vr(1,j,k)+vr(1,j,k+1)+vr(1,j-1,k)+vr(1,j-1,k+1)
      wtavg = wt(1,j,k)+wt(1,j,k+1)+wt(1,j-1,k)+wt(1,j-1,k+1)
      wavg = sqrt(vxavg*vxavg + vravg*vravg + wtavg*wtavg)
      sfor = (xfor*vxavg + rfor*vravg + tfor*wtavg)/wavg
c        
      if(ind(j).eq.0) go to 6052
c
      bforce_x(j,k) = sfor*vxavg/wavg
      bforce_t(j,k) = sfor*wtavg/wavg*ravg_cell(j,k)
      bforce_r(j,k) = sfor*vravg/wavg
      bforce_q(j,k) = bforce_t(j,k)*wrad(j)
c
 6052 continue
 6051 continue
c      
      end if
c*******************************************************************************
c   end of tflow addition
c*******************************************************************************
c
c*******************************************************************************
c       calculate fluxes for the axial-momentum equation.
c*******************************************************************************
c******************************************************************************
      do 6100 k=1,kmm1
      do 6100 j=1,jm
      do 6100 i=1,imm1
      source(i,j,k) = xforce(i,j,k)
      avgp  = peff(i,j,k)+peff(i,j,k+1)+peff(i+1,j,k+1)+peff(i+1,j,k)
      avgvx = vx(i,j,k)+vx(i,j,k+1)+vx(i+1,j,k+1)+vx(i+1,j,k)
      xflux(i,j,k) = 0.25*(avgp*aqx(i,j,k) + flowx(i,j,k)*avgvx)
 6100 continue
c
c  tflow addition
      if(im.eq.2) then
      do 6101 k=1,kmm1
      do 6101 j=2,jm
      source (1,j,k) = source(1,j,k) + bforce_x(j,k)
 6101 continue
      end if
c   end tflow addition
c
      do 6110 k=1,km
      do 6110 j=2,jm
      do 6110 i=1,imm1
      avgp  = peff(i,j,k)+peff(i,j-1,k)+peff(i+1,j-1,k)+peff(i+1,j,k)
      avgvx = vx(i,j,k)+vx(i,j-1,k)+vx(i+1,j-1,k)+vx(i+1,j,k)
      rflux(i,j,k) = 0.25*(avgp*asx(i,j,k) + flowr(i,j,k)*avgvx)
 6110 continue
c
c
      do 6120 k=1,kmm1
      do 6120 j=2,jm
      do 6120 i=1,im
      avgvx  = vx(i,j,k)+vx(i,j-1,k)+vx(i,j-1,k+1)+vx(i,j,k+1)
      avgp   = peff(i,j,k)+peff(i,j-1,k)+peff(i,j-1,k+1)+peff(i,j,k+1)
      tflux(i,j,k) = 0.25*(flowt(i,j,k)*avgvx + abx(i,j,k)*avgp)
 6120 continue
c
c     add coolant axial momentum fluxes through blades
c
      if(ncoolb.ne.0) then
      do 6650 nc=1,ncoolb
      do 6656 j = jcbs(nc)+1,jcbe(nc)
      do 6656 k = kcbs(nc),kcbe(nc)-1
      if(ic(nc).eq.1) source(1,j,k)   = source(1,j,k)    - vxcwli1(j,k)
      if(ic(nc).eq.im)source(imm1,j,k)= source(imm1,j,k) - vxcwlim(j,k)
 6656 continue
 6650 continue
      endif
c
c  q3d
      if(km.eq.2) go to 6657
c  end q3d
c     add coolant axial momentum flows through the hub and casing
c
      if(ncoolw.ne.0) then
      do 6658 nc=1,ncoolw
      do 6659 j=jcws(nc)+1,jcwe(nc)
      do 6659 i=icws(nc),icwe(nc)-1
      if(kc(nc).eq.1)  source(i,j,1)   = source(i,j,1)    - vxcwlk1(i,j)
      if(kc(nc).eq.km) source(i,j,kmm1)= source(i,j,kmm1) - vxcwlkm(i,j)
 6659 continue
 6658 continue
      endif
c
 6657 continue
c
c          balance fluxes on periodic boundaries for rovx
c
      do 6026 j=2,jm
      if(ind(j).eq.1) go to 6026
      do 6027 k=1,kmm1
      avgpcusp =  peff(1,j,k)     + peff(1,j-1,k) + peff(1,j,k+1)
     &          + peff(1,j-1,k+1) + peff(im,j,k)  + peff(im,j-1,k)
     &          + peff(im,j,k+1)  + peff(im,j-1,k+1)
      avgvxcusp = vx(1,j,k)       + vx(1,j-1,k) + vx(1,j,k+1)
     &          + vx(1,j-1,k+1)   + vx(im,j,k)  + vx(im,j-1,k)
     &          + vx(im,j,k+1)    + vx(im,j-1,k+1)      
      tflux(1,j,k) =0.125*(flowt(1,j,k)*avgvxcusp +avgpcusp*abx(1,j,k))
      tflux(im,j,k)=0.125*(flowt(im,j,k)*avgvxcusp+avgpcusp*abx(im,j,k))
 6027 continue
 6026 continue
c
c      set shroud leakage axial momentum fluxes
c
      if(ifshroud.eq.1)  call shroudflux(shroudvx)
c
c******************************************************************************
c******************************************************************************
c      call subroutine tstep to update the axial-momentum
c
      call tstep(rovx,drovx,4)
c
c******************************************************************************
c******************************************************************************
c       calculate fluxes for moment of momentum equation.
c******************************************************************************
c
      do 6200 k=1,kmm1
      do 6200 j=1,jm
      avgr  =  0.125*(r(j,k)+r(j,k+1))
      do 6200 i=1,imm1
      avgvt = vt(i,j,k)+vt(i,j,k+1)+vt(i+1,j,k)+vt(i+1,j,k+1)
      source(i,j,k) = tforce(i,j,k)
      xflux(i,j,k)  = flowx(i,j,k)*avgr*avgvt
 6200 continue
c
c  tflow addition
      if(im.eq.2) then
      do 6201 k=1,kmm1
      do 6201 j=2,jm
      source(1,j,k) = source(1,j,k) + bforce_t(j,k)
 6201 continue
      end if
c   end tflow addition
c
      do 6210 k=1,kmm1
      do 6210 j=2,jm
      avgr  = 0.25*ravg_cell(j,k)
      do 6210 i=1,im
      avgvt = vt(i,j,k)+vt(i,j-1,k)+vt(i,j-1,k+1)+vt(i,j,k+1)
      avgp  = peff(i,j,k)+peff(i,j-1,k)+peff(i,j-1,k+1)+peff(i,j,k+1)
      tflux(i,j,k) = (flowt(i,j,k)*avgvt  +  abt(j,k)*avgp)*avgr
 6210 continue
c
      do 6230 k=1,km
      do 6230 j=2,jm
      avgr=(r(j,k)+r(j-1,k))*0.125
      do 6230 i=1,imm1
      avgvt = vt(i,j,k)+vt(i,j-1,k)+vt(i+1,j-1,k)+vt(i+1,j,k)
      rflux(i,j,k) = flowr(i,j,k)*avgvt*avgr
 6230 continue
c
c      balance fluxes of angular momentum across periodic boundaries
c
      do 902 j=2,jm
      if(ind(j).eq.1) go to 902
      do 903 k=1,kmm1
      tflux(1,j,k)  = 0.5*(tflux(1,j,k) + tflux(im,j,k))
      tflux(im,j,k) = tflux(1,j,k)
  903 continue
  902 continue
c
c     add coolant tangential momentum fluxes through blades
c
      if(ncoolb.ne.0) then
      do 7750 nc = 1,ncoolb
      do 7756 j  = jcbs(nc)+1,jcbe(nc)
      do 7756 k  = kcbs(nc),kcbe(nc)-1
      if(ic(nc).eq.1)  source(1,j,k)  = source(1,j,k)    - rvtcwli1(j,k)
      if(ic(nc).eq.im) source(imm1,j,k)=source(imm1,j,k) - rvtcwlim(j,k)
 7756 continue
 7750 continue
      endif
c
c  q3d
      if(km.eq.2) go to 7757
c  end q3d
c
c     add coolant tangential momentum flows through hub and casing
c
      if(ncoolw.ne.0) then
      do 7758 nc = 1,ncoolw
      do 7759 j  = jcws(nc)+1,jcwe(nc)
      do 7759 i=icws(nc),icwe(nc)-1
      if(kc(nc).eq.1)  source(i,j,1)  = source(i,j,1)    - rvtcwlk1(i,j)
      if(kc(nc).eq.km) source(i,j,kmm1)=source(i,j,kmm1) - rvtcwlkm(i,j)
 7759 continue
 7758 continue
      endif
c
 7757 continue
c
c     set shroud leakage angular momentum fluxes
c
      if(ifshroud.eq.1) call shroudflux(shroudrvt)
c
c*******************************************************************************
c******************************************************************************
c       call subroutine tstep to update the moment of momentum
c
      call tstep(rorvt,drorvt,5)
c
c*******************************************************************************
c******************************************************************************
c       calculate fluxes for the radial momentum equation
c******************************************************************************
c
      do 6300 k=1,kmm1
      do 6300 j=1,jm
      do 6300 i=1,imm1
      avgvr = vr(i,j,k)+vr(i,j,k+1)+vr(i+1,j,k)+vr(i+1,j,k+1)
      avgp  = peff(i,j,k)+peff(i,j,k+1)+peff(i+1,j,k+1)+peff(i+1,j,k)
      xflux(i,j,k) = 0.25*(flowx(i,j,k)*avgvr + aqr(i,j,k)*avgp)
 6300 continue

      do 6310 k=1,kmm1
      do 6310 j=2,jm
      do 6310 i=1,im
      avgvr = vr(i,j,k)+vr(i,j-1,k)+vr(i,j-1,k+1)+vr(i,j,k+1)
      avgp  = peff(i,j,k)+peff(i,j-1,k)+peff(i,j-1,k+1)+peff(i,j,k+1)
      tflux(i,j,k) =.25*(flowt(i,j,k)*avgvr + avgp*abr(i,j,k))
 6310 continue
c
      do 6320 k=1,km
      do 6320 j=2,jm
      do 6320 i=1,imm1
      avgvr = vr(i,j,k)+vr(i,j-1,k)+vr(i+1,j-1,k)+vr(i+1,j,k)
      avgp  = peff(i,j,k)+peff(i,j-1,k)+peff(i+1,j-1,k)+peff(i+1,j,k)
      rflux(i,j,k) = 0.25*(flowr(i,j,k)*avgvr + asr(i,j,k)*avgp)
 6320 continue
c
c     calculate the source term. this is due to the centrifugal force and any 
c     imbalance in the radial projected areas.
c
      do 6330 k=1,kmm1
      do 6330 j=2,jm
      do 6330 i=1,imm1
      avgrvt = rovt(i,j,k)+rovt(i,j,k+1)+rovt(i,j-1,k)+rovt(i,j-1,k+1)
     & + rovt(i+1,j,k)+rovt(i+1,j,k+1)+rovt(i+1,j-1,k)+rovt(i+1,j-1,k+1)
      avgvt = 0.125*(vt(i,j,k)+vt(i,j-1,k)+vt(i,j-1,k+1)+vt(i,j,k+1)
     & + vt(i+1,j,k)+vt(i+1,j-1,k)+vt(i+1,j-1,k+1)+vt(i+1,j,k+1))
      avgp = peff(i,j,k)+ peff(i,j-1,k)+ peff(i,j-1,k+1)+ peff(i,j,k+1)
     & + peff(i+1,j,k)+peff(i+1,j-1,k)+peff(i+1,j-1,k+1)+peff(i+1,j,k+1)
      source(i,j,k) = rforce(i,j,k) - (avgrvt*avgvt + avgp)*volor(i,j,k)
 6330 continue
c
c  tflow addition
      if(im.eq.2) then
      do 6331 k=1,kmm1
      do 6331 j=2,jm
      source(1,j,k) = source(1,j,k) + bforce_r(j,k)
 6331 continue
      end if
c   end tflow addition
c
c     add coolant radial momentum fluxes through the blade surfaces.
c
      if(ncoolb.ne.0) then
      do 9650 nc=1,ncoolb
      do 9656 j = jcbs(nc)+1,jcbe(nc)
      do 9656 k = kcbs(nc),kcbe(nc)-1
      if(ic(nc).eq.1)  source(1,j,k)   = source(1,j,k)    - vrcwli1(j,k)
      if(ic(nc).eq.im) source(imm1,j,k)= source(imm1,j,k) - vrcwlim(j,k)
 9656 continue
 9650 continue
      endif
c
c     add coolant radial momentum flows through the hub and casing
c
c  q3d
      if(km.eq.2) go to 8657
c  end q3d
c
      if(ncoolw.ne.0) then
      do 8658 nc=1,ncoolw
      do 8659 j=jcws(nc)+1,jcwe(nc)
      do 8659 i=icws(nc),icwe(nc)-1
      if(kc(nc).eq.1)  source(i,j,1)   = source(i,j,1)    - vrcwlk1(i,j)
      if(kc(nc).eq.km) source(i,j,kmm1)= source(i,j,kmm1) - vrcwlkm(i,j)
 8659 continue
 8658 continue
      endif
c
 8657 continue
c
c          balance fluxes on periodic boundaries for rovr
c
      do 6260 j=2,jm
      if(ind(j).eq.1) go to 6260
      do 6261 k=1,kmm1
      avgp1  =  0.25*(peff(1,j,k)+peff(1,j-1,k)+peff(1,j,k+1)
     &              + peff(1,j-1,k+1))
      avgpm  =  0.25*(peff(im,j,k)+peff(im,j-1,k)+peff(im,j,k+1)
     &              + peff(im,j-1,k+1))
      t1 = tflux(1,j,k) - avgp1*abr(1,j,k)
      t2 = tflux(im,j,k)- avgpm*abr(im,j,k)
      tflux(1,j,k)   = 0.5*((t1+t2) + (avgp1+avgpm)*abr(1,j,k))
      tflux(im,j,k)  = 0.5*((t1+t2) + (avgp1+avgpm)*abr(im,j,k))
 6261 continue
 6260 continue
c
c      set shroud leakage radial momentum fluxes
c
      if(ifshroud.eq.1) call shroudflux(shroudvr)
c
c******************************************************************************
c******************************************************************************
c       call subroutine tstep to update the radial-momentum
c
      call tstep(rovr,drovr,6)
c
c******************************************************************************
c******************************************************************************
c      all the conservation equations have now been updated
c******************************************************************************
c******************************************************************************
c
c      check for any very steep streamwise density gradients or any very low densities..
c
      do 7200 k=1,km
      do 7200 j=3,jm-3
      do 7200 i=1,im
      roavg = 0.25*(ro(i,j-2,k)+ro(i,j-1,k)+ro(i,j+1,k)+ro(i,j+2,k))
      if(ro(i,j,k).lt.0.6*roavg) ro(i,j,k) = 0.5*(ro(i,j,k)+0.6*roavg)
      if(ro(i,j,k).gt.1.5*roavg) ro(i,j,k) = 0.5*(ro(i,j,k)+1.5*roavg)
      if(ro(i,j,k).lt.rolim) ro(i,j,k) = 0.5*(ro(i,j,k) + rolim)
 7200 continue
c
c******************************************************************************
c  check for any very high mach numbers - limit = machlim
c
      gm1      = ga - 1.0
      tratio   = 1.0 + machlim*machlim*0.5*gm1
      rfgm1    = 1.0/gm1
      roratio  = tratio**rfgm1
      do 7250 nr = 1,nrows
      jref   = jle(nr)
      treff  = t_static(imid,jref,kmid)
      wrefsq = vx(imid,jref,kmid)*vx(imid,jref,kmid) +
     &         vr(imid,jref,kmid)*vr(imid,jref,kmid) +
     &         wt(imid,jref,kmid)*wt(imid,jref,kmid)
      torel_ref   = treff + 0.5*wrefsq/cp
      tlimit      = torel_ref/tratio
      if(tlimit.gt.tcool_min) tlimit = tcool_min
      vlimit(nr)  = sqrt(2*cp*(torel_ref - tlimit) )
      rostag_rel  = ro(imid,jref,kmid) *(torel_ref/treff)**rfgm1
      rolimit(nr) = rostag_rel/roratio
 7250 continue
c
c    check and apply limits to both the relative mach number and relative mass flux .
c
      do 7100 k=1,km
      do 7100 j=1,jm
      nr       = nrow(j)
      v_limit  = vlimit(nr)
      ro_limit = rolimit(nr)
      do 7100 i=1,im
      if(ro(i,j,k).lt.ro_limit) ro(i,j,k) = ro_limit
      ronow       = ro(i,j,k)
      wrelsq      = vx(i,j,k)*vx(i,j,k) + vr(i,j,k)*vr(i,j,k)
     &            + wt(i,j,k)*wt(i,j,k)
      wrel        = sqrt(wrelsq)
      facsafe     = v_limit/wrel
      if(facsafe.lt.1.0) then
           rovx(i,j,k) = ronow*vx(i,j,k) * facsafe
           rovr(i,j,k) = ronow*vr(i,j,k) * facsafe
           rovtrel     = ronow*wt(i,j,k) * facsafe
           rowt(i,j,k) = rovtrel
           rovt(i,j,k) = rovtrel + ronow*ublade(j,k) 
           rorvt(i,j,k)= rovt(i,j,k)*r(j,k)
      end if
c
 7100 continue
c
c******************************************************************************
c******************************************************************************
c     call subroutine new_mixplan to transfer the flow across the mixing planes
c
      if(im.gt.2) then
           if(ifmix.gt.0.and.rfmix.gt.1.0e-3) call new_mixplan
      end if
c
c******************************************************************************
c******************************************************************************
c      force the trailing edge separation points if  "if_cusp" = 2
c******************************************************************************
c******************************************************************************
c
           do 7355 nr = 1,nrows
c
       if(if_cusp(nr).eq.2) then
                jtedge   = jte(nr)
                jsep_i1  = jtedge - nup_i1(nr)
                jsep_im  = jtedge - nup_im(nr)
                jthik    = max0(jsep_i1,jsep_im) 
                fdrag    = sep_drag(nr)
                nwke     = n_wake(nr)               
c
           do 7350 k=1,km
           thik_lim = sep_thik(nr)*rt_thick(jthik,k)
           ssslope = (rtheta(1,jsep_i1+1,k) - rtheta(1,jsep_i1-1,k))
     &              /(x(jsep_i1+1,k) - x(jsep_i1-1,k))
           psslope = (rtheta(im,jsep_im+1,k) - rtheta(im,jsep_im-1,k))
     &              /(x(jsep_im+1,k) - x(jsep_im-1,k))
c
           do 7351 j = jsep_i1,jtedge+nwke
           rtextrap = rtheta(1,jsep_i1,k)+ ssslope*(x(j,k)-x(jsep_i1,k))
           do 7353 i = 1,6
           if((rtextrap-rtheta(i,j,k)).gt.thik_lim) then
                rovx(i,j,k)  = fdrag*rovx(i,j,k)
                rovr(i,j,k)  = fdrag*rovr(i,j,k)
                ro_vblade  = ro(i,j,k)*ublade(j,k)
                ro_wt      = rovt(i,j,k) - ro_vblade
                ro_vt      = fdrag*ro_wt + ro_vblade                 
                rorvt(i,j,k) = ro_vt*r(j,k)
           end if
 7353      continue
 7351      continue
c
           do 7352 j = jsep_im,jtedge+nwke
           rtextrap = rtheta(im,jsep_im,k)+psslope*(x(j,k)-x(jsep_im,k))
           do 7354 i = im-6,im
           if((rtheta(i,j,k)-rtextrap).gt.thik_lim) then
                rovx(i,j,k)  = fdrag*rovx(i,j,k)
                rovr(i,j,k)  = fdrag*rovr(i,j,k)
                ro_vblade  = ro(i,j,k)*ublade(j,k)
                ro_wt      = rovt(i,j,k) - ro_vblade
                ro_vt      = fdrag*ro_wt + ro_vblade                 
                rorvt(i,j,k) = ro_vt*r(j,k)
           end if
 7354      continue
 7352      continue
c
 7350      continue
c
       end if
c
 7355      continue
c
c******************************************************************************
c*********************  end of the main time stepping loop ********************
c
c     go to 8000 to check the convergence every 5 time steps.
c
      if(mod(nstep,5).eq.0.or.(ifend.eq.1)) go to 8000
c
c******************************************************************************
c******************************************************************************
c      return to the start of the main loop for the next timestep
c
      go to 5000
c
c****************************************************************************
c      the remainder of this sub program is only executed every 5 steps
c****************************************************************************
c
 8000 continue
c
c       every five time steps check convergence and print out summary to units 4 and 6
c
c        evaluate the inlet and local mass flows.
c
c        blade_flow  is the total mass flow through the blade passages,
c        including and cooling, leakage  or bleed flows.
c
c        check_flow  is the local mass flow excluding cooling flows, leakage and bleed flows.
c        this should be conserved and is used as a check for global continuity.  
c
      do 5675 j=1,jm
      sumas = 0
      nb = nblade(j)
      do 5665 i=1,imm1
      do 5665 k=1,kmm1
      sumas = sumas  - flowx(i,j,k)
 5665 continue
      blade_flow(j) =  sumas*nb
      check_flow(j)  = (sumas + shrdflow(j))*nb+ sumbleed(j)- sumcwl(j)
 5675 continue
      econt    = 0.0
      do 5677 j=1,jm
      ratio = check_flow(j)/check_flow(1)
      emass = abs(1.-ratio)
      if(emass.gt.econt) then
          econt = emass
          jcont = j
      end if
 5677 continue
      flowrat =  blade_flow(1)
c
c******************************************************************************
c******************************************************************************
c     calculate the maximum percentage change in meridional velocity and save
c     it as store(i,j,k) for the convergence check.
c     this is only done every 5 iterations.
c
      emax=0.0
      eavg=0.0
      vref = sqrt(vrms*rijkm)
      rvef = 100./vref
      do 8500 k=1,km
c    jdd changed next line to jmm1 for if_wave =1 option.
      do 8510 j=2,jmm1
      xd = dx(j,k)
      rd = dr(j,k)
      sd = ds(j,k)
      do 8520 i=1,im
      vm_start = (xd*vx(i,j,k)  + rd*vr(i,j,k))/sd
      vm_end   = (xd*rovx(i,j,k)
     &         +  rd*rovr(i,j,k))/sd/ro(i,j,k)
      dvmer  = (vm_end - vm_start)*rvef
      store(i,j,k) = dvmer
      eavg   = eavg + abs(dvmer)
      if(abs(dvmer).gt.abs(emax)) then
           emax = dvmer
           imax = i
           jmax = j
           kmax = k
      endif
 8520 continue
 8510 continue
 8500 continue
c
      eavg=eavg/(im*jmm1*km)
c
c******************************************************************************
c******************************************************************************
c    calculate the pitchwise average values at exit for use in setting up a repeating stage condition.
c
      if(if_repeat.gt.0)     then
      if(mod(n,ninmod).eq.0) then
c
      call mix_bconds(1)
c
      call newbconds                 
c
      end if
      end if
c
c******************************************************************************
c******************************************************************************
c    call output to print out main printed output if requested, or if converged.
c    call the loss routines to update the viscous forces before printing out.
c
      if(ifprint.eq.1.or.ifend.eq.1) then
      if(ilos.eq.10)                  call loss
      if(ilos.ge.100.and.ilos.lt.200) call new_loss
      if(ilos.ge.200)                 call spal_loss
      call output
      end if
c
c******************************************************************************
c******************************************************************************
c      call stepup every to update the timestep if itimst=3
c
      if(itimst.ge.3) call stepup(0.25)
c
c******************************************************************************
c******************************************************************************
c
c      call eficool to print out mass averaged quantities
c      and machine efficiency and shroud leakage flows every 200 steps.
c
      if(mod(nstep,200).eq.0.or.ifprint.eq.1.or.ifend.eq.1) then
      call eficool(hbledtot)
c
c******************************************************************************
c******************************************************************************
c
c     write out the shroud leakage flows and frictional work. also every 200 steps.
c
      if(ifshroud.eq.1) then
c
      do 9099 nr = 1,nrows
      if(ktips(nr).lt.0) then
           write(6,*)
           write(6,*) 'row no ',nr,'percentage shroud leakage flow = ',
     &                 100.*sleak(nr)/blade_flow(1),' % '
	   write(6,*) 'work done on the shroud and hub and casing by the
     & leakage flow = ',swork(nr),'watts'
      endif
 9099 continue
c
       write(6,*)
       write(6,*) ' total frictional work done on all the shrouds = ',
     &              sworktot,' watts'
c
c    end of output for shroud leakage
      endif
c
c   end of output every 200 steps or when finished or converged.
      endif
c
c******************************************************************************
c******************************************************************************
c      write out a short output summary every 5 steps
c
      if(nstep.eq.5.or.mod(nstep,50).eq.0) write(6,9002)itimst,cfl,
     &    damp,sft,sfx,fextrap,f1,f2eff,f3,nrsmth
 9002 format(/,' istep=',i2,' cfl=',f5.2,' damp= ',f5.2,' sft,sfx= ',
     &2f6.3,' fextrap= ',f5.3,' f1=',f5.2,' f2eff=',f5.2,' f3=',f5.2,
     &' nrsmth=',i3,/)
c
      if(itimst.ge.5.and.mod(nstep,50).eq.0) write(6,*)
     &   ' artificial speed of sound   = ', vsound
c
c******************************************************************************
c******************************************************************************
c     the following timing call will differ for different computers.
c
      if(mod(nstep,50).eq.0) then
c      call clock@(fini)
      fini = float(mclock())
      runtime   = (fini - start)/50*1.0e-06
      pointime  = runtime/(im*jm*km)
      write(6,*) ' cpu time per step=', runtime, 'seconds.',
     & ' cpu time per point per step=', pointime,'seconds.'
      write(6,*)
c      call clock@(start)
      start = float(mclock())
      endif
c
c******************************************************************************
c     end of timing call
c******************************************************************************
c******************************************************************************
c     write out a summary to the screen every 5 steps
c
      if(nstep.eq.5.or.mod(nstep,50).eq.0) write(6,9001)
 9001 format(' step   emax at i   j   k  eavg    econt at j=    vref   v 
     &max  at i   j   k     in flow     out flow    rat flow')
c
      ratflow = blade_flow(jm)/flowrat
      write(6,9000) nstep,emax,imax,jmax,kmax,eavg,econt,jcont,vref,
     &              vmax,ivmax,jvmax,kvmax,flowrat,blade_flow(jm),
     &              ratflow
 9000 format(i5,f8.4,3i4,2f8.4,i6,2f8.2,3i4,3f12.4)
c
c******************************************************************************
c******************************************************************************
c   write a convergence history to unit 4 output summary every 5 steps
c
c
      write(4,9009) emax,eavg,econt,flowrat,nstep,imax,jmax,kmax
 9009 format('emax,eavg,econt',3e12.5,' flow',e12.5,'step',i5,'at',3i5)
c
c******************************************************************************
c******************************************************************************
c    call mix_bconds to write out the pitchwise average values at the mixing planes.
c
      if(ifend.eq.1) then
      call mix_bconds(1)
      end if
c
c*******************************************************************************
c*********stop if converged or maximum iterations reached**********
c
      if(ifend.eq.1) stop
c
c********************** return to start of main loop ************
c**************  if not yet converged or at last iteration ******
c
      go to 5000
c
c******************************************************************************
c******************************************************************************
c
c     end of the subroutine loop
c
      end
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c
      subroutine tstep(d,diff,ncall)
c       ====================
c
c       this routine updates the given variable every
c       timestep (i.e. ro rovr rovt rovx roe)
c       and also performs the pitchwise smoothing
c       most of the computational time is used by this subroutine
c
c       ====================
      include  'commall-open-19.2'
c
      dimension d(id,jd,kd), diff(id,jd,kd),avg(kd),
     &          b1chg(ig1,jg1,kg1),b2chg(ig2,jg2,kg2),sbchg(jd)
c
c******************************************************************************
c     set the multigrid changes to zero
c
      do 110 k=1,nkb1
      do 110 j=1,njb1+1
      do 110 i=1,nib1
      b1chg(i,j,k) = 0.0
  110 continue
      do 210 k=1,nkb2
      do 210 j=1,njb2+1
      do 210 i=1,nib2
      b2chg(i,j,k) = 0.0
  210 continue
      do 310 j = 1,nsblk
      sbchg(j) = 0.0
  310 continue
c*******************************************************************************
c     balance fluxes across the tip gap
c
      do 520 j=2,jm
      nr = nrow(j)
      if(ktips(nr).le.0) go to 520
      do 521      k = ktips(nr),ktipe(nr)-1
      tflux(1,j,k)  = 0.5*(tflux(1,j,k)+tflux(im,j,k))
      tflux(im,j,k) = tflux(1,j,k)
  521 continue
  520 continue
c
c***********************************************************
c     extrapolate the fluxes on the upstream face of the mixing planes. 
c     unless fextrap = 0.0
c 
      if(fextrap.lt.0.001) go to 4141

      do 4140 nr = 1,nrwsm1
      j = jmix(nr)
      jp1 = j+1
      jm1 = j-1
      jp2 = j+2
c
      do 4100 k = 1,kmm1
c
      flowdirn = -flowx(imid,j,k)
c
      avflxjm1 = 0.0
      avflxjp2 = 0.0
      do 4120 i = 1,imm1
      avflxjm1   = avflxjm1  + xflux(i,jm1,k)
      avflxjp2   = avflxjp2  + xflux(i,jp2,k)
 4120 continue
     
      do 4130 i=1,imm1
c
      if(flowdirn.gt.0.0) then
           dfluxjm1 = xflux(i,jm1,k) - avflxjm1*fp(i)
           dflux  = dfluxjm1*fextrap
           source(i,j,k) = source(i,j,k) - dflux
      else
           dfluxjp2 = xflux(i,jp2,k) - avflxjp2*fp(i)
           dflux = dfluxjp2*fextrap
           source(i,jp2,k) = source(i,jp2,k) + dflux
      end if
c
      source(i,jp1,k) = 0.0
c
 4130 continue
c
 4100 continue
c
 4140 continue

 4141 continue
c*******************************************************************************
c*******************************************************************************
c      sum the fluxes to form change and save it in store(i,j,k) 
c
      do 1000 k=1,kmm1
      do 1000 j=2,jm
      ratpitch = float(nblade(j-1))/nblade(j)
      do 1000 i=1,imm1
      delta        = xflux(i,j,k) - xflux(i,j-1,k)*ratpitch
     &             + rflux(i,j,k) - rflux(i,j,k+1)
     &             + tflux(i,j,k) - tflux(i+1,j,k)
     &             - source(i,j,k)
      store(i,j,k)  = f1*delta + f2*diff(i,j,k)
      diff(i,j,k)   = delta    + f3*diff(i,j,k)
 1000 continue
c
c**********************************************************************************
c**********************************************************************************
c     pitchwise average the changes at any mixing planes
c     this is only done on the downstream face of the mixing plane unless fextrap = 0
c     in which case it is done on both sides as in tblock-13.
c
      do 1750 nr = 1,nrwsm1
      j   = jmix(nr)
      do 1790 k = 1,kmm1
c
      flowdirn = -flowx(imid,j,k)
c    firstly at j = jmix + 2
      if((flowdirn.gt.0.0).or.(fextrap.lt.0.001)) then
      jmx = j + 2
      sum_store = 0.0
      do 1770 i=1,imm1
      sum_store  = sum_store  + store(i,jmx,k)
 1770 continue
      do 1780 i=1,imm1
         store(i,jmx,k) = sum_store*fp(i)
 1780 continue
      end if
c  next if fextrap = 0  at j = jmix
      if((flowdirn.lt.0.0).or.(fextrap.lt.0.001)) then
      sum_store = 0.0
      do 1771 i=1,imm1
      sum_store  = sum_store  + store(i,j,k)
 1771 continue
      do 1781 i=1,imm1
         store(i,j,k) = sum_store*fp(i)
 1781 continue
      end if
c
 1790 continue
c
 1750 continue
c
c*******************************************************************************
c*******************************************************************************
c    jump to 1020  if no multigrid. this is very unusual.
c
      if(ir.le.1.and.jr.le.1.and.kr.le.1) go to 1020
c
c*******************************************************************************
c*******************************************************************************
c      sum the element changes to form the changes for the multigrid
c      blocks. store(i,j,k) is now used as a store for the change in
c      the elements.
c
      do 700 k=1,kmm1
      k1 = kb1(k)
      k2 = kb2(k)
      do 700 j=2,jm
      jsb= jsblk(j)
      j1 = jb1(j)
      j2 = jb2(j)
      do 700 i=1,imm1
      i1 = ib1(i)
      i2 = ib2(i)
      delta = store(i,j,k)
      b1chg(i1,j1,k1) = b1chg(i1,j1,k1) + delta
      b2chg(i2,j2,k2) = b2chg(i2,j2,k2) + delta
      sbchg(jsb)      = sbchg(jsb)      + delta
  700 continue
c
 1020 continue
c
c*******************************************************************************
c*******************************************************************************
c     add the block changes to the element changes, multiply by the time step
c     jdd removed the calculation of the average change from this april 2018.
c     it is now calculated in the do 1501 loop.
c
      do 1500 k=1,kmm1
      k1 = kb1(k)
      k2 = kb2(k)
      do 1500 j=2,jm
      jsb= jsblk(j)
      j1 = jb1(j)
      j2 = jb2(j)
      do 1500 i=1,imm1
      i1 = ib1(i)
      i2 = ib2(i)
      delta = store(i,j,k)*step(i,j,k)
     &      + (b1chg(i1,j1,k1)*step1(i1,j1,k1)
     &      + b2chg(i2,j2,k2)*step2(i2,j2,k2)
     &      + sbchg(jsb)*stepsbk(jsb))*rstep(i,j,k)
      store(i,j,k) = delta
 1500 continue
c
c      end of calculating the multigrid changes
c*******************************************************************************
c*******************************************************************************
c    use residual smoothing if nrsmth > 0
c
      if(nrsmth.gt.0) then
      call smooth_resid(store,rsmth,nrsmth)
      end if
c
c   
c*******************************************************************************
c*******************************************************************************
c
c     apply the negative feedback. skip it if damp is small or large
c
      if(damp.lt.2.0.or.damp.gt.100) go to 1550
c
c     changed by jdd to use the average change per row - avg_blk(nr)- 14/02/2018
c     calculate the average changes for each row,  avg_chg(nr).
c     this is new by jdd  14/02/2018.
c
      do 1501 nr = 1,nrows
      sumchg = 0.0
      jst = jstart(nr) + 1 
      jen = jmix(nr)   - 1
      jchange = jen - jst  + 1
      nsum = imm1*kmm1*jchange
      do 1502 k=1,kmm1
      do 1502 j = jst, jen
      do 1502 i=1,imm1
      sumchg = sumchg + abs(store(i,j,k))
 1502 continue
      avg_chg(nr) = sumchg/nsum
 1501 continue
c
c     smooth the blade row changes
c
      if(nrows.eq.1) avg_blk(1) = avg_chg(1)
      if(nrows.eq.2) then
            avg_blk(2) = 0.5*(avg_chg(1) + avg_chg(2))
            avg_blk(1) = avg_blk(2)
      end if
      if(nrows.gt.2) then
      do 1503 nr = 2,nrows-1
      avg_blk(nr) = 0.25*(avg_chg(nr-1)+avg_chg(nr+1)) + 0.5*avg_chg(nr)
 1503 continue
      avg_blk(nrows) = 0.5*(avg_chg(nrows-1) + avg_chg(nrows))
      avg_blk(1)     = 0.5*(avg_chg(1) + avg_chg(2)) 
      end if 
c
c*******************************************************************************
c      apply the negative feedback to limit the maximum change.
c*******************************************************************************
c
      do 1525 k=1,kmm1
      do 1525 j=2,jm
      nr     = nrow(j)
      do 1525 i=1,imm1
      delta  = store(i,j,k)
      abschg = abs(delta)
      fdamp  = abschg/avg_blk(nr)
      store(i,j,k) = delta/(1. + fdamp/damp )
 1525 continue
c
c   end of jdd  14/02/2018  changes.
c*******************************************************************************
c*******************************************************************************
c
 1550 continue
c   end of jdd 14/02/2018 changes
c*****************************************************************************
c******************************************************************************
c   next is special treatment for the turbulent viscosity if using the sa model
c
      if(ncall.eq.2) then
c  
c     update the turbulent viscosity - rememberng that it is cell centred.
c    
      do 1600 k=1,kmm1
      do 1600 j=2,jm
      do 1600 i=1,imm1
      d(i,j,k) = d(i,j,k) + store(i,j,k)
 1600 continue
c
c   transfer the average turbulent viscosity across the mixing plane
c
      do 1560 nr = 1,nrows-1
      j = jmix(nr) 
      fmult = turbvis_damp(nr)
      do 1570 k=1,kmm1
      avg(k) = 0.0
      do 1580 i=1,imm1
      avg(k) = avg(k) + fp(i)*d(i,j,k)
 1580 continue
      js = j + 1
      je = j + 2
      do 1591 jav = js,je
      do 1590 i=1,imm1
 1590 d(i,jav,k) = fmult*avg(k)
 1591 continue
 1570 continue
 1560 continue      
c
c******************************************************************************
c   smooth the turbulent viscosity
c
      nsmth = 1
      sftvis  = fac_sfvis*sft
      call smooth_resid(d,sftvis,nsmth)
c
c     jump to 8700 for the turbulent viscosity only
      go to 8700
c
c   end of special treatment for the turbulent viscosity
c
      end if
c
c**************************************************************************************
c**************************************************************************************
c     again pitchwise average the changes at any mixing planes after including 
c     multigrid and multiplying by the timestep. 
c     this is only done on the downstream face of the mixing plane unless  fextrap = 0
c     in which case it is done on both sides as in tblock-13.
c
      do 1850 nr = 1,nrwsm1
      j   = jmix(nr)
      jp1 = j+1
      do 1890 k = 1,kmm1
c
      flowdirn = -flowx(imid,j,k)
c
      if((flowdirn.gt.0.0).or.(fextrap.lt.0.001)) then
      jmx = j + 2
      sum_store = 0.0
      do 1870 i=1,imm1
      sum_store  = sum_store  + store(i,jmx,k)
 1870 continue
      sum_store = sum_store/imm1
      do 1880 i=1,imm1
         store(i,jmx,k) = sum_store
 1880 continue
      end if
c    if fextrap = 0  also at the cells upstream of the mixing plane, j = jmix  
      if((flowdirn.lt.0.0).or.(fextrap.lt.0.001)) then
      sum_store = 0.0
      do 1871 i=1,imm1
      sum_store  = sum_store  + store(i,j,k)
 1871 continue
      sum_store = sum_store/imm1
      do 1881 i=1,imm1
         store(i,j,k) = sum_store
 1881 continue
      end if
c
 1890 continue
c
 1850 continue
c
c****************************************************************************
c****************************************************************************
c          add the changes to the old value of the variable d
c              distributing the changes to the four corners
c              with double weighting at the boundaries.
c              the factor of 1/8 is included in the 'fmi' terms.
c****************************************************************************
c
      do 1100 k=1,kmm1
      do 1100 j=2,jm
      do 1100 i=1,imm1
c
      add = store(i,j,k)
c
      d(i,j,k)      =  d(i,j,k)       + add*fbl(i,k)*facdwn(j)
      d(i+1,j,k)    =  d(i+1,j,k)     + add*fbr(i,k)*facdwn(j)
      d(i,j,k+1)    =  d(i,j,k+1)     + add*ftl(i,k)*facdwn(j)
      d(i+1,j,k+1)  =  d(i+1,j,k+1)   + add*ftr(i,k)*facdwn(j)
      d(i,j-1,k)    =  d(i,j-1,k)     + add*fbl(i,k)*facup(j)
      d(i+1,j-1,k)  =  d(i+1,j-1,k)   + add*fbr(i,k)*facup(j)
      d(i,j-1,k+1)  =  d(i,j-1,k+1)   + add*ftl(i,k)*facup(j)
      d(i+1,j-1,k+1)=  d(i+1,j-1,k+1) + add*ftr(i,k)*facup(j)
 1100 continue
c
c************************************************************************************
c************************************************************************************
c    calll  "smoothvar" to apply the smoothing (artificial viscosity) to the variable "d". 
c
      call smooth_var(d)
c
c*******************************************************************************
c*******************************************************************************
c    re enter here if ncall = 2 for the turbulent viscosity.
 8700 continue
c*******************************************************************************
c*******************************************************************************
c**********apply the periodic boundary conditions at i=1 and im ***************
c  
      ilast = im
      if(ncall.eq.2) ilast = imm1
c
      do 8000 j=1,jm
      if((ind(j).eq.1).or.(indle(j).eq.1)) go to 8000
      do 7000 k=1,km
      d(1,j,k)      = 0.5*(d(1,j,k)+d(ilast,j,k))
      d(ilast,j,k) = d(1,j,k)
 7000 continue
 8000 continue
c
c*******************************************************************************
c*******************************************************************************
c  q3d
      if(km.eq.2) go to 509
c  end q3d
c*******************************************************************************
c*******************************************************************************
c
c     apply periodicity across the tip gap.
c     note that the tip point itself,  ktip,  is not periodic.
c
      do 510 j=1,jm
      nr = nrow(j)
      if(ktips(nr).le.0) go to 510
      k1 = ktips(nr)
      k2 = ktipe(nr)
      if(k1.eq.1)  k2 = k2-1
      if(k2.eq.km) k1 = k1+1
      do 511 k=k1,k2
      d(1,j,k)     = 0.5*(d(1,j,k)+d(ilast,j,k))
      d(ilast,j,k) = d(1,j,k)
  511 continue
  510 continue
c
c*******************************************************************************
c*******************************************************************************
  509 continue
c*******************************************************************************
c*******************************************************************************
c
c    pitchwise average the variables on both faces of the mixing plane
c
      if(ifmix.eq.0) go to 1400
c
      do 100 nr = 1,nrwsm1
      j   = jmix(nr)
      jp1 = j+1
c
      do 201 k=1,km
      avgvarj    = 0.0
      avgvarjp1  = 0.0
      do 301 i=1,imm1
      avgvarj   = avgvarj   + fp(i)*(d(i,j,k)+d(i+1,j,k)) 
      avgvarjp1 = avgvarjp1 + fp(i)*(d(i,jp1,k)+d(i+1,jp1,k))  
  301 continue
      avgvarj   = 0.5*avgvarj
      avgvarjp1 = 0.5*avgvarjp1
      do 401 i=1,im
      d(i,j,k)   = avgvarj
      d(i,jp1,k) = avgvarjp1
  401 continue 
  201 continue
c 
  100 continue
c
c     end of  mixing plane treatment
c
 1400 continue
c
c*******************************************************************************
c*******************************************************************************
c
      return
      end
c*******************************************************************************
c*******************************************************************************
c
      subroutine output
c       ====================
c
c       this routine calculates flow properties and prints headings
c       before calling print to output the result arrays
c
      include  'commall-open-19.2'
c
c
      write(6,*) 'in output, writing plotting/restart file "flow_out" ' 
c
c        write  the combined restart/plotting  file   "flow_out"  if  "ifprint"  = 1
c        or  if  "ifend " = 1 when the program has converged.
c
      preff = 0.5*(p(imid,1,kmid) + p(imid,jm,kmid))
      write(7) nstep
      write(7)(((ro(i,j,k),i=1,im),j=1,jm),k=1,km)
      write(7)(((rovx(i,j,k),i=1,im),j=1,jm),k=1,km)
      write(7)(((rovr(i,j,k),i=1,im),j=1,jm),k=1,km)
      write(7)(((rovt(i,j,k),i=1,im),j=1,jm),k=1,km)
c     correct the value of roe  to the true value before sending it
c     to output file if using aretificial compressibility.
      if(itimst.eq.5.or.itimst.eq.6) then
                write(7)(((roe(i,j,k)*ro(i,j,k)/rosub(i,j,k)
     &         ,i=1,im),j=1,jm),k=1,km)
      else
                write(7)(((roe(i,j,k),i=1,im),j=1,jm),k=1,km)
      end if
c
      write(7)(((rosub(i,j,k),i=1,im),j=1,jm),k=1,km)
c
c     write out the transformed viscosity if using the sa model
c     otherwise write out the turbulent/laminar viscosity ratio for plotting only.
      if(ilos.ge.200) then
          write(7)(((trans_dyn_vis(i,j,k),i=1,im),j=1,jm),k=1,km)
      else
          write(7) (((visc_rat(i,j,k),i=1,im),j=1,jm),k=1,km)
      end if
c
c      write out the peff for the surface pressures in a throughflow calculation.
c      write(7)((( (peff(i,j,k)+preff),i=1,im),j=1,jm),k=1,km)
c      write out sparevar  which can be set to any required 3d varible.
 
      write(7) (((sparevar(i,j,k),i=1,im),j=1,jm),k=1,km)
c
      write(6,*) 'in output, plotting/restart file "flow_out" written '
c
c********************************************************************************
c********************************************************************************
c********************************************************************************
c
c    now  write out the formatted output file  "results.out" to unit 3 
c
c********************************************************************************
c********************************************************************************
c********************************************************************************
c
c     calculate some flow parameters for printing out 
c
c      temp1 = relative swirl velocity
c      temp2 = relative mach number
c      temp3 = pressure in bar
c      temp4 = absolute stagnation temperature
c     
      do 1100 i=1,im
      do 1100 k=1,km
      do 1100 j=1,jm
      eke = 0.5*(vx(i,j,k)*vx(i,j,k) + vt(i,j,k)*vt(i,j,k)
     &         + vr(i,j,k)*vr(i,j,k))
c
      if(ifgas.eq.0) then
            tstatic = (ho(i,j,k) - eke)*rcp
            if(tstatic.lt.1.) tstatic = 1.
            v_sonic = sqrt(ga*rgas*tstatic)
            temp4(i,j,k)  = ho(i,j,k)*rcp
      else
            hstag   = ho(i,j,k)
            hstat   = hstag - eke
            tstag   = tfromh(hstag,tref,href,ht1,ht2,ht3,ht4)
            tstat   = tfromh(hstat,tref,href,ht1,ht2,ht3,ht4)
            temp4(i,j,k) = tstag
            if(tstat.lt.1.) tstat = 1.
            cpnow   = cp1 + cp2*(tstat-tref) 
     &              + cp3*(tstat-tref)*(tstat-tref)
            gamnow  = cpnow/(cpnow-rgas)
            v_sonic = sqrt(gamnow*rgas*tstat)
      end if
c
      temp1(i,j,k)  = vt(i,j,k) - ublade(j,k)
      ekerel  = eke - 0.5*(vt(i,j,k)*vt(i,j,k)
     &              - temp1(i,j,k)*temp1(i,j,k))
      temp3(i,j,k)  = 0.00001*p(i,j,k)
      temp2(i,j,k)  = sqrt (2.0*ekerel)/v_sonic
 1100 continue
c
c***************************************************************************
c***************************************************************************
c    store is the percentage change in meridional velocity in the last time step.
c
      write(3,5) title
    5 format(1h ,24x,18a4/)
      if(iout(1).eq.0) go to 10
      iopt=iout(1)
      write(3,5400) nstep
 5400 format(1h ,'** percentage change in vm ** timestep =',i5,' **')
      call print (iopt,store)
c
   10 continue
c
      if(iout(2).eq.0) go to 20
      iopt=iout(2)
      write(3,5100) nstep
 5100 format(1h ,'*** axial velocity  m/sec timestep =',i5,' **')
      call print(iopt,vx)
c
   20 continue
c
      if(iout(3).eq.0) go to 30
      iopt=iout(3)
      write(3,5200) nstep
 5200 format(1h ,' radial velocity m/sec timestep =',i5,' **')
      call print(iopt,vr)
c
   30 continue
c
      if(iout(4).eq.0) go to 40
      iopt=iout(4)
      write(3,5300) nstep
 5300 format(1h0,' relative swirl  velocity m/sec timestep =',i5,' **')
      call print(iopt,temp1)
c
   40 continue
c
      if(iout(5).eq.0) go to 50
      iopt=iout(5)
      write(3,5500) nstep
 5500 format(1h0,'********** pressure in bar *** timestep =',i5,' **')
      call print(iopt,temp3(i,j,k))
c
   50 continue
c
      if(iout(6).eq.0) go to 60
      iopt=iout(6)
      write(3,5600) nstep
 5600 format(1h ,' ** relative  mach  number  ** timestep =',i5,' **')
      call print(iopt,temp2)
c
   60 continue
c
      if(iout(7).eq.0) go to 70
      write(3,5700) nstep
      iopt=iout(7)
 5700 format(1h0,'absolute stag. temperature deg k timestep= ',i5,'**')
      call print(iopt,temp4)
c
c     calculate the entropy function
c
      do 1200 i=1,im
      do 1200 k=1,km
      do 1200 j=1,jm
      d=.5*(vx(i,j,k)*vx(i,j,k)+vt(i,j,k)*vt(i,j,k)+vr(i,j,k)*vr(i,j,k))
c
      if(ifgas.eq.0) then
            tstatic = (ho(i,j,k) - d)*rcp
            if(tstatic.lt.1.) tstatic = 1.
            store(i,j,k) = p(i,j,k)/po1(kmid)*
     &      (to1(kmid)/tstatic)**(ga/(ga-1))
      else
           tstag  = temp4(i,j,k)
           hstat  = hstag - d
           tstat  = tfromh(hstat,tref,href,ht1,ht2,ht3,ht4)
           trat   = tstat/tstag
           prat   = prat_from_trat(to1(kmid),trat,alpha,beta1,beta2)
           store(i,j,k) = p(i,j,k)/(po1(kmid)*prat)
      end if
c
 1200 continue
c
   70 continue
c
      if(iout(12).eq.0) go to 75
      iopt=iout(12)
      write(3,5710) nstep
 5710 format(1h0,' p/(t**(ga/(ga-1)))/ inlet value on mid k grid line,
     &  timestep no ',i5)
      call print(iopt,store)
c
   75 continue
c*******************************************************************************
c    calculate more properties for printing out
c
c          store is meridional velocity
c          temp1 is relative swirl angle, atan(vtrel/vm)
c          temp2 is meridional pitch angle, atan(vr/vx)
c
      pin = 0.0
      do 2000 k=1,km
      do 2010 i=1,im
      do 2020 j=1,jm
      store(i,j,k) = sqrt(vx(i,j,k)*vx(i,j,k)+vr(i,j,k)*vr(i,j,k))
      temp1(i,j,k) = atan(temp1(i,j,k)/store(i,j,k)) * 57.296
      temp2(i,j,k) = asin(vr(i,j,k)/store(i,j,k)) * 57.296
 2020 continue
      pin = pin + p(i,1,k)
 2010 continue
 2000 continue
c**********************************************************************************
c
      if(iout(8).eq.0) go to 80
      iopt=iout(8)
      write(3,2030) nstep
 2030 format(1h0,' meridional velocity, m/sec, timestep=',i5 )
      call print(iopt,store)
c
   80 continue
c
      if(iout(9).eq.0) go to 90
      iopt=iout(9)
      write(3,2040) nstep
 2040 format(1h0,' swirl angle on stream surface, atan(wtrel/vm), deg,
     1 timestep=',i5)
      call print(iopt,temp1)
c
   90 continue
c
      if(iout(10).eq.0) go to 100
      iopt=iout(10)
      write(3,2050) nstep
 2050 format(1h0,' meridional pitch angle, atan(vr/vx), deg, timestep=
     1 ',i5)
      call print(iopt,temp2)
c
  100 continue
c
      if(iout(11).eq.0) go to 110
      iopt=iout(11)
      write(3,2060) nstep
 2060 format(1h0,' density kg/m**3, timestep=   ',i5)
      call print(iopt,ro)
c
  110 continue
c
      if(iout(12).eq.0) go to 120
      iopt=iout(12)
      write(3,2070) nstep
 2070 format(1h0,' artificial density kg/m**3, timestep=   ',i5)
      call print(iopt,rosub)
c
  120 continue
c
      if(iout(13).eq.0) go to 130
      iopt=iout(13)
      pin=pin/(im*km)
      do 113 k=1,km
      do 113 i=1,im
      do 113 j=1,jm
      store(i,j,k)=(p(i,j,k)-pin)/(po1(kmid)-pin)
  113 continue
      write(3,2080) nstep
 2080 format(1h0, ' pressure coefficient based on  average inlet
     & static pressure and stagnation pressure at kmid',i5)
      call print(iopt,store)
c
  130 continue
c
c    end of printed output
c
c************************************************************************************
c************************************************************************************
c
      return
      end
c
c************************************************************************************
c************************************************************************************
c
      subroutine print(iopt,f)
c
c       ====================
c       this routine print's out 3-d arrays or pitchwise average values.
c       ====================
c
      include  'commall-open-19.2'
c
c
      dimension f(id,jd,kd),sumfun(maxki)
c
      if(iopt.eq.1) go to 4000
      if(iopt.eq.3) go to 4000
      if(iopt.eq.2) go to 500
      return
c**********************************************************************
c      print out whole flow field of the input variable  "f(i,j,k)" .
c
 4000 continue
c
      do 1000 k=1,km
      if(kout(k).eq.0 )go to 1000
      write(3,5200) k
 5200 format(/' stream surface number  ',i4,/)
      do 1100 j=1,jm
      write(3,5300) j,r(j,k),x(j,k),smerid(j,k),(f(i,j,k),i=1,im)
 1100 continue
 1000 continue
 5300 format(1h ,'j= ',i4,'r= ',f10.5,' x= ',f10.5,' meridional distance
     & =',f10.5,' value i=1->im =',/,(10f11.3))
c
      return
c
c************************************************************************
c      print out pitchwise mass averaged values if requested. the values are the 
c      averaged values in each cell so only km-1 values
c
  500 continue
c
      write(3,2000)
 2000 format( /,'   fractional span of each k line ' )
      write(3,2001) (fspan(k),k=1,kmm1)
 2001 format(10f12.4)
      write(3,*)
      write(3,*)
c
      do 100 j=1,jm
      do 110 k=1,kmm1
      sumas     = 0.0
      sumfun(k) = 0.0
      do 120 i=1,imm1
      sumas     = sumas + flowx(i,j,k)
      sumfun(k) = sumfun(k)+flowx(i,j,k)*(f(i,j,k)+f(i+1,j,k)
     &          + f(i+1,j,k+1)+f(i,j,k+1))*0.25
  120 continue
      sumfun(k) = sumfun(k)/sumas
  110 continue
c
      write(3,200) j,r(j,1),r(j,km),(sumfun(k),k=1,kmm1)
  200 format( /,' j=',i5,' rhub=',f10.4,' rtip=',f10.4,' pitchwise mass 
     & average, k=1,kmm1,= ',/, (10f12.3))
c
  100 continue
c
      return
      end
c******************************************************************************
c******************************************************************************
c
      subroutine setup(ansin)
c
c**********this subroutine sets up the grid and initialises all the
c          variables used in the main program.
c
      include  'commall-open-19.2'
c
      dimension jstart1(jd),jend1(jd),jstart2(jd),jend2(jd),pdown(kd),
     &          indlete(jd)
c
      character*1 answ,ansin
c
c      set various constants and integers needed throughout the calculation
c
      write(6,1111)
 1111 format('  entering subroutine setup so data input was ok ')
c
      jlep1  = jle(1)+1
      nrwsm1 = nrows-1
      imid   = ifix(0.5*im)
      if(im.eq.2) imid = 1
      kmid   = ifix(0.5*km)
      if(km.eq.2) kmid = 1
      ifprint= 0
      ifend  = 0
c
      do j=1,jm
           indlete(j) = 0
           if(indle(j).eq.1) indlete(j) = 1
           if(indte(j).eq.1) indlete(j) = 1
      end do 
c  q3d
      if(km.eq.2) then
           kmid = 1
           kmm2 = 1
      end if
c  end q3d

      do 1000 k=1,km
      bssin(k) = sin(bs(k)*degrad)
      bscos(k) = cos(bs(k)*degrad)
      brsin(k) = sin(br(k)*degrad)
 1000 brcos(k) = cos(br(k)*degrad)
c
c******************************************************************************
c     make the fr(k) sum to unity
c
      sum = 0.0
      do 2000 k=1,kmm1
 2000 sum = sum+fr(k)
      do 2001 k=1,kmm1
 2001 fr(k) = fr(k)/sum
c
c    fspan is the fractional span at the centre of the elements.
      fspan(1) = 0.5*fr(1)
      do 1999 k=2,kmm1
 1999 fspan(k) = fspan(k-1)+0.5*(fr(k)+fr(k-1))
c
c    make the fp(i) sum to unity
      sum=0.0
      do 2002 i=1,imm1
 2002 sum = sum+fp(i)
      do 2003 i=1,imm1
 2003 fp(i) = fp(i)/sum
c 
c  the fu(i)  and fd(i)  are used in the smoothing routine.
c
      if(im.gt.2) then
c
      do 2004 i=2,imm1
      fu(i) = fp(i)/(fp(i)+fp(i-1))
 2004 fd(i) = fp(i-1)/(fp(i)+fp(i-1))
      fu(1) = fp(1)/fp(2)
      fu(im)= fp(imm1)/fp(imm2)  
c
      else
c   for tflow
      fu(1) = 0.5
      fd(1) = 0.5
      fu(2) = 0.5
      fd(2) = 0.5
c
      end if    
c
c      set fku(k)  and  fkd(k) which  are used in the smoothing routine.
c
      if(km.gt.2) then
           do 2005 k=2,kmm1
           fku(k) = fr(k)/(fr(k)+fr(k-1))
 2005      fkd(k) = fr(k-1)/(fr(k)+fr(k-1))
           fku(1) = fr(1)/fr(2)
           fku(km)= fr(kmm1)/fr(kmm2)
      else
c   q3d
           fku(1) = 0.5
           fku(2) = 0.5
      end if
c   end q3d
c******************************************************************************
c
c     set the distribution functions for use in subroutine tstep
c
      do 2006 k=1,kmm1
      do 2006 i=1,imm1
           ftl(i,k)  = 0.125
           fbl(i,k)  = 0.125
           ftr(i,k)  = 0.125
           fbr(i,k)  = 0.125
 2006 continue
      do 2007 k=1,kmm1
           ftl(1,k)    = 0.25
           fbl(1,k)    = 0.25
           ftr(imm1,k) = 0.25
           fbr(imm1,k) = 0.25
 2007 continue
      do 2008 i=1,imm1
           fbl(i,1)    = 0.25
           fbr(i,1)    = 0.25
           ftl(i,kmm1) = 0.25
           ftr(i,kmm1) = 0.25
 2008 continue
           ftl(1,kmm1)    = 0.5
           fbl(1,1)       = 0.5
           ftr(imm1,kmm1) = 0.5
           fbr(imm1,1)    = 0.5
      do j = 2,jm
           facup(j)  = 1.0
           facdwn(j) = 1.0
      end do
           facup(2)    = 2.0
           facdwn(jm)  = 2.0
c
c********************************************************************************
c********************************************************************************
c
c      call intpol to interpolate in input data to set up grid nodes on the
c      blade surfaces.
c      if km = 2 then the two input stream surfaces are the hub and casing
c      and there is no interpolation.
c  q3d
      if(km.eq.2) then
           do k=1,2
           do j=1,jm
                x(j,k) = xsurf(j,k)
                r(j,k) = rsurf(j,k)
           end do
           end do
      else
c
c  end q3d
c
      write(6,*) ' input style, ansin = ', ansin
c
      if(ansin.eq.'n'.or.ansin.eq.'n')then
          write(6,*) ' calling new_intpol from setup'
          call new_intpol
      else
          write(6,*) ' calling old_intpol from setup'          
          call old_intpol
      end if
c
      write(6,1888)
 1888 format( '  interpolation in input blade sections completed ok ')
c
      endif
c
c******************************************************************************
c******************************************************************************
c     shift the blades so they are aligned circumferentially at i=1 at mid-span
c
      shift = 0.0
      do 3022 j=3,jm
      if(indmix(j-1).eq.1) shift = (rt_upp(j-1,kmid)
     &                           -  rt_upp(j,kmid))/r(j,kmid)
      do 3022 k=1,km
      rt_upp(j,k) = rt_upp(j,k) + shift*r(j,k)
 3022 continue
c
c******************************************************************************
c     scale the blade thickness for the pinched tip model of tip clearance.
c
      do 2506 j=1,jm
      nr = nrow(j)
      do 2505 k=1,km
      rtmid = rt_upp(j,k) - 0.5*rt_thick(j,k)
      thick = fthick(nr,k)*rt_thick(j,k)
      rt_upp(j,k)    = rtmid + 0.5*thick
      rt_thick(j,k)  = thick
 2505 continue
 2506 continue
c
c*******************************************************************************
c      set various constants
c
      ga1   = ga - 1.0
      fga   = ga1/ga
      rfga  = 1.0/fga
      rcp   = 1.0/cp
      cv    = cp/ga
      rcv   = 1.0/cv
      rgas  = cp - cv
      sfex1 = 1. - sfexit
      sfexh = 0.5*sfexit
c
c*******************************************************************************
c          set temporary grid coordinates for use in working out areas.
c
c       theta .......store for theta co-ord  of grid nodes
c       rtheta ......store for rtheta co-ord of grid nodes
c       dx ..........axial extent of elements
c       ds ... ......meridional extent of elements.
c       rt_pitch ... local blade pitch.
c
      do 3000 j=1,jm
      pitch = 2.*pi/nblade(j)
      do 3021 k=1,km
      gap           = pitch*r(j,k) - rt_thick(j,k)
      rtheta(1,j,k) = rt_upp(j,k)
      theta(1,j,k)  = rt_upp(j,k)/r(j,k)
      rt_pitch(j,k) = pitch*r(j,k)
      do 3020 i=2,im
      rtheta(i,j,k) = rtheta(i-1,j,k) + fp(i-1)*gap
      theta(i,j,k)  = rtheta(i,j,k)/r(j,k)
 3020 continue
 
 3021 continue
 3000 continue
c
c*******************************************************************************
c      evaluate the meridional distance, smerid.
c
      dsmin=1.0e06
      do 3530 k=1,km
      smerid(1,k) = 0.0
      do 3540 j=2,jm
      dr(j,k)     = (r(j,k)-r(j-1,k))
      dx(j,k)     = (x(j,k)-x(j-1,k))
      ds(j,k)     = sqrt(dr(j,k)*dr(j,k) + dx(j,k)*dx(j,k))
      smerid(j,k) = smerid(j-1,k) + ds(j,k)
      if(ds(j,k).lt.dsmin) dsmin = ds(j,k)
 3540 continue
      dx(1,k) = dx(2,k)
      ds(1,k) = ds(2,k)
      dr(1,k) = dr(2,k)
 3530 continue
c
c   jdd addition august 2017
c   find the blade chords for every row
      do 3550 nr = 1,nrows
      jledge    = jle(nr)
      jtedge    = jte(nr)
      schord    = smerid(jtedge,kmid)-smerid(jledge,kmid)
      tchord    = rtheta(1,jledge,kmid) - rtheta(1,jtedge,kmid)
      chord(nr) = sqrt(schord*schord + tchord*tchord)
 3550 continue
c
c      write "smerid" to the file "globplot"
c
      write(11)    jm,1,1
      write(11)    (smerid(j,kmid), j=1,jm)
c
c*******************************************************************************
c   set ravg_cell(j,k) = average radius of an element.
c
      do 4000 j=2,jm
      do 4000 k=1,kmm1 
      ravg_cell(j,k) = 0.25*(r(j,k)+r(j-1,k)+r(j-1,k+1)+r(j,k+1))
 4000 continue
c
c     set the blade speed  ublade(j,k)
c
      do 4005 j=1,jm
      do 4005 k=1,km
      ublade(j,k) = wrad(j)*r(j,k)
 4005 continue
c
c*******************************************************************************
c   set fmup ,fmdn  for use in the smoothing routine.
c
      do 4010 k=1,km
      do 4010 j=2,jmm1
      dsup      = smerid(j,k)   - smerid(j-1,k)
      dsdwn     = smerid(j+1,k) - smerid(j,k)
      fmdn(j,k) = dsdwn/(dsup + dsdwn)
      fmup(j,k) =  dsup/(dsup + dsdwn)
 4010 continue      
c
c*******************************************************************************
c      store extrapolated values of theta and r*theta on the mixing plane
c      temporarily as temp1 and temp2
c
      do 3025 j=1,jm
      if(indmix(j).ne.1) go to  3025
      do 3024 k=1,km
      fac = ds(j+1,k)/ds(j+2,k)
      do 3024 i=1,im
      temp1(i,j,k) = theta(i,j+1,k)+ fac*(theta(i,j+1,k)-theta(i,j+2,k))
      temp2(i,j,k) = r(j,k)*temp1(i,j,k)
 3024 continue
 3025 continue
c
c******************************************************************************
c******************************************************************************
c     write the grid geometry to unit 21 for use when plotting
c
      open(unit=21,file='grid_out',form= 'unformatted' )
c
      write(21) nsteps_max
      write(21) im,jm,km
      write(21) cp,ga
      write(21) (indlete(j),j=1,jm)
      write(21) (wrad(j),j=1,jm)
      write(21) (nblade(j),j=1,jm)
c 
      do 20 j=1,jm
      do 20 k=1,km 
      write(21) x(j,k),r(j,k),(rtheta(i,j,k),i=1,im)
   20 continue
c
      close(21)
c
c******************************************************************************
c******************************************************************************
c           work out the projected areas of the blade faces
c
      write(6,*) ' starting to work out areas of faces of elements'
c
c           quasi orthogonal face first.
c
      do 3030 j=1,jm
      do 3030 i=1,imm1
      do 3015 k=1,kmm1
      x1=x(j,k+1)-x(j,k)
      x2=x1
      r1=r(j,k+1)-r(j,k)
      r2=r1
      t1=r(j,k+1)*(theta(i+1,j,k+1)-theta(i,j,k))
      t2=r(j,k+1)*(theta(i,j,k+1)-theta(i,j,k))
     &    -r(j,k)*(theta(i+1,j,k)-theta(i,j,k))
      aqx(i,j,k) = -0.5*(r2*t1-r1*t2)
      aqr(i,j,k) = -0.5*(x1*t2-x2*t1)
      aqtot(i,j,k) = sqrt(aqx(i,j,k)*aqx(i,j,k)+aqr(i,j,k)*aqr(i,j,k))
 3015 continue
      aqx(i,j,km)=aqx(i,j,kmm1)
 3030 aqr(i,j,km)=aqr(i,j,kmm1)
c
      do 3031 j=1,jm
      do 3031 k=1,km
      aqx(im,j,k) = aqx(imm1,j,k)
      aqr(im,j,k) = aqr(imm1,j,k)
 3031 continue
c
c***********next work out the areas of the bladewise face*********
c
      do 3040 j=2,jm
      do 3040 k=1,kmm1
      x1 = x(j-1,k+1)-x(j,k)
      x2 = x(j,k+1)-x(j-1,k)
      r1 = r(j-1,k+1)-r(j,k)
      r2 = r(j,k+1)-r(j-1,k)
      abt(j,k) = -0.5*(x1*r2-x2*r1)
c
      wrabt(j,k) = wrad(j)*0.25*(r(j,k)+r(j-1,k)+r(j-1,k+1)+r(j,k+1))
     &            *abt(j,k)
c
      do 3040 i=1,im
      t1 = r(j-1,k+1)*(theta(i,j-1,k+1)-theta(i,j,k))
      t2 =   r(j,k+1)*(theta(i,j,k+1)-theta(i,j,k))
     &   -   r(j-1,k)*(theta(i,j-1,k)-theta(i,j,k))
c
      if(indmix(j-1).eq.1) then
      t1 = r(j-1,k+1)*(temp1(i,j-1,k+1)-theta(i,j,k))
      t2 =   r(j,k+1)*(theta(i,j,k+1)-theta(i,j,k))
     &   -   r(j-1,k)*(temp1(i,j-1,k)-theta(i,j,k))
      endif
c
      abx(i,j,k) = -0.5*(r1*t2-r2*t1)
      abr(i,j,k) = -0.5*(x2*t1-x1*t2)
      abtot(i,j,k) = sqrt(abx(i,j,k)**2 + abr(i,j,k)**2 + abt(j,k)**2)
c
      if(i.eq.im) go to 3040
c
c   work out the volume of the elements
c
      vol(i,j,k) = abt(j,k)*0.25*(rtheta(i+1,j,k) + rtheta(i+1,j-1,k)
     &  + rtheta(i+1,j-1,k+1) + rtheta(i+1,j,k+1) - rtheta(i,j,k)
     &  - rtheta(i,j-1,k)     - rtheta(i,j-1,k+1) - rtheta(i,j,k+1))
c
      if(indmix(j-1).eq.1)
     &    vol(i,j,k) = abt(j,k)*0.25*(rtheta(i+1,j,k) + temp2(i+1,j-1,k)
     &  + temp2(i+1,j-1,k+1) + rtheta(i+1,j,k+1) - rtheta(i,j,k)
     &  - temp2(i,j-1,k)     - temp2(i,j-1,k+1)  - rtheta(i,j,k+1))
c
 3040 continue
c
c**********work out areas of streamwise face.
c
      do 3050 j=2,jm
      do 3050 i=1,imm1
      do 3050 k=1,km
      x1=x(j-1,k)-x(j,k)
      x2=x1
      r1=r(j-1,k)-r(j,k)
      r2=r1
      t1 = r(j-1,k)*(theta(i,j-1,k)  - theta(i,j,k))
     &   -   r(j,k)*(theta(i+1,j,k)  - theta(i,j,k))
      t2 = r(j-1,k)*(theta(i+1,j-1,k)- theta(i,j,k))
c
      if(indmix(j-1).eq.1) then
      t1 = r(j-1,k)*(temp1(i,j-1,k)  - theta(i,j,k))
     &   -   r(j,k)*(theta(i+1,j,k)  - theta(i,j,k))
      t2 = r(j-1,k)*(temp1(i+1,j-1,k)- theta(i,j,k))
      endif
c
      asx(i,j,k)=0.5*(r1*t2-r2*t1)
      asr(i,j,k)=0.5*(x2*t1-x1*t2)
c
 3050 continue
c
c**********************************************************************************
c     set the volume/radius term needed for the radial momentum equation.
c
      nneg = 0
      do 4045 k=1,kmm1
      do 4045 j=2,jm
      do 4045 i=1,imm1
      if(vol(i,j,k).lt.0.0) nneg = nneg + 1
      rarea  = abr(i,j,k) - abr(i+1,j,k)
     &       + aqr(i,j,k) - aqr(i,j-1,k)
     &       + asr(i,j,k) - asr(i,j,k+1)
      volor(i,j,k) = - 0.125*rarea
c      ratio = abs(rarea)/abs(asr(i,j,k))
c      check closure of the volumes in the radial direction
c      if( ratio.gt.0.00001 )then
c      write(6,*) ' rarea not well  closed at  j , k =', j,k
c      write(6,*) ' asr, rarea',  asr(i,j,k), rarea, ratio
c      end if

      xarea  = aqx(i,j,k) - aqx(i,j-1,k)
     &       + abx(i,j,k) - abx(i+1,j,k)
     &       + asx(i,j,k) - asx(i,j,k+1)
c      ratio = abs(xarea)/abs(aqx(i,j,k) )
c      check closure of the volumes in the axial direction
c      if( ratio.gt.0.00001 )then
c      write(6,*) ' xarea not well  closed at  j , k =', j,k
c      write(6,*) ' aqx, xarea',  aqx(i,j,k), xarea, ratio
c      end if

 4045 continue
c
      write(6,*)  ' areas of all faces of the elements evaluated '
c
c**********************************************************************************
c**********************************************************************************
c
c************set the initial guess of p and ro at all grid points********
c
      write(6,*)' starting to make the initial guess of the flow field.'
c
c   set the spanwise variation of the  exit pressure if ipout = 3.
      if(ipout.eq.3) go to 3212
      pd(1) = pdown_hub
      do 3211 k=2,km
      pd(k) = pd(k-1) + fr(k-1)*(pdown_tip - pdown_hub)
 3211 continue
 3212 continue
c
c**********************************************************************************
c     set values needed if using artificial compressibility
c
      if(itimst.ge.5) then
           po_ref  = po1(kmid)
           ro_ref  = po1(kmid)/rgas/to1(kmid)
           if(itimst.eq.6) ro_ref = densty
           dp_dro  = vsound*vsound
      end if
c
c*************************************************************************  
c  set the initial guess of pressure, density and temperature
c
      do 3210 k=1,km
           rpo1(k)   = 1./po1(k)
           vr(1,1,k) = vm1(k)*(r(2,k)-r(1,k))/ds(1,k)
           vinsq     = vm1(k)*vm1(k)+vtin(k)*vtin(k)
c
      if(ifgas.eq.0) then
            t1  = to1(k) - 0.5*vinsq/cp
            p1  = po1(k)*(t1/to1(k))**rfga
      else
            hstag    = hfromt(to1(k),tref,href,cp1,cp2,cp3)
            hstat    = hstag - 0.5*vinsq
            t1       = tfromh(hstat,tref,href,ht1,ht2,ht3,ht4)
            trat     = t1/to1(k)
            prat     = prat_from_trat(to1(k),trat,alpha,beta1,beta2)
            p1       = po1(k)*prat
      end if
c          
      do 3210 j=1,jm
           ps   = pguess(j)
           prat = ps/po1(k)
c
      if(ifgas.eq.0) then
           ts = (prat**fga)*to1(k)
      else
           trat = trat_from_prat(to1(k),prat,fgagas,r_alpha,
     &            balpha1,balpha2)
           ts   = to1(k)*trat
      end if
c
           ros = ps/ts/rgas
           if(itimst.eq.6) ros  = densty
c
      do 3220 i=1,im
           p(i,j,k)     = ps
           ro(i,j,k)    = ros
           if(itimst.eq.5.or.itimst.eq.6)
     &     rosub(i,j,k) = ro_ref - (po_ref - p(i,j,k))/dp_dro
           roavg_cell(i,j,k) = ros
           t_static(i,j,k)   = ts
 3220 continue
 3210 continue
c
c
c**********************************************************************************
c
      write(6,*) ' inital guess of p & ro completed '
c
c**********************************************************************************
c   set the initial guess of velocity components.  make the tangential velocity
c   vary smoothly between bladed and unbladed regions using  "facj" .
c
      kavg = 0.75*km
      do 3300 k=1,km
      do 3310 i=1,im
      vt(i,1,k) = vtin(k)
      wt(i,1,k) = vtin(k) - ublade(1,k)
c
      if(ifgas.eq.0) then
           ho(i,1,k) = cp*to1(k)
      else
           ho(i,1,k) = hfromt(to1(k),tref,href,cp1,cp2,cp3)
      end if
c
      facj = 0.9
      do 3320 j = 2,jm
c
      if(j.gt.2.and.indte(j-2).eq.1) facj = 1.0
c
      area=sqrt(aqx(imid,j,kavg)*aqx(imid,j,kavg)+aqr(imid,j,kavg)
     &         *aqr(imid,j,kavg))*nblade(j)
      if(j.eq.2) ain=area
      vs        = vm1(kavg)*ain*ro(imid,1,kavg)/(area*ro(imid,j,kavg))
      vtgrid = vs*r(j,k)*(theta(i,j,k) - theta(i,j-1,k))/ds(j,k)
     &          + ublade(j,k)
c  relax  vtheta   between the grid angle vaue and the upstream point value
      vt(i,j,k) =  (1.- facj)*vtgrid + facj*vt(i,j-1,k)*r(j-1,k)/r(j,k)
c     omit the value at the mixing plane where the grid angle value is wrong
      if(indmix(j-1).eq.1) vt(i,j,k) = vt(i,j-1,k)
c     use the inlet boundary value upstream of the first leading edge.
      if(j.lt.jle(1))  vt(i,j,k) = vtin(k)*r(1,k)/r(j,k)
c
      wt(i,j,k) = vt(i,j,k) - ublade(j,k)
      vx(i,j,k) = vs*dx(j,k)/ds(j,k)
      vr(i,j,k) = vs*dr(j,k)/ds(j,k)
      ho(i,j,k) = ho(i,j-1,k)+wrad(j)*(r(j,k)*vt(i,j,k)-r(j-1,k)
     &            *vt(i,j-1,k))
c   relax the upstream value by facj = 0.9.
      if(indle(j).eq.1) facj = 0.9
c
      if(k.eq.kmid.and.i.eq.imid) then 
          write(6,3321) j, t_static(i,j,k), p(i,j,k), ro(i,j,k),
     &                  vx(i,j,k), vr(i,j,k), vt(i,j,k) 
      end if
 3321 format('j=  ',i5,'  initial guess of: t, p ,ro, vx, vr, vt',
     &      f10.2, f12.1, 4f10.3)
c
 3320 continue
      vr(i,1,k) = vr(i,2,k)
      vx(i,1,k) = vx(i,2,k)
      wt(i,1,k) = wt(i,2,k)
 3310 continue
 3300 continue
c
c tflow 
      if(im.eq.2) then
           do 3333 k=1,km
           do 3333 j=1,jm
                vtavg = 0.5*(vt(1,j,k) + vt(2,j,k))
                vt(1,j,k) = vtavg
                vt(2,j,k) = vtavg
 3333      continue
      end if
c  end tflow 
c*********************************************************************************
c     set the mass fluxes, rovx, rovr, rorvt  and roe '
c
      preff = 0.5*(p(imid,1,kmid) + p(imid,jm,kmid))
c
      do 3500 j=1,jm
      do 3500 k=1,km
      do 3500 i=1,im
      eke  = 0.5*(vx(i,j,k)**2+vt(i,j,k)**2+vr(i,j,k)**2)
c
      if(ifgas.eq.0) then
           roe(i,j,k)    = p(i,j,k)/(ga-1) + eke*ro(i,j,k)
      else
           estat         = ho(i,j,k) - eke - p(i,j,k)/ro(i,j,k) 
           roe(i,j,k)    = ro(i,j,k)*(estat + eke)
      end if
c
      if(itimst.ge.5) roe(i,j,k)=rosub(i,j,k)*((ho(i,j,k)-eke)/ga + eke)
c
      ros           = ro(i,j,k)
      rovr(i,j,k)   = ros*vr(i,j,k)
      rovt(i,j,k)   = ros*vt(i,j,k)
      rowt(i,j,k)   = ros*wt(i,j,k)
      rorvt(i,j,k)  = rovt(i,j,k)*r(j,k)
      rovx(i,j,k)   = ros*vx(i,j,k)
      peff(i,j,k)   = p(i,j,k) - preff
      dro(i,j,k)    = 0.0
      droe(i,j,k)   = 0.0
      drovx(i,j,k)  = 0.0
      drovr(i,j,k)  = 0.0
      drorvt(i,j,k) = 0.0
      dpds_cell(i,j,k)  = 0.0
      del_dynvis(i,j,k) = 0.0
      y_plus(i,j,k) = 1000.0
      visc_rat(i,j,k) = 1.0
      trans_dyn_vis(i,j,k) = 0.0
      yplus_k1(i,j) = 25.0
      yplus_km(i,j) = 25.0
      yplus_i1(j,k) = 25.0
      yplus_im(j,k) = 25.0
      sparevar(i,j,k) = 1.0
 3500 continue
c
      write(6,*)  ' done initial guess of all flow velocities etc '
c
c**********************************************************************************
c**********************************************************************************
c      read in from restart file if "if_restart" = 1 to overwrite initial guess
c
      if(if_restart.eq.0) go to 3700
c
c******************************************************************************
      write(6,*)
      write(6,*)   ' reading in the restart file from unit 7.'
c
      read(7) nstep
      read(7)(((ro(i,j,k),i=1,im),j=1,jm),k=1,km)
      read(7)(((rovx(i,j,k),i=1,im),j=1,jm),k=1,km)
      read(7)(((rovr(i,j,k),i=1,im),j=1,jm),k=1,km)
      read(7)(((rovt(i,j,k),i=1,im),j=1,jm),k=1,km)
      read(7)(((roe(i,j,k),i=1,im),j=1,jm),k=1,km)
      read(7)(((rosub(i,j,k),i=1,im),j=1,jm),k=1,km) 
      read(7)(((trans_dyn_vis(i,j,k),i=1,im),j=1,jm),k=1,km)
c
      rewind(7)
c
      write(6,*) ' restart file read in ok '
      write(6,*)
c**********************************************************************************
c**********************************************************************************
c   setting the secondary variables from the restart file of primary variables.
c
      do 3901 k=1,km
      do 3900 j=1,jm
      do 3900 i=1,im
      vx(i,j,k)    = rovx(i,j,k)/ro(i,j,k)
      vr(i,j,k)    = rovr(i,j,k)/ro(i,j,k)
      vt(i,j,k)    = rovt(i,j,k)/ro(i,j,k)
      rorvt(i,j,k) = rovt(i,j,k)*r(j,k)
      eke = 0.5*(vx(i,j,k)*vx(i,j,k) + vr(i,j,k)*vr(i,j,k) 
     &   +       vt(i,j,k)*vt(i,j,k))
c   the value of roe sent to the plot file was the true value
c    change back to the artificial value "rosub*e" .
      if(itimst.eq.5.or.itimst.eq.6)
     &         roe(i,j,k) = roe(i,j,k)*rosub(i,j,k)/ro(i,j,k)
c
      if(ifgas.eq.0) then
           tstatic   = (roe(i,j,k)/ro(i,j,k) - eke)/cv
           t_static(i,j,k) = tstatic
           ho(i,j,k) = cp*tstatic + eke
      else
           estat     = roe(i,j,k)/ro(i,j,k)  - eke
           tstatic   = tfrome(estat,tref,eref,et1,et2,et3,et4)
           t_static(i,j,k) = tstatic
           hstat     = hfromt(tstatic,tref,href,cp1,cp2,cp3)
           ho(i,j,k) = hstat + eke
      end if
c
      if(itimst.lt.5) then
      rosub(i,j,k) = ro(i,j,k)
      p(i,j,k)     = ro(i,j,k)*rgas*tstatic
      end if
c
      if(itimst.ge.5) p(i,j,k) = po_ref - dp_dro*(ro_ref - rosub(i,j,k))
c
      peff(i,j,k)  = p(i,j,k)

 3900 continue
c
      pdown(k) = 0.0
      do 3902 i=1,imm1
      pdown(k) = pdown(k) + fp(i)*0.5*(p(i,jm,k)+p(i+1,jm,k))
 3902 continue
c
 3901 continue
c
      if(ipout.ge.1)  go to 3700
      if(ipout.eq.0)  pdiff =   pd(1)  - pdown(1)
      if(ipout.eq.-1) pdiff =   pd(km) - pdown(km) 
c         
      do 3903 k=1,km
      pd(k) = pdown(k) + pdiff
 3903 continue
c
      write(6,*) ' all flow properties set up from the restart file.'
c
c******************************************************************************
c   end of setting up the flow from the restart file.
 3700 continue
c
c******************************************************************************
c******************************************************************************
c        set the value of step(i,j,k)= timestep/volume for each element.
c        step(j,k) is the main time step multiplied by 1/volume
c        if itimst = 2  fixed non uniform time steps are taken with dt
c        proportional to ds. if itimst =3 the the time step is
c        regularly updated by subroutine stepup.
c
      write(6,*) ' setting the time step. '
c
      if(itimst.lt.5) then
           if(ifgas.eq.0) then
                vsound = sqrt(rgas*ga*to1(kmid))
           else
                cpnow = cp1 + cp2*(to1(kmid)-tref) 
     &                + cp3*(to1(kmid)-tref)*(to1(kmid)-tref)
                gamnow = cpnow/(cpnow-rgas)
                vsound = sqrt(gamnow*rgas*to1(kmid))
            end if
      end if
c
c     set the time step based on the local grid dimension and a constant speed of sound.
c     this is modified to the true speed of sound in subroutine "stepup" . 
c
      do 6000 j=2,jm
      do 6000 k=1,kmm1
      do 6000 i=1,imm1
      smin = dsmin
      attot = sqrt(abt(j,k)*abt(j,k)+abx(i,j,k)*abx(i,j,k)
     &           + abr(i,j,k)*abr(i,j,k))
      perpt = vol(i,j,k)/attot
      perpq = vol(i,j,k)/aqtot(i,j,k)
      astot = sqrt(asx(i,j,k)*asx(i,j,k)+ asr(i,j,k)*asr(i,j,k))
      perps = vol(i,j,k)/astot
      smin  = amin1(perpt,perpq,perps)
      step(i,j,k)  = smin*cfl/vol(i,j,k)/vsound
      bstep(i,j,k) = step(i,j,k)
      rstep(i,j,k) = 1.0
 6000 continue
c
c     set the time step to an average at the mixing plane
c
      do 101 nr=1,nrwsm1
      j = jmix(nr)
      do 201 k=1,kmm1
      do 201 i=1,imm1
      step(i,j+1,k)  = 0.5*(step(i,j,k)+step(i,j+2,k))
      bstep(i,j+1,k) = step(i,j+1,k)
      rstep(i,j+1,k) = 0.0
  201 continue
  101 continue
c
c     call stepup to set the locally varying time steps if itimst >= 3 .
c
      if(itimst.ge.3) then
           call stepup(1.0)
      endif
c
c*******************************************************************************
c*******************************************************************************
c
      write(6,*)  ' setting up the multigrid arrays '
c
c         starting to set up arrays for use in multigrid
c
      if(jr.eq.1.and.ir.eq.1.and.kr.eq.1) go to 7100
c
      do 7000 k=1,kmm1
      kb2(k) =  1 + (k-1)/krbb
 7000 kb1(k) =  1 + (k-1)/kr
      nkb1 = kb1(kmm1)
      nkb2 = kb2(kmm1)
c
      do 7002 i=1,imm1
      ib2(i) = 1 + (i-1)/irbb
 7002 ib1(i) = 1 + (i-1)/ir
      nib1 = ib1(imm1)
      nib2 = ib2(imm1)
c
c     sort the j values into multigrid blocks making sure that
c     no blocks overlap a mixing plane
c
      jstt = 2
      jb1(1) = 1
      jb1(2) = 1
      do 7003 j=3,jm
      mark = 0
      if(indmix(j).eq.1.or.j.eq.jm)       mark = 1
      if((mod((j-jstt),jr).eq.0).and.(mark.eq.0)) then
           jb1(j) = jb1(j-1) + 1
      else
           jb1(j) = jb1(j-1)
      endif
      if(indmix(j-1).eq.1) jb1(j) = 0
      if(indmix(j-2).eq.1) then
           jb1(j) = jb1(j-2) + 1
           jstt   = j
      endif
 7003 continue
c
      jstt = 2
      jb2(1) = 1
      jb2(2) = 1
      do 7004 j=3,jm
      jp1 = j+1
      if(j.eq.jm) jp1 = jm
      mark = 0
      if(indmix(j).eq.1.or.j.eq.jm)       mark =1
      if(indmix(jp1).eq.1.or.(j+1).eq.jm) mark =1
      if((mod((j-jstt),jrbb).eq.0).and.(mark.eq.0)) then
           jb2(j) = jb2(j-1) + 1
      else
           jb2(j) = jb2(j-1)
      endif
      if(indmix(j-1).eq.1) jb2(j) = 0
      if(indmix(j-2).eq.1) then
           jb2(j)=jb2(j-2) + 1
           jstt = j
      endif
 7004 continue
c
      njb1 = jb1(jm)
      njb2 = jb2(jm)
c
c     check that the multigrid dimensions are not too large
c
      if(nib1.gt.ig1) write(6,*) 'stopping because nib1 too large.'
      if(nkb1.gt.kg1) write(6,*) 'stopping because nkb1 too large.'
      if(njb1.gt.jg1) write(6,*) 'stopping because njb1 too large '
      if(nib1.gt.ig1.or.nkb1.gt.kg1.or.njb1.gt.jg1) stop
c
c
      if(nib2.gt.ig2) write(6,*) 'stopping because nib2 too large.'
      if(nkb2.gt.kg2) write(6,*) 'stopping because nkb2 too large.'
      if(njb2.gt.jg2) write(6,*) 'stopping because njb2 too large '
      if(nib2.gt.ig2.or.nkb2.gt.kg2.or.njb2.gt.jg2) stop

      jstart1(1)  = 2
      jstart2(1)  = 2
      do 7005 j=2,jm
      if(jb1(j).eq.0) jb1(j) = njb1+1
      if(jb2(j).eq.0) jb2(j) = njb2+1
      if(jb1(j).ne.jb1(j-1)) jstart1(jb1(j)) = j
      jend1(jb1(j)) = j
      if(jb2(j).ne.jb2(j-1)) jstart2(jb2(j)) = j
      jend2(jb2(j)) = j
 7005 continue
c
      jstart1(njb1+1) = jm
      jend1(njb1+1)   = jm
      jstart2(njb2+1) = jm
      jend2(njb2+1)   = jm
c
c
 7100 continue
c
c******************************************************************************
c        set up block sizes for the supergrid jsblk(j) is superblock
c                           indicator
c
      write(6,*)  ' setting up the super block arrays',
     &   ' there are only 4 superblocks per blade row.'
      nsb=1
      do 199 j=1,jm
      jsblk(j) = nsb
      if(indle(j).eq.1) nsb=nsb+1
      if(indmid(j).eq.1)nsb=nsb+1
      if(indte(j).eq.1) nsb=nsb+1
  199 continue
c      write(6,7102)(jsblk(j),j=1,jm)
c 7102 format( ' super block index= ',20i5)
c
      nsblk = jsblk(jm)
c
      if(nsblk.gt.jg3) then
      write(6,*)    ' stopping because the total number of super blocks,
     & nsblk, is too large.'
      stop
      endif
c
c******************************************************************************
c     set the time steps for the multigrid blocks
c     first for the level 1 blocks
c
      do 7110  i1 = 1,nib1
      istart = (i1-1)*ir + 1
      iend = istart + ir -1
      if(iend.gt.imm1) iend = imm1
      do 7110  j1 = 1,njb1
      jstrt = jstart1(j1)
      jendd = jend1(j1)
      do 7110  k1 = 1,nkb1
      kstart = (k1-1)*kr + 1
      kend = kstart + kr -1
      if(kend.gt.kmm1) kend = kmm1
      perpt = 0.0
      perpq = 0.0
      perps = 0.0
      volb  = 0.0
      do 7120 i = istart,iend
      do 7120 k = kstart,kend
      do 7120 j = jstrt,jendd
      attot = sqrt(abt(j,k)*abt(j,k)+abx(i,j,k)*abx(i,j,k)
     & + abr(i,j,k)*abr(i,j,k))
      perpt = perpt +  vol(i,j,k)/attot
      perpq = perpq + vol(i,j,k)/aqtot(i,j,k)
      astot = sqrt(asx(i,j,k)*asx(i,j,k)+ asr(i,j,k)*asr(i,j,k))
      perps = perps + vol(i,j,k)/astot
      volb  = volb + vol(i,j,k)
 7120 continue
      perps = perps/(ir*jr)
      perpt = perpt/(jr*kr)
      perpq = perpq/(ir*kr)
      perpmin = amin1(perps,perpt,perpq)
      step1(i1,j1,k1)     = cfl*fblk1*perpmin/vsound/volb
      step1(i1,njb1+1,k1) = 0.0
 7110 continue
c********************************************************************************
c     now for the level  2  blocks
c
      do 8110  i2 = 1,nib2
      istart = (i2-1)*irbb + 1
      iend   = istart + irbb - 1
      if(iend.gt.imm1) iend = imm1
      do 8110  j2 = 1,njb2
      jstrt = jstart2(j2)
      jendd = jend2(j2)
      do 8110  k2 = 1,nkb2
      kstart = (k2-1)*krbb + 1
      kend = kstart + krbb - 1
      if(kend.gt.kmm1) kend = kmm1
      perpt = 0.0
      perpq = 0.0
      perps = 0.0
      volb  = 0.0
      do 8120 i = istart,iend
      do 8120 k = kstart,kend
      do 8120 j = jstrt,jendd
      attot = sqrt(abt(j,k)*abt(j,k)+abx(i,j,k)*abx(i,j,k)
     & + abr(i,j,k)*abr(i,j,k))
      perpt = perpt +  vol(i,j,k)/attot
      perpq = perpq + vol(i,j,k)/aqtot(i,j,k)
      astot = sqrt(asx(i,j,k)*asx(i,j,k)+ asr(i,j,k)*asr(i,j,k))
      perps = perps + vol(i,j,k)/astot
      volb  = volb + vol(i,j,k)
 8120 continue
      perps = perps/(irbb*jrbb)
      perpt = perpt/(jrbb*krbb)
      perpq = perpq/(irbb*krbb)
      perpmin = amin1(perps,perpt,perpq)
      step2(i2,j2,k2)     = cfl*fblk2*perpmin/vsound/volb
      step2(i2,njb2+1,k2) = 0.0
 8110 continue
c*******************************************************************************
c     now for the super blocks
c
      perpq = 0.0
      volb  = 0.0
      do 9010 j = 2,jm
      jsb = jsblk(j)
      do 9020 i=1,imm1
      do 9020 k=1,kmm1
      perpq = perpq + vol(i,j,k)/aqtot(i,j,k)
      volb  = volb  + vol(i,j,k)
 9020 continue
c
      if(j.eq.jm) then
      perpavg      = perpq/(imm1*kmm1)
      stepsbk(jsb) = cfl*fblk3*perpavg/vsound/volb
      go to 9010
      endif
      if(jsblk(j+1).ne.jsb) then
      perpavg = perpq/(imm1*kmm1)
      stepsbk(jsb) = cfl*fblk3*perpavg/vsound/volb
      perpq = 0.0
      volb  = 0.0
      endif
 9010 continue
c
c     end of setting the multigrid time steps
c
c******************************************************************************
c*****************************************************************************
c     initialise the source terms to zero
c     note that this section has been moved for the sa model because it needs step(i,j,k)
c     to be set.
c
      do 3710 k=1,km
      do 3710 j=1,jm
c  tflow
      pblade(j,k)   = 0.0
      rovar_m1(j,k) = 0.0
      bforce_x(j,k) = 0.0
      bforce_r(j,k) = 0.0
      bforce_t(j,k) = 0.0
      bforce_q(j,k) = 0.0
c end tflow
      do 3710 i=1,im
      xforce(i,j,k)   = 0.0
      tforce(i,j,k)   = 0.0
      rforce(i,j,k)   = 0.0
      qsource(i,j,k)  = 0.0
      source(i,j,k)   = 0.0
      t_source(i,j,k) = 0.0
      sgen(i,j,k)     = 0.0
 3710 continue
c     
c
      if(nneg.gt.0) then
           write(6,*)'*******************warning***********************'
           write(6,*)   nneg,' negative volumes found ' 
           write(6,*) ' this is very likely to cause failure. '
           write(6,*) ' do you want to continue despite this ?'
           write(6,*) ' answer   y  or  n '
c           answ = 'n'
           read(1,*)    answ
           if(answ.eq.'n'.or.answ.eq.'n')  stop
      end if
      close(1)
c*******************************************************************************
c   jdd added 30/9/10. to set a limit on the vorticity
c   modified 18/9/2018 to base it on the relative velocities
      dr1 = r(1,km)  - r(1,1)
      dx1 = x(1,km)  - x(1,1)
      drm = r(jm,km) - r(jm,1)
      dxm = x(jm,km) - x(jm,1)

      vxin = vx(imid,1,kmid)
      vrin = vr(imid,1,kmid)
      wtin = vt(imid,1,jmid) - wrad(1)*r(1,kmid)
      win  = sqrt(vxin*vxin + vrin*vrin + wtin*wtin)

      vxout = vx(imid,jm,kmid)
      vrout = vr(imid,jm,kmid)
      wtout = vt(imid,jm,jmid) - wrad(jm)*r(jm,kmid)
      wout  = sqrt(vxout*vxout + vrout*vrout + wtout*wtout)
     
      spanin = sqrt(dx1*dx1 + dr1*dr1)
      spanout= sqrt(dxm*dxm + drm*drm)
  
      vortin   = win/(spanin*fr(1))
      vortout  = wout/(spanout*fr(1))
      vort_max = amax1(vortin,vortout)
c
      write(6,*)
      write(6,9021)  vort_max
 9021 format(' vorticity limit, vort_max = ',e15.3)
      write(6,*)
c
c  end of jdd addition
c******************************************************************************
c
c call set_xlength to calculate the wall distances and set the mixing lengths
c
      write(6,*)    ' calling  set_xlength to set the wall distances and 
     & mixing lengths'
c
       call set_xlength
c
      write(6,*) ' called  set_xlength '
c
c******************************************************************************
c******************************************************************************
c      call loss routines to initialise the body force terms
c
      write(6,*)  ' calling the loss routines from setup, ilos = ',ilos
c
      temp_rf_vis  = rf_vis
      rf_vis   = 1.0
      fmixup   = 1.0
      nstep    = 1 
      if(ilos.eq.10)                  call loss
      if(ilos.ge.100.and.ilos.lt.200) call new_loss
      if(ilos.ge.200)                 call spal_loss
      rf_vis   = temp_rf_vis
c
c******************************************************************************
c******************************************************************************
c  initialise the shroud flows, bleed flows and cooling flows to zero.
      do 9030 j = 1,jm
      shrdflow(j) = 0.0
      sumbleed(j) = 0.0
      sumcwl(j)   = 0.0
 9030 continue
c
c******************************************************************************
c tflow
c******************************************************************************
c
      if(im.eq.2) then
c
c  find the centre line grid angle alpha_cent and smooth its streamwise variation.
      do 9039 k=1,kmm1
      s_dist(1) = 0.0
      do 9038  j=2,jm
           abxavg     = 0.5*( abx(1,j,k)  + abx(2,j,k) )
           abravg     = 0.5*( abr(1,j,k)  + abr(2,j,k) )
           vecx       = aqx(1,j,k)/aqtot(1,j,k)
           vecr       = aqr(1,j,k)/aqtot(1,j,k)
           abmeravg   = abxavg*vecx  + abravg*vecr
           cent_angl(j) = atan(abmeravg/abt(j,k))
           xdif  = x(j,k) - x(j-1,k)
           rdif  = r(j,k) - r(j-1,k)
           s_dist(j) = s_dist(j-1) + sqrt(xdif*xdif + rdif*rdif)
 9038 continue
c
c     set the centre line angles upstream and downstream of a blade row.
      do 9036 j=1,jm
      nrw   = nrow(j)
      jlee  = jle(nrw)
      jtee  = jte(nrw) 
c      avgle = (cent_angl(jlee)+cent_angl(jlee+1)+cent_angl(jlee+2))/3
c      avgle  = 2.0*cent_angl(jlee+1) - cent_angl(jlee+2)
c      avgte = (cent_angl(jtee)+cent_angl(jtee-1)+cent_angl(jtee-2))/3
c      avgte  = 2.0*cent_angl(jtee-1) - cent_angl(jtee-2)
c
c   jdd changed this to setting the angles to equal the le and te blade angles.
c   may 2017.
      if(j.le.jlee) cent_angl(j) = cent_angl(j+1)
      if(j.ge.jtee) cent_angl(j) = cent_angl(j-1)
 9036 continue
c
c    smooth the centre line angle variation.
      do 9035 nrw = 1,nrows
      jlee  = jle(nrw)
      jtee  = jte(nrw)      
           call smooth(jlee,jtee,4,0.25,s_dist,cent_angl)
 9035 continue
c
      do 9037 j=1,jm
      alpha_cent(j,k) = cent_angl(j)
 9037 continue
c           
 9039 continue
c
c   interpolate to find the deviation angle or blade exit flow angle depending
c   on  angl_typ(nr) .
      do 9040 nr = 1,nrows
c
      jtee = jte(nr)
      qspan(1) = 0.0
      do 9041 k=2,km
           xdif = x(jtee,k) - x(jtee,k-1)
           rdif = r(jtee,k) - r(jtee,k-1)
           qspan(k) = qspan(k-1) + sqrt(xdif*xdif + rdif*rdif)
 9041 continue
           qtot = qspan(km)
      do 9042 k=1,km
           qspan(k) = qspan(k)/qtot
 9042 continue
      do 9043 k=1,nangles(nr)
           fspan(k)  = frac_span(nr,k)
           exang(k)  = exit_angl(nr,k) 
 9043 continue
      do 9044 k=1,km
      call intp(nangles(nr),fspan,exang,qspan(k),angl_out)
      if(angl_typ(nr).eq.'a') exit_angl(nr,k) = angl_out*degrad
      if(angl_typ(nr).eq.'d') devn_angl(nr,k) = angl_out*degrad
 9044 continue
c
      write(6,*)
      write(6,*) ' row number', nr
           write(6,9047) (alpha_cent(jtee,k)*raddeg,k=1,km)
      if(angl_typ(nr).eq.'a')
     &     write(6,9045) (exit_angl(nr,k)*raddeg,k=1,km)
      if(angl_typ(nr).eq.'d')
     &     write(6,9046) (devn_angl(nr,k)*raddeg,k=1,km)
 9045 format(' blade exit flow angle in degrees ',/,(10f10.3))
 9046 format(' exit flow angle deviation from grid angle, in degrees.',
     & /,(10f10.3))
 9047 format(' grid centre line angle in degrees at blade exit',
     & /,(10f10.3))
c
c   convert the exit angle to deviation angle if  "angl_typ" = a .
c
      if(angl_typ(nr).eq.'a') then 
           do 9049 k=1,kmm1
           grid_angl       = alpha_cent(jtee,k)
           devn_angl(nr,k) = grid_angl - exit_angl(nr,k)
 9049      continue
      end if      
c
 9040 continue
c
c   end of  "if im = 2" loop
      end if
c******************************************************************************
c  end tflow
c******************************************************************************
c
c************if nout(1)=0 write out and plot out the initial guess of the flow field.
      if(nout(1).eq.0) then
           write(6,*) ' writing a plot file of the initial guess '
           call output
      end if
c
c
      write(6,*)
      write(6,*)'******************************************************'
      write(6,*) ' leaving subroutine setup,  ok so far '
      write(6,*)'******************************************************'
      write(6,*)
c
      return
      end
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c
      subroutine new_intpol
c       ====================
c
c       this routine interpolates values for x  ,  r,  rt_upp  and
c       rt_thick for the blade sections on the final grid. 
c       the grid is adjusted to allow for the tip gap and number of pointrs in the gap.
c
c       ====================
c
      include  'commall-open-19.2'
c
      dimension   xint(maxki),yint(maxki),rtup_int(maxki),
     &            rttk_int(maxki),frmod(maxki),qodist(maxki)
c
      double precision sumfr,tipgap,tipspace,fac
c
c
      do 1 k=1,nsecs_in
      smerid(1,k) = 0.0
      do 2 j=2,jm
      xd = xsurf(j,k)  - xsurf(j-1,k)
      rd = rsurf(j,k)  - rsurf(j-1,k)
    2 smerid(j,k) = smerid(j-1,k) + sqrt(xd*xd+rd*rd)
    1 continue
c
c******************************************************************************
c******************************************************************************
c
      do 1000 j=1,jm
c
      nr = nrow(j)
      j1 = jstart(nr)
      jl = jle(nr)
      jt = jte(nr)
      je = jmix(nr)
c
c    set the default that   frmod(k)   = fr(k) if no tip gap .
c
      do 55 k=1,kmm1
   55 frmod(k) = fr(k)
c
c   do not modify the grid and jump to 75 if there is no tip gap.
c
      if(ktips(nr).le.0) go to 75
c
c******************************************************************************
c   set frmod(k) to modify the grid if there is a tip gap.
c   the gap is varied linearly from fractip1 at the le to fractip2  at the te.
c
          ktip   = nsecs_in
          if(ktips(nr).eq.1) ktip =1
          fac = 1.0
          if(j.ge.j1.and.j.lt.jl) then
          fac = (smerid(j,ktip)-smerid(j1,ktip))
     &         /(smerid(jl,ktip)-smerid(j1,ktip)) 
          tipgap = fractip1(nr)
          end if

          if(j.ge.jl.and.j.le.jt) then
               fac  = 1.0
               frac = (smerid(j,ktip) - smerid(jl,ktip))
     &               /(smerid(jt,ktip)- smerid(jl,ktip))
               tipgap = fractip1(nr) + frac*(fractip2(nr)-fractip1(nr)) 
          end if

          if(j.gt.jt.and.j.le.je) then
          fac = (smerid(je,ktip)-smerid(j,ktip))
     &    /(smerid(je,ktip)- smerid(jt,ktip))
          tipgap = fractip2(nr)
          end if
c
      sumfr = 0.0
      ncell    = ktipe(nr) - ktips(nr)
      tipspace = tipgap/ncell
      do 50 k=ktips(nr),ktipe(nr)-1
   50 sumfr = sumfr + fr(k)
      do 60 k=ktips(nr),ktipe(nr)-1
   60 frmod(k) =  tipspace*fac + (1.-fac)*fr(k)
      do 70 k= ks1(nr),ks2(nr)
   70 frmod(k) = fr(k)*((1.-tipgap)/(1.-sumfr)*fac + (1.-fac) )
c
c   end of setting frmod(k) for the tip gap .
c****************************************************************************
c****************************************************************************
c
   75 continue
c
c     now interpolate in the nsecs_in uniformly spaced stream surfaces to obtain
c     the coordinates on   km   stream surfaces for the final grid.
c
      qodist(1)   = 0.0
      do 10 k=1,nsecs_in
      xint(k)     = xsurf(j,k)
      yint(k)     = rsurf(j,k)
      rtup_int(k) = rt_upp(j,k)
      rttk_int(k) = rt_thick(j,k)
      if (k.gt.1) then
         rd = rsurf(j,k)  -  rsurf(j,k-1)
         xd = xsurf(j,k)  -  xsurf(j,k-1)
	 qodist(k) = qodist(k-1) + sqrt(xd*xd+rd*rd)
      end if
   10 continue
c
      xarg   = 0.0
      qospan = qodist(nsecs_in)
c
      do 30 k=1,km
      call intp(nsecs_in,qodist,xint,xarg,x(j,k))
      call intp(nsecs_in,qodist,yint,xarg,r(j,k))
      call intp(nsecs_in,qodist,rtup_int,xarg,rt_upp(j,k))
      call intp(nsecs_in,qodist,rttk_int,xarg,rt_thick(j,k))
      xarg  = xarg + frmod(k)*qospan
   30 continue
c  
 1000 continue
c
      return
      end
c
c***********************************************************************
c
      subroutine old_intpol
c       ====================
c
c       this routine interpolates values for rt_upp
c       rt_thick, x  and r for the required cross-sections
c
c       ====================
c
      include  'commall-open-19.2'
c
      dimension   xint(maxki),yint(maxki),rtup_int(maxki),
     &            rttk_int(maxki),frmod(maxki),qodist(maxki)
c
      double precision sumfr,tipgap,tipspace,fac
c
      qodist(1) = 0.0
c
      do 1 k=1,nosect
      smerid(1,k) = 0.0
      do 2 j=2,jm
      xd = xsurf(j,k)  - xsurf(j-1,k)
      rd = rsurf(j,k)  - rsurf(j-1,k)
    2 smerid(j,k) = smerid(j-1,k) + sqrt(xd*xd+rd*rd)
    1 continue
c
c******************************************************************************
c
      do 1000 j=1,jm
      nr = nrow(j)
      j1 = jstart(nr)
      jl = jle(nr)
      jt = jte(nr)
      je = jmix(nr)
c
c   do not modify the grid and jump to 75 if no tip gap.
c
      do 55 k=1,kmm1
   55 frmod(k) = fr(k)
c
      if(ktips(nr).le.0) go to 75
c
c******************************************************************************
c   modify the grid if there is a tip gap.
c   the gap is varied linearly from fractip1 at the le to fractip2  at the te.
c
          ktip   = nosect
          if(ktips(nr).eq.1) ktip =1
          fac = 1.0
          if(j.ge.j1.and.j.lt.jl) then
          fac = (smerid(j,ktip) -smerid(j1,ktip))
     &         /(smerid(jl,ktip)-smerid(j1,ktip)) 
          tipgap = fractip1(nr)
          end if

          if(j.ge.jl.and.j.le.jt) then
               fac  = 1.0
               frac = (smerid(j,ktip) - smerid(jl,ktip))
     &               /(smerid(jt,ktip)- smerid(jl,ktip))
               tipgap = fractip1(nr) + frac*(fractip2(nr)-fractip1(nr)) 
          end if

          if(j.gt.jt.and.j.le.je) then
          fac = (smerid(je,ktip)-smerid(j,ktip))
     &    /(smerid(je,ktip)- smerid(jt,ktip))
          tipgap = fractip2(nr)
          end if
c
c      write(6,*) 'j,j1,jl,jt,je, fac, tipgap, ratgap ', j,j1,jl,jt,je,
c     &            fac ,tipgap, ratgap
c
c      change fr(k) to frmod(k) to slightly adjust the blade tip 
c      position so that the tip gap is fractip.
c      first set frmod(k).
c
      sumfr = 0.0
      ncell    = ktipe(nr) - ktips(nr)
      tipspace = tipgap/ncell
      do 50 k=ktips(nr),ktipe(nr)-1
   50 sumfr = sumfr + fr(k)
      do 60 k=ktips(nr),ktipe(nr)-1
   60 frmod(k) =  tipspace*fac + (1.-fac)*fr(k)
      do 70 k= ks1(nr),ks2(nr)
   70 frmod(k) = fr(k)*((1.-tipgap)/(1.-sumfr)*fac + (1.-fac) )
c
c****************************************************************************
c
   75 continue
c
      do 10 k=1,nosect
      xint(k)     = xsurf(j,k)
      yint(k)     = rsurf(j,k)
      rtup_int(k) = rt_upp(j,k)
      rttk_int(k) = rt_thick(j,k)
      if (k.gt.1) then
         rd = rsurf(j,k)  - rsurf(j,k-1)
         xd = xsurf(j,k)  -  xsurf(j,k-1)
	 qodist(k) = qodist(k-1) + sqrt(xd*xd+rd*rd)
      end if
   10 continue
c
      dref = 0.001*qodist(nosect)
      dref = dref*dref
c
c     find the point nearest to the hub, l1
c
      l1 = 1
      flag1 = 1.
      do 20 l=1,nosect-1
      if ((rsurf(j,l)-r(j,1))*(rsurf(j,l+1)-r(j,1)).le.dref.and.
     *   (xsurf(j,l)-x(j,1))*(xsurf(j,l+1)-x(j,1)).le.dref) then
         l1=l
         goto 21
      end if
   20 continue
      flag1 = -1.
      if(insurf.ne.2) write(6,111) j
  111 format(//, ' warning !!!! the first blade section is outboard of
     &the hub',/,' the extrapolation may be very inaccurate ','j=',i3)
c
c     find the point nearest to the casing, lm
c
   21 continue
c
      lm = nosect
      flag2 = 1.
      do 22 l=nosect,2,-1
      if ((rsurf(j,l)-r(j,km))*(rsurf(j,l-1)-r(j,km)).le.dref.and.
     *   (xsurf(j,l)-x(j,km))*(xsurf(j,l-1)-x(j,km)).le.dref) then
         lm=l
         goto 25
      end if
   22 continue
c
      flag2 = -1.
      if(insurf.ne.2) write(6,112) j
  112 format(//,' warning !!!! the last blade section is inboard of the
     &casing ',/,' the extrapolation may be very inaccurate ','j= ',i3)
c
   25 continue
c
      rd = rsurf(j,l1) - r(j,1)
      xd = xsurf(j,l1)   - x(j,1)
      qdist1 = sqrt(rd*rd+xd*xd)*flag1
      rd = rsurf(j,lm) - r(j,km)
      xd = xsurf(j,lm)   - x(j,km)
      qdistm = sqrt(rd*rd+xd*xd)*flag2
      qospan = qodist(lm) - qodist(l1) - qdist1 - qdistm
c
      xarg   = qodist(l1) + qdist1
c
c      write(6,*)
c      write(6,*) 'k,   x(j,k),   r(j,k),   rt_upp(j,k),   rt_thick(j,k)'
      do 30 k=1,km
      if (k.gt.1) xarg = xarg + qospan*frmod(k-1)
      if(insurf.ne.2.and.(xarg.lt.qodist(1).or.xarg.gt.qodist(nosect)))
     &      then
c
      write(6,*) ' warning linear extrapolation at j= ',j,'k= ',k
c
      call linint(nosect,qodist,xint,xarg,x(j,k))
      call linint(nosect,qodist,yint,xarg,r(j,k))
      call linint(nosect,qodist,rtup_int,xarg,rt_upp(j,k))
      call linint(nosect,qodist,rttk_int,xarg,rt_thick(j,k))
      else
      call intp(nosect,qodist,xint,xarg,x(j,k))
      call intp(nosect,qodist,yint,xarg,r(j,k))
      call intp(nosect,qodist,rtup_int,xarg,rt_upp(j,k))
      call intp(nosect,qodist,rttk_int,xarg,rt_thick(j,k))
      endif
c
c      if(j.eq.150) write(6,66) k,x(j,k),r(j,k),rt_upp(j,k),rt_thick(j,k)
c   66 format(i5,4f10.5)
c
   30 continue
c
      do 40 k=1,km
      xsurf(j,k)   = x(j,k)
      rsurf(j,k)   = r(j,k)
   40 continue
 1000 continue
c
      return
      end
c
c*****************************************************************************
c*****************************************************************************
c
      subroutine intp(n,xn,yn,x,y)
c
c      this subroutine interpolates in the given table of yn as a
c      function of xn to find the value of y at the input value
c      of x.
c
      dimension xn(n),yn(n)
      span=xn(n)-xn(1)
      y=0.
      l=1
      nm=n
      if(n.lt.4) go to 8
      nm=4
    4 if(span.gt.0.0.and.x.lt.xn(l)) go to 5
      if(span.lt.0.0.and.x.gt.xn(l)) go to 5
      if(l.eq.n) go to 3
      l=l+1
      go to 4
    5 if(l.gt.2) go to 6
      l=1
      go to 8
    6 if(l.ne.n) go to 7
    3 l=n-3
      go to 8
    7 l=l-2
    8 do 11 l1=1,nm
      co=1
      do 10 l2=1,nm
      if(l1.eq.l2) go to 9
      temp=(x-xn(l+l2-1))/(xn(l+l1-1)-xn(l+l2-1))
      go to 10
    9 temp=1
   10 co=co*temp
   11 y=y+co*yn(l+l1-1)
      return
      end
c
c*************************************************************************************
c
       subroutine linint(npoints,x,y,xarg,yans)
c
c      this subroutine interpolates in the given table of yn as a
c      function of x to find the value of y at the input value
c      of x = xarg.
c
c      this version uses linear interpolation to avoid any possible
c      problems with overshoots or undershoots.
c
      dimension x(npoints),y(npoints)
c
      if (x(1).gt.xarg) then
      yans = y(1) + (xarg-x(1))*(y(2)-y(1))/(x(2)-x(1))
      else
      n=2
   10 continue
      if(x(n).gt.xarg) go to 20
      n=n+1
      if(n.gt.npoints) go to 30
      go to 10
   20 yans = y(n) + (xarg-x(n))*(y(n-1)-y(n))/(x(n-1)-x(n))
      go to 40
   30 yans = y(npoints) + (xarg-x(npoints))*(y(npoints)-y(npoints-1))/
     & (x(npoints)-x(npoints-1))
   40 continue
c  
      endif
c
      return
      end
c
c******************************************************************************
c 
      subroutine sumflx(blade_flow,sumpo,sumto,sumrvt,sum_entpy,
     &                  sumtstat,sumpstat) 
c
c   this subroutine mass averages some flow quantities at every  "j"  station.
c
      include  'commall-open-19.2'
c
      dimension blade_flow(jd),sum_entpy(jd),sumpo(jd),sumto(jd),
     &          sumrvt(jd),sumtstat(jd),sumpstat(jd)
c
      do 100 j=1,jm
      blade_flow(j)  = 0.0
      sumpo(j)      = 0.0
      sumto(j)      = 0.0
      sumrvt(j)     = 0.0
      sum_entpy(j)  = 0.0
      sumtstat(j)   = 0.0
      sumpstat(j)   = 0.0
      do 110 k=1,kmm1
      do 120 i=1,imm1
      dflow        = -flowx(i,j,k)*nblade(j)
      blade_flow(j) = blade_flow(j) + dflow
      dflow4       = 0.25*dflow
      sumpo(j)     = sumpo(j) + dflow4*(temp2(i,j,k)+temp2(i+1,j,k)
     &             + temp2(i+1,j,k+1)+temp2(i,j,k+1))
      sumto(j)     = sumto(j) + dflow4*(temp1(i,j,k)+temp1(i+1,j,k)
     &             + temp1(i+1,j,k+1)+temp1(i,j,k+1))
      sumrvt(j)    = sumrvt(j)+ dflow4*(temp3(i,j,k)+temp3(i+1,j,k)
     &             + temp3(i+1,j,k+1)+temp3(i,j,k+1))
      sum_entpy(j) = sum_entpy(j) + dflow4*(temp4(i,j,k)+temp4(i+1,j,k)
     &             + temp4(i+1,j,k+1)+temp4(i,j,k+1))
      sumtstat(j)  = sumtstat(j) + dflow4*(store2(i,j,k)+store2(i+1,j,k)
     &             + store2(i+1,j,k+1)+store2(i,j,k+1))
      sumpstat(j)  = sumpstat(j)  + dflow4*(p(i,j,k)+p(i+1,j,k)
     &             + p(i+1,j,k+1) + p(i,j,k+1))
  120 continue
  110 continue
  100 continue
      return
      end
c
c******************************************************************************
c
      subroutine stepup(relax)
c
c          this subroutine changes the timestep in proportion to the local
c          mach number. the changes are relaxed by the factor  "relax" .
c
      include  'commall-open-19.2'
c
c
      tlim   = 0.1*to1(kmid)
      vlim   = 0.01*vm1(kmid)*vm1(kmid)
      relax1 = 1.0 - relax
c******************************************************************************
c     calculate the local speed of sound and set the factor  c/(v+c)
c
      do 1100 k=1,km
      do 1100 j=1,jm
      do 1100 i=1,im
      vsq = vx(i,j,k)*vx(i,j,k)+vt(i,j,k)*vt(i,j,k)+vr(i,j,k)*vr(i,j,k)
      hstat  = ho(i,j,k) - 0.5*vsq
c
      if(ifgas.eq.0) then
           tstatic = hstat*rcp
           gamnow  = ga
      else
           tstatic = tfromh(hstat,tref,href,ht1,ht2,ht3,ht4)
           cpnow   = cp1 + cp2*(tstatic-tref) 
     &             + cp3*(tstatic-tref)*(tstatic-tref)
           gamnow  = cpnow/(cpnow-rgas)
      end if
c
      if(tstatic.lt.tlim) tstatic = tlim
      wtrel     = wt(i,j,k)
      wsq       = vsq - (vt(i,j,k)*vt(i,j,k) - wtrel*wtrel)
      if(wsq.lt.vlim) wsq = vlim
c
      if(itimst.ge.5) then
           v_sonic = vsound
      else
           v_sonic = sqrt(gamnow*rgas*tstatic)
      end if
c
      vplusc    =  sqrt(wsq) + v_sonic
c
 1100 temp1(i,j,k) =  vsound/vplusc
c
c       temp1 is used as a temporary store for  vsound/(w + vsound)
c       where vsound is the stagnation speed of sound and w isthe relative velocity.
c
c******************************************************************************
c  use  bstep and the local velocities to set the time step/volume   step(i,j,k).
c
      do 10 i=1,imm1
      ip1=i+1
      do 10 k=1,kmm1
      kp1=k+1
      amachp =  temp1(i,1,k)   + temp1(ip1,1,k) +
     &          temp1(i,1,kp1) + temp1(ip1,1,kp1)
      do 10 j=  2,jm
      amachl =  temp1(i,j,k)   + temp1(i,j,kp1) +
     &          temp1(ip1,j,k) + temp1(ip1,j,kp1)
      avmach =  0.125*(amachl  + amachp)
      amachp =  amachl
      stepnew      =  bstep(i,j,k)*avmach
      step(i,j,k)  =  relax*stepnew  +  relax1*step(i,j,k)
      rstep(i,j,k) =  step(i,j,k)/bstep(i,j,k)
   10 continue
c
c******************************************************************************
c      set the time step at the mixing plane, ie at j = jmix + 1..
c
c    jdd modified this  march 2015. so step(jmix+1) = step(jmix+2).
c    which gives better results
c
      do 101 nr=1,nrwsm1
      j   = jmix(nr)
      jp1 = j+1
      jp2 = j+2
      do 201 k=1,kmm1
      do 201 i=1,imm1
      step(i,jp1,k)   = step(i,jp2,k)
      rstep(i,jp1,k)  = 0.0
  201 continue
  101 continue
c
c******************************************************************************
      return
      end
c
c**********************************************************************
c
      subroutine setflo
c
c      this subroutine forces the mass flow towards an input value in in_flow=3
c      or towards the average value if in_flow=2.
c      it does this by means of a body force which generates loss if in_flow=3
c      but in_flow =2 should give improved convergence with no loss generation.
c
      include  'commall-open-19.2'
c
      do 100 j=1,jm
      flow(j)=0.0
      do 100 k=1,kmm1
      do 100 i=1,imm1
      flow(j) = flow(j)-flowx(i,j,k)*nblade(j)/nblade(1)
  100 continue
c
      avflow = flowin/nblade(1)
      if(in_flow.eq.3) go to 300
      avflow=0.0
      do 200 j=1,jm
      avflow=avflow+flow(j)
  200 continue
      avflow=avflow/jm
  300 continue
c
      do 400 j=1,jm
      rovmer=sqrt(rovx(imid,j,kmid)*rovx(imid,j,kmid)+rovr(imid,j,kmid)
     &  *rovr(imid,j,kmid))
      delta = (avflow/flow(j) -1.0)*rovmer
      delta = delta*rflow
      do 400 i=1,im
      do 400 k=1,km
      rovmer = sqrt(rovx(i,j,k)*rovx(i,j,k)+rovr(i,j,k)*rovr(i,j,k))
c
c  added   3/9/90 may not be valid for radial flow machines
c
      if(rovx(i,j,k).lt.0.0)  rovmer = -rovmer
c
      rovx(i,j,k) = rovx(i,j,k) + rovx(i,j,k)/rovmer*delta
      rovr(i,j,k) = rovr(i,j,k) + rovr(i,j,k)/rovmer*delta
      ro_wt       = rovt(i,j,k) - ublade(j,k)*ro(i,j,k)
      rovt(i,j,k) = rovt(i,j,k) + ro_wt/rovmer*delta
c
  400 continue
c
      return
      end
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c
      subroutine loss
c
c            this subroutine computes a body force based on wall functions for the
c            surface shear stress and a mixing length model of eddy viscosity.
c
c            the body force is calculated by making a thin shear layer approximation
c            to the n-s equations.
c
c            the method uses wall functions to calculate the surface shear stresses
c            if "ypluswall" < 5 , and  evaluates the surface shear stresses by assuming
c            that the wall grid point is at "ypluswall" if ypluswall > 5.
c
c            the full energy equation including heat conduction and viscous work is
c            solved insteady of assuming that they cancel as in earlier versions.
c            all solid surfaces are assumed to be adiabatic - no heat flux.
c
c
      include  'commall-open-19.2'
c
c
      dimension vxavg(maxki),vravg(maxki),wtavg(maxki),wabs(maxki),
     &          roavg(maxki),xstres(maxki),ravg(maxki),rstres(maxki),
     &          tstres(maxki),area(maxki),wbound(maxki),wtb(maxki),
     &          vtavg(maxki),tavg(maxki),wvisc(maxki),qflow(maxki),
     &          tempp(jd)
c
c
c      calculate the viscosity over the first quarter of the steps.
c      then hold it constant for the remainder of the steps.
c
      if(nstep.eq.1.or.nstep.lt.nsteps_max/4) then
c
      if(reyno.gt.100.) then
           j1=jle(1)
           j2=jte(1)
           if(jle(1).gt.jm) j1=1
           if(jte(1).gt.jm) j2=jm
           xchord = smerid(j2,kmid)-smerid(j1,kmid)
           row2   = sqrt(rovx(imid,j2,kmid)*rovx(imid,j2,kmid)
     &            + rovr(imid,j2,kmid)*rovr(imid,j2,kmid)
     &            + rowt(imid,j2,kmid)*rowt(imid,j2,kmid))
           vislam  = xchord*row2/reyno
      end if
c
      if(reyno.gt.0.0.and.reyno.lt.99.99) then
            vislam = reyno/100000.
      end if
c
      if(reyno.lt.0.0) vislam = -reyno*1.0e-5
c
      tcond  = cp*vislam/prandtl
      ftcond = tcond/vislam
c
c   end of part only used for the first quarter of the steps.
      endif
c******************************************************************************
c    save the viscosity it is only used to calculate and write out the reynolds number. 
      do nrw =1,nrows
      viscosy(nrw)  = vislam
      end do
c******************************************************************************
c******************************************************************************
c******************************************************************************
c     do up to statement 25 only on the first call to the subroutine
c
      if(nstep.gt.nlos) go to 25
c
c     evaluate some parameters needed to calculate the viscous stresses on the
c     streamwise ( k = constant) surfaces.
c
      do 27 nrw = 1,nrows
c
      sumf = 0.0
      do 26 i  = 2,imm1
      sumf          = sumf+fp(i-1)
      xlim          = xllim_i1(nrw)*(1.0-sumf) + xllim_im(nrw)*sumf
      fpitch        = sumf*(1.-sumf)
      if(fpitch.gt.xlim) fpitch = xlim
      df            = fp(i)+fp(i-1)
      filam(nrw,i)  = fp(i)/df
      fiturb(nrw,i) = 0.16*fpitch*fpitch/(df*df)
      fiwake(nrw,i) = 0.16*xllim_dwn(nrw)*xllim_dwn(nrw)/(df*df)
      fiup(nrw,i)   = 0.16*xllim_up(nrw)*xllim_up(nrw)/(df*df)
   26 continue
      df            = fp(1) + fp(imm1)
      filam(nrw,1)  = fp(1)/df
      fiwake(nrw,1) = 0.16*xllim_dwn(nrw)*xllim_dwn(nrw)/(df*df)
      fiup(nrw,1)   = 0.16*xllim_up(nrw)*xllim_up(nrw)/(df*df)
      fiturb(nrw,1) = fiturb(nrw,2)
      filam(nrw,im) = fp(imm1)/df
      fiwake(nrw,im)= 0.16*xllim_dwn(nrw)*xllim_dwn(nrw)/(df*df)
      fiup(nrw,im)  = 0.16*xllim_up(nrw)*xllim_up(nrw)/(df*df)
      fiturb(nrw,im)= fiturb(nrw,imm1)
c
c  q3d
      if(km.eq.2) go to 16
c  end q3d
c
c     find the mixing length limits on the hub and casing
c
      avgspan = 0.0
      jmid    = 0.5*(jle(nrw)+jte(nrw))
      do 13 k=2,km
      xdif    = x(jmid,k)-x(jmid,k-1)
      rdif    = r(jmid,k)-r(jmid,k-1)
      avgspan = avgspan + sqrt(xdif*xdif+rdif*rdif)
   13 continue
      avgpit  = 2*3.1415926*r(jmid,kmid)/nblade(jmid)
      xlimh   = xllim_k1(nrw)*avgpit/avgspan
      xlimt   = xllim_km(nrw)*avgpit/avgspan
c
c     evaluate parameters needed to calculate the viscous stresses on the
c     bladewise (i = constant) surfaces.
c
      sumf = 0.0
      do 15 k=2,kmm1
      sumf     = sumf + fr(k-1)
      xlim     = xlimh*(1.0-sumf) + xlimt*sumf
      fcspan   = sumf*(1.-sumf)
      if(fcspan.gt.xlim) fcspan = xlim
      df       = fr(k) + fr(k-1)
      fklam(nrw,k)  =    fr(k)/df
      fkturb(nrw,k) =    0.16*fcspan*fcspan/(df*df)
   15 continue
c
   16 continue
c
   27 continue
c
c     end of the setup used only on the first call to this subroutine.
c
   25 continue
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c    evaluate and smooth the pressure gradients if ypluswall is < -10.0.
c
      if(ypluswall.lt.-10.0) call set_pwallgrad
c
c********************************************************************************
c********************************************************************************
c      first work out the viscous stresses on the streamwise (k = constant) faces of the
c      elements in the do 50 j loop
c
c  q3d   
      if(km.eq.2) go to 555
c  end q3d
c
      do 50 j=2,jm
c
      nrw    = nrow(j)
      j1     = jstart(nrw)
      jrow   = j - j1 + 1
      jtrhub = jtran_k1(nrw)
      jtrtip = jtran_km(nrw)
      jledge = jle(nrw)
      jtedge = jte(nrw)
      wrel   = wrad(j)
c
c     evaluate the average velocities etc on the streamwise faces
c     of the elements.
      do 40 i=1,imm1
c
c
      do 30 k=1,km
      area(k)  = sqrt(asx(i,j,k)*asx(i,j,k)+asr(i,j,k)*asr(i,j,k))
      vxavg(k) = 0.25*(vx(i,j,k)+vx(i,j-1,k)+vx(i+1,j,k)+vx(i+1,j-1,k))
      vravg(k) = 0.25*(vr(i,j,k)+vr(i,j-1,k)+vr(i+1,j,k)+vr(i+1,j-1,k))
      vtavg(k) = 0.25*(vt(i,j,k)+vt(i,j-1,k)+vt(i+1,j,k)+vt(i+1,j-1,k))
      roavg(k) = 0.25*(ro(i,j,k)+ro(i,j-1,k)+ro(i+1,j,k)+ro(i+1,j-1,k))
      tavg(k)  = 0.25*(t_static(i,j,k)   + t_static(i,j-1,k)
     &               + t_static(i+1,j,k) + t_static(i+1,j-1,k))
      ravg(k)  = 0.5*(r(j,k) + r(j-1,k))
      wtavg(k) = vtavg(k) - wrel*ravg(k)
      wabsq    = vxavg(k)*vxavg(k)+vravg(k)*vravg(k)+wtavg(k)*wtavg(k)
      wabs(k)  = sqrt(wabsq)
      wtb(k)   = vtavg(k) - whub(j)*ravg(k)
      if(k.gt.kmid) wtb(k) = vtavg(k) - wtip(j)*ravg(k)
      wbound(k)= sqrt(vxavg(k)*vxavg(k)+vravg(k)*vravg(k)+wtb(k)*wtb(k))
c
   30 continue
c
c     evaluate the viscosity from a power law if reyno is negative
c
      if(reyno.lt.0.0001) then
      vislam = (abs(reyno)/100000.0) * (tavg(kmid)/288.0)**0.62
      tcond  = cp*vislam/prandtl
      ftcond = tcond/vislam
      end if
c
c    save the viscosity it is only used to calculate and write out the reynolds number.
      if(i.eq.imid.and.j.eq.jtedge) viscosy(nrw) = vislam
c
c     calculate the reynolds number based on flow conditions at the trailing edge.
c
      if(i.eq.imid.and.j.eq.jtedge) then
           rovexit(nrw)  = roavg(kmid)*wabs(kmid)
           viscosy(nrw)  = vislam
           reynolds(nrw) = chord(nrw)*rovexit(nrw)/vislam
      end if
c
c************************************************************************
c      evaluate the viscous stresses on the hub and casing.
c
c************************************************************************
c      first the hub
c**********************************************************************************
c**********************************************************************************
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then

           perpk1    =  vol(i,j,1)/area(1)
           yplus_old =  yplus_k1(i,j)
c
         call wallfun(i,j,1,1,perpk1,dpds_cell(i,j,1),roavg(1),
     &                twallk1,yplus_old,wbound(1),yplus_new)
           yplus_k1(i,j) = amin1(1000.0,yplus_new)
c
      go to 345

      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************
c
      if(ypluswall.lt.5.) then
c
      perp   =  vol(i,j,1)/area(1)
      re     =  perp*roavg(2)*wbound(2)/vislam
      relog  =  1.0/alog(re)
c
c    allow for roughness
c
      rough = rough_h(nrw)
c
      if(rough.gt.1.0e-7) then
           if(rough.gt.perp) rough = perp
           rek   = re*rough/perp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if
c
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c
c     allow for roughness if rough > 1.0e-7 .
      if(rough.gt.1.0e-7) then
           twallk1 = 0.5*cf*roavg(2)*wbound(2)*wbound(2)
           vstar   = sqrt(twallk1/roavg(2))
           plusk   = rough*vstar*roavg(2)/vislam
           relog10 = log10(perp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
      end if
c
           twallk1    = 0.5*cf*roavg(2)*wbound(2)*wbound(2)
c
      else
c
c   if ypluswall > 5 use the ypluswall value to calculate the skin friction
c
           twallk1 = cfwall*roavg(1)*wbound(1)*wbound(1)
c
      endif
c
  345 continue
c
c************************************************************************
c****************************************************************************
c    next on the casing
c
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then

           perpkm     =  vol(i,j,kmm1)/area(km)
           yplus_old  =  yplus_km(i,j)
c
        call wallfun(i,j,kmm1,km,perpkm,dpds_cell(i,j,kmm1),roavg(kmm1),
     &              twallkm,yplus_old,wbound(kmm1),yplus_new)
           yplus_km(i,j) = amin1(1000.0,yplus_new)
      go to 355
c
      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************
c
      if(ypluswall.lt.5.) then
c
      perp    = vol(i,j,kmm1)/area(km)
      re      = perp*roavg(kmm1)*wbound(kmm1)/vislam
      relog   = 1.0/alog(re)
c
c    allow for roughness
c
      rough = rough_t(nrw)
c
      if(rough.gt.1.0e-7) then
           if(rough.gt.perp) rough = perp
           rek   = re*rough/perp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if 
c
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c
c     allow for roughness if rough > 1.0e-7 .
      if(rough.gt.1.0e-7) then
           twallkm = 0.5*cf*roavg(kmm1)*wbound(kmm1)*wbound(kmm1)
           vstar   = sqrt(twallkm/roavg(kmm1))
           plusk   = rough*vstar*roavg(kmm1)/vislam
           relog10 = log10(perp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
      end if
c
           twallkm  = 0.5*cf*roavg(km m1)*wbound(kmm1)*wbound(kmm1)
c
      else
c
c   if ypluswall > 5 use the ypluswall value to calculate the skin friction
c.
           twallkm = cfwall*roavg(km)*wbound(km)*wbound(km)
c
      endif
c
  355 continue
c
c******************************************************************************
c  jdd addition. jan/14.  evaluate yplus on the endwalls.
c
      vstar     = sqrt(twallk1/roavg(1))
      do k=1,kmid
      perp          = sqrt(dwallsq(i,j,k))
      yplsk1        = vstar*perp*roavg(1)/vislam
      y_plus(i,j,k) = amin1(yplsk1,1000.0)
      end do
c
c
      vstar     = sqrt(twallkm/roavg(km))
      do k = kmid+1,kmm1
      perp          = sqrt(dwallsq(i,j,k))
      yplskm        = vstar*perp*roavg(km)/vislam
      y_plus(i,j,k) = amin1(yplskm,1000.0)
      end do
c
c   end of jan/14 addition
c******************************************************************************
c     calculate the frictional stresses on the hub
      fmult      = twallk1*area(1)/wbound(1)
      if(ibound.eq.1.or.ibound.gt.2) fmult = 0.0
      xstres(1)  = fmult*vxavg(1)
      rstres(1)  = fmult*vravg(1)
      tstres(1)  = fmult*wtb(1)
c     calculate the frictional stresses on the casing
      fmult      = -twallkm*area(km)/wbound(km)
      if(ibound.ge.2) fmult = 0.0
      xstres(km) = fmult*vxavg(km)
      rstres(km) = fmult*vravg(km)
      tstres(km) = fmult*wtb(km)
c  calculate the shear work and heat flow on the endwalls.
          wvisc(1)   = tstres(1)*whub(j)*ravg(1)
          qflow(1)   = 0.0
          wvisc(km)  = tstres(km)*wtip(j)*ravg(km)
          qflow(km)  = 0.0
c******************************************************************************
c******************************************************************************
c      evaluate the viscous stresses on the streamwise faces of the
c      elements away from the end walls.
c
c   jdd added and changed next section  jan/14.
c    to reduce the turbulent viscosity in the sublayer and buffer region.
      do 35 k=2,kmm1
      if(y_plus(i,j,k).le.yplam) fyplus = 0.0
      if(y_plus(i,j,k).gt.yplam)  then
               xfac = (y_plus(i,j,k) - yplam)/(ypturb- yplam)
               if(xfac.gt.1.0) xfac = 1.0
               fyplus = xfac*xfac*(3.0 - 2.0*xfac)
      end if
c    
      dperp   = vol(i,j,k)/area(k)
      vlam(k) = fklam(nrw,k)*vislam/dperp
      vturb(k)= fkturb(nrw,k)*roavg(k)*abs(wabs(k+1)-wabs(k-1))*fmixup
      vturb(k)= vturb(k)*fyplus
   35 continue
c     end jdd jan/14 addition
c******************************************************************************
c******************************************************************************
c     check for transition
c     first on the hub.
      if(jrow.lt.jtrhub) go to 34
      rmax = 0.0
      do 37 k=2,kmid
      ratvis = vturb(k)/vlam(k)
      if(ratvis.lt.rmax) go to 37
      rmax = ratvis
   37 continue
      if(rmax.gt.ftrans) go to 38
   34 continue
      do 39 k = 2,kmid
   39 vturb(k) = 0.0
   38 continue
c    next on the casing
      if(jrow.lt.jtrtip) go to 46
      rmax = 0.0
      do 47 k = kmid,kmm1
      ratvis = vturb(k)/vlam(k)
      if(ratvis.lt.rmax) go to 47
      rmax = ratvis
   47 continue
      if(rmax.gt.ftrans) go to 49
   46 continue
      do 48 k=kmid,kmm1
   48 vturb(k) = 0.0
   49 continue
c******************************************************************************
c     calculate the stresses on the streamwise faces of the elements.
      do 41 k=2,kmm1
      fmult     = area(k)*(vlam(k)+vturb(k))
      xstres(k) = fmult*(vxavg(k+1)-vxavg(k-1))
      rstres(k) = fmult*(vravg(k+1)-vravg(k-1))
      vistot    = vislam*(1.0 + vturb(k)/vlam(k))
      visc_rat(i,j,k) = vturb(k)/vlam(k)
      tstres(k) = fmult*(vtavg(k+1)-vtavg(k-1))
     &          - vistot*area(k)*vtavg(k)/ravg(k)
      qflow(k)  =  fmult*ftcond*(tavg(k+1)-tavg(k-1))
      wvisc(k)  =  xstres(k)*vxavg(k) + tstres(k)*vtavg(k)
     &          +  rstres(k)*vravg(k)
   41 continue 
c   
c******************************************************************************
c      form the component of viscous force due to the difference of the
c      stresses on adjacent streamwise surfaces.
c
      rf_vis1 = 1.-rf_vis
      do 45 k=1,kmm1
      xforce(i,j,k)=rf_vis1*xforce(i,j,k)+rf_vis*(xstres(k)-xstres(k+1))
      rforce(i,j,k)=rf_vis1*rforce(i,j,k)+rf_vis*(rstres(k)-rstres(k+1))
      tforce(i,j,k)=rf_vis1*tforce(i,j,k) +
     &              rf_vis*(tstres(k)*ravg(k)-tstres(k+1)*ravg(k+1))
      qsource(i,j,k) = rf_vis1*qsource(i,j,k) 
     &     + rf_vis*(qflow(k)  - qflow(k+1) + wvisc(k) - wvisc(k+1))
c
   45 continue
c
c************************************************************************
c     work out the entropy generation rate per unit volume only if ifend = 1 or ifprint = 1.
c     first on the streamwise surfaces.
c
      if(ifprint.eq.1.or.ifend.eq.1) then
      do 57 k = 1,kmm1
      xgen  = (xstres(k)+xstres(k+1))*(vxavg(k+1) -vxavg(k))
      rgen  = (rstres(k)+rstres(k+1))*(vravg(k+1) -vravg(k))
      tgen  = (tstres(k)+tstres(k+1))*(wtavg(k+1) -wtavg(k))
      tavgg = (tavg(k) + tavg(k+1))
      qgen  = 2.0*(qflow(k) + qflow(k+1))*(tavg(k+1) - tavg(k))/tavgg
      sgen(i,j,k) = (xgen + rgen + tgen + qgen)/vol(i,j,k)/tavgg
   57 continue
      end if
c     end of entropy generation rate calculatiion on streamwise surfaces
c**************************************************************************
c     end of i loop
   40 continue
c     end of j loop
   50 continue
c   end of calculating the viscous stresses on the streamwise surfaces
  555 continue
c
c
c********************************************************************************
c********************************************************************************
c********************************************************************************
c      now work out the viscous stresses on the blades and bladewise
c      surfaces in the do 100 j  loop.
c********************************************************************************
c********************************************************************************
c********************************************************************************
c
      do 100 j=2,jm
c
      nrw    = nrow(j)
      j1     = jstart(nrw)
      jrow   = j - j1 + 1
      jtrlow = jtran_i1(nrw)
      jtrup  = jtran_im(nrw)
      jledge = jle(nrw)
      jtedge = jte(nrw)
      wrel   = wrad(j)
c
c      first evaluate the average velocities etc on the bladewise surfaces.
c
      do 90 k=1,kmm1
      ravg(k)  = ravg_cell(j,k)
      do 60 i=1,im
      area(i)  = sqrt(abx(i,j,k)*abx(i,j,k) + abr(i,j,k)*abr(i,j,k)
     &              +abt(j,k)*abt(j,k))
      vxavg(i) = 0.25*(vx(i,j,k)+vx(i,j-1,k)+vx(i,j-1,k+1)+vx(i,j,k+1))
      vravg(i) = 0.25*(vr(i,j,k)+vr(i,j-1,k)+vr(i,j-1,k+1)+vr(i,j,k+1))
      vtavg(i) = 0.25*(vt(i,j,k)+vt(i,j-1,k)+vt(i,j-1,k+1)+vt(i,j,k+1))      
      wtavg(i) = 0.25*(wt(i,j,k)+wt(i,j-1,k)+wt(i,j-1,k+1)+wt(i,j,k+1))
      roavg(i) = 0.25*(ro(i,j,k)+ro(i,j-1,k)+ro(i,j-1,k+1)+ro(i,j,k+1))
      tavg(i)  = 0.25*(t_static(i,j,k)     + t_static(i,j-1,k)
     &               + t_static(i,j-1,k+1) + t_static(i,j,k+1))
      wabsq    = vxavg(i)*vxavg(i)+vravg(i)*vravg(i)+wtavg(i)*wtavg(i)
      wabs(i)  = sqrt(wabsq)
      wbound(i) = wabs(i)
   60 continue
c
c**********************************************************************************
c     evaluate the viscosity from a power law if reyno is negative
c
      if(reyno.lt.0.0001) then
      vislam = (abs(reyno)/100000.0) * (tavg(imid)/288.0)**0.62
      tcond  = cp*vislam/prandtl
      ftcond = tcond/vislam
      end if
c
c**********************************************************************************
c      go to 65 to work out viscous forces before the leading edge and in the wakes.
c    
      if(j.le.jledge.or.j.gt.jtedge) go to 65
c
c**********************************************************************************
c**********************************************************************************
c      now work out the stresses on the blade surfaces.
c
c      skip the tip gap
      if((ktips(nrw).gt.0).and.(k.ge.ktips(nrw)).and.(k.lt.ktipe(nrw))) 
     &    go to 51
c
c**********************************************************************************
c**********************************************************************************
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then

           perpi1    =  vol(1,j,k)/area(1)
           yplus_old =  yplus_i1(j,k)
c
         call wallfun(1,j,k,1,perpi1,dpds_cell(1,j,k),roavg(1),
     &                twalli1,yplus_old,wbound(2),yplus_new)
          yplus_i1(j,k) = amin1(1000.0,yplus_new)
c
      go to 365
c
      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************
      if(ypluswall.lt.5.) then
c
      perp   =  vol(1,j,k)/area(1)
      re     =  roavg(2)*wabs(2)*perp/vislam
      relog  =  1.0/alog(re)
c
c    allow for roughness
c
      rough  = rough_l(nrw)
      if(rough.gt.1.0e-7) then
           if(rough.gt.perp) rough = perp
           rek   = re*rough/perp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if
c
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c
c     allow for roughness if rough > 1.0e-7
      if(rough.gt.1.0e-7) then
           twalli1    = 0.5*cf*roavg(2)*wbound(2)*wbound(2)
           vstar   = sqrt(twalli1/roavg(2))
           plusk   = rough*vstar*roavg(2)/vislam
           relog10 = log10(perp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
      end if 
c
      twalli1    = 0.5*cf*roavg(2)*wbound(2)*wbound(2)
c
      else
c   if ypluswall > 5 use the ypluswall value to calculate the skin friction
      twalli1    = cfwall*roavg(1)*wabs(1)*wabs(1)
c   end of ypluswall < 5 loop
      endif
c
  365 continue
c
c******************************************************************************
c******************************************************************************
c  next the i = im , upper, surface
c
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then

           perpim    =  vol(imm1,j,k)/area(im)
           yplus_old =  yplus_im(j,k)
c
           call wallfun(imm1,j,k,im,perpim,dpds_cell(imm1,j,k),
     &             roavg(imm1),twallim,yplus_old,wbound(imm1),yplus_new)
           yplus_im(j,k) = amin1(1000.0,yplus_new)
c
      go to 375
c
      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************
c
      if(ypluswall.lt.5.) then
c
      perp   =  vol(imm1,j,k)/area(im)
      re     =  roavg(imm1)*wabs(imm1)*perp/vislam
      relog  =  1.0/alog(re)
c
c    allow for roughness
c
       rough = rough_u(nrw)
       if(rough.gt.1.0e-7) then
           if(rough.gt.perp) rough = perp
           rek   = re*rough/perp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if
c
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c    allow for roughness if  rough > 1.0e-7 .
      if(rough.gt.1.0e-7) then
           twallim = 0.5*cf*roavg(imm1)*wbound(imm1)*wbound(imm1)
           vstar   = sqrt(twallim/roavg(imm1))
           plusk   = rough*vstar*roavg(imm1)/vislam
           relog10 = log10(perp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
      end if
c
           twallim   = 0.5*cf*roavg(imm1)*wbound(imm1)*wbound(imm1)
c
      else
c   if ypluswall > 5 use the ypluswall value to calculate the skin friction
           twallim   = cfwall*roavg(im)*wabs(im)*wabs(im)
c    end of ypluswall < 5 loop
      endif
c
  375 continue
c
c******************************************************************************
c     jdd addition jan/14  . evaluate yplus on the blade surfaces.
c
      vstar     = sqrt(twalli1/roavg(1))
      do i=1,imid
      yplsi1        = vstar*sqrt(dwallsq(i,j,k))*roavg(1)/vislam
      yplsp         = y_plus(i,j,k)
      y_plus(i,j,k) = amin1(yplsi1,yplsp)
      end do
c
c
      vstar     = sqrt(twallim/roavg(im))
      do i = imid+1,imm1
      yplsim        = vstar*sqrt(dwallsq(i,j,k))*roavg(im)/vislam
      yplsp         = y_plus(i,j,k)
      y_plus(i,j,k) = amin1(yplsim,yplsp)
      end do
c   end of jan/14 addition
c
      fmult     = twalli1*area(1)/wabs(1)
      xstres(1) = fmult*vxavg(1)
      rstres(1) = fmult*vravg(1)
      tstres(1) = fmult*wtavg(1)
c
      fmult     = -twallim*area(im)/wabs(im)
      xstres(im)= fmult*vxavg(im)
      rstres(im)= fmult*vravg(im)
      tstres(im)= fmult*wtavg(im)
      qflow(1)  = 0.0
      wvisc(1)  = tstres(1)*wrel*ravg(k)
      qflow(im) = 0.0
      wvisc(im) = tstres(im)*wrel*ravg(k)
c
   51 continue
c
c*********************************************************************************
c*********************************************************************************
c     now work out the stresses on the bladewise surfaces of the main flow.
c     within the blade passage
c
c   jdd added and changed the next section  jan/14.
c    reduce the turbulent viscosity in the sublayer and buffer region.
      do 55 i=2,imm1     
      if(y_plus(i,j,k).le.yplam) fyplus = 0.0
      if(y_plus(i,j,k).gt.yplam)  then
               xfac = (y_plus(i,j,k) - yplam)/(ypturb- yplam)
               if(xfac.gt.1.0) xfac = 1.0
               fyplus = xfac*xfac*(3.0 - 2.0*xfac)
      end if

      dperp   = vol(i,j,k)/area(i)
      vlam(i) = filam(nrw,i)*vislam/dperp
      vturb(i)= fiturb(nrw,i)*roavg(i)*abs(wabs(i+1)-wabs(i-1))*fmixup
      vturb(i)= vturb(i)*fyplus
   55 continue
c   end of jdd addition jan/14
c*********************************************************************************
c*********************************************************************************
c     check for transition
c
      if(jrow.lt.jtrlow) go to 70
      rmax=0.0
      imax=2
      do 56 i=2,imid
      ratvis = vturb(i)/vlam(i)
      if(ratvis.lt.rmax) go to 56
      rmax = ratvis
      imax = i
   56 continue
      if(rmax.gt.ftrans)  go to 71
   70 continue
      do 72 i = 2,imid
   72 vturb(i)= 0.0
   71 continue
c
      if(jrow.lt.jtrup) go to 69
      rmax = 0.0
      imax = imm1
      do 73 i = imid,imm1
      ratvis  = vturb(i)/vlam(i)
      if(ratvis.lt.rmax) go to 73
      rmax = ratvis
      imax = i
   73 continue
      if(rmax.gt.ftrans) go to  77
   69 continue
      do 74 i  = imid,imm1
   74 vturb(i) = 0.0
   77 continue
c
c*********************************************************************************
c     form the viscous stresses away from the blade surfaces
      do 76 i=2,imm1
      fmult        = (vlam(i)+vturb(i))*area(i)
      xstres(i)    = fmult*(vxavg(i+1)-vxavg(i-1))
      rstres(i)    = fmult*(vravg(i+1)-vravg(i-1))
      tstres(i)    = fmult*(wtavg(i+1)-wtavg(i-1))
      visc_rat(i,j,k) = visc_rat(i,j,k) + vturb(i)/vlam(i)
      qflow(i)     =  fmult*ftcond*(tavg(i+1)-tavg(i-1))
      wvisc(i)     =  xstres(i)*vxavg(i) + rstres(i)*vravg(i) 
     &             +  tstres(i)*vtavg(i)
   76 continue
c*********************************************************************************
c     set the stresses on the periodic boundaries above the tip
c  q3d
      if(km.eq.2) go to 78
c  end q3d
c
      if((ktips(nrw).gt.0).and.(k.ge.ktips(nrw)).and.(k.lt.ktipe(nrw))) 
     &    then
c
      xstres(1)  = 0.5*(xstres(2)+xstres(imm1))
      rstres(1)  = 0.5*(rstres(2)+rstres(imm1))
      tstres(1)  = 0.5*(tstres(2)+tstres(imm1))
      xstres(im) = xstres(1)
      tstres(im) = tstres(1)
      rstres(im) = rstres(1)
      qflow(1)   = 0.5*(qflow(2) + qflow(imm1))
      wvisc(1)   = 0.5*(wvisc(2) + wvisc(imm1))
      qflow(im)  = qflow(1)
      wvisc(im)  = wvisc(1)
c
      endif
c
   78 continue
c
c     end of the calculation of viscous forces within the blade row
c*********************************************************************************
c*********************************************************************************
c
      go to 75
c
c      form the viscous stresses on the bladewise surfaces
c      outside a blade row where the mixing length is taken to be
c      fracpw or fracpup * the blade pitch
c
   65 continue
c
c*********************************************************************************
c*********************************************************************************
c     form the viscous stresses in the wake. blend them to the blade values using facj
c
      if(j.gt.jtedge) then
c
      facj = 0.1*(j-jtedge)
      if(facj.gt.1.)  facj = 1.
c
      do 66 i=1,im
      isub = i
      if(i.eq.im) isub=imm1
      im1 = i-1
      ip1 = i+1
      if(i.eq.1)  im1 = imm1
      if(i.eq.im) ip1 = 2
      fblend       = (facj*fiwake(nrw,i)+(1.-facj)*fiturb(nrw,i))*fmixup
      dperp        = vol(isub,j,k)/area(i)
      vlam(i)      = filam(nrw,i)*vislam/dperp
      vturb(i)     = fblend*roavg(i)*abs(wabs(ip1)-wabs(im1))
      visc_rat(i,j,k) = visc_rat(i,j,k) + vturb(i)/vlam(i)
      fmult        = area(i)*(vlam(i)+vturb(i))
      xstres(i)    = fmult*(vxavg(ip1)-vxavg(im1))
      rstres(i)    = fmult*(vravg(ip1)-vravg(im1))
      tstres(i)    = fmult*(wtavg(ip1)-wtavg(im1))
      qflow(i)     = fmult*ftcond*(tavg(ip1) - tavg(im1))
      wvisc(i)     = xstres(i)*vxavg(i) + tstres(i)*vtavg(i)
     &             + rstres(i)*vravg(i)
   66 continue
c
      endif
c*********************************************************************************
c     form the viscous stresses upstream of the leading edge.
c
      if(j.le.jledge) then
c
      do 67 i=1,im
      isub = i
      if(i.eq.im) isub=imm1
      im1 = i-1
      ip1 = i+1
      if(i.eq.1)  im1 = imm1
      if(i.eq.im) ip1 = 2
      dperp       = vol(isub,j,k)/area(i)
      vlam(i)     = filam(nrw,i)*vislam/dperp
      vturb(i)    = fmixup*fiup(nrw,i)*roavg(i)*abs(wabs(ip1)-wabs(im1))
      visc_rat(i,j,k)= visc_rat(i,j,k) + vturb(i)/vlam(i)
      fmult       = area(i)*(vlam(i)+vturb(i))
      xstres(i)   = fmult*(vxavg(ip1)-vxavg(im1))
      rstres(i)   = fmult*(vravg(ip1)-vravg(im1))
      tstres(i)   = fmult*(wtavg(ip1)-wtavg(im1))
      qflow(i)    = fmult*ftcond*(tavg(ip1) - tavg(im1))
      wvisc(i)    = xstres(i)*vxavg(i) + tstres(i)*vtavg(i)
     &            + rstres(i)*vravg(i)
   67 continue
c
      endif
c
c   make the viscous stresses periodic upstream and downstream of the blade.
c

      xstres(1)  = 0.5*(xstres(im)+xstres(1))
      rstres(1)  = 0.5*(rstres(im)+rstres(1))
      tstres(1)  = 0.5*(tstres(im)+tstres(1))
      xstres(im) = xstres(1)
      rstres(im) = rstres(1)
      tstres(im) = tstres(1)
      qflow(1)   = 0.5*(qflow(im) + qflow(1))
      wvisc(1)   = 0.5*(wvisc(im) + wvisc(1))
      qflow(im)  = qflow(1)
      wvisc(im)  = wvisc(1)
c
c*********************************************************************************
c*********************************************************************************
c   end of calculation of viscous forces on the bladewise (i) surfaces
c
   75 continue
c
c      complete the viscous force terms by adding the difference of the
c      stresses on the bladewise faces of the elements. relaxing the changes by rfvis .
c
c   q3d
      if(km.ne.2) then
c   end q3d
      do 80 i=1,imm1
      xforce(i,j,k)  = xforce(i,j,k)+rf_vis*(xstres(i)-xstres(i+1))
      rforce(i,j,k)  = rforce(i,j,k)+rf_vis*(rstres(i)-rstres(i+1))
      tforce(i,j,k)=tforce(i,j,k)+rf_vis*(tstres(i)-tstres(i+1))*ravg(k)
      qsource(i,j,k) = qsource(i,j,k) + rf_vis*(qflow(i) - qflow(i+1)
     &               + wvisc(i) - wvisc(i+1))    
   80 continue
c
      else
c   q3d  .  if km = 2
      do 81 i=1,imm1
      xforce(i,j,k)  = rf_vis1*xforce(i,j,k)
     &               + rf_vis*(xstres(i)-xstres(i+1))
      rforce(i,j,k)  = rf_vis1*rforce(i,j,k)
     &               + rf_vis*(rstres(i)-rstres(i+1))
      tforce(i,j,k)  = rf_vis1*tforce(i,j,k)
     &               + rf_vis*(tstres(i)-tstres(i+1))*ravg(k)
      qsource(i,j,k) = rf_vis1*qsource(i,j,k)
     &   + rf_vis*(qflow(i) - qflow(i+1) + wvisc(i) - wvisc(i+1)) 
   81 continue
c  end q3d
c
      end if
c
c
c***************************************************************************
c     complete the entropy generation rate per unit volume if ifend = 1 or ifprint = 1.
c     by adding the generation due to the bladewise surfaces.
c     non dimensionalise the entropy generation rate by:
c     0.001*(one blade row pressure change) *(blade exit velocity)/(blade pitch).
c
      if(ifprint.eq.1.or.ifend.eq.1) then
      delp_blade = abs((po1(kmid) - 0.5*(pdown_hub + pdown_tip)))/nrows
      tau_ref    = 0.001*delp_blade
      rho_ref    =  po1(kmid)/rgas/to1(kmid)
      vel_ref    = sqrt(2*delp_blade/rho_ref)
      pitch      = 2*3.14159*r(1,kmid)/nblade(1)
      sref       =  tau_ref*vel_ref/pitch
      do 91 i = 1,imm1
      xgen  = (xstres(i)+xstres(i+1))*(vxavg(i+1) -vxavg(i))
      rgen  = (rstres(i)+rstres(i+1))*(vravg(i+1) -vravg(i))
      tgen  = (tstres(i)+tstres(i+1))*(wtavg(i+1) -wtavg(i))
      tavgg = (tavg(i) + tavg(i+1))
      qgen  = 2.0*(qflow(i) + qflow(i+1))*(tavg(i+1) - tavg(i))/tavgg
      sgen_now     =  (xgen + rgen + tgen + qgen)/vol(i,j,k)/tavgg
      temp4(i,j,k) =  (sgen(i,j,k) + sgen_now)/sref
   91 continue
c
      end if
c
c     end of entropy generation rate calculatiion on bladewise surfaces
c******************************************************************************8
c
c     end of k loop for viscous calculations on the streamwise surfaces.
   90 continue
c
c   end of j loop
  100 continue
c*********************************************************************************
c*********************************************************************************
c    distribute entropy generation rate to the cell corners.
c
      if(ifprint.eq.1.or.ifend.eq.1) then
c
           do k=1,km-1
           do i=1,im-1
           temp4(i,1,k) = 0.0
           end do
           end do
c
           call cell_to_node(temp4,sgen)
c
c*********************************************************************************
cc  jdd added  24/9/2018  so  visc_rat can be plotted  
c   end 24/9/18 addition
      do 28 k=1,km,kmm1
      do 28 j=1,jm
      do 28 i=1,im
      visc_rat(i,j,k) =  0.0
   28 continue
      do 29 k=1,km
      do 29 j=1,jm
      do 29 i=1,im,imm1
      visc_rat(i,j,k) =  0.0
   29 continue
c     end 24/9/2018 addition
c**********************************************************************************
c**********************************************************************************
      end if
c
c*********************************************************************************
c    end of subroutine loss 
c*********************************************************************************
c
      return
      end
c
c******************************************************************************
c******************************************************************************
c*****************************************************************************
c******************************************************************************
c
      subroutine inpint
c
      include  'commall-open-19.2'
c
      dimension  sumfn(maxki),sumf_in(maxki),ans(maxki)
c
c     this subroutine interpolates in the inlet flow data at kin points
c     spaced by fr_in(k) to produce new data at km points spaced via fr(k).
c
      write(6,*) '  entered  inpint '
c
c    make  fr(k)  and fr_in(k)  both sum to 1.0 .
c
      sumfn(1) = 0.0
      do 20 k=2,km
           sumfn(k) = sumfn(k-1) + fr(k-1)
   20 continue
      do 30 k=2,km
           sumfn(k) = sumfn(k)/sumfn(km)
   30 continue
c
      sumf_in(1)=0.0
      do 40 k=2,kin
           sumf_in(k) = sumf_in(k-1) + fr_in(k-1)
   40 continue
      do 50 k=2,kin
           sumf_in(k) = sumf_in(k)/sumf_in(kin)
   50 continue
c
c
      write(6,*) ' sum fr for actual grid points '
      write(6,10)(sumfn(k),k=1,km)
      write(6,*)
      write(6,*) ' sum fr  for boundary condition points'
      write(6,10)(sumf_in(k),k=1,kin)
      write(6,*)
   10 format(8f10.5)
c
c    interpolate for the stagnation pressures at the grid points.
      do 100 k=1,km
           call intp(kin,sumf_in,po1,sumfn(k),ans(k))
  100 continue
      do 110 k=1,km
           po1(k)=ans(k)
  110 continue
c
      write(6,*) 'po1 at grid points on inlet boundary.'
      write(6,11)(po1(k),k=1,km)
   11 format(8f10.1)
c
c    interpolate for the stagnation temperatures at the grid points.
      do 120 k=1,km
           call intp(kin,sumf_in,to1,sumfn(k),ans(k))
  120 continue
      do 130 k=1,km
           to1(k)=ans(k)
  130 continue
c
      write(6,*)  'to1 at grid points on inlet boundary.'
      write(6,10)(to1(k),k=1,km)
c
c    interpolate for the tangential velocity at the grid points.
      do 140 k=1,km
           call intp(kin,sumf_in,vtin,sumfn(k),ans(k))
  140 continue
      do 150 k=1,km
           vtin(k)=ans(k)
  150 continue
c
      write(6,*)  'vtin at grid points on inlet boundary.'
      write(6,10)(vtin(k),k=1,km)
c
c    interpolate for the meridional velocity at the grid points.
      do 160 k=1,km
           call intp(kin,sumf_in,vm1,sumfn(k),ans(k))
  160 continue
      do 170 k=1,km
           vm1(k)=ans(k)
  170 continue
c
      write(6,*)  'vm at grid points on inlet boundary.'
      write(6,10)(vm1(k),k=1,km)
c
      if(ipout.eq.3) then
c
c    interpolate for the static pressures at the grid points.
      do 180 k=1,km
           call intp(kin,sumf_in,pd,sumfn(k),ans(k))
  180 continue
      do 190 k=1,km
           pd(k)=ans(k)
  190 continue
c
      write(6,*)  'p static at grid points on exit boundary.'
      write(6,11)(pd(k),k=1,km)
c
      end if
c
c    interpolate for the meridional pitch angle at the grid points.
      do 200 k=1,km
           call intp(kin,sumf_in,br,sumfn(k),ans(k))
  200 continue
      do 210 k=1,km
           br(k)=ans(k)
  210 continue
c
      write(6,*) 'pitch angle at grid points on inlet boundary.'
      write(6,10)(br(k),k=1,km)
c
c    interpolate for the yaw angle at the grid points.
      do 220 k=1,km
           call intp(kin,sumf_in,bs,sumfn(k),ans(k))
  220 continue
      do 230 k=1,km
           bs(k)=ans(k)
  230 continue
c
      write(6,*)  'yaw angle at grid points on inlet boundary.'
      write(6,10)(bs(k),k=1,km)
c
c
      return
      end
c
c******************************************************************************c
c******************************************************************************c
c******************************************************************************c
c
      subroutine grid_down(k,j1,j2,s1,s2,ngap,sdist,xint,rint,nextrap,
     &  beta_dwn1,beta_dwn2,ifcusp,icusp,lcusp,lcuspup,ifangles)
c
c     this subroutine sets the grid downstream of the trailing edge
c     down to the mixing plane or downstream boundary.
c     it also fits a cusp at the trailing edge and sets the grid angles
c     downstream of the trailing edge
c
      include  'commall-open-19.2'
c
      dimension  sdist(jd),xint(jd),rint(jd),snew(jd),ansx(jd),
     &           ansr(jd),th_upp(jd),th_thick(jd),th_mid(jd)
c
c******************************************************************************c
c    calculate the meridional grid spacings downstream of the blade. using
c    a geometric progression with ratio  "rat"  which is set by "solve".
c******************************************************************************c
c     set the new meridional grid using the geometric progression found by  "solve" .
c
      snew(j1) = s1
      snew(j2) = s2
      sdiff    = s2-s1
      nint = j2-j1
      r1   = smerid(j1,k)-smerid(j1-1,k)
      call solve(nint,sdiff,r1,rat,k)
      xd=r1
      do 10 j = j1+1,j2-1
      snew(j) = snew(j-1) + xd
      xd      = xd*rat
   10 continue
c
c******************************************************************************c
c     interpolate to find xsurf  and rsurf at the new grid points
c
      do 30 j = j1,j2
      call intp(ngap,sdist,xint,snew(j),ansx(j))
      call intp(ngap,sdist,rint,snew(j),ansr(j))
   30 continue
      do 31 j = j1,j2
      xsurf(j,k)   = ansx(j)
      rsurf(j,k)   = ansr(j)
   31 continue
c
c******************************************************************************c
c******************************************************************************c
c     maintain the cusp (if any) set in the input data if ifcusp = 0.
c
      if(ifcusp.eq.0) then
c

      if(ifangles.eq.0) then
c     reset the downstream grid angles if the input value of if_angles  is zero.
c
      jstrt      = j1 - nextrap
      do 20    j = jstrt,j2
      xd          = xsurf(j,k) - xsurf(j-1,k)
      rd          = rsurf(j,k) - rsurf(j-1,k)
      snew(j)     = snew(j-1) + sqrt(xd*xd + rd*rd)
      th_upp(j)   = rt_upp(j,k)/rsurf(j,k) 
      th_thick(j) = rt_thick(j,k)/rsurf(j,k)
      th_mid(j)   = th_upp(j) - 0.5*th_thick(j)
   20 continue
c
      dsmer      = snew(j1)  - snew(jstrt)
      slopes     = (th_upp(j1)   - th_upp(jstrt))/dsmer
      slopetk    = (th_thick(j1) - th_thick(jstrt))/dsmer
      slopep     = slopes - slopetk
      slopecent  = 0.5*(slopes+slopep)
c   set betadwn1  from the slope of the blade center line at the trailing edge.
c   if if_angles  was zero.
      betadown1  = slopecent
      betadown2  = betadown1
c
      else
c
c   use the input values of betadwn1 and betadwn2   if if_angles  was not zero.
      betadown1 = tan(beta_dwn1*degrad)/rsurf(j1,k)
      betadown2 = tan(beta_dwn2*degrad)/rsurf(j2,k)
c
      end if
c
c******************************************************************************c
c     find the end of any existing cusp, where the thickmess = 0 .
      js = j1
      do 35 j=j1,j2
           if(rt_thick(j,k).lt.1.0e-4*sdiff) then
                js = j
                go to 36 
           end if
   35 continue
c
   36 continue
c
c     js  is the end point of any existing cusp.
c     do not change the existing cusp upstream of js .
c     change from r*theta to theta to get better extrapolation when the radius changes.
c     set the grid angle, r_theta, and thickness downstream of the existing cusp.
c
      do 40 j = js+1,j2
      fdown     = (snew(j) - snew(j1))/(snew(j2) - snew(j1))
      betadown  = betadown1 + fdown*(betadown2 - betadown1)
      dsdist    = snew(j) - snew(j-1)
      rt_upp(j,k)   = rt_upp(j-1,k) + rsurf(j,k)*betadown*dsdist
      rt_thick(j,k) = 0.0
   40 continue
c
      return
c
c   end of ifcusp = 0 option. 
      end if
c
c*******************************************************************************
c*******************************************************************************
c     now generate a new cusp  if ifcusp is not zero
c
c     change from r*theta to theta to get better extrapolation when the radius changes.
c     reset snew  for those points on the blade and downstream which might be used later.
      lup   = 20
      jstrt = j1 - lup
      snew(jstrt-1) = 0.0
      do 100    j = jstrt,j2
      xd          = xsurf(j,k) - xsurf(j-1,k)
      rd          = rsurf(j,k) - rsurf(j-1,k)
      snew(j)     = snew(j-1) + sqrt(xd*xd + rd*rd)
      th_upp(j)   = rt_upp(j,k)/rsurf(j,k) 
      th_thick(j) = rt_thick(j,k)/rsurf(j,k)
      th_mid(j)   = th_upp(j) - 0.5*th_thick(j)
  100 continue
c
      js = jstrt + 1
      if(js.gt.j1)  js = j1
c
c  search for the start of the trailing edge thinning.
      do 102 j = js,j1
           rthik   = th_thick(j)/th_thick(j-1)
           if(rthik.lt.0.9) go to 103
  102 continue
c
  103 continue
c
c   jcs  is the start of the trailing edge thinning
      jcs = j-1 
c
c     extapolate the blade surfaces from upstream, j = jcs, to the trailing edge 
c     so that the te thinning is removed.
      tkgrad = (th_thick(jcs) - th_thick(jcs-2))/(snew(jcs)-snew(jcs-2))
      ygrad  = (th_mid(jcs)   - th_mid(jcs-2))/(snew(jcs) - snew(jcs-2))     
      do 104 j = jcs,j1
           th_mid(j)   = th_mid(jcs)   +  ygrad*(snew(j) - snew(jcs))
           th_thick(j) = th_thick(jcs) + tkgrad*(snew(j) - snew(jcs))
           if(th_thick(j).lt.0.0) th_thick(j) = 0.0
           th_upp(j)     = th_mid(j)  + 0.5*th_thick(j)
           rt_upp(j,k)   = rsurf(j,k)*th_upp(j)
           rt_thick(j,k) = rsurf(j,k)*th_thick(j)
  104 continue
c
c   end of adjusting the blade upstream of the trailing edge
c
c**********************************************************************************
c     set the slope of the blade surface at the new trailing edge.
c     which is at  jstrt = j1 - lcuspup
c
      jstrt      = j1 - lcuspup
      dsmer      = snew(jstrt)  - snew(jstrt-2)
      slopes     = (th_upp(jstrt)      - th_upp(jstrt-2))/dsmer
      slopetk    = (th_thick(jstrt)    - th_thick(jstrt-2))/dsmer
      slopep     = slopes - slopetk
      slopecent  = 0.5*(slopes+slopep)
c 
      if(ifangles.eq.0) then
c     use the extrapolated blade centre line angle if  if_angles = zero.
	   betadown1 = slopecent
           betadown2 = slopecent
      else
c     use the input downstream grid angles if if_angles > zero .
	   betadown1 = tan(beta_dwn1*degrad)/rsurf(jstrt,k)
	   betadown2 = tan(beta_dwn2*degrad)/rsurf(jstrt,k)
      endif
c
c**********************************************************************************
c     calculate the trailing edge thickness and centreline rtheta value at the.
c     new trailing edge, which is at  j = j1 - lcuspup
c
      ythickte = rt_thick(jstrt,k)
      if(ythickte.lt.0.0)  ythickte = 0.0
      ymidte   = rt_upp(jstrt,k) - 0.5*ythickte
c
c     form the cusp of length  lcusp, starting lcuspup points before the.
c     trailing edge.
c     the cusp is centred on the blade centre line if icusp = 0.
c     the cusp makes the i=1 surface continuous on the cusp if  icusp =  1.
c     the cusp makes the i=im surface continuous on the cusp if icusp = -1.
c     the cusp extends lcuspup upstream onto the solid part of the blade.
c
      do 500 j = jstrt+1,j2
c
      fdown     = (snew(j) - snew(jstrt))/(snew(j2) - snew(jstrt))
      betadown  = betadown1 + fdown*(betadown2 - betadown1)
      dsdist    = snew(j) - snew(j-1)
      ymid      = ymidte  + rsurf(j,k)*slopecent*(snew(j) - snew(jstrt))
      clength   = snew(jstrt+lcusp) - snew(jstrt)
      plength   = snew(j) - snew(jstrt)
      if(plength.gt.clength) plength = clength
c
c     set the downstream grid slope according to  "icusp" and "lcusp".
c     note that the local radius, rsurf, is used to change back from theta to rtheta.
c
      jcuspend  = jstrt+lcusp
      if((jcuspend-j).ge.0)   then
c     if on the cusp
	   ytnew         = ythickte*(1. -plength/clength)
	   rt_thick(j,k) = ytnew
	   if(icusp.eq.0) rt_upp(j,k) = ymid + 0.5*ytnew
	   if(icusp.eq.1) rt_upp(j,k) = rt_upp(j-1,k)
     &                                + rsurf(j,k)*slopes*dsdist
	   if(icusp.eq.-1) rt_upp(j,k) = rt_upp(j-1,k) - rt_thick(j-1,k) 
     &                   + rsurf(j,k)*slopep*dsdist  + ytnew
      else
c     if past end of cusp
	   rt_upp(j,k)   = rt_upp(j-1,k) + rsurf(j,k)*betadown*dsdist
	   rt_thick(j,k) = 0.0
      endif
c
  500 continue
c
c     end of new cusp generation. november 2015 .  
c
      return
      end
c
c**********************************************************************
c
      subroutine grid_up(k,j1,j2,s1,s2,ngap,sdist,xint,rint,nextrap,
     &                 beta_up,ifangles)
c
      include  'commall-open-19.2'
c
c     this subroutine sets the grid upstream of the leading edge,
c     up to the mixing plane or upstream boundary.
c
      dimension sdist(jd),xint(jd),rint(jd),snew(jd),ansx(jd),
     &          ansr(jd),th_upp(jd),th_thick(jd)
c
c      form the new meridional grid spacings "snew" using a geometric progression
c      with expansion ratio "rat" which is calculated by subroutine  "solve" .
c
      snew(j1) = s1
      snew(j2) = s2
      sdiff    = s2 - s1
      nint = j2 - j1
      r1 = smerid(j2+1,k)-smerid(j2,k)
      call solve(nint,sdiff,r1,rat,k)
      xd = r1
      do 10 j = j1+1,j2-1
      jsub    = j2 + j1 - j
      snew(jsub) = snew(jsub+1) - xd
      xd = xd*rat
   10 continue
c
c     set the grid points at the upstream mixing plane to be very closely,
c     but not exactly, coincident, except for the first row.
c
      if(j1.ne.1) snew(j1) = snew(j1) + 0.01*(snew(j1+1)-snew(j1))
c
c******************************************************************************c
c  interpolate to find  xsurf   and  rsurf   at the new grid points.
c
      do 30 j = j1,j2
      call intp(ngap,sdist,xint,snew(j),ansx(j))
      call intp(ngap,sdist,rint,snew(j),ansr(j))
   30 continue
      do 35 j = j1,j2
      xsurf(j,k)   = ansx(j)
      rsurf(j,k)   = ansr(j)
   35 continue
c
c******************************************************************************c
c******************************************************************************c
c     change from r*theta to theta to get a better extrapolation for radial blades.
c
      do 100   j = j1,j2 + nextrap
      th_upp(j)   = rt_upp(j,k)/rsurf(j,k)
      th_thick(j) = rt_thick(j,k)/rsurf(j,k)
  100 continue
c
c     set the slope of the upstream grid as  dtheta/ds, or use the input value 
c     of  beta_up  if  if_angles  is greater than zero
c
      if(ifangles.eq.0) then
      slopet = (th_upp(j2) - 0.5*th_thick(j2) - th_upp(j2+nextrap)
     & + 0.5*th_thick(j2+nextrap))/(smerid(j2,k) - smerid(j2+nextrap,k))
      else
      slopet = tan(beta_up*degrad)/rsurf(j2,k)
      endif
c
c     use the calculated slope to set  rtheta  upstream of the leading edge.
c     multiply by the local radius, rsurf, to change back from   theta  to  rtheta.
c
      rtheta_le = rt_upp(j2,k) - 0.5*rt_thick(j2,k)
      do 20 j = j1,j2-1
      rt_upp(j,k) = rtheta_le + rsurf(j,k)*slopet*(snew(j)-snew(j2))
   20 continue
c
c******************************************************************************c
c******************************************************************************c
c
      return
      end
c
c******************************************************************************
c******************************************************************************c
c
      subroutine solve(npup,xup,dx1,rat,k)
c
c     find the ratio of spacings for the upstream and downstream grids 
c 
c  
      rup = xup/dx1 
      rat = 1.2
      if(rup.lt.npup) rat = 0.9
      m = 1 
  50  rhs = (rat-1.)*rup 
      if(rup.ge.npup) rnew = (rhs+1.)**(1./npup)
      if(rup.lt.npup) rnew = (rat**npup -1.)/rup  + 1.
      diff =abs(rnew/rat-1.0)
c      may be more stable if use  "rat = 0.5*(rnew+rat) "
      rat = rnew
      m = m+1
      if(diff.lt.0.00001) go to 51
      if(m.gt.250) go to 51
      go to 50 
   51 continue
c
      if(m.ge.250) write(6,*) ' warning !!!    gridup/down iteration not
     & converged.'
c
      write(6,1) k, m, rat
    1 format(' k= ',i5,'grid_up/down its= ',i5,'grid expansion ratio= ',
     &       f10.5)
c
c
      return
      end
c
c******************************************************************************c
c******************************************************************************
c
      subroutine newgrid(jmrow,j1,j2,jlerow,jterow,jroths,jrothe,
     & jrotts,jrotte)
c
      include  'commall-open-19.2'
c
c
      dimension upf(jd),onf(jd),downf(jd),sold(jd),snew(jd),
     & xbuf(jd),ybuf(jd),ysnew(jd), xfracup(jd),relspup(jd),
     & xfracon(jd),relspon(jd),xfracdwn(jd),relspdwn(jd)
c
c      this suboutine interpolates in the input data to set up a new grid spacing in the streamwise (j)
c      direction.
c
c      nup  -   is number of upstream points including the l.e. point
c      		i.e. it is the number of upstream elements + 1.
c      non  -   is number of points on blade including the t.e. point but not
c      		the l.e. point. it = the number of elements on the blade.
c      ndown  - is the number of points downstream not including the t.e. point
c     		it = the number of elements downstream of the blade.
c
c
c      the relative spacings are all ordered in the direction of
c      increasing   j  . i.e. upf(1) will tend to be large and upf(nup)
c      will tend to be small.
c
c       upf(nup) is not used as there are nup-1 intervals between nup points
c
      read(5,*)     dummy_input
      write(6,*)    dummy_input
      write(6,*) 'starting to generate a new grid in subroutine newgrid'
c
      read(5,*)     nup,non,ndown
      write(6,1600) nup,non,ndown
 1600 format('  nup, non, ndown = ',3i10)
c
c******************************************************************************
c     input only a few relative grid spacings
c     and generate the grid spacings automatically if  nup, non or ndown = 0.
c
c     read in  the relative grid spacings, upf(j) , onf(j)  and  downfj) ,
c     if  nup, non or ndown > 0 .
c
c******************************************************************************
c******************************************************************************
c
      if(nup.eq.0) then
           read(5,*) nup
           ninup = 0
 880       continue
           ninup = ninup+1
           read(5,*) xfracup(ninup),relspup(ninup)
           if(xfracup(ninup).gt.0.9999) go to 881
           go to 880
 881       continue
           do 886 n=1,nup
                xx = float(n-1)/(nup-1)
                call intp(ninup,xfracup,relspup,xx,upf(n))
 886       continue
c
      else
           read(5,*)     (upf(j),j=1,nup)
      endif
c
c******************************************************************************
      if(non.eq.0) then
           read(5,*) non
           ninon = 0
 882       continue
           ninon = ninon+1
           read(5,*) xfracon(ninon),relspon(ninon)
           if(xfracon(ninon).gt.0.9999) go to 883
           go to 882
 883       continue
           do 887 n=1,non
                xx = float(n-1)/(non-1)
           call intp(ninon,xfracon,relspon,xx,onf(n))
 887  conti   nue
c
      else
           read(5,*)     (onf(j),j=1,non)
      endif
c
c******************************************************************************
      if(ndown.eq.0) then
           read(5,*) ndown
           nindwn = 0
 884       continue
           nindwn = nindwn+1
           read(5,*) xfracdwn(nindwn),relspdwn(nindwn)
           if(xfracdwn(nindwn).gt.0.9999) go to 885
           go to 884
 885       continue
           do 888 n=1,ndown
                xx = float(n-1)/(ndown-1)
                call intp(nindwn,xfracdwn,relspdwn,xx,downf(n))
 888       continue
c
      else
           read(5,*)     (downf(j),j=1,ndown)
      endif
c
c******************************************************************************
c******************************************************************************
c     upf(j), onf(j) and downf(j) are now set, print them out .
c
      write(6,1704) (upf(j),j=1,nup)
      write(6,1704) (onf(j),j=1,non)
      write(6,1704) (downf(j),j=1,ndown)
1704  format(1h ,12f5.2)
c
c******************************************************************************
c      read the change in the upstream and downstream extent of the grid
c      these multiply the original upstream and downstream
c      lengths of the grid by upext and dwnext.
c
      read(5,*)     upext,dwnext
      write(6,1704) upext,dwnext
c
c     end of input data for generating the new grid
c******************************************************************************
c
      jmnew  = nup + non + ndown
      jlenew = nup
      jtenew = nup + non
c
      if(jroths.gt.jmrow) jroths = jmrow
      if(jrothe.gt.jmrow) jrothe = jmrow
      if(jrotts.gt.jmrow) jrotts = jmrow
      if(jrotte.gt.jmrow) jrotte = jmrow
c
c******************************************************************************
c******************************************************************************
c   start the 'do 5000' loop over the number of input sections
c
      nloop = nsecs_in
      if(nsecs_in.eq.1) nloop = 2
c
      do 5000 k=1,nloop
c
c     calculate s coords. and put in array sold
c
      sold(1) = 0
      do 601 j=2,jmrow
      js = j -1 + j1
      sold(j) = sold(j-1) + sqrt((rsurf(js,k)-rsurf(js-1,k))**2
     &                         + (xsurf(js,k)-xsurf(js-1,k))**2)
601   continue
c
c      use snew to store s-coords of new grid lines
c
c     first calculate uptot, ontot, and downtot:
c
      uptot = 0
      do 605 j=1,nup-1
      uptot = uptot + upf(j)
605   continue
c
      ontot = 0
      do 607 j=1,non
      ontot = ontot + onf(j)
607   continue
c
      downtot = 0
      do 609 j=1,ndown
      downtot = downtot + downf(j)
609   continue
c
c     find the streamwise coordinate where the endwall rotations start and end.
c
      if(k.eq.1) then
      rhstart= sold(jroths)*1.0001
      rhend  = sold(jrothe)*0.9999
      endif
      if(k.eq.nsecs_in) then
      rtstart= sold(jrotts)*1.0001
      rtend  = sold(jrotte)*0.9999
      endif
c
c     find the streamwise coordinates of the new grid
c
      snew(1) = sold(jlerow)*(1.- upext)
      do 610 j=2,nup
      snew(j) = snew(j-1) + upf(j-1)/uptot*(sold(jlerow) - snew(1))
610   continue
c
      do 620 j=nup+1,nup+non
      snew(j) =snew(j-1) + onf(j-nup)/ontot*(sold(jterow)-sold(jlerow))
620   continue
c
      do 625 j=nup+non+1,nup+non+ndown
      snew(j) = snew(j-1) + downf(j-nup-non)/downtot*(sold(jmrow)
     &        - sold(jterow))*dwnext
625   continue
c
c     find the coordinates where endwall rotation starts and ends on the new grid
c
      do 626 j=1,jmnew
      if(k.eq.1) then
      if(snew(j).lt.rhstart) jroths = j
      if(snew(j).lt.rhend)   jrothe = j + 1
      endif
      if(k.eq.nsecs_in) then
      if(snew(j).lt.rtstart) jrotts = j
      if(snew(j).lt.rtend)   jrotte = j + 1
      endif
  626 continue
c
c******************************************************************************
c******************************************************************************
c     now find x,r,ysuct..,thick.. using 4 point interpolation
c
c     set up buffers xbuf and ybuf to contain the old values of  s  and x .
c
      do 632 jj=1,jmrow
      js = jj-1+j1
      xbuf(jj) = sold(jj)
      ybuf(jj) = xsurf(js,k)
c      write(6,*) 'jrow,jtot,xbuf,ybuf,snew ',
c     &            jj,js,xbuf(jj),ybuf(jj),snew(jj)
632   continue
c
c     find interpolated x-vals, and put them in xin:
c
      do 634 j=1,jmnew
      js = j-1+j1
      call intp(jmrow,xbuf,ybuf,snew(j),xsurf(js,k))
634   continue
c
c     find new r:
c
      do 636 jj=1,jmrow
      js = jj-1+j1
      ybuf(jj) = rsurf(js,k)
636   continue 
c
      do 638 j=1,jmnew
      js = j-1+j1
      call intp(jmrow,xbuf,ybuf,snew(j),rsurf(js,k))
638   continue
c
c     now  the blade upper surface  rt_upp .
c
c     first upstream of the leading edge
c
      do 640 jj=1,jlerow
      js = jj-1+j1
      ybuf(jj) = rt_upp(js,k)
640   continue
c
      do 642 j=1,jlenew
      js = j-1+j1
      call intp(jlerow,xbuf,ybuf,snew(j),ysnew(js))
642   continue
c
c    next on the blade row
c
      do 643 j=jlerow,jterow
      js =  j-1+j1
      jj =  j+1-jlerow
      ybuf(jj) = rt_upp(js,k)
      xbuf(jj) = sold(j)
643   continue
c    
      npoints= jterow - jlerow + 1
      do 644 j=jlenew,jtenew
      js = j-1+j1
      call intp(npoints,xbuf,ybuf,snew(j),ysnew(js))
644   continue
c
c  next downstream of the trailing edge
c
      do 645 j=jterow,jmrow
      js = j-1+j1
      jj =  j+1-jterow
      ybuf(jj) = rt_upp(js,k)
      xbuf(jj) = sold(j)
645   continue
c
      npoints= jmrow - jterow + 1
      do 646 j=jtenew,jmnew
      js = j-1+j1
      call intp(npoints,xbuf,ybuf,snew(j),ysnew(js))
646   continue
c
      do 670 j=1,jmnew
      js = j-1+j1
  670 rt_upp(js,k) = ysnew(js)
c
c     now the blade thickness  rt_thick
c     first on the blade surface
c
      do 647 j=jlerow,jterow
      js =  j-1+j1
      jj =  j+1-jlerow
      ybuf(jj) = rt_thick(js,k)
      xbuf(jj) = sold(j)
647   continue
c
      npoints= jterow - jlerow + 1
      do 648 j=jlenew,jtenew
      js = j-1+j1
      call intp(npoints,xbuf,ybuf,snew(j),rt_thick(js,k))
648   continue
c
c    set the thickness to zero upstream of the leading edge
      do 649 j=1,jlenew
      js = j-1+j1
      rt_thick(js,k) = 0.0
649   continue
c
c  set the thickness to zero downstream of the trailing edge
c  leaving 2 points for a cusp.
      do 651 j=jtenew+2,jmnew
      js = j-1+j1
      rt_thick(js,k)=0.0
651   continue
c
c     end of new grid generation on one stream surface.
c
      write(6,995) k
  995 format(//,' new interpolated grid points for k =',i5,//)
      write(6,997)
  997 format('         j     x            y           r          tk ')
      do 999 j=1,jmnew
      js = j - 1 + j1
  999 write(6,998) js,xsurf(js,k),rt_upp(js,k),rsurf(js,k),
     &             rt_thick(js,k)
  998 format( i10,4f12.5)
c
c******************************************************************************
c******************************************************************************
c   return to start a new stream surface
c
 5000 continue
c
c******************************************************************************
c******************************************************************************
c    reset the j indices for the new grid
c
      jmrow  =  jmnew
      jlerow = jlenew
      jterow = jtenew
      if(jroths.gt.jmnew) jroths = jmnew
      if(jrothe.gt.jmnew) jrothe = jmnew
      if(jrotts.gt.jmnew) jrotts = jmnew
      if(jrotte.gt.jmnew) jrotte = jmnew
c
      write(6,*) ' returning from subroutine newgrid '
      write(6,*)
      write(6,*) 'new values of jm,jle,jte relative to the start of this
     & blade row. ' 
      write(6,*)  jmrow,jlerow,jterow
      write(6,*)
      write(6,*) 'new values of jroths,jrothe,jrotts,jrotte relative to 
     &the start of this blade row.'
      write(6,*)  jroths,jrothe,jrotts,jrotte
c
c******************************************************************************
c******************************************************************************
c
      return
      end
c
c******************************************************************************
c******************************************************************************
c
      subroutine shroudflow
c
      include  'commall-open-19.2'
c
c
c     read in data for the shroud leakage model if ktips(nr) < 1 .
c
      rfshrd = 0.1
c
      if(nstep.eq.1) then
c
      do 50 nr = 1,nrows
c
      kshroud(nr) = 0
c
      if(ktips(nr).ge.0) go to 75
c 
      write(6,*)
      write(6,*)'shroud leakage data for row number', nr
      read(5,*)  kshroud(nr),jleaks(nr),jleake(nr),jlkins(nr),jlkine(nr)
      write(6,*)'kshroud, jleaks, jleake, jleakins, jleakine',
     &           kshroud(nr),jleaks(nr),jleake(nr),jlkins(nr),jlkine(nr) 
      read(5,*)  sealgap(nr),nseal(nr),cfshroud(nr),cfcasing(nr)
      write(6,*)'sealgap, nseal, cfshroud, cfcasing',
     &           sealgap(nr),nseal(nr),cfshroud(nr),cfcasing(nr)
      read(5,*)  wcase(nr),pitchin(nr)
      write(6,*)'rpm casing, inlet pitch angle', wcase(nr),pitchin(nr)
      write(6,*)
c
      sealgap(nr) = sealgap(nr)/sqrt(float(nseal(nr)))
      if(abs(pitchin(nr)).lt.10.0) pitchin(nr)= 10.0
      pitchin(nr) = pitchin(nr)*degrad
      wcase(nr)   = wcase(nr)*pi/30
c
   75 continue
c
      do 99 j=2,jm
   99 shrdflow(j) = 0.0
c
   50 continue
c
c     end of shroud leakage input data
c
c     set rfshrd = 1 on the first iteration only
      rfshrd = 1.
c
c    end of first iteration only part
c
      endif
c
c
      sworktot = 0.0
c
c     start do 100 loop over all blade rows
c
      do 100 nr =1,nrows
c
      if(ktips(nr).ge.0) go to 100
c
      js = jstart(nr)
      k  = kshroud(nr)
      j1 = jleaks(nr) + js -1
      j2 = jleake(nr) + js -1
      j3 = jlkins(nr) + js -1
      j4 = jlkine(nr) + js -1
      jshroud = 0.5*(j2 + j3)
      nbld    = nblade(jshroud)
c
c     find the area average pressure downstream of the seal
c
      poutavg = 0.0
      astot   = 0.0
      do 200 i = 1,imm1
      do 200 j = j3+1,j4
      asabs    = sqrt(asx(i,j,k)*asx(i,j,k) + asr(i,j,k)*asr(i,j,k))
      astot    = astot + asabs
      poutavg  = poutavg +
     &           asabs*(p(i,j,k)+p(i+1,j,k)+p(i,j-1,k)+p(i+1,j-1,k))
      as(i,j)  = asabs
  200 continue
      poutavg  = 0.25*poutavg/astot
      aleakin  = astot
c
c     find the average quantities for the flow leaking over the shroud
c     and calculate poax the po neglecting swirl velocity
c
      do 250 j=j1,j2
      do 250 i=1,im
      vmersq    = vx(i,j,k)*vx(i,j,k) + vr(i,j,k)*vr(i,j,k)
      eke       = vmersq + vt(i,j,k)*vt(i,j,k)
c
      if(ifgas.eq.0) then
           tstat     = (ho(i,j,k) - 0.5*eke)/cp
           torel     = tstat + 0.5*vmersq/cp
           poax(i,j) = p(i,j,k)*(torel/tstat)**rfga
      else
           hstag     = ho(i,j,k)
           hstat     = hstag - 0.5*eke
           horel     = hstat + 0.5*vmersq
           tstat     = tfromh(hstat,tref,href,ht1,ht2,ht3,ht4)
           torel     = tfromh(horel,tref,href,ht1,ht2,ht3,ht4)
           trat      = torel/tstat
           prat      = prat_from_trat(tstat,trat,alpha,beta1,beta2)
           poax(i,j) = p(i,j,k)*prat
      end if
c
      toax(i,j) = torel
  250 continue
c
      poaxavg = 0.0
      hoavg   = 0.0
      toaxavg = 0.0
      vtavg   = 0.0
      vxavg   = 0.0
      vravg   = 0.0
      astot   = 0.0
      do 300 i=1,imm1
      do 300 j=j1+1,j2
      asabs  = sqrt(asx(i,j,k)*asx(i,j,k) + asr(i,j,k)*asr(i,j,k))
      astot  = astot + asabs
      poaxavg = poaxavg +
     &        asabs*(poax(i,j)+poax(i+1,j)+poax(i,j-1)+poax(i+1,j-1))
      toaxavg = toaxavg +
     &        asabs*(toax(i,j)+toax(i+1,j)+toax(i,j-1)+toax(i+1,j-1))
      hoavg = hoavg +
     &        asabs*(ho(i,j,k)+ho(i+1,j,k)+ho(i,j-1,k)+ho(i+1,j-1,k))
      vtavg = vtavg +
     &        asabs*(vt(i,j,k)+vt(i+1,j,k)+vt(i,j-1,k)+vt(i+1,j-1,k))
      vravg = vravg +
     &        asabs*(vr(i,j,k)+vr(i+1,j,k)+vr(i,j-1,k)+vr(i+1,j-1,k))
      vxavg = vxavg +
     &        asabs*(vx(i,j,k)+vx(i+1,j,k)+vx(i,j-1,k)+vx(i+1,j-1,k))
      as(i,j) = asabs
 300  continue
c
      aleakout  = astot
      poaxavg   = 0.25*poaxavg/astot
      holeakout = 0.25*hoavg/astot 
      vtleakout = 0.25*vtavg/astot
      vxleakout = 0.25*vxavg/astot
      vrleakout = 0.25*vravg/astot
      toaxavg   = 0.25*toaxavg/astot

c
c     calculate the leakage mass flow rate.
c
      prat  = poutavg/poaxavg
      if(ifgas.eq.0) then
            tjet       = toaxavg*prat**fga
            if(tjet.gt.toaxavg)  tjet = 0.99999*toaxavg
            vjet       = sqrt(2*cp*(toaxavg - tjet))
      else
            trat     = trat_from_prat(toaxavg,prat,fgagas,r_alpha,
     &                 balpha1,balpha2)
            tjet     = trat*toaxavg
            hjet     = hfromt(tjet,tref,href,cp1,cp2,cp3)
            hojet    = hfromt(toaxavg,tref,href,cp1,cp2,cp3)
            if(hjet.gt.hojet) hjet = 0.99999*hojet
            vjet     = sqrt(2.0*(hojet-hjet))
      end if
c            
      rojet      = poutavg/(rgas*tjet)
      roshroud   = poutavg/(rgas*toaxavg)
      contract   = 0.6
      aleak      = 2*pi*r(jshroud,k)*sealgap(nr)*contract
      shroudmass = aleak*vjet*rojet
      shroudflbl = shroudmass/nbld
c
c  calculate the shroud surface area, assumed from le to te
c  the casing area is taken as this multilpled by a factor of 1.25
c
      rle   = r(jle(nr),k)
      xle   = x(jle(nr),k)
      rte   = r(jte(nr),k)
      xte   = x(jte(nr),k)
      difr  = rte - rle
      difx  = xte - xle
      dmer  = sqrt(difr*difr + difx*difx)
      sarea = pi*(rle+rte)*dmer
      sarea = sarea * 1.1
      carea = sarea * 1.25
      cfs   = cfshroud(nr)
      cfc   = cfcasing(nr)
      ravg  = 0.5*(rle+rte)
      vtshroud    = wrad(jshroud)*ravg
      vtcase      = wcase(nr)*ravg
      vtflow      = vtleakout
      shroudwork  = 0.0
      casework    = 0.0
      vlim        = 2.*amax1(abs(vtflow),abs(vtshroud),abs(vtcase))
c
      n_intervals = 20
c
      facs    =    0.5*cfs*sarea*roshroud/n_intervals
      facc    =    0.5*cfc*carea*roshroud/n_intervals
      do ns   = 1,n_intervals
      vtrelshrd   = vtflow - vtshroud
      vtrelcase   = vtflow - vtcase
      shroudforce = facs*vtrelshrd*abs(vtrelshrd)
      caseforce   = facc*vtrelcase*abs(vtrelcase)
      vtflow      = vtflow - (shroudforce + caseforce)/shroudmass
      if(abs(vtflow).gt.vlim) vtflow = vlim*vtflow/abs(vtflow)
      shroudwork  = shroudwork + shroudforce*vtshroud
      casework    = casework   + caseforce*vtcase
      end do
c
      shroudwork  = shroudwork + casework
      vtleakin    = vtflow
      holeakin    = holeakout  - shroudwork/shroudmass
      hlim = 1.25*holeakout
      if(holeakin.gt.hlim) holeakin = hlim
      sworktot    = sworktot   + shroudwork
c
c     note "sworktot" is the frictional work done on the shroud by
c     the leakage flow.
c
c   calculate the properties of the rentry flow
c
      rovnin =  shroudflbl/(0.25*aleakin)
      vnin   =  rovnin/roshroud
      if(vnin.gt.vjet)  vnin = vjet
      vsin   =  vnin/abs(tan(pitchin(nr)))
      if(k.eq.km) then
      vxin = vsin*difx/dmer + vnin*difr/dmer
      vrin = vsin*difr/dmer - vnin*difx/dmer
      endif
      if(k.eq.1) then
      vxin = vsin*difx/dmer - vnin*difr/dmer
      vrin = vsin*difr/dmer + vnin*difx/dmer
      endif
c
c    set the wall fluxes for leakage and re-entry flow
c
      fk  = 1.0
      rfshrd1 = 1.-rfshrd
      if(k.eq.1) fk = -1.0
      facgr    = fk*shroudflbl/aleakout
      rin      = 0.5*(r(j4,k) + r(j3,k))
      rout     = 0.5*(r(j1,k) + r(j2,k))
      rvtlkin  = rin*vtleakin
      rvtlkout = rout*vtleakout
      do 410 j = j1+1,j2
      sflow    = 0.0
      do 400 i = 1,imm1
      grout         = facgr*as(i,j)
      rfgrout       = rfshrd*grout
      sflow         = sflow + abs(grout)
      shroudgr(i,j) = rfshrd1*shroudgr(i,j)  + rfgrout
      shroudho(i,j) = rfshrd1*shroudho(i,j)  + rfgrout*holeakout
      shroudvx(i,j) = rfshrd1*shroudvx(i,j)  + rfgrout*vxleakout
      shroudrvt(i,j)= rfshrd1*shroudrvt(i,j) + rfgrout*rvtlkout
      shroudvr(i,j) = rfshrd1*shroudvr(i,j)  + rfgrout*vrleakout
  400 continue
      shrdflow(j)   = shrdflow(j-1) + sflow
  410 continue
c
c     set the shroud flow above/below the shrouded part of the blade.
c
      if(j3.gt.j2) then
      do 450 j = j2+1,j3
      shrdflow(j) = shroudflbl
  450 continue
      endif
      if(j4.lt.j1) then
      do 455 j = j4+1,j1
      shrdflow(j) = -shroudflbl
  455 continue
      endif
c
      facgr    = -fk*shroudflbl/aleakin
      do 510 j = j3+1,j4
      sflow    = 0.0
      do 500 i=1,imm1
      grout         = facgr*as(i,j)
      rfgrout       = rfshrd*grout
      sflow         = sflow + abs(grout)
      shroudgr(i,j) = rfshrd1*shroudgr(i,j)  + rfgrout
      shroudho(i,j) = rfshrd1*shroudho(i,j)  + rfgrout*holeakin
      shroudvx(i,j) = rfshrd1*shroudvx(i,j)  + rfgrout*vxin
      shroudrvt(i,j)= rfshrd1*shroudrvt(i,j) + rfgrout*rvtlkin
      shroudvr(i,j) = rfshrd1*shroudvr(i,j)  + rfgrout*vrin
  500 continue
      shrdflow(j)   = shrdflow(j-1) - sflow
  510 continue
c
      sleak(nr) = shroudmass
      swork(nr) = shroudwork
c
c    end of loop over all blade rows
c
  100 continue
c
      return
      end
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c
      subroutine shroudflux(sflux)
c
      include  'commall-open-19.2'
c
c
      dimension sflux(id,jd)
c
      do 9555 nr =1,nrows
      if(ktips(nr).ge.0) go to 9555
      js = jstart(nr)
      k  = kshroud(nr)
      j1 = jleaks(nr) + js -1
      j2 = jleake(nr) + js -1
      j3 = jlkins(nr) + js -1
      j4 = jlkine(nr) + js -1
      do 9554 j=j1,j2
      do 9554 i=1,imm1
 9554 rflux(i,j,k) = rflux(i,j,k) + sflux(i,j)
      do 9553 j=j3,j4
      do 9553 i=1,imm1
 9553 rflux(i,j,k) = rflux(i,j,k) + sflux(i,j)
 9555 continue
c
      return
      end
c
c******************************************************************************
c******************************************************************************
c
      subroutine coolin_1
c
c  this subroutine sets up fixed cooling flows which are only set once and not
c  changed during the iterations. the mach mumber of the ejected flow is specified
c  and the stagnation pressure of the coolant is only used to calculate the efficiency.
c
      include  'commall-open-19.2'
c
      write(6,*) ' entered coolin_1'
c
      do n = 1,nstages
           wpump(n) = 0.0
      end do
c
c*************************************************************************
c*************************************************************************
c
c      set initial cooling flows to zero and
c      work out area magnitudes for coolant flows
c
      do 400 j=2,jm
      do 400 k=1,kmm1
      cflowi1(j,k) = 0.0
      cflowim(j,k) = 0.0
      atoti1(j,k) = sqrt(abx(1,j,k)*abx(1,j,k) + abr(1,j,k)*abr(1,j,k)
     &                  + abt(j,k)*abt(j,k))
      atotim(j,k) = sqrt(abx(im,j,k)*abx(im,j,k)+abr(im,j,k)*abr(im,j,k)
     &                  + abt(j,k)*abt(j,k))
  400 continue
c
      do 500 j=2,jm
      do 500 i=1,imm1
      cflowk1(i,j) = 0.0
      cflowkm(i,j) = 0.0
      atotk1(i,j)= sqrt(asx(i,j,1)*asx(i,j,1)  + asr(i,j,1)*asr(i,j,1))
      atotkm(i,j)= sqrt(asx(i,j,km)*asx(i,j,km)+asr(i,j,km)*asr(i,j,km))
  500 continue
c
c***********************************************************************
c***********************************************************************
c
c      now evaluate the cooling mass flows through the blade surfaces.
c      and the associated energy and momentum fluxes.
c
c
      do 600 ncb = 1,ncoolb
c
      acoolb(ncb) = 0.0
      do 610 j = jcbs(ncb)+1,jcbe(ncb)
      do 610 k = kcbs(ncb),kcbe(ncb)-1
      if(ic(ncb).eq.1)  acoolb(ncb) = acoolb(ncb) + atoti1(j,k)
      if(ic(ncb).eq.im) acoolb(ncb) = acoolb(ncb) + atotim(j,k)
  610 continue
c
      do 620 j = jcbs(ncb)+1,jcbe(ncb)
c
             n_row    = nrow(j)
             n_stage  = 1 + (n_row-1)/2
c
      do 620 k = kcbs(ncb),kcbe(ncb)-1
c
      if(ic(ncb).eq.1)  then
           cflowin  = cflowb(ncb)*atoti1(j,k)/(acoolb(ncb)*nblade(j))
      else
           cflowin  = cflowb(ncb)*atotim(j,k)/(acoolb(ncb)*nblade(j))
      end if
c
      toin     = tocoolb(ncb)
      poin     = pocoolb(ncb)
      ravg     = ravg_cell(j,k)
      vblade   = wrad_coolb(ncb)*ravg
c  
c     pre swirl included 4/12/06  jdd
c
      w_preswirl     = wrad_coolb(ncb)*rvt_in_b(ncb)
      wpump(n_stage) = wpump(n_stage) + nblade(j)*
     &                 cflowin*(vblade*vblade - w_preswirl)
c
      if(ifgas.eq.0) then
           torel       = toin  + (0.5*vblade*vblade - w_preswirl)/cp
           po_rel      = poin*(torel/toin)**rfga
           toabs       = torel + 0.5*vblade*vblade/cp
           tin         = torel/(1. + 0.5*(ga-1.)*machcool(ncb)*
     &                   machcool(ncb))
           vin         = sqrt(2.*cp*(torel - tin))
           hoabs       = cp*toabs
      else
           delt        = toin-tref
           cpgas       = cp1 + cp2*delt + cp3*delt*delt
           gagas       = cpgas/(cpgas - rgas)
           rfgagas     = gagas/(gagas - 1.)
           ho_in       = hfromt(toin,tref,href,cp1,cp2,cp3)
           horel       = ho_in  + 0.5*vblade*vblade - w_preswirl
           hoabs       = horel + 0.5*vblade*vblade
           torel       = tfromh(horel,tref,href,ht1,ht2,ht3,ht4)
           po_rel      = poin*(torel/toin)**rfgagas
           tin         = torel/(1. + 0.5*(gagas-1.)*machcool(ncb)*
     &                   machcool(ncb))
           hstat       = hfromt(tin,tref,href,cp1,cp2,cp3) 
           if(hstat.gt.horel) hstat = 0.99999*horel
           vin         = sqrt(2.0*(horel - hstat))
      end if         
c
      torelb(ncb) = torel
      tstat_ejectb(ncb) = tin
      po_ejectb(ncb)    = po_rel
      vnin     = vin*sin(sangleb(ncb))
      vsin     = vin*cos(sangleb(ncb))
c
      if(ic(ncb).eq.1)  then      
           anorm    = atoti1(j,k)
           xnorm    = abx(1,j,k)/anorm
           tnorm    = abt(j,k)/anorm
           rnorm    = abr(1,j,k)/anorm
      else
           anorm    = atotim(j,k)
           xnorm    = abx(im,j,k)/anorm
           tnorm    = abt(j,k)/anorm
           rnorm    = abr(im,j,k)/anorm
      end if
c
      ang      =  xangleb(ncb)
      xtnorm   =  sqrt(xnorm*xnorm + tnorm*tnorm)
      sx       =  (tnorm*cos(ang) - rnorm*xnorm*sin(ang))/xtnorm
      st       = -(xnorm*cos(ang) + rnorm*tnorm*sin(ang))/xtnorm
      sr       =  sin(ang)*xtnorm
c
      if(ic(ncb).eq.1)  then
           xvel          =  vnin*xnorm + vsin*sx
           rvel          =  vnin*rnorm + vsin*sr
           tvel          =  vnin*tnorm + vsin*st
           cflowi1(j,k)  = cflowin
           hocwli1(j,k)  = cflowin*(hoabs + vblade*tvel)
           vxcwli1(j,k)  = cflowin*xvel
           vrcwli1(j,k)  = cflowin*rvel
           rvtcwli1(j,k) = cflowin*(vblade  +  tvel)*ravg
      else
           xvel          = -vnin*xnorm + vsin*sx
           rvel          = -vnin*rnorm + vsin*sr
           tvel          = -vnin*tnorm + vsin*st
           cflowim(j,k)  = cflowin
           hocwlim(j,k)  = cflowin*(hoabs + vblade*tvel)
           vxcwlim(j,k)  = cflowin*xvel
           vrcwlim(j,k)  = cflowin*rvel
           rvtcwlim(j,k) = cflowin*(vblade  +  tvel)*ravg
      end if
c
  620 continue
c
c     end of blade surface cooling setup
c
  600 continue
c
c**************************************************************************
c**************************************************************************
c
c      now evaluate the cooling mass flows through the hub and casing.
c
      do 700 ncw = 1,ncoolw
c
      acoolw(ncw) = 0.0
      do 710 j =  jcws(ncw)+1,jcwe(ncw)
      do 710 i =  icws(ncw),icwe(ncw)-1
      if(kc(ncw).eq.1)  acoolw(ncw) = acoolw(ncw) + atotk1(i,j)
      if(kc(ncw).eq.km) acoolw(ncw) = acoolw(ncw) + atotkm(i,j)
  710 continue
c
      do 720 j = jcws(ncw)+1,jcwe(ncw)
c
           n_row    = nrow(j)
           n_stage  = nstage(n_row)
c
      do 720 i = icws(ncw),icwe(ncw)-1
c
      if(kc(ncw).eq.1) then
           cflowin  = cfloww(ncw)*atotk1(i,j)/(acoolw(ncw)*nblade(j))
           ravg     = 0.5*(r(j,1)+r(j-1,1))
      else
           cflowin  = cfloww(ncw)*atotkm(i,j)/(acoolw(ncw)*nblade(j))
           ravg     = 0.5*(r(j,km)+r(j-1,km))
      end if 
c
      toin   = tocoolw(ncw)
      poin   = pocoolw(ncw)
      vwall  = wrad_coolw(ncw)*ravg
c
c     pre swirl included 4/12/06  jdd
c
      w_preswirl     = wrad_coolw(ncw)*rvt_in_w(ncw)
      wpump(n_stage) = wpump(n_stage) + nblade(j)*
     &                 cflowin*(vwall*vwall - w_preswirl)
c
      if(ifgas.eq.0) then
           torel  = toin  + (0.5*vwall*vwall - w_preswirl)/cp
           po_rel = poin*(torel/toin)**rfga
           toabs  = torel + 0.5*vwall*vwall/cp
           tin    = torel/(1. + 0.5*(ga-1.)*machcool(ncw)*machcool(ncw))
           vin    = sqrt(2.*cp*(torel - tin))
           hoabs  = cp*toabs
      else
           delt        = toin-tref
           cpgas       = cp1 + cp2*delt + cp3*delt*delt
           gagas       = cpgas/(cpgas - rgas)
           rfgagas     = gagas/(gagas - 1.)
           ho_in       = hfromt(toin,tref,href,cp1,cp2,cp3)
           horel       = ho_in  + 0.5*vwall*vwall - w_preswirl
           hoabs       = horel + 0.5*vwall*vwall
           torel       = tfromh(horel,tref,href,ht1,ht2,ht3,ht4)
           po_rel      = poin*(torel/toin)**rfgagas
           tin         = torel/(1. + 0.5*(gagas-1.)*machcool(ncb)*
     &                   machcool(ncb)) 
           hstat       = hfromt(tin,tref,href,cp1,cp2,cp3) 
           if(hstat.gt.horel) hstat = 0.99999*horel
           vin         = sqrt(2.0*(horel - hstat))
      end if         
c
      torelw(ncw) = torel
      tstat_ejectw(ncw) = tin
      po_ejectw(ncw)    = po_rel
      vnin   = vin*sin(sanglew(ncw))
      vsin   = vin*cos(sanglew(ncw))
c
      if(kc(ncw).eq.1) then
           anorm  = atotk1(i,j)
           xnorm  = asx(i,j,1)/anorm
           rnorm  = asr(i,j,1)/anorm
      else
           anorm  = atotkm(i,j)
           xnorm  = asx(i,j,km)/anorm
           rnorm  = asr(i,j,km)/anorm
      end if
c
      ang    =  tanglew(ncw)
      xrnorm =  sqrt(xnorm*xnorm + rnorm*rnorm)
      sx     =  cos(ang)*rnorm/xrnorm
      st     =  sin(ang)
      sr     = -cos(ang)*xnorm/xrnorm

c
      if(kc(ncw).eq.1) then
           xvel   =  vnin*xnorm + vsin*sx
           rvel   =  vnin*rnorm + vsin*sr
           tvel   =  vsin*st
      else
           xvel   =  -vnin*xnorm + vsin*sx
           rvel   =  -vnin*rnorm + vsin*sr
           tvel   =   vsin*st
      end if
c
      if(kc(ncw).eq.1) then    
           cflowk1(i,j)  = cflowin
           hocwlk1(i,j)  = cflowin*(hoabs + vwall*tvel)
           vxcwlk1(i,j)  = cflowin*xvel
           vrcwlk1(i,j)  = cflowin*rvel
           rvtcwlk1(i,j) = cflowin*(vwall + tvel)*ravg
      else
           cflowkm(i,j)  = cflowin
           hocwlkm(i,j)  = cflowin*(hoabs + vwall*tvel)
           vxcwlkm(i,j)  = cflowin*xvel
           vrcwlkm(i,j)  = cflowin*rvel
           rvtcwlkm(i,j) = cflowin*(vwall + tvel)*ravg
      endif
c
  720 continue
c
  700 continue
c
c     find the sum of the cooling flows added up to each j station
c
      sumcwl(1) = 0.0
      do 20 j=2,jm
      cooladd   = 0.0
      do 25 k=1,kmm1
   25 cooladd   = cooladd  + (cflowi1(j,k) + cflowim(j,k))*nblade(j)
      do 30 i=1,imm1
   30 cooladd   = cooladd  + (cflowk1(i,j) + cflowkm(i,j))*nblade(j)
      sumcwl(j) = sumcwl(j-1) + cooladd
   20 continue
c
c     end of cooling flow setup
c
      return
      end
c
c**********************************************************************
c
      subroutine bleedout(ibleed,kbleed)
c
c
      include  'commall-open-19.2'
c
c
      dimension 
     & kbleeds(nbled),kbleede(nbled),jbleeds(nbled),jbleede(nbled),
     & ibleeds(nbled),ibleede(nbled),mbleed(nbled)
c
      real massbled,mbleed
c
   99 format(a72)
      nbleedtot = 0
      ibleed    = 0
      kbleed    = 0
c
      do 100 n_row = 1,nrows
c
      j1 = jstart(n_row)
      read(5,99)  dummy_input
      write(6,99) dummy_input
      read(5,*)   nbleed
c
c    read in the bleed flow details. the j values in the input are
c    defined for the current blade row. ie j = 1 at the upstream mixing plane.
c
      do 10 nbld =1,nbleed
      nbleedtot = nbleedtot + 1
      read(5,*) iblds,iblde,jblds,jblde,kblds,kblde,massbled
      kbleeds(nbleedtot)  = kblds
      kbleede(nbleedtot)  = kblde
      ibleeds(nbleedtot)  = iblds
      ibleede(nbleedtot)  = iblde
      jbleeds(nbleedtot)  = j1 + jblds
      jbleede(nbleedtot)  = j1 + jblde
      mbleed(nbleedtot)   = massbled
   10 continue
c
  100 continue
c
c      work out area magnitudes for bleed flows
c
      do 150 j=2,jm
      do 150 k=1,kmm1
      blflowi1(j,k) = 0.0
      blflowim(j,k) = 0.0
      atoti1(j,k) = sqrt(abx(1,j,k)*abx(1,j,k) + abr(1,j,k)*abr(1,j,k)
     &                 + abt(j,k)*abt(j,k))
      atotim(j,k) = sqrt(abx(im,j,k)*abx(im,j,k)+abr(im,j,k)*abr(im,j,k)
     &                 + abt(j,k)*abt(j,k))
  150 continue
c
      do 500 j=2,jm
      do 500 i=1,imm1
      blflowk1(i,j) = 0.0
      blflowkm(i,j) = 0.0
      atotk1(i,j) =sqrt(asx(i,j,1)*asx(i,j,1) + asr(i,j,1)*asr(i,j,1))
      atotkm(i,j) =sqrt(asx(i,j,km)*asx(i,j,km)+asr(i,j,km)*asr(i,j,km))
  500 continue
c
c
      do 200 nbld =1,nbleedtot
c
      ks = kbleeds(nbld)
      ke = kbleede(nbld)
      js = jbleeds(nbld)
      je = jbleede(nbld)
      is = ibleeds(nbld)
      ie = ibleede(nbld)
c
      if(ks.eq.ke.and. (ks.ne.0) ) then
c
c     if hub or casing bleed
c
      k = ks
      kbleed = 1
      ableed = 0.0
      do 210 j = js+1,je
      do 210 i=is,ie-1
      if(k.eq.1)  ableed = ableed + atotk1(i,j)
      if(k.eq.km) ableed = ableed + atotkm(i,j)
  210 continue
      do 300 j = js+1,je
      nblad = nblade(j)
      do 300 i = is,ie-1
      if(k.eq.1)   blflowk1(i,j) = mbleed(nbld)*atotk1(i,j)/ableed/nblad
      if(k.eq.km)  blflowkm(i,j) = mbleed(nbld)*atotkm(i,j)/ableed/nblad
  300 continue
c
      endif
c
      if(is.eq.ie.and. (is.ne.0) ) then
c
c     if blade surface bleed
c
      i = is
      ibleed = 1
      ableed = 0.0
      do 250 j = js+1,je
      do 250 k = ks,ke-1
      if(i.eq.1)  ableed = ableed + atoti1(j,k)
      if(i.eq.im) ableed = ableed + atotim(j,k)
  250 continue
      do 270 j = js+1,je
      nblad = nblade(j)
      do 270 k = ks,ke-1
      if(i.eq.1)  blflowi1(j,k) = mbleed(nbld)*atoti1(j,k)/ableed/nblad
      if(i.eq.im) blflowim(j,k) = mbleed(nbld)*atotim(j,k)/ableed/nblad
  270 continue
c
      endif
c
c     obtain the total bled flow for use in the continuity check
c
      sumbleed(1) = 0.0
      do 400 j=2,jm
      nblad  = nblade(j)
      bleedj = 0.0
      do 410 i=1,imm1
  410 bleedj       = bleedj + nblad*(blflowk1(i,j) + blflowkm(i,j))
      do 420 k=1,kmm1
  420 bleedj       = bleedj + nblad*(blflowi1(j,k) + blflowim(j,k))
      sumbleed(j)  = sumbleed(j-1) + bleedj
  400 continue
c
c
  200 continue
c
      return
      end
c
c**********************************************************************
c
      subroutine eficool(hbledtot)
c
c     this subroutine calculates the mass averaged quantities and machine
c     efficiency and writes them out to unit 6.
c
      include  'commall-open-19.2'
c
c
      dimension blade_flow(jd),sum_entpy(jd),sumpo(jd),sumto(jd),
     &          sumrvt(jd),sumtstat(jd), eta_loss(jd),sumpstat(jd),
     &          check_flow(jd)
c
c      sum fluxes of mass,stagnation presure  etc every 100 time steps
c
c      set temp1, temp2, temp3 ,temp4 ,store2 to stagn temp, stagn press, rvt, entropy,
c      and static temperature. for use in subroutine sumflx.
c
c    temp3  stores rvt
c    temp2  stores the stagnation pressure ratio
c    temp1  stores the stagnation temperature
c    temp4 stores the entropy.
c    store2 stores the static temperature
c
      enthin = 0.0
c
      do 5610 k=1,km
      do 5610 j=1,jm
      do 5610 i=1,im
      eke  = 0.5*(vx(i,j,k)*vx(i,j,k) + vt(i,j,k)*vt(i,j,k)
     &          + vr(i,j,k)*vr(i,j,k))
c
      if(ifgas.eq.0) then
           ts    = (ho(i,j,k) - eke)*rcp
           tstag = ho(i,j,k)*rcp
      else
           hstat  = ho(i,j,k) - eke
           ts     = tfromh(hstat,tref,href,ht1,ht2,ht3,ht4)
           hstag  = ho(i,j,k)
           tstag  = tfromh(hstag,tref,href,ht1,ht2,ht3,ht4)
      end if
c
      if(tstag.lt.0.1*to1(kmid))  tstag = 0.1*to1(kmid)
c
      if(ts.lt.0.1*to1(kmid))      ts   = 0.1*to1(kmid)
c
      pstat = p(i,j,k)
      if(pstat.lt.0.02*po1(kmid)) pstat=0.02*po1(kmid)
c
      if(ifgas.eq.0) then
            pstag = pstat*(tstag/ts)**rfga
      else
            trat  = tstag/ts
            prat  = prat_from_trat(ts,trat,alpha,beta1,beta2)
            pstag = pstat*prat
      end if
c
      trat    = ts/to1(kmid)
      prat    = pstat/po1(kmid)
      entrpy  = 0.0
      if(prat.gt.0.99999.and.prat.lt.1.00001) go to 5605
      if(trat.gt.0.99999.and.trat.lt.1.00001) go to 5605
c
      if(ifgas.eq.0) then
           entrpy = cp*log(trat)- rgas*log(prat)
      else
           prat   = prat_from_trat(to1(kmid),trat,alpha,beta1,beta2)
           entrpy = pstat/(po1(kmid)*prat)
      end if
c
 5605 continue
c
      store2(i,j,k)= ts
      temp1(i,j,k) = tstag
      temp2(i,j,k) = pstag/po1(kmid)
      temp3(i,j,k) = rorvt(i,j,k)/ro(i,j,k)
      temp4(i,j,k) = entrpy
c
 5610 continue
c
      call sumflx(blade_flow,sumpo,sumto,sumrvt,sum_entpy,sumtstat,
     &            sumpstat) 
c
c******************************************************************************
c
      if(ifend.eq.1)  then
c
c     write output for globplot only at end of whole calculation.
c
      write(11)   1,1
      write(11)   (blade_flow(j), j=1,jm)
      write(11)   (sumpo(j), j=1,jm)
      write(11)   (sumto(j), j=1,jm)
      write(11)   (sumrvt(j),j=1,jm)
      write(11)   (sum_entpy(j),j=1,jm)

      endif
c********************************************************
c  mass average the 1d values
c
      do 5611 j=1,jm
      tflow         = blade_flow(j)
      sumpo(j)      = sumpo(j)/tflow
      sumto(j)      = sumto(j)/tflow
      sumrvt(j)     = sumrvt(j)/tflow
      sum_entpy(j)  = sum_entpy(j)/tflow
      sumtstat(j)   = sumtstat(j)/tflow
      sumpstat(j)   = sumpstat(j)/tflow
 5611 continue
c
      dtref  = abs( sumto(1) - sumto(jm) )
      ifout  = 0
      if(dtref.lt.0.005*sumto(1)) then
           dtin = sumto(1)  - sumtstat(1)
           dtex = sumto(jm) - sumtstat(jm)
           dtref = amax1(dtin,dtex)
           ifout = 1
      end if
      dhref = cp*dtref
      do j = 1,jm
      eta_loss(j) = sumto(jm)*(sum_entpy(j)-sum_entpy(1))/dhref
      end do
c
c ********************************************************
c    write a file to plot the entropy loss coefficient or lost efficiency'
c  
      if(ifend.eq.1) then   
      open(unit=23,file='loss-co.plt')
      write(23,*) ' plotting output for lost efficiency '
      write(23,*) ' number of output points ', jm
      write(23,*) ' meridional distance '
      write(23,5612) (smerid(j,kmid), j=1,jm)
      write(23,*) ' lost efficiency '
      write(23,5612) (eta_loss(j),   j=1,jm)
 5612 format(10f10.5)
      end if
c
c******************************************************************************
c
c     write out the main output to unit 6
c
      if(ifout.eq.0) then
           write(6,5645)
           write(6,5646)
      else
           write(6,5647)
           write(6,5648)
      end if
c
 5645 format( '    j           mass          local/inlet         lost
     &        stagnation      stagnation       mass avg     ')
 5646 format( '  value      flow. kg/s          flow          efficiency
     &        pressure        temperature        r*vt       ')
c
 5647 format( '    j           mass          local/inlet       entropy
     &        stagnation      stagnation       mass avg     ')
 5648 format( '  value      flow. kg/s          flow         loss coeff.
     &        pressure        temperature        r*vt       ')
c
c
      do 5640  j=1,jm
      check_flow(j) =  blade_flow(j) + shrdflow(j)*nblade(j)
     &               +  sumbleed(j) - sumcwl(j)
      ratio          =  check_flow(j)/check_flow(1)
      write(6,5630) j,check_flow(j),ratio,eta_loss(j),sumpo(j),
     &                sumto(j),sumrvt(j)
 5630 format( i5, 6f17.5)
 5640 continue
c
c************************************************************************************
c
c     write out the reynolds number for each blade row
c
      write(6,*)
      write(6,*)
      do nrw = 1,nrows
           jledge  = jle(nrw)
           jtedge  = jte(nrw)
           rowref = sqrt(rovx(imid,jtedge,kmid)*rovx(imid,jtedge,kmid) +
     &                  rovr(imid,jtedge,kmid)*rovr(imid,jtedge,kmid) 
     &                + rowt(imid,jtedge,kmid)*rowt(imid,jtedge,kmid)  )
           reynolds(nrw) = chord(nrw)*rowref/viscosy(nrw)
      write(6,*) ' blade row number ',nrw
      write(6,*) ' chord = ',chord(nrw),' rov-exit =',rowref,
     &           ' viscosity = ',viscosy(nrw)
      write(6,*) ' reynolds number based on blade exit conditions = ',
     &             reynolds(nrw)
      end do
      write(6,*)
c
c***************************************************************************************************
c***************************************************************************************************
c
c     work out the efficiency of the whole machine if there are cooling flows present
c
      if(ifcool.ne.0) then
c
c     find the total work done in pumping the coolant 
c
      wpump_tot = 0.0
      do 145 n  = 1,nstages
      wpump_tot = wpump_tot + wpump(n)
  145 continue
      wpump(nstages+1) = wpump_tot
c
      nstgp1  = nstages + 1
      if(nstages.eq.1) nstgp1 = 1
c 
c***********************************************************************************************
c    now write out the stage performances.
c    note:  the nstgp1 stage is the whole machine
c
      do 3333  nstg = 1,nstgp1
c
      write(6,*)
      write(6,*) '******************************************************
     &******'
      if(nstg.eq.nstgp1) then
           write(6,*) ' the following output is for the whole machine '
      else
           write(6,*) ' the following output is for stage number ',nstg
      end if
      write(6,*) '******************************************************
     &******'
c
      if(nstg.le.nstages) then
              j1      = jstg_start(nstg)
              j2      = jstg_end(nstg) + 1
              j2m1    = j2
c
c     jdd change  2/7/07  to include the mixing loss in the efficiency of the upstream blade row.
c     for the last stage there is no mixing downstream. base the effy on the last but one point.
c
              if(nstg.eq.nstages) j2m1 = jm -1 
c
c     end of jdd addition
c
c     jmix_stage is only used to find the stage reaction.
              n_row      = nrow(j1+1)
              jmix_stage = jmix(n_row) 
      else
c
c   for the whole machine there is no mixing downstream. base the effy on the last but one point, j2m1.
c
              j1   = 2
              j2   = jm
              j2m1 = jm - 1
              jmix_stage = (j1+j2)/2
      end if
c
      write(6,*)
      write(6,*) ' efficiencies based on mass averaged values at j = ',
     &             j1,j2m1
      write(6,*)
c
c     find the stage inlet and outlet stagnation conditions
c
      poin    = sumpo(j1)*po1(kmid)
      po2     = sumpo(j2m1)*po1(kmid)
      porat   = poin/po2
      to_in   = sumto(j1)
      to_out  = sumto(j2m1)
      psin    = sumpstat(j1)
      psout   = sumpstat(j2m1)
c
c     if there are any cooling flows, sum the cooling mass flows and enthalpies
c     and find the isentropic work available from the cooling flows.
c
      enthin        = 0.0
      wcool_is      = 0.0
      wcool_at_hole = 0.0
      cflow_stg     = 0.0
c
      do 150 ncb = 1,ncoolb
c
      sum_pstat  = 0.0
      sum_mas    = 0.0

      if(jcbs(ncb).lt.j1.and.jcbe(ncb).le.j1) go to 150
      if(jcbs(ncb).ge.j2.and.jcbe(ncb).gt.j2) go to 150

      do 151 j = jcbs(ncb)+1,jcbe(ncb)

      if(j.le.j1.or.j.gt.j2) go to 151

      do 152 k = kcbs(ncb),kcbe(ncb)-1
      if(ic(ncb).eq.1) then
           sum_pstat  = sum_pstat + cflowi1(j,k)*
     &     (p(1,j-1,k)+p(1,j,k)+p(1,j-1,k+1)+p(1,j,k+1))
           sum_mas = sum_mas + cflowi1(j,k)
      else
           sum_pstat  = sum_pstat + cflowim(j,k)*
     &     (p(im,j-1,k)+p(im,j,k)+p(im,j-1,k+1)+p(im,j,k+1))
           sum_mas = sum_mas + cflowim(j,k)
      endif
  152 continue

  151 continue
c
c     "ps_eject"  is the average static pressure on the surface of the cooling patch.
c     "po_at_hole"  is the absolute stagnation pressure at which the cooling 
c     flow would need to be supplied to the blade row if there was no loss in
c     the internal cooling passages.
c
      ps_eject  = 0.25*sum_pstat/sum_mas
c
      if(ifgas.eq.0) then
           gagas    = ga
           cpgas    = cp
      else
           delt   = torelb(ncb) - tref
           cpgas  = cp1 + cp2*delt + cp3*delt*delt
           gagas  = cpgas/(cpgas - rgas)
      end if
c
c      ts_eject  = torelb(ncb)
c     &          /(1.+.5*(gagas-1)*machcool(ncb)*machcool(ncb))
c     jdd changed to use stored value of ts_eject. may 2018 .
      ts_eject = tstat_ejectb(ncb)
c
      if(ifgas.eq.0) then
           po_at_hole= ps_eject*(tocoolb(ncb)/ts_eject)**rfga
      else
           trat = tocoolb(ncb)/ts_eject
           prat = prat_from_trat(ts_eject,trat,alpha,beta1,beta2)
           po_at_hole= ps_eject*prat
      end if
c
      pocool     = pocoolb(ncb)
      tocool     = tocoolb(ncb)
      sum_mas    = sum_mas*nblade(jcbs(ncb)+1)
      cflow_stg  = cflow_stg + sum_mas
c
      if(ifgas.eq.0) then
           to2_is         = tocool*(po2/pocool)**fga
           wcool_is       = wcool_is + sum_mas*cp*(tocool-to2_is)
           to2is_at_hole  = tocool*(po2/po_at_hole)**fga
           wcool_at_hole  = wcool_at_hole 
     &                    + sum_mas*cp*(tocool - to2is_at_hole)
           enthin         = enthin   + sum_mas*cp*tocool
      else
           prat       = po2/pocool
           trat       = trat_from_prat(tocool,prat,fgagas,r_alpha,
     &                  balpha1,balpha2)
           to2_is     = tocool*trat
           hocool_is  = hfromt(to2_is,tref,href,cp1,cp2,cp3)
           hocool     = hfromt(tocool,tref,href,cp1,cp2,cp3)
           wcool_is   = wcool_is + sum_mas*(hocool - hocool_is)
           prat       = po2/po_at_hole
           trat       = trat_from_prat(tocool,prat,fgagas,r_alpha,
     &                  balpha1,balpha2)
           to2is_at_hole = tocool*trat
           ho2is_at_hole = hfromt(to2is_at_hole,tref,href,cp1,cp2,cp3)
           wcool_at_hole = wcool_at_hole 
     &                   + sum_mas*(hocool - ho2is_at_hole)
           enthin        = enthin + sum_mas*hocool
      end if
c
      write(6,164) ncb,pocool,po_ejectb(ncb),po_at_hole
  164 format(' blade cooling patch number',i3,' specified supply stagnat 
     &ion pressure = ',t67,f10.1,/,' isentropic stagnation pressure with
     &in the blade = ',t67,f10.1,/,' coolant stagnation pressure at cool
     &ing hole exit= ',t67,f10.1)  
      write(6,*)
c
  150 continue
c
      do 160 ncw = 1,ncoolw
c
      if(jcws(ncw).lt.j1.and.jcwe(ncw).le.j1) go to 160
      if(jcws(ncw).ge.j2.and.jcwe(ncw).gt.j2) go to 160
c
      sum_pstat  = 0.0
      sum_mas    = 0.0
      do 161 j = jcws(ncw)+1,jcwe(ncw)

      if(j.le.j1.or.j.gt.j2) go to 161

      do 162 i=icws(ncw),icwe(ncw)-1    
      if(kc(ncw).eq.1) then
           sum_pstat  = sum_pstat + cflowk1(i,j)*
     &     (p(i,j-1,1)+p(i,j,1)+p(i+1,j-1,1)+p(i+1,j,1))
           sum_mas = sum_mas + cflowk1(i,j)
      else
           sum_pstat  = sum_pstat + cflowkm(i,j)*
     &     (p(i,j-1,km)+p(i,j,km)+p(i,j-1,km)+p(im,j,km))
           sum_mas = sum_mas + cflowkm(i,j)
      endif
  162 continue

  161 continue
c
c     "ps_eject"  is the average static pressure on the surface of the cooling patch.
c     "po_at_hole"  is the absolute stagnation pressure at which the cooling 
c     flow would need to be supplied to the blade row if there was no loss in
c     the internal cooling passages.
c
      ps_eject  = 0.25*sum_pstat/sum_mas
c      ts_eject  = torelw(ncw)/(1.+.5*(ga-1)*machcool(ncw)*machcool(ncw))
c  jdd changed to use stored value of ts_eject. may 2018 .
      ts_eject  = tstat_ejectw(ncw)
c
      if(ifgas.eq.0) then
           po_at_hole = ps_eject*(tocoolw(ncw)/ts_eject)**rfga
      else
           trat       = tocoolw(ncw)/ts_eject
           prat       = prat_from_trat(ts_eject,trat,alpha,beta1,beta2)
           po_at_hole = ps_eject*prat
      end if
c
      pocool     = pocoolw(ncw)
      tocool     = tocoolw(ncw)
      sum_mas    = sum_mas*nblade(jcws(ncw)+1)
      cflow_stg  = cflow_stg + sum_mas
c
c      if(ifcool.eq.2) po_at_hole = pocool
c
      if(ifgas.eq.0) then
           to2_is        = tocool*(po2/pocool)**fga
           wcool_is      = wcool_is + sum_mas*cp*(tocool - to2_is)
           to2is_at_hole = tocool*(po2/po_at_hole)**fga
           wcool_at_hole = wcool_at_hole 
     &                   + sum_mas*cp*(tocool - to2is_at_hole)
           enthin        = enthin   + sum_mas*cp*tocool
      else
           prat          = po2/pocool
           trat          = trat_from_prat(tocool,prat,fgagas,r_alpha,
     &                     balpha1,balpha2)
           to2_is        = tocool*trat
           hocool_is     = hfromt(to2_is,tref,href,cp1,cp2,cp3)
           hocool        = hfromt(tocool,tref,href,cp1,cp2,cp3)
           wcool_is      = wcool_is + sum_mas*(hocool - hocool_is)
           prat          = po2/po_at_hole
           trat          = trat_from_prat(tocool,prat,fgagas,r_alpha,
     &                     balpha1,balpha2)
           to2is_at_hole = tocool*trat
           ho2is_at_hole = hfromt(to2is_at_hole,tref,href,cp1,cp2,cp3)
           wcool_at_hole = wcool_at_hole 
     &                  + sum_mas*(hocool - ho2is_at_hole)
           enthin       = enthin + sum_mas*hocool
      end if
c
      write(6,163) ncw,pocool,po_ejectw(ncw),po_at_hole
  163 format(' wall cooling patch number',i3,' specified supply stagnati 
     &on pressure  = ',t67,f10.1,/,' isentropic stagnation pressure with
     &in the blade = ',t67,f10.1,/,' coolant stagnation pressure at cool
     &ing hole exit= ',t67,f10.1)  
      write(6,*)
c
  160 continue
c
c
c******************************************************************************
c******************************************************************************
c
c     mod 9/5/2002. to allow for errors in continuity it is better to base the actual
c     work on the inlet flow, added cooling flow, and the actual temperature change.
c     otherwise small continuity errors cause significant efficiency errors.
c
      flow_out =   blade_flow(j1) +  cflow_stg
c
      if(ifgas.eq.0) then
           wnet     = cp*(blade_flow(j1)*to_in - flow_out*to_out)
     &              + enthin
           to2_is   = to_in*(po2/poin)**fga
           wmain_is = blade_flow(j1)*cp*(to_in - to2_is)
      else
           ho_in    = hfromt(to_in,tref,href,cp1,cp2,cp3)
           ho_out   = hfromt(to_out,tref,href,cp1,cp2,cp3)
           prat     = po2/poin
           trat     = trat_from_prat(to_in,prat,fgagas,r_alpha,
     &                balpha1,balpha2)
           to2_is   = to_in*trat
           ho2_is   = hfromt(to2_is,tref,href,cp1,cp2,cp3)
           wmain_is = blade_flow(j1)*(ho_in - ho2_is)
           wnet     = blade_flow(j1)*ho_in  + enthin - flow_out*ho_out
      end if
c
      wtot_is  =   wmain_is + wcool_is
      pccool   =   100.*cflow_stg/blade_flow(j1)
c
      if(wnet.ge.0.0) etatt      = wnet/wtot_is
      if(wnet.ge.0.0) etatt_nop  = (wnet + wpump(nstg))/wtot_is
      if(wnet.lt.0.0) etatt      = wtot_is/wnet
      if(wnet.lt.0.0) etatt_nop  = wtot_is/(wnet + wpump(nstg))
c
      wtot_at_hole = wmain_is + wcool_at_hole
      if(wnet.ge.0.0) etatt_at_hole =  wnet/wtot_at_hole
      if(wnet.ge.0.0) eta_ath_nop   =  (wnet + wpump(nstg))/wtot_at_hole
      if(wnet.lt.0.0) etatt_at_hole =  wtot_at_hole/wnet
      if(wnet.lt.0.0) eta_ath_nop   =  wtot_at_hole/(wnet + wpump(nstg))
c
c     write the ouput for the stage or the whole machine
c
      write(6,*)
      if(wnet.gt.0.0) then
           write(6,*) ' this stage is a turbine '
           write(6,*)
           write(6,*) ' pressure ratio, poinlet/poexit             = ',
     &                  porat 
      else
           write(6,*) ' this stage is a compressor '
           write(6,*)
           write(6,*) ' stagnation pressure ratio, poexit/poinlet  = ',
     &                  1.0/porat 
      end if
c
c     work out the gas properties at stage inlet
c
      if(ifgas.ne.0) then
           delt  = to_in - tref
           cp = cp1 + cp2*delt + cp3*delt*delt
           ga = cp/(cp - rgas)
      end if
c
c
      write(6,*)      ' gas properties at stage inlet, cp = ',cp,
     &                ' gamma = ',ga,' r = ',rgas
      write(6,*)      ' inlet and exit stagnation pressures        = ',
     & poin, po2
      write(6,*)      ' inlet and exit static pressures            = ',
     &                  psin, psout
      write(6,*)      ' inlet and exit stagnation temperatures     = ',
     &                  to_in, to_out
      write(6,*)      ' inlet and outlet mass flow rates           = ',
     &                  blade_flow(j1), flow_out, ' kg/sec .'           
      write(6,*)      ' net power output                           = ',
     &                  wnet/1000.,' kw. this includes any work done on 
     &the coolant by the rotor- i.e. the pumping work.'
      write(6,*)      ' neglecting pumping power, net power output = ', 
     &                  (wnet + wpump(nstg))/1000. ,' kw.'
      write(6,*)      ' coolant pumping power                      = ',
     &                  wpump(nstg)/1000., ' kw.'
      write(6,*)      ' total coolant flow added                   = ',
     &                  cflow_stg, ' kg/sec.'
      write(6,*)      ' percentage coolant flow added              = ',
     &                  pccool
c
      if(nstg.le.nstages) then
           reaction = (sumtstat(jmix_stage) - sumtstat(j2m1))
     &               /(sumtstat(j1)         - sumtstat(j2m1))
           if(wnet.lt.0.0) reaction = 1.0 - reaction
           write(6,*) ' mean stage reaction                        = ',
     &     reaction
           write(6,*)
      end if
c
      write(6,*) '*****************************************************'
      write(6,*)
      write(6,*) ' using the values of coolant stagnation pressure input
     & as data'
      write(6,*) '- then:-  '
      write(6,*)
      write(6,*) ' total-to-total isentropic efficiency, allowing for th
     &e potential work '
      write(6,*) ' of all cooling flows and the pumping power          =        
     &               ',  etatt
      write(6,*)
      write(6,*) ' neglecting the pumping power, total-to-total isentrop
     &ic efficiency= ', etatt_nop
      write(6,*)
      write(6,*) '*****************************************************'
c
      if(ifcool.eq.1) then
      write(6,*) '*****************************************************'
      write(6,*)
      write(6,*) ' assuming no stagnation pressure loss or gain in the c
     &coolant supply passages'
      write(6,*) ' and using the static pressure at exit from the ejecti
     &on holes to estimate'
      write(6,*) ' the coolant supply stagnation  pressure'
      write(6,*) ' -then:- '
      write(6,*)
      write(6,*) ' total-to-total isentropic efficiency, allowing for th
     &e potential work '
      write(6,*) ' of all cooling flows and the pumping power          =        
     &               ',  etatt_at_hole
      write(6,*)
      write(6,*) ' neglecting the pumping power, total-to-total isentrop
     &ic efficiency= ', eta_ath_nop
      write(6,*)
      write(6,*) '*****************************************************'
      end if
c
 3333 continue
c
      end if
c
c****************************************************************************
c****************************************************************************
c
c     now calculate the machine performance if there are no cooling flows
c
      if(ifcool.eq.0) then
c
c******************************************************************************
c******************************************************************************
c
c      if no cooling flows work out the  isentropic efficiency for each stage
c      this is based on the mass average quantities at inlet and at jmm1
c      first for every stage
c
      nstgp1    = nstages + 1
      if(nstages.eq.1) nstgp1 = 1
c
c     nstgp1 is used to calculate the performance of the whole machine
c
c*****************************************************************************************
c
      do 2222  nstg = 1,nstgp1
c
      if(nstg.le.nstages) then
c
              j1      = jstg_start(nstg)
              j2      = jstg_end(nstg) + 1
              j2m1    = j2
c
c     jdd change  2/7/07  to include the mixing loss in the efficiency of the upstream blade row.
c     for the last stage there is no mixing downstream. base the effy on the last but one point.
c
              if(nstg.eq.nstages) j2m1 = jm -1 
c
c     end of jdd addition
c
c     jmix_stage is only used to find the stage reaction.
              n_row      = nrow(j1+1)
              jmix_stage = jmix(n_row) 
      else
c
c   for the whole machine there is no mixing downstream. base the effy on the last but one point, j2m1.
c
              j1   = 2
              j2   = jm
              j2m1 = jm - 1
              jmix_stage = (j1+j2)/2
      end if
c
      if(nstgp1.gt.1.and.nstg.eq.nstgp1) go to 123
c
      kmid = km/2
      phub_in = 0.0
      do  i = 1,imm1
      phub_in  = phub_in + 0.5*(p(i,j1,1)+p(i+1,j1,1))*fp(i)
      end do
c
      pmid_in = 0.0
      do  i = 1,imm1
      pmid_in  = pmid_in + 0.5*(p(i,j1,kmid)+p(i+1,j1,kmid))*fp(i)
      end do
c
c
      ptip_in = 0.0
      do  i = 1,imm1
      ptip_in  = ptip_in + 0.5*(p(i,j1,km)+p(i+1,j1,km))*fp(i)
      end do
c
c
      phub_mid = 0.0
      do  i = 1,imm1
      phub_mid  = phub_mid 
     &          + 0.5*(p(i,jmix_stage,1)+p(i+1,jmix_stage,1))*fp(i)
      end do
c
      pmid_mid = 0.0
      do  i = 1,imm1
      pmid_mid  = pmid_mid
     &        + 0.5*(p(i,jmix_stage,kmid)+p(i+1,jmix_stage,kmid))*fp(i)
      end do
c
      ptip_mid = 0.0
      do  i = 1,imm1
      ptip_mid  = ptip_mid 
     &          + 0.5*(p(i,jmix_stage,km)+p(i+1,jmix_stage,km))*fp(i)
      end do
c
c
      phub_out = 0.0
      do  i = 1,imm1
      phub_out  = phub_out 
     &          + 0.5*(p(i,j2m1,1)+p(i+1,j2m1,1))*fp(i)
      end do
c
      pmid_out = 0.0
      do  i = 1,imm1
      pmid_out  = pmid_out
     &          + 0.5*(p(i,j2m1,kmid)+p(i+1,j2m1,kmid))*fp(i)
      end do
c
      ptip_out = 0.0
      do  i = 1,imm1
      ptip_out  = ptip_out 
     &          + 0.5*(p(i,j2m1,km)+p(i+1,j2m1,km))*fp(i)
      end do
c
c     calculate the reactions as for a turbine
      reac_hub  = (phub_mid - phub_out)/(phub_in - phub_out)
      reac_mid  = (pmid_mid - pmid_out)/(pmid_in - pmid_out)
      reac_tip  = (ptip_mid - ptip_out)/(ptip_in - ptip_out)
c
c     change the reactions if a compressor
      if(pmid_out.gt.pmid_in) then
           reac_hub  = 1.0 - reac_hub
           reac_mid  = 1.0 - reac_mid
           reac_tip  = 1.0 - reac_tip
      end if
c
  123 continue
c 
      poin    = sumpo(j1)*po1(kmid)
      po2     = sumpo(j2m1)*po1(kmid)
      psin    = sumpstat(j1)
      psout   = sumpstat(j2m1)
      to_in   = sumto(j1)
      to_out  = sumto(j2m1)
      delto   = to_out - to_in
c
      if(ifgas.eq.0) then
           torat   = to_out/to_in
           porat   = po2/poin
           psrat   = psout/poin
           if(porat.lt.0.01)  porat   = 0.01
           to2is     = porat**fga *to_in
           tstat2_is = psrat**fga *to_in
           deltois = to2is  - to_in
           delt_tot_stat = tstat2_is - to_in
           etais   = 1.0
           if(delto.gt.0.0) etais = deltois/delto
           if(delto.lt.0.0) etais = delto/deltois
           epoly   = alog(torat)/alog(porat)/fga
           if(torat.gt.1.0) epoly = 1.0/epoly
c
           etats   = 1.0
           if(delto.gt.0.0) etats = delt_tot_stat/delto
           if(delto.lt.0.0) etats = delto/delt_tot_stat
c
      else
           porat   = po2/poin
           trat   = trat_from_prat(to_in,porat,fgagas,r_alpha,
     &              balpha1,balpha2)
           to_out_is = to_in*trat
           ho_out_is = hfromt(to_out_is,tref,href,cp1,cp2,cp3)
           ho_in     = hfromt(to_in,tref,href,cp1,cp2,cp3)
           ho_out    = hfromt(to_out,tref,href,cp1,cp2,cp3)
           delho     = ho_out    - ho_in
           delho_is  = ho_out_is - ho_in
           etais     = 1.0
           if(delho.gt.0.0) etais = delho_is/delho
           if(delho.lt.0.0) etais = delho/delho_is
           epoly     = 1.0
           delt  = to_in - tref
           cp = cp1 + cp2*delt + cp3*delt*delt
           ga = cp/(cp - rgas)
c
      end if
c
      write(6,*)'*******************************************************
     &*****************************************************************'
      if(nstg.le.nstages) then
           write(6,*) ' results for stage number ',nstg 
      else
           write(6,*) ' the following results are for the whole machine'
      end if
c
      write(6,*)'*******************************************************
     &*****************************************************************'
c
      write(6,*)
      if(delto.lt.0.0) then
           write(6,*) ' ********** this is a turbine ********** '
           write(6,*)
           write(6,*) ' pressure ratio, poinlet/poexit            = ',
     &                  1.0/porat 
      else
           write(6,*) ' ***********this is a compressor ********** '
           write(6,*)
           write(6,*) ' pressure ratio, poexit/poinlet            = ',
     &                  porat 
      end if
c
      write(6,*)      ' gas properties at stage inlet, cp = ',cp,
     &                ' gamma = ',ga,' r = ',rgas
      write(6,*)
      write(6,*)
      write(6,*) ' efficiencies based on values at j = ', j1,j2m1
      write(6,*) ' jmix  = ', jmix_stage
      write(6,*)
      write(6,*)      ' inlet and exit stagnation pressures   = ',
     &                  poin, po2
      write(6,*)      ' inlet and exit static pressures       = ',
     &                   psin,psout
      write(6,*)      ' inlet and exit stagnation temperatures= ',
     &                  to_in, to_out
      write(6,*)      ' inlet and outlet mass flow rates      = ',
     &                  blade_flow(j1),blade_flow(j2m1)
c
      write(6,*)
      write(6,*) ' using the mass averaged stagnation pressures and '
      write(6,*) ' temperatures at inlet and jm-1 .'
      write(6,*)
      write(6,*) 
     & ' total to total isentropic efficiency       = ', etais 
      if(ifgas.eq.0)  write(6,*) 
     & ' total to static isentropic efficiency      = ', etats 
      if(ifgas.eq.0)  write(6,*)
     & ' total to total polytropic efficiency       = ', epoly
      write(6,*)
c
      if(nstg.le.nstages) then
c
      write(6,*)      ' inlet, mid and exit static temperatures   = ',
     &                 sumtstat(j1),sumtstat(jmix_stage),sumtstat(j2m1)
c
                       reaction = (sumtstat(jmix_stage)- sumtstat(j2m1))
     &                 /(sumtstat(j1)         - sumtstat(j2m1))
                       if(psout.gt.psin) reaction = 1.0 - reaction
c     
      write(6,*)  
      write(6,*) 'average stage reaction based on mass averaged temperat
     &ures = ',   reaction
      write(6,*)
c
      write(6,*)
      write(6,*) 'reactions for stage number ', nstg
      write(6,*) 'hub reaction based on pressure changes     ',reac_hub
      write(6,*) 'mid span reaction based on pressure changes',reac_mid
      write(6,*) 'tip reaction based on pressure changes     ',reac_tip
      write(6,*)
c
      if(ifgas.eq.0) then
           work = (cp*(blade_flow(j1)*to_in
     &          - blade_flow(j2m1)*to_out))/1000.
      else
           work = (blade_flow(j1)*ho_in
     &          -  blade_flow(j2m1)*ho_out)/1000.
      end if
c
      write(6,*)
      write(6,*)  ' stage work neglecting any bleed or cooling flows = '
     &     , work,' kilowatts.'
      write(6,*)
c
c
      else
c
c********************************************************************************
c
c     now for the whole machine, nstg = nstages + 1
c
      if(ifbleed.eq.0) hbledtot = 0.0
c
      if(ifgas.eq.0) then
           work = (cp*(blade_flow(1)*to_in-blade_flow(jmm1)*to_out)
     &          +  enthin  - hbledtot  )/1000.
      else
           work = (blade_flow(1)*ho_in  - blade_flow(jmm1)*ho_out 
     &          + enthin  - hbledtot)/1000.
      end if
c   
      write(6,*)
      write(6,*) ' inlet enthalpy flux,  kw  = ',
     &             cp*blade_flow(1)*to_in/1000.
      write(6,*) ' outlet enthalpy flux, kw  = ',
     &             cp*blade_flow(jmm1)*to_out/1000.
      write(6,*) ' enthalpy bled off,    kw  = ', hbledtot/1000.
c
      write(6,*)
      write(6,*) ' for the whole machine '
      write(6,*) ' overall power output             = ',work,
     &                ' kilowatts.'
      write(6,*)
      write(6,*) ' inlet and outlet mass flow rates = ',
     &             blade_flow(1),blade_flow(jmm1),' kg/sec.'   
      write(6,*)
c        
      if(ifbleed.ne.0) then
           write(6,*)
           write(6,*) ' total mass flow bled off = ',sumbleed(jm),
     &                ' kg/sec'
           write(6,*)
           write(6,*) ' the power output allows for the bleed flows and 
     &cooling flows.'
           write(6,*) ' but the efficiencies only relate the mass averag
     & ed states of the fluid at machine entry and exit.'
           write(6,*)
      end if
c
      end if
c
 2222 continue
c
      end if
c
c
      return
      end
c
c********************************************************************************
c********************************************************************************
c
      subroutine smooth(n1,n2,nsmooth,fsmooth,frac,var)
c
      include 'commall-open-19.2'
c 
      dimension frac(jd),var(jd),temp(jd)
c
      do 10 its = 1,nsmooth
c
      do n = n1,n2
      temp(n) = var(n)
      end do
c
      do n = n1+1, n2-1
      fleft  = frac(n)   - frac(n-1)
      fright = frac(n+1) - frac(n)
      avg = ( fleft*temp(n+1) + fright*temp(n-1) )/(fleft + fright)
      var(n) = (1.0-fsmooth)*var(n) + fsmooth*avg
      end do
c
   10 continue
c
      return
      end
c******************************************************************************
c
      subroutine set_coeffs
c
c      this subroutine sets coefficients needed for the non-perfect gas property functions    
c
      include 'commall-open-19.2'
c
      href = cp1*tref
      cv1  = cp1 - rgas
      eref = cv1*tref
c
      write(6,*) 'tref, href,cv1, eref',tref, href,cv1, eref
c
      ht1 = 1.0/cp1   
      ht2 = -0.5*cp2/cp1**3
      ht3 = (3.0*cp2**2 - 2.0*cp3*cp1)/6/cp1**5
      ht4 = (-15*cp2**3/cp1**7 + 20*cp2*cp3/cp1**6  )/24
      write(6,*) 'ht1,ht2,ht3,ht4',ht1,ht2,ht3,ht4
c
      et1 = 1.0/cv1
      et2 = -0.5*cp2/cv1**3
      et3 = (3.0*cp2**2 - 2.0*cp3*cv1)/6/cv1**5
      et4 = (-15*cp2**3/cv1**7 + 20*cp2*cp3/cv1**6  )/24
      write(6,*) 'et1,et2,et3,et4',et1,et2,et3,et4
c
      alpha   = (cp1 - tref*cp2 + tref**2*cp3)/rgas
      r_alpha = 1.0/alpha
      beta1   = (cp2 -2.*cp3*tref)/rgas
      beta2   = 0.5*cp3/rgas
      balpha1 = -beta1/alpha
      balpha2 = -beta2/alpha 
c
      eps0    = 1.0/(alpha - 1.0)
      eps1    = -beta1*eps0
      eps2    = -beta2*eps0
      write(6,*) ' eps0,eps1,eps2 ',eps0,eps1,eps2
c
      gagas   = cp1/cv1
      fgagas  = (gagas-1.0)/gagas
      ga1gas  = gagas - 1.0
      write(6,*) ' gagas,fgagas,ga1gas ',gagas,fgagas,ga1gas
c
      return
      end
c
c******************************************************************************
c
      function hfromt(tin,tref,href,cp1,cp2,cp3)
c
      hfromt = cp1*tin + cp2*(tin-tref)*(tin-tref)/2 + 
     &         cp3*(tin-tref)*(tin-tref)*(tin-tref)/3
c
      return
      end  
c
c******************************************************************************
c
      function tfromh(hin,tref,href,ht1,ht2,ht3,ht4)
c
      dh     = hin - href
      dh2    = dh*dh
      tfromh = tref + ht1*dh + ht2*dh2 + ht3*dh*dh2 + ht4*dh2*dh2
c
      return
      end     
c
c******************************************************************************
c
      function tfrome(ein,tref,eref,et1,et2,et3,et4)
c
      de     = ein - eref
      de2    = de*de
      tfrome = tref + et1*de + et2*de2 + et3*de*de2 + et4*de2*de2
c
      return
      end  
c
c******************************************************************************
c
      function prat_from_trat(tstart,trat,alpha,beta1,beta2)
c
      t2 = tstart*trat
      beta  = beta1*(t2-tstart) + beta2*(t2*t2 - tstart*tstart)
      prat_from_trat  = trat**alpha   *   exp(beta)
c
      return
      end  
c
c******************************************************************************
c
c
      function trat_from_prat(tstart,prat,fgagas,r_alpha,
     &         balpha1,balpha2)
c
      palpha = prat**r_alpha
c
      trat   = prat**fgagas
      t2     = tstart*trat
      beta   = balpha1*(t2-tstart) + balpha2*(t2*t2-tstart*tstart)
      trat   = palpha  *  exp(beta)
      t2     = 0.5*(t2 + tstart*trat)
c
      beta   = balpha1*(t2-tstart) + balpha2*(t2*t2-tstart*tstart)
      trat   = palpha  *  exp(beta)
      t2     = 0.5*(t2 + tstart*trat)
c
      beta   = balpha1*(t2-tstart) + balpha2*(t2*t2-tstart*tstart)
      trat   = palpha  *  exp(beta)
      t2     = 0.5*(t2 + tstart*trat)
c
      beta   = balpha1*(t2-tstart) + balpha2*(t2*t2-tstart*tstart)
      trat   = palpha  *  exp(beta) 
      t2     = 0.5*(t2 + tstart*trat)
c
      trat_from_prat     =  t2/tstart
      return
      end
c
c******************************************************************************
c
c
      function trat_from_rorat(tstart,rorat,ga1gas,eps0,
     &         eps1,eps2)
c
      talpha = rorat**eps0
      trat   = rorat**ga1gas
      t2     = tstart*trat
      beta   = eps1*(t2-tstart) + eps2*(t2*t2-tstart*tstart)
      trat   = talpha  *  exp(beta)
      t2     = 0.5*(t2 + tstart*trat)
      beta   = eps1*(t2-tstart) + eps2*(t2*t2-tstart*tstart)
      trat   = talpha  *  exp(beta)
      t2     = 0.5*(t2 + tstart*trat)
      beta   = eps1*(t2-tstart) + eps2*(t2*t2-tstart*tstart)
      trat   = talpha  *  exp(beta)
      t2     = 0.5*(t2 + tstart*trat)
      beta   = eps1*(t2-tstart) + eps2*(t2*t2-tstart*tstart)
      trat   = talpha  *  exp(beta)
      t2     = 0.5*(t2 + tstart*trat)
      trat_from_rorat     =  t2/tstart
      return
      end
c
c**********************************************************************
c
c      
      subroutine mass_avg(im,jm,km,fcn,flowx,sumfun,avgall,id,jd,kd)
c
      dimension fcn(id,jd,kd),flowx(id,jd,kd), sumfun(kd)
c
      j = jm
      sumall   = 0.0
      flowall  = 0.0
c
      do 110 k=1,km-1
      sumas    = 0.0
      sumfun(k)= 0.0
      do 120 i=1,im-1
      dflow = -flowx(i,j,k)
      if(dflow.lt.0.0) dflow = 0.0
      sumas     = sumas+dflow
      sumfun(k) = sumfun(k) + dflow*(fcn(i,j,k)+fcn(i+1,j,k)
     &          + fcn(i+1,j,k+1)+fcn(i,j,k+1))*0.25  
  120 continue
      sumall    = sumall  + sumfun(k)
      flowall   = flowall + sumas
      sumfun(k) = sumfun(k)/sumas
  110 continue
c
      do 200 k=2,km-1
      sumfun(k) = 0.5*(sumfun(k) + sumfun(k-1))
  200 continue
c
      sumfun(1)  = 2.0*sumfun(2)    - sumfun(3)
      sumfun(km) = 2.0*sumfun(km-1) - sumfun(km-2)
      avgall     = sumall/flowall
c
      return
      end
c
c************************************************************************************
c************************************************************************************
c
      subroutine restagg(nsec,j1,j2,jlerow,jterow,rotate,fracx_rot)
c
c************************************************************************************
c     this subroutine restaggers a blade section.
c     it assumes negligible change of radius of the stream surfaceand so is not usable
c     for radial flow machines.
c
c     it is not used in  version 16.3  and above. use subroutine "restagger" instead .
c************************************************************************************
      include 'commall-open-19.2'
c
      dimension  xin(jd),yup(jd),ylow(jd),ymid(jd),
     &           xup(jd),xmid(jd),xlow(jd)
c
      rotate     = 0.0
      fracx_rot  = 0.5
      read(5,*) dummy_input
      read(5,*,err=1502) rotate, fracx_rot
 1502 continue
      write(6,*) ' input for restaggering option'
      write(6,*) ' rotation =', rotate, 'about ',fracx_rot 
      write(6,*)
c
      rad_rot = rotate*degrad
c
c     find the y coordinates of the upper and lower blade surfaces.
c
      nj  = j2 - j1 + 1
      do 10 j  = j1,j2
      jloc  = j - j1 + 1
      xin(jloc)   = xsurf(j,nsec)
      yup(jloc)   = rt_upp(j,nsec)
      ylow(jloc)  = yup(jloc) - rt_thick(j,nsec)
      ymid(jloc)  = 0.5*(yup(jloc) + ylow(jloc))      
   10 continue
c
c     find the coordinates of the axis of rotation
c
      xrot = xin(jlerow) + fracx_rot*(xin(jterow) - xin(jlerow))
      call intp(nj,xin,ymid,xrot,yrot)
c
c     rotate all the points clockwise by rad_rot radians.
c
      do 20 j = 1,nj
c
      xrel      = xin(j) - xrot
      yuprel    = yup(j) - yrot
      ymidrel   = ymid(j)- yrot
      ylowrel   = ylow(j)- yrot
      anglup    = atan2(yuprel,xrel) - rad_rot
      anglmid   = atan2(ymidrel,xrel)- rad_rot
      angllow   = atan2(ylowrel,xrel)- rad_rot
      radup     = sqrt(xrel*xrel + yuprel*yuprel)
      radmid    = sqrt(xrel*xrel + ymidrel*ymidrel)
      radlow    = sqrt(xrel*xrel + ylowrel*ylowrel)
      xup(j)    = radup*cos(anglup)   + xrot
      xmid(j)   = radmid*cos(anglmid) + xrot
      xlow(j)   = radlow*cos(angllow) + xrot
      yup(j)    = radup*sin(anglup)   + yrot
      ymid(j)   = radmid*sin(anglmid) + yrot
      ylow(j)   = radlow*sin(angllow) + yrot
c      
c      write(6,100) j,xmid(j),yup(j),ymid(j),ylow(j)
   20 continue

c
c    interpolate to find new  y values on the blade surfaces at the xmid values of x.
c    scale the movement by facmove so that it is reduced upstream and downstream of the
c    blade such that first and last grid points do not move unless they are the le and te points.
c
      write(6,*)
      write(6,*) ' blade restaggered, new values of x,yup,ythick,fmove=' 
      do 30 j = 1,nj
      jall = j1 + j - 1
      facmove = 1.0
      if(j.lt.jlerow) facmove=(xmid(j)-xmid(1))/(xmid(jlerow)-xmid(1))
      if(j.gt.jterow) facmove=(xmid(nj)-xmid(j))/(xmid(nj)-xmid(jterow)) 
      fm1 = 1.0 - facmove
      xarg = xmid(j)
      call intp(nj,xup,yup,xarg,yupnew)
      call intp(nj,xlow,ylow,xarg,ylownew)
      xsurf(jall,nsec)    = fm1*xsurf(jall,nsec)  + facmove*xarg
      rt_upp(jall,nsec)   = yupnew
      rt_thick(jall,nsec) = (yupnew - ylownew)
      write(6,100) jall,xsurf(jall,nsec),rt_upp(jall,nsec),
     &             rt_thick(jall,nsec),facmove
   30 continue
c
  100 format('j = ',i5,4f15.7)
c
      return
      end
c
c******************************************************************************
c******************************************************************************
c
      subroutine lean(k,j1,j2,anglean)
c
      include 'commall-open-19.2'
c******************************************************************************
c******************************************************************************
c     lean the blade by anglean if if_lean is greater than zero.
c     if "anglean" is positive the hub is held fixed and the other sections are 
c     moved in the positive theta direction .
c
      read(5,*) dummy_input
c
      anglean = 0.0
      read(5,*,err=20) anglean
   20 continue
c
      write(6,*) ' leaning the blade section by anglean = ',
     &             anglean,' degrees relative to the hub section.'
c
      anglean = anglean*degrad
c
      do 10 j= j1,j2
      jm1 = j-1
      jp1 = j+1
      if(j.eq.j1) jm1 = 1
      if(j.eq.j2) jp1 = j2
      xdif = xsurf(jp1,1) - xsurf(jm1,1)
      rdif = rsurf(jm1,1) - rsurf(jp1,1)
      sdif = sqrt(xdif*xdif + rdif*rdif)
      rnorm   =   xdif/sdif
      xnorm   =  -rdif/sdif
      qdistr  =   rsurf(j,k) - rsurf(j,1)
      qdistx  =   xsurf(j,k) - xsurf(j,1)
      qdistn  =   qdistr*rnorm + qdistx*xnorm         
      t_shift =   anglean*qdistn
      rt_upp(j,k)   = t_shift + rt_upp(j,k)
   10 continue
c
      return
      end
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c
      subroutine mix_bconds(ifout)
c
c***********************************************************************
c     write the mass average values at exit to a file 'outbconds'
c***********************************************************************
c
      include 'commall-open-19.2'
c
      dimension favg(kd)
c
      open(unit = 12,file ='mixbconds')
c
c     calculate and write out the exit stagnation pressure.
c
      do 7777 nr = 1, nrows + 1
c
      write(12,*)
      write(12,*) ' blade row number ', nr, ' jmix = ',jmix(nr) 
      write(12,*)
c
      if(nr.eq.1) then 
c
      j = 1
      if(ifout.eq.1) then
               write(12,*)
               write(12,*) ' inlet boundary conditions to row number 1.'
               write(12,*)
      end if
c
      else
c
      j = jmix(nr-1)
      if(ifout.eq.1) then
      write(12,*)  ' row number ', nr-1,' jmix  =  ', j
      write(12,*) ' exit flow conditions from this blade row'
      write(12,*) ' these can be used as the inlet boundary conditions 
     &for the next stage.'
      write(12,*)
      end if
c
      end if
c
c
      do 1110 i=1,im
      do 1110 k=1,km
      d=.5*(vx(i,j,k)*vx(i,j,k)+vt(i,j,k)*vt(i,j,k)+vr(i,j,k)*vr(i,j,k))
      tstag     = ho(i,j,k)/cp
      ts        = tstag - d/cp
      if(ts.lt.1.) ts = 1.
      store(i,j,k) = p(i,j,k)*(tstag/ts)**(ga/(ga-1))
 1110 continue
c
      call mass_avg(im,j,km,store,flowx,po_out_avg,poavgjm,id,jd,kd)
      if(ifout.eq.1) then
             write(12,*) ' stagnation pressure '
             write(12,1111) (po_out_avg(k),k=1,km)
      end if
 1111 format(8f10.1)
c
c     calculate and write out the exit static pressure
c
      call mass_avg(im,j,km,p,flowx,favg,psavgjm,id,jd,kd)
      if(ifout.eq.1) then
            write(12,*) ' static pressure '
            write(12,1111)(favg(k),k=1,km)
      end if
c
c
c     calculate and write out the exit stagnation temp.
c
      do 1014 i=1,im
      do 1014 k=1,km
      store(i,j,k) =  ho(i,j,k)/cp
 1014 continue
c
      call mass_avg(im,j,km,store,flowx,to_out_avg,toavgjm,id,jd,kd)
      if(ifout.eq.1) then
            write(12,*) ' stagnation temperature '
            write(12,1111) (to_out_avg(k),k=1,km)
      end if
 1112 format(8f10.4)
c
c
c     calculate and write out the exit static temperature.
c
      do 1119 i = 1, im
      do 1119 k = 1, km
      d=.5*(vx(i,j,k)*vx(i,j,k)+vt(i,j,k)*vt(i,j,k)+vr(i,j,k)*vr(i,j,k))
      store(i,j,k) = (ho(i,j,k) - d)/cp
 1119 continue
c
      call mass_avg(im,j,km,store,flowx,favg,tsavgjm,id,jd,kd)
      if(ifout.eq.1) then
            write(12,*) ' static temperature '
            write(12,1111) (favg(k),k=1,km)
      end if
c
c     write out the exit tangential velocity.
c
      call mass_avg(im,j,km,vt,flowx,vt_out_avg,vtavgjm,id,jd,kd)
      if(ifout.eq.1) then
            write(12,*) ' absolute tangential velocity '
            write(12,1112) (vt_out_avg(k),k=1,km)
      end if
c
c     calculate and write out the exit meridional velocity
c
      do 1120 i=1,im
      do 1120 k=1,km
      store(i,j,k) = sqrt(vx(i,j,k)*vx(i,j,k) + vr(i,j,k)*vr(i,j,k))
 1120 continue
c
      call mass_avg(im,j,km,store,flowx,favg,vmavgjm,id,jd,kd)
      if(ifout.eq.1) then
            write(12,*) ' meridional velocity '
            write(12,1112) (favg(k),k=1,km)
      end if
c
c    calculate and write out the yaw angle at exit
c
      do 1130 i=1,im
      do 1130 k=1,km
      vmer = sqrt(vx(i,j,k)*vx(i,j,k) + vr(i,j,k)*vr(i,j,k))
      store(i,j,k) = atan2(vt(i,j,k),vmer)*raddeg
 1130 continue
c
      call mass_avg(im,j,km,store,flowx,yaw_out_avg,yawavgjm,id,jd,kd)
      if(ifout.eq.1) then
            write(12,*) ' yaw angle '
            write(12,1112) (yaw_out_avg(k),k=1,km)
      end if 
c
c    calculate and write out the pitch angle at exit
c
      do 1140 i=1,im
      do 1140 k=1,km
      vmer = sqrt(vx(i,j,k)*vx(i,j,k) + vr(i,j,k)*vr(i,j,k))
      store(i,j,k) = atan2(vr(i,j,k),vmer)*raddeg
 1140 continue
c
      call mass_avg(im,j,km,store,flowx,pitch_out_avg,
     &              pitchavgjm,id,jd,kd)
      if(ifout.eq.1) then
            write(12,*) ' pitch angle '
            write(12,1112) (pitch_out_avg(k),k=1,km)
      end if
c
c      write out   fr  and   fp  .
c
      if(ifout.eq.1) then
            write(12,*) ' spanwise grid spacing '
            write(12,1113) (fr(k),k=1,km-1)
            write(12,*) ' pitchwise grid spacing '
            write(12,1113) (fp(i),i=1,im-1)
            fspan(1) = 0.0
            do k=2,km
            fspan(k) = fspan(k-1) + fr(k-1)
            end do
            write(12,*) ' spanwise grid positions'
            write(12,1113) (fspan(k),k=1,km)            
      end if
c
 1113 format(8f10.6)
c
c     end of one blade row, return for next row
 7777 continue
c
c      end of writing file to unit 12
c
      close(12)
c  
c*********************************************************************
c*********************************************************************
      return
      end

c
c******************************************************************************
c
      subroutine newbconds 
c
      include 'commall-open-19.2'
c
c     this subroutine changes the inlet boundary conditions to simulate a repeating stage
c
      rfinbc1  =  1.0 - rfinbc
c
      do 10 k = 1,km
      po1(k)  = rfinbc1*po1(k) 
     &        + rfinbc*(po_in_mid  + po_out_avg(k)   - po_out_avg(kmid))
c      to1(k)  = rfinbc1*to1(k)
c     &        + rfinbc*(to_in_mid  + to_out_avg(k)   - to_out_avg(kmid))
c      bs(k)   = rfinbc1*bs(k)
c     &        + rfinbc*(yaw_in_mid + yaw_out_avg(k) - yaw_out_avg(kmid))
c
       bs(k)   = rfinbc1*bs(k) + rfinbc*yaw_out_avg(k)
c
c      br(k)   = rfinbc1*br(k)
c     &      + rfinbc*(pitch_in_mid+pitch_out_avg(k)-pitch_out_avg(kmid))
      vtin(k) = rfinbc1*vtin(k)
     &        + rfinbc*(vt_in_mid  + vt_out_avg(k)   - vt_out_avg(kmid))
      rpo1(k) = 1./po1(k)
      bssin(k)= sin(bs(k)*degrad)
      bscos(k)= cos(bs(k)*degrad)
      brsin(k)= sin(br(k)*degrad)
      brcos(k)= cos(br(k)*degrad)
   10 continue
c
      return
      end
c*****************************************************************************      
c
      subroutine cell_to_node(vcell,vnode)
c
      include 'commall-open-19.2'
c
      dimension  vcell(id,jd,kd), vnode(id,jd,kd) 
c    
c     distribute the cell centred variable vcell 
c     to the node stored variable vnode
c*****************************************************************************
c
c     first average along the i direction
c
      if(im.gt.2) then

      do 1100 k = 1,kmm1
      do 1100 j = 1,jmm1
      do 1101 i = 2,imm1
      temp1(i,j,k)  = (vcell(i-1,j,k) + vcell(i,j,k))
 1101 continue
      temp1(1,j,k)  = 3.0*vcell(1,j,k)    - vcell(2,j,k)
      temp1(im,j,k) = 3.0*vcell(imm1,j,k) - vcell(imm2,j,k)
 1100 continue

      else

      do 1102 k=1,kmm1
      do 1102 j=1,jmm1
      temp1(1,j,k) = vcell(1,j,k)
      temp1(2,j,k) = vcell(1,j,k)
 1102 continue

      end if
c
c     next average along the k direction
c  q3d
c
      if(km.gt.2) then
c
      do 1120 j = 1,jmm1
      do 1120 i = 1,im
      do 1121 k = 2,kmm1
      temp2(i,j,k) = (temp1(i,j,k-1) + temp1(i,j,k))
 1121 continue
      temp2(i,j,1)   = 3.0*temp1(i,j,1)    - temp1(i,j,2)
      temp2(i,j,km)  = 3.0*temp1(i,j,kmm1) - temp1(i,j,kmm2)
 1120 continue
c
      else
c
      do 1122 j = 1,jmm1
      do 1122 i = 1,im
      do 1122 k=1,2
      temp2(i,j,k) = temp1(i,j,1) + temp1(i,j,2)
 1122 continue
c
      end if
c
c   end q3d
c
c      now average along the j direction
c
      do 1125 k  = 1,km
      do 1125 i  = 1,im
      do 1126 j  = 2,jmm1
      vnode(i,j,k)   = 0.125*(temp2(i,j-1,k) + temp2(i,j,k))
 1126 continue
      vnode(i,1,k)   = 0.375*temp2(i,1,k)    - 0.125*temp2(i,2,k)
      vnode(i,jm,k)  = 0.375*temp2(i,jmm1,k) - 0.125*temp2(i,jmm2,k)
 1125 continue
c
c
      return
      end
c
c**************************************************************************************
c
      subroutine set_xlength
c
c********************************************************************************
c********************************************************************************
c    this subroutine calculates the distance from every grid point to the nearest
c    solid wall and uses this to set the mixing length.
c********************************************************************************
c********************************************************************************
      include 'commall-open-19.2'
c
c   the nearest wall is sought over the range  j  +/- jrange,  k  +/- krange
c
      write(6,*) ' starting set_xlength. '
c
c     jrange and krange are the range of points over which the nearest 
c     wall point is sought
      jrange = 20
      krange = 5
c
c
      do 1000 k=1,km
      do 1000 j=1,jm
      nr       = nrow(j)
      pitch    = 6.283185*r(j,k)/nblade(j)
      dist_ref = 0.5*pitch
c
c********************************************************************************
c********************************************************************************
      do 1000 i=1,im

c     find the distance to the hub or casing
c
      iwall = imid
      j1 = j - jrange
      if(j1.lt.2)  j1 = 2
      j2 = j + jrange
      if(j2.gt.jm) j2 = jm
c
      hubdist = 1.0e10
      tipdist = 1.0e10
      dmin    = 1.0e10
c
c   find the nearest point on the hub
      do 10 jwall = j1,j2
      xdif  = x(j,k) - x(jwall,1)
      rdif  = r(j,k) - r(jwall,1)
      distsq = xdif*xdif + rdif*rdif
      if(distsq.lt.dmin) then
           dmin = distsq
           jmin = jwall
      end if
   10 continue
c
c    find the perpendicular distance to the nearest point on the hub
      atot  = sqrt(asx(iwall,jmin,1)*asx(iwall,jmin,1)
     &      + asr(iwall,jmin,1)*asr(iwall,jmin,1))
      xnorm = asx(iwall,jmin,1)/atot
      rnorm = asr(iwall,jmin,1)/atot
      xdif  = x(j,k) - x(jmin,1)
      rdif  = r(j,k) - r(jmin,1)
      hubdist = abs(xdif*xnorm + rdif*rnorm)
      if(ibound.eq.1.or.ibound.gt.2) hubdist = dist_ref
c
c   find the nearest point on the casing
      dmin    = 1.0e10
      do 15 jwall = j1,j2
      xdif  = x(j,k) - x(jwall,km)
      rdif  = r(j,k) - r(jwall,km)
      distsq = xdif*xdif + rdif*rdif
      if(distsq.lt.dmin) then
           dmin = distsq
           jmin = jwall
      end if
   15 continue
c
c    find the perpendicular distance to the nearest point on the casing
      atot  = sqrt(asx(iwall,jmin,km)*asx(iwall,jmin,km)
     &      + asr(iwall,jmin,km)*asr(iwall,jmin,km))
      xnorm = asx(iwall,jmin,km)/atot
      rnorm = asr(iwall,jmin,km)/atot
      xdif  = x(j,k) - x(jmin,km)
      rdif  = r(j,k) - r(jmin,km)
      tipdist = abs(xdif*xnorm + rdif*rnorm)
      if(ibound.ge.2) tipdist = dist_ref
c
c   store the distance to the nearest end wall
c
       endwall_dist = hubdist*tipdist/(hubdist + tipdist)
c
c*************************************************************************
c*************************************************************************
c     next find the distance to the nearest blade surface
c
      j1 = j - jrange
      if(j1.lt.2) j1 = 2
      if(j1.lt.jstart(nr)) j1 = jstart(nr)
      if(j1.ge.jte(nr))    j1 = jte(nr) - 1
      j2 = j + jrange
      if(j2.gt.jm) j2 = jm
      if(j2.gt.jmix(nr)) j2 = jmix(nr)
      if(j2.le.jle(nr))  j2 = jle(nr) + 1
      k1 = k - krange
      if(k1.lt.1) k1 = 1
      k2 = k + krange
      if(k2.gt.km) k2 = km
c
c    find the nearest point on the  i = 1 blade surface
c
      dmin = 1.0e10
      if_found = 0
      do 20 jsurf = j1,j2
      do 25 ksurf = k1,k2  
      xdif  = x(j,k)  - x(jsurf,ksurf)
      rdif  = r(j,k)  - r(jsurf,ksurf)
c
c      tdif  = rtheta(i,j,k) - rtheta(1,jsurf,ksurf)
c      changed by ws/lx  17 may, 2017  
      tdifsq = 4.* r(j,k)*r(jsurf,ksurf)*sin((theta(i,j,k)
     &       - theta(1,jsurf,ksurf))/2.)**2
      distsq = xdif*xdif + rdif*rdif + tdifsq
c
      if(distsq.lt.dmin) then
           dmin = distsq
           jmin = jsurf
           kmin = ksurf
           if_found = 1
      end if
   25 continue
   20 continue
c
      if(if_found.eq.1) then
c     find the perpendicular distance to the nearest point on the i=1 blade surface.
      knorm = kmin
      if(kmin.eq.km) knorm = km-1
      atot  = sqrt(abx(1,jmin,knorm)*abx(1,jmin,knorm)  
     &      + abr(1,jmin,knorm)*abr(1,jmin,knorm)
     &      + abt(jmin,knorm)*abt(jmin,knorm))
      xnorm  = abx(1,jmin,knorm)/atot
      rnorm  = abr(1,jmin,knorm)/atot
      tnorm  = abt(jmin,knorm)/atot
      xdif   = x(j,k)  - x(jmin,kmin)
      rdif   = r(j,k)  - r(jmin,kmin)
c
c      tdif   = rtheta(i,j,k) - rtheta(1,jmin,kmin)
c      changed by ws/lx  17 may, 2017  
      tdif = 2.* sqrt(r(j,k)*r(jmin,kmin))*sin((theta(i,j,k)
     &       - theta(1,jmin,kmin))/2.) 
c
      ssdist = abs(xdif*xnorm + rdif*rnorm + tdif*tnorm)
c   jdd added 11/13 to improve treatment of thin leading and trailing edges.
      if(jmin.le.jle(nr).or.jmin.ge.jte(nr))
     &       ssdist = sqrt(xdif*xdif + rdif*rdif+ tdif*tdif)
      else
             ssdist = dist_ref
      end if
c
c    find the nearest point on the  i = im blade surface
c
      dmin = 1.0e10
      if_found = 0
      do 30 jsurf = j1,j2
      do 35 ksurf = k1,k2  
      xdif  = x(j,k)  - x(jsurf,ksurf)
      rdif  = r(j,k)  - r(jsurf,ksurf)
c
c      tdif  = rtheta(i,j,k) - rtheta(im,jsurf,ksurf)  
c      changed by ws/lx  17 may, 2017  
      tdifsq = 4.* r(j,k)*r(jsurf,ksurf)*sin((theta(i,j,k)
     &       - theta(im,jsurf,ksurf))/2.)**2
      distsq = xdif*xdif + rdif*rdif + tdifsq
c
      if(distsq.lt.dmin) then
           dmin = distsq
           jmin = jsurf
           kmin = ksurf
           if_found = 1
      end if
   35 continue
   30 continue
c
c     find the perpendicular distance to the nearest point on the i=im blade surface.
      knorm = kmin
      if(if_found.eq.1) then
      if(kmin.eq.km) knorm = km-1
      atot = sqrt(abx(im,jmin,knorm)*abx(im,jmin,knorm)  
     &     + abr(im,jmin,knorm)*abr(im,jmin,knorm)
     &     + abt(jmin,knorm)*abt(jmin,knorm))
      xnorm = abx(im,jmin,knorm)/atot
      rnorm = abr(im,jmin,knorm)/atot
      tnorm = abt(jmin,knorm)/atot
      xdif  = x(j,k)  - x(jmin,kmin)
      rdif  = r(j,k)  - r(jmin,kmin)
c
c      tdif  = rtheta(i,j,k) - rtheta(im,jmin,kmin)
c      changed by ws/lx  17 may, 2017  
      tdif = 2.* sqrt(r(j,k)*r(jmin,kmin))*sin((theta(i,j,k)
     &       - theta(im,jmin,kmin))/2.)
c
      psdist = abs(xdif*xnorm + rdif*rnorm + tdif*tnorm)
c   jdd added 11/13 to improve treatment of thin leading and trailing edges.
      if(jmin.le.jle(nr).or.jmin.ge.jte(nr))
     &      psdist = sqrt(xdif*xdif + rdif*rdif+ tdif*tdif)
      else
            psdist = dist_ref
      end if
c
c     store the distance to the nearest blade surface.
c
      blade_dist = ssdist*psdist/(ssdist + psdist)
c
c*******************************************************************************  
c*******************************************************************************    
c   set the distance from the nearest solid surface as  xdist*ydist/sqrt(xdist**2 + ydist**2)
c
      if(k.gt.1.and.k.lt.km) then
      distsqrt  = sqrt(endwall_dist*endwall_dist
     &          + blade_dist*blade_dist)
      if(distsqrt.lt.1.0e-10) distsqrt = 1.0e-10
      dist_min(i,j,k) = blade_dist*endwall_dist/distsqrt
      else
      dist_min(i,j,k) = 0.0
      end if
c
c  q3d
      if(km.eq.2) dist_min(i,j,k) = blade_dist
      if(im.eq.2) dist_min(i,j,k) = endwall_dist
c  end q3d
c*******************************************************************************  
c*******************************************************************************
c   set the mixing length limit, xlimit(j) varying with meridional distance.
c
      if(i.eq.imid.and.k.eq.kmid) then
            if(j.le.jle(nr))  xllim = xllim_in(nr) +
     &      (smerid(j,k) - smerid(jstart(nr),k))
     &     /(smerid(jle(nr),k) - smerid(jstart(nr),k))
     &     *(xllim_le(nr) - xllim_in(nr))
            if(j.gt.jle(nr).and.j.le.jte(nr)) xllim = xllim_le(nr) +
     &      (smerid(j,k) - smerid(jle(nr),k))
     &     /(smerid(jte(nr),k) - smerid(jle(nr),k))
     &     *(xllim_te(nr) - xllim_le(nr))
            if(j.gt.jte(nr))  xllim = xllim_te(nr) +
     &      (smerid(j,k) - smerid(jte(nr),k))
     &     /(smerid(jmix(nr),k) - smerid(jte(nr),k))
     &     *(xllim_dn(nr) - xllim_te(nr))
c
c     use a factor of 2 on xlimit so that is roughly the throat width.
c
            xlimit(j) = 2.0*dist_min(i,j,k)*xllim
c
      end if 
c   end of loop over all  i,j, k,  points .
 1000 continue
c
c*******************************************************************************  
c*******************************************************************************  
c*******************************************************************************  
c     blend the linear region of the mixing length and the constant region at 
c     ratio = dblend.
c     check if the mixing length is greater than the mixing length limit
c     then store the mixing length for each point as store(i,j,k)
c
      dblend = 0.6666666
      expon  = dblend/(1.0 - dblend)
      fblend = (1.0 - dblend)*dblend**expon
c     
      do 2100 j = 1,jm
      do 2000 k = 1,km
      do 2000 i = 1,im
      x_length  = dist_min(i,j,k)
      ratio   = x_length/xlimit(j)
      if(ratio.gt.dblend)
     &           x_length = xlimit(j)*(1.0 - fblend*ratio**(-expon))         
      store(i,j,k) = x_length
c     temp3  is used in the next section to damp the free stream turbulence near a wall.
      temp3(i,j,k) = (x_length/xlimit(j))**4
 2000 continue
 2100 continue
c
c*******************************************************************************      
c    now average the mixing length for a cell and store the average 
c    multiplied by the von karmen constant  (= 0.41) all squared.
c
      do 2200 k = 1,kmm1
      do 2200 j = 2,jm
      nr = nrow(j)
      do 2200 i = 1,imm1
      avg_xlength = 0.125*(store(i,j,k) + store(i,j,k+1)
     &     + store(i,j-1,k) + store(i,j-1,k+1) + store(i+1,j,k)
     &     + store(i+1,j,k+1) + store(i+1,j-1,k) + store(i+1,j-1,k+1))
      xlength(i,j,k) = 0.41*0.41*avg_xlength*avg_xlength
c
c     set the mixing length to zero for the first cell on a blade row.
      if(j.eq.jstart(nr)) xlength(i,j,k) = 0.0 
c   set the function used to damp the free stream turbulence near a wall.
      avg_fstrat = 0.125*(temp3(i,j,k) + temp3(i,j,k+1)
     &     + temp3(i,j-1,k) + temp3(i,j-1,k+1) + temp3(i+1,j,k)
     &     + temp3(i+1,j,k+1) + temp3(i+1,j-1,k) + temp3(i+1,j-1,k+1))
      fst_rat(i,j,k) = avg_fstrat*fsturb(nr)
c     set the free stream turbulence to zero for the first cell on a blade row.
      if(j.eq.jstart(nr)) fst_rat(i,j,k) = 0.0 
c
 2200 continue
c
c*******************************************************************************  
c     set the average wall distance for a cell , squared .
c
      do 2300 k=1,kmm1
      do 2300 j=2,jm
      do 2300 i=1,imm1
      avgdist = 0.125*(dist_min(i,j,k) + dist_min(i,j,k+1)
     &   + dist_min(i,j-1,k)   + dist_min(i,j-1,k+1) + dist_min(i+1,j,k)
     &   + dist_min(i+1,j,k+1) + dist_min(i+1,j-1,k)
     &   + dist_min(i+1,j-1,k+1))
      dwallsq(i,j,k) = avgdist*avgdist
 2300 continue
c
c*******************************************************************************  
c
      write(6,*)
      write(6,*) ' leaving subroutine    set_xlength. '
      write(6,*)
c
      return
      end
c******************************************************************************
c*************************************************************************************
c******************************************************************************
c
      subroutine new_loss
c
c******************************************************************************
c************this subroutine computes a body force based on wall functions for the
c            surface shear stress and a mixing length model of eddy viscosity.
c
c            new model based on tblock model. january 2010.
c
c
      include  'commall-open-19.2'
c
      common/bkstress/  txx(id,jd,kd),txr(id,jd,kd),
     &                  txt(id,jd,kd),trx(id,jd,kd),trr(id,jd,kd),
     &                  trt(id,jd,kd),ttx(id,jd,kd),ttr(id,jd,kd),
     &                  ttt(id,jd,kd),qxx(id,jd,kd),qrr(id,jd,kd),
     &                  qtt(id,jd,kd)
c
      dimension  forcex(id,jd,kd),forcer(id,jd,kd),forcet(id,jd,kd),
     &           esource(id,jd,kd),tempp(jd)
c
c      calculate the viscosity over the first quarter of the steps.
c      then hold it constant for the remainder of the steps.
c
      if(nstep.eq.1.or.nstep.lt.nsteps_max/4) then
c
      if(reyno.gt.100.) then
           j1=jle(1)
           j2=jte(1)
           if(jle(1).gt.jm) j1=1
           if(jte(1).gt.jm) j2=jm
           xchord  = smerid(j2,kmid)-smerid(j1,kmid)
           row2    = sqrt(rovx(imid,j2,kmid)*rovx(imid,j2,kmid)
     &             + rovr(imid,j2,kmid)*rovr(imid,j2,kmid) 
     &             + rowt(imid,j2,kmid)*rowt(imid,j2,kmid) )
           vislam  = xchord*row2/reyno
      end if
c
      if(reyno.gt.0.0.and.reyno.lt.99.99) then
            vislam = reyno/100000.
      end if
c
      if(reyno.lt.0.0) vislam = -reyno*1.0e-5
c
      thcond  = cp*vislam/prandtl
      ftcond = thcond/vislam
c
c   end of part only called over first 1/4 of the steps.
      endif
c
c********************************************************************************
c     set the viscous forces and energy sources to zero
c
      do 100 k=1,kmm1
      do 100 j=2,jm
      do 100 i=1,imm1
      forcex(i,j,k)  = 0.0
      forcer(i,j,k)  = 0.0
      forcet(i,j,k)  = 0.0
      esource(i,j,k) = 0.0
  100 continue
c
c*******************************************************************************************
c*******************************************************************************************
c     start to evaluate the turbulent viscosity for each cell.
c
      do 150 k=1,kmm1
      kp1 = k+1
      do 150 j=2,jm
      jm1 = j-1
      nrw = nrow(j)
      jtedge = jte(nrw)
      jledge = jle(nrw)
      do 150 i=1,imm1
      ip1 = i+1
c
c    average the conditions for the cells. note these are true averages.
c
      wtavg = 0.125*(wt(i,j,k)+wt(ip1,j,k)+wt(ip1,j,kp1)+wt(i,j,kp1)
     &      + wt(i,jm1,k)+wt(ip1,jm1,k)+wt(ip1,jm1,kp1)+wt(i,jm1,kp1))
      roavg = roavg_cell(i,j,k)
      tsavg = 0.125*(t_static(i,j,k)+t_static(ip1,j,k)
     &        +t_static(ip1,j,kp1)+t_static(i,j,kp1)+t_static(i,jm1,k)
     &        +t_static(ip1,jm1,k)+t_static(ip1,jm1,kp1)
     &        +t_static(i,jm1,kp1))
      ravg  = ravg_cell(j,k)
c
c     calculate the derivatives of the velocity components and temperature.
c     note the vorticity is based on the relative velocity,  wt.
c     this was changed from the absolute velocity by jdd august 2017.
      call gradvel(i,j,k,vx,dvxdx,dvxdr,dvxdt)
      call gradvel(i,j,k,vr,dvrdx,dvrdr,dvrdt)
      call gradvel(i,j,k,wt,dwtdx,dwtdr,dwtdt)
      call gradvel(i,j,k,t_static,dtempdx,dtempdr,dtempdt)
c
c     calculate the rates of strain for calculating the viscous stresses.
c     these are not yet the stresses as the final viscosity is not yet known.
c     
      txr(i,j,k) = dvxdr + dvrdx
      txx(i,j,k) = dvxdx
      txt(i,j,k) = dwtdx + dvxdt
      trr(i,j,k) = dvrdr
      trt(i,j,k) = dwtdr + dvrdt - wtavg/ravg
      ttt(i,j,k) = dwtdt
      qxx(i,j,k) = dtempdx
      qrr(i,j,k) = dtempdr
      qtt(i,j,k) = dtempdt
c
c        use the vorticity to determine the turbulent viscosity for each cell.
c
      vortx   = dvrdt - dwtdr - wtavg/ravg 
      vortr   = dwtdx - dvxdt
      vortt   = dvxdr - dvrdx
      vort    = sqrt(vortx*vortx + vortr*vortr + vortt*vortt)
c  set a limit to the vorticity
      if(vort.gt.vort_max) vort = vort_max
c
      abs_vort(i,j,k) = vort
c
c   calculate the turbulent viscosity in the boundary layers.
c
      turbvis_wall = roavg*xlength(i,j,k)*vort*fmixup
c.
c    use fyplus to reduce the turbulent viscosity in the buffer region. 
c  
      if(y_plus(i,j,k).le.yplam) fyplus = 0.0
      if(y_plus(i,j,k).gt.yplam)  then
               xfac = (y_plus(i,j,k) - yplam)/(ypturb- yplam)
               if(xfac.gt.1.0) xfac = 1.0
               fyplus = xfac*xfac*(3.0 - 2.0*xfac)
      end if
      turbvis_wall = turbvis_wall*fyplus
c
c     set the laminar viscosity if it varies with temperature.
c
      if(reyno.lt.0.0001) then
            vislam = (abs(reyno)/100000.0) * (tsavg/288.0)**0.62
      end if
c
c    save the viscosity it is only used to calculate and write out the reynolds number.
      if(i.eq.imid.and.k.eq.kmid.and.j.eq.jtedge) viscosy(nrw) = vislam
c
c     calculate the reynolds number at mid passage at the trailing edge.
c
      if(i.eq.imid.and.k.eq.kmid.and.j.eq.jtedge) then
           wabs       = sqrt(vx(i,j,k)*vx(i,j,k) + vr(i,j,k)*vr(i,j,k)
     &                + wt(i,j,k)*wt(i,j,k) )
           rovexit(nrw)  = ro(i,j,k)*wabs
           viscosy(nrw)  = vislam
           reynolds(nrw) = chord(nrw)*rovexit(nrw)/vislam
      end if
c
c  set the total viscosity including that due to free stream turbulence.
c
      vistot  = vislam*(1 + fst_rat(i,j,k) ) + turbvis_wall 
c
c    set a limit to the turbulent viscosity
c
      vislim = turbvis_lim*vislam
      if(vistot.gt.vislim) vistot = vislim
c
c     save the turbulent/laminar viscosity ratio as "visc_rat"
c
       visc_rat(i,j,k)   = vistot/vislam
       visc_lam(i,j,k)   = vislam
c
  150 continue
c
c    end of setting the turbulent viscosity
c****************************************************************************
c****************************************************************************
c     next  check for transition.
c
      do 160 j = 2,jm
c
      nrw    = nrow(j)
      j1     = jstart(nrw)
      jrow   = j - j1 + 1
      jtrhub = jtran_k1(nrw)
      jtrtip = jtran_km(nrw)
      jtrlow = jtran_i1(nrw)
      jtrup  = jtran_im(nrw)
c
c  q3d
      if(km.eq.2) go to 171
c  end q3d
c
c   first on the hub
c
      do 170 i = 1,imm1
      if(jrow.lt.jtrhub) go to 34
      rmax = 0.0
      do 37 k=1,kmid
      ratvis = visc_rat(i,j,k)
      if(ratvis.lt.rmax) go to 37
      rmax = ratvis
   37 continue
      if(rmax.gt.ftrans) go to 38
   34 continue
      do 39 k = 1,kmid
   39 visc_rat(i,j,k) = 1.0
   38 continue
c
c  next on the casing
c
      if(jrow.lt.jtrtip) go to 46
      rmax = 0.0
      do 47 k = kmid,kmm1
      ratvis = visc_rat(i,j,k)
      if(ratvis.lt.rmax) go to 47
      rmax = ratvis
   47 continue
      if(rmax.gt.ftrans) go to 49
   46 continue
      do 48 k=kmid,kmm1
   48 visc_rat(i,j,k) = 1.0
   49 continue
c
  170 continue
c
  171 continue
c
c    next on the lower - i=1, blade surface
c
      do 180  k=1,kmm1
      if(jrow.lt.jtrlow) go to 70
      rmax=0.0
      do 56 i=1,imid
      ratvis = visc_rat(i,j,k)
      if(ratvis.lt.rmax) go to 56
      rmax = ratvis
   56 continue
      if(rmax.gt.ftrans)  go to 71
   70 continue
      do 72 i = 1,imid
   72 visc_rat(i,j,k)= 1.0
   71 continue
c
c    next on the upper, i = im, blade surface
c
      if(jrow.lt.jtrup) go to 69
      rmax = 0.0
      do 73 i = imid,imm1
      ratvis  = visc_rat(i,j,k)
      if(ratvis.lt.rmax) go to 73
      rmax = ratvis
   73 continue
      if(rmax.gt.ftrans) go to  77
   69 continue
      do 74 i  = imid,imm1
   74 visc_rat(i,j,k) = 1.0
   77 continue
c
  180 continue
c
  160 continue
c
c********************************************************************************
c********************************************************************************
c
c     now the final viscosity is fixed set the stresses and heat flows in all elements
c     jdd warning the div v term in the normal stresses is not yet included.
c
      do 190 k = 1,kmm1
      do 190 j = 2,jm
      do 190 i = 1,imm1
c
      vislam  = visc_lam(i,j,k)
      vistot  = visc_rat(i,j,k)*vislam
      thcond  = cp*vistot/prandtl
      vistot2 = 2.0*vistot
c     
      txr(i,j,k) = vistot*txr(i,j,k)
      txx(i,j,k) = vistot2*txx(i,j,k)
      txt(i,j,k) = vistot*txt(i,j,k)
      trx(i,j,k) = txr(i,j,k)
      trr(i,j,k) = vistot2*trr(i,j,k)
      trt(i,j,k) = vistot*trt(i,j,k)
      ttx(i,j,k) = txt(i,j,k)
      ttt(i,j,k) = vistot2*ttt(i,j,k)
      ttr(i,j,k) = trt(i,j,k)
      qxx(i,j,k) = -thcond*qxx(i,j,k)
      qrr(i,j,k) = -thcond*qrr(i,j,k)
      qtt(i,j,k) = -thcond*qtt(i,j,k)
c
  190 continue
c
c******************************************************************************
c******************************************************************************
c    evaluate and smooth the pressure gradients if ypluswall is < -10.0.
c
      if(ypluswall.lt.-10.0) call set_pwallgrad
c
c***************************************************************************************
c****************************************************************************************
c     calculate the wall shear stresses on all the solid surfaces.
c
      vislam16  = vislam*16.0
c
c******************first for blade surfaces.****************************
c
      do 200 k=1,kmm1
      kp1 = k+1
      do 200 j=2,jm
c
c   skip if not on a blade surface
      if(ind(j).eq.0)  go to 220
c  
      jm1 = j-1
      nrw = nrow(j)
c 
c   also skip for cells above the tip gap
      if( (ktips(nrw).gt.0).and. 
     &    (k.ge.ktips(nrw)).and.(k.lt.ktipe(nrw)) ) 
     &    go to 220
c
c*********************first for the i = 1  blade surface.*******************************
c                  calculate the wall shear stress.
c
      vxavg1 = vx(2,j,k)+vx(2,jm1,k)+vx(2,jm1,kp1)+vx(2,j,kp1)
      vravg1 = vr(2,j,k)+vr(2,jm1,k)+vr(2,jm1,kp1)+vr(2,j,kp1)
      wtavg1 = wt(2,j,k)+wt(2,jm1,k)+wt(2,jm1,kp1)+wt(2,j,kp1)
      roavg1 = ro(2,j,k)+ro(2,jm1,k)+ro(2,jm1,kp1)+ro(2,j,kp1)
      wavg   = sqrt(wtavg1*wtavg1 + vxavg1*vxavg1 + vravg1*vravg1)
      area1  = sqrt(abx(1,j,k)*abx(1,j,k) + abr(1,j,k)*abr(1,j,k)
     &            + abt(j,k)*abt(j,k)) 
      dperp  = vol(1,j,k)/area1
c
c**********************************************************************************
c**********************************************************************************
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then
c
         yplus_old =  yplus_i1(j,k)
         density   = 0.25*roavg1
         wslip     = 0.25*wavg
         call wallfun(1,j,k,1,dperp,dpds_cell(1,j,k),density,
     &                twalli1,yplus_old,wslip,yplus_new)
          yplus_i1(j,k) = amin1(1000.0,yplus_new)
c
      go to 365
c
      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************
c
      re     = dperp*roavg1*wavg/vislam16
      relog  = 1./alog(re)
c******************************************************************************
c    allow for roughness
c
      rough  = rough_l(nrw)
c
      if(rough.gt.1.0e-7) then
           if(rough.gt.dperp) rough = dperp
           rek   = re*rough/dperp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if
c******************************************************************************
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c
c  check for fully rough surfaces
c
      if(rough.gt.1.0e-7) then
           twalli1    = cf*roavg1*wavg*wavg/128
           vstar   = sqrt(twalli1/(roavg1/4))
           plusk   = rough*vstar*roavg1/vislam/4
           relog10 = log10(dperp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
      end if
c
c******************************************************************************
c   set the shear stress on the i = 1 blade surface.
c
      twalli1    = cf*roavg1*wavg*wavg/128
c
c******************************************************************************
c
  365 continue
c
c******************************************************************************
c    calculate  yplus for use later.
      do i = 1,imid
      vstar         = sqrt(twalli1/(0.25*roavg1))
      yplsi1        = vstar*sqrt(dwallsq(i,j,k))*(0.25*roavg1)/vislam
      y_plus(i,j,k) = amin1(1000.0,yplsi1)
      end do  
c    
      fmult   = -twalli1*area1/wavg
c
      xforce1 =  fmult*vxavg1
      rforce1 =  fmult*vravg1
      tforce1 =  fmult*wtavg1
      wvisc1  =  tforce1*wrad(j)*ravg_cell(j,k)
c
      forcex(1,j,k)    = forcex(1,j,k)    + xforce1
      forcer(1,j,k)    = forcer(1,j,k)    + rforce1
      forcet(1,j,k)    = forcet(1,j,k)    + tforce1
      esource(1,j,k)   = esource(1,j,k)   + wvisc1
c
c*********************now for the i = im blade surface************************
c                     calculate the wall shear stress 
c
      vxavgim = vx(imm1,j,k)+vx(imm1,jm1,k)
     &        + vx(imm1,jm1,kp1)+vx(imm1,j,kp1)
      vravgim = vr(imm1,j,k)+vr(imm1,jm1,k)
     &        + vr(imm1,jm1,kp1)+vr(imm1,j,kp1)
      wtavgim = wt(imm1,j,k)+wt(imm1,jm1,k)
     &        + wt(imm1,jm1,kp1)+wt(imm1,j,kp1)
      roavgim = ro(imm1,j,k)+ro(imm1,jm1,k)
     &        + ro(imm1,jm1,kp1)+ro(imm1,j,kp1)
      wavg    = sqrt(wtavgim*wtavgim+vxavgim*vxavgim+vravgim*vravgim)
c 
      areaim  = sqrt(abx(im,j,k)*abx(im,j,k) + abr(im,j,k)*abr(im,j,k)
     &              +abt(j,k)*abt(j,k)) 
      dperp    = vol(imm1,j,k)/areaim
c
c******************************************************************************
c******************************************************************************
c  next the i = im , upper, surface
c
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then

           yplus_old =  yplus_im(j,k)
           density   = 0.25*roavgim
           wslip     = 0.25*wavg
           call wallfun(imm1,j,k,im,dperp,dpds_cell(imm1,j,k),
     &             density,twallim,yplus_old,wslip,yplus_new)
           yplus_im(j,k) = amin1(1000.0,yplus_new)
c
      go to 375
c
      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************
      re       = dperp*roavgim*wavg/vislam16
      relog    = 1./alog(re)
c******************************************************************************
c    allow for roughness
c
      rough  = rough_u(nrw)
c
      if(rough.gt.1.0e-7) then
           if(rough.gt.dperp) rough = dperp
           rek   = re*rough/dperp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if
c******************************************************************************
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c
c  check for fully rough surfaces
c
      if(rough.gt.1.0e-7) then
           twallim = cf*roavgim*wavg*wavg/128
           vstar   = sqrt(twallim/(roavgim/4))
           plusk   = rough*vstar*roavgim/vislam/4
           relog10 = log10(dperp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
      end if
c
c******************************************************************************
c   set the shear stress on the i = im blade surface
c
      twallim    = cf*roavgim*wavg*wavg/128
c
c******************************************************************************
c******************************************************************************
c 
  375 continue
c
c*****************************************************************************
c*****************************************************************************   
c     calculate  yplus for use later.
      do i = imid+1,imm1
      vstar         = sqrt(twallim/(0.25*roavgim))
      yplsim        = vstar*sqrt(dwallsq(i,j,k))*(0.25*roavgim)/vislam
      y_plus(i,j,k) = amin1(1000.0,yplsim)
      end do
c
      fmult    = -twallim*areaim/wavg
c 
      xforceim =  fmult*vxavgim
      rforceim =  fmult*vravgim
      tforceim =  fmult*wtavgim
      wviscim  =  tforceim*wrad(j)*ravg_cell(j,k)
c
      forcex(imm1,j,k)  = forcex(imm1,j,k)   + xforceim
      forcer(imm1,j,k)  = forcer(imm1,j,k)   + rforceim
      forcet(imm1,j,k)  = forcet(imm1,j,k)   + tforceim
      esource(imm1,j,k) = esource(imm1,j,k)  + wviscim
c
  220 continue
c
c*****************************************************************************
c
c     now set the forces on the periodic cells, i=1 and i= im-1,upstream of the
c     le and downstream of the te .
c     and also for cells in the tip gap.
c*****************************************************************************
c     do not treat cells on the solid blade surfaces as periodic.
c
       if( (ind(j).eq.1).and.(k.lt.ktips(nrw)).or.(k.gt.ktipe(nrw)) )
     &      go to 225
c
      i   = 1
      im1 = im-1
      xfor = 0.5*((txr(i,j,k)+txr(im1,j,k))*abr(i,j,k)
     &    +       (txt(i,j,k)+txt(im1,j,k))*abt(j,k)
     &    +       (txx(i,j,k)+txx(im1,j,k))*abx(i,j,k))
      forcex(im1,j,k) = forcex(im1,j,k) + xfor
      forcex(i,j,k)   = forcex(i,j,k)   - xfor
      avg = 0.5*(forcex(i,j,k) + forcex(im1,j,k))
      forcex(i,j,k)   = avg
      forcex(im1,j,k) = avg

      rfor = 0.5*((trr(i,j,k)+trr(im1,j,k))*abr(i,j,k)
     &    +       (trt(i,j,k)+trt(im1,j,k))*abt(j,k)
     &    +       (trx(i,j,k)+trx(im1,j,k))*abx(i,j,k))
      forcer(im1,j,k)  = forcer(im1,j,k) + rfor
      forcer(i,j,k)    = forcer(i,j,k)   - rfor
      avg = 0.5*(forcer(i,j,k) + forcer(im1,j,k))
      forcer(i,j,k)   = avg
      forcer(im1,j,k) = avg

      tfor = 0.5*((ttr(i,j,k)+ttr(im1,j,k))*abr(i,j,k)
     &    +       (ttt(i,j,k)+ttt(im1,j,k))*abt(j,k)
     &    +       (ttx(i,j,k)+ttx(im1,j,k))*abx(i,j,k))
      forcet(im1,j,k)  = forcet(im1,j,k) + tfor
      forcet(i,j,k)    = forcet(i,j,k)   - tfor
      avg = 0.5*(forcet(i,j,k) + forcet(im1,j,k))
      forcet(i,j,k)   = avg
      forcet(im1,j,k) = avg
c
c     set the heat flow and viscous work
c
      qflow =  0.5*((qxx(i,j,k) + qxx(im1,j,k))*abx(i,j,k)
     &      +       (qrr(i,j,k) + qrr(im1,j,k))*abr(i,j,k)
     &      +       (qtt(i,j,k) + qtt(im1,j,k))*abt(j,k))
c
      vxavg = vx(i,j,k)+vx(i,jm1,k)+vx(i,jm1,kp1)+vx(i,j,kp1)
      vravg = vr(i,j,k)+vr(i,jm1,k)+vr(i,jm1,kp1)+vr(i,j,kp1)
      vtavg = vt(i,j,k)+vt(i,jm1,k)+vt(i,jm1,kp1)+vt(i,j,kp1)
c
      wvisc = -0.25*(xfor*vxavg + rfor*vravg + tfor*vtavg)
c
      esource(i,j,k)  = esource(i,j,k)   + qflow + wvisc
      esource(im1,j,k)= esource(im1,j,k) - qflow - wvisc 
      avg = 0.5*(esource(i,j,k) + esource(im1,j,k))
      esource(i,j,k)   = avg
      esource(im1,j,k) = avg
c
  225 continue
c
c******************************************************************************
c******************************************************************************
c
c     now sum the stresses on the other i = constant faces 
c
      do 250 i=2,imm1
      im1  = i-1
      xfor = 0.5*((txr(i,j,k)+txr(im1,j,k))*abr(i,j,k)
     &    +       (txt(i,j,k)+txt(im1,j,k))*abt(j,k)
     &    +       (txx(i,j,k)+txx(im1,j,k))*abx(i,j,k))
      forcex(im1,j,k) = forcex(im1,j,k) + xfor
      forcex(i,j,k)   = forcex(i,j,k)   - xfor

      rfor = 0.5*((trr(i,j,k)+trr(im1,j,k))*abr(i,j,k)
     &    +       (trt(i,j,k)+trt(im1,j,k))*abt(j,k)
     &    +       (trx(i,j,k)+trx(im1,j,k))*abx(i,j,k))
      forcer(im1,j,k)  = forcer(im1,j,k) + rfor
      forcer(i,j,k)    = forcer(i,j,k)   - rfor

      tfor = 0.5*((ttr(i,j,k)+ttr(im1,j,k))*abr(i,j,k)
     &    +       (ttt(i,j,k)+ttt(im1,j,k))*abt(j,k)
     &    +       (ttx(i,j,k)+ttx(im1,j,k))*abx(i,j,k))
      forcet(im1,j,k)  = forcet(im1,j,k) + tfor
      forcet(i,j,k)    = forcet(i,j,k)   - tfor
c
c     set the heat flow and viscous work
c
      qflow =  0.5*((qxx(i,j,k) + qxx(im1,j,k))*abx(i,j,k)
     &      +       (qrr(i,j,k) + qrr(im1,j,k))*abr(i,j,k)
     &      +       (qtt(i,j,k) + qtt(im1,j,k))*abt(j,k))
c
      vxavg = vx(i,j,k)+vx(i,jm1,k)+vx(i,jm1,kp1)+vx(i,j,kp1)
      vravg = vr(i,j,k)+vr(i,jm1,k)+vr(i,jm1,kp1)+vr(i,j,kp1)
      vtavg = vt(i,j,k)+vt(i,jm1,k)+vt(i,jm1,kp1)+vt(i,j,kp1)
c
      wvisc = -0.25*(xfor*vxavg + rfor*vravg + tfor*vtavg)
c
      esource(i,j,k)  = esource(i,j,k)   + qflow + wvisc
      esource(im1,j,k)= esource(im1,j,k) - qflow - wvisc 

  250 continue

  200 continue
c
c**************************************************************************************
c***************************************************************************************
c**********calculate the wall shear stresses on the hub and casing surfaces**********
c
c   q3d
      if(km.eq.2) go to 401
c   end q3d
c
      do 400 j=2,jm
      jm1  = j-1
      nrw  = nrow(j)
      do 400 i=1,imm1
      ip1 = i+1
c
c     first for the k = 1,  hub, wall.
c
c     calculate the wall shear stress.
c
      vxavg1 = vx(i,j,2)+vx(ip1,j,2)+vx(ip1,jm1,2)+vx(i,jm1,2)
      vravg1 = vr(i,j,2)+vr(ip1,j,2)+vr(ip1,jm1,2)+vr(i,jm1,2)
      vtavg1 = vt(i,j,2)+vt(ip1,j,2)+vt(ip1,jm1,2)+vt(i,jm1,2)
      roavg1 = ro(i,j,2)+ro(ip1,j,2)+ro(ip1,jm1,2)+ro(i,jm1,2)
      ravg1  = r(j,2)   + r(j,2)    + r(jm1,2)    + r(jm1,2)
      wtavg1 = vtavg1 - whub(j)*ravg1
      wavg   = sqrt(vxavg1*vxavg1 + vravg1*vravg1 + wtavg1*wtavg1)
      area1  = sqrt(asx(i,j,1)*asx(i,j,1) + asr(i,j,1)*asr(i,j,1))
      dperp  = vol(i,j,1)/area1
c
c**********************************************************************************
c**********************************************************************************
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then

           yplus_old =  yplus_k1(i,j)
           density   = 0.25*roavg1
           wslip     = 0.25*wavg
         call wallfun(i,j,1,1,dperp,dpds_cell(i,j,1),density,
     &                twallk1,yplus_old,wslip,yplus_new)
           yplus_k1(i,j) = amin1(1000.0,yplus_new)
c
      go to 345

      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************
      re     = dperp*roavg1*wavg/vislam16
      relog  = 1./alog(re)
c******************************************************************************
c    allow for roughness
c
      rough  = rough_h(nrw)
c
      if(rough.gt.1.0e-7) then
           if(rough.gt.dperp) rough = dperp
           rek   = re*rough/dperp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if
c
c******************************************************************************
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c
c  check for fully rough surfaces
c
      if(rough.gt.1.0e-7) then
           twallk1     = cf*roavg1*wavg*wavg/128
           vstar   = sqrt(twallk1/(roavg1/4))
           plusk   = rough*vstar*roavg1/vislam/4
           relog10 = log10(dperp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
      end if
c
c******************************************************************************
c  set the shear stress on the k=1, (hub) , endwall .
c
      twallk1    = cf*roavg1*wavg*wavg/128
c
c******************************************************************************
c
  345 continue
c
c******************************************************************************
c   set yplus for use later
      if(ibound.eq.0.or.ibound.eq.2) then
      do k = 1,kmid
      vstar         = sqrt(twallk1/(0.25*roavg1))
      yplsk1        = vstar*sqrt(dwallsq(i,j,k))*(0.25*roavg1)/vislam
      yplsp         = y_plus(i,j,k)
      y_plus(i,j,k) = amin1(yplsk1,yplsp)
      end do
      end if
c
c******************************************************************************
c
      fmult   = -twallk1*area1/wavg
      if(ibound.eq.1.or.ibound.gt.2) fmult = 0.0
c
      xforce1 = fmult*vxavg1
      rforce1 = fmult*vravg1
      tforce1 = fmult*wtavg1
      wvisc1  = tforce1*0.25*whub(j)*ravg1
c
      forcex(i,j,1)    = forcex(i,j,1)    + xforce1
      forcer(i,j,1)    = forcer(i,j,1)    + rforce1
      forcet(i,j,1)    = forcet(i,j,1)    + tforce1
      esource(i,j,1)   = esource(i,j,1)   + wvisc1
c
c
c**************now for the k = km (casing) end wall**********************************
c
c     calculate the wall shear stress.
c
      vxavgkm = vx(i,j,kmm1)+vx(ip1,j,kmm1)
     &        + vx(ip1,jm1,kmm1)+vx(i,jm1,kmm1)
      vravgkm = vr(i,j,kmm1)+vr(ip1,j,kmm1)
     &        + vr(ip1,jm1,kmm1)+vr(i,jm1,kmm1)
      vtavgkm = vt(i,j,kmm1)+vt(ip1,j,kmm1)
     &        + vt(ip1,jm1,kmm1)+vt(i,jm1,kmm1)
      roavgkm = ro(i,j,kmm1)+ro(ip1,j,kmm1)
     &        + ro(ip1,jm1,kmm1)+ro(i,jm1,kmm1)
      ravgkm  = r(j,kmm1) + r(j,kmm1) + r(jm1,kmm1) + r(jm1,kmm1)
      wtavgkm = vtavgkm   - wtip(j)*ravgkm
      wavg    = sqrt(vxavgkm*vxavgkm+vravgkm*vravgkm+wtavgkm*wtavgkm)
      areakm  = sqrt(asx(i,j,km)*asx(i,j,km) + asr(i,j,km)*asr(i,j,km))
      dperp   = vol(i,j,kmm1)/areakm
c************************************************************************
c****************************************************************************
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then

           yplus_old  =  yplus_km(i,j)
           density    = 0.25*roavgkm
           wslip      = 0.25*wavg
        call wallfun(i,j,kmm1,km,dperp,dpds_cell(i,j,kmm1),density,
     &              twallkm,yplus_old,wslip,yplus_new)
           yplus_km(i,j) = amin1(1000.0,yplus_new)
c
      go to 355
c
      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************

      re      = dperp*roavgkm*wavg/vislam16
      relog   = 1./alog(re)
c******************************************************************************
c    allow for roughness
c
      rough  = rough_t(nrw)
c
      if(rough.gt.1.0e-7) then
           if(rough.gt.dperp) rough = dperp
           rek   = re*rough/dperp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if
c
c******************************************************************************
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c
c  check for fully rough surfaces
c
      if(rough.gt.1.0e-7) then
           twallkm = cf*roavgkm*wavg*wavg/128
           vstar   = sqrt(twallkm/(roavgkm/4))
           plusk   = rough*vstar*roavgkm/vislam/4
           relog10 = log10(dperp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
      end if
c
c******************************************************************************
c  set the shear stress on the casing .
c
      twallkm    = cf*roavgkm*wavg*wavg/128
c
c******************************************************************************
c
  355 continue
c
c******************************************************************************
c     set yplus for use later
      if(ibound.le.1) then
      do k = kmid+1,kmm1
      vstar         = sqrt(twallkm/(0.25*roavgkm))
      yplskm        = vstar*sqrt(dwallsq(i,j,k))*(0.25*roavgkm)/vislam
      yplsp         = y_plus(i,j,k)
      y_plus(i,j,k) = amin1(yplskm,yplsp)
      end do
      end if
c 
      fmult   = -twallkm*areakm/wavg
      if(ibound.ge.2) fmult = 0.0 
c
      xforcekm = fmult*vxavgkm
      rforcekm = fmult*vravgkm
      tforcekm = fmult*wtavgkm
      wvisckm  =  tforcekm*0.25*wtip(j)*ravgkm
c
      forcex(i,j,kmm1) = forcex(i,j,kmm1) + xforcekm
      forcer(i,j,kmm1) = forcer(i,j,kmm1) + rforcekm
      forcet(i,j,kmm1) = forcet(i,j,kmm1) + tforcekm
      esource(i,j,kmm1)= esource(i,j,kmm1)+ wvisckm
c
c     now sum the stresses on the other k = constant faces
c
      do 450 k=2,kmm1
      km1  = k-1
      xfor = 0.5*((txr(i,j,k)+txr(i,j,km1))*asr(i,j,k)
     &    +       (txx(i,j,k)+txx(i,j,km1))*asx(i,j,k))
      forcex(i,j,km1) = forcex(i,j,km1) + xfor
      forcex(i,j,k)   = forcex(i,j,k)   - xfor

      rfor = 0.5*((trr(i,j,k)+trr(i,j,km1))*asr(i,j,k)
     &    +       (trx(i,j,k)+trx(i,j,km1))*asx(i,j,k))
      forcer(i,j,km1)  = forcer(i,j,km1) + rfor
      forcer(i,j,k)    = forcer(i,j,k)   - rfor

      tfor = 0.5*((ttr(i,j,k)+ttr(i,j,km1))*asr(i,j,k)
     &    +       (ttx(i,j,k)+ttx(i,j,km1))*asx(i,j,k))
      forcet(i,j,km1)  = forcet(i,j,km1) + tfor
      forcet(i,j,k)    = forcet(i,j,k)   - tfor
c
c     set the heat flow and viscous work
c
      qflow =  0.5*((qxx(i,j,k) + qxx(i,j,km1))*asx(i,j,k)
     &      +       (qrr(i,j,k) + qrr(i,j,km1))*asr(i,j,k))
c
      vxavg = vx(i,j,k)+vx(ip1,j,k)+vx(ip1,jm1,k)+vx(i,jm1,k)
      vravg = vr(i,j,k)+vr(ip1,j,k)+vr(ip1,jm1,k)+vr(i,jm1,k)
      vtavg = vt(i,j,k)+vt(ip1,j,k)+vt(ip1,jm1,k)+vt(i,jm1,k)

      wvisc = -0.25*(xfor*vxavg + rfor*vravg + tfor*vtavg)

      esource(i,j,k)  = esource(i,j,k)   + qflow + wvisc
      esource(i,j,km1)= esource(i,j,km1) - qflow - wvisc 
c
  450 continue
c
  400 continue
c
  401 continue
c
c******************************************************************************
c
c     all the viscous forces are now set. use them to update the
c     global body force terms.
c
c******************************************************************************
c     note the negative sign because the source term is subtracted in subroutine tstep.
c 
      rf_vis1 = 1.0 - rf_vis
      do 500 k=1,kmm1
      do 500 j=2,jm
      do 500 i=1,imm1
      xforce(i,j,k)   = rf_vis1*xforce(i,j,k)   - rf_vis*forcex(i,j,k)
      tforce(i,j,k)   = rf_vis1*tforce(i,j,k)   - rf_vis*forcet(i,j,k)
     &                 *ravg_cell(j,k)
      rforce(i,j,k)   = rf_vis1*rforce(i,j,k)   - rf_vis*forcer(i,j,k)
      qsource(i,j,k)  = rf_vis1*qsource(i,j,k)  - rf_vis*esource(i,j,k)
  500 continue
c
c******************************************************************************
c******************************************************************************
c
      return
      end
c
c******************************************************************************
c
c*************************************************************************************
c
      subroutine gradvel(i,j,k, v, dvdx, dvdr, dvdt)
c
c******************************************************************************************
c     this subroutine calculates the derivatives of the input function  "v " 
c     for a cell with corner storage using gauss' theorem .
c
c     because the areas of the faces are evaluated relative to the radial direction 
c     at the middle of the face  an extra  "v/r"  term must be subtracted from the 
c     radial derivative.
c******************************************************************************************
c
      include 'commall-open-19.2'
c
      dimension  v(id,jd,kd)
c
      ip1 = i+1
      jm1 = j-1
      kp1 = k+1
c  note that the j index of a cell is the j value on its right hand side. 
c  all areas are positive in the direction of an inwards normal.    
      rvol    = 0.25/vol(i,j,k)
c
      vavgj1  = v(i,jm1,k)+ v(ip1,jm1,k)+ v(ip1,jm1,kp1)+ v(i,jm1,kp1)
      vavgj2  = v(i,j,k)  + v(ip1,j,k)  + v(ip1,j,kp1)  + v(i,j,kp1)
      vavgi1  = v(i,j,k)  + v(i,jm1,k)  + v(i,jm1,kp1)  + v(i,j,kp1)
      vavgi2  = v(ip1,j,k)+ v(ip1,jm1,k)+ v(ip1,jm1,kp1)+ v(ip1,j,kp1)
      vavgk1  = v(i,j,k)  + v(ip1,j,k)  + v(ip1,jm1,k)  + v(i,jm1,k)
      vavgk2  = v(i,j,kp1)+ v(ip1,j,kp1)+ v(ip1,jm1,kp1)+ v(i,jm1,kp1)
      vavgcell= 0.125*(vavgj1 + vavgj2)
c
      dvdx = -(vavgj2*aqx(i,j,k)   - vavgj1*aqx(i,jm1,k)
     &       + vavgi1*abx(i,j,k)   - vavgi2*abx(ip1,j,k)
     &       + vavgk1*asx(i,j,k)   - vavgk2*asx(i,j,kp1))*rvol
c
      dvdr = -(vavgj2*aqr(i,j,k)   - vavgj1*aqr(i,jm1,k)
     &       + vavgi1*abr(i,j,k)   - vavgi2*abr(ip1,j,k)
     &       + vavgk1*asr(i,j,k)   - vavgk2*asr(i,j,kp1))*rvol
     &       - vavgcell/ravg_cell(j,k)
c
      dvdt = -(vavgi1  - vavgi2)*abt(j,k)*rvol
c
      return
      end
c*************************************************************************************
c*************************************************************************************
c
c*************************************************************************************
c
      subroutine gradcell(i,j,k, v, dvdx, dvdr, dvdt)
c
c******************************************************************************************
c     this subroutine calculates the derivatives of the input function  "v" 
c     for a cell with cell centre storage using gauss' theorem, where  "v" is the average
c     value for the cell .
c
c     because the areas of the faces are evaluated relative to the radial direction 
c     at the middle of the face  an extra  "v/r"  term must be subtracted from the 
c     radial derivative.
c******************************************************************************************
c******************************************************************************************
c
      include 'commall-open-19.2'
c
      dimension  v(id,jd,kd)
c
      ip1 = i+1
      jp1 = j+1
      kp1 = k+1
      im1 = i-1
      jm1 = j-1
      km1 = k-1
      if(i.eq.1)    im1 = imm1
      if(i.eq.imm1) ip1 = 1
      if(k.eq.1)    km1 = 1
      if(k.eq.kmm1) kp1 = kmm1
c  note that the j index of a cell is the j value on its right hand side. 
c  all areas are positive in the direction of an inwards normal.  
c   1/rvol is x 0.5  to allow for averaging.
c  
      rvol    = 0.5/vol(i,j,k)
c
      vavgj1  = v(i,j,k)  + v(i,jm1,k)
      vavgj2  = v(i,j,k)  + v(i,jp1,k)
      vavgi1  = v(i,j,k)  + v(im1,j,k)
      vavgi2  = v(i,j,k)  + v(ip1,j,k)
      vavgk1  = v(i,j,k)  + v(i,j,km1)
      vavgk2  = v(i,j,k)  + v(i,j,kp1)
c   set zero value on blade surfaces
      if(ind(j).eq.1) then
           if(i.eq.1)    vavgi1 = 0.0
           if(i.eq.imm1) vavgi2 = 0.0
      end if
c   set zero value on hub and casing
c  jdd changed to allow for no shear on hub or casing. jan/14                
      if(k.eq.1.and.(ibound.eq.0.or.ibound.eq.2)) vavgk1 = 0.0
      if(k.eq.kmm1.and.ibound.le.1)               vavgk2 = 0.0
c
      dvdx = -(vavgj2*aqx(i,j,k)   - vavgj1*aqx(i,jm1,k)
     &       + vavgi1*abx(i,j,k)   - vavgi2*abx(ip1,j,k)
     &       + vavgk1*asx(i,j,k)   - vavgk2*asx(i,j,kp1))*rvol
c
      dvdr = -(vavgj2*aqr(i,j,k)   - vavgj1*aqr(i,jm1,k)
     &       + vavgi1*abr(i,j,k)   - vavgi2*abr(ip1,j,k)
     &       + vavgk1*asr(i,j,k)   - vavgk2*asr(i,j,kp1))*rvol
     &       - v(i,j,k)/ravg_cell(j,k)
c
      dvdt = -(vavgi1  - vavgi2)*abt(j,k)*rvol
c
c*************************************************************************************
c*************************************************************************************
c
      return
      end
c
c*************************************************************************************
c*************************************************************************************
c*************************************************************************************
c
      subroutine spal_loss
c
c*****************************************************************************************
c     this subroutine computes a viscous body force based on wall functions for the
c     surface shear stress and the spalart - allmaras turbulence model.          
c******************************************************************************************
c
      include  'commall-open-19.2'
c
      common/bkstress/  txx(id,jd,kd),txr(id,jd,kd),
     &                  txt(id,jd,kd),trx(id,jd,kd),trr(id,jd,kd),
     &                  trt(id,jd,kd),ttx(id,jd,kd),ttr(id,jd,kd),
     &                  ttt(id,jd,kd),qxx(id,jd,kd),qrr(id,jd,kd),
     &                  qtt(id,jd,kd)
c
      dimension  forcex(id,jd,kd),forcer(id,jd,kd),forcet(id,jd,kd),
     &           esource(id,jd,kd)
c
      rf_vis1 = 1.0 - rf_vis
      rf_vis4 = 0.25*rf_vis
      ndamp   = 100000
c
c      calculate the viscosity over the first quarter of the steps.
c      then hold it constant for the remainder of the steps.
c
      if(nstep.eq.1.or.nstep.lt.nsteps_max/4) then
c
      if(reyno.gt.100.) then
           j1=jle(1)
           j2=jte(1)
           if(jle(1).gt.jm) j1=1
           if(jte(1).gt.jm) j2=jm
           xchord = smerid(j2,kmid)-smerid(j1,kmid)
           row2   = sqrt(rovx(imid,j2,kmid)*rovx(imid,j2,kmid)
     &            + rovr(imid,j2,kmid)*rovr(imid,j2,kmid)
     &            + rowt(imid,j2,kmid)*rowt(imid,j2,kmid) )
           vislam = xchord*row2/reyno
      end if
c
      if(reyno.gt.0.0.and.reyno.lt.99.99) then
            vislam = reyno/100000.
      end if
      if(reyno.lt.0.0) vislam = -reyno*1.0e-5
c
c     energy eqn addition
c
      thcond  = cp*vislam/prandtl
      ftcond = thcond/vislam
c
      endif
c
c    end of setting the viscosity, etc over the first quarter of the steps.
c
c********************************************************************************
c     initially set the viscous forces and energy sources to zero
c
      do 100 k=1,kmm1
      do 100 j=2,jm
      do 100 i=1,imm1
      forcex(i,j,k)  = 0.0
      forcer(i,j,k)  = 0.0
      forcet(i,j,k)  = 0.0
      esource(i,j,k) = 0.0
  100 continue
c
c*******************************************************************************************
c*******************************************************************************************
c     evaluate the turbulent viscosity for each cell in the do 150 loop.
c     only use the pressure gradient term if fac_pgrad is greater than zero.
c
      if(fac_pgrad.gt.0.001) then
          call set_pwallgrad
      end if
c
      do 1000 k=1,kmm1
      kp1 = k+1
      do 1000 j=2,jm
      jm1 = j-1
      nrw = nrow(j)
      jtedge = jte(nrw)
      jledge = jle(nrw)
      do 1000 i=1,imm1
      ip1 = i+1
c
c    average the conditions for the cells. note these are true averages.
c
      wtavg = 0.125*(wt(i,j,k)+wt(ip1,j,k)+wt(ip1,j,kp1)+wt(i,j,kp1)
     &      + wt(i,jm1,k)+wt(ip1,jm1,k)+wt(ip1,jm1,kp1)+wt(i,jm1,kp1))
      roavg = roavg_cell(i,j,k)
      tsavg = 0.125*(t_static(i,j,k)+t_static(ip1,j,k)
     &        +t_static(ip1,j,kp1)+t_static(i,j,kp1)+t_static(i,jm1,k)
     &        +t_static(ip1,jm1,k)+t_static(ip1,jm1,kp1)
     &        +t_static(i,jm1,kp1))
      ravg  = ravg_cell(j,k)
c
c     calculate the derivatives of the velocity components and temperature.
c     note the vorticity is based on the relative velocity,  wt.
c     this was changed from the absolute velocity by jdd august 2017.
      call gradvel(i,j,k,vx,dvxdx,dvxdr,dvxdt)
      call gradvel(i,j,k,vr,dvrdx,dvrdr,dvrdt)
      call gradvel(i,j,k,wt,dwtdx,dwtdr,dwtdt)
c
      call gradvel(i,j,k,t_static,dtempdx,dtempdr,dtempdt)
c
c     calculate the rates of strain for use when calculating the viscous stresses.
c     these are not yet the stresses as the final viscosity is not yet known.
c     
      txr(i,j,k) = dvxdr + dvrdx
      txx(i,j,k) = dvxdx
      txt(i,j,k) = dwtdx + dvxdt
      trr(i,j,k) = dvrdr
      trt(i,j,k) = dwtdr + dvrdt - wtavg/ravg
      ttt(i,j,k) = dwtdt
      qxx(i,j,k) = dtempdx
      qrr(i,j,k) = dtempdr
      qtt(i,j,k) = dtempdt
c
c        use the vorticity to determine the turbulent viscosity for each cell.
c
      vortx   = dvrdt - dwtdr - wtavg/ravg 
      vortr   = dwtdx - dvxdt
      vortt   = dvxdr - dvrdx
      vort    = sqrt(vortx*vortx + vortr*vortr + vortt*vortt)
      if(j.eq.jstart(nrw)) vort = 0.0
c
c   set a limit to the vorticity
      if(vort.gt.vort_max) vort = vort_max
c
      if(indmix(j).eq.1) then
            vort = abs_vort(i,j-1,k)
      end if
c
      abs_vort(i,j,k) = vort
c******************************************************************************
c******************************************************************************
c
c     calculate the factor to increase the source term due to streamwise vorticity.
c     as used by lee, wilson & vahdati.
c
      if(fac_vort.gt.0.001.or.fac_pgrad.gt.0.001) then
c
      wtavg = (wt(i,j,k)+wt(ip1,j,k)+wt(ip1,j,kp1)+wt(i,j,kp1)
     &      + wt(i,jm1,k)+wt(ip1,jm1,k)+wt(ip1,jm1,kp1)+wt(i,jm1,kp1))
      vxavg = (vx(i,j,k)+vx(ip1,j,k)+vx(ip1,j,kp1)+vx(i,j,kp1)
     &      + vx(i,jm1,k)+vx(ip1,jm1,k)+vx(ip1,jm1,kp1)+vx(i,jm1,kp1))
      vravg =(vr(i,j,k)+vr(ip1,j,k)+vr(ip1,j,kp1)+vr(i,j,kp1)
     &      + vr(i,jm1,k)+vr(ip1,jm1,k)+vr(ip1,jm1,kp1)+vr(i,jm1,kp1)) 
      wrel  = sqrt(vxavg*vxavg + vravg*vravg + wtavg*wtavg)
c
      end if
c
      if(fac_vort.gt.0.001) then
c  
      helicity  = (vxavg*vortx + vravg*vortr + wtavg*vortt)
     &            /(wrel*vort)
      helfac    = abs(helicity)
c
c    the next commented out line is the value of fac_vort used by lee, wislon & vahdati.
c    vorticity_fac  = fac_vort*(1.0 + 0.9191*tanh(3.0*helfac*helfac))
      vorticity_fac  = (1.0 + fac_vort*tanh(3.0*helfac*helfac))
      if(vort.lt.0.001*vort_max) vorticity_fac = 1.0
      sparevar(i,j,k) = vorticity_fac
c    
      end if
c
c******************************************************************************
c******************************************************************************
c
c     now set the pressure gradient term as used by lee, wislon & vahdati .
c
      if(fac_pgrad.gt.0.001) then
c
c  set a limiting relative velocity = 0.1 x the local mid point value.
      wlim = 0.1*sqrt(vx(imid,j,k)*vx(imid,j,k)
     &              + vr(imid,j,k)*vr(imid,j,k)
     &              + wt(imid,j,k)*wt(imid,j,k))
c   jdd correctin 20/9/18.  the factor of 0.125  should be used before checking against wlin.
      wrel  = 0.125*wrel
      if(wrel.lt.wlim) wrel = wlim
      pgrad = dpds_cell(i,j,k)
c
c   set sparevar to plot out the pressure gradient field.
c      sparevar(i,j,k) = pgrad/1.0e6
c
c    set the pressure gradient term to zero for favourable pressure gradients.
      if(pgrad.lt.0.0) pgrad = 0.0
c
      pterm = pgrad*vislam/(roavg*roavg*wrel*wrel*wrel)
c
c   note the pressure gradient term is scaled by reynolds number so that it
c   becomes independent of viscosity. this is different to lee, wislon & vahdati .
      pterm = pterm*reynolds(nrw)
c
c  the next coommented out line is the vaue of fac_pgrad used by lee, wilson & vahdati.
c      pgrad_fac = fac_pgrad*(1.0 + 0.6565*tanh(pterm*pterm) )
      pgrad_fac = (1.0 + fac_pgrad*tanh(pterm*pterm) )
      if(fac_vort.lt.0.001)  sparevar(i,j,k) = pgrad_fac
      if(fac_vort.gt.0.001) sparevar(i,j,k)  = sparevar(i,j,k)*pgrad_fac
c
      end if     
       
c******************************************************************************
c     set a turbulent viscosity using the mixing length model
c
      turbvis_wall = roavg*xlength(i,j,k)*vort*fmixup
c
c******************************************************************************
c     set the laminar viscosity if it is to vary with temperature
      if(reyno.lt.0.0001) then
            vislam = (abs(reyno)/100000.0) * (tsavg/288.0)**0.62 
      end if
c
c******************************************************************************
c    save the viscosity it is only used to calculate and write out the reynolds number.
      if(i.eq.imid.and.k.eq.kmid.and.j.eq.jtedge) viscosy(nrw)  = vislam
c******************************************************************************
c     calculate the reynolds number at mid passage at the trailing edge.
c
      if(i.eq.imid.and.k.eq.kmid.and.j.eq.jtedge) then
           wabs       = sqrt(vx(i,j,k)*vx(i,j,k) + vr(i,j,k)*vr(i,j,k)
     &                     + wt(i,j,k)*wt(i,j,k) )
           rovexit(nrw)  = ro(i,j,k)*wabs
           viscosy(nrw)  = vislam
           reynolds(nrw) = chord(nrw)*rovexit(nrw)/vislam
      end if
c******************************************************************************
c     add the free stream turbulent viscosity
c
      vistot  = vislam*(1 + fst_rat(i,j,k) ) + turbvis_wall 
c
c******************************************************************************
c    set a limit to the turbulent viscosity
c
      vislim = turbvis_lim*vislam
      if(vistot.gt.vislim) vistot = vislim
      visklam  = vislam/roavg
c
c************************************************************
c     save the mixing length turbulent dynamic viscosity. 
c
      dyn_vis_ml(i,j,k) = (vistot - vislam) 
c
c******************************************************************************
c******************************************************************************
c     now set the initial guess of viscosity for the s-a model when nstep = 1.
c
      if(if_restart.eq.0.and.nstep.eq.1)then
c
           if(fsturb(1).lt.1.0) then 
           transvisin  = 4.35*vislam*fsturb(1)**0.25
           else
           transvisin  = vislam*(3.5 + fsturb(1)*0.85)
           end if
c
          trans_dyn_vis(i,j,k) = transvisin
c
      end if 
c
      if(nstep.eq.1) trans_kvis(i,j,k) = trans_dyn_vis(i,j,k)/roavg
c
c******************************************************************************
c  now set the final value of turbulent viscosity
c
      transkvis        = trans_kvis(i,j,k)
      xi               = transkvis/visklam
c
      if(xi.lt.1.1) xi = 1.1
      if(xi.gt.turbvis_lim) xi = turbvis_lim
c
      xi3              = xi*xi*xi
      fv1              = xi3/(xi3 + 357.9)
      transvislim      = vislim/roavg/fv1
      if(transkvis.gt.transvislim) transkvis = transvislim
c
c*******************************************************************************************
c*******************************************************************************************
c   this is the main result of the sa model. vistru is the true turbulent dynamic viscosity.
c
      vistru = transkvis*fv1*roavg
c
c*******************************************************************************************
c
      fv2              = 1.0 - xi/(1.0 + xi*fv1)
c 
      if(xi.lt.2.0.and.fv2.gt.0.0) fv2 = 0.0
c
c      trans_vort(i,j,k) = abs_vort(i,j,k) 
c     &                  + 5.9488*fv2*transkvis/dwallsq(i,j,k)
c     
c   jdd changed this jan/14.  use the true vorticity rather than the transformed one.
c   this is different to the original sa model which is commented out above.
c
      trans_vort(i,j,k) = abs_vort(i,j,k)
c      
      vistot           = vistru + vislam
c
      vislim = turbvis_lim*vislam
      if(vistot.gt.vislim) vistot = vislim
c
c     save the turbulent/laminar viscosity ratio as "visc_rat" and the laminar viscosity as "visc_lam".
c     these and  vistot are the dynamic viscosities. 
c
       visc_rat(i,j,k)   = vistot/vislam
c      save the viscosity ratio for the output file.
       if(fac_vort.lt.0.001.and.fac_pgrad.lt.0.001)
     &    sparevar(i,j,k) = visc_rat(i,j,k)
c
       visc_lam(i,j,k)   = vislam
       visc_turb(i,j,k)  = vistot - vislam
c
 1000 continue
c
c******************************************************************************
c******************************************************************************
c   form the grad of turbulent viscosity for the source terms.
c
      do 1500,k=1,kmm1
      do 1500 j=2,jmm1
      nrw = nrow(j)
      do 1500,i=1,imm1
      call gradcell(i,j,k,trans_kvis,d_tvisdx,d_tvisdr,d_tvisdt)
      grad_tvis(i,j,k) = d_tvisdx*d_tvisdx + d_tvisdr*d_tvisdr
     &                 + d_tvisdt*d_tvisdt
      visceff          = trans_kvis(i,j,k) + visc_lam(i,j,k)
      tvisgradx(i,j,k) = visceff*d_tvisdx
      tvisgradr(i,j,k) = visceff*d_tvisdr
      tvisgradt(i,j,k) = visceff*d_tvisdt
      if(j.eq.jstart(nrw).or.j.eq.jstart(nrw)-1) then
      grad_tvis(i,j,k) = 0.0
      tvisgradx(i,j,k) = 0.0
      tvisgradr(i,j,k) = 0.0
      tvisgradt(i,j,k) = 0.0
      end if
 1500 continue
c******************************************************************************
c
c   
c    calculate the source terms
c
      cw2    = 0.3
      sigma  = 0.6667
      cb2    = 0.622
      cb1    = 0.1355
      cw1    = 0.80607 + 2.433
      fst0   = fac_st0*cb1
      fst1   = fac_st1/sigma
      fst2   = fac_st2*cb2/sigma
      fst3   = fac_st3*cw1
      fstmix = fac_stmix/ndamp
c
      do 1600 k=1,kmm1
      do 1600 j=2,jm
      do 1600 i=1,imm1
c
c   jdd jan/14 set fyplus to reduce the source term in the sublayer and buffer regions.
c
      if(y_plus(i,j,k).le.yplam) fyplus = 0.0
      if(y_plus(i,j,k).gt.yplam)  then
               xfac = (y_plus(i,j,k) - yplam)/(ypturb- yplam)
               if(xfac.gt.1.0) xfac = 1.0
               fyplus = xfac*xfac*(3.0 - 2.0*xfac)
      end if
c
c  end jan/14 addition.

      call gradcell(i,j,k,tvisgradx,tvisddx,dumy,dummy2) 
      call gradcell(i,j,k,tvisgradr,dumy,tvisddr,dummy2) 
      call gradcell(i,j,k,tvisgradt,dumy,dummy2,tvisddt) 
c
      roavg   = roavg_cell(i,j,k)
c
      tvis    = trans_kvis(i,j,k)

      visklam = visc_lam(i,j,k)/roavg
c
      st0     = fst0*trans_vort(i,j,k)*tvis
      if(fac_vort.gt.0.001)  st0 = st0*vorticity_fac
      if(fac_pgrad.gt.0.001) st0 = st0*pgrad_fac
c
c    jdd addition jan/14.  reduce the source term in the sublayer and buffer layer.
      st0     = st0*fyplus
c  end jan/14 addition
c
      st1     = fst1*(tvisddx + tvisddr + tvisddt)
c
      st2     = fst2*grad_tvis(i,j,k)
c
      tau     = 5.9488*tvis/trans_vort(i,j,k)/dwallsq(i,j,k)
c
c   set some limits on tau
c
      if(tau.gt.1.2)   tau = 1.2
      if(tau.lt.0.001) tau = 0.001
c
      g       = tau + cw2*(tau*tau*tau*tau*tau*tau - tau)
      g6cw6   =  g*g*g*g*g*g/64.0
      fw      =  g*((1 + 0.015625)/(1 + g6cw6))**0.166667
c
      st3     =  fst3*fw*tvis*tvis/dwallsq(i,j,k)
c
c     set an extra source term which drives the turbulent viscosity towards
c     the mixing length value
c
      st_mixl = fstmix*(dyn_vis_ml(i,j,k)-visc_turb(i,j,k))/step(i,j,k)
c
c******************************************************************************
c  set the combined source term , tsource  .
c
      t_source(i,j,k) = rf_vis1*t_source(i,j,k) + 
     & rf_vis*(roavg*(st0 + st1 + st2 - st3 )*vol(i,j,k) + st_mixl)
c
c******************************************************************************
c  set the source term to zero for the dummy cell at the mixing plane.
c
      if(j.gt.2.and.indmix(j-1).eq.1) t_source(i,j,k) = 0.0
      if(j.gt.2) then
          if(indmix(j-2).eq.1)        t_source(i,j,k) = 0.0
      end if
c
 1600 continue
c
c
c    end of setting the turbulent viscosity
c****************************************************************************
c****************************************************************************
c     next  check for transition and reduce the source term in laminar regions
c
      iend1 = imid/2
      iendm = imm1 - iend1    
      kend1 = kmid/2
      kendm = kmm1 - kend1
c
      do 160 j = 2,jm
c
      nrw    = nrow(j)
      j1     = jstart(nrw)
      jrow   = j - j1 + 1
      jtrhub = jtran_k1(nrw)
      jtrtip = jtran_km(nrw)
      jtrlow = jtran_i1(nrw)
      jtrup  = jtran_im(nrw)
      jlerow = jle(nrw)
c
c   first on the hub
c
c   q3d
      if(km.eq.2)  go to 171
c   end q3d
      do 170 i = 1,imm1
      if(jtrhub.le.1.and.ftrans.lt.0.1) go to 38
      if(jrow.lt.jtrhub) go to 34
      rmax = 0.0
      do 37 k=1,kend1
      ratvis = visc_rat(i,j,k)
      if(ratvis.lt.rmax) go to 37
      rmax = ratvis
   37 continue
      if(rmax.gt.ftrans) go to 38
   34 continue
c
c      tvedge = visc_turb(i,kend1,j)
c      changed by ws/lx   17 may, 2017 
      tvedge = visc_turb(i,j,kend1)
c
      do 39 k = 1, kend1
      fk = float(kend1 + 1 -k)/kend1
      fk  = fk*fk*(3 - 2*fk)
      tvspec = fk*tvedge
   39 t_source(i,j,k) = t_source(i,j,k)
     &    -(visc_turb(i,j,k)-tvspec)/step(i,j,k)/ndamp
c
   38 continue

c  next on the casing

      if(jtrtip.le.1.and.ftrans.lt.0.1) go to 49
      if(jrow.lt.jtrtip) go to 46
      rmax = 0.0
      do 47 k = kendm, kmm1
      ratvis = visc_rat(i,j,k)
      if(ratvis.lt.rmax) go to 47
      rmax = ratvis
   47 continue
      if(rmax.gt.ftrans) go to 49
   46 continue
c
c      tvedge = visc_turb(i,kendm,j)
c      changed by ws/lx  17 may, 2017  
      tvedge = visc_turb(i,j,kendm)
c
      do 48 k = kendm, kmm1
      fk = float(k + 1 -kendm)/(km - kendm)
      fk  = fk*fk*(3 - 2*fk)
      tvspec  = fk*tvedge
   48 t_source(i,j,k) = t_source(i,j,k)
     &    -(visc_turb(i,j,k)-tvspec)/step(i,j,k)/ndamp
c
   49 continue
c
  170 continue
c
  171 continue
c
c    next on the lower - i=1, blade surface

      do 180  k=1,kmm1 
c
      if(jtrlow.le.1.and.ftrans.lt.0.1) go to 71 
      if(j.le.jlerow) go to 71
      if(jrow.lt.jtrlow.or.j.lt.jlerow) go to 70
      rmax = 0.0
      do 56 i=1,iend1
      ratvis = visc_rat(i,j,k)
      if(ratvis.lt.rmax) go to 56
      rmax = ratvis
   56 continue
      if(rmax.gt.ftrans)  go to 71
   70 continue
      tvedge = visc_turb(iend1,j,k)
      do 72 i = 1,iend1
      fi = float(iend1 + 1 -i)/iend1
      fi  = fi*fi*(3 - 2*fi)
      tvspec = fi*tvedge
   72 t_source(i,j,k) = t_source(i,j,k)
     &    -(visc_turb(i,j,k)-tvspec)/step(i,j,k)/ndamp
c
   71 continue

c    next on the upper, i = im, blade surface

      if(jtrup.le.1.and.ftrans.lt.0.1) go to 77 
      if(j.le.jlerow) go to 77
      if(jrow.lt.jtrup) go to 69
      rmax = 0.0
      do 73 i = iendm,imm1
      ratvis  = visc_rat(i,j,k)
      if(ratvis.lt.rmax) go to 73
      rmax = ratvis
   73 continue
      if(rmax.gt.ftrans) go to  77
   69 continue
      tvedge = visc_turb(iendm,j,k)
      do 74 i  = iendm,imm1
      fi = float(i + 1 -iendm)/(im - iendm)
      fi  = fi*fi*(3 - 2*fi)
      tvspec = fi*tvedge
   74 t_source(i,j,k) = t_source(i,j,k)
     &    -(visc_turb(i,j,k)-tvspec)/step(i,j,k)/ndamp

   77 continue
c
  180 continue
c
  160 continue
c
c********************************************************************************
c********************************************************************************
c
c     now the final viscosity is fixed set the stresses and heat flows in the element
c     jdd warning the div v term in the normal stresses is not yet included.
c
      do 190 k = 1,kmm1
      do 190 j = 2,jm
      do 190 i = 1,imm1
c
      vislam  = visc_lam(i,j,k)
      vistot  = visc_rat(i,j,k)*vislam
      thcond  = cp*vistot/prandtl
      vistot2 = 2.0*vistot
c     
      txr(i,j,k) = vistot*txr(i,j,k)
      txx(i,j,k) = vistot2*txx(i,j,k)
      txt(i,j,k) = vistot*txt(i,j,k)
      trx(i,j,k) = txr(i,j,k)
      trr(i,j,k) = vistot2*trr(i,j,k)
      trt(i,j,k) = vistot*trt(i,j,k)
      ttx(i,j,k) = txt(i,j,k)
      ttt(i,j,k) = vistot2*ttt(i,j,k)
      ttr(i,j,k) = trt(i,j,k)
      qxx(i,j,k) = -thcond*qxx(i,j,k)
      qrr(i,j,k) = -thcond*qrr(i,j,k)
      qtt(i,j,k) = -thcond*qtt(i,j,k)
c
  190 continue
c
c
c***************************************************************************************
c****************************************************************************************
c     calculate the wall shear stresses on all the solid surfaces.
c
      vislam16  = vislam*16.0
c
c******************first for blade surfaces.****************************
c******************************************************************************
c******************************************************************************
c    evaluate and smooth the pressure gradients if ypluswall is < -10.0
c
      if(ypluswall.lt.-10.0) call set_pwallgrad
c
c***************************************************************************************
c****************************************************************************************
      do 200 k=1,kmm1
      kp1 = k+1
      do 200 j=2,jm
c
c   skip if not on a blade surface
      if(ind(j).eq.0)  go to 220
c   
      jm1 = j-1
      nrw = nrow(j)
c   
c   also skip for cells above the tip gap
      if( (ktips(nrw).gt.0).and. 
     &    (k.ge.ktips(nrw)).and.(k.lt.ktipe(nrw)) ) 
     &    go to 220
c
c******************first for the i = 1  blade surface.*******************************
c                  calculate the wall shear stress.
c
      vxavg1 = vx(2,j,k)+vx(2,jm1,k)+vx(2,jm1,kp1)+vx(2,j,kp1)
      vravg1 = vr(2,j,k)+vr(2,jm1,k)+vr(2,jm1,kp1)+vr(2,j,kp1)
      wtavg1 = wt(2,j,k)+wt(2,jm1,k)+wt(2,jm1,kp1)+wt(2,j,kp1)
      roavg1 = ro(2,j,k)+ro(2,jm1,k)+ro(2,jm1,kp1)+ro(2,j,kp1)
      wavg   = sqrt(wtavg1*wtavg1 + vxavg1*vxavg1 + vravg1*vravg1)
      area1  = sqrt(abx(1,j,k)*abx(1,j,k) + abr(1,j,k)*abr(1,j,k)
     &              +abt(j,k)*abt(j,k)) 
      dperp  = vol(1,j,k)/area1
c
c**********************************************************************************
c**********************************************************************************
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then
c
         yplus_old =  yplus_i1(j,k)
         density   = 0.25*roavg1
         wslip     = 0.25*wavg
         call wallfun(1,j,k,1,dperp,dpds_cell(1,j,k),density,
     &                twalli1,yplus_old,wslip,yplus_new)
          yplus_i1(j,k) = amin1(1000.0,yplus_new)
c
      go to 365
c
      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************

      re     = dperp*roavg1*wavg/vislam16
      relog  = 1./alog(re)
c
c*******************************************
c    allow for roughness
c
      rough  = rough_l(nrw)
c
      if(rough.gt.1.0e-7) then
           if(rough.gt.dperp) rough = dperp
           rek   = re*rough/dperp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if
c
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c
c  check for fully rough surfaces
c
      if(rough.gt.1.0e-7) then
           twalli1    = cf*roavg1*wavg*wavg/128
c
c          vstar   = sqrt(twalli1/roavg1/4)  
c      changed by ws/lx  17 may, 2017 
           vstar   = sqrt(twalli1/(0.25*roavg1))
c
           plusk   = rough*vstar*roavg1/vislam/4
           relog10 = log10(dperp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
      end if
c*********************************************************************************
c     set the shear stress on the  i = 1 blade surface
c
      twalli1 = cf*roavg1*wavg*wavg/128
c
c*********************************************************************************
c
  365 continue
c
c*********************************************************************************
c    calculate yplus for use in setting  fyplus.
c
      vstar   = sqrt(twalli1/(0.25*roavg1))
      do i = 1,imid
c      ypls = ro(i,j,k)*vstar*sqrt(dwallsq(i,j,k))/vislam 
c      changed by ws/lx  17 may, 2017 
      yplsi1        = vstar*sqrt(dwallsq(i,j,k))*(0.25*roavg1)/vislam
      y_plus(i,j,k) = amin1(1000.0,yplsi1)
      end do
c        
      fmult   = -twalli1*area1/wavg
c
      xforce1 =  fmult*vxavg1
      rforce1 =  fmult*vravg1
      tforce1 =  fmult*wtavg1
      wvisc1  =  tforce1*wrad(j)*ravg_cell(j,k)
c
      forcex(1,j,k)    = forcex(1,j,k)    + xforce1
      forcer(1,j,k)    = forcer(1,j,k)    + rforce1
      forcet(1,j,k)    = forcet(1,j,k)    + tforce1
      esource(1,j,k)   = esource(1,j,k)   + wvisc1
c
c***********************************************************************************
c***************************now for the i = im blade surface************************
c                     calculate the wall shear stress.
c
      vxavgim = vx(imm1,j,k)+vx(imm1,jm1,k)
     &        + vx(imm1,jm1,kp1)+vx(imm1,j,kp1)
      vravgim = vr(imm1,j,k)+vr(imm1,jm1,k)
     &        + vr(imm1,jm1,kp1)+vr(imm1,j,kp1)
      wtavgim = wt(imm1,j,k)+wt(imm1,jm1,k)
     &        + wt(imm1,jm1,kp1)+wt(imm1,j,kp1)
      roavgim = ro(imm1,j,k)+ro(imm1,jm1,k)
     &        + ro(imm1,jm1,kp1)+ro(imm1,j,kp1)
      wavg    = sqrt(wtavgim*wtavgim+vxavgim*vxavgim+vravgim*vravgim)
      areaim  = sqrt(abx(im,j,k)*abx(im,j,k) + abr(im,j,k)*abr(im,j,k)
     &              +abt(j,k)*abt(j,k)) 
      dperp   = vol(imm1,j,k)/areaim
c******************************************************************************
c******************************************************************************
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then

           yplus_old =  yplus_im(j,k)
           density   = 0.25*roavgim
           wslip     = 0.25*wavg
           call wallfun(imm1,j,k,im,dperp,dpds_cell(imm1,j,k),
     &             density,twallim,yplus_old,wslip,yplus_new)
           yplus_im(j,k) = amin1(1000.0,yplus_new)
c
      go to 375
c
      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************
      re      = dperp*roavgim*wavg/vislam16
      relog   = 1./alog(re)
c
c*******************************************
c    allow for roughness
c
      rough  = rough_u(nrw)
c
      if(rough.gt.1.0e-7) then
           if(rough.gt.dperp) rough = dperp
           rek   = re*rough/dperp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if
c
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c
c  check for fully rough surfaces
c
      if(rough.gt.1.0e-7) then
           twallim     = cf*roavgim*wavg*wavg/128
c
c          vstar   = sqrt(twallim/roavgim/4)  
c      changed by ws/lx  17 may, 2017  
           vstar   = sqrt(twallim/(0.25*roavgim))
c
           plusk   = rough*vstar*roavgim/vislam/4
           relog10 = log10(dperp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
      end if
c
c*********************************************************************************
c    set the shear stress on the  i = im blade surface
c
      twallim  = cf*roavgim*wavg*wavg/128
c
c*********************************************************************************
c
  375 continue
c
c*********************************************************************************
c     calculate yplus for use in setting  fyplus.
c
      vstar   = sqrt(twallim/(0.25*roavgim))
      istrt = imid + 1
      do i = istrt,im
c      ypls = ro(i,j,k)*vstar*sqrt(dwallsq(i,j,k))/vislam
c      changed by ws/lx  17 may, 2017  
      yplsim        = vstar*sqrt(dwallsq(i,j,k))*(0.25*roavgim)/vislam
      y_plus(i,j,k) = amin1(1000.,yplsim)
      end do
c
      fmult    = -twallim*areaim/wavg
c 
      xforceim =  fmult*vxavgim
      rforceim =  fmult*vravgim
      tforceim =  fmult*wtavgim
      wviscim  =  tforceim*wrad(j)*ravg_cell(j,k) 
c
      forcex(imm1,j,k)  = forcex(imm1,j,k)   + xforceim
      forcer(imm1,j,k)  = forcer(imm1,j,k)   + rforceim
      forcet(imm1,j,k)  = forcet(imm1,j,k)   + tforceim
      esource(imm1,j,k) = esource(imm1,j,k)  + wviscim
c
  220 continue
c
c*****************************************************************************
c
c     now set the forces on the periodic cells, i=1 and i= im-1,upstream of the 
c     le and downstream of the te and also for cells in the tip gap.
c     do not treat cells on the solid blade surfaces as periodic.
c*********************************************************************************
c
       if( (ind(j).eq.1).and.(k.lt.ktips(nrw)).or.(k.gt.ktipe(nrw)) )
     &      go to 225
c
      i   = 1
      im1 = imm1
      xfor = 0.5*((txr(i,j,k)+txr(im1,j,k))*abr(i,j,k)
     &    +       (txt(i,j,k)+txt(im1,j,k))*abt(j,k)
     &    +       (txx(i,j,k)+txx(im1,j,k))*abx(i,j,k))
      forcex(im1,j,k) = forcex(im1,j,k) + xfor
      forcex(i,j,k)   = forcex(i,j,k)   - xfor
      avg = 0.5*(forcex(i,j,k) + forcex(im1,j,k))
      forcex(i,j,k)   = avg
      forcex(im1,j,k) = avg

      rfor = 0.5*((trr(i,j,k)+trr(im1,j,k))*abr(i,j,k)
     &    +       (trt(i,j,k)+trt(im1,j,k))*abt(j,k)
     &    +       (trx(i,j,k)+trx(im1,j,k))*abx(i,j,k))
      forcer(im1,j,k)  = forcer(im1,j,k) + rfor
      forcer(i,j,k)    = forcer(i,j,k)   - rfor
      avg = 0.5*(forcer(i,j,k) + forcer(im1,j,k))
      forcer(i,j,k)   = avg
      forcer(im1,j,k) = avg

      tfor = 0.5*((ttr(i,j,k)+ttr(im1,j,k))*abr(i,j,k)
     &    +       (ttt(i,j,k)+ttt(im1,j,k))*abt(j,k)
     &    +       (ttx(i,j,k)+ttx(im1,j,k))*abx(i,j,k))
      forcet(im1,j,k)  = forcet(im1,j,k) + tfor
      forcet(i,j,k)    = forcet(i,j,k)   - tfor
      avg = 0.5*(forcet(i,j,k) + forcet(im1,j,k))
      forcet(i,j,k)   = avg
      forcet(im1,j,k) = avg
c
c     set the heat flow and viscous work
c
      qflow =  0.5*((qxx(i,j,k) + qxx(im1,j,k))*abx(i,j,k)
     &      +       (qrr(i,j,k) + qrr(im1,j,k))*abr(i,j,k)
     &      +       (qtt(i,j,k) + qtt(im1,j,k))*abt(j,k))
c
      vxavg = vx(i,j,k)+vx(i,jm1,k)+vx(i,jm1,kp1)+vx(i,j,kp1)
      vravg = vr(i,j,k)+vr(i,jm1,k)+vr(i,jm1,kp1)+vr(i,j,kp1)
      vtavg = vt(i,j,k)+vt(i,jm1,k)+vt(i,jm1,kp1)+vt(i,j,kp1)
c
      wvisc = -0.25*(xfor*vxavg + rfor*vravg + tfor*vtavg)
c
      esource(i,j,k)  = esource(i,j,k)   + qflow + wvisc
      esource(im1,j,k)= esource(im1,j,k) - qflow - wvisc 
      avg = 0.5*(esource(i,j,k) + esource(im1,j,k))
      esource(i,j,k)   = avg
      esource(im1,j,k) = avg
c
  225 continue
c
c*********************************************************************************
c*********************************************************************************
c     now sum the stresses on the other i = constant faces 
c
      do 250 i=2,imm1
      im1  = i-1
      xfor = 0.5*((txr(i,j,k)+txr(im1,j,k))*abr(i,j,k)
     &    +       (txt(i,j,k)+txt(im1,j,k))*abt(j,k)
     &    +       (txx(i,j,k)+txx(im1,j,k))*abx(i,j,k))
      forcex(im1,j,k) = forcex(im1,j,k) + xfor
      forcex(i,j,k)   = forcex(i,j,k)   - xfor

      rfor = 0.5*((trr(i,j,k)+trr(im1,j,k))*abr(i,j,k)
     &    +       (trt(i,j,k)+trt(im1,j,k))*abt(j,k)
     &    +       (trx(i,j,k)+trx(im1,j,k))*abx(i,j,k))
      forcer(im1,j,k)  = forcer(im1,j,k) + rfor
      forcer(i,j,k)    = forcer(i,j,k)   - rfor

      tfor = 0.5*((ttr(i,j,k)+ttr(im1,j,k))*abr(i,j,k)
     &    +       (ttt(i,j,k)+ttt(im1,j,k))*abt(j,k)
     &    +       (ttx(i,j,k)+ttx(im1,j,k))*abx(i,j,k))
      forcet(im1,j,k)  = forcet(im1,j,k) + tfor
      forcet(i,j,k)    = forcet(i,j,k)   - tfor
c
c     set the heat flow and viscous work
c
      qflow =  0.5*((qxx(i,j,k) + qxx(im1,j,k))*abx(i,j,k)
     &      +       (qrr(i,j,k) + qrr(im1,j,k))*abr(i,j,k)
     &      +       (qtt(i,j,k) + qtt(im1,j,k))*abt(j,k))
c
      vxavg = vx(i,j,k)+vx(i,jm1,k)+vx(i,jm1,kp1)+vx(i,j,kp1)
      vravg = vr(i,j,k)+vr(i,jm1,k)+vr(i,jm1,kp1)+vr(i,j,kp1)
      vtavg = vt(i,j,k)+vt(i,jm1,k)+vt(i,jm1,kp1)+vt(i,j,kp1)
c
      wvisc = -0.25*(xfor*vxavg + rfor*vravg + tfor*vtavg)
c
      esource(i,j,k)  = esource(i,j,k)   + qflow + wvisc
      esource(im1,j,k)= esource(im1,j,k) - qflow - wvisc 

  250 continue

  200 continue
c
c**************************************************************************************
c***************************************************************************************
c**********calculate the wall shear stresses on the hub and casing surfaces**********
c
c   q3d
      if(km.eq.2) go to 401
c   end q3d
      do 400 j=2,jm
      jm1  = j-1
      do 400 i=1,imm1
      ip1 = i+1
c
c     first calculate the wall shear stress.for the k = 1,  hub, wall.   
c
      vxavg1 = vx(i,j,2)+vx(ip1,j,2)+vx(ip1,jm1,2)+vx(i,jm1,2)
      vravg1 = vr(i,j,2)+vr(ip1,j,2)+vr(ip1,jm1,2)+vr(i,jm1,2)
      vtavg1 = vt(i,j,2)+vt(ip1,j,2)+vt(ip1,jm1,2)+vt(i,jm1,2)
      roavg1 = ro(i,j,2)+ro(ip1,j,2)+ro(ip1,jm1,2)+ro(i,jm1,2)
      ravg1  = r(j,2)   + r(j,2)    + r(jm1,2)    + r(jm1,2)
      wtavg1 = vtavg1   - whub(j)*ravg1
      wavg   = sqrt(vxavg1*vxavg1 + vravg1*vravg1 + wtavg1*wtavg1)
      area1  = sqrt(asx(i,j,1)*asx(i,j,1) + asr(i,j,1)*asr(i,j,1))
      dperp  = vol(i,j,1)/area1
c
c**********************************************************************************
c**********************************************************************************
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then

           yplus_old =  yplus_k1(i,j)
           density   = 0.25*roavg1
           wslip     = 0.25*wavg
         call wallfun(i,j,1,1,dperp,dpds_cell(i,j,1),density,
     &                twallk1,yplus_old,wslip,yplus_new)
           yplus_k1(i,j) = amin1(1000.0,yplus_new)
c
      go to 345

      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************
      re     = dperp*roavg1*wavg/vislam16
      relog  = 1./alog(re)
c
c**************************************************************************************
c
      rough  = rough_h(nrw)
c
      if(rough.gt.1.0e-7) then
           if(rough.gt.dperp) rough = dperp
           rek   = re*rough/dperp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if
c
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c
c  check for fully rough surfaces
c
      if(rough.gt.1.0e-7) then
           twallk1     = cf*roavg1*wavg*wavg/128
c
c		   vstar   = sqrt(twallk1/roavg1/4)  
c      changed by ws/lx  17 may, 2017  
           vstar   = sqrt(twallk1/(0.25*roavg1))
c
           plusk   = rough*vstar*roavg1/vislam/4
           relog10 = log10(dperp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
c
      end if
c**************************************************************************************
c    set the shear stress on the hub surface, k=1.
c
      twallk1    = cf*roavg1*wavg*wavg/128
c
c**************************************************************************************
c
  345 continue
c
c**************************************************************************************
c     calculate yplus for use in setting  fyplus
      if(ibound.eq.0.or.ibound.eq.2) then
      do k = 1,kmid
      vstar         = sqrt(twallk1/(0.25*roavg1))
      yplsk1        = vstar*sqrt(dwallsq(i,j,k))*(0.25*roavg1)/vislam
      yplsp         = y_plus(i,j,k)
      y_plus(i,j,k) = amin1(yplsk1,yplsp)
      end do
      end if
c
      fmult   = -twallk1*area1/wavg
      if(ibound.eq.1.or.ibound.gt.2) fmult = 0.0
c
      xforce1 = fmult*vxavg1
      rforce1 = fmult*vravg1
      tforce1 = fmult*wtavg1
      wvisc1  = tforce1*0.25*whub(j)*ravg1
c
      forcex(i,j,1)    = forcex(i,j,1)    + xforce1
      forcer(i,j,1)    = forcer(i,j,1)    + rforce1
      forcet(i,j,1)    = forcet(i,j,1)    + tforce1
      esource(i,j,1)   = esource(i,j,1)   + wvisc1
c
c**************************************************************************************
c**************now set the shear stress on  the casing, k = km ************************
c
      vxavgkm = vx(i,j,kmm1)+vx(ip1,j,kmm1)
     &        + vx(ip1,jm1,kmm1)+vx(i,jm1,kmm1)
      vravgkm = vr(i,j,kmm1)+vr(ip1,j,kmm1)
     &        + vr(ip1,jm1,kmm1)+vr(i,jm1,kmm1)
      vtavgkm = vt(i,j,kmm1)+vt(ip1,j,kmm1)
     &        + vt(ip1,jm1,kmm1)+vt(i,jm1,kmm1)
      roavgkm = ro(i,j,kmm1)+ro(ip1,j,kmm1)
     &        + ro(ip1,jm1,kmm1)+ro(i,jm1,kmm1)
      ravgkm  = r(j,kmm1) + r(j,kmm1) + r(jm1,kmm1) + r(jm1,kmm1)
      wtavgkm = vtavgkm   - wtip(j)*ravgkm
      wavg    = sqrt(vxavgkm*vxavgkm+vravgkm*vravgkm+wtavgkm*wtavgkm)
      areakm  = sqrt(asx(i,j,km)*asx(i,j,km) + asr(i,j,km)*asr(i,j,km))
      dperp   = vol(i,j,kmm1)/areakm
c
c************************************************************************
c****************************************************************************
c    use the shih et al  wall functions if ypluswall is negative.
      if(ypluswall.lt.-0.001) then

           yplus_old  =  yplus_km(i,j)
           density    = 0.25*roavgkm
           wslip      = 0.25*wavg
        call wallfun(i,j,kmm1,km,dperp,dpds_cell(i,j,kmm1),density,
     &              twallkm,yplus_old,wslip,yplus_new)
           yplus_km(i,j) = amin1(1000.0,yplus_new)
c
      go to 355
c
      end if
c    end of shih et al wallfunctions
c**********************************************************************************
c**********************************************************************************
c
      re      = dperp*roavgkm*wavg/vislam16
      relog   = 1./alog(re)
c
c**************************************************************************************
c    allow for roughness
c
      rough  = rough_t(nrw)
c
      if(rough.gt.1.0e-7) then
           if(rough.gt.dperp) rough = dperp
           rek   = re*rough/dperp
           rek   = rek - 80.0
           if(rek.lt.0.0)  rek = 0.0
           reksq = rek*rek
           a1 = -.00178493 + .0000814923*rek + .000000150445*reksq
           a2 = .029072 - .001584*rek - .00000225194*reksq
           a3=.270313+.0091409*rek+.00000451537*reksq +
     &     .00000000464767*reksq*rek
           cf    = a1 + a2*relog + a3*relog*relog
      else
c
c    end roughness, next eqn for smooth surfaces.
c
           cf = -0.00178493 + 0.029072*relog + 0.270313*relog*relog
c
      end if
c
c   take cf as the max of the laminar and turbulent values.
c
      cflam = 2.0/re
      cf    = amax1(cf,cflam)
      if(re.lt.125.) cf = cflam
c
c  check for fully rough surfaces
c
      if(rough.gt.1.0e-7) then
           twallkm    = cf*roavgkm*wavg*wavg/128
c
c		   vstar   = sqrt(twallkm/roavgkm/4)  
c      changed by ws/lx  17 may, 2017  
           vstar   = sqrt(twallkm/(0.25*roavgkm))
c
           plusk   = rough*vstar*roavgkm/vislam/4
           relog10 = log10(dperp/rough)
           funcn   = 5.75*relog10 + 8.5
           cf_full_rough = 2.0/(funcn*funcn)
           if(plusk.gt.45.and.cf.gt.cf_full_rough) cf = cf_full_rough
      end if
c
c**************************************************************************************
c   set the shear stress on the casing , k = km, surface
c
      twallkm    = cf*roavgkm*wavg*wavg/128
c
c**************************************************************************************
c
  355 continue
c
c**************************************************************************************
c
c     calculate yplus for use in setting fyplus.
      if(ibound.eq.0.or.ibound.eq.1) then
      do k = kmid+1,kmm1
      vstar         = sqrt(twallkm/(0.25*roavgkm))
      yplskm        = vstar*sqrt(dwallsq(i,j,k))*(0.25*roavgkm)/vislam
      yplsp         = y_plus(i,j,k)
      y_plus(i,j,k) = amin1(yplskm,yplsp)
      end do
      end if
c
      fmult   = -twallkm*areakm/wavg
      if(ibound.ge.2) fmult = 0.0 
c
      xforcekm = fmult*vxavgkm
      rforcekm = fmult*vravgkm
      tforcekm = fmult*wtavgkm
      wvisckm  =  tforcekm*0.25*wtip(j)*ravgkm
c
      forcex(i,j,kmm1) = forcex(i,j,kmm1) + xforcekm
      forcer(i,j,kmm1) = forcer(i,j,kmm1) + rforcekm
      forcet(i,j,kmm1) = forcet(i,j,kmm1) + tforcekm
      esource(i,j,kmm1)= esource(i,j,kmm1)+ wvisckm
c
c     now sum the stresses on the other k = constant faces
c
      do 450 k=2,kmm1
      km1  = k-1
      xfor = 0.5*((txr(i,j,k)+txr(i,j,km1))*asr(i,j,k)
     &    +       (txx(i,j,k)+txx(i,j,km1))*asx(i,j,k))
      forcex(i,j,km1) = forcex(i,j,km1) + xfor
      forcex(i,j,k)   = forcex(i,j,k)   - xfor

      rfor = 0.5*((trr(i,j,k)+trr(i,j,km1))*asr(i,j,k)
     &    +       (trx(i,j,k)+trx(i,j,km1))*asx(i,j,k))
      forcer(i,j,km1)  = forcer(i,j,km1) + rfor
      forcer(i,j,k)    = forcer(i,j,k)   - rfor

      tfor = 0.5*((ttr(i,j,k)+ttr(i,j,km1))*asr(i,j,k)
     &    +       (ttx(i,j,k)+ttx(i,j,km1))*asx(i,j,k))
      forcet(i,j,km1)  = forcet(i,j,km1) + tfor
      forcet(i,j,k)    = forcet(i,j,k)   - tfor
c
c     set the heat flow and viscous work
c
      qflow =  0.5*((qxx(i,j,k) + qxx(i,j,km1))*asx(i,j,k)
     &      +       (qrr(i,j,k) + qrr(i,j,km1))*asr(i,j,k))
c
      vxavg = vx(i,j,k)+vx(ip1,j,k)+vx(ip1,jm1,k)+vx(i,jm1,k)
      vravg = vr(i,j,k)+vr(ip1,j,k)+vr(ip1,jm1,k)+vr(i,jm1,k)
      vtavg = vt(i,j,k)+vt(ip1,j,k)+vt(ip1,jm1,k)+vt(i,jm1,k)

      wvisc = -0.25*(xfor*vxavg + rfor*vravg + tfor*vtavg)

      esource(i,j,k)  = esource(i,j,k)   + qflow + wvisc
      esource(i,j,km1)= esource(i,j,km1) - qflow - wvisc 
c
  450 continue
c
  400 continue
c  re enter here if q3d
  401 continue
c   end q3d 
c**************************************************************************************
c**************************************************************************************
c     all the viscous forces are now set. use them to update the global body force terms.
c     note the negative sign because the source term is subtracted in subroutine tstep.
c 
      do 500 k=1,kmm1
      do 500 j=2,jm
      do 500 i=1,imm1
      xforce(i,j,k)   = rf_vis1*xforce(i,j,k)   - rf_vis*forcex(i,j,k)
      tforce(i,j,k)   = rf_vis1*tforce(i,j,k)   - rf_vis*forcet(i,j,k)
     &                                           *ravg_cell(j,k)
      rforce(i,j,k)   = rf_vis1*rforce(i,j,k)   - rf_vis*forcer(i,j,k)
      qsource(i,j,k)  = rf_vis1*qsource(i,j,k)  - rf_vis*esource(i,j,k)
  500 continue
c
c**************************************************************************************
c**************************************************************************************
c
      return
      end
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c    
      subroutine smooth_resid(var,smfac,nsmth)
c
c******************************************************************************************
c     this subroutine applies a simple linear smoothing to the
c     variable var(i,j,k). a zero gradient condition is applied at the boundaries.
c     with no smoothing across the mixing plane.
c******************************************************************************************
c
      include 'commall-open-19.2'
c
      dimension tempvar(jd),var(id,jd,kd)
c
      sfac1 = 1.0 -smfac
      sfach = 0.5*smfac
c 
      do 200, nsm=1,nsmth
c
c     smooth in the j direction
c
      do 100 nr = 1,nrows
      j1   = jstart(nr)
      j1p1 = j1+1
      j1p2 = j1+2
      j2   = jmix(nr)
      j2m1 = j2-1
c
      do 40 k=1,kmm1
      do 40 i=1,imm1
      do 50 j=j1p2,j2m1
      tempvar(j) = sfac1*var(i,j,k)
     &           + sfach*(var(i,j+1,k)+var(i,j-1,k))
   50 continue
      tempvar(j1p1) = var(i,j1p1,k)+ sfach*(var(i,j1p2,k)-var(i,j1p1,k))
      tempvar(j2)   = var(i,j2,k)  + sfach*(var(i,j2m1,k)-var(i,j2,k))
      do 60 j = j1p1,j2
      var(i,j,k) = tempvar(j)
   60 continue
   40 continue
c
c     smooth in the i direction
c     tflow addition
      if(im.eq.2) go to 11

      do 10 k=1,kmm1
      do 10 j=j1p1,j2
      do 20 i=2,imm2
      tempvar(i) = sfac1*var(i,j,k)
     &           + sfach*(var(i+1,j,k)+var(i-1,j,k))
   20 continue
      tempvar(1)   = var(1,j,k)    + sfach*(var(2,j,k)   -var(1,j,k))
      tempvar(imm1)= var(imm1,j,k) + sfach*(var(imm2,j,k)-var(imm1,j,k))
      do 30 i=1,imm1
      var(i,j,k) = tempvar(i)
   30 continue
   10 continue

   11 continue
c
c  q3d
      if(km.eq.2) go to 100
c  end q3d
c
c     smooth in the k direction
c
      do 70 i=1,imm1
      do 70 j=j1p1,j2
      do 80 k=2,kmm2
      tempvar(k) = sfac1*var(i,j,k) 
     &           + sfach*(var(i,j,k+1)+var(i,j,k-1))
   80 continue
      tempvar(1)   = var(i,j,1)    + sfach*(var(i,j,2)   -var(i,j,1))
      tempvar(kmm1)= var(i,j,kmm1) + sfach*(var(i,j,kmm2)-var(i,j,kmm1))
      do 90 k=1,kmm1
      var(i,j,k) = tempvar(k)
   90 continue
   70 continue
c
  100 continue
c
  200 continue
c
c******************************************************************************************
c******************************************************************************************
      return
      end
c
c*****************************************************************************
c*****************************************************************************
c
      subroutine new_mixplan
c
c******************************************************************************************
c     this subroutine relaxes the properties downstream of the mixing plane towards pitchwise
c     uniform isentropic values.and also allows the flow angles to be extrapolated from downstream.
c     the flow can cross the mixing plane in either direction.
c******************************************************************************************
c
      include  'commall-open-19.2'
c
      dimension  alphaa(id),gama(id),wref_jp2(id),roref_jp2(id),
     &           hstat(id),pisent(id),roisent(id),
     &           wmach_jp2(id),tisent(id),wis_jp2(id),ro_sub_jp2(id)
c
      rfmix1   = 1.0 - rfmix
c
      do 500 nr = 1,nrwsm1  
c
      do 400  k = 1,km
c
      ksub = k
      if(k.eq.km) ksub = km-1
c
c   check the flow direction at the mixing plane.
c
           j     = jmix(nr)
      if(flowx(imid,j,ksub).lt.0.0) then 
           jp1   = j+1
           jp2   = j+2  
           jp3   = j+3 
           jp4   = j+4
           f_angle  = fangle
      else
           j   = jmix(nr) + 1
           jp1 = j-1
           jp2 = j-2
           jp3 = j-3
           jp4 = j-4
           f_angle  = 0.0
      end if
c
           f_angle1 = 1.0 - f_angle
c
c*************************************************************************************************
c   calculate the mixed out flow condition at the mixing plane. the flow there is pitchwise uniform
c   so we only need to use the mid-pitch values.
c
      rjp1     = r(jp1,k)
      ujp1     = ublade(jp1,k)
      romix    = ro(imid,jp1,k)
      vxmix    = rovx(imid,jp1,k)/romix
      vrmix    = rovr(imid,jp1,k)/romix
      vtmix    = rorvt(imid,jp1,k)/romix/rjp1
      roemix   = roe(imid,jp1,k)
      homix    = ho(imid,jp1,k)   
      vmsq     = vxmix*vxmix + vrmix*vrmix
      vmmix    = sqrt(vmsq)
      vabsq    = vmsq + vtmix*vtmix
      wtmix    = vtmix - ujp1
      wmixsq   = vmsq  + wtmix*wtmix
      wmix     = sqrt(wmixsq)
      alphamix = asin(wtmix/wmix)
      gamamix  = atan2(vrmix,vxmix)

c
      if(ifgas.eq.0) then
           tmix      = (roemix/romix - 0.5*vabsq)/cv
c
           if(itimst.eq.5.or.itimst.eq.6) 
     &      tmix = (roemix/rosub(imid,jp1,k) - 0.5*vabsq)/cv
c
           rothalpy  = cp*tmix + 0.5*(wmixsq - ujp1*ujp1)
           gagas     = ga
      else
            emix     = roemix/romix - 0.5*vabsq
            tmix     = tfrome(emix,tref,eref,et1,et2,et3,et4)
            hstatic  = hfromt(tmix,tref,href,cp1,cp2,cp3)
            rothalpy = hstatic + 0.5*(wmixsq - ujp1*ujp1)
            delt     =  tmix - tref
            cpgas    = cp1 + cp2*delt + cp3*delt*delt
            gagas    = cpgas/(cpgas-rgas)
      end if
c
      pmix    = romix*rgas*tmix
c
      if(itimst.eq.5.or.itimst.eq.6) 
     &   pmix = po_ref + dp_dro*(rosub(imid,jp1,k) - ro_ref)
c
      if(itimst.eq.6) porelmix = pmix + 0.5*densty*wmixsq
c
c
c*********************************************************************************************
c
      nmach = 0
      do 350 i = 1,im
c
c     calculate the velocities at jmix + 2
c
      rjp2       = r(jp2,k)
      unowjp2    = ublade(jp2,k)
      rojp2      = ro(i,jp2,k)
      xveljp2    = rovx(i,jp2,k)/rojp2
      rveljp2    = rovr(i,jp2,k)/rojp2
      tveljp2    = rorvt(i,jp2,k)/rjp2/rojp2
      wtjp2      = tveljp2 - unowjp2
      vmsqjp2    = xveljp2*xveljp2 + rveljp2*rveljp2
      velsqjp2   = vmsqjp2 + tveljp2*tveljp2
      wact_sq_jp2= vmsqjp2 + wtjp2*wtjp2
      wact_jp2   = sqrt(wact_sq_jp2)
      ro_act_jp2 = rojp2
      t_act_jp2  = (roe(i,jp2,k)/rojp2 - 0.5*velsqjp2)/cv
c
      if(itimst.eq.5.or.itimst.eq.6) 
     &       t_act_jp2 = (roe(i,jp2,k)/rosub(i,jp2,k) - 0.5*velsqjp2)/cv
c
c
c**********************************************************************************************
c    obtain velocities and angles at jmix + 3. only the angles are used.
c
      rjp3       = r(jp3,k)
      unowjp3    = ublade(jp3,k)
      rojp3      = ro(i,jp3,k)
      rovxjp3    = rovx(i,jp3,k)
      rovrjp3    = rovr(i,jp3,k)
      rovtjp3    = rorvt(i,jp3,k)/rjp3
      rowtjp3    = rovtjp3 - unowjp3*rojp3
      rovmsqjp3  = rovxjp3*rovxjp3 + rovrjp3*rovrjp3
      rowrelsqjp3= rovmsqjp3 + rowtjp3*rowtjp3
      rowreljp3  = sqrt(rowrelsqjp3)
      alphajp3   = asin(rowtjp3/rowreljp3)
      gamajp3    = atan2(rovrjp3,rovxjp3)
c
c*********************************************************************************************
c    obtain velocities and angles at jmix + 4. only the angles are used.
c
      rjp4       = r(jp4,k)
      unowjp4    = ublade(jp4,k)
      rojp4      = ro(i,jp4,k)
      rovxjp4    = rovx(i,jp4,k)
      rovrjp4    = rovr(i,jp4,k)
      rovtjp4    = rorvt(i,jp4,k)/rjp4
      rowtjp4    = rovtjp4 - unowjp4*rojp4
      rovmsqjp4    = rovxjp4*rovxjp4 + rovrjp4*rovrjp4
      rowrelsqjp4  = rovmsqjp4 + rowtjp4*rowtjp4
      rowreljp4    = sqrt(rowrelsqjp4)
      alphajp4     = asin(rowtjp4/rowreljp4)
      gamajp4      = atan2(rovrjp4,rovxjp4)
c
c*********************************************************************************************
c
c   next the isentropic values at jmix+2 assuming an isentropic expansion from the mixing plane conditions
c   to  pstat  at jmix+2 .
c
      if(ifgas.eq.0) then
           pstat      = rojp2*rgas*t_act_jp2

           if(itimst.eq.5.or.itimst.eq.6)
     &     pstat   = po_ref  + dp_dro*(rosub(i,jp2,k) - ro_ref)

           tisent_jp2 = tmix*(pstat/pmix)**fga
           wis_sq_jp2 = 2*(rothalpy- cp*tisent_jp2+ 0.5*unowjp2*unowjp2)
c
           if(itimst.eq.6) wis_sq_jp2 = 2.0*(porelmix - pstat)/densty
c
      else
           estat  = roe(i,jp2,k)/rojp2 - 0.5*velsqjp2
           tstat  = tfrome(estat,tref,eref,et1,et2,et3,et4)
           t_act_jp2 = tstat
           pstat  = rojp2*rgas*tstat
           prat   = pstat/pmix
           trat   = trat_from_prat
     &     (tmix,prat,fgagas,r_alpha,balpha1,balpha2)
           tisent_jp2 = tmix*trat
           hstat(i)  = hfromt(tisent_jp2,tref,href,cp1,cp2,cp3)
           wis_sq_jp2  = 2*(rothalpy - hstat(i) + 0.5*unowjp2*unowjp2)
      end if
c
      roisent_jp2  = pstat/rgas/tisent_jp2
c
      if(itimst.eq.6) roisent_jp2 = densty
c
      if(itimst.eq.5.or.itimst.eq.6) 
     &       ro_sub_jp2(i)= ro_ref + (pstat - po_ref)/dp_dro

      roisent(i)  = roisent_jp2
      pisent(i)   = pstat
      tisent(i)   = tisent_jp2
c
      if(wis_sq_jp2.lt.1.0) wis_sq_jp2 = 1.0
      wis_jp2(i)     = sqrt(wis_sq_jp2)
c
c**********************************************************************************************
c     set the current values at jmix+2 for use in finding the mach number and angle variation.
c
      wref_jp2(i)    =  wact_jp2 
      roref_jp2(i)   =  ro_act_jp2
      wmach_jp2(i)   =  wact_jp2/sqrt(gagas*rgas*t_act_jp2)
      if(wmach_jp2(i).gt.1.0) nmach = nmach + 1
c
c***********************************************************************************************
c    set the flow angles at jmix+2 weighting the flow direction according to f_angle
c
c      alphaa(i) = f_angle1*alphamix  + f_angle*(2*alphajp3 - alphajp4)
c      gama(i)   = f_angle1*gamamix   + f_angle*(2*gamajp3  - gamajp4)
      alphaa(i) = f_angle1*alphamix  + f_angle*alphajp3
      gama(i)   = f_angle1*gamamix   + f_angle*gamajp3
c      alphaa(i) = f_angle1*alphamix  + f_angle*(alphajp3 +alphajp4)/2
c      gama(i)   = gamamix
c
  350 continue
c
c*********************************************************************************
c  
      if(itimst.ge.5) go to 250
      if(nmach.eq.0)  go to 250
c
c********************************************************************************
c     apply characteristics relationship between the yaw angle and mach no. if
c     the relative flow is supersonic
c
c**********************************************************************************************
c    set frotn to allow for left or right running mach waves
c      
      frotn  =  1.0
      if(wtmix.lt.0.0) frotn = -1.0
c
c   if the relative flow is supersonic,calculate the pitchwise variation in angle needed
c   to satisfy the prandtl-meyer relationship.
c
      do 211 i=imid+1,im 
      im1 = i-1  
      if(wmach_jp2(i).gt.1.001.and.wmach_jp2(im1).gt.1.001) then      
      delp     = pisent(i) - pisent(im1)   
      deltai   = sqrt(wmach_jp2(i)*wmach_jp2(i) - 1 )/
     &           (roref_jp2(i)*wref_jp2(i)*wref_jp2(i))     
      deltaim1 = sqrt(wmach_jp2(im1)*wmach_jp2(im1) -1 )/
     &           (roref_jp2(im1)*wref_jp2(im1)*wref_jp2(im1))
      alphaa(i) = alphaa(im1) + frotn*delp*0.5*(deltai + deltaim1)
      end if       
  211 continue  
c
      do 212 i=imid-1,1,-1   
      ip1  = i + 1
      if(wmach_jp2(i).gt.1.001.and.wmach_jp2(ip1).gt.1.001) then      
      delp     = pisent(i) - pisent(ip1)   
      deltai   = sqrt(wmach_jp2(i)*wmach_jp2(i) - 1 )/
     &           (roref_jp2(i)*wref_jp2(i)*wref_jp2(i))     
      deltaip1 = sqrt(wmach_jp2(ip1)*wmach_jp2(ip1) -1 )/
     &           (roref_jp2(ip1)*wref_jp2(ip1)*wref_jp2(ip1))
      alphaa(i) = alphaa(ip1) + frotn*delp*0.5*(deltai + deltaip1)
      end if       
  212 continue        
c
c    enforce periodicity of angle
c  
      rotn  = alphaa(im) - alphaa(1)
      sumfp = 0.0
      do 230 i=1,im
      alphaa(i) = alphaa(i) - rotn*(sumfp - 0.5)
      sumfp = sumfp + fp(i)
  230 continue
c
c******************************************************************************
  250 continue
c******************************************************************************
c 
c    store the value of pitchwise flow angle on the first time step
c
c      if(nstep.eq.1) then
c      do 240 i=1,im
c      alpha_s(i,k,nr) = alphaa(i)
c  240 continue
c      end if
c****************************************************************************************
c****************************************************************************************
c
      do 300 i = 1,im
c
c     relax changes to the value of alphaa by rfmix
c
c      alphaa(i)   = rfmix1*alpha_s(i,k,nr) + rfmix*alphaa(i)
c      alpha_s(i,k,nr) = alphaa(i)
c
c     set the isentropic velocities at jmix+2. based on the calculated relative velocity magnitude
c     and the relative flow direction.
c
      ro_jp2    = roisent(i)      
      wt_jp2    = wis_jp2(i)*sin(alphaa(i))
      vm_jp2    = wis_jp2(i)*cos(alphaa(i))
      vx_jp2    = vm_jp2*cos(gama(i))
      vr_jp2    = vm_jp2*sin(gama(i))
      vm_jp2sq  = vx_jp2*vx_jp2 + vr_jp2*vr_jp2
      vt_jp2    = wt_jp2 + unowjp2
      v_jp2_sq  = vm_jp2*vm_jp2 + vt_jp2*vt_jp2
c
      if(ifgas.eq.0) then
              roe_jp2   = roisent(i)*(cv*tisent(i)  + 0.5*v_jp2_sq)
c
              if(itimst.eq.5.or.itimst.eq.6)
     &        roe_jp2   = ro_sub_jp2(i)*(cv*tisent(i)  + 0.5*v_jp2_sq)
c
      else
              eisent    = hstat(i) - rgas*tisent(i) + 0.5*v_jp2_sq
              roe_jp2   = eisent*roisent(i)
      end if
c
c*****************************************************************************************
c*****************************************************************************************  
c     relax the primary variables at j = jmix + 2 towards the "isentropic"  values 
c
      ro(i,jp2,k)    = rfmix1*ro(i,jp2,k)   + rfmix*ro_jp2
      rovx(i,jp2,k)  = rfmix1*rovx(i,jp2,k) + rfmix*ro_jp2*vx_jp2
      rovr(i,jp2,k)  = rfmix1*rovr(i,jp2,k) + rfmix*ro_jp2*vr_jp2
      rorvt(i,jp2,k) = rfmix1*rorvt(i,jp2,k)+ rfmix*ro_jp2*vt_jp2*rjp2
      roe(i,jp2,k)   = rfmix1*roe(i,jp2,k)  + rfmix*roe_jp2
c
  300 continue      
c
c******************************************************************************************
c
  400 continue
c
c******************************************************************************************
c
  500 continue
c
c******************************************************************************************
c
      return
      end
c
c******************************************************************************************
c******************************************************************************************
c******************************************************************************************
c
      subroutine insect(jd,maxki,iin,xi,yi,jin,xj,yj,xinsect,yinsect)
c
c   this subroutine finds the intersection point of two lines defined by their coordinates
c   xi,yi   and  xj,yj .
c
      dimension   xi(jd),yi(jd),xj(maxki),yj(maxki),dperp(500)
c
c
c   first find the two closest points on the two curves.
c
      dminov = 1.0e12
      do i=1,iin
      dmin_i= 1.0e12
c
      do j=1,jin
      xdif = xj(j) - xi(i)
      ydif = yj(j) - yi(i)
      dist = xdif*xdif + ydif*ydif
      if(dist.lt.dmin_i) then
           jmin   = j
           dmin_i = dist
      end if
      end do
c
      if(dmin_i.lt.dminov) then
           iminov = i
           jminov = jmin
           dminov = dmin_i
      end if 
      end do
c
c      write(6,*) ' the closest distance between points =', sqrt(dminov)
c      write(6,*) ' at  i = ',iminov, ' j= ', jminov
c
c    form a simple estimate of the slope of the curve at the closest point.
c
      i = iminov
      im1 = i-1
      if(i.eq.1)   im1 = 1
      ip1 = i+1
      if(i.eq.iin) ip1 = iin
      dx = xi(ip1) - xi(im1)
      dy = yi(ip1) - yi(im1)
      dydx = dy/dx
      ds   = sqrt(dx*dx + dy*dy)
      vecx =  dx/ds
      vecy =  dy/ds
c
c   decide which is the next to closest point
c   and form a better estimate of slope of the curve at the closest point.
c
c
      proj = (xj(jminov)-xi(i))*vecx + (yj(jminov)-yi(i))*vecy
c
      if(proj.gt.0.0) then
               dx = xi(ip1) - xi(i)
               dy = yi(ip1) - yi(i)
      else
               dx = xi(i) - xi(im1)
               dy = yi(i) - yi(im1)
      end if
c
c    now find a unit normal vector to the "i"  line at the closest point.
c
      dydx = dy/dx
      ds   = sqrt(dx*dx + dy*dy)
      vecx =  dx/ds
      vecy =  dy/ds
      vnormx = -vecy
      vnormy =  vecx
c
c    now find the perpendicular distance between the closest point on the "i" line
c    and all points on the "j" line
c
      do 20 j=1,jin
      xdif = xj(j) - xi(i)
      ydif = yj(j) - yi(i)
      dperp(j) = xdif*vnormx + ydif*vnormy
   20 continue
c
c     interpolate to find the point where the perpendicular distance = zero.
c
      call intp(jin,dperp,xj,0.0,xinsect)
      call intp(jin,dperp,yj,0.0,yinsect)
c
c      write(6,*) ' xinsect, yinsect = ', xinsect,yinsect
c
      return
c
      end
c
c******************************************************************************************
c******************************************************************************************
c******************************************************************************************
c  
      subroutine smooth_var(d)   
c
      include 'commall-open-19.2'
c
      dimension d(id,jd,kd),avg(jd),curve(jd),scurve(jd)
      num3  = 3
c
c******************************************************************************************
c
c     apply combined 2nd and 4th order streamwise smoothing to the variable d
c
c******************************************************************************************
c
      if(sfx.lt.0.0001) go to 4001
c
c      streamwise smoothing with constant coefficient, sfx
c
      do 3990 nr = 1,nrows
           j1 = jstart(nr)
           j2 = jmix(nr)
      do 3900 j = j1+3,j2-3
      do 3900 k=1,km
      do 3900 i=1,im
           averg = 0.5*(d(i,j+1,k)+d(i,j-1,k))
           curv  = averg - 0.5*d(i,j,k) - 0.25*(d(i,j-2,k)+d(i,j+2,k))
           store(i,j,k) = sfx1*d(i,j,k) + sfx*(averg + fac_4th*curv)
 3900 continue
c******************************************************************************************
c******************************************************************************************
c      extra smoothing at the upstream, downstream and and mixing plane boundaries.
c      this modified by jdd  march 2015. new version gives better results for supersonic
c      flows at the mixing plane
c
      do 3905 k=1,km
      do 3905 i=1,im
      store(i,j2,k)   = d(i,j2,k)
      store(i,j1,k)   = d(i,j1,k)
 3905 continue
c
      do 3910 k=1,km
      do 3910 i=1,im
c
      store(i,j1+2,k) =   sfxb1*d(i,j1+2,k) 
     &                +  sfxbh*(d(i,j1+3,k) + d(i,j1+1,k)  )
c
      store(i,j1+1,k) = sfxb1*d(i,j1+1,k) + sfxb*(1.25*d(i,j1+2,k)
     &                 + 0.5*d(i,j1+3,k) -  0.75*d(i,j1+4,k) )
c 
      store(i,j2-2,k) =  sfxb1*d(i,j2-2,k)
     &                + sfxbh*(d(i,j2-1,k)  + d(i,j2-3,k)  )
c
      store(i,j2-1,k) = sfxb1*d(i,j2-1,k) + sfxb*(1.25*d(i,j2-2,k)
     &                  + 0.5*d(i,j2-3,k) - 0.75*d(i,j2-4,k) )
c
 3910 continue
c
 3990 continue
c
c     end of streamwise smoothing
c
c*******************************************************************************
c*******************************************************************************
c
c      now apply special smoothing around the leading edge
c
      do 4006 nr = 1,nrows
      js = jled(nr)-2
      je = js + 1
      do 4004 j = js,je
      do 4004 k=1,km
      store(1,j,k)  = sfxhm1*d(1,j,k) + sfx14*
     &                         (d(2,j,k) + d(imm1,j,k))
      store(im,j,k) = sfxhm1*d(im,j,k) + sfx14*
     &                         (d(2,j,k)+d(imm1,j,k))
 4004 continue
      j = jled(nr)
      do 4007 k=1,km
      store(1,j,k)  = sfxhm1*d(1,j,k)+sfx14*
     &                      (d(1,j+1,k)+d(im,j,k))
      store(im,j,k) = sfxhm1*d(im,j,k)+sfx14*
     &                       (d(im,j+1,k)+d(1,j,k))
 4007 continue
 4006 continue
c*******************************************************************************
c     re-set the variable   "d"  to the new smoothed value..
c
      do 4015 nr = 1,nrows
      j1 = jstart(nr)
      j2 = jmix(nr)
      do 4015 j  = j1,j2
      do 4015 k  = 1,km
      do 4015 i  = 1,im
      d(i,j,k)   = store(i,j,k)
 4015 continue
c
c*******************************************************************************
c*******************************************************************************
c     re enter here if no streamwise smoothing.
c
 4001 continue
c
c*******************************************************************************
c*******************************************************************************
c
      if(sft.lt.0.0001) go to 9500
c      
c          now perform the pitchwise and spanwise smoothing.
c
      do 9000 j=2,jm
c  tflow addition
      if(im.eq.2) go to 9101
c
      do 9100 k=1,km
      do 9110 i=2,imm1
      avg(i)    = fd(i)*d(i+1,j,k)+fu(i)*d(i-1,j,k)
      curve(i)  = d(i,j,k)-avg(i)
 9110 continue
c
      if(ind(j).eq.1.or.indle(j).eq.1) go to 9112
c    set values on the periodic boundaries
      avg(1)    = 0.5*(d(2,j,k)+d(imm1,j,k))
      avg(im)   = avg(1)
      curve(1)  = d(1,j,k)-avg(1)
      curve(im) = curve(1)
      scurve(1) = 0.5*(curve(2)+curve(imm1))
      scurve(im)= scurve(1)  
      go to 9113
c
 9112 continue
c     set values on the blade surfaces
      avg(1)    = d(2,j,k)    + fu(1)*(d(2,j,k)-d(num3,j,k))
      avg(im)   = d(imm1,j,k) + fu(im)*(d(imm1,j,k)-d(imm2,j,k))
      curve(1)  = 0.0
      curve(im) = 0.0
      scurve(1) = 0.0
      scurve(im)= 0.0
c
 9113 continue
c   smooth the second derivative, curve, to form scurve .
      do 9120 i=2,imm1
      scurve(i) = fu(i)*curve(i-1)+fd(i)*curve(i+1)
 9120 continue
c
c     apply the final pitchwise smoothing.
      do 9130 i=1,im
      d(i,j,k) = sft1*d(i,j,k) +sft*(avg(i) + fac_4th*scurve(i))
 9130 continue
c
 9100 continue
c  tflow addition
 9101 continue
c*******************************************************************************
c     start to apply the spanwise smoothing
c  q3d  
      if(km.eq.2) go to 9301
c  end q3d
c
      do 9300 i=1,im
      do 9200 k=2,kmm1
      avg(k)   = fkd(k)*d(i,j,k+1) + fku(k)*d(i,j,k-1)
      curve(k) = d(i,j,k)-avg(k)
 9200 continue
c     form values on the end walls
      avg(1)   = d(i,j,2) +0.5*fku(1)*(d(i,j,2)-d(i,j,3))
      avg(km)  = d(i,j,kmm1)+0.5* fku(km)*(d(i,j,kmm1)-d(i,j,kmm2))
      curve(1) = 0.0
      curve(km)= 0.0
      scurve(1)= 0.0
      scurve(km)=0.0
c     smooth the second derivative
      do 9210 k=2,kmm1
      scurve(k) = fku(k)*curve(k-1) + fkd(k)*curve(k+1)
 9210 continue
c
c    apply the final spanwise smoothing
      do 9230 k=1,km
      d(i,j,k)  = sft1*d(i,j,k) + sft*(avg(k)+fac_4th*scurve(k))
 9230 continue
c
 9300 continue
c
 9301 continue
c
c*******************************************************************************
c   smoothing of the corner points using sft .
c
      d(1,j,1)  = sft1*d(1,j,1)  + sfth*(d(2,j,1)+d(1,j,2))
      d(im,j,1) = sft1*d(im,j,1) + sfth*(d(imm1,j,1)+d(im,j,2))
      d(1,j,km) = sft1*d(1,j,km) + sfth*(d(1,j,kmm1)+d(2,j,km))
      d(im,j,km)= sft1*d(im,j,km)+ sfth*(d(imm1,j,km)+d(im,j,kmm1))
 9000 continue
c
c*******************************************************************************
c
c      re enter here if no pitchwise or spanwise smoothing
 9500 continue
c
c*******************************************************************************
c*******************************************************************************
c      smooth the exit flow if required, nsfexit points are smoothed by sfexit.
c
      if(sfexit.lt.0.0001) go to 8700
c      
      jsstart = jm + 1 - nsfexit
      do 8502 j= jsstart,jm
c   smooth in the "i" direction
c
      if(im.eq.2) go to 8505
c
      do 8500 k=1,km
      do 8501 i=2,imm1
      d(i,j,k) = sfex1*d(i,j,k)+sfexh*(d(i+1,j,k)+d(i-1,j,k))
 8501 continue
      d(1,j,k) = d(1,j,k) +sfexit*(2.*d(2,j,k)-d(num3,j,k)-d(1,j,k))
      d(im,j,k)= d(im,j,k)+sfexit*(2.*d(imm1,j,k)-d(imm2,j,k)-d(im,j,k))
 8500 continue
c
 8505 continue
c
c  q3d
      if(km.eq.2) go to 8602
c  end q3d    
c     smooth in the "k" direction
      do 8600 i=1,im
      do 8601 k=2,kmm1
      d(i,j,k) = sfex1*d(i,j,k)+sfexh*(d(i,j,k-1) +d(i,j,k+1))
 8601 continue
      d(i,j,1) = d(i,j,1) +sfexit*(2.*d(i,j,2)-d(i,j,3)-d(i,j,1))
      d(i,j,km)= d(i,j,km)+sfexit*(2.*d(i,j,kmm1)-d(i,j,kmm2)-d(i,j,km))
 8600 continue
c
 8602 continue
c
 8502 continue
c
c
 8700 continue
c
c******************************************************************************
c     end of exit flow smoothing
c******************************************************************************
      return
      end
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c
      subroutine re_design(nr,k,j1,j2,jlerow,jterow) 
c
      include 'commall-open-19.2'
c
      dimension  smer(jd),fracnew(jd),betanew(jd),
     &           slope(jd),thickup(jd),thicklow(jd),tkup(jd),tklow(jd),
     &		 xss(jd),rss(jd),s_ss(jd),relspace(jd),s_rel(jd),
     &           space(jd),frac_chord(jd),xnew(jd),rnew(jd),
     &           ssdist(jd)
c
c*******************************************************************************
c     subroutine to design a new blade section on the stream surface input below.
c*******************************************************************************
c
      write(6,*)
      write(6,*) ' using the blade redesign option on section number ',k
      write(6,*) ' j1, j2,  jlerow, jterow = ', j1,j2,jlerow,jterow
c
      jlee   = j1+jlerow-1
      jtee   = j1+jterow-1
      jmrow  = j2 - j1 + 1
      nxtrap = 4
c
      write(6,*) ' absolute values of jlee and jtee= ',jlee,jtee
      write(6,*) ' jm for the row being designed is  ', jmrow
c
c      input the new stream surface coordinates and relative grid spacing.
c      mark the position of the blade leading edge and the trailing edge.
c
      read(5,*)    dummy_input
      read(5,*) n_ss, n_le, n_te
c
      read(5,*)   dummy_input
      do n = 1,n_ss
      read(5,*) xss(n),rss(n),relspace(n)
      end do
c  
c
c   input the number of points at which a new blade camber line and thickness will de input.
c     also the number of times they will be smoothed.
c
      read(5,*)    dummy_input
      read(5,*)    nnew, nsmooth
      write(6,*)
c
c   input the new camber line slope and upper and lower tangential thickness.
c   as fractions of the meridional chord. 
c
      read(5,*)    dummy_input
      do  nn = 1,nnew
      read(5,*)     fracnew(nn),betanew(nn),thickup(nn),thicklow(nn)
      write(6,*)  ' fracnew(nn),betanew(nn),thickup(nn),thicklow(nn) ',
     &              fracnew(nn),betanew(nn),thickup(nn),thicklow(nn)
      end do
      write(6,*)
c
      read(5,*)   dummy_input
      read(5,*)   frac_chord_up, frac_chord_dwn, rtheta_mid
      write(6,*) 'frac_chord_up, frac_chord_dwn, rtheta_mid',
     &            frac_chord_up, frac_chord_dwn, rtheta_mid
      write(6,*)
c
c    smooth the input values of blade camber line angle and thickness.      
      fsmooth = 0.25
      call smooth(1,nnew,nsmooth,fsmooth,fracnew,betanew)
      call smooth(1,nnew,nsmooth,fsmooth,fracnew,thickup)
      call smooth(1,nnew,nsmooth,fsmooth,fracnew,thicklow)
c
c    find the meridional distance of the points on the stream surface.
      s_ss(1) = 0.0
      do n = 2,n_ss
      xdif = xss(n) - xss(n-1)
      rdif = rss(n) - rss(n-1)
      s_ss(n) = s_ss(n-1) + sqrt(xdif*xdif + rdif*rdif)
      end do
c
c   set the meridional chord
      s_merid  = s_ss(n_te) - s_ss(n_le)
c
c     scale the thickness by the meridional chord
      do nn = 1,nnew 
           tkup(nn)  = thickup(nn)*s_merid
           tklow(nn) = thicklow(nn)*s_merid
      end do
c
c    non-dimensionalise the meridional distance by the meridional chord
c    s_rel(n_le) = 0, s_rel(n_te) = 1.
      do n = 1,n_ss
      s_rel(n) = ( s_ss(n) - s_ss(n_le) )/s_merid
      end do
c
c    interpolate to find the relative grid spacings at the final grid points.
c    this is for all grid points.
      do j = 1,jmrow
      fracj = float(j-jlerow)/float(jterow - jlerow)
      call intp(n_ss,s_rel,relspace,fracj,space(j))
      end do
c
c    find the fraction of meridional chord at the stream surface points input
c    on the blade surface.
      frac_chord(jlerow) = 0.0
      do j = jlerow+1,jterow
      frac_chord(j) = frac_chord(j-1) + 0.5*(space(j)+space(j-1))
      end do
c     make sure that frac_chord = 0 at the  le  and  = 1 at the te.
      do j = jlerow,jterow
      frac_chord(j) = frac_chord(j)/frac_chord(jterow)
      end do
c
c     set the final meridional positions of the grid points on the blade surface
      do j = jlerow,jterow
      ssdist(j) = s_ss(n_le) + s_merid*frac_chord(j)
      end do
c
c    find the fraction of meridional chord at the stream surface points input
c    upstream of the le.
      frac_chord(1) = 0.0
      do j = 2,jlerow
      frac_chord(j) = frac_chord(j-1) + 0.5*(space(j)+space(j-1))
      end do

c     make sure that frac_chord = 0 at the  le  and  = 1 at the leading edge.
      do j = 1,jlerow
      frac_chord(j) = frac_chord(j)/frac_chord(jlerow)
      end do
c
c     set the final meridional positions of the grid points upstream of the le.
      do j = 1,jlerow-1
      ssdist(j) = s_ss(n_le) + frac_chord_up*s_merid*(frac_chord(j)-1.0)
      end do
c
c    find the fraction of meridional chord at the stream surface points input
c    downstream of the trailing edge.
      frac_chord(jterow) = 0.0
      do j = jterow+1,jmrow
      frac_chord(j) = frac_chord(j-1) + 0.5*(space(j)+space(j-1))
      end do

c     make sure that frac_chord = 0 at the  te  and  = 1 at the exit.
      sum_frac = frac_chord(jmrow) - frac_chord(jterow)
      do j = jterow,jmrow
      frac_chord(j) = (frac_chord(j) - frac_chord(jterow))/sum_frac
      end do
c
c     set the final meridional positions of the grid points downstream of the te.
      do j = jterow+1,jmrow
      ssdist(j) = ssdist(jterow) + frac_chord_dwn*s_merid*frac_chord(j)
      end do
c
c   interpolate to find the x  and  r  coordinates of the final grid.
      do j = 1,jmrow
      jall  = j1 + j - 1
      call intp(n_ss,s_ss,xss,ssdist(j),xsurf(jall,k))
      call intp(n_ss,s_ss,rss,ssdist(j),rsurf(jall,k))
      end do 
c
c   interpolate in the new blade to obtain the coordinates at the grid points.
c 
      do j = jlerow,jterow
      jall = j + j1 - 1
      fracs = (ssdist(j)-ssdist(jlerow))/(ssdist(jterow)-ssdist(jlerow))
      call intp(nnew,fracnew, betanew,fracs, slope(jall))
      call intp(nnew,fracnew, tkup, fracs, thickup(jall))
      call intp(nnew,fracnew, tklow,fracs, thicklow(jall))
      end do
c
c     form the grid upstream of the leading edge
      dthdx  = tan(slope(jlee+nxtrap)*degrad)/rsurf(jlee+nxtrap,k)
      do j = j1,jlee
      distup        = ssdist(j) - ssdist(jlerow)
      rt_upp(j,k)   = distup*dthdx*rsurf(j,k)
      rt_thick(j,k) = 0.0
      end do
c
c    form the grid on the blade
      dthdxp = tan(slope(jlee)*degrad)/rsurf(jlee,k)
      th_mid_p = rt_upp(jlee,k)/rsurf(jlee,k)
      do j   = jlee + 1,jtee
      dthdx  = tan(slope(j)*degrad)/rsurf(j,k)
      th_mid = th_mid_p + 0.5*(dthdx+dthdxp)*(ssdist(j) - ssdist(j-1) )
      rt_upp(j,k)   = th_mid*rsurf(j,k)  + thickup(j)
      rt_thick(j,k) = thickup(j) + thicklow(j)
      dthdxp        = dthdx
      th_mid_p      = th_mid
      end do
c
c     form the new grid downstream of the trailing edge
      dthdx  = tan(slope(jtee-nxtrap)*degrad)/rsurf(jtee-nxtrap,k)
      do j = jtee + 1, j2
      distdwn       = ssdist(j) - ssdist(jterow)
      rt_upp(j,k)   = rt_upp(jtee,k)  + distdwn*dthdx*rsurf(j,k)
      rt_thick(j,k) = 0.0
      end do                 
c
c     set rtheta = rtheta_mid  at the mid j point of this section
c
      jmdd  = (jlee + jtee)/2
      thshift = (rtheta_mid
     &       - (rt_upp(jmdd,k)- 0.5*rt_thick(jmdd,k)))/rsurf(jmdd,k)
      do j = j1,j2
      rt_upp(j,k) = rt_upp(j,k) + thshift*rsurf(j,k)
      end do
c
c    write out the new blade coordinates.
c
      write(6,*)
      write(6,*)'new blade coordinates, j   x , rt_upp, rt_thick, rsurf'
      do j=j1,j2
      write(6,666) j, xsurf(j,k),rt_upp(j,k),rt_thick(j,k),rsurf(j,k)
      end do
  666 format(i5,4f20.5)
c
c*********************************************************************************** 
c     end of blade section redesign
c***********************************************************************************     
c 
      return
      end
c
c*****************************************************************************
c
c*****************************************************************************
      subroutine set_ssthick(j1,j2)
c
c  this subroutinr modified to allow for multiple blade rows. august 2018.
c
      include 'commall-open-19.2'
c
      dimension fracss(20),tkss(20),tk_ss(jd),sdist(jd)
      save      srange1,tkss_last
c
      write(6,*) '  inputting the data for a q3d calculation '
      write(6,*)
      read(5,*)    dummy_input
      read(5,*)    q3dforce
      write(6,*) ' q3dforce = ', q3dforce
      read(5,*)    nss
      write(6,*) ' stream surface data '
      read(5,*)   (fracss(n), n=1,nss)
      write(6,*)  'fracss',( fracss(n), n=1,nss)
      read(5,*)   (tkss(n),n=1,nss)
      write(6,*)  '  tkss',( tkss(n),n=1,nss)
      write(6,*)
c
      sdist(1) = 0.0
      do 10 j = j1+1, j2
           jloc = j+1- j1
           xd = xsurf(j,1)  - xsurf(j-1,1)
           rd = rsurf(j,1)  - rsurf(j-1,1)
           sdist(jloc) = sdist(jloc-1) + sqrt(xd*xd+rd*rd)
   10 continue
c
      nval = j2+1-j1
c
      srange   = sdist(nval) - sdist(1)
      if(j1.eq.1) srange1 = srange
c
      do 20  j = 1, nval
      sdist(j) = (sdist(j)-sdist(j1))/srange
   20 continue
c
c  set the inlet ss thickness = 0.05 x the meridional length ogf the first row.
      tkss1 = tkss(1)
      do 25 n=1,nss
      tkss(n) = 0.05*tkss(n)/tkss1*srange1
   25 continue
c
cc   make the ss thickness continuous at a mixing plane
           
      if(j1.gt.1) then
           tkss1 = tkss(1)
           do 26 n=1,nss
           tkss(n) = tkss(n)*tkss_last/tkss1
   26      continue
      end if
      tkss_last = tkss(nss)
c   
      do 30 j = 1,nval 
      call intp(nss,fracss,tkss,sdist(j),tk_ss(j))
   30 continue
c
c  smooth the stream surface thickness
c
      call smooth(1,nval,4,0.25,sdist,tk_ss)
c
      tk_ss(1)     = tk_ss(2)
      tk_ss(nval)  = tk_ss(nval-1) 
c
      do 40 j= j1,j2
      jp1 = j+1
      if(j.eq.j2) jp1 = j2
      jm1 = j-1
      if(j.eq.j1) jm1 = j1
      jloc   = j+1-j1
      jlocp1 = jloc +1
      jlocm1 = jloc -1
      if(j.eq.1)  jlocm1 = 1
      if(j.eq.j2) jlocp1 = jloc
      drds = (rsurf(jp1,1) -  rsurf(jm1,1))
     &     /( srange*(sdist(jlocp1) - sdist(jlocm1)) )
      dxds = (xsurf(jp1,1) -  xsurf(jm1,1))
     &     /( srange*(sdist(jlocp1) - sdist(jlocm1)) )
      rsurf(j,2)    = rsurf(j,1) + tk_ss(jloc)*dxds
      xsurf(j,2)    = xsurf(j,1) - tk_ss(jloc)*drds
      rt_upp(j,2)   = rt_upp(j,1)*rsurf(j,2)/rsurf(j,1)
      rt_thick(j,2) = rt_thick(j,1)
   40 continue  
c      
      return
      end
c
c************************************************************************************
c************************************************************************************
c
      subroutine restagger(k,j1,j2,jlerow,jterow,rotate,fracx_rot)
c
c************************************************************************************
c     this subroutine restaggers a blade section.
c     the rotation is applied on a stream surface so that it is applicable
c     to both radial flow and axial flow machines.
c************************************************************************************
      include 'commall-open-19.2'
c
      dimension  sdis_t(jd),thet_a(jd),slope(jd)
c
      rotate     = 0.0
      fracx_rot  = 0.5
      read(5,*) dummy_input
      read(5,*,err=1502) rotate, fracx_rot
 1502 continue
      write(6,*) ' input for restaggering option'
      write(6,*) ' rotation =', rotate, 'about ',fracx_rot 
      write(6,*)
c
      rad_rot = rotate*degrad 
c
c  calculate the meridional distances and angles.
c
      sdis_t(j1) = 0.0  
      thet_a(j1) =  (rt_upp(j1,k) - 0.5*rt_thick(j1,k))/rsurf(j1,k)   
      do 10 j = j1+1,j2
      xdif = xsurf(j,k) - xsurf(j-1,k)
      rdif = rsurf(j,k) - rsurf(j-1,k)
      sdis_t(j) = sdis_t(j-1) + sqrt(xdif*xdif + rdif*rdif)
      thet_a(j) = (rt_upp(j,k)- 0.5*rt_thick(j,k))/rsurf(j,k) 
   10 continue
c
c     find the j value at the centre of rotation, jrot .
c
      jrot = 1
      jlee = j1+jlerow-1
      jtee = j1+jterow-1
      s_merid = sdis_t(j2)
      do 20 j = j1,j2
      fracrd = (sdis_t(j) - sdis_t(jlee))/ (sdis_t(jtee) - sdis_t(jlee))
      if(fracrd.gt.fracx_rot)  go to 30
   20 continue
   30 continue
      jrot = j
      rtrot = rt_upp(jrot,k)- 0.5*rt_thick(jrot,k)
c
c     calculate the centre line slopes on the stream surface
c     defined as as  " r d(theta)/ds  "
c       
      do 40 j = j1+1,j2-1
      r_dt_ds = rsurf(j,k)*(thet_a(j+1) - thet_a(j-1) )
     &         / (sdis_t(j+1) - sdis_t(j-1) )
      slope(j) = atan(r_dt_ds)
   40 continue
      slope(j1) = slope(j1+1)
      slope(j2) = slope(j2-1)
c
c     rotate the blade on the stream surface so that all centre line
c     slopes change by "rotate" degrees in a clockwise sense.
c     maintain the same perpendicular thickness.
c
      thmid   = thet_a(j1)
      do 50 j = j1,j2
      jm1 = j-1
      if(j.eq.1) jm1 = 1
      avgslope = 0.5*(slope(j) + slope(jm1)) - rad_rot
      angnew   = slope(j) - rad_rot
      rt_thick(j,k) = rt_thick(j,k)*cos(slope(j))/cos(angnew)
      thmid = thmid 
     &      + (sdis_t(j) - sdis_t(jm1))*tan(avgslope)/rsurf(jm1,k)
      rt_upp(j,k) = rsurf(j,k)*thmid + 0.5*rt_thick(j,k)
   50 continue
c
c      rotate as a solid body so that the centre of rotation, jrot, does not move.
c 
      thshift = (rtrot
     &        - (rt_upp(jrot,k)- 0.5*rt_thick(jrot,k)))/rsurf(jrot,k)
      do 60 j = j1,j2
      rt_upp(j,k) = rt_upp(j,k) + thshift*rsurf(j,k)
   60 continue
c
      return
      end
c
c*********************************************************************************
c*********************************************************************************
c*********************************************************************************
      subroutine wallfun(i,j,k,nwall,perp,dpds,density,twall,
     &                   yplus_old,wrel_wall,yplus_new)
c
c   this subroutine applies the wall functions proposed by shih et al. nasa/tm-1999-209398 
c
      include  'commall-open-19.2'
c
c    do not use the pressure gradient term unless ypluswall is less than -10.
      if(ypluswall.gt.-10.0) go to 100
c
      if(abs(dpds).lt.1.0) dpds = 1.0
      sgn_dpds = abs(dpds)/dpds
c
      vpres    = (vislam/density/density* abs(dpds))**0.33333
      vstar    = yplus_old*vislam/(density*perp)
c
c     set limits to the pressure gradient velocity term.
      if(sgn_dpds.lt.0.0.and.vpres.gt.0.33*vstar) vpres = 0.33*vstar
      if(sgn_dpds.gt.0.0.and.vpres.gt.2.0*vstar)  vpres = 2.0*vstar
c 
      yplus_pres = vpres*perp*density/vislam
c  
c     set the pressure gradient term from yplus_pres.
      yp1     = yplus_pres
      yp2     = yp1*yp1
      yp3     = yp2*yp1
      yp4     = yp2*yp2
      if(yp1.le.4.0)                  fyplusp = 0.5*yp2 - 7.31e-3*yp3
      if(yp1.gt.4.0.and.yp1.le.15.0)  fyplusp = -15.138 + 8.4688*yp1
     &                 - 0.81976*yp2 + 3.7292e-2*yp3 - 6.3866e-4*yp4
      if(yp1.gt.15.0.and.yp1.le.30.0) fyplusp = 11.925 + 0.934*yp1
     &             - 2.7805e-2*yp2 + 4.6262e-4*yp3 - 3.1442e-6*yp4
      if(yp1.gt.30.0)                 fyplusp = 5*alog(yp1) + 8.0
c
  100 continue
c
c    re-enter here if ypluswall is between -5 and -10.
c    set the wall distance-velocity  term from yplus_old
      yp1      = yplus_old
      yp2      = yp1*yp1
      yp3      = yp2*yp1
      yp4      = yp2*yp2 
      if(yp1.le.5.0) fyplus = yp1 + 1.0e-2*yp2 -2.9e-03*yp3
      if(yp1.gt.5.0.and.yp1.le.30.0) fyplus = -0.872 + 1.465*yp1
     &                 -7.02e-2*yp2 + 1.66e-03*yp3 - 1.495e-5*yp4
      if(yp1.gt.30.and.yp1.le.140.0)  fyplus = 8.6 + 0.1864*yp1
     &             - 2.006e-3*yp2 + 1.144e-5*yp3 -2.551e-8*yp4
      if(yp1.gt.140.0)  fyplus   = 2.439*alog(yp1) + 5.0
c
c
      if(ypluswall.gt.-10.0) then
           fracp = 0.0
      else
           fracp = sgn_dpds*(vpres/vstar) * fyplusp/fyplus
      end if
c
      vplus     = fyplus*(1.0 + fracp)
c
      vstar     = wrel_wall/vplus
c      
      twall     = density*vstar*vstar
c
c  set the new value of yplus. this is used on the next call of this subroutine
c  and the iteration  should converge to the final value.
      yplus_new = density*vstar*perp/vislam
c
      return
      end
c
c********************************************************************************* 
c*********************************************************************************
c      
      subroutine set_pwallgrad
c
c    this subroutine evaluatres the pressure gradients on thw walls for use with
c    the shih et al wall functions
c
      include  'commall-open-19.2'
c
      dimension  tempp(jd)
c
c   the limit on pressure gradient is if the average pressure change per blade row
c   occurs over 10% of the chord. 
c    jdd corrected error in dpds_lim  20/9/18 . 
      dpds_lim = 20.0*(po1(kmid) - 0.5*(pdown_hub+pdown_tip))/nrows/
     &                (chord(1)  + chord(nrows))
      dpds_lim =  abs(dpds_lim)
c
c     calculate the pressure gradient along a streamline
      do 226 j = 2,jm
      do 226 k = 1,kmm1
      do 226 i = 1,imm1
      call gradvel(i,j,k,peff,dpdx,dpdr,dpdt)
      vxcell = vx(i,j,k)+vx(i,j-1,k)+vx(i+1,j,k)+vx(i+1,j-1,k) 
     &       + vx(i,j,k+1)+vx(i,j-1,k+1)+vx(i+1,j,k+1)+vx(i+1,j-1,k+1)
      vrcell = vr(i,j,k)+vr(i,j-1,k)+vr(i+1,j,k)+vr(i+1,j-1,k)
     &       + vr(i,j,k+1)+vr(i,j-1,k+1)+vr(i+1,j,k+1)+vr(i+1,j-1,k+1)
      wtcell = wt(i,j,k)+wt(i,j-1,k)+wt(i+1,j,k)+wt(i+1,j-1,k)
     &       + wt(i,j,k+1)+wt(i,j-1,k+1)+wt(i+1,j,k+1)+wt(i+1,j-1,k+1)
      wrel_cell = sqrt(vxcell*vxcell + vrcell*vrcell + wtcell*wtcell)
      dpds      = (dpdx*vxcell+dpdr*vrcell+dpdt*wtcell)/wrel_cell
      if(abs(dpds).gt.dpds_lim) dpds = dpds_lim*dpds/abs(dpds)
      dpds_cell(i,j,k) = dpds      
  226 continue
c 
c   smooth the pressure gradient along a grid line.
      do 250 nr = 1,nrows
c
      j1 = jstart(nr)
      j2 = jmix(nr)
      do 252 k=1,kmm1
      do 252 i=1,imm1
c
      dpds_cell(i,j2+1,k) = 0.0
      dpds_cell(i,j2+2,k) = 2*dpds_cell(i,j2+3,k) - dpds_cell(i,j2+4,k)
      dpds_cell(i,j2,k)   = 2*dpds_cell(i,j2-1,k) - dpds_cell(i,j2-2,k)
c
      do 251 j = j1+2,j2-1
      tempp(j) = 0.5*dpds_cell(i,j,k)
     &         + 0.25*(dpds_cell(i,j-1,k) + dpds_cell(i,j+1,k))  
  251 continue 
      tempp(j2+1) = 0.0
      do 253 j = j1+3,j2-2
      dpds_cell(i,j,k)= 0.5*tempp(j) + 0.25*(tempp(j-1)+tempp(j+1))
  253 continue
c
      dpds_cell(i,j2+1,k) = 0.0
      dpds_cell(i,j2+2,k) = 2*dpds_cell(i,j2+3,k) - dpds_cell(i,j2+4,k)
      dpds_cell(i,j2,k)   = 2*dpds_cell(i,j2-1,k) - dpds_cell(i,j2-2,k)
c
  252 continue
c 
  250 continue 
c
      return
      end
c******************************************************************************
c****************************************************************************** 
c
      subroutine cool_input
c
c  this subroutine reads in the data forany cooling flows.
c
      include  'commall-open-19.2'
c
  99  format(a72)
c
      ncoolb     = 0
      ncoolw     = 0
      ncwlbladep = 0
      ncwlwallp  = 0
      pi_180     = 3.1415926/180
      pi_30      = 3.1415926/30
c******************************************************************************
c******************************************************************************
      do 100 n_row = 1,nrows
c
      j1      =  jstart(n_row)
c
      write(6,*) ' reading the cooling flow data for row no.', n_row
      read(5,99)  dummy_input
      write(6,99) dummy_input
c
      read(5,*)  ncwlblade, ncwlwall
      write(6,*) ncwlblade, ncwlwall
c
      if(ncwlblade.eq.0.and.ncwlwall.eq.0) go to 100
c
      ncoolb  = ncwlbladep + ncwlblade
      ncoolw  = ncwlwallp  + ncwlwall
c
c******************************************************************************
c******************************************************************************
c     read in the cooling flow details for the blade surfaces.
c
      if(ncwlblade.ne.0) then
c
      write(6,*) ' reading the blade surface cooling data for row no.',
     &             n_row
c
      read(5,99)  dummy_input
      write(6,99) dummy_input
c
      do 200 nbsurf = ncwlbladep+1,ncoolb
      read(5,*)  ic(nbsurf),jcbs(nbsurf),jcbe(nbsurf),kcbs(nbsurf),
     &           kcbe(nbsurf)
      read(5,*)  cflowb(nbsurf),tocoolb(nbsurf),pocoolb(nbsurf),
     &           machcool(nbsurf),sangleb(nbsurf),xangleb(nbsurf),
     &           rvt_in_b(nbsurf),rpm_cool
      write(6,*) ic(nbsurf),jcbs(nbsurf),jcbe(nbsurf),kcbs(nbsurf),
     &           kcbe(nbsurf)
      write(6,*) cflowb(nbsurf),tocoolb(nbsurf),pocoolb(nbsurf),
     &           machcool(nbsurf),sangleb(nbsurf),xangleb(nbsurf),
     &           rvt_in_b(nbsurf),rpm_cool
c
      sangleb(nbsurf)    = sangleb(nbsurf)*pi_180
      xangleb(nbsurf)    = xangleb(nbsurf)*pi_180
      wrad_coolb(nbsurf) = rpm_cool*pi_30
      jcbs(nbsurf)       = jcbs(nbsurf) + j1 -1
      jcbe(nbsurf)       = jcbe(nbsurf) + j1 -1
c
      if(tocoolb(nbsurf).lt.tcool_min) tcool_min = tocoolb(nbsurf)
c
  200 continue
c
      end if
c
c      jcbs & jcbe are the j values where cooling starts and ends for
c      this patch, the j indices are defined relative to the start
c      of the current blade row. ie the upstream mixing plane is j = 1.
c
c      cflowb  is the coolant mass flow through this patch in kg/s.
c
c      tocoolb is the absolute stagnation temperature of the coolant when
c      it is first supplied to the blade row. (note absolute).
c
c      machcool is the relative mach number of the coolant as it leaves the
c      blade and enters the mainstream flow. (note relative).
c
c      sangleb is the angle between the coolant jet and the blade surface.
c
c      xangleb is the angle between the projection of the cooling jet onto
c      the blade surface and the intersection of the blade surface with
c      a surface of constant radius (cylindrical surface).
c
c      rvt_in_b  is the angular momentum with which the cooling flow is supplied
c      disk chamber before entering the blades.
c
c      rpm_cool is the rotational speed of the surface through which the
c      coolant is ejected. the coolant velocity and directions are relative
c      to this rotational speed.
c******************************************************************************
c******************************************************************************
c     read in the cooling flow details for the endwall surfaces.
c******************************************************************************
c******************************************************************************
c
      if(ncwlwall.ne.0) then
c
      write(6,*) ' reading the endwall cooling data for row no.', n_row
c
      read(5,99)  dummy_input
      write(6,99) dummy_input
c
      do 300 nwsurf = ncwlwallp+1,ncoolw
      read(5,*)  kc(nwsurf),jcws(nwsurf),jcwe(nwsurf),icws(nwsurf),
     &           icwe(nwsurf)
      read(5,*)  cfloww(nwsurf),tocoolw(nwsurf),pocoolw(nwsurf),
     &           machcool(nwsurf),sanglew(nwsurf),tanglew(nwsurf),
     &           rvt_in_w(nwsurf),rpm_cool
      write(6,*) kc(nwsurf),jcws(nwsurf),jcwe(nwsurf),icws(nwsurf),
     &           icwe(nwsurf)
      write(6,*) cfloww(nwsurf),tocoolw(nwsurf),pocoolw(nwsurf),
     &           machcool(nwsurf),sanglew(nwsurf),tanglew(nwsurf),
     &           rvt_in_w(nwsurf),rpm_cool
c
      sanglew(nwsurf)    = sanglew(nwsurf)*pi_180
      tanglew(nwsurf)    = tanglew(nwsurf)*pi_180
      wrad_coolw(nwsurf) = rpm_cool*pi_30
      jcws(nwsurf)       = jcws(nwsurf) + j1 -1
      jcwe(nwsurf)       = jcwe(nwsurf) + j1 -1
c
      if(tocoolw(nwsurf).lt.tcool_min) tcool_min = tocoolw(nwsurf)
c
  300 continue
c
      end if
c
c     assume the maxi<um coolant ejection mach number = 1.0 .
      tcool_min = tcool_min/(1.+0.5*(ga-1.0))
c   end of endwall cooling flow data input.

c
c      jcws & jcwe are the j values where cooling starts and ends for
c      this patch, the j indices are defined relative to the start
c      of the current blade row. ie the upstream mixing plane is j = 1.
c
c      cfloww is the coolant mass flow through this patch in kg/s.
c
c      tocoolw is the absolute stagnation temperature of the coolant when
c      it is first supplied to the blade row. (note absolute).
c
c      machcool is the relative mach number of the coolant as it leaves the
c      endwall and enters the mainstream flow. (note relative).
c
c      sanglew is the angle between the coolant jet and the endwall surface.
c
c      tanglew is the angle between the projection of the cooling jet onto
c      the endwall surface and the intersection of the endwall surface with
c      the axial- radial plane ( i.e with the plane theta = constant).
c
c      rvt_in_w  is the angular momentum with which the cooling flow is supplied
c      disk chamber before entering the blades.
c
c      rpm_cool is the rotational speed of the surface through which the
c      coolant is ejected. the coolant velocity and directions are relative
c      to this rotational speed.
c
c

c
c******************************************************************************
c******************************************************************************
c
      ncwlbladep = ncoolb
      ncwlwallp  = ncoolw
c
      write(6,*)'blade row number ', n_row
      write(6,*)'number of cooling patches this row,on blade,on wall =',
     &           ncwlblade, ncwlwall
      write(6,*)'total cooling patches up to now, ncoolb,ncoolw = ',
     &           ncoolb, ncoolw
c
  100 continue
c
c    end of coolant data input loop.
c*************************************************************************
c*************************************************************************
      write(6,*) ' cooling flow data input complete. '
      write(6,*) ' total number of blade cooling patches   = ',ncoolb
      write(6,*) ' total number of endwall cooling patches = ',ncoolw
c
c*************************************************************************
c*************************************************************************
      return
      end
c*************************************************************************
c*************************************************************************
c*************************************************************************
c
      subroutine coolin_2
c
c  this subroutine sets the cooling flows based on the input coolant stagnation temperature
c  and pressure and the local static pressure at the point of ejection. the cooling flows
c  therefore change during convergence and must be reset every few iterations. the coolant
c  ejection mach number is input but not used.

      include  'commall-open-19.2'
c
      dimension ho_abs(jd,maxki),vnin(jd,maxki),vsin(jd,maxki),
     &          cell_flow(jd,maxki)
c
      do n = 1,nstages
           wpump(n) = 0.0
      end do
c
c*************************************************************************
c*************************************************************************
c
c      set the initial cooling flows to zero and
c      work out area magnitudes for coolant flows
c
      if(ncoolb.eq.0) go to 401
      do 400 j=2,jm
      do 400 k=1,kmm1
      cflowi1(j,k) = 0.0
      cflowim(j,k) = 0.0
      atoti1(j,k) = sqrt(abx(1,j,k)*abx(1,j,k) + abr(1,j,k)*abr(1,j,k)
     &                  + abt(j,k)*abt(j,k))
      atotim(j,k) = sqrt(abx(im,j,k)*abx(im,j,k)+abr(im,j,k)*abr(im,j,k)
     &                  + abt(j,k)*abt(j,k))
  400 continue
  401 continue
c
c
      if(ncoolw.eq.0) go to 501
      do 500 j=2,jm
      do 500 i=1,imm1
      cflowk1(i,j) = 0.0
      cflowkm(i,j) = 0.0
      atotk1(i,j)= sqrt(asx(i,j,1)*asx(i,j,1)  + asr(i,j,1)*asr(i,j,1))
      atotkm(i,j)= sqrt(asx(i,j,km)*asx(i,j,km)+asr(i,j,km)*asr(i,j,km))
  500 continue
  501 continue
c
c******************************************************************************
c******************************************************************************
c
c      now evaluate the cooling mass flows through the blade surfaces.
c      and the associated energy and momentum fluxes.
c
c******************************************************************************
c******************************************************************************
      do 600 ncb = 1,ncoolb
c
      if(ifgas.eq.0) then
           cpgas       = cp
           gagas       = ga
           fgagas      = fga
      else
           delt        = toin - tref
           cpgas       = cp1  + cp2*delt + cp3*delt*delt
           gagas       = cpgas/(cpgas - rgas)
           fgagas      = (gagas - 1.0)/gagas
      end if 
c
c
      toin        = tocoolb(ncb)
      poin        = pocoolb(ncb)
      w_preswirl  = wrad_coolb(ncb)*rvt_in_b(ncb)
      ssangleb    = sin(sangleb(ncb))
      csangleb    = cos(sangleb(ncb))
      patch_flow  = 0.0
c
c***********************************************************************
c***********************************************************************
c
      do 620 j = jcbs(ncb)+1,jcbe(ncb)
             n_row    = nrow(j)
             n_stage  = nstage(n_row)
      do 620 k = kcbs(ncb),kcbe(ncb)-1
c 
      if(ic(ncb).eq.1)
     &     pavg_in  = 0.25*(p(1,j,k)+p(1,j-1,k)+p(1,j,k+1)+p(1,j-1,k+1))
      if(ic(ncb).eq.im)
     &     pavg_in  = 0.25*(p(im,j,k)+p(im,j-1,k)+p(im,j,k+1)
     &                +p(im,j-1,k+1))
c
      ravg        = ravg_cell(j,k)
      vblade      = wrad_coolb(ncb)*ravg
      torel       = toin  + (0.5*vblade*vblade - w_preswirl)/cpgas
      toabs       = torel + 0.5*vblade*vblade/cpgas
      po_rel      = poin*(torel/toin)**(1./fgagas)
      po_ejectb(ncb) = po_rel
      pcool_rat    = pavg_in/po_rel
      if(pcool_rat.lt.0.5)    pcool_rat = 0.5
      if(pcool_rat.gt.0.9999) pcool_rat = 0.9999
      pavg_in     = pcool_rat*po_rel
      tcool_rat   = pcool_rat**fgagas
      tin         = torel*tcool_rat
      roin        = pavg_in/rgas/tin           
      vin         = sqrt(2.*cpgas*(torel - tin))
      ho_abs(j,k) = cpgas*toabs
      torelb(ncb) = torel
      tstat_ejectb(ncb) = tin
c
c   evaluate the normal and tangential coolant relative velocities
      vnin(j,k)     = vin*ssangleb
      vsin(j,k)     = vin*csangleb
c
c      evaluate the  mass flow through the patch if it fills the whole cell surface area.
c
      if(ic(ncb).eq.1)  cell_flow(j,k) = vnin(j,k)*roin*atoti1(j,k)
      if(ic(ncb).eq.im) cell_flow(j,k) = vnin(j,k)*roin*atotim(j,k)
c   calculate the total flow through the patch if the full area is used.
      patch_flow = patch_flow + cell_flow(j,k)*nblade(j) 
c
c   end of  j  , k  loops
  620 continue
c
c    evaluate the fraction of cell area used by the cooling flow in order to obtain the specified flow.
       frac_flow = cflowb(ncb)/patch_flow
c******************************************************************************
c******************************************************************************
c
c   now recalculate the coolant flows so that they sum to the specified total flow
      do 650 j = jcbs(ncb)+1,jcbe(ncb)
      do 650 k = kcbs(ncb),kcbe(ncb)-1
c
c   set the coolant flow through the cell faces so that they sum to the specified flow.
      cflowin        = cell_flow(j,k) * frac_flow
c
      ravg           = ravg_cell(j,k)
      vblade         = wrad_coolb(ncb)*ravg
      wpump(n_stage) = wpump(n_stage) + nblade(j)*
     &                 cflowin*(vblade*vblade - w_preswirl)
c
      if(ic(ncb).eq.1)  then      
           anorm    = atoti1(j,k)
           xnorm    = abx(1,j,k)/anorm
           tnorm    = abt(j,k)/anorm
           rnorm    = abr(1,j,k)/anorm
      endif
      if(ic(ncb).eq.im) then
           anorm    = atotim(j,k)
           xnorm    = abx(im,j,k)/anorm
           tnorm    = abt(j,k)/anorm
           rnorm    = abr(im,j,k)/anorm
      end if
c
      ang      =  xangleb(ncb)
      xtnorm   =  sqrt(xnorm*xnorm + tnorm*tnorm)
      sx       =  (tnorm*cos(ang) - rnorm*xnorm*sin(ang))/xtnorm
      st       = -(xnorm*cos(ang) + rnorm*tnorm*sin(ang))/xtnorm
      sr       =  sin(ang)*xtnorm
c
      if(ic(ncb).eq.1)  then
           xvel          =  vnin(j,k)*xnorm + vsin(j,k)*sx
           rvel          =  vnin(j,k)*rnorm + vsin(j,k)*sr
           tvel          =  vnin(j,k)*tnorm + vsin(j,k)*st
           cflowi1(j,k)  = cflowin
           hocwli1(j,k)  = cflowin*(ho_abs(j,k) + vblade*tvel)
           vxcwli1(j,k)  = cflowin*xvel
           vrcwli1(j,k)  = cflowin*rvel
           rvtcwli1(j,k) = cflowin*(vblade  +  tvel)*ravg
      else
           xvel          = -vnin(j,k)*xnorm + vsin(j,k)*sx
           rvel          = -vnin(j,k)*rnorm + vsin(j,k)*sr
           tvel          = -vnin(j,k)*tnorm + vsin(j,k)*st
           cflowim(j,k)  = cflowin
           hocwlim(j,k)  = cflowin*(ho_abs(j,k) + vblade*tvel)
           vxcwlim(j,k)  = cflowin*xvel
           vrcwlim(j,k)  = cflowin*rvel
           rvtcwlim(j,k) = cflowin*(vblade  +  tvel)*ravg
      end if
c   end of the  j,k  loop
  650 continue
c
c******************************************************************************
c     end of setting the blade surface cooling fluxes.
  600 continue
c
c**************************************************************************
c**************************************************************************
c**************************************************************************
c**************************************************************************
c**************************************************************************
c
c      now evaluate the cooling mass flows through the endwall surfaces.
c      and the associated energy and momentum fluxes.
c
      do 700 ncw = 1,ncoolw
c
      if(ifgas.eq.0) then
           cpgas       = cp
           gagas       = ga
           fgagas      = fga
      else
           delt        = toin-tref
           cpgas       = cp1 + cp2*delt + cp3*delt*delt
           gagas       = cpgas/(cpgas - rgas)
           fgagas      = (gagas - 1.0)/gagas
      end if 
c
c
      toin        = tocoolw(ncw)
      poin        = pocoolw(ncw)
      w_preswirl  = wrad_coolw(ncw)*rvt_in_w(ncw)
      ssanglew    = sin(sanglew(ncw))
      csanglew    = cos(sanglew(ncw))
      patch_flow  = 0.0
c
c***********************************************************************
c***********************************************************************
c
      do 720 j = jcws(ncw)+1,jcwe(ncw)
             n_row    = nrow(j)
             n_stage  = nstage(n_row)
      do 720 i = icws(ncw),icwe(ncw)-1
c 
      if(kc(ncw).eq.1) then
           pavg_in  = 0.25*(p(i,j,1)+p(i,j-1,1)+p(i+1,j,1)+p(i+1,j-1,1))
           ravg     = 0.5*(r(j,1)+r(j-1,1))
      end if
      if(kc(ncw).eq.km) then
           pavg_in  = 0.25*(p(i,j,km)+p(i,j-1,km)+p(i+1,j,km)
     &                    +p(i+1,j-1,km))
           ravg     = 0.5*(r(j,km) + r(j-1,km))
      end if

      vwall       = wrad_coolw(ncw)*ravg
      torel       = toin  + (0.5*vwall*vwall - w_preswirl)/cpgas
      toabs       = torel + 0.5*vwall*vwall/cpgas
      po_rel      = poin*(torel/toin)**(1.0/fgagas)
      po_ejectw(ncw) = po_rel
      pcool_rat   = pavg_in/po_rel
      if(pcool_rat.lt.0.5)    pcool_rat = 0.5
      if(pcool_rat.gt.0.9999) pcool_rat = 0.9999
      pavg_in     = pcool_rat*po_rel
      tcool_rat   = pcool_rat**fgagas
      tin         = torel*tcool_rat
      roin        = pavg_in/rgas/tin           
      vin         = sqrt(2.*cpgas*(torel - tin))
      ho_abs(j,i) = cpgas*toabs
c
      torelw(ncw)       = torel
      tstat_ejectw(ncw) = tin
c   evaluate the normal and tangential coolant relative velocities
      vnin(j,i)     = vin*ssanglew
      vsin(j,i)     = vin*csanglew
c
c      evaluate the  mass flow through the patch if it fills the whole cell surface area.
c      this is not the final coolant flow.
      if(kc(ncw).eq.1)  cell_flow(j,i) = vnin(j,i)*roin*atotk1(i,j)
      if(kc(ncw).eq.km) cell_flow(j,i) = vnin(j,i)*roin*atotkm(i,j)
c   calculate the total flow through the patch if the full area is used.
      patch_flow = patch_flow + cell_flow(j,i)*nblade(j) 
c
c   end of  i, j   loops
  720 continue
c
c    evaluate the fraction of cell area used by the cooling flow in order to obtain
c    the specified coolant flowrate.
       frac_flow = cfloww(ncw)/patch_flow
c******************************************************************************
c******************************************************************************
c
c   now recalculate the coolant flows so that they sum to the specified total flow
      do 750 j = jcws(ncw)+1,jcwe(ncw)
      do 750 i = icws(ncw),icwe(ncw)-1
c
c   set the coolant flow through the cell faces so that they sum to the specified flow.
      cflowin        = cell_flow(j,i) * frac_flow
c
      if(kc(ncw).eq.1)  ravg = 0.5*(r(j,1)  + r(j-1,1))
      if(kc(ncw).eq.km) ravg = 0.5*(r(j,km) + r(j-1,km))
      vwall          = wrad_coolw(ncw)*ravg
      wpump(n_stage) = wpump(n_stage) + nblade(j)*
     &                 cflowin*(vwall*vwall - w_preswirl)
c
      if(kc(ncw).eq.1) then
           anorm  = atotk1(i,j)
           xnorm  = asx(i,j,1)/anorm
           rnorm  = asr(i,j,1)/anorm
      else
           anorm  = atotkm(i,j)
           xnorm  = asx(i,j,km)/anorm
           rnorm  = asr(i,j,km)/anorm
      end if
c
      ang    =  tanglew(ncw)
      xrnorm =  sqrt(xnorm*xnorm + rnorm*rnorm)
      sx     =  cos(ang)*rnorm/xrnorm
      st     =  sin(ang)
      sr     = -cos(ang)*xnorm/xrnorm
c
      if(kc(ncw).eq.1)  then
           xvel          =  vnin(j,i)*xnorm + vsin(j,i)*sx
           rvel          =  vnin(j,i)*rnorm + vsin(j,i)*sr
           tvel          =  vsin(j,i)*st
           cflowk1(i,j)  = cflowin
           hocwlk1(i,j)  = cflowin*(ho_abs(j,i) + vwall*tvel)
           vxcwlk1(i,j)  = cflowin*xvel
           vrcwlk1(i,j)  = cflowin*rvel
           rvtcwlk1(i,j) = cflowin*(vwall  +  tvel)*ravg
      else
           xvel          =  -vnin(j,i)*xnorm + vsin(j,i)*sx
           rvel          =  -vnin(j,i)*rnorm + vsin(j,i)*sr
           tvel          =   vsin(j,i)*st
           cflowkm(i,j)  = cflowin
           hocwlkm(i,j)  = cflowin*(ho_abs(j,i) + vwall*tvel)
           vxcwlkm(i,j)  = cflowin*xvel
           vrcwlkm(i,j)  = cflowin*rvel
           rvtcwlkm(i,j) = cflowin*(vwall  +  tvel)*ravg
      end if
c   end of the  i, j  loop
  750 continue
c
c******************************************************************************
c     end of setting the blade surface cooling fluxes.
  700 continue
c
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c
c     find the sum of the cooling flows added up to each j station
c
      sumcwl(1) = 0.0
      do 20 j=2,jm
      cooladd   = 0.0
      do 25 k=1,kmm1
   25 cooladd   = cooladd  + (cflowi1(j,k) + cflowim(j,k))*nblade(j)
      do 30 i=1,imm1
   30 cooladd   = cooladd  + (cflowk1(i,j) + cflowkm(i,j))*nblade(j)
      sumcwl(j) = sumcwl(j-1) + cooladd
   20 continue
c
c     end of cooling flow calculation.
c******************************************************************************
c******************************************************************************
c
      return
      end
c
c**********************************************************************
