c ======================================================================
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,
     &  ddsddt,drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,
     &  dpred,cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,
     &  pnewdt,celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c ======================================================================
c UMAT: simple linear-elastic material
c ----------------------------------------------------------------------
      implicit none
      character(80) cmname
      integer(4) ntens,ndi, nshr, nstatv, nprops, noel, npt,
     & layer, kspt, kstep, kinc
      real(8) stress(ntens), statev(nstatv), ddsdde(ntens,ntens),
     &  ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     &  time(2), predef(1), dpred(1), props(nprops), coords(3),
     &  drot(3,3), dfgrd0(3,3), dfgrd1(3,3), sse, spd, scd, rpl,
     &  drpldt, dtime, temp, dtemp, pnewdt, celent

      integer i, j, k
      real(8) E, possion
      
c elastic properties
      E  = props(1)      ! young's modulus (Pa)
      possion = props(2) ! possion's ratio

      if (ntens.eq.3) then
        ddsdde = 0.0d0
        ddsdde(1,1) = 1.0d0
        ddsdde(2,2) = 1.0d0
        ddsdde(1,2) = possion/(1.0d0-possion)
        ddsdde(2,1) = possion/(1.0d0-possion)
        ddsdde(3,3) = (1.0d0-2.0d0*possion)/(2.0d0*(1.0d0-possion))
        ddsdde = ddsdde*E*(1.0d0-possion)/
     &        ((1.0d0+possion)*(1.0d0-2.0d0*possion)) 
      elseif (ntens.eq.6) then
        ddsdde = 0.0d0
        ddsdde(1,1) = 1.0d0
        ddsdde(2,2) = 1.0d0
        ddsdde(3,3) = 1.0d0
        ddsdde(1,2) = possion/(1-possion)
        ddsdde(2,1) = possion/(1-possion)
        ddsdde(2,3) = possion/(1-possion)
        ddsdde(3,2) = possion/(1-possion)
        ddsdde(1,3) = possion/(1-possion)
        ddsdde(3,1) = possion/(1-possion)
        ddsdde(4,4) = (1-2*possion)/(2*(1-possion))
        ddsdde(5,5) = (1-2*possion)/(2*(1-possion))
        ddsdde(6,6) = (1-2*possion)/(2*(1-possion))
        ddsdde = ddsdde*E*(1-possion)/((1+possion)*(1-2*possion))
      endif

      do i=1,ntens
        do j=1,ntens
          stress(i) = stress(i) + ddsdde(i,j)*dstran(j)
        enddo           
      enddo      

      return
      end subroutine umat