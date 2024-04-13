c **********************************************************************
c written by Zexu Fan, Tongji Univeristy
c email: 2110049@tongji.edu.cn
c **********************************************************************

        subroutine UEL(rhs,amatrx,svars,energy,ndofel,nrhs,
     &     nsvars,props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,
     &     dtime,kstep,kinc,jelem,params,ndload,jdltyp,adlmag,
     &     predef,npredf,lflags,mlvarx,ddlmag,mdload,pnewdt,
     &     jprops,njprop,period)

      implicit none

      integer(4) ndofel, nrhs, nsvars, nprops, njprop, mcrd, nnode,
     &           jtype, kstep, kinc, jelem, ndload, npredf, mlvarx,
     &           mdload, jprops(njprop), lflags(*), jdltyp(mdload,*)
      real(8) dtime, pnewdt, period
      real(8) rhs(mlvarx,*), amatrx(ndofel,ndofel), svars(nsvars),
     &        energy(8), props(nprops), coords(mcrd,nnode), u(ndofel),
     &        du(mlvarx,*), v(ndofel), a(ndofel), time(2), params(*),
     &        adlmag(mdload,*), ddlmag(mdload,*), predef(2,npredf,nnode)

      character(80) cmname
      integer(4) ntens, ndi, nshr, nstatv, nprops_umat, npt, layer, kspt
      real(8) stress(6), ddsdde(6,6), ddsddt(6), drplde(6), stran(6), 
     &  dstran(6), predefUmat(1), dpred(1), coordsUmat(3), drot(3,3), 
     &  dfgrd0(3,3), dfgrd1(3,3), sse, spd, scd, rpl, drpldt, temp, 
     &  dtemp, celent
      real(8), allocatable :: statev(:), props_umat(:)
      
      integer i, j, nprops_uel
      real(8) Nu(2,8), Bu(3,8), m(3,1), matStiff0(3,3), matStiff(3,3),
     &        coordVertex(2,4), coordGussian(2,4), invmat_kbar(2,2), 
     &        matMBase(8,8), matCBase(8,8), matkBar(2,2), 
     &        matD1uu(8,8), matD4uu(8,8), matK1uu(8,8), matK4uu(8,8), 
     &        matK1uw(8,8), matK4uw(8,8), matK1ww(8,8), matK4ww(8,8), 
     &        matKK(16,16), matKK0(16,16), matMM(16,16), matCC(16,16),
     &        dhdx, dhdy, Fu_u(3,2), Lu(2,8), 
     &        h(4,1), b1(4,1), b2(4,1), vecGamma(4,1), vecUu(8,1), 
     &        vecUw(8,1), vecU(16,1), vecV(16,1), vecA(16,1), 
     &        resid1(16,1), resid0(16,1), resForce(16,1), vecRhs(16,1), 
     &        vecbAcc(2,1), strain0(6), tmpScal1, tmpScal2, 
     &        tmpMat3m8(3,8), tmpMat1m8(1,8), tmpMat8m2(8,2), 
     &        tmpMat3m2(3,2), tmpMat2m8(2,8), tmpMat4m1(4,1), 
     &        tmpMat1m1(1,1), tmpMat3m1(3,1)
      real(8) rhoS, rhoF, poro, e0, aa, Qb, Kf, kBar_x, kBar_y, bAcc_x, 
     &        bAcc_y, pp, void
      real(8) xi, eta, detJacobian, area, alpha, beta, gamma,
     &        dadu, dvdu

      nprops_uel = 9    
      rhoS    = props(1)
      rhoF    = props(2)
      e0      = props(3)
      aa      = props(4)
      Kf      = props(5)
      kBar_x  = props(6)
      kBar_y  = props(7)
      bAcc_x  = props(8)
      bAcc_y  = props(9)
      
      nprops_umat = 2 
      nstatv      = 5 
      allocate(props_umat(nprops_umat))
      do i=1,nprops_umat
        props_umat(i) = props(nprops_uel+i)
      enddo

      ntens  = 6
      ndi    = 3
      nshr   = 3

      drot   = 0.0d0
      dfgrd0 = 0.0d0
      dfgrd1 = 0.0d0
      do i=1,3
        drot(i,i) = 1.0d0
        dfgrd0(i,i) = 1.0d0
        dfgrd1(i,i) = 1.0d0
      enddo

      allocate(statev(nstatv))

      matKK     = 0.0d0
      matKK0    = 0.0d0
      matMM     = 0.0d0
      matCC     = 0.0d0

      matD1uu   = 0.0d0
      matD4uu   = 0.0d0
      matK1uu   = 0.0d0
      matK4uu   = 0.0d0
      matK1uw   = 0.0d0
      matK4uw   = 0.0d0
      matK1ww   = 0.0d0
      matK4ww   = 0.0d0

      matMBase  = 0.0d0
      matCBase  = 0.0d0

      resForce  = 0.0d0
      m(:,1)    = (/1.0d0, 1.0d0, 0.0d0/)

      vecU(:,1)  = u
      vecV(:,1)  = v
      vecA(:,1)  = a
      vecUu(:,1) = u(1:8)
      vecUw(:,1) = u(9:16)

      matkBar    = 0.0d0
      invmat_kbar = 0.0d0
      matkBar(1,1) = kBar_x
      matkBar(2,2) = kBar_y
      invmat_kbar(1,1) = 1.0d0/kBar_x
      invmat_kbar(2,2) = 1.0d0/kBar_y

      vecbAcc(:,1) = (/bAcc_x, bAcc_y/)   

      coordVertex(:,1) = (/-1.0d0, -1.0d0/)
      coordVertex(:,2) = (/ 1.0d0, -1.0d0/)
      coordVertex(:,3) = (/ 1.0d0,  1.0d0/)
      coordVertex(:,4) = (/-1.0d0,  1.0d0/)

      coordGussian = coordVertex * 0.5773502692

      area=0.5d0*((coords(1,3)-coords(1,1))*(coords(2,4)-coords(2,2)))
     &    +0.5d0*((coords(1,2)-coords(1,4))*(coords(2,3)-coords(2,1)))
      b1(1,1) = 0.5d0/area*(coords(2,2)-coords(2,4))
      b1(2,1) = 0.5d0/area*(coords(2,3)-coords(2,1))
      b1(3,1) = 0.5d0/area*(coords(2,4)-coords(2,2))
      b1(4,1) = 0.5d0/area*(coords(2,1)-coords(2,3))
      b2(1,1) = 0.5d0/area*(coords(1,4)-coords(1,2))
      b2(2,1) = 0.5d0/area*(coords(1,1)-coords(1,3))
      b2(3,1) = 0.5d0/area*(coords(1,2)-coords(1,4))
      b2(4,1) = 0.5d0/area*(coords(1,3)-coords(1,1))

      h(:,1) = (/1.0d0,-1.0d0,1.0d0,-1.0d0/)
      tmpMat4m1(:,1) = coords(1,:)
      tmpMat1m1 = MATMUL(TRANSPOSE(h),tmpMat4m1)
      tmpScal1 = tmpMat1m1(1,1)
      tmpMat4m1(:,1) = coords(2,:)
      tmpMat1m1 = MATMUL(TRANSPOSE(h),tmpMat4m1)
      tmpScal2 = tmpMat1m1(1,1)        
      vecGamma = 0.25d0 * (h - tmpScal1*b1 - tmpScal2*b2)
      Lu = 0.0d0
      do i=1,4
        Lu(:,2*i-1)=(/vecGamma(i,1),0.0d0/)
        Lu(:,2*i-0)=(/0.0d0,vecGamma(i,1)/)
      enddo

      do npt = 1,1

        xi    = 0.0d0
        eta   = 0.0d0

        detJacobian = 0.125d0*
     &     (((coords(1,2)-coords(1,4))*(coords(2,3)-coords(2,1))+
     &       (coords(1,3)-coords(1,1))*(coords(2,4)-coords(2,2)))+
     &   xi*((coords(1,2)-coords(1,1))*(coords(2,3)-coords(2,4))+
     &       (coords(1,3)-coords(1,4))*(coords(2,1)-coords(2,2)))+
     &  eta*((coords(1,1)-coords(1,4))*(coords(2,3)-coords(2,2))+
     &       (coords(1,3)-coords(1,2))*(coords(2,4)-coords(2,1))))

        if (detJacobian.le.0.0d0) then
          write(*,*) "UEL ERROR: detJacobian <= 0.0"
          call XIT
        endif

        Nu = 0.0d0
        do i = 1,2
          do j = 1,4
            Nu(i,2*(j-1)+i) = 0.25d0 * (1.0d0+coordVertex(1,j)*xi) *
     &                                 (1.0d0+coordVertex(2,j)*eta)
          enddo
        enddo

        Bu = 0.0d0
        do i=1,4
          Bu(:,2*i-1) = (/b1(i,1), 0.0d0, b2(i,1)/)
          Bu(:,2*i-0) = (/0.0d0, b2(i,1), b1(i,1)/)
        enddo

        stress     = svars(1:6)
        strain0    = svars(7:12)
        statev     = svars(41:40+nstatv)
        tmpMat3m1  = MATMUL(Bu,vecUu)
        stran      = 0.0d0
        stran(1)   = tmpMat3m1(1,1)
        stran(2)   = tmpMat3m1(2,1)
        stran(4)   = tmpMat3m1(3,1)
        dstran     = stran - strain0

        call umat(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,
     &  predefUmat,dpred,cmname,ndi,nshr,ntens,nstatv,props_umat,
     &  nprops_umat,coordsUmat,drot,pnewdt,celent,dfgrd0,dfgrd1,
     &  jelem,npt,layer,kspt,kstep,kinc)

        void  = svars(14)
        void  = void + (1+void)*(dstran(1)+dstran(2))
        poro = e0 / (1+e0)
        Qb   = Kf / poro

        tmpMat3m8 = Bu
        tmpMat1m8 = MATMUL(TRANSPOSE(m),tmpMat3m8)
        tmpMat1m1 = MATMUL(tmpMat1m8,vecUw) + aa*MATMUL(tmpMat1m8,vecUu)               
        pp = Qb*tmpMat1m1(1,1)

        svars(1:6)          = stress
        svars(7:12)         = stran
        svars(13)           = pp
        svars(14)           = void
        svars(41:40+nstatv) = statev

        matStiff(1,1) = ddsdde(1,1)
        matStiff(2,1) = ddsdde(2,1)
        matStiff(3,1) = ddsdde(4,1)
        matStiff(1,2) = ddsdde(1,2)
        matStiff(2,2) = ddsdde(2,2)
        matStiff(3,2) = ddsdde(4,2)
        matStiff(1,3) = ddsdde(1,4)
        matStiff(2,3) = ddsdde(2,4)
        matStiff(3,3) = ddsdde(4,4)

        matStiff0 = matStiff

        tmpMat3m8 = MATMUL(matStiff,Bu) 
        matD1uu = matD1uu + MATMUL(TRANSPOSE(Bu),tmpMat3m8)
     &                    *detJacobian*4.0d0

        tmpMat1m8 = MATMUL(TRANSPOSE(m),Bu)
        matK1ww = matK1ww + MATMUL(TRANSPOSE(tmpMat1m8),tmpMat1m8)
     &                    *detJacobian*4.0d0

        matK1uu = matK1ww

        matK1uw = matK1ww

        tmpMat3m1(1,1) = stress(1)
        tmpMat3m1(2,1) = stress(2)
        tmpMat3m1(3,1) = stress(4)
        resForce(1:8,:) = resForce(1:8,:)+ 
     &             MATMUL(TRANSPOSE(Bu),tmpMat3m1)*detJacobian*4.0d0

      enddo     ! end of the loop of 1*1 integration points

      do npt = 1,4

        xi    = coordGussian(1,npt)
        eta   = coordGussian(2,npt)

        detJacobian = 0.125d0*
     &     (((coords(1,2)-coords(1,4))*(coords(2,3)-coords(2,1))+
     &       (coords(1,3)-coords(1,1))*(coords(2,4)-coords(2,2)))+
     &   xi*((coords(1,2)-coords(1,1))*(coords(2,3)-coords(2,4))+
     &       (coords(1,3)-coords(1,4))*(coords(2,1)-coords(2,2)))+
     &  eta*((coords(1,1)-coords(1,4))*(coords(2,3)-coords(2,2))+
     &       (coords(1,3)-coords(1,2))*(coords(2,4)-coords(2,1))))

        if (detJacobian.le.0.0d0) then
           write(*,*) "UEL ERROR: detJacobian <= 0.0"
           call XIT
        endif

        Nu = 0.0d0
        do i = 1,2
          do j = 1,4
            Nu(i,2*(j-1)+i) = 0.25d0 * (1.0d0+coordVertex(1,j)*xi) *
     &                                 (1.0d0+coordVertex(2,j)*eta)
          enddo
        enddo

        dhdx = 0.25d0/detJacobian*
     &   (xi*((+coords(2,1)-coords(2,2)-coords(2,3)+coords(2,4)))+
     &   eta*((-coords(2,1)-coords(2,2)+coords(2,3)+coords(2,4))))
        dhdy = -0.25d0/detJacobian*
     &   (xi*((+coords(1,1)-coords(1,2)-coords(1,3)+coords(1,4)))+
     &   eta*((-coords(1,1)-coords(1,2)+coords(1,3)+coords(1,4))))
        Fu_u = 0.0d0
        Fu_u(:,1) = 1.0d0 * (/dhdx, -dhdx, 0.0d0/)
        Fu_u(:,2) = 1.0d0 * (/-dhdy, dhdy, 0.0d0/)

        tmpMat3m2 = MATMUL(matStiff0,Fu_u)
        tmpMat3m8 = MATMUL(tmpMat3m2,Lu)
        tmpMat2m8 = MATMUL(TRANSPOSE(Fu_u),tmpMat3m8)
        matD4uu = matD4uu + MATMUL(TRANSPOSE(Lu),tmpMat2m8)
     &                    *detJacobian

        tmpMat3m8 = MATMUL(Fu_u,Lu)
        tmpMat1m8 = MATMUL(TRANSPOSE(m),tmpMat3m8)
        matK4uu = matK4uu + MATMUL(TRANSPOSE(tmpMat1m8),tmpMat1m8)
     &                    *detJacobian
     
        matK4ww = matK4uu

        matK4uw = matK4uu

        matMBase = matMBase + MATMUL(TRANSPOSE(Nu),Nu)*detJacobian

        tmpMat8m2 = MATMUL(TRANSPOSE(Nu),invmat_kbar)
        matCBase  = matCBase + MATMUL(tmpMat8m2,Nu)*detJacobian

        resForce(1:8,:)  = resForce(1:8,:)
     &                     - MATMUL(TRANSPOSE(Nu),vecbAcc)         
     &                     * ((1.0d0-poro)*rhoS+poro*rhoF)
     &                     * detJacobian
        resForce(9:16,:) = resForce(9:16,:)
     &                     - MATMUL(TRANSPOSE(Nu),vecbAcc)*(rhoF)                             
     &                     * detJacobian

      enddo     ! end of the loop of 2*2 integration points


      matMM(1:8,1:8)   = matMBase*((1.0d0-poro)*rhoS+poro*rhoF)
      matMM(1:8,9:16)  = matMBase*rhoF
      matMM(9:16,1:8)  = matMBase*rhoF
      matMM(9:16,9:16) = matMBase*rhoF/poro

      matCC(9:16,9:16) = matCBase

      matKK(1:8,1:8)   = matD1uu+matD4uu + 
     &                   (matK1uu+matK4uu)*aa**2*Qb
      matKK(1:8,9:16)  = (matK1uw+matK4uw)*aa*Qb
      matKK(9:16,1:8)  = TRANSPOSE(matKK(1:8,9:16))
      matKK(9:16,9:16) = (matK1ww+matK4ww)*Qb

      matKK0(1:8,1:8)   = matD4uu + 
     &                    (matK1uu+matK4uu)*aa**2*Qb
      matKK0(1:8,9:16)  = (matK1uw+matK4uw)*aa*Qb
      matKK0(9:16,1:8)  = TRANSPOSE(matKK0(1:8,9:16))
      matKK0(9:16,9:16) = (matK1ww+matK4ww)*Qb


      if (lflags(1).eq.11 .or. lflags(1).eq.12) then
        
        alpha=params(1)
        beta=params(2)
        gamma=params(3)
        dadu=1.0d0/(beta*dtime**2)
        dvdu=gamma/(beta*dtime)     

        resid1 = MATMUL(matCC,vecV) + MATMUL(matKK0,vecU) + resForce

        resid0(:,1) = svars(41+nstatv:56+nstatv)
        svars(41+nstatv:56+nstatv) = resid1(:,1)

        if (lflags(3).eq.1) then
          amatrx = matMM*dadu + matCC*(1.0d0+alpha)*dvdu +
     &             matKK*(1.0d0+alpha)
          vecRhs = - MATMUL(matMM,vecA) - (1.0d0+alpha)*resid1
     &             + alpha*resid0
          rhs(1:16,1) = vecRhs(:,1)
        elseif (lflags(3).eq.4) then
          amatrx = matMM
        elseif (lflags(3).eq.6) then
          amatrx = matMM
          vecRhs = -resid1
          rhs(1:16,1) = vecRhs(:,1)
        endif          
       
      elseif (lflags(1).eq.1 .or. lflags(1).eq.2) then

        resid1 = MATMUL(matKK0,vecU) + resForce
        amatrx = matKK
        svars(41+nstatv:56+nstatv) = resid1(:,1)
        vecRhs = -resid1
        rhs(1:16,1) = vecRhs(:,1)

      else

        write(*,*) "Step type not supported by CPE4RUW now"
        call xit

      endif

      deallocate(statev)
      deallocate(props_umat)

      return
      end subroutine UEL