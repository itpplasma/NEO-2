
C .coinas
      SUBROUTINE COINas(RT,R0,RC,NCI,NPH,MMAI,NMAI,NMAXI,
     *cbf,L1I,NOPRI1,NOPRI2)
      implicit double precision (a-h,o-z),integer(i-n)
      save
      parameter (kmp=3,klp=6)
      dimension a1c(kmp,klp),a2c(kmp,klp),a1lc(kmp)
      dimension rc(1)
      logical prop,nopri1,nopri2
      COMMON/cgbas/a1c,a2c,a1lc,a101c,a202c,a103c,a204c,a105c,a206c,
     *rtc,cbfc
      COMMON/cgbasi/mpc
      COMMON/cgbasl/prop
      PROP=.FALSE.
      mpc=l1i
      rtc=rt
      cbfc=cbf
      read(17,*) a101c,a202c,a103c,a204c,a105c,a206c
      read(17,*) a1lc
      read(17,*) a1c
      read(17,*) a2c

c      read(17,*)(kM(ML),NM(ML),H(1,ML),ML=1,BMMA)
c      read(17,*)(NAI(ML),AI(1,ML),ML=1,15)

      IF(NOPRI1) return
!      PRINT 200,a101c,a202c,a103c
!      PRINT 200,a204c,a105c,a206c
!      PRINT 200,a1lc
  200 FORMAT(4X,3(1pE13.5,3X))
      DO 210 ML=1,klp
!      PRINT 201,(a1c(mm,ml),mm=1,kmp)
!      PRINT 202,(a2c(mm,ml),mm=1,kmp)
  201 FORMAT(4X,3(1pE13.5,3X),'a1')
  202 FORMAT(4X,3(1pE13.5,3X),'a2')
  210 CONTINUE

      RETURN
      END
C .coinas

      SUBROUTINE GBas_1(RIN,FIN,ZIN,BROFF,BFOFF,BZOFF,
     *BRROF,BRFOF,BRZOF,BFROF,BFFOF,BFZOF,BZROF,BZFOF,BZZOF)
      implicit double precision (a-h,o-z),integer(i-n)
      save
      parameter (kmp=3,klp=6)
      dimension a1c(kmp,klp),a2c(kmp,klp),a1lc(kmp)
      dimension a1(kmp,klp),a2(kmp,klp),a1l0(kmp)
      dimension cn10(kmp),cn20(kmp),
     *an11(kmp),an10(kmp),
     *bn11(kmp),bn10(kmp),
     *an22(kmp),an21(kmp),an20(kmp),
     *bn22(kmp),bn21(kmp),bn20(kmp)
      dimension cd10(kmp),cd20(kmp),cd30(kmp),
     *ad11(kmp),ad10(kmp),
     *bd11(kmp),bd10(kmp),
     *ad22(kmp),ad21(kmp),ad20(kmp),
     *bd22(kmp),bd21(kmp),bd20(kmp),
     *ad33(kmp),ad32(kmp),ad31(kmp),ad30(kmp),
     *bd33(kmp),bd32(kmp),bd31(kmp),bd30(kmp)
      LOGICAL PROP
      COMMON/cgbas/a1c,a2c,a1lc,a101c,a202c,a103c,a204c,a105c,a206c,
     *rtc,cbfc
      COMMON/cgbasi/mpc
      COMMON/cgbasl/prop

      IF(PROP) GO TO 30
      PROP=.TRUE.

      mp=mpc
      rt=rtc
      cbf=cbfc
      o2=1d0/2
      o3=1d0/3
      o4=1d0/4
      o6=1d0/6
      o8=1d0/8
      o9=1d0/9
      o10=1d0/10
      o15=1d0/15
      o16=1d0/16
      o4m3=3d0/4
      o8m3=3d0/8
      o9m19=19d0/9

      a101=a101c
      a202=a202c
      a103=a103c
      a204=a204c
      a105=a105c
      a206=a206c

      do 20 km=1,kmp
      a1l0(km)=a1lc(km)
      mre=mp*km
      om=1d0/mre
      omp1=1d0/(mre+1)
      omm1=1d0/(mre-1)
      omp2=1d0/(mre+2)
      omm2=1d0/(mre-2)
      omp3=1d0/(mre+3)
      omm3=1d0/(mre-3)

      cn10(km)=o8*om*om
      an11(km)=-mre*omp1
      an10(km)=mre*omm1
      bn11(km)=-mre*omm1
      bn10(km)=mre*omp1

      cn20(km)=cn10(km)*o4*om
      an22(km)=-an11(km)*mre*o2*omp2
      an21(km)=-an10(km)*mre*omp1
      an20(km)=an10(km)+an11(km)*o2-bn10(km)*omm1-bn11(km)*omm2

      bn22(km)=bn11(km)*mre*o2*omm2
      bn21(km)=bn10(km)*mre*omm1
      bn20(km)=-(bn10(km)+bn11(km)*o2+an10(km)*omp1+an11(km)*omp2)

      cd10(km)=o8*om
      ad11(km)=-mre*omp1
      ad10(km)=(mre-2)*omm1
      bd11(km)=mre*omm1
      bd10(km)=-(mre+2)*omp1

      cd20(km)=cd10(km)*o4*om
      ad22(km)=-ad11(km)*mre*o2*omp2
      ad21(km)=-ad10(km)*mre*omp1
      ad20(km)=ad10(km)+ad11(km)*o2-bd10(km)*omm1-bd11(km)*omm2

      bd22(km)=bd11(km)*mre*o2*omm2
      bd21(km)=bd10(km)*mre*omm1
      bd20(km)=-(bd10(km)+bd11(km)*o2+ad10(km)*omp1+ad11(km)*omp2)

      cd30(km)=cd20(km)*o4*om
      ad33(km)=-ad22(km)*mre*o3*omp3
      ad32(km)=-ad21(km)*mre*o2*omp2
      ad31(km)=-ad20(km)*mre*omp1
      ad30(km)=ad20(km)+ad21(km)*o2+ad22(km)*o3-
     *bd20(km)*omm1-bd21(km)*omm2-bd22(km)*omm3

      bd33(km)=bd22(km)*mre*o3*omm3
      bd32(km)=bd21(km)*mre*o2*omm2
      bd31(km)=bd20(km)*mre*omm1
      bd30(km)=-(bd20(km)+bd21(km)*o2+bd22(km)*o3+
     *ad20(km)*omp1+ad21(km)*omp2+ad22(km)*omp3)
      do 20 kl=1,klp
      a1(km,kl)=a1c(km,kl)
      a2(km,kl)=a2c(km,kl)

   20 continue
   30 rotn=rin/rt
      zotn=zin/rt
      fi=fin
      or=1d0/rotn
      r2=rotn*rotn
      r4=r2*r2
      r6=r2*r4
      z2=zotn*zotn
      z2d2=z2*o2

      cn00=dlog(rotn)
      cn01=-o4*((r2+1d0)*cn00+1d0-r2)
      cn02=o16*((o4*r4+r2+o4)*cn00+o8m3*(1d0-r4))
      cd00=1d0
      cd01=o2*cn00+o4*(1d0-r2)
      cd02=o2*cn01+o16*(cn00+o4*r4-r2+o4m3)
c      cd03=o2*cn02+o16*(cn01+o6*cn00-o16*(o9*r6-r4+3d0*r2-o9m19))

      cn00r=or
      cn01r=-o4*(2d0*r2*cn00+1d0-r2)*or
      cn02r=o16*((2d0+r2)*r2*cn00+o4+r2-5d0*o4*r4)*or

      cd00r=0d0
      cd01r=o2*(cn00r-rotn)
      cd02r=o2*cn01r+((r2-2d0)*rotn+cn00r)*o16

      cn00rr=-or*or
      cn01rr=-cn00-cn01r*or
      cn02rr=-cn01-cn02r*or

      cd00rr=0d0
      cd01rr=-cd00-cd01r*or
      cd02rr=-cd01-cd02r*or

      vn00=cn00
      vd00=cd00

      vn01=zotn*cn00
      vd01=zotn*cd00

      vn02=z2d2*cn00+cn01
      vd02=z2d2*cd00+cd01

      vn03=(o3*z2d2*cn00+cn01)*zotn
      vd03=(o3*z2d2*cd00+cd01)*zotn

      vn04=(o6*z2d2*cn00+cn01)*z2d2+cn02
      vd04=(o6*z2d2*cd00+cd01)*z2d2+cd02

      vn01z=vn00
      vd01z=vd00

      vn03z=vn02
      vd03z=vd02

      vn05z=vn04
      vd05z=vd04

      vn00r=cn00r
      vd00r=cd00r

      vn01r=zotn*cn00r
      vd01r=zotn*cd00r

      vn02r=z2d2*cn00r+cn01r
      vd02r=z2d2*cd00r+cd01r

      vn03r=(o3*z2d2*cn00r+cn01r)*zotn
      vd03r=(o3*z2d2*cd00r+cd01r)*zotn

      vn04r=(o6*z2d2*cn00r+cn01r)*z2d2+cn02r
      vd04r=(o6*z2d2*cd00r+cd01r)*z2d2+cd02r

      vn05r=((o10*z2d2*cn00r+cn01r)*z2d2*o3+cn02r)*zotn
      vd05r=((o10*z2d2*cd00r+cd01r)*z2d2*o3+cd02r)*zotn

      vn01rr=zotn*cn00rr
      vd01rr=zotn*cd00rr

      vn03rr=(o3*z2d2*cn00rr+cn01rr)*zotn
      vd03rr=(o3*z2d2*cd00rr+cd01rr)*zotn

      vn05rr=((o10*z2d2*cn00rr+cn01rr)*z2d2*o3+cn02rr)*zotn
      vd05rr=((o10*z2d2*cd00rr+cd01rr)*z2d2*o3+cd02rr)*zotn

      vn01zz=0d0
      vd01zz=0d0

      vn03zz=vn01
      vd03zz=vd01

      vn05zz=vn03
      vd05zz=vd03

      vn01zr=vn00r
      vd01zr=vd00r

      vn03zr=vn02r
      vd03zr=vd02r

      vn05zr=vn04r
      vd05zr=vd04r

      v0z=a206*vn05z+a105*vd05z+a204*vn03z+a103*vd03z+
     *a202*vn01z+a101*vd01z

      v0r=a206*vn05r+a105*vd05r+a204*vn03r+a103*vd03r+
     *a202*vn01r+a101*vd01r

      v0rr=a206*vn05rr+a105*vd05rr+a204*vn03rr+a103*vd03rr+
     *a202*vn01rr+a101*vd01rr

      v0zz=a206*vn05zz+a105*vd05zz+a204*vn03zz+a103*vd03zz+
     *a202*vn01zz+a101*vd01zz

      v0zr=a206*vn05zr+a105*vd05zr+a204*vn03zr+a103*vd03zr+
     *a202*vn01zr+a101*vd01zr

      vrb=v0r
      vfb=0d0
      vzb=v0z
      vrrb=v0rr
      vzzb=v0zz
      vzrb=v0zr
      vrfb=0d0
      vzfb=0d0
      vffb=0d0

      do 40 km=1,kmp
      mre=km*mp
      mre2=mre*mre
      fim=mre*fi
      sinmf=dsin(fim)
      cosmf=dcos(fim)
      rm=rotn**mre
      orm=1d0/rm
      cnm0=o2*(rm-orm)/mre
      cnm1=cn10(km)*(rm*(an11(km)*r2+an10(km))+orm*(bn11(km)*r2+
     *bn10(km)))
      cnm2=cn20(km)*(rm*(an22(km)*r4+an21(km)*r2+an20(km))+
     *orm*(bn22(km)*r4+bn21(km)*r2+bn20(km)))

      cdm0=o2*(rm+orm)
      cdm1=cd10(km)*(rm*(ad11(km)*r2+ad10(km))+orm*(bd11(km)*r2+
     *bd10(km)))
      cdm2=cd20(km)*(rm*(ad22(km)*r4+ad21(km)*r2+ad20(km))+
     *orm*(bd22(km)*r4+bd21(km)*r2+bd20(km)))
      cdm3=cd30(km)*(rm*(ad33(km)*r6+ad32(km)*r4+ad31(km)*r2+
     *ad30(km))+orm*(bd33(km)*r6+bd32(km)*r4+bd31(km)*r2+bd30(km)))

      an22r=an22(km)*(mre+4)
      an21r=an21(km)*(mre+2)
      an20r=an20(km)*mre

      an11r=an11(km)*(mre+2)
      an10r=an10(km)*mre

      bn22r=-bn22(km)*(mre-4)
      bn21r=-bn21(km)*(mre-2)
      bn20r=-bn20(km)*mre

      bn11r=-bn11(km)*(mre-2)
      bn10r=-bn10(km)*mre

      ad33r=ad33(km)*(mre+6)
      ad32r=ad32(km)*(mre+4)
      ad31r=ad31(km)*(mre+2)
      ad30r=ad30(km)*mre

      ad22r=ad22(km)*(mre+4)
      ad21r=ad21(km)*(mre+2)
      ad20r=ad20(km)*mre

      ad11r=ad11(km)*(mre+2)
      ad10r=ad10(km)*mre

      bd33r=-bd33(km)*(mre-6)
      bd32r=-bd32(km)*(mre-4)
      bd31r=-bd31(km)*(mre-2)
      bd30r=-bd30(km)*mre

      bd22r=-bd22(km)*(mre-4)
      bd21r=-bd21(km)*(mre-2)
      bd20r=-bd20(km)*mre

      bd11r=-bd11(km)*(mre-2)
      bd10r=-bd10(km)*mre

      cnm0r=cdm0*or
      cnm1r=cn10(km)*(rm*(an11r*r2+an10r)+orm*(bn11r*r2+
     *bn10r))*or
      cnm2r=cn20(km)*(rm*(an22r*r4+an21r*r2+an20r)+
     *orm*(bn22r*r4+bn21r*r2+bn20r))*or

      cdm0r=mre*mre*cnm0*or
      cdm1r=cd10(km)*(rm*(ad11r*r2+ad10r)+orm*(bd11r*r2+
     *bd10r))*or
      cdm2r=cd20(km)*(rm*(ad22r*r4+ad21r*r2+ad20r)+
     *orm*(bd22r*r4+bd21r*r2+bd20r))*or
      cdm3r=cd30(km)*(rm*(ad33r*r6+ad32r*r4+ad31r*r2+
     *ad30r)+orm*(bd33r*r6+bd32r*r4+bd31r*r2+bd30r))*or

      cnm0rr=(mre2*cnm0*or-cnm0r)*or
      cdm0rr=(mre2*cdm0*or-cdm0r)*or

      cnm1rr=(mre2*cnm1*or-cnm0*rotn-cnm1r)*or
      cdm1rr=(mre2*cdm1*or-cdm0*rotn-cdm1r)*or

      cnm2rr=(mre2*cnm2*or-cnm1*rotn-cnm2r)*or
      cdm2rr=(mre2*cdm2*or-cdm1*rotn-cdm2r)*or

      cdm3rr=(mre2*cdm3*or-cdm2*rotn-cdm3r)*or

      vnm0=cnm0
      vnm1=zotn*cnm0
      vnm2=z2d2*cnm0+cnm1
      vnm3=(o3*z2d2*cnm0+cnm1)*zotn
      vnm4=(o6*z2d2*cnm0+cnm1)*z2d2+cnm2
      vnm5=((o10*z2d2*cnm0+cnm1)*z2d2*o3+cnm2)*zotn

      vdm0=cdm0
      vdm1=zotn*cdm0
      vdm2=z2d2*cdm0+cdm1
      vdm3=(o3*z2d2*cdm0+cdm1)*zotn
      vdm4=(o6*z2d2*cdm0+cdm1)*z2d2+cdm2
      vdm5=((o10*z2d2*cdm0+cdm1)*z2d2*o3+cdm2)*zotn
      vdm6=((o15*z2d2*cdm0+cdm1)*z2d2*o6+cdm2)*z2d2+cdm3

      vnm0r=cnm0r
      vnm1r=zotn*cnm0r
      vnm2r=z2d2*cnm0r+cnm1r
      vnm3r=(o3*z2d2*cnm0r+cnm1r)*zotn
      vnm4r=(o6*z2d2*cnm0r+cnm1r)*z2d2+cnm2r
      vnm5r=((o10*z2d2*cnm0r+cnm1r)*z2d2*o3+cnm2r)*zotn

      vdm0r=cdm0r
      vdm1r=zotn*cdm0r
      vdm2r=z2d2*cdm0r+cdm1r
      vdm3r=(o3*z2d2*cdm0r+cdm1r)*zotn
      vdm4r=(o6*z2d2*cdm0r+cdm1r)*z2d2+cdm2r
      vdm5r=((o10*z2d2*cdm0r+cdm1r)*z2d2*o3+cdm2r)*zotn
      vdm6r=((o15*z2d2*cdm0r+cdm1r)*z2d2*o6+cdm2r)*z2d2+cdm3r

      vnm0rr=cnm0rr
      vnm1rr=zotn*cnm0rr
      vnm2rr=z2d2*cnm0rr+cnm1rr
      vnm3rr=(o3*z2d2*cnm0rr+cnm1rr)*zotn
      vnm4rr=(o6*z2d2*cnm0rr+cnm1rr)*z2d2+cnm2rr
      vnm5rr=((o10*z2d2*cnm0rr+cnm1rr)*z2d2*o3+cnm2rr)*zotn

      vdm0rr=cdm0rr
      vdm1rr=zotn*cdm0rr
      vdm2rr=z2d2*cdm0rr+cdm1rr
      vdm3rr=(o3*z2d2*cdm0rr+cdm1rr)*zotn
      vdm4rr=(o6*z2d2*cdm0rr+cdm1rr)*z2d2+cdm2rr
      vdm5rr=((o10*z2d2*cdm0rr+cdm1rr)*z2d2*o3+cdm2rr)*zotn
      vdm6rr=((o15*z2d2*cdm0rr+cdm1rr)*z2d2*o6+cdm2rr)*z2d2+cdm3rr

      vnm1z=vnm0
      vnm2z=vnm1
      vnm3z=vnm2
      vnm4z=vnm3
      vnm5z=vnm4

      vdm1z=vdm0
      vdm2z=vdm1
      vdm3z=vdm2
      vdm4z=vdm3
      vdm5z=vdm4
      vdm6z=vdm5

      vnm1zz=0d0
      vnm2zz=vnm0
      vnm3zz=vnm1
      vnm4zz=vnm2
      vnm5zz=vnm3

      vdm1zz=0d0
      vdm2zz=vdm0
      vdm3zz=vdm1
      vdm4zz=vdm2
      vdm5zz=vdm3
      vdm6zz=vdm4

      vnm1zr=vnm0r
      vnm2zr=vnm1r
      vnm3zr=vnm2r
      vnm4zr=vnm3r
      vnm5zr=vnm4r

      vdm1zr=vdm0r
      vdm2zr=vdm1r
      vdm3zr=vdm2r
      vdm4zr=vdm3r
      vdm5zr=vdm4r
      vdm6zr=vdm5r

      vmc=a2(km,6)*vnm5+a1(km,5)*vdm5+a2(km,4)*vnm3+
     *a1(km,3)*vdm3+a2(km,2)*vnm1+a1(km,1)*vdm1
      vms=a1(km,6)*vdm6+a2(km,5)*vnm4+a1(km,4)*vdm4+
     *a2(km,3)*vnm2+a1(km,2)*vdm2+a2(km,1)*vnm0+a1l0(km)*vdm0

      vmcz=a2(km,6)*vnm5z+a1(km,5)*vdm5z+a2(km,4)*vnm3z+
     *a1(km,3)*vdm3z+a2(km,2)*vnm1z+a1(km,1)*vdm1z
      vmsz=a1(km,6)*vdm6z+a2(km,5)*vnm4z+a1(km,4)*vdm4z+
     *a2(km,3)*vnm2z+a1(km,2)*vdm2z

      vmcr=a2(km,6)*vnm5r+a1(km,5)*vdm5r+a2(km,4)*vnm3r+
     *a1(km,3)*vdm3r+a2(km,2)*vnm1r+a1(km,1)*vdm1r
      vmsr=a1(km,6)*vdm6r+a2(km,5)*vnm4r+a1(km,4)*vdm4r+
     *a2(km,3)*vnm2r+a1(km,2)*vdm2r+a2(km,1)*vnm0r+a1l0(km)*vdm0r

      vmcrr=a2(km,6)*vnm5rr+a1(km,5)*vdm5rr+a2(km,4)*vnm3rr+
     *a1(km,3)*vdm3rr+a2(km,2)*vnm1rr+a1(km,1)*vdm1rr
      vmsrr=a1(km,6)*vdm6rr+a2(km,5)*vnm4rr+a1(km,4)*vdm4rr+
     *a2(km,3)*vnm2rr+a1(km,2)*vdm2rr+a2(km,1)*vnm0rr+a1l0(km)*vdm0rr

      vmczz=a2(km,6)*vnm5zz+a1(km,5)*vdm5zz+a2(km,4)*vnm3zz+
     *a1(km,3)*vdm3zz+a2(km,2)*vnm1zz+a1(km,1)*vdm1zz
      vmszz=a1(km,6)*vdm6zz+a2(km,5)*vnm4zz+a1(km,4)*vdm4zz+
     *a2(km,3)*vnm2zz+a1(km,2)*vdm2zz

      vmczr=a2(km,6)*vnm5zr+a1(km,5)*vdm5zr+a2(km,4)*vnm3zr+
     *a1(km,3)*vdm3zr+a2(km,2)*vnm1zr+a1(km,1)*vdm1zr
      vmszr=a1(km,6)*vdm6zr+a2(km,5)*vnm4zr+a1(km,4)*vdm4zr+
     *a2(km,3)*vnm2zr+a1(km,2)*vdm2zr

      vmf=mre*(vms*cosmf-vmc*sinmf)
      vmz=vmcz*cosmf+vmsz*sinmf
      vmr=vmcr*cosmf+vmsr*sinmf

      vmrr=vmcrr*cosmf+vmsrr*sinmf
      vmrf=mre*(vmsr*cosmf-vmcr*sinmf)
      vmrz=vmczr*cosmf+vmszr*sinmf

      vmzr=vmrz
      vmzf=mre*(vmsz*cosmf-vmcz*sinmf)
      vmzz=vmczz*cosmf+vmszz*sinmf

      vmfr=vmrf
      vmff=-mre2*(vmc*cosmf+vms*sinmf)
      vmfz=vmzf

      vrb=vrb+vmr
      vfb=vfb+vmf
      vzb=vzb+vmz

      vrrb=vrrb+vmrr
      vrfb=vrfb+vmrf
      vzrb=vzrb+vmzr
      vzzb=vzzb+vmzz
      vzfb=vzfb+vmzf
      vffb=vffb+vmff

   40 continue

      BROFF=vrb
      BZOFF=vzb
      bcf=(vfb+cbf)*or
      BFOFF=bcf
      BRROF=vrrb/rt
      BRZOF=vzrb/rt
      BRFOF=vrfb

      BZROF=vzrb/rt
      BZZOF=vzzb/rt
      BZFOF=vzfb
      BFROF=(vrfb-bcf)*or/rt
      BFZOF=vzfb*or/rt
      BFFOF=vffb*or
      RETURN
      END
C .GBRFZY


      SUBROUTINE stevvo_1(RT0,R0i,L1i,cbfi,BY0i,bf0)
cccccccc input of W7-AS parameters ccccccccccccccccc
      INTEGER BMMA,BNC
      LOGICAL NOPRI1,NOPRI2
      double precision RC(6),cbf,rt,r0,by0
      double precision cbfi,rt0,r0i,by0i,bf0
      BNC=1
      NOPRI1=.false.
      NOPRI2=.FALSE.

      print *,'this is W7AS'  
      open(17,form='FORMATTED',file='pcas1.dat',
     ,     status='old')
      READ(17,340)nopri1,nopri2
  340 FORMAT(7X,L5,6X,7X,L5)

      READ(17,120)NMAX,MMA,NMA,RT,R0,L1,NC
  120 FORMAT(5X,I2,2X,4X,I3,2X,5X,I2,2X,2(3X,E10.4,2X),2(3X,I3,2X))
      IF(NC.GT.BNC) NC=BNC
!      PRINT 130,NMAX,BMMA,MMA,NMA,L1,RT,R0,NC
  130 FORMAT(//5HNMAX=,I3,4X,4X,5HBMMA=,I3,4X//
     *4HMMA=,I3,4X,4HNMA=,I3,4X,3HL1=,I3,4X,
     *3HRT=,E14.7,4X,3HR0=,E14.7,4X,3HNC=,I3/)
      READ(17,420)cbf,RC
  420 FORMAT(5X,E13.7/3(3X,E13.7)/3(3X,E13.7))
      READ(17,180)BY0
  180 FORMAT(4X,E14.7)
!      PRINT 250,BY0,cbf
  250 FORMAT(/4HBY0=,E14.7,3X,4Hcbf=,E14.7/)
!      PRINT 480,RC
  480 FORMAT(/3HRC=,6(E14.7,2X))
      CALL COINas(RT,R0,RC,NC,BMMA,MMA,NMA,NMAX,
     *cbf,L1,NOPRI1,NOPRI2)
      close(17)
      RT0=rt
!      R0i=r0
      R0i=37.d0
      L1i=l1
      cbfi=cbf
      BY0i=by0
      bf0=cbf
      return
      END


