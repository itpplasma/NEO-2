
C .coinhs
      SUBROUTINE COINhs_l(RT,R0,RC,NCI,NPH,MMAI,NMAI,NMAXI,
     *KBZ,L1I,NOPRI1,NOPRI2)

      save /brfzsl/,/brfzse/,/brfzsi/
      parameter (nvez=30,nh=360,nap=15,nnc=6)

      double precision h,ai,ch0,q,q1
      double precision rt,r0,rc,ch02,bzrt

c      REAL MR(6),MR2(6),KBZ
      double precision MR(nnc),MR2(nnc),KBZ
c      double precision MR(6),MR2(6),KBZ
      LOGICAL PROP,NOPRI1,NOPRI2
      INTEGER BMMA
      DIMENSION NAI(nap),RC(1)
      COMMON/BRFZsl/prop
      COMMON/BRFZse/MR,MR2,H(1,nh),AI(1,nap),
     *CH02,BZRT
      COMMON/BRFZsi/KM(nh),NM(nh),
     *MAM(nh),NC,NMAX,MMA,NMA,L1,nai
      PROP=.FALSE.
      MC1=1
      MC2=2
      MC3=3

      MC4=4
      MC5=5
      MC6=6
c      BMMA=NPH
c      MMA=MMAI
c      NMA=NMAI
c      NC=NCI
c      NMAX=NMAXI
      BMMA=nh
      MMA=MMAI
      NMA=NMAI
      NC=NCI
      NMAX=NMAXI
      L1=L1I
      BZRT=-KBZ*L1

c      IF(BMMA.GT.360) BMMA=360
      IF(MMA.GT.BMMA) MMA=BMMA
      IF(NMA.GT.nap) NMA=nap
      IF(NMAX.GT.nvez) NMAX=nvez
      IF(NC.GT.nnc) NC=nnc
      DO 440 MC=1,NC
      MR2(MC)=RT*RT-RC(MC)**2
      MR(MC)=dSQRT(MR2(MC))
  440 CONTINUE

      read(17,*)(kM(ML),NM(ML),H(1,ML),ML=1,BMMA)
      read(17,*)(NAI(ML),AI(1,ML),ML=1,nap)

c     DO 150 ML=1,BMMA

c     READ 140,KM(ML),NM(ML),(H(K,ML),K=1,6)
c 140 FORMAT(I3,I3,1X,6(E12.6))
c 150 CONTINUE
c     DO 170 ML=1,15
c     READ 160,NAI(ML),(AI(K,ML),K=1,6)
c 160 FORMAT(I3,4X,6(E12.6))
c 170 CONTINUE
      IF(NMA.LT.1) GO TO 185
!      PRINT 180,BMMA,MMA,NC,NC,NMAX,(NAI(K),K=1,NMA)
  180 format(/' coinhs: vvedeno',i4,'+5 k; cohs:',
     *i4,'*',i1,' co.hnmq, sloev ',i1,
     *',chlenov gip.ryad.',i3,'garm ANQ ',15(i2,','))
      GO TO 188
  185 continue
!  185 PRINT 187,BMMA,MMA,NC,NC,NMAX
  187 format(/' co-hs: vvedeno',i4,'+5 k; coytc:',
     *i4,'*',i1,' ko.hnmq, sloev',i1,
     *',chl.giper.rya',i3,',garm. ANQ net')

  188 IF(NOPRI1) GO TO 250
!      PRINT 190,MC1,MC2,MC3,MC4,MC5,MC6
  190 format(' coinhs: pech.vved.koeff.',/
     *' HNMQ',6(5x,'mc=',i1,7x),2x,'nm',3x,'k'/)

      DO 210 ML=1,BMMA
!      PRINT 200,H(1,ML),NM(ML),KM(ML)
  200 FORMAT(4X,1(E13.6,3X),I4,I4)
  210 CONTINUE
!      PRINT 220,MC1,MC2,MC3,MC4,MC5,MC6
  220 FORMAT(/' ANQ',6(5X,'MC=',I1,7X),2X,'NA'/)
      DO 240 ML=1,nap
!      PRINT 230,AI(1,ML),NAI(ML)
  230 FORMAT(4X,1(E13.6,3X),I4)

  240 CONTINUE
!      PRINT 245
  245 format(' coinhs: konec pech. vved. koeff')


  250 CH0=RT/R0
      CH02=CH0**2
      IF(MMA.LT.1) GO TO 265
      DO 260 ML=1,MMA
      MAM(ML)=IABS(L1*KM(ML))
      CALL LEGQKB_l(NM(ML),MAM(ML),CH0,Q,Q1,1d-17,NP)
c     DO 255 MC=1,NC

      H(1,ML)=H(1,ML)/Q
c 255 CONTINUE
  260 CONTINUE
  265 IF(NMA.LT.1) GO TO 275
      DO 270 ML=1,NMA
      CALL LEGQKB_l(nAi(ML),1,CH0,Q,Q1,1d-17,NP)
c     DO 267 MC=1,NC
      AI(1,ML)=AI(1,ML)/Q
c 267 CONTINUE
  270 CONTINUE

  275 IF(NOPRI2) RETURN
!      PRINT 280,MC1,MC2,MC3,MC4,MC5,MC6
  280 format(' coinhs: pech. peresch. koeff.',/
     *' HNM',6(5X,'MC=',I1,7X),2X,'NM',3X,'K'/)

      IF(MMA.LT.1) GO TO 305
      DO 300 ML=1,MMA
!      PRINT 290,H(1,ML),NM(ML),KM(ML)
  290 FORMAT(4X,1(E13.6,3X),I4,I4)
  300 CONTINUE
  305 continue
!  305 PRINT 310,MC1,MC2,MC3,MC4,MC5,MC6

  310 FORMAT(/'  AN',6(5X,'MC=',I1,7X),2X,'NA'/)
      IF(NMA.LT.1) GO TO 335
      DO 330 ML=1,NMA
!      PRINT 320,AI(1,ML),NAI(ML)
  320 FORMAT(4X,1(E13.6,3X),I4)
  330 CONTINUE
  335 continue
!  335 PRINT 340
  340 format(' coinhs: konec pech.')


      RETURN
      END
C .coiny

C .GBRFZY
      SUBROUTINE GBhs_l(RIN,FIN,ZIN,BROFF,BFOFF,BZOFF,
     *BRROF,BRFOF,BRZOF,BFROF,BFFOF,BFZOF,BZROF,BZFOF,BZZOF)
c      implicit real*8(a-h,o-z),integer(i-n)
      implicit double precision(a-h,o-z),integer(i-n)
      parameter (nvez=30,nh=360,nap=15,nnc=6)
      parameter (nmal=nvez*nh, nmal1=nvez*nap)

c     REAL*8 MAL(2400),MBL(2400),MAL1(270),MBL1(270),MR(6),MR2(6),
c    *MFI,

      save

c      double precision MAL(10800),MBL(10800),MAL1(450),MBL1(450),
c     *MR(6),MR2(6),mfi,moa(6)
      double precision MAL(nmal),MBL(nmal),MAL1(nmal1),MBL1(nmal1),
     *MR(nnc),MR2(nnc),mfi,moa(nnc)
c     *MFI,
c     *MOA(6)
c      common/comal/mal
c      common/combl/mbl

c      DIMENSION VEZ(30),H(360),AI(15),KM(360),NM(360),MAM(360)
      DIMENSION VEZ(nvez),H(nh),AI(nap),KM(nh),NM(nh),MAM(nh)
      LOGICAL PROP
      COMMON/BRFZsl/prop
      COMMON/BRFZse/RI,R2I,HNM,AN,
     *CH02I,BZRT
      COMMON/BRFZsi/KMI(nh),
     *NMI(nh),MAMI(nh),NCI,NMAXI,MMAI,NMAI,L1I,nai(nap)
c      real RI(6),R2I(6),HNM(1,360),AN(1,15),
c      real RI(6),R2I(6),
c     *CH02I,BZRT

c      dimension HNM(1,360),AN(1,15),RI(6),R2I(6)
      dimension HNM(1,nh),AN(1,nap),RI(nnc),R2I(nnc)

      IF(PROP) GO TO 10
      PROP=.TRUE.
c     NC=NCI
      NC=1
      NMAX=NMAXI
      if(nmax.gt.nvez) nmax=nvez
      NMAX1=NMAX-1
      L1=L1I
      MMA=MMAI
      NMA=NMAI

      if(mma.gt.nh) mma=nh
      if(nma.gt.nap) nma=nap

!      print 25,nvez,nh,nap,nnc,nmal,nmal1,nmax,mma,nma,nc
   25 format(/'   GBhs ',' nvez=',i3,'   nh=',i4,'  nap=',i3,
     *'  nnc=',i2,'  nmal=',i6,'  nmal1=',i4/8x,' nmax=',i3,
     *'  mma=',i4,'  nma=',i3,'   nc=',i2/)

      CH02=CH02I
      CH05=0.5d0*dSQRT(CH02)
      DO 5 MC=1,NC
      MR(MC)=RI(MC)
      MR2(MC)=R2I(MC)
      MOA(MC)=1.0d0/MR(MC)
      DO 2 ML=1,MMA
      H(ML)=HNM(MC,ML)
    2 CONTINUE
      DO 4 ML=1,NMA
      AI(ML)=AN(MC,ML)
    4 CONTINUE
    5 CONTINUE

      IF(MMA.LT.1) GO TO 505
      DO 6 ML=1,MMA
      KM(ML)=KMI(ML)
      NM(ML)=NMI(ML)
      MAM(ML)=MAMI(ML)
    6 CONTINUE
      KML=1
      KMAX=NMAX1
      DO 500 ML=1,MMA
      MM=MAM(ML)

      AL=0.5d0*(NM(ML)+MM+1.5d0)
      BL=0.5d0*(NM(ML)+MM+0.5d0)
      CL=NM(ML)+1
      AL1=AL+0.5d0
      BL1=BL+0.5d0
      MAL(KML)=1.0d0
      MBL(KML)=1.0d0
      DO 490 K=KML,KMAX
      K1=K-KML
      BC1=(K1+1)*CH02

      BC1=1.0d0/((CL+K1)*BC1)
      MAL(K+1)=BC1*MAL(K)*(AL+K1)*(BL+K1)
      MBL(K+1)=BC1*MBL(K)*(AL1+K1)*(BL1+K1)
  490 CONTINUE
      KML=KML+NMAX
      KMAX=KMAX+NMAX
  500 CONTINUE
  505 KML=1
      KMAX=NMAX1
      IF(NMA.LT.1) GO TO 10

      DO 520 ML1=1,NMA
      MM=1
c     ML=ML1-1
      ML=nai(ml1)
      AL=0.5d0*(ML+MM+1.5d0)
      BL=0.5d0*(ML+MM+0.5d0)
      CL=ML+1
      AL1=AL+0.5d0
      BL1=BL+0.5d0
      MAL1(KML)=1.0d0
      MBL1(KML)=1.0d0

      DO 510 K=KML,KMAX
      K1=K-KML
      BC1=(K1+1)*CH02
      BC1=1.0d0/((CL+K1)*BC1)
      MAL1(K+1)=BC1*MAL1(K)*(AL+K1)*(BL+K1)
      MBL1(K+1)=BC1*MBL1(K)*(AL1+K1)*(BL1+K1)
  510 CONTINUE
      KML=KML+NMAX
      KMAX=KMAX+NMAX
  520 CONTINUE

   10 X=RIN
      Y=ZIN
      MFI=FIN*L1

      X2=X*X
      OX=1.0d0/X
      Y2=Y*Y
      BXB=0.0d0
      BYB=0.0d0
      BZB=0.0d0
      BRRB=0.0d0

      BRZB=0.0d0
      BRFB=0.0d0
      BZRB=0.0d0
      BZZB=0.0d0
      BZFB=0.0d0
      BFRB=0.0d0
      BFZB=0.0d0
      BFFB=0.0d0
      DO 90 MC=1,NC
      OA=MOA(MC)

      R=MR(MC)
      R2=MR2(MC)
      BZ1=2.0d0*R*Y
      XPA=X+R
      XMA=X-R
      BX=XPA*XMA
      BY=BX+Y2
      FX=dATAN2(BZ1,BY)
      BY=R2+X2+Y2
      COSV=dCOS(FX)

      SINV=dSIN(FX)
      CH2=BY/((XPA*XPA+Y2)*(XMA*XMA+Y2))
      CH2=CH2*BY
      CH=dSQRT(CH2)
      SH2=CH2-1.0d0
      OSH2=1.0d0/SH2
      SH=dSQRT(SH2)
      CHC=CH-COSV
      SCHC=dSQRT(CHC)
      BKO2=0.5d0*OA

      BKO1=BKO2*Y
      BKO2=BKO2*X
      BKO5=BKO1*BKO1
      BKO3=0.5d0*COSV/CHC-BKO5
      UR=(1.0d0-CH*COSV)*OA
      UZ=-OA*SH*SINV
      VR=UZ
      VZ=-UR
      CHS=CH*SINV
      SHC=SH*COSV

      URR=(CHS*VR-SHC*UR)*OA
      URZ=(CHS*VZ-SHC*UZ)*OA
      UZZ=-URR
      VRR=URZ
      VRZ=UZZ
      VZZ=-URZ
      BKO4=-X*UR*OSH2
      UR2=UR*UR
      UZ2=VR*VR
      URUZ=UR*UZ

      OCH=1.0d0/CH
      BQ=-SH*OCH
      BQ2=CH05*OCH
      SBQ2=dSQRT(BQ2)
      ZL=CH02/CH2
      VEZ(1)=1.0d0
      DO 15 K=2,NMAX
      VEZ(K)=VEZ(K-1)*ZL
   15 CONTINUE
      FIU1=0.0d0

      FIU2=0.0d0
      FIV2=0.0d0
      FIF1=0.0d0
      FIU3=0.0d0
      FIU4=0.0d0
      FIU5=0.0d0
      FIU6=0.0d0
      FIU7=0.0d0
      KMO=0
      DO 30 ML=1,MMA

      T=NM(ML)*FX+KM(ML)*MFI
      MM=MAM(ML)
      NL=NM(ML)
      COST=dCOS(T)
      SINT=dSIN(T)
      F=0.0d0
      F1=0.0d0
      DO 20 K=1,NMAX
      F=F+MAL(K+KMO)*VEZ(K)
      F1=F1+MBL(K+KMO)*VEZ(K)

   20 CONTINUE
      BQ1=BQ**MM*BQ2**NL*SBQ2
      Q2=F*BQ1
      Q=Q2
      Q1=(NL+MM+0.5d0)*F1*BQ1*BQ-MM*Q2/BQ
      BC1=H(ML)*SINT
      BC2=H(ML)*COST
      BC6=BC1*Q
      FIU1=FIU1+BC6
      FIU2=BC1*Q1+FIU2
      BC3=BC2*Q
      FIV2=BC3*NL+FIV2
      FIF1=KM(ML)*BC3+FIF1
      BC4=BC2*Q1
      BC5=BC6*NL
      FIU3=BC5*NL+FIU3
      FIU4=BC6*MM*MM+FIU4
      FIU5=BC4*NL+FIU5
      FIU6=BC4*KM(ML)+FIU6
      FIU7=BC5*KM(ML)+FIU7

      KMO=KMO+NMAX
   30 CONTINUE
      AU1=0.0d0
      AU2=0.0d0
      AU3=0.0d0
      AV2=0.0d0
      AV3=0.0d0
      KMO=0
      IF(NMA.LT.1) GO TO 85
      DO 80 ML1=1,NMA

c     ML=ML1-1
      ML=nai(ml1)
      IF(ML-1) 40,50,60
   40 SINNV=0.0d0
      COSNV=1.0d0
      GO TO 65
   50 SINNV=SINV
      COSNV=COSV
      GO TO 65
   60 T=ML*FX
      COSNV=dCOS(T)

      SINNV=dSIN(T)
   65 F=0.0d0
      F1=0.0d0
      DO 70 K=1,NMAX
      F=F+MAL1(K+KMO)*VEZ(K)
      F1=F1+MBL1(K+KMO)*VEZ(K)
   70 CONTINUE
      BQ1=BQ*BQ2**ML*SBQ2
      Q2=F*BQ1
      Q=Q2

      Q1=(ML+1.5d0)*F1*BQ1*BQ-Q2/BQ
      BC1=AI(ML1)*COSNV
      BC2=AI(ML1)*SINNV*ML
      BC3=BC1*Q
      AU1=AU1+BC3
      AU2=BC1*Q1+AU2
      AU3=BC3*ML*ML+AU3
      AV2=BC2*Q+AV2
      AV3=BC2*Q1+AV3


      KMO=KMO+NMAX

   80 CONTINUE
   85 FIU=(BKO2*FIU1+FIU2)*SCHC
      FIV=(BKO1*FIU1+FIV2)*SCHC
      FIFI=L1*SCHC*FIF1
      FIUU=(BKO4*FIU2+BKO5*FIU1+FIU3+FIU4*OSH2)*SCHC
      FIUV=((FIU2-BKO2*FIU1)*BKO1+FIU5+FIV2*BKO2)*SCHC
      FIVV=(FIU1*BKO3-FIU3+FIV2*2.0d0*BKO1)*SCHC
      FIUF=(FIF1*BKO2+FIU6)*SCHC*L1
      FIVF=(FIF1*BKO1-FIU7)*SCHC*L1
      FIFF=-SCHC*FIU4

      BX=UR*FIU+VR*FIV
      BY=UZ*FIU+VZ*FIV
      BF=FIFI*OX
      BZB=BZB+BF
      BC1=2.0d0*FIUV*URUZ
      BC2=UZ2-UR2
      BRR=FIUU*UR2+BC1+FIVV*UZ2+FIU*URR+FIV*VRR
      BRZ=(FIUU-FIVV)*URUZ+FIUV*BC2+FIU*URZ+FIV*VRZ
      BRF=FIUF*UR+FIVF*VR
      BZR=BRZ

      BZZ=FIUU*UZ2-BC1+FIVV*UR2+FIU*UZZ+FIV*VZZ
      BZF=FIUF*UZ+FIVF*VZ
      BFR=(BRF-BF)*OX
      BFZ=OX*BZF
      BFF=OX*FIFF
      AU=(BKO2*AU1+AU2)*SCHC
      AV=(BKO1*AU1-AV2)*SCHC
      AUU=((OSH2+BKO5)*AU1+AU3+AU2*BKO4)*SCHC
      AUV=-((AU1*BKO2-AU2)*BKO1+AV3+AV2*BKO2)*SCHC
      AVV=(AU1*BKO3-AU3-2.0d0*BKO1*AV2)*SCHC

      BC1=2.0d0*AUV*URUZ
      AR=AU*UR+AV*VR
      AZ=AU*UZ+AV*VZ
      AF=AU1*SCHC
      ARR=AUU*UR2+BC1+AVV*UZ2+AU*URR+AV*VRR
      ARZ=(AUU-AVV)*URUZ+AUV*BC2+AU*URZ+AV*VRZ
      AZZ=AUU*UZ2-BC1+AVV*UR2+AU*UZZ+AV*VZZ
      BRR=BRR-ARZ
      BRZ=BRZ-AZZ
      BZR=(AR-AF*OX)*OX+ARR+BZR

      BZZ=AZ*OX+ARZ+BZZ
      BXB=BXB+BX-AZ
      BYB=BYB+BY+(AR+AF*OX)
      BRRB=BRRB+BRR
      BRZB=BRZB+BRZ
      BRFB=BRFB+BRF
      BZRB=BZRB+BZR
      BZZB=BZZB+BZZ
      BZFB=BZFB+BZF
      BFRB=BFRB+BFR

      BFZB=BFZB+BFZ
      BFFB=BFFB+BFF

   90 CONTINUE
      BROFF=BXB
      BZOFF=BYB
      BC1=BZRT*OX
      BFOFF=BZB+BC1
      BRROF=BRRB
      BRZOF=BRZB
      BRFOF=BRFB

      BZROF=BZRB
      BZZOF=BZZB
      BZFOF=BZFB
      BFROF=BFRB-BC1*OX
      BFZOF=BFZB
      BFFOF=BFFB

      RETURN
      END
C .GBRFZY
C .GBRZY
      SUBROUTINE GBRZd_l(RI,ZI,BRI,BZI,BRRI,BRZI,BZRI,BZZI)
      implicit double precision(a-h,o-z),integer(i-n)
      save
      dimension ARZ(16),ARK(16),ARC(16)
      dimension ARZi(16),ARKi(16),ARCi(16)
c      REAL ARZi(16),ARKi(16),ARCi(16)

      INTEGER SK,SMA
      double precision KAP,KAP2,KAP2R,KAP2Z

      LOGICAL PROP
      DATA PROP/.FALSE./
      IF(PROP) GO TO 20
      PROP=.TRUE.
      EE0=2.0d0*dATAN(1.0d0)
      CALL XBRZd_l(NK,ARCi,ARKi,ARZi)
!      PRINT 15
   15 format(2x,'gbrzd variant XEK'/)
      SMA=NK
      IF(SMA.LE.16) GO TO 21

!      PRINT 10,NK
   10 FORMAT(2X,'NK>16',3X,'NK=',I3/)
      STOP

   21 CONTINUE
      DO 90 SK=1,SMA
      aRK(sk)=ARKi(SK)
      arz(sk)=ARZi(SK)
      arc(sk)=arci(sk)
   90 CONTINUE
   20 CONTINUE
      R=RI
      Z=ZI
      R2=R*R
      BR=0.0d0
      BZP=0.0d0
      BZR=0.0d0
      Brr=0.0d0

      BZZP=0.0d0
      DO 100 SK=1,SMA
      RK=ARK(SK)
      ZPR=Z-ARZ(SK)
      ZPR2=ZPR*ZPR
      RKPR=RK+R
      RKMR=RK-R
      OZH12=1.0d0/(RKPR**2+ZPR2)
      OZH1=dSQRT(OZH12)
      OZH22=1.0d0/(RKMR**2+ZPR2)

      KAP2=4.0d0*RK*R*OZH12
      XEK=1.0d0-KAP2
      FEEA=(((0.01736506451d0*XEK+0.04757383546d0)*XEK+
     *0.06260601220d0)*XEK+0.44325141463d0)*XEK+1.0d0
      FEEB=(((0.00526449639d0*XEK+0.04069697526d0)*XEK+
     *0.09200180037d0)*XEK+0.24998368310d0)*XEK
      FEKA=(((0.01451196212d0*XEK+0.03742563713d0)*XEK+
     *0.03590092383d0)*XEK+0.09666344259d0)*XEK+1.38629436112d0
      FEKB=(((0.00441787012d0*XEK+0.03328355346d0)*XEK+
     *0.06880248576d0)*XEK+0.12498593597d0)*XEK+0.5d0

      ALNX=-dLOG(XEK)
      EE=ALNX*FEEB+FEEA
      EK=ALNX*FEKB+FEKA
      BC1=(RKMR*RKPR-ZPR2)*OZH22
      KAP2R=(RKMR*RKPR+ZPR2)*4.0d0*RK*OZH12*OZH12
      KAP2Z=-KAP2*2.0d0*ZPR*OZH12
      O2K2=0.5d0/KAP2
      EEK=(EE-EK)*O2K2
      EKK=(EE/XEK-EK)*O2K2
      BC2=BC1*EEK+EKK

      BC3=2.0d0*RK*EE*OZH22*OZH22
      BC4=BC1*EE+EK
      BC5=BC4*OZH12
      BC6=2.0d0*OZH1*ARC(SK)
      BZP=BC4*BC6+BZP
      BZR=((RKMR*RKMR-ZPR2)*BC3+BC2*KAP2R-BC5*RKPR)*BC6+BZR
      BZZP=(BC2*KAP2Z-2.0d0*BC3*ZPR*RKMR-BC5*ZPR)*BC6+BZZP
      IF(R2.EQ.0.0d0) GO TO 40
      OR=1.0d0/R
      BR=((RK*RK+R2+ZPR2)*EE*OZH22-EK)*BC6*ZPR*OR+BR

      BRR=-BR*OR-BZZP
      GO TO 100
   40 BRR=-0.5d0*BZZP
  100 CONTINUE
      BRI=BR
      BZI=BZP
      BRRI=BRR
      BRZI=BZR
      BZRI=BZR
      BZZI=BZZP

      RETURN
      END
C .GBRZY

C .legqkb
      SUBROUTINE LEGQKB_l(N,M,X,Q,Q1,EL,KL)
      implicit double precision(a-h,o-z),integer(i-n)
c      REAL Q2,Q3,Z,A,B,F,F1,BQ,BQ1
      INTEGER C
      E=dABS(EL)
      Z=1d0/(X*X)
      A=0.5d0*(N+M+1.5d0)
      B=0.5d0*(N+M+0.5d0)
      C=N+1
      IF(Z.LT.1d0) GO TO 20
!      PRINT 10

   10 FORMAT(/' APYMEHT LEGQF < 1'//)
!      Q=dSQRT(-1d0)
      RETURN
   20 F=1d0
      F1=1d0
      BQ=1d0
      BQ1=1d0
      K=-1
   30 K=K+1
      Q2=Z/((C+K)*(K+1))

      Q3=K+0.5d0
      BQ=BQ*(A+K)*(B+K)*Q2
      BQ1=BQ1*(A+Q3)*(B+Q3)*Q2
      F=F+BQ
      F1=F1+BQ1
      IF(dABS(BQ).GT.E) GO TO 30
      BQ1=1d0/X
      BQ=-SQRT(X*X-1.)*BQ1
      BQ1=0.5d0*BQ1
      BQ1=BQ**M*0.5d0**N*dSQRT(0.5d0)

      Q2=F*BQ1
      Q=Q2
      Q1=(N+M+0.5d0)*F1*BQ1*BQ-M*Q2/BQ
      KL=K
      RETURN
      END
C .legqkb

c .rkin1
      SUBROUTINE RKIN1d_l(KZ3,L1,XGA,OBRA,RT,R0,BY0,IK,ZK,RK,SMA)

      save /crki1l/,/crki1i/,/crki1e/,/cxbrza/,/cxbrzi/

      INTEGER SMA,SMAC
      LOGICAL OBRAC,OBRA,PROP
      double precision iKC(16),ZKC(16),RKC(16),IK(1),ZK(1),RK(1)
      double precision rt,r0,by0,xga,rtc,by0c,r0c,xgac,pi
      COMMON/CRKI1l/ OBRAC,prop
      COMMON/CRKI1i/ KZ3C,L1C,smac
      COMMON/CRKI1e/ PI,XGAC,RTC,BY0C,R0C
      COMMON/CXBRZa/IKC,RKC,ZKC
      COMMON/CXBRZi/NKC
      PROP=.FALSE.
      SMAC=SMA
      NKC=SMA
      KZ3C=KZ3
      OBRAC=OBRA
      XGAC=XGA
      RTC=RT
      L1C=L1
      BY0C=BY0
      R0C=R0
      PI=4d0*dATAN(1d0)
      IF(SMAC.GT.16) SMAC=16
      DO 10 M=1,SMAC
      IKC(M)=IK(M)
      ZKC(M)=ZK(M)
   10 RKC(M)=RK(M)
      RETURN
      END
C .rkin1
c .xbrz
      SUBROUTINE XBRZd_l(NK,ARC,ARK,ARZ)

      save /cxbrza/,/cxbrzi/

      double precision ARCC(16),ARKC(16),ARZC(16)
      double precision ARC(1),ARK(1),ARZ(1)
      COMMON/CXBRZa/ARCC,ARKC,ARZC
      COMMON/CXBRZi/NKC

      NPK2=NKC
      NK=NKC
      IF(NPK2.LE.16) GO TO 5
!      PRINT 3,NKC
    3 FORMAT(/5X,'NK>16',3X,'NK=',I3/)
      STOP
    5 CONTINUE

      DO 10 K=1,NPK2
      ARC(K)=ARCC(K)
      ARK(K)=ARKC(K)
      ARZ(K)=ARZC(K)
   10 CONTINUE
!      PRINT 60,NK

   60 format(/5x,'  XBRZd  prorabotala   ',
     *2x,'sma=',i3/)
      RETURN
      END
C .xbrz
      SUBROUTINE stevvo_l(RT0,R0i,L1i,cbfi,BY0i,bf0)

      double precision IK(16),ZK(16),RK(16),rt0,r0i,cbfi,by0i,bf0
      INTEGER sma

      INTEGER BMMA,BNC
      LOGICAL NOPRI1,NOPRI2
      double precision RC(6),cbf,rt,r0,by0
      BNC=1
      print *,'this is QHS'
      open(17,form='FORMATTED',file='sdat.dat',status='old')
      READ(17,340)nopri1,nopri2
  340 FORMAT(7X,L5,6X,7X,L5)

      READ(17,120)NMAX,MMA,NMA,RT,R0,L1,NC
  120 FORMAT(5X,I2,2X,4X,I3,2X,5X,I2,2X,2(3X,E10.4,2X),2(3X,I3,2X))
      IF(NC.GT.BNC) NC=BNC
!      PRINT 130,NMAX,BMMA,MMA,NMA,L1,RT,R0,NC
  130 FORMAT(//'NMAX=',I3,4X,4X,'BMMA=',I3,4X//
     *'MMA=',I3,4X,'NMA=',I3,4X,3HL1=,I3,4X,
     *'RT=',E14.7,4X,'R0=',E14.7,4X,'NC=',I3/)
      READ(17,140)cbf,RC
  140 FORMAT(5X,E13.7/3(3X,E13.7)/3(3X,E13.7))
      READ(17,180)BY0
  180 FORMAT(4X,E14.7)
!      PRINT 250,BY0,cbf
  250 FORMAT(/'BY0=',E14.7,3X,'cbf=',E14.7/)
!      PRINT 480,RC
  480 FORMAT(/'RC=',6(E14.7,2X))
      CALL COINhs_l(RT,R0,RC,NC,BMMA,MMA,NMA,NMAX,
     *cbf,L1,NOPRI1,NOPRI2)
      READ(17,420)SMA
  420 FORMAT(5X,I3)
      READ(17,400)IK
  400 FORMAT(4(3X,E13.7))
      READ(17,400)ZK
      READ(17,400)RK
!      PRINT 450
  450 FORMAT('  IK')
!      PRINT 460,IK
  460 FORMAT(8(2X,E14.7))
!      PRINT 470
  470 FORMAT('  ZK')
!      PRINT 460,ZK
!      PRINT 475
  475 FORMAT('  RK')
!      PRINT 460,RK
      CALL RKIN1d_l(40,L1,0d0,.false.,RT,RK(16),BY0,IK,ZK,RK,SMA)
      close(17)
      RT0=rt
!      print *, 'RT0 = ',RT0
      R0i=r0
      L1i=l1
      cbfi=cbf
      BY0i=by0
      bf0=-cbf*L1/rt
      return
      END
