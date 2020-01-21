cfil  The output  : 	ect is DFIRE, ect2 is GOAP_AG
c              GOAP=DFIRE+GOAP_AG
c	all inputs are at /gpfs1/active/hzhou2/fold/dplanedfire/ 
c       which was set in variable base
c	 basic flow:  (1) read in structure list 
c                     (2) read in potential data
c	              (3) loop over each structure
c	       inside the loop: (a) read in structure
c	                        (b) calculate local system vectors for each atom
c	                        (c) loop over all pairs and evaluate
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8(a-h,o-z)
      integer maxa,nwt,newwt,maxp,natom,numdel,
     &        ntyp,maxit,maxres,lengthtemp
      logical iwrit,ioutwt,slmut,slmut1
      real*8 cbox,cx,cy,cz,vdwradi(5),pi,pi3,theta,rcut,sgcut,rd
      parameter(maxa=200000,cbox=18.856,maxit=10000,
     &          ntyp=20,natyp=8,maxres=50000,ibin=20)
      integer idel(maxa),resnum(maxa),restyp(maxa),pol(20),
     &        ihflg(maxa),rnum,ncel,iused(ntyp),imut(maxres),
     &        itgt(maxres),iatyp(maxa),iu(natyp,natyp),
     &        mutmap(ntyp,12),ind1(maxa),ind2(maxa),
     &        map(500),ibk(maxres,15),ianum(20),ib0(maxa)
      real*8 xp(maxa),yp(maxa),zp(maxa),ddg(maxres),resdep(maxres),  
     &        hd(20),asa(20),qq(maxa),rco(maxres),ele(maxa),aa
      real*8 xpn(maxa),ypn(maxa),zpn(maxa),atmdep(maxa),hvnum(ntyp),
     &       icnt(50),pot(20,15,20,15,ibin),rrme(ibin),rghist(1000),
     &       scalco(ibin),ddr(ibin),pb(20,15,20,15,ibin)
        
      real*8 sx,sy,sz,tx,ty,tz,xm(3,4),dep(maxa,maxit),avedep,
     &       entropy(ntyp),volume(ntyp),aenv(maxres),para(ntyp),
     &       para1(ntyp),alf(natyp,natyp),bkbn(natyp,natyp),
     &       chg(ntyp,15),rmsd(3000),cntnum(20,15,20,15,ibin),
     &       ttlcnt(ibin),xmol(20,15),expnum(20,15,20,15,ibin),
     &       unk1(ibin,ibin),unk2(ibin),hnum,ha1,ha2,dh,hnumex,
     &       ymol(20,15),pbtl(ibin),xn(maxa,3),xd(maxa,3),
     &       xn2(3),xd2(3),xn1(3),xd1(3),cnttheta(20,15,20,15,ibin,12,5)
     &      ,vt(3),vt2(3),cnttheta0(12,5),av0(5),sd0(5),ppx(12)
     &      ,xxh(4),yyh(4),zzh(4),dih,dih2,etcrecord(maxa,2)
      character resname*4,atmname*4,ctmp*4,resn(ntyp)*3,rname(maxa)*4,
     &    aname(maxa)*4,line1*80,cct*2,chain*1,ares0*5,ares1*5
     &  ,ctmp1*4,ctmp2*2,ctmp3*4,tempetc1*8,tempetc2*8
       character resnm*3,atnm*3,cind(ntyp,20)*3,base*500

       data qq/maxa*0.0/,ele/maxa*0.0/,aa/0.0/
       data map/500*0/
       data cind/400*'XOP'/
       common /fff/ips,base

       data asa/360.,425.,429.,402.,399.,375.,470.,447.,347.,
     &  323,370.,352.,408.,379.,405.,378.,432.,468.,445.,364./

       data hd/4.19,4.80,6.03,5.87,5.89,4.51,7.14,4.90,2.03,1.01,
     &     2.53,1.59,2.48,1.40,1.78,1.09,2.95,1.78,2.01,3.51/



      data hvnum/2,4,7,4,4,3,10,8,1,0,
     &           3,2,5,4,5,4,6,7,5,3/               

      data pol/1,0,0,0,0,0,0,1,0,0,1,1,1,1,1,1,1,1,1,0/

      data volume/45.,104.,130.,101.,101.,75.,168.,133.,26.,0.,
     &     56.,30.0,86.0,64.,77.,53.,96.0,129.,106.,59./



      data entropy/-0.55,-1.61,-0.58,-0.89,-0.78,-0.51,-0.97,-0.98,
     & 0.0,0.0,-1.63,-1.71,-2.11,-1.57,-1.80,-1.25,-0.96,-2.03,
     & -1.94,0.0/            

      data vdwradi,pi/2.0,1.85,1.70,2.0,0.,3.1415926/

      data resn/'CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR',
     &          'ALA','GLY','THR','SER','GLN','ASN','GLU','ASP',
     &          'HIS','ARG','LYS','PRO'/
   

      common /ddd/icnt


      character*80 whole, filename*100,whole1(maxa),fil1*40,fil2*40,
     &              afil(200000)*80,recordwang(maxa),temppath

          read(*,'(a500)') base

c updating by xiao wang: xiaowang20140001@gmail.com
        i=1 
c        open(unit=20,file='pdbfile.list',status='old')
 2001   read(*,*,end=2000) afil(i)
        i=i+1
        goto 2001

 2000   continue
        nfil=i-1
c        close(20)


            do i=1,500
            map(i)=0
            enddo
            hcut=3.5
            dh=2./ibin

c         base='/home/hzhou2/fold/dplanedfire '
c         base='/gpfs1/active/hzhou2/fold/dplanedfire '
        ips=index(base,' ')-1
c        open(unit=20,file=base(1:ips)//'/fort.21_1.61_10',
c        print *,base
        open(unit=20,file=base(1:ips)//'/fort.21_1.61_2',
     &               status='old')
         read(20,*)
         read(20,*) ibinme,mapnum
         if(ibinme.gt.ibin) then
         write(*,*) 'bin exceeds limit:',ibin,'!'
         stop
         endif
         do i=1,mapnum
         read(20,*) map(i)
         enddo

 3001   read(20,*,end=3000) ctmp1,ctmp2,ctmp3,m,pp,i,j,k,l
c            if(pp.gt.9) pp=0.0
        if(i.eq.21.and.k.eq.21) then
c        hdpot(j,l)=pp
        else
        pot(i,j,k,l,m)=pp
        pot(k,l,i,j,m)=pp
        endif
        goto 3001
 3000   continue
        close(20)

              
          do l1=1,20
          do m1=1,15
          do l2=1,20
          do m2=1,15
          do i=1,ibin
          do l=1,12
          do m=1,5
          cnttheta(l1,m1,l2,m2,i,l,m)=0.0
          enddo
          enddo
          enddo
          enddo
          enddo
          enddo
          enddo

        open(unit=21,file=base(1:ips)//'/fort.31_g72_noshift5_new',
     &               status='old')
         read(21,*)
         read(21,*) ibinme,mapnum,ig 
         if(ibinme.gt.ibin) then
         write(*,*) 'bin exceeds limit:',ibin,'!'
         stop
         endif
         do i=1,mapnum
         read(21,*) map(i)
         enddo
c m:distance bins? i:index for 1st amino acid j: item type index, k:2nd amino acid, l:2nd type index, mx:1-5 angle index lx:angle bins from 1-12 corresponds to 0 degree to 360 degree
 4001   read(21,*,end=4000) ctmp1,ctmp2,ctmp3,m,pp,i,j,k,l
c            if(pp.gt.9) pp=0.0
         do mm=1,5
        read(21,*,end=4000) ctmp1,ctmp2,ctmp3,m,mx,
     &  (ppx(lx),lx=1,12),ppy,i,j,k,l
        if(i.eq.21.and.k.eq.21) then
        else
        do lx=1,12
        cnttheta(i,j,k,l,m,lx,mx)=ppx(lx)
        enddo 
        endif
         enddo 
       
        do lx=1,12
        cnttheta(k,l,i,j,m,lx,1)=cnttheta(i,j,k,l,m,lx,3)
        cnttheta(k,l,i,j,m,lx,2)=cnttheta(i,j,k,l,m,lx,4)
        cnttheta(k,l,i,j,m,lx,3)=cnttheta(i,j,k,l,m,lx,1)
        cnttheta(k,l,i,j,m,lx,4)=cnttheta(i,j,k,l,m,lx,2)
        cnttheta(k,l,i,j,m,lx,5)=cnttheta(i,j,k,l,m,lx,5)
        enddo 


        goto 4001
 4000   continue
        close(21)


          
	open(unit=10,file=base(1:ips)//'/charge_inp.dat',
     &                     status='old')
        read(10,*) 
        id=1
5        read (10,'(a3,1x,i2)',end=100) resnm,mm
        resn(id)=resnm
        ianum(id)=mm
        do i=1,mm
        read(10,*) atnm,xx
        cind(id,i)=atnm
        chg(id,i)=xx        
        enddo
        id=id+1
        goto 5
100     continue
        close(10)
     



         do ii=1,nfil
c            filename='/home/hzhou2/sparks/1010db/'//afil(ii)
            filename=afil(ii)

            if(chain.eq.'0') chain=' '
       ect=0.0
       ect2=0.0
 1    open(11,file=filename,status='old')
	     do ttt=1,maxa
	     do mmm=1,2
			etcrecord(ttt,mmm)=0
		    enddo
	     enddo
          ime=-1 
       id=1
       cx=0
       cy=0
       cz=0
        ires0=0
        ires=0
        ares0='     '
        rnum=0
        resnum(1)=0
 15    line = line + 1
c      print *,'reading finished for  line :',line
      read(11,'(a80)',err=10,end=10) whole1(id) 
c here load all data finish and then goto 10,notice here i should put all the data in recordwang, but i will modify some to install the some atoms' energy to it.
c      if(whole1(id)(22:22).ne.chain.and.ime.eq.1) goto 10 
		
      if(whole1(id)(1:4).ne.'ATOM') goto 15
c      if(whole1(id)(1:4).ne.'ATOM'.or.whole1(id)(22:22).ne.chain)
c     &                     goto 15
         read(whole1(id)(18:21),'(a4)') resname 
             iusd=-1
          do i=1,20
          if(resn(i)(1:3).eq.resname(1:3)) iusd=1
          enddo
          if(iusd.lt.0) goto 15

      ime=1

      if(whole1(id)(17:17).ne.' ') 
     &  READ (whole1(id)(57:60),'(f4.2)') prob
        
      if(whole1(id)(17:17).ne.' '.and.prob.lt.0.5)   goto 15
      if(whole1(id)(17:17).eq.'B'.and.abs(prob-0.5).lt.1.d-4)
     &    goto 15

      READ (whole1(id)(23:27),'(a5)') ares1
        if(ares1.ne.ares0) then
           ires=ires+1
           ares0=ares1
           do i=1,15
           ibk(ires,i)=-1
           enddo
        endif            



      natom=id
      rnum=ires
  

        resnum(rnum+1)=id           
c here resnum denotes the id, specially correspond to atom, so that we can use this to denote our every atom's energy       
c           write(*,'(a80)') whole1(id)
c      print "(i10)",id
c      print *,whole1(id)
      READ (whole1(id)(31:54),'(3F8.3)') Xp(id),Yp(id),Zp(id)
       cx=cx+xp(id)
       cy=cy+yp(id)
       cz=cz+zp(id)




         read(whole1(id)(14:17),'(a4)') atmname 
         if(atmname(1:2).eq.'OT') atmname='O   '   
         if(atmname(1:3).eq.'OXT') atmname='O   '   
         if(resname(1:3).eq.'ILE'.and.atmname(1:3).eq.'CD1') 
     &           atmname='CD  '   
         rname(rnum)=resname
         aname(id)=atmname
          restyp(rnum)=0

         i=id


            ihflg(id)=-1
            ind2(i)=-1
c         if(atmname(1:1).ne.'H'.and.atmname(1:1).ne.'A') then
         if(atmname(1:1).ne.'H')then
            ihflg(id)=1
         iok=-1
         do j1=1,20
         if(resn(j1).eq.resname(1:3)) then
            restyp(rnum)=j1
            do j2=1,15 
            if(cind(j1,j2).eq.atmname(1:3)) then
c            qq(i)=chg(j1,j2)
            ind1(i)=j1
            ind2(i)=j2
            ibk(rnum,j2)=id
            iok=1
            endif
            enddo
         endif
         enddo
          
         if(iok.lt.0) then
c         write(*,*) afil(ii),ires,resname(1:3),atmname(1:3), 
c     &            'not found !'
            ind1(i)=10
            ind2(i)=1
            ibk(rnum,1)=id
c         stop
         endif
         endif 


      natom=id
	  recordwang(id)=whole1(id)
      id=id+1
c      print *,'record wang info',recordwang(id)
      goto 15

c      print *,'reading file finished with residue',rnum

 10   continue
c      print *,'reading file finished with residue',rnum
           close(11)

c setting xn,xd
          do k1=2,rnum
         if((resnum(k1+1)-resnum(k1)).lt.ianum(restyp(k1))) goto 404
           do k=resnum(k1)+1,resnum(k1+1)   
           if(ind2(k).gt.0) then
           call caldplane(xp,yp,zp,k1,k,ind1(k),ind2(k),ibk,xn1,xd1,ibb)
           ib0(k)=ibb
           do l=1,3
           xn(k,l)=xn1(l)
           xd(k,l)=xd1(l)
           enddo 
           endif
           enddo
      
404	continue
          enddo
		 
c         print *,xn
         do k1=2,rnum-1
         if((resnum(k1+1)-resnum(k1)).lt.ianum(restyp(k1))) goto 405
           do k=resnum(k1)+1,resnum(k1+1)   
           cct=aname(k)(1:2)
           if(ind2(k).gt.0) then
           do l=1,3
           xn1(l)=xn(k,l)
           xd1(l)=xd(k,l)
           enddo


             do j1=k1+1,rnum
         if((resnum(j1+1)-resnum(j1)).lt.ianum(restyp(j1))) goto 406
             do j=resnum(j1)+1,resnum(j1+1)     
             cct=aname(j)(1:2)
           if(ind2(j).gt.0.) then
           do l=1,3
           xn2(l)=xn(j,l)
           xd2(l)=xd(j,l)
           enddo

              rd=sqrt((xp(k)-xp(j))**2+(yp(k)-yp(j))**2+
     &                (zp(k)-zp(j))**2)



c             if(jj.gt.ibin) jj=ibin 

             jj=map(int(rd*2))

             if(jj.le.ibin.and.jj.gt.0.1.and.rd.lt.30) then 

             ee=pot(ind1(k),ind2(k),ind1(j),ind2(j),jj)
             ect=ect+ee
			 etcrecord(k,1)=etcrecord(k,1)+ee
			 etcrecord(j,1)=etcrecord(j,1)+ee
c      print *,'starting the final calculation',k1,j1       
      if((j1-k1).ge.ig) then
             xdd=1./rd
             vt(1)=(xp(j)-xp(k))*xdd   ! 1->2
             vt(2)=(yp(j)-yp(k))*xdd   ! 1->2
             vt(3)=(zp(j)-zp(k))*xdd   ! 1->2

             vt2(1)=-vt(1)
             vt2(2)=-vt(2)
             vt2(3)=-vt(3)
c            print *,xn1,xn2 
c           print *,'call calang2' 
             call calang2(xn1,vt,ang1,cs1)
             mm1=int((cs1+1.001)*6.0)+1
             if(mm1.gt.12) mm1=12
c            print *,'call calphi'
             call calphi(xn1,xd1,vt,phi1)
             mm2=int((phi1+180.001)/30.0)+1
             if(mm2.gt.12) mm2=12
c            print *,'call calang2 again'
             call calang2(xn2,vt2,ang2,cs2)
             mm3=int((cs2+1.001)*6.0)+1
             if(mm3.gt.12) mm3=12
c             print *,'cs2',cs2
c           print *,'call calphi again'
             call calphi(xn2,xd2,vt2,phi2)
             mm4=int((phi2+180.001)/30.0)+1
             if(mm4.gt.12) mm4=12

             xxh(1)=xp(j)+xn2(1)
             yyh(1)=yp(j)+xn2(2)
             zzh(1)=zp(j)+xn2(3)
             xxh(2)=xp(j)
             yyh(2)=yp(j)
             zzh(2)=zp(j)
             xxh(3)=xp(k)
             yyh(3)=yp(k)
             zzh(3)=zp(k)
             xxh(4)=xp(k)+xn1(1)
             yyh(4)=yp(k)+xn1(2)
             zzh(4)=zp(k)+xn1(3)
c            print *,'call dihedral'
            call dihedral(xxh,yyh,zzh,dih)

             mm5=min(int((dih+180.001)/30.0)+1,12)


c            print *,'fill all the records data',k,j
c            print *,ind1(k),ind2(k),ind1(j),ind2(j),jj,mm1
c	      print *,cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm1,1)
          temp0=0.0
          
          temp0=temp0+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm1,1)
c          print *,'temp0 calculate finished1',temp0
          temp0=temp0+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm2,2)
c          print *,'temp0 calculate finished2',temp0
c           print *,mm3
          temp0=temp0+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm3,3)
c          print *,'temp0 calculate finished3',temp0
          temp0=temp0+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm4,4)
c		  print *,'temp0 calculate finished4',temp0
          temp0=temp0+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm5,5)
c  		  print *,'temp0 calculate finished',temp0
          flush(6)
          etcrecord(k,2)=temp0+etcrecord(k,2)
		  ect2=ect2+temp0
		  etcrecord(j,2)=temp0+etcrecord(j,2)
c          print *,'etcrecord finished calculation'
           flush(6)
             endif ! .ge.ig

             endif

            endif
c          print *,'calculation call finished',k1,j1
1234              continue
            enddo
406	continue
            enddo

           endif
          enddo


405	continue
         enddo
	  lengthtemp=len_trim(afil(ii))
	  temppath=afil(ii)(9:lengthtemp-4)//'_goap.pdb'
	  write(*,555) ii,' ',temppath,ect+ect2,ect,ect2
	  open(unit=38,file=temppath,status='new')
      do uuu=1,id-1
	    write(tempetc1,'(f8.2)')etcrecord(uuu,1)
		write(tempetc2,'(f8.2)')etcrecord(uuu,2)
		recordwang(uuu)(63:69)=tempetc1
		recordwang(uuu)(71:78)=tempetc2
		write(38,*)recordwang(uuu)
	     enddo
      write(38,*)'RESULT ',ect+ect2
	     close(38)
        write(*,555) ii,' ',afil(ii),ect+ect2,ect,ect2
	   
	   
	   
	   
	   
403	continue
      enddo

555	format(1x,i5,a1,1x,a30,1x,f12.2,1x,2f12.2)

     

        stop
      end

c#################################################
c Uniformlu Distributed Random Number Generator 
c Output:
c     Double Precision RN from 0 to 1:
c Input (changed by Han 2/13/92):
c     Iseed --- random seed for initialization
c     the range should be 0<= Iseed <=31328
c
      FUNCTION RAN(ISEED)
      REAL*8 RAN,RAN1
      INTEGER*4 ISEED
      SAVE INIT
      DATA INIT /1/
      IF (INIT.EQ.1) THEN
        INIT=0
        CALL RMARIN(ISEED,23456)
      END IF

  10  CALL RANMAR(RAN1)
      IF (RAN1.LT.1D-16) GOTO 10
      RAN=RAN1
	RETURN
      END

      FUNCTION RN()
      REAL*8 RN,RAN1
      INTEGER*4 ISEED
      COMMON/SEED/ISEED
      SAVE INIT
      DATA INIT /1/
      IF (INIT.EQ.1) THEN
        INIT=0
        CALL RMARIN(ISEED,23456)
      END IF

  10  CALL RANMAR(RAN1)
      IF (RAN1.LT.1D-16) GOTO 10
      RN=RAN1
	RETURN
      END
C	##############################
      SUBROUTINE RANMAR(RVEC)
*     -----------------
* Universal random number generator proposed by Marsaglia and Zaman
* in report FSU-SCRI-87-50
* In this version RVEC is a double precision variable.
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
      COMMON/ RASET2 / IRANMR,JRANMR
C      SAVE /RASET1/,/RASET2/
      UNI = RANU(IRANMR) - RANU(JRANMR)
      IF(UNI .LT. 0D0) UNI = UNI + 1D0
      RANU(IRANMR) = UNI
      IRANMR = IRANMR - 1
      JRANMR = JRANMR - 1
      IF(IRANMR .EQ. 0) IRANMR = 97
      IF(JRANMR .EQ. 0) JRANMR = 97
      RANC = RANC - RANCD
      IF(RANC .LT. 0D0) RANC = RANC + RANCM
      UNI = UNI - RANC
      IF(UNI .LT. 0D0) UNI = UNI + 1D0
      RVEC = UNI
      END
 
      SUBROUTINE RMARIN(IJ,KL)
*     -----------------
* Initializing routine for RANMAR, must be called before generating
* any pseudorandom numbers with RANMAR. The input values should be in
* the ranges 0<=ij<=31328 ; 0<=kl<=30081
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
      COMMON/ RASET2 / IRANMR,JRANMR
      SAVE /RASET1/,/RASET2/
* This shows correspondence between the simplified input seeds IJ, KL
* and the original Marsaglia-Zaman seeds I,J,K,L.
* To get the standard values in the Marsaglia-Zaman paper (i=12,j=34
* k=56,l=78) put ij=1802, kl=9373
      I = MOD( IJ/177 , 177 ) + 2
      J = MOD( IJ     , 177 ) + 2
      K = MOD( KL/169 , 178 ) + 1
      L = MOD( KL     , 169 )
      DO 300 II = 1 , 97
        S =  0D0
        T = .5D0
        DO 200 JJ = 1 , 24
          M = MOD( MOD(I*J,179)*K , 179 )
          I = J
          J = K
          K = M
          L = MOD( 53*L+1 , 169 )
          IF(MOD(L*M,64) .GE. 32) S = S + T
          T = .5D0*T
  200   CONTINUE
        RANU(II) = S
  300 CONTINUE
      RANC  =   362436D0 / 16777216D0
      RANCD =  7654321D0 / 16777216D0
      RANCM = 16777213D0 / 16777216D0
      IRANMR = 97
      JRANMR = 33
      END
c       ###################
      subroutine calexp(natom1,rg,map)
      implicit real*8 (a-h,o-z)
      parameter(maxa=200000,coe=1.17)
      real*8 x(maxa),y(maxa),z(maxa),icnt(50)      
      integer map(*) 
      common /ddd/icnt
       
        xr=coe*rg
        rho=natom1/(4./3.*3.1416*xr**3)
        dd=1./(rho**(1./3.))
        d2=dd*0.5 

        ir=int(xr/dd)
        mm=mm+1 

         ia=0
   
         do i=-ir,ir
          do j=-ir,ir
           do k=-ir,ir
           rd=sqrt(float(i*i+j*j+k*k))*dd
           if(rd.le.xr) then
           ia=ia+1
           x(ia)=i*dd+d2*rn()
           y(ia)=j*dd+d2*rn()
           z(ia)=k*dd+d2*rn()
           endif
         enddo
         enddo
         enddo
         natom=ia
c       write(*,*) mm,xr,rho,natom,natom1
       do i=1,natom-1
        do j=i+1,natom                            
         rd=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+
     &           (z(i)-z(j))**2)
c        jj=int(rd-1.5)+1
        jj=map(int(rd*2))
        if(jj.gt.0.1.and.jj.lt.50) icnt(jj)=icnt(jj)+1.0
        enddo
       enddo


      return
       end

c##################################################
         subroutine findh(k1,j1,ibk,xp,yp,zp,ha1,ha2)  
         integer k1,j1,ibk(2000,15)
         real*8 xp(*),yp(*),zp(*),rd,ha1,ha2
         real*8 hnum,x,y,z,xh,yh,zh,dx,dy,dz,dr,dx1,dy1,dz1,dr1
c     k1 donor j1 acceptor
         data rd0/3.0/ 
c     determine h(k1):
         if(k1.le.1) return ! cann't determine h(k1)
         k0=k1-1
         l0=ibk(k0,1)  ! C
         l1=ibk(k1,2)  ! CA
         l3=ibk(k1,3)  ! N
         m1=ibk(j1,1)  ! C
         m2=ibk(j1,4)  ! O


         x=(xp(l0)+xp(l1))*0.5        
         y=(yp(l0)+yp(l1))*0.5        
         z=(zp(l0)+zp(l1))*0.5                  
         dx=xp(l3)-x
         dy=yp(l3)-y
         dz=zp(l3)-z
         dr=1./sqrt(dx**2+dy**2+dz**2)
         xh=xp(l3)+dx*dr        
         yh=yp(l3)+dy*dr        
         zh=zp(l3)+dz*dr        
         dx1=xp(m2)-xh         
         dy1=yp(m2)-yh         
         dz1=zp(m2)-zh         
         dr1=1./sqrt(dx1**2+dy1**2+dz1**2)         
         ha1=-(dx*dx1+dy*dy1+dz*dz1)*dr*dr1

         dx=xp(m1)-xp(m2)  
         dy=yp(m1)-yp(m2)  
         dz=zp(m1)-zp(m2)  
         dr=1./sqrt(dx**2+dy**2+dz**2)
         dx1=xh-xp(m2)  
         dy1=yh-yp(m2)  
         dz1=zh-zp(m2)  
         dr1=1./sqrt(dx1**2+dy1**2+dz1**2)
         ha2=(dx*dx1+dy*dy1+dz*dz1)*dr*dr1        
                 
c                write(*,*) ha1,ha2,k1,j1,l0,l1,l3,m1,m2

         return
         end
c	##########################
        subroutine caldplane(xp,yp,zp,k1,k,id1,id2,ibk,xn,xd,ibb)
        real*8 xp(*),yp(*),zp(*),xn(*),xd(*),xn1(3),xd1(3)
        real*8 v1(3),v2(3),v3(3),sidegeo(20,4,15)
         character ctmp*3,base*70
        parameter(maxres=50000)
        integer ibk(maxres,15),sidelb(20,3,15)
         data init/0/,sidelb/900*0/
       common /fff/ips,base

         if(init.eq.0) then
             init=1

       open(unit=10,file=base(1:ips)//'/side_geometry.dat',
     & status='old')
        do i=1,20
       read(10,'(1x,a3,1x,i2)') ctmp,l
       read(10,'(1x,i2,1x,10(f9.3,1x,i2))') k0,
     &   (sidegeo(i,1,j),sidelb(i,1,j),j=1,l)
       read(10,'(1x,i2,1x,10(f9.3,1x,i2))') k0,
     &   (sidegeo(i,2,j),sidelb(i,2,j),j=1,l)
       read(10,'(1x,i2,1x,10(f9.3,1x,i2))') k0,
     &   (sidegeo(i,3,j),sidelb(i,3,j),j=1,l)
       read(10,'(1x,i2,1x,10(f9.3,1x,i2))') k0,
     &   (sidegeo(i,4,j),sidelb(i,3,j),j=1,l)
        enddo

         endif 

         i1=k 
         n0=-1
         if(id2.eq.1) then ! N
           i2=ibk(k1-1,3)        ! k1-1's C
           i3=ibk(k1,2)          ! k1's CA
           n0=1
         elseif(id2.eq.2) then !  CA
           i2=ibk(k1,1)         !k1  N
           i3=ibk(k1,3)         !k1  C
           n0=1
         elseif(id2.eq.3) then !  C
           i2=ibk(k1,2)         !k1  CA
           i3=ibk(k1,4)         !k1  O
           n0=1
         elseif(id2.eq.4) then !  O
           i2=ibk(k1,3)         !k1  C
           i3=ibk(k1,2)         !k1  CA
        elseif(id2.gt.4) then
           i2=ibk(k1,sidelb(id1,1,id2-4))         !k1  CB = sidelb(id1,1,2)
           i3=-1
           if(sidelb(id1,1,id2-3).gt.0) then
              do n=id2+1,15
              if(sidelb(id1,1,n-4).eq.id2) then 
              i3=ibk(k1,n)
              n0=1
              endif 
              enddo
           endif
           if(i3.lt.0) i3=ibk(k1,sidelb(id1,2,id2-4))
        endif
           ibb=i2

c           v1(1)=xp(i2)-xp(i1)
c           v1(2)=yp(i2)-yp(i1)
c           v1(3)=zp(i2)-zp(i1)
c           v2(1)=xp(i3)-xp(i1)
c           v2(2)=yp(i3)-yp(i1)
c           v2(3)=zp(i3)-zp(i1)
c           call xdot(v1,v2,xn)
c           xd(1)=v1(1)+v2(1)
c           xd(2)=v1(2)+v2(2)
c           xd(3)=v1(3)+v2(3)


           if(n0.lt.0) then
           xn(1)=xp(i2)-xp(i1)
           xn(2)=yp(i2)-yp(i1)
           xn(3)=zp(i2)-zp(i1)
           else
           xn(1)=xp(i2)-xp(i1)+xp(i3)-xp(i1)
           xn(2)=yp(i2)-yp(i1)+yp(i3)-yp(i1)
           xn(3)=zp(i2)-zp(i1)+zp(i3)-zp(i1)
           endif

           v2(1)=xp(i3)-xp(i1)
           v2(2)=yp(i3)-yp(i1)
           v2(3)=zp(i3)-zp(i1)

           call xdot(xn,v2,v1)
           call xdot(v1,xn,xd)


c	normalize
        xx1=1./sqrt(xd(1)**2+xd(2)**2+xd(3)**2)
        xx2=1./sqrt(xn(1)**2+xn(2)**2+xn(3)**2)

        do l=1,3
        xd(l)=xd(l)*xx1
        xn(l)=xn(l)*xx2
        enddo



        return
        end
       




c	#############################
        subroutine xdot(x,y,z)
         real*8 x(*),y(*),z(*)

         z(1)=x(2)*y(3)-x(3)*y(2)
         z(2)=x(3)*y(1)-x(1)*y(3)
         z(3)=x(1)*y(2)-x(2)*y(1)

         return
         end
c       ############################
       subroutine calang(xn,xh1,xo,ang,cs)
       real*8 xn(3),xh1(3),xo(3),ang,cs,r1,r2,xt1,xt2,co
       data co/57.296/

         cs=0.0
         r1=0.0
         r2=0.0
       do i=1,3
          xt1=xn(i)-xh1(i)
          xt2=xo(i)-xh1(i)
          cs=cs+xt1*xt2
          r1=r1+xt1**2
          r2=r2+xt2**2
       enddo

       cs=cs/sqrt(r1*r2)

       ang=acos(cs)*co

       return
       end

c     #######################################
       subroutine calang2(xn,xh1,ang,cs)
       real*8 xn(3),xh1(3),xo(3),ang,cs,r1,r2,xt1,xt2,co
       data co/57.296/

         cs=0.0
         r1=0.0
         r2=0.0
c       print *,xn,xh1
       do i=1,3
          xt1=xn(i)
          xt2=xh1(i)
          cs=cs+xt1*xt2
          r1=r1+xt1**2
          r2=r2+xt2**2
       enddo
c       print *,'r1*r2 result',sqrt(r1*r2)
       cs=cs/sqrt(r1*r2)

       ang=acos(cs)*co

       return
       end


c          #######################################
       subroutine calphi(xn,xd,xh1,ang)
       real*8 xn(3),xh1(3),xo(3),ang,cs,r1,r2,xt1,xt2,co
     &  ,xd(3),xy(3)
       data co/57.296/

         cs=0.0
         r1=0.0
         r2=0.0
       do i=1,3
          xt1=xn(i)
          xt2=xh1(i)
          cs=cs+xt1*xt2
          r1=r1+xt1**2
          r2=r2+xt2**2
       enddo
       cs=cs/sqrt(r1*r2)
       ss=sqrt(1-cs*cs+0.0001)
        xo(1)=xh1(1)-cs*xn(1)
        xo(2)=xh1(2)-cs*xn(2)
        xo(3)=xh1(3)-cs*xn(3)

       cs2=xo(1)*xd(1)+xo(2)*xd(2)+xo(3)*xd(3)
       cs2=cs2/sqrt(xo(1)**2+xo(2)**2+xo(3)**2)
       call xdot(xn,xd,xy)
       cs3=xy(1)*xo(1)+xy(2)*xo(2)+xy(3)*xo(3)
       if(cs2.lt.-1.0) cs2=-0.9999
       if(cs2.gt.1.0) cs2=0.9999
       ang=acos(cs2)*co
       if(cs3.lt.0) ang=-ang
c       write(*,*) ang,' YYY ',cs2,cs3

       return
       end

c	################
      subroutine dihedral (x,y,z,dih)
      parameter (Nn=4)
      real*8  x(Nn),y(Nn), z(Nn),dih
      real*8 mconst
c
c bond angle distribution
c
      pi = 4.d0*datan(1.d0)
         i  = 1
         ipi = i
         jpi = i+1
         kpi = i+2
         lpi = i+3
c
c dihedral angle
         dihedral_angle = 0.d0
c
         rijx = x(ipi) - x(jpi)
         rijy = y(ipi) - y(jpi)
         rijz = z(ipi) - z(jpi)
         rjkx = x(jpi) - x(kpi)
         rjky = y(jpi) - y(kpi)
         rjkz = z(jpi) - z(kpi)
         rklx = x(kpi) - x(lpi)
         rkly = y(kpi) - y(lpi)
         rklz = z(kpi) - z(lpi)
c
         ax = rijy*rjkz - rijz*rjky
         ay = rijz*rjkx - rijx*rjkz
         az = rijx*rjky - rijy*rjkx
         bx = rjky*rklz - rjkz*rkly
         by = rjkz*rklx - rjkx*rklz
         bz = rjkx*rkly - rjky*rklx
C
C     set to MCONST if smaller than MCONST
         mconst=1.0d-10
         RG=SQRT(MAX(MCONST,RJKX*RJKX+RJKY*RJKY+RJKZ*RJKZ))
         RGR=1.d0/RG
         RA2R=1.d0/MAX(MCONST,AX*AX+AY*AY+AZ*AZ)
         RB2R=1.d0/MAX(MCONST,BX*BX+BY*BY+BZ*BZ)
         RABR=SQRT(RA2R*RB2R)
CP=COS(PHI)
         CP=RABR*(AX*BX+AY*BY+AZ*BZ)
C SP=SIN(PHI)
C which can be simplify to sin(phi)=|G|H.A/(|A|.|B|)
Cab...B950603 06/29/95
Cab        SP=RG*RABR*(AX*HX+AY*HY+AZ*HZ)
         SP=-RG*RABR*(AX*RKLX+AY*RKLY+AZ*RKLZ)
Cab...
C
         if(cp.gt.1.d0) then
c            write(99,*) 'dh:', cp, ipi,jpi,kpi,lpi
            cp = 1.d0
            theta  = 0.d0
         else if(cp.lt.-1.d0) then
c            write(99,*) 'dh:', cp, ipi,jpi,kpi,lpi
            cp = -1.d0
            theta  = 180.d0
         else             
            theta    = acos(cp)
            theta    = theta*180/pi
         endif
         if(sp.lt.0.0) theta = -theta          
         dih = theta
         return
      end












