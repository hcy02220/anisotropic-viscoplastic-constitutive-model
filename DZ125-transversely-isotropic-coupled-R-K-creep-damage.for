C#########################################################################
C#########################################################################
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
      implicit none
C
C=====================================================================
c=====================================================================
C
      character*16 cmname
      
      integer ndi, nshr, ntens, nstatv, nprops, noel, npt, layer, kspt,
     1 kinc, kstep, k1, k2, temp_num, param_num,i,j,o_max_plane
      integer o_yield_cri(12), c_yield_cri(6)
      logical time_match,crossed
      
      double precision stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3),
     4 param(props(2)), euler(3), trans_a(6,6), trans_b(6,6),
     5 trans_c(6,6), trans_d(6,6), trans_m(3,3), p_stran_r(6),
     6 delta_stress(6), delta_p_stran(6), o_stress(12), o_gama(12),
     7 stress_loc(6), o_schmid(6,12), d_o_gama(12),o_slip_plane(3,12),
     8 delta_stress_d(6),dev_delta_stress_d(6)
      double precision o_x_stress_r(12), o_x_stress(12), c_x_stress(6),
     1 c_x_stress_r(6), o_iso_hard(12), o_iso_hard_r(12), c_iso_hard(6),
     2 c_iso_hard_r(6), o_effect_stress(12), c_effect_stress(6),
     3 o_gama_r(12), c_gama_r(6),max_stress_n_temp(3),stress_temp(3,3)
	
      
      double precision sse, spd, scd, rpl, temp, dtemp, pnewdt, celent,
     1 drpldt,dtime, emod, enu, ebulk3, eg2, eg, elam, err_cri,tempv,
     2 subinc, error, fai, theta, sigma,max_stress_n,max_gama_range,
     3 PERIOD,ct,steptime,res,newcyc,cyc_num,target_time,start_timepoint
	integer  yie(12)
	double precision E11,E33,G23,MIU12,MIU23,K,N,K0,A1,A2,C1,C2,M1,M2,
	1 PHAI_S,OMEGA,B,Q,M11,M33,M55,EFF_S,J2_BACK,R,EFF_P_R,j2_s_d,f,D_r,
	1 EFF_P,PHAI_P,R_R,DELTA,D,r1,r2,beta1,beta2,pa,pr,pk,j2_back1,j2_back2
	double precision S_BACK(6),DEV_STRESS(6),DEV_BACK(6),MM(6,6),
	1 MM_T(6,6),EFF_S_T(6),S_BACK1(6),S_BACK2(6),S_BACK_R1(6),
     2 S_BACK_R2(6),DDSDDET(6,6),p_stran(6),EFF_P_R_T(6),S(6),
     3 stress_d(6),dev_stress_d(6)
C--------------------------------------------------------------------
C      在线调试语句部分
C--------------------------------------------------------------------
C      logical, save :: FirstCall = .true.
C      integer :: dummyVar
C      if (FirstCall==.true.) then
C          FirsrCall = .false.
C          read(*,*) dummyVar
C      end if
C      dummyVar = 1234
C---------------------------------------------------------------------   
      

C读取材料参数      
C     误差下限(建议值1e-5，杨晓光书)
      err_cri = 2e-5
C     err_U = 1e-4,误差上限
      period = 0.0125
      start_timepoint = 0.009375
      E11 = props(1)
      E33 = props(2)
      G23 = props(3)
      MIU12 = props(4)
      MIU23 = props(5)
      omega = PROPS(6) 
      M33 = PROPS(7) 
      M55 = PROPS(8) 
      PA = PROPS(9) 
      PR = PROPS(10)
      PK = PROPS(11)
      k = props(12)
	n = PROPS(13)
      k0 = PROPS(14)
      b = PROPS(15)
      q = PROPS(16)
      m11 = PROPS(17)
      a1 = PROPS(18)
      a2 = PROPS(19) 
      c1 = PROPS(20) 
      c2 = PROPS(21) 
      m1 = props(22)
      m2 = props(23)
      phai_s = props(24)
      r1 = props(25)
      r2 = props(26)
      beta1 = props(27)
      beta2 = props(28)
      
C读取状态变量,背应力，等效塑性应变，背应力分量，各向同性硬化
	do i = 1,6
		s_back(i) = statev(i)
          s_back1(i) = statev(6+i)
          s_back2(i) = statev(12+i)
          p_stran(i) = statev(18+i)
          p_stran_r(i) = statev(24+i)
      enddo
      
	eff_p = statev(31)
	R = statev(32)
      D = statev(33)
      D_r = statev(79)
      subinc = statev(34)

      do i = 1,12
          o_gama(i) = statev(34+i)
      enddo 
C组装弹性刚度矩阵
      DDSDDE = 0.0              
      delta = (1.0+miu12)*(1.0-miu12-2.0*miu23**2*E33/E11)/E11/E11/E33
	DDSDDE(1,1) = (1.0-MIU23**2.0*E33/E11)/E11/E33/DELTA
      DDSDDE(1,2) = (MIU12+MIU23**2*E33/E11)/E11/E33/DELTA
      DDSDDE(1,3) = (MIU23+MIU12*MIU23)*E33/E11/E11/E33/DELTA
      DDSDDE(2,1) = DDSDDE(1,2)
      DDSDDE(2,2) = DDSDDE(1,1)
      DDSDDE(2,3) = DDSDDE(1,3)
      DDSDDE(3,1) = DDSDDE(1,3)
      DDSDDE(3,2) = DDSDDE(2,3)
      DDSDDE(3,3) = (1.0-MIU12**2)/E11/E11/DELTA
      DDSDDE(4,4) = E11/2.0/(1+MIU12)
      DDSDDE(6,6) = G23
      DDSDDE(5,5) = G23
C	反映各向异性屈服面的张量
	MM = 0.0
	MM(1,1) = M11
	MM(2,2) = M11
	MM(3,3) = M33
	MM(5,5) = M55
	MM(6,6) = M55
	MM(1,2) = 0.5*M33-M11
      MM(2,1) = MM(1,2)
	MM(1,3) = -0.5*M33
      MM(3,1) = MM(1,3)
	MM(2,3) = MM(1,3)
	MM(3,2) = MM(2,3)
	MM(4,4) = 2.0*M11-0.5*M33

C	更新应力
	call stress_strain(ddsdde, stress, dstran, p_stran_r,ntens, 
     1                            dtime, delta_stress,delta_p_stran)
      p_stran =  p_stran+delta_p_stran
C     损伤应力
      do i = 1,6
          stress_d(i) = stress(i)/(1.0-D)
      enddo
      
C	应力和背应力偏量
	call deviation_form(stress_d,dev_stress)
	call deviation_form(s_back,dev_back)

C	有效应力
      s = dev_stress - dev_back
C     屈服函数
      eff_s = 0.0	
      do i = 1,6
          eff_s_t(i) = 0.0
          DO J = 1,6
              eff_s_t(i) = eff_s_t(i)+s(j)*MM(j,i)
          enddo
          if (i<=3)then
              eff_s = eff_s + eff_s_t(i)*s(i)
          else
              eff_s = eff_s + 2.0*eff_s_t(i)*s(i)
          endif
      enddo
      eff_s = sqrt(1.5*eff_s)
      f = eff_s - R - k0
      
C	晶体坐标与整体坐标重合。计算各滑移系的切应变率
      call shear_stran_transfer(dstran,d_o_gama,o_slip_plane)   
      o_gama = o_gama + d_o_gama


C	各滑移系最小塑性应变和最大塑性应变，状态变量初始值为0！    
      do k1 = 1,12
	    if ( o_gama(k1) - statev(46+k1) < 0.0 )then
		    statev(46+k1) = o_gama(k1)
	    endif
		
	    if ( o_gama(k1) - statev(58+k1) > 0.0 ) then
		    statev(58+k1) = o_gama(k1)
	    endif
	enddo
      max_gama_range = 0.0
	o_max_plane = 1
      do i =1,12
		if (statev(58+i)-statev(46+i)-max_gama_range > 0.0 )then
			max_gama_range = statev(58+i)-statev(46+i)
			o_max_plane = i
		endif
      enddo
C     max shear strain amp of current cycle
      if (statev(71) < max_gama_range)then
          statev(71) = max_gama_range
      endif
          
      do k1 = 1,3
		stress_temp(k1,k1) = stress(k1)
	enddo
	stress_temp(1,2) = stress(4)
	stress_temp(2,1) = stress(4)
	stress_temp(1,3) = stress(5)
	stress_temp(3,1) = stress(5)
	stress_temp(2,3) = stress(6)
	stress_temp(3,2) = stress(6)
      
C	等效塑性应变增量
	do k1 = 1,6
		if (k1 < 4)then
      		statev(72) = statev(72) + delta_p_stran(k1)**2
          else
              statev(72) = statev(72) + 2.0*delta_p_stran(k1)**2
          endif
	enddo
	statev(73) = (2.0/3.0*statev(72))**0.5
      
      do k1 = 1,3
		stress_temp(k1,k1) = stress(k1)
	enddo
	stress_temp(1,2) = stress(4)
	stress_temp(2,1) = stress(4)
	stress_temp(1,3) = stress(5)
	stress_temp(3,1) = stress(5)
	stress_temp(2,3) = stress(6)
	stress_temp(3,2) = stress(6)
      do i = 1,3
		max_stress_n_temp(i) = 0.0
		do j = 1,3
			max_stress_n_temp(i) = max_stress_n_temp(i) + 
	1 o_slip_plane(j,o_max_plane)*stress_temp(j,i)
		enddo
	enddo
	
	max_stress_n = 0.0
	do i = 1,3
		max_stress_n = max_stress_n + 
	1 max_stress_n_temp(i)*o_slip_plane(i,o_max_plane)
      enddo
      
      if (max_stress_n - statev(74) > 1E-20)then
          statev(74) = max_stress_n
      endif
      
C     fatemi damage
      statev(75) =  statev(71)/2.0*
	1	(1.0+statev(74)/905.0)
C	accumulated plastic strain energy density
	do i = 1,6
		statev(76) = statev(76) + stress(i)*p_stran_r(i)*dtime
	enddo
	
C	accumulated shear strain
	do i = 1,12
		statev(77) = statev(77) + abs(d_o_gama(i))
      enddo
      
C	背应力的第二主不变量
	j2_back1 = 0.0
      j2_back2 = 0.0
      
	call inv2(s_back1,j2_back1)
      call inv2(s_back2,j2_back2)

C	根据当前的误差确定下一步步长
      call error_control(delta_stress, delta_p_stran, err_cri, 
     1 subinc, pnewdt,error,E33)

C	在增量步结尾判断是否屈服，如果屈服，计算出相应的内变量变化率
	if (f > 1e-16)then
		p_stran_r = 1.5*(f/k)**n*matmul(MM,s)/eff_s/(1.0-D)
		eff_p_r = (f/K)**n
		eff_p = eff_p + eff_p_r*dtime
		phai_p = phai_s+(1-phai_s)*exp(-omega*eff_p)
		s_back_r1 = 2.0/3.0*c1*a1*p_stran_r-c1*phai_p*
	1	abs(j2_back1/a1*phai_p)**m1*s_back1*eff_p_r
     2     - beta1*j2_back1**(r1-1)*s_back1
		s_back_r2 = 2.0/3.0*c2*a2*p_stran_r-c2*phai_p*
	1	abs(j2_back2/a2*phai_p)**m2*s_back2*eff_p_r
     2     - beta2*j2_back2**(r2-1)*s_back2
		s_back1 = s_back1+s_back_r1*dtime
		s_back2 = s_back2+s_back_r2*dtime
		s_back = s_back1+s_back2
		r_r = b*(q-r)*eff_p_r
		r = r+r_r*dtime
      endif
      j2_s_d = 0.0
      call inv2(dev_stress,j2_s_d)
C     注意这里算出来的蠕变断裂时间单位是h，所以除3600
      D_r = (j2_s_d/pa)**pr*(1.0-D)**(-1.0*pk)/3600.0
      D = D + D_r*dtime
      statev(79) = D_r
      
C     判断是否进入新的循环，当前时间等于分析步开始时间＋时间增量
      cyc_num = statev(80)
      ct = time(2)+dtime
      steptime = time(1)+dtime
c      target_time = period * cyc_num + start_timepoint
c      time_match =  ABS(steptime - target_time) < 0.9 E-10
C      crossed = (time(1)<target_time) .and. (steptime>=target_time) 
      
	if (steptime - cyc_num * period - start_timepoint > 1e-20)then
c          if (statev(80)<cyc_num+0.5)then
            cyc_num = cyc_num + 1.0
            do k1 = 1,12
                statev(46+k1) = 0.0
                statev(58+k1) = 0.0
            enddo
            statev(81) = statev(81) + statev(75)
            statev(71) = 0.0
            statev(74) = 0.0
            statev(75) = 0.0
c          endif
      endif
      do i = 1,6
          statev(i) = s_back(i)
          statev(6+i) = s_back1(i)
          statev(12+i) = s_back2(i)
          statev(18+i) = p_stran(i)  
	    statev(24+i) = p_stran_r(i)
      enddo
      
      statev(31) = eff_p
	statev(32) = r 
      statev(33) = D
      statev(34) = subinc
      do i = 1,12
           statev(34+i) = o_gama(i)
      enddo 
      statev(80) = cyc_num
      
      return
      end
      
C#######################################################################
      subroutine stress_strain(ddsdde, stress, dstran, p_stran_r, 
     1 ntens,dtime, delta_stress,delta_p_stran)
C     功能：在每个增量步内按照应力-应变关系，更新应力
C     流程：应力增量等于弹性矩阵*（总应变增量-弹性应变增量）
C     输入：弹性矩阵，应变增量
C     输出：增量步内的应力增量
      implicit none
      integer ntens
      double precision ddsdde(ntens,ntens),stress(ntens), dstran(ntens),
     1 p_stran_r(ntens), delta_stress(ntens),delta_p_stran(ntens)
      
      double precision dtime
      
      integer k1,k2
C---------------------------------------------------------------------
C     STRESS，基于全局坐标的计算
      
      do k1 = 1, ntens
          delta_p_stran(k1) = 0.0
          delta_p_stran(k1) = p_stran_r(k1) * dtime
      end do
      
      do k1=1, ntens
          delta_stress(k1) = 0.0
          do k2=1, ntens
              delta_stress(k1) = delta_stress(k1) + ddsdde(k1, k2)
     1                               *(dstran(k2) - delta_p_stran(k2))
          end do
          stress(k1) = stress(k1) + delta_stress(k1)
      end do
C---------------------------------------------------------------------
      return
      end
C#######################################################################
C子程序stress-strain结束
C#######################################################################
      subroutine inv2(s,j2)
      implicit none
      double precision s(6)
      double precision j2
      integer i
      do i = 1,3
          j2 = j2 + s(i)**2
      enddo    
      do i = 4,6
          j2 = j2 + 2.0*s(i)**2
      enddo
      j2 = sqrt(1.5*j2)
      return 
      end
C######################################################################
      subroutine deviation_form(state, state_dev)
C     功能：求应力或背应力的偏量
C     调用：yield_cri子程序
C     输入：应力张量或背应力张量
C     输出：应力偏张量或背应力偏张量
      implicit none
      double precision state(6), state_dev(6)
      
      double precision hydro_pres
      integer k1
      
      hydro_pres = 1.0/3.0*(state(1) + state(2) + state(3))
      
      do k1 = 1, 6
          if (k1 <=3) then
              state_dev(k1) = state(k1) - hydro_pres
          else
              state_dev(k1) = state(k1)
          end if
      end do
      
      return
      end
      
C#######################################################################
C子程序deviation_form结束
C#######################################################################
      
      
c#######################################################################
      subroutine error_control(delta_stress, delta_p_stran,  
     1 err_cri,subinc, pnewdt, error,E33)
C     功能：通过判断误差控制步长
C     调用：主程序
C     输入：应力增量、塑形应变增量、步长控制状态变量、设定误差
C     输出：步长控制变量 pnewdt
      implicit none
      double precision delta_stress(6), delta_p_stran(6)
      double precision err_cri, subinc, pnewdt, error, 
     1 E33, delta_p_stran_c(6)
      
      double precision delta_stress_dev(6), error_stress
      integer k1
      
      call deviation_form(delta_stress, delta_stress_dev)
      error = 0.0
      do k1 =1,6
          error = error + delta_p_stran(k1)**2
      end do
      error = sqrt(2.0/3.0*error)
      error_stress = 0.0
      do k1 =1,6
          error_stress = error_stress + delta_stress_dev(k1)**2
      end do
      error_stress = sqrt(1.5*error_stress)/2.0/E33
      error = error_stress+error
      if(subinc <=4 ) then
          if(error <= err_cri ) then
              pnewdt = 2.0                                      
              subinc = 0.0
          elseif(error >= err_cri .and. error <= 10*err_cri) then
              pnewdt = 1.0
              subinc = 0.0
          else
              pnewdt = 0.5
              subinc = subinc + 1.0
          end if
              
      else
	     pnewdt = 1.0
           subinc = subinc + 1.0
      end if
      
      return
      end
C#######################################################################
C子程序error_control结束
C#######################################################################


      subroutine shear_stran_transfer(dstran,d_o_gama,o_slip_plane)
C     功能：求不同滑移系上的分解剪切应力
C     调用：主程序
C     输入：应力、转换矩阵
C     输出：八面体滑移系与立方对称滑移系下的分解剪切应力
      implicit none
      
      double precision dstran(6),d_o_gama(12)
      integer k1,k2
      double precision o_slip_direction(3,12), o_slip_plane(3,12),
     1  o_schmid(6,12),o_gama(12)

C#######################################################################
C指定八面体主滑移系滑移面法向
      do k1 = 1,12
          
          if (k1 <= 3) then
              o_slip_plane(1, k1) = sqrt(3.0)/3
              o_slip_plane(2, k1) = sqrt(3.0)/3
              o_slip_plane(3, k1) = sqrt(3.0)/3
              
          elseif (k1 <= 6) then
              o_slip_plane(1, k1) = -sqrt(3.0)/3
              o_slip_plane(2, k1) = sqrt(3.0)/3
              o_slip_plane(3, k1) = sqrt(3.0)/3
              
          elseif (k1 <= 9) then
              o_slip_plane(1, k1) = sqrt(3.0)/3
              o_slip_plane(2, k1) = -sqrt(3.0)/3
              o_slip_plane(3, k1) = sqrt(3.0)/3
              
          elseif (k1 <= 12) then
              o_slip_plane(1, k1) = sqrt(3.0)/3
              o_slip_plane(2, k1) = sqrt(3.0)/3
              o_slip_plane(3, k1) = -sqrt(3.0)/3
              
          end if
          
      end do
C八面体主滑移系滑移面法向指定完毕
      
C#######################################################################
C指定八面体主滑移系滑移方向      
      do k1 = 1, 12
          
          if(k1 == 1 .or. k1 == 6) then
              o_slip_direction(1, k1) = 0.0
              o_slip_direction(2, k1) = -sqrt(2.0)/2
              o_slip_direction(3, k1) = sqrt(2.0)/2
              
          elseif(k1 == 2 .or. k1 == 9) then
              o_slip_direction(1, k1) = sqrt(2.0)/2
              o_slip_direction(2, k1) = 0.0
              o_slip_direction(3, k1) = -sqrt(2.0)/2
          
          elseif(k1 == 3 .or. k1 == 12) then
              o_slip_direction(1, k1) = -sqrt(2.0)/2
              o_slip_direction(2, k1) = sqrt(2.0)/2 
              o_slip_direction(3, k1) = 0.0
          
          elseif(k1 == 4 .or. k1 == 11) then
              o_slip_direction(1, k1) = sqrt(2.0)/2
              o_slip_direction(2, k1) = 0.0
              o_slip_direction(3, k1) = sqrt(2.0)/2
              
          elseif(k1 == 5 .or. k1 == 8) then
              o_slip_direction(1, k1) = sqrt(2.0)/2
              o_slip_direction(2, k1) = sqrt(2.0)/2 
              o_slip_direction(3, k1) = 0.0
          
          elseif(k1 == 7 .or. k1 == 10) then
              o_slip_direction(1, k1) = 0.0
              o_slip_direction(2, k1) = sqrt(2.0)/2
              o_slip_direction(3, k1) = sqrt(2.0)/2 
              
          end if
      end do
C八面体主滑移系滑移方向指定完毕
      
C############################################################
C计算八面体滑移系Schmid因子
      do k1 = 1, 12
C公式勘误
C前三个分量与正应力、正应变对应
          o_schmid(1, k1) = o_slip_direction(1, k1) * o_slip_plane(1,k1)
          o_schmid(2, k1) = o_slip_direction(2, k1) * o_slip_plane(2,k1)
          o_schmid(3, k1) = o_slip_direction(3, k1) * o_slip_plane(3,k1)
C后三个分量与切应力、切应变对应
          o_schmid(4, k1)=0.5*(o_slip_direction(1,k1)*o_slip_plane(2,k1)
     1               + o_slip_direction(2, k1) * o_slip_plane(1,k1))
          o_schmid(5, k1)=0.5*(o_slip_direction(1,k1)*o_slip_plane(3,k1)
     1               + o_slip_direction(3, k1) * o_slip_plane(1,k1))
          o_schmid(6, k1)=0.5*(o_slip_direction(2,k1)*o_slip_plane(3,k1)
     1               + o_slip_direction(3, k1) * o_slip_plane(2,k1))
      end do
      
      do k1 = 1,12
          d_o_gama(k1) = 0.0
          do k2 = 1,6
                d_o_gama(k1) = d_o_gama(k1)+o_schmid(k2,k1)*dstran(k2)
          enddo
      enddo
      return
      end

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA) 

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(82)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2) 
      DIMENSION ARRAY(82),JARRAY(82),JMAC(*),JMATYP(*),COORD(*)

C      Error counter:
      JERROR = 0
C      Stress tensor:
      CALL GETVRM('SDV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)
C         #FATEMI,累积塑性应变能，累积切应变，累积塑性应变，蠕变损伤 
      JERROR = JERROR + JRCD
      UVAR(1) = ARRAY(75)
      UVAR(2) = ARRAY(81)
      UVAR(3) = ARRAY(80)
      UVAR(4) = ARRAY(31)
      UVAR(5) = ARRAY(33)
      RETURN
      END
