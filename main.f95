! main.f95
! Date: Jan. 12, 2018

!************************************************************************
!************************ Variables and functions ***********************
!************************************************************************

module init

  use vectors
  implicit none
  save

  real(dp), parameter :: pi = 4 * atan(1.0_dp)
  real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
  logical, parameter  :: debug = .false.
  integer,  parameter :: Nradius=512 ! num of grids on disk radius
  integer, parameter :: Nphi=256 ! devide half a disk into Nphi parts polarly
  integer,  parameter :: Ntheta0max=100, Ncos=100, Ndegree=90
  real(dp), dimension(Nradius) :: radius, height

end module init


program main

  use init
  implicit none

  integer  :: shape = 1
  real(dp) :: diskRadius=4000.0, para1=0.0, h0=0.0, r0=1.0
  real(dp) :: rNS = huge(one)
  real(dp) :: rg
  character(len=10) :: los = 'degree'
  real(dp) :: directAni, diskAni, reflectedAni
  real(dp) :: theta0, dtheta0, costheta0, radia_tot
  integer  :: i
  logical  :: lambertian=.true. ! If true, use Lambert's cosine law
  logical  :: GR=.false.
  real(dp), dimension(Nradius) :: mid, delh, gap, hypotenuse, &
      height_m, radia_delphi
  real(dp), dimension(Nradius) :: radia_flux

  namelist / GEOMETRY / shape, para1, diskRadius, r0, h0, los
  namelist / PHYSICS / lambertian, GR, rNS

  open (5, file='./input.nml') !, status='old')
  read (5, nml=GEOMETRY)
  read (5, nml=PHYSICS)
  close (5)

  rg = one / rNS ! normalize length units to rNS

  write (*,*) 'Running with the following parameters:'
  write (*,3) '     shape = ', shape
  write (*,9) '     para1 = ', para1
  write (*,4) 'lambertian = ', lambertian
  write (*,4) '        GR = ', GR
  if (GR) write (*,8) "rg = ", rg
  write (*,9) ' radius(inner) = ', r0
  write (*,9) 'radius(outer) = ', diskRadius
  write (*,9) ' height(inner) = ', h0
8 format (1X, A20, F20.10, " rNS")
9 format (1X, A20, F20.10)
3 format (1X, A20, I20)
4 format (1X, A20, L20)

  directAni = zero
  diskAni = zero
  reflectedAni = zero

  write (*,*)
  write (*,*) "Calculating radiation onto the disk..."
  call radiation(shape, diskRadius, para1, h0, r0, rg, lambertian, GR, &
      radia_delphi, radia_flux, radia_tot, mid, delh, gap, hypotenuse, &
      height_m)

  write (*,9) 'RF_intrinsic = ', radia_tot / pi

  write (*,*)
  write (*,*) "Calculating anisotropy factors..."
  if (los == 'cos') then
    write (*,10) 'cos(theta)', 'direct ani', 'reflected ani', 'reflected/direct'
    do i = 0, Ncos
      costheta0 = real(i, dp) / real(Ncos, dp)
      theta0 = acos(costheta0)
      call direct(r0, h0, radius(Nradius), height(Nradius), theta0, lambertian, &
          GR, rg, directAni)
      call disk(mid, delh, gap, hypotenuse, height_m, radia_delphi, reflectedAni, &
          theta0, lambertian, GR, rg)
      !      radia_delphi = zero
      !      call disk(mid, delh, gap, hypotenuse, height_m, radia_delphi, diskAni, &
      !          theta0_rad, Ntheta0)
      write (*,6) costheta0, directAni, reflectedAni, reflectedAni / directAni
    end do
  else
    dtheta0 = pi / 2.0_dp / Ndegree
    theta0 = zero
    write (*,10) 'degree', 'direct ani', 'reflected ani', 'reflected/direct'
    do i = 0, Ndegree
      call direct(r0, h0, radius(Nradius), height(Nradius), theta0, lambertian, &
          GR, rg, directAni)
      call disk(mid, delh, gap, hypotenuse, height_m, radia_delphi, reflectedAni, &
          theta0, lambertian, GR, rg)
!      radia_delphi = zero
!      call disk(mid, delh, gap, hypotenuse, height_m, radia_delphi, diskAni, &
!          theta0_rad, Ntheta0)
      write (*,7) i, directAni, reflectedAni, reflectedAni / directAni
      theta0 = theta0 + dtheta0
    end do
  end if

!  if (los == 'cos') then
!    write (*,10) 'radian', 'direct ani', 'reflected ani', 'reflected/direct'
!    do i = 0, Ntheta0
!      write (*,6) theta0_rad(i), directAni(i), reflectedAni(i), &
!          reflectedAni(i) / directAni(i)
!    end do
!  else
!    write (*,10) 'degree', 'direct ani', 'reflected ani', 'reflected/direct'
!    do i = 0, Ntheta0
!      write (*,7) i, directAni(i), reflectedAni(i), &
!          reflectedAni(i) / directAni(i)
!    end do
!  end if

6 format (1X, F14.2, F18.12, F18.12, F18.12)
7 format (1X, I14, F18.12, F18.12, F18.12)
10 format (1X, A14, A18, A18, A18)

end program main


module functions
  use init
  implicit none

contains

  ! Here we define functions of disk_height and radiative tranfer

  !***************************** disk_height *****************************
  elemental real(dp) function disk_height(r, shape, para1, r0)
    ! Input a list of radial coordinates
    ! output a list of coresponding disk heights
    implicit none
    integer,  intent(in) :: shape
    real(dp), intent(in) :: para1, r0
    real(dp), intent(in) :: r
    ! real(dp), dimension(size(r)), intent(out) :: disk_height

    select case(shape)
    case(1)
      disk_height = zero
    case(2)
      disk_height = para1 * r
    case(3)
      ! disk_height = para1 * (r - 1)
      disk_height = para1 * (r - r0)
    end select
  end function disk_height

  !*************************** Distri_star *******************************
  elemental real(dp) function distrib_star(cosBeta, lambertian)
    implicit none
    real(dp), intent(in) :: cosBeta
    logical,  intent(in) :: lambertian
    ! Normalization factor of (1 + 2.06 nu)
    real(dp), parameter :: norm_factor = 7.45604656451977595_dp

    if (lambertian) then
       distrib_star = cosBeta / pi
    else
       distrib_star = (one + 2.06 * cosBeta) * cosBeta
       distrib_star = distrib_star / norm_factor  ! Normalization
    end if
  end function distrib_star

  !*************************** Light bending **************************
  elemental real(dp) function getCosAlpha(mu, R, rg)
    implicit none
    real(dp), intent(in) :: mu, R, rg

    !    getCosAlpha = 1.0_dp - (1.0_dp - mu) * (1.0 - rg / R)
    getCosAlpha = mu * (1.0_dp - rg / R) + rg / R
    !    cosAlpha = mu
  end function getCosAlpha

  elemental real(dp) function dFdS(cosAlpha, R, lambertian, rg)
    implicit none
    real(dp), intent(in) :: cosAlpha, R, rg
    logical, intent(in) :: lambertian

    dFdS = (one - rg/R)**2 * distrib_star(cosAlpha, lambertian)
    ! dFdS = (1 - rg/R)**2 * cosAlpha / pi
  end function dFdS

  elemental real(dp) function rPsi(mu, b, rg)
    ! returns r(psi) in units of r_g
    implicit none

    real(dp), intent(in) :: mu, b, rg
    real(dp) :: part1, part2

    part1 = rg**2 * (1.0_dp - mu)**2 / (4 * (1.0_dp + mu)**2) + b**2 / &
         (1.0_dp - mu**2)
    part2 = rg * (1.0_dp - mu) / (2 * (1.0_dp + mu))
    rPsi = sqrt(part1) - part2
  end function rPsi

end module functions


subroutine radiation(shape, diskRadius, para1, h0, r0, rg, lambertian, GR, &
    radia_delphi, radia_flux, radia_tot, mid, delh, gap, hypotenuse, &
    height_m)
  ! Calculate the radiation distribution onto the disk surface

  use functions
  implicit none

  ! Data dictionary: declare calling parameter types & definitions
  integer, intent(in)  :: shape ! =1 for flat disk; =2 for inclined disk
  real(dp), intent(in) :: diskRadius, para1, h0, r0, rg
  logical, intent(in) :: lambertian, GR
  real(dp), intent(out), dimension(Nradius) :: radia_delphi, radia_flux
  real(dp), intent(out), dimension(Nradius) :: mid, delh, gap, hypotenuse, &
      height_m
  real(dp), intent(out) :: radia_tot

  ! Data dictionary: local variables
  real(dp), parameter :: del_phi = pi / 2 / Nphi ! On disk plane
  real(dp) :: delTheta, del_epsilon
  integer :: Nepsilon = 256 ! epsilon is the azimuzal angle on stellar surface
  integer :: Ntheta = 256 ! theta is the polar angle on stellar surface
  real(dp) :: logradius, logRadiusI
  real(dp), dimension(Nradius) :: area, surfArea ! Area inside a column,
  ! Angles on the equatorial plane or the disk surface
  real(dp)  ::  theta, phi, thetaMax, sint, cost, del_omega, cos_p
  real(dp), dimension(Nradius) :: clt, TQ, solid_an, flux
  real(dp), dimension(Nradius) :: cosBeta, cosGamma
  integer :: i

  if (shape < 1 .or. shape >3) then
    write (*,*) 'Disk shape not defined'
    stop                  ! Need to improve
  end if

  delTheta = pi / 2.0_dp / Ntheta
  del_epsilon = pi / 2.0_dp / Nepsilon

  ! Make a radius list
  logradius = log10(diskRadius)
  logRadiusI = log10(r0)
  do i = 1, Nradius
    radius(i) = logRadiusI + real(i - 1, dp) / real(Nradius, dp) * &
        (logradius - logRadiusI)
  end do
  radius = 10.**radius

  ! Others paras
  gap(1:Nradius-1) = radius(2:Nradius) - radius(1:Nradius-1)
  gap(Nradius) = 2 * gap(Nradius-1) - gap(Nradius-2)

  mid = radius + gap / 2.0_dp
  area = mid * gap * sin(del_phi) ! area of trapezoids
  mid = mid * cos(del_phi / 2.0_dp)

  height = disk_height(radius, shape, para1, r0) ! h
  height_m = disk_height(mid, shape, para1, r0)

  delh(:Nradius - 1) = height(2:) - height(:Nradius - 1)
  delh(Nradius) = 2.0_dp * delh(Nradius - 1) - delh(Nradius - 2)
  hypotenuse = sqrt(gap*gap + delh*delh)
  surfArea = area * hypotenuse / gap

  thetaMax = acos(1 / max(1.0_dp, height_m(Nradius)))  ! Eqn(15)

!  if (GR) then
!    write (*,*) "Current version does not support GR in calculating radiation"
!    write (*,*) "onto disk. Using non-GR scheme. The GR scheme for directed"
!    write (*,*) "and reflected anisotropy is working properly."
!  end if

  clt = zero
  theta = - thetaMax + delTheta / 2
  thetaloop: do
    sint = sin(theta)
    cost = cos(theta)
    ! area on the surface
    del_omega = sint * delTheta * del_epsilon
    phi = del_phi / 2
    do
      ! do phi = del_phi / 2, pi / 2, del_phi
      ! see figure 4 in tex/Essay/accretion_disk Equations are from
      ! tex/Essay/accretion_disk/numerical_approach_to_RF.pdf
      cos_p = cos(phi)
      TQ = sqrt(mid**2 - 2 * sint * cos_p * mid + &
          height_m**2 - 2 * height_m * cost + 1)
      ! Eqn(6). middle point Eqn(7). middle point
      cosBeta = sint * cos_p * mid + height_m * cost - 1
      cosBeta = cosBeta / TQ
      ! Eqn(14)
      cosGamma = abs(delh * (mid - cos_p * sint) - &
          gap * (height_m - cost)) / (hypotenuse * TQ)
      ! Eqn(15) solid angle of the area on the disk
      solid_an = surfArea * cosGamma / TQ**2

      where (mid**2 + height_m**2 < 1 + TQ**2) solid_an = zero
      where (-delh * (mid - cos_p * sint) + gap * (height_m - cost) > 0) &
          solid_an = zero

      flux = del_omega * distrib_star(cosBeta, lambertian) * solid_an
      clt = clt + flux
      ! Eqn(5) flux onto every small area
      ! pi is there because int_{0}^{pi/2} cos(theta) sin(theta)
      ! d theta d phi = pi
      ! clt += flux

      phi = phi + del_phi
      if (phi > pi / 2) exit
    end do
    theta = theta + delTheta
    if (theta > pi / 2) exit
  end do thetaloop

  radia_delphi = clt * 2 * Nepsilon / Nphi ! radia in r*dphi*dr
  radia_tot = sum(radia_delphi) * 4 * Nphi ! the intrinsic reflection fraction
  radia_flux = radia_delphi / area ! flux as a function of radius

end subroutine radiation

!************************************************************************
!***************************** Main program *****************************
!************************************************************************


!*************************** Anisotropy factors *******************************

subroutine direct(r0, h0, rEnd, hEnd, theta0, lambertian, GR, rg, directAni)

  ! Dependancies: This subroutine uses another subroutine functions, which
  ! uses subroutine init, which uses subroutine vectors

  ! Inputs:
  ! -------
  ! r0: real(dp) The radius of disk inner edge
  ! h0: real(dp) The height of disk inner edge
  ! rEnd: real(dp) The radius of disk outer edge
  ! hEnd: real(dp) The height of disk outer edge
  ! theta0_rad: real(dp)(0:Ntheta0max)
  ! labertian: logical
  ! GR: logical
  ! rg: real(dp)

  ! Outputs:
  ! --------
  ! directAni: real(dp)(0:Ntheta0max)

  use functions
  implicit none

  !  real(dp), parameter :: pi = 4 * atan(1.0_dp)
  integer, parameter :: num = 100
  real(dp), intent(in) :: rEnd, hEnd ! the radius and height of disk outer edge
  real(dp), intent(in) :: h0, r0, rg
  real(dp), intent(in) :: theta0 ! theta0_rad(0:Ntheta0max)
  real(dp), intent(out) :: directAni
  logical, intent(in) :: lambertian, GR

  integer :: i, j, k
  real(dp), dimension(num*2) :: cosThetas, sinThetas, cos_psi
  real(dp), dimension(num*2) :: localFlux, blocked, blockedI
  real(dp) :: totflux, dCosTheta, dPhi, phi, sinPhi, cosPhi, dS
  !       real(dp) :: directAni
  real(dp) :: cost0, sint0, tant0
  real(dp), dimension(num*2) :: cA, cB, cC
  real(dp), dimension(num*2) :: delta, cosPhi1, sinPhi1, cosPsi1, sinPsi1
  real(dp), dimension(num*2) :: cosAlpha, sinAlpha1, b1, r1
  real(dp) :: edge2center, triProduct
  type (vector) :: vec0, vecP, vecQ, vec0crossP, vec0crossQ

  dPhi   = pi / num

  ! cosThetas  = (/ (i, i = 1, num) /) / real(num, dp)
  cosThetas  = (/ (i, i = -num, num - 1) /) / real(num, dp) + &
      0.5 / real(num, dp)
  sinThetas  = sqrt(1 - cosThetas**2)
  dCosTheta = 1.0_dp / num
  edge2center = sqrt(r0**2 + h0**2)

!  ! loop through all los angles
!  l1: do j = 0, Ntheta0
!    theta0 = theta0_rad(j)
    cost0  = cos(theta0)
    sint0  = sin(theta0)
    tant0  = tan(theta0)

    vec0 = (/ sint0, zero, cost0 /) ! los unit vector

    totflux = zero
    phi = dPhi / 2.0
    l2:      do
      sinPhi = sin(phi)
      cosPhi = cos(phi)

      ! <surface normal, los>
      cos_psi = sinThetas * cosPhi * sint0 + cosThetas * cost0
      dS = dCosTheta * dPhi
      if (.not. GR) then
        localFlux = dS * distrib_star(cos_psi, lambertian)
        where (cos_psi < zero) localFlux = zero
      else
        cosAlpha = getCosAlpha(cos_psi, one, rg)
        localFlux = dS * dFdS(cosAlpha, one, lambertian, rg)
        where (cosAlpha < zero) localFlux = zero
      end if

      if (.not. GR) then

        ! blocked by the outer disk
        blocked = sinThetas*cosPhi - (cosThetas - hEnd) * tant0 - &
            sqrt(rEnd**2 - (sinThetas * sinPhi)**2)
        where (blocked > zero) localFlux = zero

        ! blocked by the inner disk
        blockedI = sinThetas * cosPhi - (cosThetas - h0) * &
            tant0 - sqrt(r0**2 - (sinThetas * sinPhi)**2)
        where (blockedI > zero) localFlux = zero

      else

        ! <p1 Check blocking by inner disk>
        ! <p2 Solve for intersections with disk inner edge>
        cA = r0 * sinThetas * sinPhi * cost0
        cB = r0 * (-sinThetas * cosPhi * cost0 + cosThetas * sint0)
        cC = -h0 * sinThetas * sinPhi * sint0
        delta = cB**2 * (cA**2 + cB**2 - cC**2) ! delta of quadratic eq
        cosPhi1 = (-cA * cC + sqrt(delta)) / (cA**2 + cB**2) ! on the ring
        sinPhi1 = sqrt(1.0_dp - cosPhi1**2)
        cosPsi1 = sint0 * r0 * cosPhi1 + cost0 * h0 ! <los, ring>
        sinPsi1 = sqrt(1 - cosPsi1**2)
        sinAlpha1 = sqrt(1 - cosAlpha**2)
        b1 = one * sinAlpha1 / sqrt(one - rg/one) ! impact factor
        r1 = rPsi(cosPsi1, b1, rg) ! r(psi) at the direction of the ring
        ! </p2>

        do k = 1, 2*num

          ! <p2 If no real root, the los is totally blocked by the disk>
          if (delta(k) < zero) then
            localFlux(k) = zero
            cycle
          end if
          ! </p2>

          ! <p2 Otherwise, check whether or not the los is blocked>
          vecP = (/ sinThetas(k) * cosPhi, sinThetas(k) * sinPhi, &
              cosThetas(k) /)
          vecQ = (/ r0 * cosPhi1(k), r0 * sinPhi1(k), h0 /) ! one case
          vec0crossP = vec0 * vecP ! cross product
          vec0crossQ = vec0 * vecQ ! cross product
          triProduct = vec0crossP .dot. vecQ

          ! If vec0, vecP, vecQ are not in one plane, flip Q_y
          if (abs(triProduct) > real(1.0e-9, dp)) then
            vecQ%y = -vecQ%y
            vec0crossQ = vec0 * vecQ
            triProduct = vec0crossP .dot. vecQ
            if (abs(triProduct) > real(1.0e-9, dp)) then
              print *, "O, P, Q fails to be in one plane. Algorithm error?"
              call exit(-1)
            end if
          end if

          ! Now that vec0, vecP, and vecQ are in one plane, if P and Q
          ! are in different side of vec0, it's not blocked
          triProduct = vec0crossP .dot. vec0crossQ
          if (triProduct < zero) cycle

          ! Now that P and Q are in same side of vec0, check whether
          ! the photon path is blocked by the disk
          if (r1(k) > edge2center) localFlux(k) = zero
          ! </p2>

        end do
        ! </p1>

        ! <p1 Checking blocking by outer disk>
        ! TODO
        ! </p1>
      end if

      totflux = totflux + sum(localFlux)

      phi = phi + dPhi
      if (phi > pi) exit
    end do l2

    directAni = totflux * 2
!    directAni(j) = totflux * 2
!  end do l1

end subroutine direct


subroutine disk(mid, delh, gap, hypotenuse, height_m, radia_delphi, diskAni, &
    theta0, lambertian, GR, rg)

  use functions
  implicit none

  ! Data dictionary: input and output parameters
  real(dp), intent(in), dimension(Nradius) :: mid, delh, gap, hypotenuse, &
      height_m, radia_delphi
  real(dp), intent(in) :: theta0, rg
  logical, intent(in) :: lambertian, GR
  real(dp), intent(out) :: diskAni

  real(dp), parameter :: del_phi = pi / 2 / Nphi ! On disk plane
  real(dp), dimension(Nradius) :: cos_psi2, disk_blocked !, localFlux
  real(dp), dimension(Nradius) :: unblocked_by_star
  real(dp), dimension(Nradius) :: cosPsi, sinPsi, r0, cosAlpha, sinAlpha, b
  real(dp) :: phi0 = zero, cost0, sint0, tant0, cosPhi, sinPhi
  real(dp) :: totflux, phi
  !      type (vector) :: vec0, vecP, vecQ, vec0crossP, vec0crossQ
  real(dp) :: bCritical, epsilon, cosZeta
  integer :: j
  type(vector) :: vec0, vecP, vecPPrime, vecE, vecn
  real(dp), dimension(Nradius) :: sinPsiMinusAlpha, cosPsiMinusAlpha
  integer :: i

  bCritical = sqrt(one / (one - rg)) ! (b/Rc)^2 = 1 / (1 - rg / Rc)
  r0 = sqrt(mid**2 + height_m**2)

!  l1 :  do j = 0, Ntheta0
!    theta0 = theta0_rad(j)
    cost0  = cos(theta0)
    sint0  = sin(theta0)
    tant0  = tan(theta0)

    vec0 = (/ sint0, zero, cost0 /) ! los unit vector

    totflux = zero
    phi = del_phi / 2
    l1 : do
      cosPhi = cos(phi)
      sinPhi = sin(phi)

      ! surface normal
      epsilon = zero
      vecn = (/ -sin(epsilon) * cosPhi, &
          -sin(epsilon) * sinPhi, &
          cos(epsilon) /)

      if (debug .and. j == 45 .and. phi > 3.0) then
        print *, 'debug'
      end if

      if (.not. GR) then
        ! <los, surface normal>
        cos_psi2 = (-(cosPhi * sint0 * cos(phi0) + sinPhi * sint0 * &
            sin(phi0)) * delh + gap * cost0) / hypotenuse

        disk_blocked = mid**2 - 2 * mid * cosPhi * (height_m - &
            height_m(Nradius)) * tant0 + (height_m - &
            height_m(Nradius))**2 * tant0**2 - mid(Nradius)**2
        unblocked_by_star = mid**2 * (cosPhi**2 * cost0**2 + sinPhi**2) - &
            2 *  mid * cosPhi * cost0 * sint0 * height_m + &
            height_m**2 * sint0**2 - 1
        do i = 1, Nradius
          if (.not. disk_blocked(i) < zero) cycle
          if (.not. cos_psi2(i) > zero) cycle
          if (phi > pi / 2) then
            if (.not. unblocked_by_star(i) > zero) cycle
          end if
          totflux = totflux + radia_delphi(i) * distrib_star(&
              cos_psi2(i), lambertian)
        end do
      else
        !            if (.not. shape == 1) then
        if (.not. all(abs(height_m) < 1e-14)) then
          print *, ""
          print *, "Fail! GR = .true. only supports flat disk"
          call exit(-1)
        end if
        ! psi and alpha
        cosPsi = (sint0 * mid * cosPhi + cost0 * height_m) / r0
        sinPsi = sqrt(1 - cosPsi**2)
        cosAlpha = getCosAlpha(cosPsi, r0, rg)
        sinAlpha = sqrt(1 - cosAlpha**2)
        cosPsiMinusAlpha = cosPsi * cosAlpha + sinPsi * sinAlpha
        sinPsiMinusAlpha = sinPsi * cosAlpha - cosPsi * sinAlpha
        b = r0 * sinAlpha / sqrt(one - rg / r0) ! impact factor

        do i = 1, Nradius
          ! blocked by star
          if (cosAlpha(i) < zero .and. b(i) < bCritical) cycle
          ! a mesh grid on disk surface
          vecP = (/ mid(i) * cosPhi, mid(i) * sinPhi, height_m(i) /)
          vecPPrime = (vec0 * vecP) * vec0
          vecE = cosPsiMinusAlpha(i) * vec0 + sinPsiMinusAlpha(i) / &
              sqrt(vecPPrime .dot. vecPPrime) * vecPPrime
          cosZeta = (vecE .dot. vecn) / sqrt(vecE .dot. vecE)
          if (cosZeta < zero) cycle ! obtuse emergent angle
          totflux = totflux + radia_delphi(i) * dFdS(cosZeta, r0(i), &
              lambertian, rg)
        end do

      end if

      phi = phi + del_phi
      if (phi > pi) exit
    end do l1

    ! do
    !   ! # account for the shadow of both the disk and the star.
    !   ! # For phi form pi/2 to pi
    !   cosPhi = cos(phi)
    !   sinPhi = sin(phi)

    !   ! shadow of the disk
    !   cos_psi2 = (-(cosPhi * sint0 * cos(phi0) + sinPhi * sint0 * &
    !       sin(phi0)) * delh + gap * cost0) / hypotenuse

    !   blocked_by_disk = mid**2 - 2 * mid * cosPhi * (height_m - &
    !       hEnd) * tant0 + (height_m - hEnd)**3 * &
    !       tant0**2 - r1**2

    !   ! F_mid < zero is visible

    !   ! shadow of the star
    !   unblocked_by_star = mid**2 * (cosPhi**2 * cost0**2 + sinPhi**2) - &
    !       2 *  mid * cosPhi * cost0 * sint0 * height_m + &
    !       height_m**2 * sint0**2 - 1
    !   ! G_mid > zero is visible
    !   do i = 1, Nradius
    !     if (blocked_by_disk(i) < zero .and. unblocked_by_star(i) > zero &
    !         .and. cos_psi2(i) > zero) then
    !       totflux = totflux + radia_delphi(i) * &
    !           distrib_star(cos_psi2(i), lambertian)
    !     end if
    !   end do
    !   phi = phi + del_phi
    !   if (phi > pi) exit
    ! end do

    diskAni = totflux * 2 ! diskAni must be inout
!  end do l1

end subroutine disk
