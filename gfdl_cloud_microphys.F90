!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Cloud Microphysics.
!*
!* The GFDL Cloud Microphysics is free software: you can
!* redistribute it and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The GFDL Cloud Microphysics is distributed in the hope it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the GFDL Cloud Microphysics.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'gfdl_cloud_microphys' contains the full GFDL cloud
!! microphysics \cite chen2013seasonal.
!>@details The module is paired with 'fv_cmp', which performs the "fast"
!! processes
!>author Shian-Jiann Lin, Linjiong Zhou

! =======================================================================
! cloud micro - physics package for gfdl global cloud resolving model
! the algorithms are originally derived from lin et al 1983. most of the
! key elements have been simplified / improved. this code at this stage
! bears little to no similarity to the original lin mp in zetac.
! therefore, it is best to be called gfdl micro - physics (gfdl mp) .
! developer: shian-jiann lin, linjiong zhou
! =======================================================================

module gfdl2_cloud_microphys_mod

  use omp_lib

  implicit none
  private

  public gfdl_cloud_microphys_driver, gfdl_cloud_microphys_init

  real :: missing_value = - 1.e10

  logical :: module_is_initialized = .false.
  logical :: qsmith_tables_initialized = .false.

  character (len = 17) :: mod_name = 'gfdl_cloud_microphys'

  real, parameter :: MAPL_PI = 4.0 * atan(1.0)
  real, parameter :: grav = 9.80665 !< gfs: acceleration due to gravity
  real, parameter :: rdgas = 287.05 !< gfs: gas constant for dry air
  real, parameter :: rvgas = 461.50 !< gfs: gas constant for water vapor
  real, parameter :: cp_air = 1004.6 !< gfs: heat capacity of dry air at constant pressure
  real, parameter :: hlv = 2.5e6 !< gfs: latent heat of evaporation
  real, parameter :: hlf = 3.3358e5 !< gfs: latent heat of fusion
  real, parameter :: pi = 3.1415926535897931 !< gfs: ratio of circle circumference to diameter
  ! real, parameter :: cp_air = rdgas * 7. / 2. ! 1004.675, heat capacity of dry air at constant pressure
  real, parameter :: cp_vap = 4.0 * rvgas !< 1846.0, heat capacity of water vapore at constnat pressure
  ! real, parameter :: cv_air = 717.56 ! satoh value
  real, parameter :: cv_air = cp_air - rdgas !< 717.55, heat capacity of dry air at constant volume
  ! real, parameter :: cv_vap = 1410.0 ! emanuel value
  real, parameter :: cv_vap = 3.0 * rvgas !< 1384.5, heat capacity of water vapor at constant volume

  ! the following two are from emanuel's book "atmospheric convection"
  ! real, parameter :: c_ice = 2106.0 ! heat capacity of ice at 0 deg c: c = c_ice + 7.3 * (t - tice)
  ! real, parameter :: c_liq = 4190.0 ! heat capacity of water at 0 deg c

  real, parameter :: c_ice = 1972.0 !< gfdl: heat capacity of ice at - 15 deg c
  real, parameter :: c_liq = 4185.5 !< gfdl: heat capacity of water at 15 deg c
  ! real, parameter :: c_liq = 4218.0 ! ifs: heat capacity of liquid at 0 deg c

  real, parameter :: eps = rdgas / rvgas ! 0.6219934995
  real, parameter :: zvir = rvgas / rdgas - 1. !< 0.6077338443

  real, parameter :: t_ice = 273.16 !< freezing temperature
  real, parameter :: table_ice = 273.16 !< freezing point for qs table

  ! real, parameter :: e00 = 610.71 ! gfdl: saturation vapor pressure at 0 deg c
  real, parameter :: e00 = 611.21 !< ifs: saturation vapor pressure at 0 deg c

  real, parameter :: dc_vap = cp_vap - c_liq !< - 2339.5, isobaric heating / cooling
  real, parameter :: dc_ice = c_liq - c_ice !< 2213.5, isobaric heating / colling

  real, parameter :: hlv0 = hlv !< gfs: evaporation latent heat coefficient at 0 deg c
  ! real, parameter :: hlv0 = 2.501e6 ! emanuel appendix - 2
  real, parameter :: hlf0 = hlf !< gfs: fussion latent heat coefficient at 0 deg c
  ! real, parameter :: hlf0 = 3.337e5 ! emanuel

  real, parameter :: lv0 = hlv0 - dc_vap * t_ice!< 3.13905782e6, evaporation latent heat coefficient at 0 deg k
  real, parameter :: li00 = hlf0 - dc_ice * t_ice!< - 2.7105966e5, fusion latent heat coefficient at 0 deg k

  real, parameter :: d2ice = dc_vap + dc_ice !< - 126, isobaric heating / cooling
  real, parameter :: li2 = lv0 + li00 !< 2.86799816e6, sublimation latent heat coefficient at 0 deg k

  real, parameter :: qpmin = 1.e-8  !< min value for suspended rain/snow/liquid/ice precip
  real, parameter :: qvmin = 1.e-20 !< min value for water vapor (treated as zero)
  real, parameter :: qcmin = 1.e-12 !< min value for cloud condensates

  real, parameter :: vr_min = 1.e-3 !< min fall speed for rain
  real, parameter :: vf_min = 1.e-5 !< min fall speed for cloud ice, snow, graupel

  real, parameter :: dz_min = 1.e-2 ! use for correcting flipped height

  real, parameter :: sfcrho = 1.2 !< surface air density
  real, parameter :: rhor = 1.e3 !< density of rain water, lin83

  real, parameter :: rc = (4. / 3.) * pi * rhor

  integer, parameter :: TABLE_LENGTH = 2621

  real :: cracs, csacr, cgacr, cgacs, csacw, craci, csaci, cgacw, cgaci, cracw !< constants for accretions
  real :: acco (3, 4) !< constants for accretions
  real :: cssub (5), cgsub (5), crevp (5), cgfr (2), csmlt (5), cgmlt (5)

  real :: es0, ces0
  real :: pie, rgrav
  real :: c_air, c_vap

  real :: lati, latv, lats, lat2, lcp, icp, tcp !< used in bigg mechanism and wet bulk

  real :: d0_vap !< the same as dc_vap, except that cp_vap can be cp_vap or cv_vap
  real :: lv00 !< the same as lv0, except that cp_vap can be cp_vap or cv_vap

  ! cloud microphysics switchers

  integer :: icloud_f = 0 !< cloud scheme
  integer :: irain_f = 0 !< cloud water to rain auto conversion scheme

  logical :: de_ice = .false. !< to prevent excessive build - up of cloud ice from external sources
  logical :: sedi_transport = .false. !< transport of momentum in sedimentation
  logical :: do_sedi_w = .false. !< transport of vertical motion in sedimentation
  logical :: do_sedi_heat = .false. !< transport of heat in sedimentation
  logical :: prog_ccn = .false. !< do prognostic ccn (yi ming's method)
  logical :: do_bigg = .false. !< do bigg mechanism freezing of supercooled liquid on aerosol nuclei
  logical :: do_evap = .false. !< do evaporation
  logical :: do_subl = .false. !< do sublimation
  logical :: do_qa = .false. !< do inline cloud fraction (WMP: in FV3 dynamics)
  logical :: preciprad = .true. !< consider precipitates in cloud fraciton calculation
  logical :: fix_negative = .false. !< fix negative water species
  logical :: do_setup = .true. !< setup constants and parameters
  logical :: p_nonhydro = .false. !< perform hydrosatic adjustment on air density

  real, dimension(TABLE_LENGTH) :: table, table2, table3, tablew
  real, dimension(TABLE_LENGTH) :: des, des2, des3, desw

  logical :: tables_are_initialized = .false.

  ! logical :: root_proc
  ! integer :: id_rh, id_vtr, id_vts, id_vtg, id_vti, id_rain, id_snow, id_graupel, &
  ! id_ice, id_prec, id_cond, id_var, id_droplets
  ! integer :: gfdl_mp_clock ! clock for timing of driver routine

  real :: dt_fr = 8. !< epsilon on homogeneous freezing of cloud water at t_wfr + dt_fr
  ! minimum temperature water can exist (moore & molinero nov. 2011, nature)
  ! dt_fr can be considered as the error bar

  real, parameter :: p_min = 100. !< minimum pressure (pascal) for mp to operate

  ! slj, the following parameters are for cloud - resolving resolution: 1 - 5 km

  ! qi0_crt = 0.8e-4
  ! qs0_crt = 0.6e-3
  ! c_psaci = 0.1
  ! c_pgacs = 0.1
  ! c_pgaci = 0.05

  ! -----------------------------------------------------------------------
  !> namelist parameters
  ! -----------------------------------------------------------------------

  real :: cld_min = 0.05 !< minimum cloud fraction
  real :: tice = 273.16 !< set tice = 165. to trun off ice - phase phys (kessler emulator)

  real :: log_10 = log (10.)
  real :: tice0 = 273.16 - 0.01
  real, parameter :: t_wfr = 273.16 - 40.0 ! supercooled water can exist down to - 40 c, which is the "absolute"

  real :: t_min = 178. !< min temp to freeze - dry all water vapor
  real :: t_sub = 184. !< min temp for sublimation of cloud ice
  real :: mp_time = 150. !< maximum micro - physics time step (sec)

  ! relative humidity increment

  real :: rh_inc = 0.25 !< rh increment for complete evaporation of cloud water and cloud ice
  real :: rh_inr = 0.25 !< rh increment for minimum evaporation of rain
  real :: rh_ins = 0.25 !< rh increment for sublimation of snow

  ! conversion time scale

  real :: tau_r2g = 900. !< rain freezing during fast_sat
  real :: tau_smlt = 900. !< snow melting
  real :: tau_g2r = 600. !< graupel melting to rain
  real :: tau_imlt = 600. !< cloud ice melting
  real :: tau_i2s = 1000. !< cloud ice to snow auto - conversion
  real :: tau_l2r = 900. !< cloud water to rain auto - conversion
  real :: tau_v2l = 150. !< water vapor to cloud water (condensation)
  real :: tau_l2v = 300. !< cloud water to water vapor (evaporation)
  real :: tau_i2v = 300. !< cloud ice to water vapor (sublimation)
  real :: tau_s2v = 600. !< snow sublimation
  real :: tau_v2s = 21600. !< snow deposition -- make it a slow process
  real :: tau_g2v = 900. !< graupel sublimation
  real :: tau_v2g = 21600. !< graupel deposition -- make it a slow process
  real :: tau_revp = 600. !< rain re-evaporation
  real :: tau_frz = 450. !, timescale for liquid-ice freezing
  ! horizontal subgrid variability

  real :: dw_land = 0.20 !< base value for subgrid deviation / variability over land
  real :: dw_ocean = 0.10 !< base value for ocean

  ! prescribed ccn

  real :: ccn_o = 90. !< ccn over ocean (cm^ - 3)
  real :: ccn_l = 270. !< ccn over land (cm^ - 3)

  real :: rthreshu =  7.0e-6 !< critical cloud drop radius (micro m)
  real :: rthreshs = 10.0e-6 !< critical cloud drop radius (micro m)

  real :: sat_adj0 = 0.90 !< adjustment factor (0: no, 1: full) during fast_sat_adj

  real :: qc_crt = 5.0e-8 !< mini condensate mixing ratio to allow partial cloudiness

  real :: qi_lim = 1. !< cloud ice limiter to prevent large ice build up

  real :: ql_mlt = 2.0e-3 !< max value of cloud water allowed from melted cloud ice
  real :: qs_mlt = 1.0e-6 !< max cloud water due to snow melt

  real :: ql_gen = 1.0e-3 !< max cloud water generation
  real :: qi_gen = 9.82679e-5 !< max cloud ice generation at -40 C

  ! cloud condensate upper bounds: "safety valves" for ql & qi

  real :: ql0_max = 2.0e-3 !< max cloud water value (auto converted to rain)
  real :: qi0_max = 1.0e-4 !< max cloud ice value (by other sources)

  real :: qi0_crt = 1.0e-4 !< cloud ice to snow autoconversion threshold (was 1.e-4)
  !! qi0_crt is highly dependent on horizontal resolution
  real :: qr0_crt = 1.0e-4 !< rain to snow or graupel / hail threshold
  !! lfo used * mixing ratio * = 1.e-4 (hail in lfo)
  real :: qs0_crt = 1.0e-3 !< snow to graupel density threshold (0.6e-3 in purdue lin scheme)

  real :: c_paut = 0.55 !< autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
  real :: c_psaci = 0.02 !< accretion: cloud ice to snow (was 0.1 in zetac)
  real :: c_piacr = 5.0 !< accretion: rain to ice:
  real :: c_cracw = 0.9 !< rain accretion efficiency
  real :: c_pgacs = 2.0e-3 !< snow to graupel "accretion" eff. (was 0.1 in zetac)
  real :: c_pgaci = 0.05   !<  ice to graupel "accretion" eff.

  ! decreasing clin to reduce csacw (so as to reduce cloud water --- > snow)

  real :: alin = 842.0 !< "a" in lin1983
  real :: clin = 4.8 !< "c" in lin 1983, 4.8 -- > 6. (to ehance ql -- > qs)

  ! fall velocity tuning constants:

  logical :: const_vi = .false. !< if .t. the constants are specified by v * _fac
  logical :: const_vs = .false. !< if .t. the constants are specified by v * _fac
  logical :: const_vg = .false. !< if .t. the constants are specified by v * _fac
  logical :: const_vr = .false. !< if .t. the constants are specified by v * _fac

  ! good values:

  real :: vi_fac = 1. !< if const_vi: 1 / 3
  real :: vs_fac = 1. !< if const_vs: 1.
  real :: vg_fac = 1. !< if const_vg: 2.
  real :: vr_fac = 1. !< if const_vr: 4.

  ! upper bounds of fall speed (with variable speed option)

  real :: vi_max = 1.0 !< max fall speed for ice
  real :: vs_max = 2.0 !< max fall speed for snow
  real :: vg_max = 12. !< max fall speed for graupel
  real :: vr_max = 12. !< max fall speed for rain

  ! cloud microphysics switchers

  logical :: fast_sat_adj = .false. !< has fast saturation adjustments
  logical :: z_slope_liq = .true. !< use linear mono slope for autocconversions
  logical :: z_slope_ice = .false. !< use linear mono slope for autocconversions
  logical :: use_ccn = .false. !< use input ccn when .T. else use ccn_o/ccn_l
  logical :: use_ppm = .false. !< use ppm fall scheme
  logical :: mono_prof = .true. !< perform terminal fall with mono ppm scheme
  logical :: mp_print = .false. !< cloud microphysics debugging printout

  ! real :: global_area = - 1.

  ! -----------------------------------------------------------------------
  ! namelist
  ! -----------------------------------------------------------------------

  namelist / gfdl_cloud_microphysics_nml /                                  &
       mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, &
       vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max,  &
       vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max,    &
       qi0_crt, qr0_crt, fast_sat_adj, rh_inc, rh_ins, rh_inr, const_vi,     &
       const_vs, const_vg, const_vr, use_ccn, rthreshu, rthreshs, ccn_l, ccn_o, qc_crt, &
       tau_g2v, tau_v2g, tau_s2v, tau_v2s, &
       tau_revp, tau_frz, do_bigg, do_evap, do_subl, &
       sat_adj0, c_piacr, tau_imlt, tau_v2l, tau_l2v, tau_i2v, &
       tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, c_pgaci,  &
       z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin,              &
       preciprad, cld_min, use_ppm, mono_prof,         &
       do_sedi_heat, sedi_transport, do_sedi_w, dt_fr, de_ice, icloud_f, irain_f, mp_print

  public                                                                    &
       mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, &
       vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max,  &
       vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max,    &
       qi0_crt, qr0_crt, fast_sat_adj, rh_inc, rh_ins, rh_inr, const_vi,     &
       const_vs, const_vg, const_vr, use_ccn, rthreshu, rthreshs, ccn_l, ccn_o, qc_crt, &
       tau_g2v, tau_v2g, tau_s2v, tau_v2s, &
       tau_revp, tau_frz, do_bigg, do_evap, do_subl, &
       sat_adj0, c_piacr, tau_imlt, tau_v2l, tau_l2v, tau_i2v, &
       tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, c_pgaci,  &
       z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin,              &
       preciprad, cld_min, use_ppm, mono_prof,         &
       do_sedi_heat, sedi_transport, do_sedi_w, dt_fr, de_ice, icloud_f, irain_f, mp_print

  !$omp declare target( &
  !$omp   des2, desw, table2, tablew, &

  !$omp   d0_vap, lv00, c_vap, c_air, tau_revp, &
  !$omp   tau_v2l, tau_l2v, tau_i2v, tau_s2v, tau_v2s, tau_g2v, &
  !$omp   tau_v2g, tau_frz, tau_imlt, tau_smlt, tau_i2s, tau_g2r, &
  !$omp   tice, tice0, rh_inc, rh_inr, t_min, do_qa, t_sub, do_evap, &
  !$omp   do_bigg, qi_lim, do_subl, preciprad, icloud_f, qc_crt, lat2, z_slope_ice, &
  !$omp   c_paut, prog_ccn, fix_negative, p_nonhydro, sedi_transport, ql_mlt, qs_mlt, qi0_crt, qs0_crt, &
  !$omp   const_vi, vi_fac, vi_max, const_vs, vs_fac, vs_max, const_vg, vg_fac, vg_max, const_vr, vr_fac, vr_max, &
  !$omp   do_sedi_w, use_ppm, mono_prof, rthreshs, rthreshu, irain_f, z_slope_liq, do_sedi_heat, &
  !$omp   ql0_max, dt_fr, sat_adj0, dw_land, dw_ocean, c_psaci, c_pgacs, &
  !$omp   ccn_l, ccn_o, c_cracw, use_ccn, de_ice, mp_time, &

  !$omp   ces0, cracs, cracw, &
  !$omp   csaci, csacr, csacw, cgaci, cgacr, cgacs, cgacw, &
  !$omp   cssub, crevp, csmlt, cgmlt, cgfr, acco)

contains

  ! -----------------------------------------------------------------------
  ! the driver of the gfdl cloud microphysics
  ! -----------------------------------------------------------------------

  !>@brief The subroutine 'gfdl_cloud_microphys_driver' executes the full GFDL
  !! cloud microphysics.
  subroutine gfdl_cloud_microphys_driver ( &
       qv, ql, qr, qi, qs, qg, qa, qn, &
       qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, &
       pt_dt, pt, w, &
       uin, vin, udt, vdt, &
       dz, delp, area, dt_in, &
       land, cnv_fraction, srf_type, eis,                                &
       rhcrit, anv_icefall, lsc_icefall,                                 &
       revap, isubl,                                                     &
       rain, snow, ice, graupel, &
       m2_rain, m2_sol, &
       hydrostatic, phys_hydrostatic, &
       iis, iie, jjs, jje, kks, kke, ktop, kbot)

    implicit none

    logical, intent (in) :: hydrostatic, phys_hydrostatic
    integer, intent (in) :: iis, iie, jjs, jje !< physics window
    integer, intent (in) :: kks, kke !< vertical dimension
    integer, intent (in) :: ktop, kbot !< vertical compute domain

    real, intent (in) :: dt_in !< physics time step

    real, intent (in), dimension (:, :) :: area !< cell area
    real, intent (in), dimension (:, :) :: land !< land fraction
    real, intent (in), dimension (:, :) :: cnv_fraction !< diagnosed convective fraction
    real, intent (in), dimension (:, :) :: srf_type
    real, intent (in), dimension (:, :) :: eis  !< estimated inversion strength
    real, intent (in), dimension (:, :, :) :: rhcrit

    real, intent (in) :: anv_icefall, lsc_icefall

    real, intent (in), dimension (:, :, :) :: delp, dz, uin, vin
    real, intent (in), dimension (:, :, :) :: pt, qv, ql, qr, qg, qa, qn

    real, intent (inout), dimension (:, :, :) :: qi, qs
    real, intent (inout), dimension (:, :, :) :: pt_dt, qa_dt, udt, vdt, w
    real, intent (inout), dimension (:, :, :) :: qv_dt, ql_dt, qr_dt
    real, intent (inout), dimension (:, :, :) :: qi_dt, qs_dt, qg_dt

    real, intent (out), dimension (:, :) :: rain, snow, ice, graupel
    real, intent (out), dimension (:, :, :) :: m2_rain, m2_sol ! Rain and Ice fluxes (Pa kg/kg)
    real, intent (out), dimension (:, :, :) :: revap ! Rain evaporation
    real, intent (out), dimension (:, :, :) :: isubl ! Ice sublimation

    ! logical :: used

    real :: mpdt, rdt, dts, convt, tot_prec

    integer :: i, j, k
    integer :: is, ie, js, je !< physics window
    integer :: ks, ke !< vertical dimension
    integer :: days, ntimes

    real, dimension (iie - iis + 1, jje - jjs + 1) :: prec_mp, prec1, cond, w_var, rh0

    real, dimension (iie - iis + 1, jje - jjs + 1, kke - kks + 1) :: vt_r, vt_s, vt_g, vt_i, qn2

    real :: allmax

    is = 1
    js = 1
    ks = 1
    ie = iie - iis + 1
    je = jje - jjs + 1
    ke = kke - kks + 1

    ! call mpp_clock_begin (gfdl_mp_clock)

    ! -----------------------------------------------------------------------
    ! define heat capacity of dry air and water vapor based on hydrostatical property
    ! -----------------------------------------------------------------------

    if (phys_hydrostatic .or. hydrostatic) then
       c_air = cp_air
       c_vap = cp_vap
       p_nonhydro = .false.
    else
       c_air = cv_air
       c_vap = cv_vap
       p_nonhydro = .true.
    endif

    d0_vap = c_vap - c_liq
    lv00 = hlv0 - d0_vap * t_ice

    if (hydrostatic) do_sedi_w = .false.

    ! -----------------------------------------------------------------------
    ! define latent heat coefficient used in wet bulb and bigg mechanism
    ! -----------------------------------------------------------------------

    latv = hlv
    lati = hlf
    lats = latv + lati
    lat2 = lats * lats

    !$omp target update to(c_air, c_vap, p_nonhydro, d0_vap, lv00, do_sedi_w, lat2)

    lcp = latv / cp_air
    icp = lati / cp_air
    tcp = (latv + lati) / cp_air

    ! tendency zero out for am moist processes should be done outside the driver

    ! -----------------------------------------------------------------------
    ! define cloud microphysics sub time step
    ! -----------------------------------------------------------------------

    mpdt   = min (dt_in, mp_time)
    rdt    = 1. / dt_in
    ntimes = nint (dt_in / mpdt)

    ! small time step:
    dts = dt_in / real (ntimes)

    ! call get_time (time, seconds, days)

    ! -----------------------------------------------------------------------
    ! initialize precipitation
    ! -----------------------------------------------------------------------

    do j = js, je
       do i = is, ie
          graupel (i, j) = 0.
          rain (i, j) = 0.
          snow (i, j) = 0.
          ice (i, j) = 0.
          cond (i, j) = 0.
       enddo
    enddo

    ! -----------------------------------------------------------------------
    ! major cloud microphysics
    ! -----------------------------------------------------------------------

    ! print *, 'gfdl_cloud_microphys_driver - calling mpdrv'
    call mpdrv ( &
         hydrostatic, uin, vin, w, delp, pt, qv, ql, qr, qi, qs, qg, &
         qa, qn, dz, is, ie, js, je, ks, ke, ktop, kbot, dt_in, ntimes, &
         rain(:, js:je), snow(:, js:je), graupel(:, js:je), ice(:, js:je), m2_rain, &
         m2_sol, cond(:, js:je), area(:, js:je), &
         land(:, js:je), cnv_fraction(:, js:je), srf_type(:, js:je), eis(:, js:je), &
         rhcrit, anv_icefall, lsc_icefall, &
         revap, isubl, &
         udt, vdt, pt_dt, &
         qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, w_var, vt_r, &
         vt_s, vt_g, vt_i, qn2)
    ! print *, 'gfdl_cloud_microphys_driver - completed mpdrv'

    ! -----------------------------------------------------------------------
    ! no clouds allowed above ktop
    ! -----------------------------------------------------------------------

    if (ks < ktop) then
       do k = ks, ktop
          do j = js, je
             do i = is, ie
                qa_dt (i, j, k) = 0.
             enddo
          enddo
       enddo
    endif

    ! -----------------------------------------------------------------------
    ! diagnostic output
    ! -----------------------------------------------------------------------

    ! if (id_vtr > 0) then
    ! used = send_data (id_vtr, vt_r, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_vts > 0) then
    ! used = send_data (id_vts, vt_s, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_vtg > 0) then
    ! used = send_data (id_vtg, vt_g, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_vti > 0) then
    ! used = send_data (id_vti, vt_i, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_droplets > 0) then
    ! used = send_data (id_droplets, qn2, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_var > 0) then
    ! used = send_data (id_var, w_var, time, is_in = iis, js_in = jjs)
    ! endif

    ! convert to mm / day

    convt = 86400. * rdt * rgrav
    do j = js, je
       do i = is, ie
          rain (i, j) = rain (i, j) * convt
          snow (i, j) = snow (i, j) * convt
          ice (i, j) = ice (i, j) * convt
          graupel (i, j) = graupel (i, j) * convt
          prec_mp (i, j) = rain (i, j) + snow (i, j) + ice (i, j) + graupel (i, j)
       enddo
    enddo

    ! if (id_cond > 0) then
    ! do j = js, je
    ! do i = is, ie
    ! cond (i, j) = cond (i, j) * rgrav
    ! enddo
    ! enddo
    ! used = send_data (id_cond, cond, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_snow > 0) then
    ! used = send_data (id_snow, snow, time, iis, jjs)
    ! used = send_data (id_snow, snow, time, is_in = iis, js_in = jjs)
    ! if (mp_print .and. seconds == 0) then
    ! tot_prec = g_sum (snow, is, ie, js, je, area, 1)
    ! if (root_proc) write (*, *) 'mean snow = ', tot_prec
    ! endif
    ! endif
    !
    ! if (id_graupel > 0) then
    ! used = send_data (id_graupel, graupel, time, iis, jjs)
    ! used = send_data (id_graupel, graupel, time, is_in = iis, js_in = jjs)
    ! if (mp_print .and. seconds == 0) then
    ! tot_prec = g_sum (graupel, is, ie, js, je, area, 1)
    ! if (root_proc) write (*, *) 'mean graupel = ', tot_prec
    ! endif
    ! endif
    !
    ! if (id_ice > 0) then
    ! used = send_data (id_ice, ice, time, iis, jjs)
    ! used = send_data (id_ice, ice, time, is_in = iis, js_in = jjs)
    ! if (mp_print .and. seconds == 0) then
    ! tot_prec = g_sum (ice, is, ie, js, je, area, 1)
    ! if (root_proc) write (*, *) 'mean ice_mp = ', tot_prec
    ! endif
    ! endif
    !
    ! if (id_rain > 0) then
    ! used = send_data (id_rain, rain, time, iis, jjs)
    ! used = send_data (id_rain, rain, time, is_in = iis, js_in = jjs)
    ! if (mp_print .and. seconds == 0) then
    ! tot_prec = g_sum (rain, is, ie, js, je, area, 1)
    ! if (root_proc) write (*, *) 'mean rain = ', tot_prec
    ! endif
    ! endif
    !
    ! if (id_rh > 0) then !not used?
    ! used = send_data (id_rh, rh0, time, iis, jjs)
    ! used = send_data (id_rh, rh0, time, is_in = iis, js_in = jjs)
    ! endif
    !
    !
    ! if (id_prec > 0) then
    ! used = send_data (id_prec, prec_mp, time, iis, jjs)
    ! used = send_data (id_prec, prec_mp, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (mp_print) then
    ! prec1 (:, :) = prec1 (:, :) + prec_mp (:, :)
    ! if (seconds == 0) then
    ! prec1 (:, :) = prec1 (:, :) * dt_in / 86400.
    ! tot_prec = g_sum (prec1, is, ie, js, je, area, 1)
    ! if (root_proc) write (*, *) 'daily prec_mp = ', tot_prec
    ! prec1 (:, :) = 0.
    ! endif
    ! endif

    ! call mpp_clock_end (gfdl_mp_clock)

  end subroutine gfdl_cloud_microphys_driver

  ! -----------------------------------------------------------------------
  !>@brief gfdl cloud microphysics, major program
  !>@details lin et al., 1983, jam, 1065 - 1092, and
  !! rutledge and hobbs, 1984, jas, 2949 - 2972
  !! terminal fall is handled lagrangianly by conservative fv algorithm
  !>@param pt: temperature (k)
  !>@param 6 water species:
  !>@param 1) qv: water vapor (kg / kg)
  !>@param 2) ql: cloud water (kg / kg)
  !>@param 3) qr: rain (kg / kg)
  !>@param 4) qi: cloud ice (kg / kg)
  !>@param 5) qs: snow (kg / kg)
  !>@param 6) qg: graupel (kg / kg)
  ! -----------------------------------------------------------------------
  subroutine mpdrv (hydrostatic, uin, vin, w, delp, pt, qv, ql, qr, qi, qs, &
       qg, qa, qn, dz, is, ie, js, je, ks, ke, ktop, kbot, dt_in, ntimes, &
       rain, snow, graupel, ice, m2_rain, m2_sol, cond, area1, land, &
       cnv_fraction, srf_type, eis, rhcrit, anv_icefall, lsc_icefall, revap, isubl, &
       u_dt, v_dt, pt_dt, qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, &
       w_var, vt_r, vt_s, vt_g, vt_i, qn2)

    implicit none

    logical, intent (in) :: hydrostatic

    integer, intent (in) :: is, ie, js, je, ks, ke
    integer, intent (in) :: ntimes, ktop, kbot

    real, intent (in) :: dt_in

    real, intent (in), dimension (is:, js:) :: area1, land
    real, intent (in), dimension (is:, js:) :: cnv_fraction
    real, intent (in), dimension (is:, js:) :: srf_type
    real, intent (in), dimension (is:, js:) :: eis

    real, intent (in), dimension (is:, js:, ks:) :: rhcrit

    real, intent (in) :: anv_icefall, lsc_icefall

    real, intent (in), dimension (is:, js:, ks:) :: uin, vin, delp, pt, dz
    real, intent (in), dimension (is:, js:, ks:) :: qv, qi, ql, qr, qs, qg, qa, qn

    real, intent (inout), dimension (is:, js:, ks:) :: u_dt, v_dt, w, pt_dt, qa_dt
    real, intent (inout), dimension (is:, js:, ks:) :: qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt
    real, intent (  out), dimension (is:, js:, ks:) :: revap, isubl

    real, intent (inout), dimension (is:, js:) :: rain, snow, ice, graupel, cond

    real, intent (out), dimension (is:, js:) :: w_var

    real, intent (out), dimension (is:, js:, ks:) :: vt_r, vt_s, vt_g, vt_i, qn2

    real, intent (out), dimension (is:, js:, ks:) :: m2_rain, m2_sol

    real, dimension (is:ie, js:je, ktop:kbot) :: h_var1d
    real, dimension (is:ie, js:je, ktop:kbot) :: qvz, qlz, qrz, qiz, qsz, qgz, qaz
    real, dimension (is:ie, js:je, ktop:kbot) :: vtiz, vtsz, vtgz, vtrz
    real, dimension (is:ie, js:je, ktop:kbot) :: dp1, dz1
    real, dimension (is:ie, js:je, ktop:kbot) :: qv0, ql0, qr0, qi0, qs0, qg0
    real, dimension (is:ie, js:je, ktop:kbot) :: den, tz, p1, denfac
    real, dimension (is:ie, js:je, ktop:kbot) :: ccn, c_praut, m1_rain, m1_sol, m1, evap1, subl1
    real, dimension (is:ie, js:je, ktop:kbot) :: w1

    real :: cpaut, t0, dts, rdt, den0
    real, dimension(is:ie, js:je) :: r1, s1, i1, g1
    real :: cvm, omq
    real :: u1_k, u1_km1, v1_k, v1_km1

    integer :: i, j, k, n

    integer :: num_devices, nteams, nthreads
    logical :: initial_device

    dts = dt_in / real (ntimes)
    rdt = 1. / dt_in

    ! -----------------------------------------------------------------------
    ! calculate cloud condensation nuclei (ccn)
    ! the following is based on klein eq. 15
    ! -----------------------------------------------------------------------

    !$omp target data &
    ! IN
    !$omp   map(to: &
    !$omp     area1, land, cnv_fraction, srf_type, eis, &
    !$omp     rhcrit, anv_icefall, lsc_icefall, &
    !$omp     uin, vin, delp, pt, dz, &
    !$omp     qv, qi, ql, qr, qs, qg, qa, qn) &

    ! IN/OUT
    !$omp   map(tofrom: &
    !$omp     u_dt, v_dt, w, pt_dt, qa_dt, &
    !$omp     qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, &
    !$omp     rain, snow, ice, graupel, cond) &

    ! OUT
    !$omp   map(from: &
    !$omp     revap, isubl, w_var, &
    !$omp     vt_r, vt_s, vt_g, vt_i, qn2, m2_rain, m2_sol) &

    ! LOCAL
    !$omp   map(alloc: &
    !$omp     h_var1d, &
    !$omp     qvz, qlz, qiz, qrz, qsz, qgz, &
    !$omp     vtiz, vtsz, vtgz, vtrz, &
    !$omp     dp1, dz1, &
    !$omp     qv0, ql0, qr0, qi0, qs0, qg0, &
    !$omp     den, tz, p1, denfac, &
    !$omp     ccn, c_praut, m1_rain, m1_sol, m1, evap1, subl1, w1, &
    !$omp     r1, i1, s1, g1)

    !$omp target teams distribute parallel do simd collapse(3) private(t0)
    do k = ktop, kbot
       do j = js, je
          do i = is, ie
             ! Initialize
             m2_rain (i, j, k) = 0.
             m2_sol (i, j, k) = 0.
             revap (i, j, k) = 0.
             isubl (i, j, k) = 0.

             t0 = pt (i, j, k)
             tz (i, j, k) = t0
             dp1 (i, j, k) = delp (i, j, k)
             ! dp0 (k) = dp1 (k) ! moist air mass * grav

             ! import horizontal subgrid variability with pressure dependence
             ! total water subgrid deviation in horizontal direction
             ! default area dependent form: use dx ~ 100 km as the base
             h_var1d (i, j, k) = min(0.30,1.0 - rhcrit(i,j,k)) ! restricted to 70%
#ifdef __GFORTRAN_TEST__
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! Convert moist mixing ratios to dry mixing ratios

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do j = js, je
          do i = is, ie
#endif
             qvz (i, j, k) = qv (i, j, k)
             qlz (i, j, k) = ql (i, j, k)
             qiz (i, j, k) = qi (i, j, k)
             qrz (i, j, k) = qr (i, j, k)
             qsz (i, j, k) = qs (i, j, k)
             qgz (i, j, k) = qg (i, j, k)
#ifdef __GFORTRAN_TEST__
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do j = js, je
          do i = is, ie
#endif
             qaz (i, j, k) = qa (i, j, k)
#ifdef __GFORTRAN_TEST__
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    !$omp target teams distribute parallel do simd collapse(3) default(shared) private(omq, den0)
    do k = ktop, kbot
       do j = js, je
          do i = is, ie
#endif
             ! dp1: dry air_mass
             ! dp1 (k) = dp1 (k) * (1. - (qvz (k) + qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)))
             dp1 (i, j, k) = dp1 (i, j, k) * (1. - qvz (i, j, k)) ! gfs
             omq = delp (i, j, k) / dp1 (i, j, k)
             qvz (i, j, k) = qvz (i, j, k) * omq
             qlz (i, j, k) = qlz (i, j, k) * omq
             qrz (i, j, k) = qrz (i, j, k) * omq
             qiz (i, j, k) = qiz (i, j, k) * omq
             qsz (i, j, k) = qsz (i, j, k) * omq
             qgz (i, j, k) = qgz (i, j, k) * omq

             den0 = - dp1 (i, j, k) / (grav * dz (i, j, k)) ! density of dry air
             p1 (i, j, k) = den0 * rdgas * t0 ! dry air pressure
#ifdef __GFORTRAN_TEST__
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    !$omp target teams distribute parallel do simd collapse(3) default(shared) private(omq, den0)
    do k = ktop, kbot
       do j = js, je
          do i = is, ie
#endif
             ! -----------------------------------------------------------------------
             ! save a copy of old value for computing tendencies
             ! -----------------------------------------------------------------------

             qv0 (i, j, k) = qvz (i, j, k)
             ql0 (i, j, k) = qlz (i, j, k)
             qr0 (i, j, k) = qrz (i, j, k)
             qi0 (i, j, k) = qiz (i, j, k)
             qs0 (i, j, k) = qsz (i, j, k)
             qg0 (i, j, k) = qgz (i, j, k)
#ifdef __GFORTRAN_TEST__
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! for sedi_momentum

    !$omp target teams distribute parallel do simd collapse(3) default(shared) private(cpaut)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
#endif
             m1 (i, j, k) = 0.
             if (do_sedi_w) w1 (i, j, k) = w (i, j, k)
#ifdef __GFORTRAN_TEST__
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    !$omp target teams distribute parallel do simd collapse(3) default(shared) private(cpaut)
    do k = ktop, kbot
       do j = js, je
          do i = is, ie
#endif
             cpaut = c_paut * 0.104 * grav / 1.717e-5
             ! ccn needs units #/m^3
             if (prog_ccn) then
                ! qn has units # / m^3
                ccn (i, j, k) = qn (i, j, k)
                c_praut (i, j, k) = cpaut * (ccn (i, j, k) * rhor) ** (- 1. / 3.)
             else
                ! qn has units # / m^3
                ccn (i, j, k) = qn (i, j, k)
                !!! use GEOS ccn: ccn (i, j, k) = (ccn_l * land (i) + ccn_o * (1. - land (i))) * 1.e6
                c_praut (i, j, k) = cpaut * (ccn (i, j, k) * rhor) ** (- 1. / 3.)
             endif
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    !$omp end target data

  end subroutine mpdrv

  ! =======================================================================
  !>@brief The subroutine 'setup'm' sets up
  !! gfdl cloud microphysics parameters.
  ! =======================================================================

  subroutine setupm

    implicit none

    real :: gcon, cd, scm3, pisq, act (8)
    real :: vdifu, tcond
    real :: visk
    real :: ch2o, hltf
    real :: hlts, hltc, ri50

    real, parameter :: gam263 = 1.456943, gam275 = 1.608355, gam290 = 1.827363, &
         gam325 = 2.54925, gam350 = 3.323363, gam380 = 4.694155, &
         gam425 = 8.285063, gam450 = 11.631769, gam480 = 17.837789, &
         gam625 = 184.860962, gam680 = 496.604067

    ! intercept parameters

    real, parameter :: rnzr = 8.0e6 ! lin83
    real, parameter :: rnzs = 3.0e6 ! lin83
    real, parameter :: rnzg = 4.0e6 ! rh84

    ! density parameters

    real, parameter :: rhos = 0.1e3 !< lin83 (snow density; 1 / 10 of water)
    real, parameter :: rhog = 0.4e3 !< rh84 (graupel density)
    real, parameter :: acc (3) = (/ 5.0, 2.0, 0.5 /)

    integer :: i, k

    pie = 4. * atan (1.0)

    vdifu = 2.11e-5
    tcond = 2.36e-2

    visk = 1.259e-5
    hlts = 2.8336e6
    hltc = 2.5e6
    hltf = 3.336e5

    ch2o = 4.1855e3
    ri50 = 1.e-4

    pisq = pie * pie
    scm3 = (visk / vdifu) ** (1. / 3.)

    cracs = pisq * rnzr * rnzs * rhos
    csacr = pisq * rnzr * rnzs * rhor
    cgacr = pisq * rnzr * rnzg * rhor
    cgacs = pisq * rnzg * rnzs * rhos
    cgacs = cgacs * c_pgacs

    ! act: 1 - 2:racs (s - r) ; 3 - 4:sacr (r - s) ;
    ! 5 - 6:gacr (r - g) ; 7 - 8:gacs (s - g)

    act (1) = pie * rnzs * rhos
    act (2) = pie * rnzr * rhor
    act (6) = pie * rnzg * rhog
    act (3) = act (2)
    act (4) = act (1)
    act (5) = act (2)
    act (7) = act (1)
    act (8) = act (6)

    do i = 1, 3
       do k = 1, 4
          acco (i, k) = acc (i) / (act (2 * k - 1) ** ((7 - i) * 0.25) * act (2 * k) ** (i * 0.25))
       enddo
    enddo

    gcon = 40.74 * sqrt (sfcrho) ! 44.628

    csacw = pie * rnzs * clin * gam325 / (4. * act (1) ** 0.8125)
    ! decreasing csacw to reduce cloud water --- > snow

    craci = pie * rnzr * alin * gam380 / (4. * act (2) ** 0.95)
    csaci = csacw * c_psaci

    cgacw = pie * rnzg * gam350 * gcon / (4. * act (6) ** 0.875)

    cgaci = cgacw * c_pgaci

    cracw = craci ! cracw = 3.27206196043822
    cracw = c_cracw * cracw

    ! subl and revp: five constants for three separate processes

    cssub (1) = 2. * pie * vdifu * tcond * rvgas * rnzs
    cgsub (1) = 2. * pie * vdifu * tcond * rvgas * rnzg
    crevp (1) = 2. * pie * vdifu * tcond * rvgas * rnzr
    cssub (2) = 0.78 / sqrt (act (1))
    cgsub (2) = 0.78 / sqrt (act (6))
    crevp (2) = 0.78 / sqrt (act (2))
    cssub (3) = 0.31 * scm3 * gam263 * sqrt (clin / visk) / act (1) ** 0.65625
    cgsub (3) = 0.31 * scm3 * gam275 * sqrt (gcon / visk) / act (6) ** 0.6875
    crevp (3) = 0.31 * scm3 * gam290 * sqrt (alin / visk) / act (2) ** 0.725
    cssub (4) = tcond * rvgas
    cssub (5) = hlts ** 2 * vdifu
    cgsub (4) = cssub (4)
    crevp (4) = cssub (4)
    cgsub (5) = cssub (5)
    crevp (5) = hltc ** 2 * vdifu

    cgfr (1) = 20.e2 * pisq * rnzr * rhor / act (2) ** 1.75
    cgfr (2) = 0.66

    ! smlt: five constants (lin et al. 1983)

    csmlt (1) = 2. * pie * tcond * rnzs / hltf
    csmlt (2) = 2. * pie * vdifu * rnzs * hltc / hltf
    csmlt (3) = cssub (2)
    csmlt (4) = cssub (3)
    csmlt (5) = ch2o / hltf

    ! gmlt: five constants

    cgmlt (1) = 2. * pie * tcond * rnzg / hltf
    cgmlt (2) = 2. * pie * vdifu * rnzg * hltc / hltf
    cgmlt (3) = cgsub (2)
    cgmlt (4) = cgsub (3)
    cgmlt (5) = ch2o / hltf

    es0 = 6.107799961e2 ! ~6.1 mb
    ces0 = eps * es0

    !$omp target update to( &
    !$omp     ces0, cracs, cracw, &
    !$omp     csaci, csacr, csacw, cgaci, cgacr, cgacs, cgacw, &
    !$omp     cssub(:), crevp(:), csmlt(:), cgmlt(:), cgfr(:), acco(:,:))

  end subroutine setupm

  ! =======================================================================
  ! initialization of gfdl cloud microphysics
  !>@brief The subroutine 'gfdl_cloud_microphys_init' initializes the GFDL
  !! cloud microphysics.
  ! =======================================================================

  subroutine gfdl_cloud_microphys_init ()

    implicit none

    integer :: file_handle, rc
    character (len = 64) :: file_name = 'input-data/input.nml'
    logical :: exists

    inquire (file = trim (file_name), exist = exists)
    if (.not. exists) then
       write (6, *) 'gfdl - mp :: namelist file: ', trim (file_name), ' does not exist'
       stop
    else
       open(newunit = file_handle, file = file_name, status = 'old')
       read(nml = gfdl_cloud_microphysics_nml, unit = file_handle, iostat=rc)
       if (rc /=0) error stop "Could not read input namelist file"
       close(file_handle)
    endif

    !$omp target update to( &
    !$omp     tau_revp, tau_v2l, tau_l2v, tau_i2v, tau_s2v, tau_v2s, tau_g2v, &
    !$omp     tau_v2g, tau_frz, tau_imlt, tau_smlt, tau_i2s, tau_g2r, &
    !$omp     tice, tice0, rh_inc, rh_inr, t_min, do_qa, t_sub, do_evap, &
    !$omp     do_bigg, qi_lim, do_subl, preciprad, icloud_f, qc_crt, z_slope_ice, &
    !$omp     c_paut, prog_ccn, fix_negative, sedi_transport, ql_mlt, qs_mlt, qi0_crt, qs0_crt, &
    !$omp     const_vi, vi_fac, vi_max, const_vs, vs_fac, vs_max, const_vg, vg_fac, vg_max, const_vr, vr_fac, vr_max, &
    !$omp     use_ppm, mono_prof, rthreshs, rthreshu, irain_f, z_slope_liq, do_sedi_heat, &
    !$omp     ql0_max, dt_fr, sat_adj0, dw_land, dw_ocean, c_psaci, c_pgacs, &
    !$omp     ccn_l, ccn_o, c_cracw, use_ccn, de_ice, mp_time)

    if (do_setup) then
       call setup_con
       call setupm
       do_setup = .false.
    endif

    module_is_initialized = .true.

  end subroutine gfdl_cloud_microphys_init

  ! =======================================================================
  ! qsmith table initialization
  !>@brief The subroutine 'setup_con' sets up constants and calls 'qsmith_init'.
  ! =======================================================================

  subroutine setup_con

    implicit none

    ! root_proc = (mpp_pe () .eq.mpp_root_pe ())

    rgrav = 1. / grav

    if (.not. qsmith_tables_initialized) call qsmith_init

    qsmith_tables_initialized = .true.

  end subroutine setup_con

  ! =======================================================================
  ! initialization
  ! prepare saturation water vapor pressure tables
  ! =======================================================================
  !>@brief The subroutine 'qsmith_init' initializes lookup tables for saturation
  !! water vapor pressure for the following utility routines that are designed
  !! to return qs consistent with the assumptions in FV3.
  !>@details The calculations are highly accurate values based on the Clausius-Clapeyron
  !! equation.
  ! =======================================================================
  subroutine qsmith_init

    implicit none

    integer, parameter :: length = TABLE_LENGTH

    integer :: i

    if (.not. tables_are_initialized) then

       ! root_proc = (mpp_pe () .eq. mpp_root_pe ())
       ! if (root_proc) print *, ' gfdl mp: initializing qs tables'

       ! debug code
       ! print *, mpp_pe (), allocated (table), allocated (table2), &
       ! allocated (table3), allocated (tablew), allocated (des), &
       ! allocated (des2), allocated (des3), allocated (desw)
       ! end debug code

       ! generate es table (dt = 0.1 deg. c)

       call qs_table (length)
       call qs_table2 (length)
       call qs_table3 (length)
       call qs_tablew (length)

       do i = 1, length - 1
          des (i) = max (0., table (i + 1) - table (i))
          des2 (i) = max (0., table2 (i + 1) - table2 (i))
          des3 (i) = max (0., table3 (i + 1) - table3 (i))
          desw (i) = max (0., tablew (i + 1) - tablew (i))
       enddo
       des (length) = des (length - 1)
       des2 (length) = des2 (length - 1)
       des3 (length) = des3 (length - 1)
       desw (length) = desw (length - 1)

       tables_are_initialized = .true.

    !$omp target update to(table2(:), des2(:), tablew(:), desw(:))

    endif

  end subroutine qsmith_init

  ! =======================================================================
  !>@brief saturation water vapor pressure table ii
  ! 1 - phase table
  ! =======================================================================

  subroutine qs_tablew (n)

    implicit none

    integer, intent (in) :: n

    real :: delt = 0.1
    real :: tmin, tem, fac0, fac1, fac2

    integer :: i

    tmin = table_ice - 160.

    ! -----------------------------------------------------------------------
    ! compute es over water
    ! -----------------------------------------------------------------------

    do i = 1, n
       tem = tmin + delt * real (i - 1)
       fac0 = (tem - t_ice) / (tem * t_ice)
       fac1 = fac0 * lv0
       fac2 = (dc_vap * log (tem / t_ice) + fac1) / rvgas
       tablew (i) = e00 * exp (fac2)
    enddo

  end subroutine qs_tablew

  ! =======================================================================
  !>@brief saturation water vapor pressure table iii
  ! 2 - phase table
  ! =======================================================================

  subroutine qs_table2 (n)

    implicit none

    integer, intent (in) :: n

    real :: delt = 0.1
    real :: tmin, tem0, tem1, fac0, fac1, fac2

    integer :: i, i0, i1

    tmin = table_ice - 160.

    do i = 1, n
       tem0 = tmin + delt * real (i - 1)
       fac0 = (tem0 - t_ice) / (tem0 * t_ice)
       if (i <= 1600) then
          ! -----------------------------------------------------------------------
          ! compute es over ice between - 160 deg c and 0 deg c.
          ! -----------------------------------------------------------------------
          fac1 = fac0 * li2
          fac2 = (d2ice * log (tem0 / t_ice) + fac1) / rvgas
       else
          ! -----------------------------------------------------------------------
          ! compute es over water between 0 deg c and 102 deg c.
          ! -----------------------------------------------------------------------
          fac1 = fac0 * lv0
          fac2 = (dc_vap * log (tem0 / t_ice) + fac1) / rvgas
       endif
       table2 (i) = e00 * exp (fac2)
    enddo

    ! -----------------------------------------------------------------------
    ! smoother around 0 deg c
    ! -----------------------------------------------------------------------

    i0 = 1600
    i1 = 1601
    tem0 = 0.25 * (table2 (i0 - 1) + 2. * table (i0) + table2 (i0 + 1))
    tem1 = 0.25 * (table2 (i1 - 1) + 2. * table (i1) + table2 (i1 + 1))
    table2 (i0) = tem0
    table2 (i1) = tem1

  end subroutine qs_table2

  ! =======================================================================
  !>@brief saturation water vapor pressure table iv
  ! 2 - phase table with " - 2 c" as the transition point
  ! =======================================================================

  subroutine qs_table3 (n)

    implicit none

    integer, intent (in) :: n

    real :: delt = 0.1
    real :: esbasw, tbasw, esbasi, tmin, tem, aa, b, c, d, e
    real :: tem0, tem1

    integer :: i, i0, i1

    esbasw = 1013246.0
    tbasw = table_ice + 100.
    esbasi = 6107.1
    tmin = table_ice - 160.

    do i = 1, n
       tem = tmin + delt * real (i - 1)
       ! if (i <= 1600) then
       if (i <= 1580) then ! change to - 2 c
          ! -----------------------------------------------------------------------
          ! compute es over ice between - 160 deg c and 0 deg c.
          ! see smithsonian meteorological tables page 350.
          ! -----------------------------------------------------------------------
          aa = - 9.09718 * (table_ice / tem - 1.)
          b = - 3.56654 * alog10 (table_ice / tem)
          c = 0.876793 * (1. - tem / table_ice)
          e = alog10 (esbasi)
          table3 (i) = 0.1 * 10 ** (aa + b + c + e)
       else
          ! -----------------------------------------------------------------------
          ! compute es over water between - 2 deg c and 102 deg c.
          ! see smithsonian meteorological tables page 350.
          ! -----------------------------------------------------------------------
          aa = - 7.90298 * (tbasw / tem - 1.)
          b = 5.02808 * alog10 (tbasw / tem)
          c = - 1.3816e-7 * (10 ** ((1. - tem / tbasw) * 11.344) - 1.)
          d = 8.1328e-3 * (10 ** ((tbasw / tem - 1.) * (- 3.49149)) - 1.)
          e = alog10 (esbasw)
          table3 (i) = 0.1 * 10 ** (aa + b + c + d + e)
       endif
    enddo

    ! -----------------------------------------------------------------------
    ! smoother around - 2 deg c
    ! -----------------------------------------------------------------------

    i0 = 1580
    i1 = 1581
    tem0 = 0.25 * (table3 (i0 - 1) + 2. * table (i0) + table3 (i0 + 1))
    tem1 = 0.25 * (table3 (i1 - 1) + 2. * table (i1) + table3 (i1 + 1))
    table3 (i0) = tem0
    table3 (i1) = tem1

  end subroutine qs_table3

  ! =======================================================================
  !>@brief saturation water vapor pressure table i
  ! 3 - phase table
  ! =======================================================================

  subroutine qs_table (n)

    implicit none

    integer, intent (in) :: n

    real :: delt = 0.1
    real :: tmin, tem, esh40
    real :: wice, wh2o, fac0, fac1, fac2
    real :: esupc (400)

    integer :: i
    real :: tc

    tmin = table_ice - 160.

    ! -----------------------------------------------------------------------
    ! compute es over ice between - 160 deg c and 0 deg c.
    ! -----------------------------------------------------------------------

    do i = 1, 1600
       tem = tmin + delt * real (i - 1)
       fac0 = (tem - t_ice) / (tem * t_ice)
       fac1 = fac0 * li2
       fac2 = (d2ice * log (tem / t_ice) + fac1) / rvgas
       table (i) = e00 * exp (fac2)
    enddo

    ! -----------------------------------------------------------------------
    ! compute es over water between - 40 deg c and 102 deg c.
    ! -----------------------------------------------------------------------

    do i = 1, 1421
       tem = 233.16 + delt * real (i - 1)
       fac0 = (tem - t_ice) / (tem * t_ice)
       fac1 = fac0 * lv0
       fac2 = (dc_vap * log (tem / t_ice) + fac1) / rvgas
       esh40 = e00 * exp (fac2)
       if (i <= 400) then
          esupc (i) = esh40
       else
          table (i + 1200) = esh40
       endif
    enddo

    ! -----------------------------------------------------------------------
    ! derive blended es over ice and supercooled water between - 40 deg c and 0 deg c
    ! -----------------------------------------------------------------------

    do i = 1, 400
       tem = 233.16 + delt * real (i - 1)
       ! wice = 0.05 * (table_ice - tem)
       ! wh2o = 0.05 * (tem - 253.16)
       ! GEOS ! WMP impose CALIPSO ice polynomial from 0 C to -40 C
       wice = ice_fraction(tem,0.0,0.0)
       wh2o = 1.0 - wice
       table (i + 1200) = wice * table (i + 1200) + wh2o * esupc (i)
    enddo

  end subroutine qs_table

  function ICE_FRACTION (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
    !$acc routine seq
    !$omp declare target
    real, intent(in) :: TEMP,CNV_FRACTION,SRF_TYPE
    real             :: ICEFRCT
    real             :: tc, ptc
    real             :: ICEFRCT_C, ICEFRCT_M

    ! In anvil/convective clouds
    real, parameter :: aT_ICE_ALL = 252.16
    real, parameter :: aT_ICE_MAX = 268.16
    real, parameter :: aICEFRPWR  = 2.0
    ! Over snow/ice SRF_TYPE = 2
    real, parameter :: iT_ICE_ALL = 236.16
    real, parameter :: iT_ICE_MAX = 261.16
    real, parameter :: iICEFRPWR  = 6.0
    ! Over Land     SRF_TYPE = 1
    real, parameter :: lT_ICE_ALL = 239.16
    real, parameter :: lT_ICE_MAX = 261.16
    real, parameter :: lICEFRPWR  = 2.0
    ! Over Oceans   SRF_TYPE = 0
    real, parameter :: oT_ICE_ALL = 238.16
    real, parameter :: oT_ICE_MAX = 263.16
    real, parameter :: oICEFRPWR  = 4.0

    ! Anvil clouds
    ! Anvil-Convective sigmoidal function like figure 6(right)
    ! Sigmoidal functions Hu et al 2010, doi:10.1029/2009JD012384
    ICEFRCT_C  = 0.00
    if ( TEMP <= aT_ICE_ALL ) then
       ICEFRCT_C = 1.000
    else if ( (TEMP > aT_ICE_ALL) .AND. (TEMP <= aT_ICE_MAX) ) then
       ICEFRCT_C = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - aT_ICE_ALL ) / ( aT_ICE_MAX - aT_ICE_ALL ) ) )
    end if
    ICEFRCT_C = MIN(ICEFRCT_C,1.00)
    ICEFRCT_C = MAX(ICEFRCT_C,0.00)
    ICEFRCT_C = ICEFRCT_C**aICEFRPWR
#ifdef MODIS_ICE_POLY
    ! Use MODIS polynomial from Hu et al, DOI: (10.1029/2009JD012384)
    tc = MAX(-46.0,MIN(TEMP-MAPL_TICE,46.0)) ! convert to celcius and limit range from -46:46 C
    ptc = 7.6725 + 1.0118*tc + 0.1422*tc**2 + 0.0106*tc**3 + 0.000339*tc**4 + 0.00000395*tc**5
    ICEFRCT_M = 1.0 - (1.0/(1.0 + exp(-1*ptc)))
#else
    ! Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
    if (SRF_TYPE == 2.0) then
       ! Over snow/ice
       ICEFRCT_M  = 0.00
       if ( TEMP <= iT_ICE_ALL ) then
          ICEFRCT_M = 1.000
       else if ( (TEMP > iT_ICE_ALL) .AND. (TEMP <= iT_ICE_MAX) ) then
          ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - iT_ICE_ALL ) / ( iT_ICE_MAX - iT_ICE_ALL ) ) )
       end if
       ICEFRCT_M = MIN(ICEFRCT_M,1.00)
       ICEFRCT_M = MAX(ICEFRCT_M,0.00)
       ICEFRCT_M = ICEFRCT_M**iICEFRPWR
    else if (SRF_TYPE > 1.0) then
       ! Over Land
       ICEFRCT_M  = 0.00
       if ( TEMP <= lT_ICE_ALL ) then
          ICEFRCT_M = 1.000
       else if ( (TEMP > lT_ICE_ALL) .AND. (TEMP <= lT_ICE_MAX) ) then
          ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - lT_ICE_ALL ) / ( lT_ICE_MAX - lT_ICE_ALL ) ) )
       end if
       ICEFRCT_M = MIN(ICEFRCT_M,1.00)
       ICEFRCT_M = MAX(ICEFRCT_M,0.00)
       ICEFRCT_M = ICEFRCT_M**lICEFRPWR
    else
       ! Over Oceans
       ICEFRCT_M  = 0.00
       if ( TEMP <= oT_ICE_ALL ) then
          ICEFRCT_M = 1.000
       else if ( (TEMP > oT_ICE_ALL) .AND. (TEMP <= oT_ICE_MAX) ) then
          ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - oT_ICE_ALL ) / ( oT_ICE_MAX - oT_ICE_ALL ) ) )
       end if
       ICEFRCT_M = MIN(ICEFRCT_M,1.00)
       ICEFRCT_M = MAX(ICEFRCT_M,0.00)
       ICEFRCT_M = ICEFRCT_M**oICEFRPWR
    endif
#endif
    ! Combine the Convective and MODIS functions
    ICEFRCT  = ICEFRCT_M*(1.0-CNV_FRACTION) + ICEFRCT_C*(CNV_FRACTION)

  end function ICE_FRACTION

end module gfdl2_cloud_microphys_mod
