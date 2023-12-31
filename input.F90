module input_mod

  implicit none

  private

  public InputScalars_T, InputArrays_T, get_data_from_file, write_difference

  character(len=*), parameter :: fmt_diff = '(1x, a10, 1x, a1, 1x, e15.9)'

  type InputScalars_T
     logical :: hydrostatic, phys_hydrostatic
     integer :: iis, iie, jjs, jje !< physics window
     integer :: kks, kke !< vertical dimension
     integer :: ktop, kbot !< vertical compute domain
     real :: dt_in !< physics time step
     real :: anv_icefall, lsc_icefall
   contains
     procedure :: write_scalars
  end type InputScalars_T

  type InputArrays_T
     real, allocatable, dimension (:, :) :: area, land, cnv_fraction, srf_type, eis
     real, allocatable, dimension (:, :, :) :: rhcrit
     real, allocatable, dimension (:, :, :) :: delp, dz, uin, vin, pt
     real, allocatable, dimension (:, :, :) :: qv, ql, qr, qg, qa, qn
     real, allocatable, dimension (:, :, :) :: qi, qs
     real, allocatable, dimension (:, :, :) :: pt_dt, qa_dt, udt, vdt, w
     real, allocatable, dimension (:, :, :) :: qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt
     ! real, allocatable, dimension (:, :) :: rain, snow, ice, graupel
     ! real, allocatable, dimension (:, :, :) :: m2_rain, m2_sol ! Rain and Ice fluxes (Pa kg/kg)
     ! real, allocatable, dimension (:, :, :) :: revap ! Rain evaporation
     ! real, allocatable, dimension (:, :, :) :: isubl ! Ice sublimation
   contains
     procedure :: write_arrays
  end type InputArrays_T

contains

  subroutine write_scalars(self)

    ! Arguments
    class(InputScalars_T), intent(in) :: self

    ! Start
    print *, 'iis/iie/jjs/jje: ', self%iis, self%iie, self%jjs, self%jje
    print *, '/kks/kke/ktop/kbot: ', self%kks, self%kke, self%ktop, self%kbot
    print *, 'dt_in: ', self%dt_in
    print *, 'anv/lsc_icefall: ', self%anv_icefall, self%lsc_icefall
    print *, 'hydrostatic/phys_hydrostatic: ', self%hydrostatic, self%phys_hydrostatic

  end subroutine write_scalars
  
  subroutine write_arrays(self)

    ! Arguments
    class(InputArrays_T), intent(in) :: self

    ! Start
    print *, 'area (min/max/sum): ', minval(self%area), maxval(self%area), sum(self%area)
    print *, 'delp (min/max/sum): ', minval(self%delp), maxval(self%delp), sum(self%delp)
    print *, 'qv (min/max/sum): ', minval(self%qv), maxval(self%qv), sum(self%qv)
    print *, 'qi (min/max/sum): ', minval(self%qi), maxval(self%qi), sum(self%qi)
    print *, 'pt_dt (min/max/sum): ', minval(self%pt_dt), maxval(self%pt_dt), sum(self%pt_dt)
    print *, 'qv_dt (min/max/sum): ', minval(self%qv_dt), maxval(self%qv_dt), sum(self%qv_dt)

  end subroutine write_arrays

  subroutine get_data_from_file(file_name, sclr, arr)

    ! Arguments
    character(len=*), intent(in) :: file_name
    type(InputScalars_T), intent(out) :: sclr
    type(InputArrays_T), intent(out) :: arr

    ! Locals
    integer :: file_handle
    integer :: iis, iie, jjs, jje, kks, kke

    ! Start
    ! -Read-input-data-
    open(newunit = file_handle, file = file_name, form = 'unformatted', status = 'old')
    ! -scalars-
    read(file_handle) &
         iis, iie, jjs, jje, &
         kks, kke, sclr%ktop, sclr%kbot, &
         sclr%dt_in, &
         sclr%anv_icefall, sclr%lsc_icefall, &
         sclr%hydrostatic, sclr%phys_hydrostatic
    sclr%iis = iis
    sclr%iie = iie
    sclr%jjs = jjs
    sclr%jje = jje
    sclr%kks = kks
    sclr%kke = kke
    
    ! -arrays-
    ! --first, allocate memory--
    allocate(arr%area(iis:iie, jjs:jje))
    allocate(arr%land, arr%cnv_fraction, arr%srf_type, arr%eis, mold=arr%area)
    allocate(arr%rhcrit(iis:iie, jjs:jje, kks:kke))
    allocate(&
         arr%delp, arr%dz, arr%uin, arr%vin, arr%pt, &
         arr%qv, arr%ql, arr%qr, arr%qg, arr%qa, arr%qn, &
         arr%qi, arr%qs, &
         arr%pt_dt, arr%qa_dt, arr%udt, arr%vdt, arr%w, &
         arr%qv_dt, arr%ql_dt, arr%qr_dt, arr%qi_dt, arr%qs_dt, arr%qg_dt, &
         mold=arr%rhcrit)
    ! allocate(arr%rain, arr%snow, arr%ice, arr%graupel, mold=arr%area)
    ! allocate(arr%m2_rain, arr%m2_sol, arr%revap, arr%isubl, mold=arr%rhcrit)
    ! --now, read--
    read(file_handle) &
         ! intent(in)
         arr%area, arr%land, arr%cnv_fraction, arr%srf_type, arr%eis, arr%rhcrit, &
         arr%delp, arr%dz, arr%uin, arr%vin, arr%pt, &
         arr%qv, arr%ql, arr%qr, arr%qg, arr%qa, arr%qn, &
         ! intent(inout)
         arr%qi, arr%qs, &
         arr%pt_dt, arr%qa_dt, arr%udt, arr%vdt, arr%w, &
         arr%qv_dt, arr%ql_dt, arr%qr_dt, arr%qi_dt, arr%qs_dt, arr%qg_dt
    close(file_handle)
    
  end subroutine get_data_from_file

  subroutine write_difference(arr1, arr2)

    ! Arguments
    type(InputArrays_T), intent(in) :: arr1, arr2

    ! Start
    print *, ''
    print *, '-----------|-----------------'
    print *, ' inout var |  abs error'
    print *, '-----------|-----------------'
    write(*, fmt_diff) 'qi', '|', norm2(arr2%qi-arr1%qi)
    write(*, fmt_diff) 'qs', '|', norm2(arr2%qs-arr1%qs)
    write(*, fmt_diff) 'qv_dt', '|', norm2(arr2%qv_dt-arr1%qv_dt)
    write(*, fmt_diff) 'ql_dt', '|', norm2(arr2%ql_dt-arr1%ql_dt)
    write(*, fmt_diff) 'qr_dt', '|', norm2(arr2%qr_dt-arr1%qr_dt)
    write(*, fmt_diff) 'qi_dt', '|', norm2(arr2%qi_dt-arr1%qi_dt)
    write(*, fmt_diff) 'qs_dt', '|', norm2(arr2%qs_dt-arr1%qs_dt)
    write(*, fmt_diff) 'qg_dt', '|', norm2(arr2%qg_dt-arr1%qg_dt)
    write(*, fmt_diff) 'qa_dt', '|', norm2(arr2%qa_dt-arr1%qa_dt)
    write(*, fmt_diff) 'w', '|', norm2(arr2%w-arr1%w)
    write(*, fmt_diff) 'udt', '|', norm2(arr2%udt-arr1%udt)
    write(*, fmt_diff) 'vdt', '|', norm2(arr2%vdt-arr1%vdt)
    print *, '-----------|-----------------'

  end subroutine write_difference

end module input_mod
