program MicrophysicsDriver

  implicit none

  call serial_driver(3)

contains

  subroutine serial_driver(irank)

    use input_mod, only: InputScalars_T, InputArrays_T, get_data_from_file
    use input_mod, only: write_inout_difference => write_difference
    use output_mod, only: OutputArrays_T, write_output_difference => write_difference
    use gfdl2_cloud_microphys_mod, only: gfdl_cloud_microphys_init, gfdl_cloud_microphys_driver

    implicit none

    ! Arguments
    integer, intent(in):: irank

    ! Locals
    character(len=*), parameter :: fmt = '(1x, a1, i2, a1, 1x, a, f10.7, 1x, a1)'
    character(len=256) :: file_name
    integer :: file_handle
    real :: start, finish, cpu_time_, gpu_time_
    type(InputScalars_T) :: sclr1, sclr2
    type(InputArrays_T) :: inarr1, inarr2
    type(OutputArrays_T) :: outarr1, outarr2

    ! Input file
    write(file_name, '(a26, i2.2, a4)') 'input-data/microphys_data.', irank, '.bin'

    ! Read data and call new (gpu) version
    call get_data_from_file(file_name, sclr2, inarr2)
    ! print *, ''
    ! print *, 'GPU:'
    ! call sclr2%write_scalars()
    ! call inarr2%write_arrays()
    outarr2 = OutputArrays_T(sclr2%iis, sclr2%iie, sclr2%jjs, sclr2%jje, sclr2%kks, sclr2%kke)

    call gfdl_cloud_microphys_init()

    call cpu_time(start)
    call gfdl_cloud_microphys_driver ( &
         ! intent (in)
         inarr2%qv, inarr2%ql, inarr2%qr, &
         ! intent (inout)
         inarr2%qi, inarr2%qs, &
         ! intent (in)
         inarr2%qg, inarr2%qa, inarr2%qn, &
         ! intent (inout)
         inarr2%qv_dt, inarr2%ql_dt, inarr2%qr_dt, inarr2%qi_dt, &
         inarr2%qs_dt, inarr2%qg_dt, inarr2%qa_dt, inarr2%pt_dt, &
         ! intent (in)
         inarr2%pt, &
         ! intent (inout)
         inarr2%w, &
         ! intent (in)
         inarr2%uin, inarr2%vin, &
         ! intent (inout)
         inarr2%udt, inarr2%vdt, &
         ! intent (in)
         inarr2%dz, inarr2%delp, inarr2%area, sclr2%dt_in, inarr2%land, inarr2%cnv_fraction, &
         inarr2%srf_type, inarr2%eis, inarr2%rhcrit, sclr2%anv_icefall, sclr2%lsc_icefall, &
         ! intent (out)
         outarr2%revap, outarr2%isubl, &
         outarr2%rain, outarr2%snow, outarr2%ice, outarr2%graupel, &
         outarr2%m2_rain, outarr2%m2_sol, &
         ! intent (in)
         sclr2%hydrostatic, sclr2%phys_hydrostatic, &
         sclr2%iis, sclr2%iie, sclr2%jjs, sclr2%jje, sclr2%kks, sclr2%kke, sclr2%ktop, sclr2%kbot)
    call cpu_time(finish)
    gpu_time_ = finish - start
    write(*, fmt) '[', irank, ']', 'Time taken (gpu):', gpu_time_, 's'
    ! call outarr2%write_arrays()

  end subroutine serial_driver

end program MicrophysicsDriver
