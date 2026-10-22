program test_phi1_density
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phi1_mod, only: phi1_density
    implicit none
    real(dp) :: pitch(5)
    complex(dp) :: y(10,2), expected_density, expected_energy, value
    integer :: unit
    character(1024) :: fixture
    namelist /density_fixture/ pitch, y, expected_density, expected_energy
    call get_command_argument(1,fixture)
    open(newunit=unit,file=trim(fixture),status='old')
    read(unit,nml=density_fixture)
    close(unit)
    value = phi1_density(y,pitch,[.25_dp,.375_dp],.005_dp)
    if (abs(value-expected_density) > 2e-12_dp) error stop 'D8 density'
    value = phi1_density(y,pitch,[.375_dp,.9375_dp],.005_dp)
    if (abs(value-expected_energy) > 2e-12_dp) error stop 'D8 energy'
end program
