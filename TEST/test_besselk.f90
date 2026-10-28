program test_besselk
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use gsl_specialfunctions_mod, only: besselk

  implicit none

  real(dp), parameter :: rtol = 1.0e-12_dp

  call assert_close(besselk(0, 0.5_dp), 9.2441907122766565e-01_dp, "K_0(0.5)")
  call assert_close(besselk(0, 1.0_dp), 4.2102443824070834e-01_dp, "K_0(1.0)")
  call assert_close(besselk(0, 2.0_dp), 1.1389387274953341e-01_dp, "K_0(2.0)")
  call assert_close(besselk(0, 5.0_dp), 3.6910983340425942e-03_dp, "K_0(5.0)")
  call assert_close(besselk(1, 0.5_dp), 1.6564411200033007e+00_dp, "K_1(0.5)")
  call assert_close(besselk(1, 1.0_dp), 6.0190723019723458e-01_dp, "K_1(1.0)")
  call assert_close(besselk(1, 2.0_dp), 1.3986588181652246e-01_dp, "K_1(2.0)")
  call assert_close(besselk(1, 5.0_dp), 4.0446134454521637e-03_dp, "K_1(5.0)")
  call assert_close(besselk(2, 0.5_dp), 7.5501835512408686e+00_dp, "K_2(0.5)")
  call assert_close(besselk(2, 1.0_dp), 1.6248388986351774e+00_dp, "K_2(1.0)")
  call assert_close(besselk(2, 2.0_dp), 2.5375975456605587e-01_dp, "K_2(2.0)")
  call assert_close(besselk(2, 5.0_dp), 5.3089437122234599e-03_dp, "K_2(5.0)")

  print *, "All tests passed!"

contains

  subroutine assert_close(actual, expected, label)
    real(dp), intent(in) :: actual, expected
    character(*), intent(in) :: label
    real(dp) :: rel_err

    rel_err = abs(actual - expected)/abs(expected)
    if (rel_err > rtol) then
      print *, "FAIL: ", label, " expected ", expected, " got ", actual, &
           " rel_err ", rel_err
      error stop 1
    end if
  end subroutine assert_close

end program test_besselk
