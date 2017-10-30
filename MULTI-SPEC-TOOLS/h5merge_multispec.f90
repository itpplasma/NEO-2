PROGRAM h5merge_multispec

  !************************************
  ! HDF5
  !************************************
  USE hdf5_tools
  USE hdf5_tools_f2003

  IMPLICIT NONE

  INTEGER :: funit = 10
  INTEGER :: stat
  CHARACTER(256) :: surfname, collectionname, cmd

  INTEGER(HID_T) :: h5id_final, h5id_collection

  collectionname = "final_neo2_multispecies_out.h5"
  
  cmd = "if [ -e " // TRIM(ADJUSTL(collectionname)) // " ]; then rm "
  cmd = TRIM(ADJUSTL(cmd)) // " " // TRIM(ADJUSTL(collectionname)) // "; fi"
  !print *,cmd
  !stop
  CALL execute_command_LINE(cmd)
  
  cmd = "find . -name neo2_multispecies_out.h5  -type f -printf '%h\n' | "
  cmd = TRIM(ADJUSTL(cmd)) // "awk '{ c=$0; gsub(""\.\/"","""",c); print c }' > jobs_list.txt"
  !print *,cmd
  !stop
  CALL execute_command_LINE(cmd)
  
  CALL h5_init()
  

  CALL h5_open_rw(collectionname, h5id_collection)
  
  OPEN(funit, file='jobs_list.txt', status='old', action='read')
  DO
     READ(funit, *, iostat=stat) surfname
     IF (stat /= 0) EXIT

     WRITE (*,*) "Processing ", TRIM(surfname)
     CALL h5_open(TRIM(surfname) // '/neo2_multispecies_out.h5', h5id_final)
     CALL h5_copy(h5id_final, '.' , h5id_collection, '/' // TRIM(surfname))
     CALL h5_close(h5id_final)
  END DO

  CALL h5_close(h5id_collection)

  CLOSE(funit)

  CALL h5_deinit()

END PROGRAM h5merge_multispec
