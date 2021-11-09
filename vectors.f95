MODULE vectors

   IMPLICIT NONE
   integer,  parameter :: dp = kind(0.d0) ! double precision

   ! Declare vector data type:
   TYPE :: vector
      REAL(dp) :: x
      REAL(dp) :: y
      REAL(dp) :: z
   END TYPE

   INTERFACE ASSIGNMENT (=)
      MODULE PROCEDURE array_to_vector
      MODULE PROCEDURE vector_to_array
   END INTERFACE ASSIGNMENT (=)

   INTERFACE OPERATOR (+)
      MODULE PROCEDURE vector_add
   END INTERFACE

   INTERFACE OPERATOR (*)
!      MODULE PROCEDURE vector_times_real
      MODULE PROCEDURE real_times_vector
!      MODULE PROCEDURE vector_times_int
!      MODULE PROCEDURE int_times_vector
      MODULE PROCEDURE cross_product
   END INTERFACE

   INTERFACE OPERATOR (.DOT.)
      MODULE PROCEDURE dot_product
   END INTERFACE OPERATOR (.DOT.)

CONTAINS

   SUBROUTINE array_to_vector(vec_result, array)
      TYPE (vector), INTENT(OUT) :: vec_result
      REAL(dp), DIMENSION(3), INTENT(IN) :: array
      vec_result%x = array(1)
      vec_result%y = array(2)
      vec_result%z = array(3)
   END SUBROUTINE array_to_vector

   SUBROUTINE vector_to_array(array_result, vec_1)
      REAL(dp), DIMENSION(3), INTENT(OUT) :: array_result
      TYPE (vector), INTENT(IN) :: vec_1
      array_result(1) = vec_1%x
      array_result(2) = vec_1%y
      array_result(3) = vec_1%z
   END SUBROUTINE vector_to_array

   FUNCTION vector_add(vec_1, vec_2)
      TYPE (vector) :: vector_add
      TYPE (vector), INTENT(IN) :: vec_1, vec_2
      vector_add%x = vec_1%x + vec_2%x
      vector_add%y = vec_1%y + vec_2%y
      vector_add%z = vec_1%z + vec_2%z
   END FUNCTION vector_add

   FUNCTION real_times_vector(scaler, vec_1)
      REAL(dp), INTENT(IN) :: scaler
      TYPE (vector), INTENT(IN) :: vec_1
      TYPE (vector) :: real_times_vector
      real_times_vector%x = scaler * vec_1%x
      real_times_vector%y = scaler * vec_1%y
      real_times_vector%z = scaler * vec_1%z
   END FUNCTION real_times_vector

   FUNCTION cross_product(vec_1, vec_2)
      TYPE (vector), INTENT(IN) :: vec_1, vec_2
      TYPE (vector) :: cross_product
      cross_product%x = vec_1%y * vec_2%z - vec_1%z * vec_2%y
      cross_product%y = vec_1%z * vec_2%x - vec_1%x * vec_2%z
      cross_product%z = vec_1%x * vec_2%y - vec_1%y * vec_2%x
   END FUNCTION cross_product

   FUNCTION dot_product(vec_1, vec_2)
      REAL(dp) :: dot_product
      TYPE (vector), INTENT(IN) :: vec_1, vec_2
      dot_product = vec_1%x*vec_2%x + vec_1%y*vec_2%y &
                  + vec_1%z*vec_2%z 
   END FUNCTION dot_product

END MODULE vectors
