module sparse_mod
    !> Facade module providing unified interface to modular sparse matrix functionality
    !>
    !> This module acts as a unified interface to the refactored sparse matrix modules,
    !> maintaining backward compatibility while leveraging the improved modular design.
    !>
    !> The functionality is distributed across specialized modules:
    !> - sparse_types_mod: Core type definitions
    !> - sparse_conversion_mod: Format conversion utilities  
    !> - sparse_io_mod: I/O operations
    !> - sparse_arithmetic_mod: Basic arithmetic operations (includes sparse_talk)
    !> - sparse_solvers_mod: Advanced solver algorithms (includes sparse_solve_method)
    !> - sparse_utils_mod: Utility functions
    
    ! Re-export all functionality from the modular sparse modules
    use sparse_types_mod
    use sparse_conversion_mod
    use sparse_io_mod
    use sparse_arithmetic_mod
    use sparse_solvers_mod
    use sparse_utils_mod
    
    implicit none
    
    ! Re-export all public entities from the modular modules
    public
    
end module sparse_mod