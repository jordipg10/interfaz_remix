!subroutine compute_rech_mat_PDE(this)
!!> $E_ii=r_i$
!    use PDE_m, only: PDE_1D_c
!    use transport_m, only: transport_1D_c, diffusion_1D_c
!    use transport_transient_m, only: transport_1D_transient_c, diffusion_1D_transient_c
!    implicit none
!    class(PDE_1D_c) :: this
!    
!    !allocate(this%source_term_PDE(this%spatial_discr%Num_targets))
!    !this%source_term_PDE=0d0 !> $g=0$ chapuza
!     select type (this)
!     type is (transport_1D_c)
!         this%rech_mat%diag=this%tpt_props_heterog%source_term
!     type is (transport_1D_transient_c)
!         this%rech_mat%diag=this%tpt_props_heterog%source_term
!     type is (diffusion_1D_c)
!         this%rech_mat%diag=this%diff_props_heterog%source_term
!     type is (diffusion_1D_transient_c)
!         this%rech_mat%diag=this%diff%diff_props_heterog%source_term
!     end select
!
!end subroutine