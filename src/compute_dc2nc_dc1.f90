!> Computes Jacobian of secondary variable activity species concentrations
!> \f$\partial \mathbf{c}_{2,v} / \partial \mathbf{c}_1\f$ with respect to the primary species concentrations for
!> a NON-ideal aqueous solution: the exact derivative of the map \f$\mathbf{c}_1 \mapsto \mathbf{c}_{2,v}\f$
!> evaluated by compute_c2v_from_c1_Picard (mass action law with activity coefficients depending on the
!> ionic strength through the extended Debye-Hückel model of aq_phase_c%compute_log_act_coeffs_aq_phase).
!>
!> With \f$p=\log_{10}c_1\f$, \f$q=\log_{10}c_{2,v}\f$, \f$lg=\log_{10}\gamma\f$ the mass action law reads
!> \f$ q = S^*_{e,1}(p + lg_1) + \log_{10}K^* - lg_2\f$ and the aqueous \f$lg\f$ depend on
!> \f$ I = \tfrac12\sum_k z_k^2 c_k\f$ (all aqueous species). With \f$D_k = \partial lg_k/\partial I\f$:
!> \f[
!>   \frac{\partial q_i}{\partial p_j} = S^*_{ij} + W_i\frac{\partial I}{\partial p_j},\qquad
!>   W_i = \sum_k S^*_{ik}D_{1,k} - D_{2,i},\qquad
!>   \frac{\partial I}{\partial p_j} = \frac{\ln 10}{2}\Bigl[z_{1,j}^2c_{1,j} + \sum_l z_{2,l}^2c_{2,l}\frac{\partial q_l}{\partial p_j}\Bigr]
!> \f]
!> i.e. the linear system for \f$L = \partial q/\partial p\f$ (\f$n_e \times n_p\f$)
!> \f[
!>   \Bigl[\mathbf{I} - \tfrac{\ln 10}{2}\,W\,(z_2^2\circ c_2)^T\Bigr] L = S^*_{e,1} + \tfrac{\ln 10}{2}\,W\,(z_1^2\circ c_1)^T,
!>   \qquad \frac{\partial c_{2,i}}{\partial c_{1,j}} = \frac{c_{2,i}\,L_{ij}}{c_{1,j}},
!> \f]
!> solved with one LU factorisation (LAPACK dgesv, \f$n_p\f$ right-hand sides).
!>
!> Species ordering: c1, c2v and the rows/columns of Se_nc_1_star are in the speciation order of the
!> reactive zone (aqueous primaries, non-aqueous primaries | aqueous secondaries, non-aqueous
!> secondaries). The aqueous species are reached through this water's index maps (ind_prim_species,
!> ind_var_act_species, indices_aq_phase), NOT by position in the aqueous phase, so the result does not
!> depend on where water sits in the aqueous phase nor on a basis swap of the zone. Non-aqueous
!> variable activity species (free exchange site, variable activity minerals and gases, surface
!> complexes) and water have \f$D=0\f$ and \f$z^2=0\f$: their activity coefficients do not depend on I
!> and they do not contribute to it (water's log gamma is Raoult bookkeeping, fixed in the Picard
!> iteration). The ionic strength is the object's, i.e. that of the Picard-converged state at (c1,c2v).
subroutine compute_dc2v_dc1(this,c1,c2v,dc2v_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object, speciated at (c1,c2v): its ionic strength is used
    real(kind=8), intent(in) :: c1(:) !> primary species concentrations, speciation order (n_p)
    real(kind=8), intent(in) :: c2v(:) !> secondary variable activity species concentrations, speciation order (n_e)
    real(kind=8), intent(out) :: dc2v_dc1(:,:) !> Jacobian dc2v/dc1 (n_e x n_p, already allocated)
!> Variables
    integer(kind=4) :: i,j,n_p,n_p_aq,n_e,n_v2_aq,info
    integer(kind=4), allocatable :: ipiv(:) !> pivots of dgesv (n_e)
    real(kind=8) :: fac !> ln(10)/2
    real(kind=8), allocatable :: d_lg_dI_1(:) !> D of the primaries, speciation order, 0 for the non-aqueous ones (n_p)
    real(kind=8), allocatable :: d_lg_dI_2(:) !> D of the secondaries, speciation order, 0 for the non-aqueous ones (n_e)
    real(kind=8), allocatable :: z2c1(:) !> z^2*c of the primaries, 0 for the non-aqueous ones (n_p)
    real(kind=8), allocatable :: z2c2(:) !> z^2*c of the secondaries, 0 for the non-aqueous ones (n_e)
    real(kind=8), allocatable :: W(:) !> W_i = sum_k Se*_ik D_1k - D_2i (n_e)
    real(kind=8), allocatable :: mat_lin_syst(:,:) !> I - (ln10/2) W (z2^2 o c2)^T (n_e x n_e), overwritten by dgesv
    real(kind=8), allocatable :: L(:,:) !> Se* + (ln10/2) W (z1^2 o c1)^T on entry, L = dq/dp on exit (n_e x n_p)
    real(kind=8), parameter :: ln10=log(1d1)
!> Pre-process: ionic strength coupling terms of each species in speciation order
    associate(sa => this%solid_chemistry%reactive_zone%speciation_alg, &
              Se => this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star)
    n_p=sa%num_prim_species
    n_p_aq=sa%num_aq_prim_species
    n_e=sa%num_eq_reactions
    n_v2_aq=sa%num_aq_sec_var_act_species
    allocate(d_lg_dI_1(n_p),z2c1(n_p),d_lg_dI_2(n_e),z2c2(n_e),W(n_e),mat_lin_syst(n_e,n_e),L(n_e,n_p),ipiv(n_e))
    d_lg_dI_1=0d0; z2c1=0d0; d_lg_dI_2=0d0; z2c2=0d0 !> non-aqueous species: no ionic strength coupling
    do j=1,n_p_aq !> aqueous primaries (first n_p_aq columns of Se_nc_1_star)
        call species_terms(this%ind_prim_species(j),c1(j),d_lg_dI_1(j),z2c1(j))
    end do
    do i=1,n_v2_aq !> aqueous secondaries (first n_v2_aq rows of Se_nc_1_star)
        call species_terms(this%ind_var_act_species(n_p_aq+i),c2v(i),d_lg_dI_2(i),z2c2(i))
    end do
!> Process: linear system for L = dq/dp and conversion to dc2v/dc1
    W=matmul(Se,d_lg_dI_1)-d_lg_dI_2
    fac=0.5d0*ln10
    do i=1,n_e
        mat_lin_syst(i,:)=-fac*W(i)*z2c2
        mat_lin_syst(i,i)=mat_lin_syst(i,i)+1d0
        L(i,:)=Se(i,:)+fac*W(i)*z2c1
    end do
    call dgesv(n_e,n_p,mat_lin_syst,n_e,ipiv,L,n_e,info)
    if (info/=0) error stop "compute_dc2v_dc1: singular ionic strength coupling matrix"
    do j=1,n_p
        dc2v_dc1(:,j)=c2v*L(:,j)/c1(j)
    end do
    end associate
    deallocate(d_lg_dI_1,z2c1,d_lg_dI_2,z2c2,W,mat_lin_syst,L,ipiv)
contains
    !> D = d(log10 gamma)/dI and z^2*c of the aqueous species with concentration index k in this water:
    !> log gamma = -alpha*A*z^2*sqrt(I)/(1+beta*sqrt(I)) + gamma*I (aq_phase_c%compute_log_act_coeffs_aq_phase),
    !> so D = -alpha*A*z^2/(2*sqrt(I)*(1+beta*sqrt(I))^2) + gamma. Water: 0 (Raoult, fixed in the Picard iteration).
    subroutine species_terms(k,c,dlg,z2c)
        integer(kind=4), intent(in) :: k !> concentration index of the species in this water
        real(kind=8), intent(in) :: c !> its concentration
        real(kind=8), intent(out) :: dlg !> d(log10 gamma)/dI
        real(kind=8), intent(out) :: z2c !> z^2*c
        integer(kind=4) :: ka !> index in the aqueous phase
        real(kind=8) :: z2,sqrt_I
        ka=this%indices_aq_phase(k)
        if (ka==this%aq_phase%ind_wat) then
            dlg=0d0
            z2c=0d0
            return
        end if
        z2=dble(this%aq_phase%aq_species(ka)%valence**2)
        z2c=z2*c
        sqrt_I=sqrt(max(this%ionic_strength,0d0))
        associate(pac => this%aq_phase%aq_species(ka)%params_act_coeff)
        if (sqrt_I>0d0) then
            dlg=-pac%alpha*(this%params_aq_sol%A*z2)/(2d0*sqrt_I*(1d0+pac%get_beta(this%params_aq_sol)*sqrt_I)**2)+ &
                pac%get_gamma(this%params_aq_sol)
        else
            dlg=pac%get_gamma(this%params_aq_sol) !> I = 0: only uncharged species present, the Debye-Hückel term is killed by z2c = 0
        end if
        end associate
    end subroutine species_terms
end subroutine