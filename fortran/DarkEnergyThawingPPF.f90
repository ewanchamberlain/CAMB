module DarkEnergyThawingPPF
    use DarkEnergyPPF
    use DarkEnergyInterface
    use classes
    implicit none

    private

    type, extends(TDarkEnergyPPF) :: TDarkEnergyThawingPPF
    
    contains
        procedure, nopass :: PythonClass => TDarkEnergyThawingPPF_PythonClass
        procedure, nopass :: SelfPointer => TDarkEnergyThawingPPF_SelfPointer
        procedure :: w_de => TDarkEnergyThawingPPF_w_de
        procedure :: grho_de => TDarkEnergyThawingPPF_grho_de
    end type TDarkEnergyThawingPPF

    public TDarkEnergyThawingPPF
    contains

    subroutine TDarkEnergyThawingPPF_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TDarkEnergyThawingPPF), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType) 
    P => PType

    end subroutine TDarkEnergyThawingPPF_SelfPointer

    function TDarkEnergyThawingPPF_PythonClass()
        character(LEN=:), allocatable :: TDarkEnergyThawingPPF_PythonClass

        TDarkEnergyThawingPPF_PythonClass = 'DarkEnergyThawingPPF'
    end function TDarkEnergyThawingPPF_PythonClass

    function TDarkEnergyThawingPPF_w_de(this, a)
    class(TDarkEnergyThawingPPF) :: this
    real(dl) :: TDarkEnergyThawingPPF_w_de, al
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w) then
        TDarkEnergyThawingPPF_w_de= max(-1., this%w_lam+ this%wa*(1._dl-a))
    else
        al=dlog(a)
        if(al <= this%equation_of_state%Xmin_interp) then
            TDarkEnergyThawingPPF_w_de= max(-1._dl, this%equation_of_state%F(1))
        elseif(al >= this%equation_of_state%Xmax_interp) then
            TDarkEnergyThawingPPF_w_de= max(-1._dl, this%equation_of_state%F(this%equation_of_state%n))
        else
            TDarkEnergyThawingPPF_w_de = max(-1._dl, this%equation_of_state%Value(al))
        endif
    endif

    end function TDarkEnergyThawingPPF_w_de  ! equation of state of the PPF DE
    
    function TDarkEnergyThawingPPF_grho_de(this, a) result(grho_de) !relative density (8 pi G a^4 rho_de /grhov)
    class(TDarkEnergyThawingPPF) :: this
    real(dl) :: grho_de, al, fint, a0
    real(dl), intent(IN) :: a

    ! TODO: implement case where DE leaves CPL and enters Lambda
    if(.not. this%use_tabulated_w) then
        a0 = (1._dl + this%w_lam + this%wa)/this%wa ! scale factor at which w=-1
        if ((a0 >= 1) .or. (a0 < 0)) then
            if ((this%w_lam + this%wa) < -1) then
                grho_de = a**4
            else
                grho_de = a**4 * (a) ** (-3._dl - 3. * this%w_lam - 3 * this%wa)
                if (this%w_lam/=0) grho_de=grho_de*exp(-3. * this%wa * (1._dl - a))
            endif
        else if ((this%w_lam + this%wa) <= -1) then
            if (a <= a0) then
                grho_de = a**4 * a0 ** (-3._dl - 3. * this%w_lam - 3. * this%wa) * exp(-3. * this%wa * (1._dl - a0))
            else
                grho_de = a ** (1._dl - 3. * this%w_lam - 3. * this%wa)
                if (this%wa/=0) grho_de=grho_de*exp(-3. * this%wa * (1._dl - a))
            endif
        else
            if (a <= a0) then
                grho_de = a**4 * (a/a0) ** (-3._dl - 3. * this%w_lam - 3. * this%wa)
                if (this%wa/=0) grho_de=grho_de*exp(-3. * this%wa * (a0 - a))
            else
                grho_de = a**4
            endif
        endif
    else
        if(a == 0.d0)then
            grho_de = 0.d0      !assume rho_de*a^4-->0, when a-->0, OK if w_de always <0.
        else
            if (a>=1) then
                fint= 1
            else
                al = dlog(a)
                if(al <= this%logdensity%X(1)) then
                    ! assume here w=w_de(a_min)
                    fint = exp(this%logdensity%F(1) + (1. - 3. * this%equation_of_state%F(1))*(al - this%logdensity%X(1)))
                else
                    fint = exp(this%logdensity%Value(al))
                endif
            end if
            grho_de = fint
        endif
    endif

    end function TDarkEnergyThawingPPF_grho_de
end module DarkEnergyThawingPPF    


