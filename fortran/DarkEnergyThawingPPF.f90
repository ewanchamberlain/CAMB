module DarkEnergyThawingPPF
    use DarkEnergyPPF
    use classes
    implicit none

    private

    type, extends(TDarkEnergyPPF) :: TDarkEnergyThawingPPF
        real(dl) :: w_lam = -1.0_dl
        real(dl) :: wa = 0.0_dl
    
    contains
        procedure :: ReadParams => TDarkEnergyThawingPPF_ReadParams
        procedure, nopass :: PythonClass => TDarkEnergyThawingPPF_PythonClass
        procedure :: Init => TDarkEnergyThawingPPF_Init
        procedure :: TDarkEnergyEqnOfState_w_de => TDarkEnergyThawingPPF_w_de
        procudure, nopass :: SelfPointer => TDarkEnergyThawingPPF_SelfPointer
    end type TDarkEnergyThawingPPF

    public TDarkEnergyThawingPPF
    contains

    subroutine TDarkEnergyThawingPPF_ReadParams(this, Ini)
        use IniObjects
        class(TDarkEnergyThawingPPF) :: this
        class(TIniFile), intent(in) :: Ini

        call this%TDarkEnergyPPF%ReadParams(Ini)
        this%w_lam = Ini%Read_Double('w_lam', -1.d0)
        this%wa = Ini%Read_Double('wa', 0.d0)
        call this%setcgammappf

    end subroutine TDarkEnergyThawingPPF_ReadParams

    subroutine TDarkEnergyThawingPPF_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TDarkEnergyPPF), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType

    end subroutine TDarkEnergyThawingPPF_SelfPointer

    function TDarkEnergyThawingPPF_PythonClass()
        character(LEN=:), allocatable :: TDarkEnergyThawingPPF_PythonClass

        TDarkEnergyThawingPPF_PythonClass = 'DarkEnergyThawingPPF'
    end function TDarkEnergyThawingPPF_PythonClass


    function TDarkEnergyThawingPPF_w_de(this, a)
    class(TDarkEnergyEqnOfState) :: this
    real(dl) :: TDarkEnergyEqnOfState_w_de, al
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w) then
        TDarkEnergyEqnOfState_w_de= max(-1._dl, this%w_lam+ this%wa*(1._dl-a))
    else
        al=dlog(a)
        if(al <= this%equation_of_state%Xmin_interp) then
            TDarkEnergyEqnOfState_w_de= this%equation_of_state%F(1)
        elseif(al >= this%equation_of_state%Xmax_interp) then
            TDarkEnergyEqnOfState_w_de= this%equation_of_state%F(this%equation_of_state%n)
        else
            TDarkEnergyEqnOfState_w_de = this%equation_of_state%Value(al)
        endif
    endif

    end function TDarkEnergyEqnOfState_w_de  ! equation of state of the PPF DE
    
    


