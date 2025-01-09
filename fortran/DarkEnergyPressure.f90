    module DarkEnergyPressure
    use precision
    use DarkEnergyInterface
    use interpolation
    use classes
    implicit none

    private
    type, extends(TDarkEnergyModel) :: TDarkEnergyDensityAndPressure
        real(dl) :: P_lam = -1._dl !pressure of DE
        real(dl) :: rho_lam = 1._dl !density of DE
        real(dl) :: cs2_lam = 1_dl !rest-frame sound speed of DE
        logical :: no_perturbations = .false. !Don't change this, no perturbations is unphysical
        Type(TCubicSpline) :: density, pressure
    contains
    procedure :: ReadParams => TDarkEnergyDensityAndPressure_ReadParams
    procedure :: Init => TDarkEnergyDensityAndPressure_Init
    procedure :: SetrhoTable => TDarkEnergyDensityAndPressure_SetrhoTable
    procedure :: SetPTable => TDarkEnergyDensityAndPressure_SetPTable
    procedure :: PrintFeedback => TDarkEnergyDensityAndPressure_PrintFeedback
    procedure :: grho_de => TDarkEnergyDensityAndPressure_grho_de
    procedure :: P_de => TDarkEnergyDensityAndPressure_P_de
    procedure :: Effective_w_wa => TDarkEnergyDensityAndPressure_Effective_w_wa
    end type TDarkEnergyDensityAndPressure

    public TDarkEnergyDensityAndPressure
    contains

    subroutine TDarkEnergyDensityAndPressure_SetrhoTable(this, a, rho, n)
    class(TDarkEnergyDensityAndPressure) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n), rho(n)

    if (abs(a(size(a)) -1) > 1e-5) error stop 'rho table must end at a=1'

    call this%density%Init(log(a), rho)
    this%rho_lam = rho(size(a))
    end subroutine TDarkEnergyDensityAndPressure_SetrhoTable

    subroutine TDarkEnergyDensityAndPressure_SetPTable(this, a, P, n)
    class(TDarkEnergyDensityAndPressure) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n), P(n)

    if (abs(a(size(a)) -1) > 1e-5) error stop 'P table must end at a=1'
    
    call this%pressure%Init(log(a), P)

    this%P_lam = P(size(a))

    end subroutine TDarkEnergyDensityAndPressure_SetPTable


    function TDarkEnergyDensityAndPressure_P_de(this, a) result(P_de)
    class(TDarkEnergyDensityAndPressure) :: this
    real(dl) :: P_de, al, al_larger, al_smaller, w_rhop
    real(dl), intent(IN) :: a
    
    al=dlog(a)
    if(this%density%X(1) > this%pressure%X(1)) then
        al_larger = this%density%X(1)
    else
        al_larger = this%pressure%X(1)
    endif
    if(this%density%X(this%density%n) < this%pressure%X(this%pressure%n)) then
        al_smaller = this%density%X(this%density%n)
    else
        al_smaller = this%pressure%X(this%pressure%n) 
    endif
    if(al <= al_larger) then
        w_rhop = this%pressure%Value(al_larger)/this%density%Value(al_larger)
        P_de= this%grho_de(a)/(a ** 4)*w_rhop
    elseif(al >= al_smaller) then
        w_rhop = this%pressure%Value(al_smaller)/this%density%Value(al_smaller)
        P_de= this%grho_de(a)/(a ** 4)*w_rhop
    else
        P_de = this%pressure%Value(al)
    endif
    end function TDarkEnergyDensityAndPressure_P_de  ! pressure of the PPF DE

    function TDarkEnergyDensityAndPressure_grho_de(this, a) result(grho_de)
    class(TDarkEnergyDensityAndPressure) :: this
    real(dl) :: grho_de, al, w_rhop, al_larger, al_smaller, fint
    real(dl), intent(IN) :: a
    
    if(a == 0.d0) then
        grho_de = 0.d0
    else
        if (a>=1) then
            fint = 1
        else
            al=dlog(a)
            if(this%density%X(1) > this%pressure%X(1)) then
                al_larger = this%density%X(1)
            else
                al_larger = this%pressure%X(1)
            endif
            if(this%density%X(this%density%n) < this%pressure%X(this%pressure%n)) then
                al_smaller = this%density%X(this%density%n)
            else
                al_smaller = this%pressure%X(this%pressure%n) 
            endif
            if(al < al_larger) then
                w_rhop = this%pressure%Value(al_larger)/this%density%Value(al_larger)
                !the exp(al_larger) **4 is required because the density can be seen as a_larger^-3(1+w) and you need an extra
                !a_larger^4 to make it such that it fits with the other parts.
                fint = this%density%F(1) * exp(al_larger) **4 *exp((1. - 3. * w_rhop) * (al - al_larger))
            elseif(al > al_smaller) then
                w_rhop = this%pressure%Value(al_smaller)/this%density%Value(al_smaller)
                fint = this%density%F(this%density%n) * exp(al_smaller) **4 * exp((1. - 3. * w_rhop) * (al - al_smaller))
            else
                fint = this%density%Value(al) * a ** 4._dl
            endif
        endif
        grho_de = fint
    endif
    end function TDarkEnergyDensityAndPressure_grho_de  ! pressure of the PPF DE

    subroutine TDarkEnergyDensityAndPressure_Effective_w_wa(this, w, wa)
    class(TDarkEnergyDensityAndPressure), intent(inout) :: this
    real(dl), intent(out) :: w, wa
    w = this%P_lam/this%rho_lam
    wa = - (this%pressure%Derivative(0._dl)*this%rho_lam - this%density%Derivative(0._dl)*this%P_lam)/(this%rho_lam**2)
    print *, w, wa
    end subroutine TDarkEnergyDensityAndPressure_Effective_w_wa

   ! function TDarkEnergyDensityAndPressure_w_de(this, a) result(w_de)
   ! class(TDarkEnergyDensityAndPressure) :: this
   ! real(dl) :: w_de, al, fint
   ! real(dl), intent(IN) :: a

   ! al = dlog(a)
   ! if(al <= this%density%Xmin_interp) then
   !     fint = this%density%F(1)
    
    
    subroutine TDarkEnergyDensityAndPressure_PrintFeedback(this, FeedbackLevel)
    class(TDarkEnergyDensityAndPressure) :: this
    integer, intent(in) :: FeedbackLevel

    end subroutine TDarkEnergyDensityAndPressure_PrintFeedback

    subroutine TDarkEnergyDensityAndPressure_ReadParams(this, Ini)
    use IniObjects
    use FileUtils
    class(TDarkEnergyDensityAndPressure) :: this
    class(TIniFile), intent(in) :: Ini
    real(dl), allocatable :: table(:,:)
    
    call File%LoadTxt(Ini%Read_String('rhoafile'), table)
    call this%SetrhoTable(table(:,1),table(:,2), size(table(:,1)))

    call File%LoadTxt(Ini%Read_String('Pafile'), table)
    
    call this%SetPTable(table(:,1),table(:,2), size(table(:,1)))

    end subroutine TDarkEnergyDensityAndPressure_ReadParams


    subroutine TDarkEnergyDensityAndPressure_Init(this, State)
    use classes
    class(TDarkEnergyDensityAndPressure), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    this%is_cosmological_constant = .false.
    this%is_no_mod_P = .false.

    end subroutine TDarkEnergyDensityAndPressure_Init

    end module DarkEnergyPressure
