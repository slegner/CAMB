    module DarkEnergyPressure
    use precision
    use interpolation
    use classes
    implicit none

    private
    type, extends(TDarkEnergyModel) :: TDarkEnergyDensityAndPressure
        !Type supporting general rho(z) and P(z) table
        logical :: use_tabulated_rho  !Use interpolated table; note this is quite slow.
        logical :: use_tabulated_P 
        !Interpolations.
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

    public TDarkEnergyModel, TDarkEnergyDensityAndPressure
    contains


    subroutine TDarkEnergyDensityAndPressure_SetrhoTable(this, a, rho, n)
    class(TDarkEnergyDensityAndPressure) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n), rho(n)

    if (abs(a(size(a)) -1) > 1e-5) error stop 'rho table must end at a=1'

    call this%density%Init(log(a), rho)
    end subroutine TDarkEnergyDensityAndPressure_SetrhoTable

    subroutine TDarkEnergyDensityAndPressure_SetPTable(this, a, P, n)
    class(TDarkEnergyDensityAndPressure) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n), P(n)

    if (abs(a(size(a)) -1) > 1e-5) error stop 'P table must end at a=1'

    call this%pressure%Init(log(a), rho)
    end subroutine TDarkEnergyDensityAndPressure_SetPTable


    function TDarkEnergyDensityAndPressure_P_de(this, a)
    class(TDarkEnergyDensityAndPressure) :: this
    real(dl) :: TDarkEnergyDensityAndPressure_P_de, al
    real(dl), intent(IN) :: a

    al=dlog(a)
    if(al <= this%pressure%Xmin_interp) then
        TDarkEnergyDensityAndPressure_P_de= this%pressure%F(1)
    elseif(al >= this%pressure%Xmax_interp) then
        TDarkEnergyDensityAndPressure_P_de= this%pressure%F(this%pressure%n)
    else
        TDarkEnergyDensityAndPressure_P_de = this%pressure%Value(al)
    endif

    end function TDarkEnergyDensityAndPressure_P_de  ! pressure of the PPF DE

    function TDarkEnergyDensityAndPressure_grho_de(this, a)
    class(TDarkEnergyDensityAndPressure) :: this
    real(dl) :: TDarkEnergyDensityAndPressure_grho_de, al
    real(dl), intent(IN) :: a

    al=dlog(a)
    if(al <= this%density%Xmin_interp) then
        TDarkEnergyDensityAndPressure_grho_de= this%density%F(1)
    elseif(al >= this%density%Xmax_interp) then
        TDarkEnergyDensityAndPressure_grho_de= this%density%F(this%density%n)
    else
        TDarkEnergyDensityAndPressure_grho_de = this%density%Value(al)
    endif

    end function TDarkEnergyDensityAndPressure_grho_de  ! pressure of the PPF DE

    subroutine TDarkEnergyDensityAndPressure_PrintFeedback(this, FeedbackLevel)
    class(TDarkEnergyDensityAndPressure) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel >0) write(*,'("(w0, wa) = (", f8.5,", ", f8.5, ")")') &
        &   this%w_lam, this%wa

    end subroutine TDarkEnergyDensityAndPressure_PrintFeedback

    subroutine TDarkEnergyDensityAndPressure_ReadParams(this, Ini)
    use IniObjects
    use FileUtils
    class(TDarkEnergyDensityAndPressure) :: this
    class(TIniFile), intent(in) :: Ini
    real(dl), allocatable :: table(:,:)

    call File%LoadTxt(Ini%Read_String('rhoafile'), table)
    call File%LoadTxt(Ini%Read_String('Pafile'), table)

    call this%SetPTable(table(:,1),table(:,2), size(table(:,1)))
    call this%SetrhoTable(table(:,1),table(:,2), size(table(:,1)))

    end subroutine TDarkEnergyDensityAndPressure_ReadParams


    subroutine TDarkEnergyDensityAndPressure_Init(this, State)
    use classes
    class(TDarkEnergyDensityAndPressure), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    this%is_cosmological_constant = .false.

    end subroutine TDarkEnergyDensityAndPressure_Init


    end module DarkEnergyDensityAndPressure
