# ---------------------------------------------------------
# This example can be found in
# [Chapter 9, in Finite Element Procedures, K-J Bathe].
# ---------------------------------------------------------
function oscillator()
    # analytic solution (cf. Example 9.7)
    A = [1/√3  (1/2)*√(2/3);
         1/√3      -√(2/3)]
    x₁(t) = (5 / √3) * (1 - cos(t*√2))
    x₂(t) = (2 * √(2/3)) * (-1 + cos(t*√5))
    U(t) = A * [x₁(t), x₂(t)]

    U12 = U(12*0.28) # ≈ [1.157, 2.489]

    # problem formulation
    M = [2 0; 0 1.]
    K = [6 -2; -2 4.]
    C = zeros(2, 2)
    R = [0, 10.]
    SecondOrderAffineContinuousSystem(M, C, K, R)
end

# ---------------------------------------------------------
# Heat transfer problem
# ---------------------------------------------------------
function heat_transfer()
    rho = 1. ;
    csh = 1. ;
    kco = 4. ;
    nelem = 3 ;
    lelem = 1.0 / nelem ;
    Area  = .25 ;
    C = rho * csh * lelem * Area * ( [0.5+0.5+0.0 0; 0  0.0+0.5+0.5] )
    K = kco *       lelem * Area * ( [1.0+1.0 -1.; -1. 1.0+1.0] * 1.0 / lelem^2 )
    M = zeros(2, 2)
    R = zeros(2)
    SecondOrderAffineContinuousSystem(M, C, K, R)
end
