import numpy as np
# x : cAMP concentration (pmol/mg protein)

def update_pka(x):

    x = np.array(x, dtype=float)

    den = 300.0 + 15.0 * x + x**2

    sqrt_arg = (
        -1643644331734924815.0 * x**6
        - 19572532067205238230.0 * x**5
        - 393461126229536661375.0 * x**4
        - 80438302732946971500.0 * x**3
        - 805328393113024740000.0 * x**2
        - 22688778805320600000.0 * x
        - 151258525368804000000.0
    )

    sqrt_term = np.sqrt(sqrt_arg + 0j)  

    A = (
        -468891428493892500.0 * x**2
        - 46425728249607375.0 * x**3
        - 14479420376562975.0 * x**4
        - 666054712166085.0 * x**5
        - 47282992935913.0 * x**6
        - 741463359651000000.0
        - 111219503947650000.0 * x
        + 9000000.0 * x * sqrt_term
        + 450000.0 * x**2 * sqrt_term
        + 30000.0 * x**3 * sqrt_term
    )

    cube_A = A ** (1.0 / 3.0)

    coeff = 3.3333333333333333e-05  # 0.3333...e-4

    term1 = coeff * cube_A / den

    term2_num = (
        349407388425.0 * x**2
        + 17367968670.0 * x**3
        + 1548762289.0 * x**4
        + 819206010000.0
        + 81920601000.0 * x
    )
    term2_den = den * cube_A
    term2 = coeff * term2_num / term2_den

    term3_num = 905100.0 + 45255.0 * x + 23017.0 * x**2
    term3 = coeff * term3_num / den

    PKA = term1 + term2 - term3

    PKA = np.abs(PKA)

    return np.real_if_close(PKA)
