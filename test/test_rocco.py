import math

Sqrtln2 = math.sqrt(math.log(2))
OneminSqrtPIln2 = 1 - math.sqrt(math.pi * math.log(2))

# simple erfcx
def erfcx(x):
    return math.exp(x**2) * math.erfc(x)

def PseudoRoccoVoigt(alpha, gamma):
    for y in [0.1, 1.0, 5.0, 10.0, 20.0, 50.0]:
        exponent = -0.6055*y + 0.0718*y**2 - 0.0049*y**3 + 0.000136*y**4
        print(f"y: {y}, exponent: {exponent}")
        try:
            bhalfy = y + Sqrtln2 * math.exp(exponent)
            Vy = bhalfy * erfcx(y)
            eta = (Vy - Sqrtln2) / (Vy * OneminSqrtPIln2)
            print(f"  eta: {eta}")
        except Exception as e:
            print(f"  ERROR: {e}")

PseudoRoccoVoigt(0.003, 0.05)
