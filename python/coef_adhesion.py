import numpy as np
from scipy.integrate import quad

# Définir la fonction à intégrer
def integrand(x, h):
    return ((- 4 * (x**2) / h + 6 * x - 2 * h)**(1/4)) * x

# Définir la valeur de h (vous pouvez la changer pour vos besoins spécifiques)
s = 0.05
h = 1.2*s

# Effectuer l'intégration numérique
result, error = quad(integrand, h/2, h, args=(h,))

print(f"Résultat de l'intégrale: {1/(result*2*np.pi)}")
print(f"Erreur estimée: {error}")


print(0.007/(np.power(h,3.25)))