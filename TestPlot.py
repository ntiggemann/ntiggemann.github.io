import numpy as np
import matplotlib.pyplot as plt

# Bereich für x definieren
x = np.linspace(-2, 2, 400)

# Funktion f(x) = x^2 definieren
y = x**2

# Plotten
plt.plot(x, y, label="f(x) = x^2")
plt.xlim(-2, 2)  # Begrenzung der x-Achse
plt.ylim(0, 4)   # Begrenzung der y-Achse
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Plot der Funktion f(x) = x^2")
plt.grid(True)
plt.legend()
plt.show()
