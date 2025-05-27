import matplotlib.pyplot as plt

matrix_sizes = [100, 200, 400, 800, 1200, 1600, 2000, 4000, 10000]
seq_avg = [0.0003, 0.0010, 0.0041, 0.0138, 0.0336, 0.0552, 0.0856, 0.2890, 1.8595]
omp_avg = [0.0013, 0.0013, 0.0024, 0.0047, 0.0099, 0.0160, 0.0223, 0.0760, 0.3960]
speedup = [0.24, 0.75, 1.72, 2.95, 3.40, 3.44, 3.84, 3.80, 4.70]


# Kreiranje grafa
plt.figure(figsize=(10, 6))

# Vremena izvršavanja
plt.subplot(2, 1, 1)
plt.plot(matrix_sizes, seq_avg, marker="o", label="Sekvencijalno vrijeme")
plt.plot(matrix_sizes, omp_avg, marker="o", label="OMP vrijeme")
plt.title("Vremena izvršavanja u zavisnosti od veličine matrice")
plt.xlabel("Veličina matrice")
plt.ylabel("Prosječno vrijeme (s)")
plt.legend()
plt.grid(True)

# Speedup graf
plt.subplot(2, 1, 2)
plt.plot(matrix_sizes, speedup, marker="o", color="green")
plt.title("Ubrzanje (Speedup) u zavisnosti od veličine matrice")
plt.xlabel("Veličina matrice")
plt.ylabel("Speedup (x)")
plt.grid(True)

plt.tight_layout()

# Čuvanje grafa kao sliku
plt.savefig("graf_vremena_speedup.png")

# Prikazivanje grafa
plt.show()
