# Power Iteration algorithm(algoritam stepena iteracije): Paralelna Implementacija i upoređivanje sa sekvencijalnom

## 🧾 Opis projekta

Ovaj projekat predstavlja implementaciju **Power Iteration algoritma**, odnosno algoritma stepene iteracije za pronalaženje **dominantne sopstvene vrijednosti i sopstvenog vektora** kvadratne matrice. Algoritam je implementiran u dvije varijante:

- ✅ **Sekvencijalna verzija** – koristi osnovne petlje i jedno jezgro
- ⚡ **Paralelna verzija** – koristi OpenMP za višenedno (multi-threaded) izvršavanje

Cilj projekta je da se uporede performanse ovih pristupa na matricama različitih dimenzija.

---

## ⚙️ Kako pokrenuti projekat

### 1. 📦 Zahtjevi

- GNU `gcc` sa podrškom za OpenMP (npr. `gcc-14`)
- **macOS korisnici**: preporučuje se instalacija GCC preko **Homebrew-a**:

```bash
brew install gcc
```

- **Windows korisnici**: MinGW + MSYS2 (Preporučeno)
- Preuzmi i instaliraj MSYS2.
- Otvori MSYS2 MSYS terminal i izvrši:

```bash
pacman -Syu       # Ažuriranje paketa
pacman -S mingw-w64-ucrt-x86_64-gcc
```

### 2. 🔨 Kompajliranje
Otvorite terminal u direktorijumu projekta i izvršite sledeću komandu:
```bash
gcc-13 -fopenmp Implementacija.c -o power_benchmark
```

3. ▶️ Pokretanje
Izvršite sledeću komandu:
```bash
./power_benchmark
```
## 📊 Izlaz

Program će izvršiti više testova nad matricama različitih dimenzija i ispisati prosječna vremena i ubrzanje (speedup)