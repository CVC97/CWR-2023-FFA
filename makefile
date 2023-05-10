# Flags
FLAGS := -lm -lgsl -lgslcblas -Wall -std=c99


# Executable
# EXECUTABLE := $(find . -type f -executable)

# Argumente, mit dem das Programm ausgeführt wird,
# falls das Programm keine Argumente aus der Konsole übernimmt, kann das leer gelassen werden
# ARGS :=


# EXECUTABLE: EXECUTABLE.o cvc_numerics.o
# 	gcc EXECUTABLE.o cvc_numerics.o -o EXECUTABLE $(FLAGS)

# EXECUTABLE.o: EXECUTABLE.c
# 	gcc -c EXECUTABLE.c

# cvc_numerics.o: cvc_numerics.c
# 	gcc -c cvc_numerics.c


# A01
A01_Einfache_Integration: A01_Einfache_Integration.c
	gcc A01_Einfache_Integration.c -o A01_Einfache_Integration $(FLAGS)


# A02
A02_Modularisieren: A02_Modularisieren.c
	gcc A02_Modularisieren.c -o A02_Modularisieren $(FLAGS)


# A03
A03_Fehlerfunktion: A03_Fehlerfunktion.o cvc_numerics.o
	gcc A03_Fehlerfunktion.o cvc_numerics.o -o A03_Fehlerfunktion $(FLAGS)
	rm *.o

A03_Fehlerfunktion.o: A03_Fehlerfunktion.c
	gcc -c A03_Fehlerfunktion.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c


# A05
A05_Numerische_Ableitungen: A05_Numerische_Ableitungen.c
	gcc A05_Numerische_Ableitungen.c -o A05_Numerische_Ableitungen $(FLAGS)


# A06
A06_NewtonRaphson: A06_NewtonRaphson.o cvc_numerics.o
	gcc A06_NewtonRaphson.o cvc_numerics.o -o A06_NewtonRaphson $(FLAGS)
	rm *.o

A06_NewtonRaphson.o: A06_NewtonRaphson.c
	gcc -c A06_NewtonRaphson.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c


# A07
A07_Fourier_Integral: A07_Fourier_Integral.o cvc_numerics.o
	gcc A07_Fourier_Integral.o cvc_numerics.o -o A07_Fourier_Integral $(FLAGS)
	rm *.o

A07_Fourier_Integral.o: A07_Fourier_Integral.c
	gcc -c A07_Fourier_Integral.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c


# A08
A08_Gleitkommagenauigkeit: A08_Gleitkommagenauigkeit.c
	gcc A08_Gleitkommagenauigkeit.c -o A08_Gleitkommagenauigkeit $(FLAGS)


# A09
A09_Epidemie_Modellierung: A09_Epidemie_Modellierung.o cvc_numerics.o
	gcc A09_Epidemie_Modellierung.o cvc_numerics.o -o A09_Epidemie_Modellierung $(FLAGS)
	rm *.o

A09_Epidemie_Modellierung.o: A09_Epidemie_Modellierung.c
	gcc -c A09_Epidemie_Modellierung.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c



# A10
A10_Arrays_und_Pointer: A10_Arrays_und_Pointer.c
	gcc A10_Arrays_und_Pointer.c -o A10_Arrays_und_Pointer $(FLAGS)

# A11
A11_Populationsdynamik: A11_Populationsdynamik.c
	gcc A11_Populationsdynamik.c -o A11_Populationsdynamik $(FLAGS)


# A12
A12_pendulums: A12_pendulums.o cvc_numerics.o
	gcc A16_pendulums.o cvc_numerics.o -o A12_pendulums $(FLAGS)
	rm *.o

A12_pendulums.o: A12_pendulums.c
	gcc -c A16_pendulums.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c


# A13
A13_Euler_Integration: A13_Euler_Integration.o cvc_numerics.o
	gcc A13_Euler_Integration.o cvc_numerics.o -o A13_Euler_Integration $(FLAGS)
	rm *.o

A13_Euler_Integration.o: A13_Euler_Integration.c
	gcc -c A13_Euler_Integration.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c


# A14
A14_Schwingungssensor: A14_Schwingungssensor.o cvc_numerics.o
	gcc A14_Schwingungssensor.o cvc_numerics.o -o A14_Schwingungssensor $(FLAGS)
	rm *.o

A14_Schwingungssensor.o: A14_Schwingungssensor.c
	gcc -c A14_Schwingungssensor.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c



# A15
A15_2D_MC_Integration: A15_2D_MC_Integration.o cvc_numerics.o
	gcc A15_2D_MC_Integration.o cvc_numerics.o -o A15_2D_MC_Integration $(FLAGS)
	rm *.o

A15_2D_MC_Integration.o: A15_2D_MC_Integration.c
	gcc -c A15_2D_MC_Integration.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c



# A18
A18_RK2_Pendel: A18_RK2_Pendel.o cvc_numerics.o
	gcc A18_RK2_Pendel.o cvc_numerics.o -o A18_RK2_Pendel $(FLAGS)
	rm *.o

A18_RK2_Pendel.o: A18_RK2_Pendel.c
	gcc -c A18_RK2_Pendel.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c


# A19
A19_Konvergenz_RK2: A19_Konvergenz_RK2.o cvc_numerics.o
	gcc A19_Konvergenz_RK2.o cvc_numerics.o -o A19_Konvergenz_RK2 $(FLAGS)
	rm *.o

A19_Konvergenz_RK2.o: A19_Konvergenz_RK2.c
	gcc -c A19_Konvergenz_RK2.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c


# A21
A21_Verlet_Integrator: A21_Verlet_Integrator.o cvc_numerics.o
	gcc A21_Verlet_Integrator.o cvc_numerics.o -o A21_Verlet_Integrator $(FLAGS)
	rm *.o

A21_Verlet_Integrator.o: A21_Verlet_Integrator.c
	gcc -c A21_Verlet_Integrator.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c


# (old) A12
old_A12_Polarmethode_normalverteilter_Zufallszahlen: old_A12_Polarmethode_normalverteilter_Zufallszahlen.c
	gcc old_A12_Polarmethode_normalverteilter_Zufallszahlen.c -o old_A12_Polarmethode_normalverteilter_Zufallszahlen $(FLAGS)


# A13
old_A13_Histogramme: old_A13_Histogramme.o cvc_numerics.o
	gcc old_A13_Histogramme.o cvc_numerics.o -o old_A13_Histogramme $(FLAGS)
	rm *.o

old_A13_Histogramme.o: old_A13_Histogramme.c
	gcc -c old_A13_Histogramme.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c


# A14
old_A14_Monte_Carlo_Integration: old_A14_Monte_Carlo_Integration.o cvc_numerics.o
	gcc A14_Monte_Carlo_Integration.o cvc_numerics.o -o old_A14_Monte_Carlo_Integration $(FLAGS)
	rm *.o

old_A14_Monte_Carlo_Integration.o: old_A14_Monte_Carlo_Integration.c
	gcc -c Aold_14_Monte_Carlo_Integration.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c


# old_A20
old_A20_Reibung: old_A20_Reibung.o cvc_numerics.o
	gcc old_A20_Reibung.o cvc_numerics.o -o old_A20_Reibung $(FLAGS)
	rm *.o

old_A20_Reibung.o: old_A20_Reibung.c
	gcc -c old_A20_Reibung.c

cvc_numerics.o: cvc_numerics.c
	gcc -c cvc_numerics.c


# Commands
o_delete:
	rm *.o


# Update die Numerik-Bibliotheken im submissions-Folder
update_numerics:
	cp cvc_numerics.c submissions/cvc_numerics.c
	cp cvc_numerics.h submissions/cvc_numerics.h
