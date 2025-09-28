#!/bin/bash

# Questo script compila il progetto utilizzando g++, un compilatore C++ moderno
# e ampiamente disponibile. È usato come sostituto del Borland Classic C++ Compiler (bcc32),
# che non è disponibile in questo ambiente. Il codice è scritto in C++98 per
# massimizzare la compatibilità.

echo "Compilando il progetto con g++..."

# Comando per compilare e linkare i file sorgente con g++
g++ -std=c++98 -Wall -o encoder_test main.cpp convolutional_encoder.cpp

# Controlla se la compilazione è riuscita
if [ $? -eq 0 ]; then
    echo "Compilazione completata con successo."
    echo "Eseguibile creato: encoder_test"
else
    echo "Compilazione fallita."
fi