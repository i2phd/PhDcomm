#!/bin/bash

echo "Compilando il progetto Viterbi K=11 (versione finale e corretta) con g++..."

# Comando per compilare e linkare i file sorgente
g++ -o viterbi_test main.cpp viterbi_k11.cpp

# Controlla se la compilazione è riuscita
if [ $? -eq 0 ]; then
    echo "Compilazione completata con successo."
    echo "Eseguibile creato: viterbi_test"
else
    echo "Compilazione fallita."
fi