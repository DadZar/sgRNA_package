import numpy as np
import pandas as pd
import joblib

import tensorflow as tf

# Cargar el modelo entrenado
try:
    rf_model = joblib.load("modelo_sgRNA.pkl")
    print("Modelo RandomForest cargado correctamente")
except Exception as e:
    print(f"Error cargando modelo RandomForest: {e}")

try:
    nn_model = tf.keras.models.load_model("modelo_sgRNA_nn.h5", compile=False)
    nn_model.compile(optimizer='adam', loss=tf.keras.losses.MeanSquaredError(), metrics=['mae'])
    print("Modelo TensorFlow cargado correctamente")
except Exception as e:
    print(f"Error cargando modelo TensorFlow: {e}")

try:
    xgb_model = joblib.load("modelo_sgRNA_xgb.pkl")
    print("Modelo XGBoost cargado correctamente")
except Exception as e:
    print(f"Error cargando modelo XGBoost: {e}")
    


# Función para obtener la secuencia complemento inversa
def complemento_inverso(seq):
    """Devuelve la secuencia complemento inversa de una secuencia de ADN."""
    complementos = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complementos[base] for base in reversed(seq))

# Función para codificar la secuencia
def one_hot_encode(seq):
    mapping = {'A': [1,0,0,0], 'C': [0,1,0,0], 'G': [0,0,1,0], 'T': [0,0,0,1]}
    return np.array([mapping[nuc] for nuc in seq]).flatten()

# Función para obtener la secuencia complemento inversa
def complemento_inverso(seq):
    """Devuelve la secuencia complemento inversa de una secuencia de ADN."""
    complementos = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complementos[base] for base in reversed(seq))

# Función para predecir la eficiencia de una guía sgRNA
def predecir_eficiencia_guia(secuencia_20nt):
    # 1. Verificar longitud
    if len(secuencia_20nt) != 20:
        raise ValueError("La secuencia de la guía debe tener 20 nucleótidos")

    # 2. Características simples de la secuencia
    gc = (secuencia_20nt.count('G') + secuencia_20nt.count('C')) / 20.0  # Fracción GC
    ultima_base_G = 1 if secuencia_20nt[-1] == 'G' else 0

    # 3. Inicializar puntaje de eficiencia
    puntaje = 1.0  # Aumento el puntaje inicial para evitar valores negativos

    # 4. Aplicar reglas heurísticas
    # Penalización si el GC está fuera del rango óptimo (40-80%)
    if gc < 0.40 or gc > 0.80:
        puntaje -= 0.3  # Penalización más suave

    # Sumar puntos si la última base es G (factor favorable)
    if ultima_base_G:
        puntaje += 0.5

    # 5. Incorporar modelo posicional (ejemplo basado en Rule Set 1)
    peso_posicional = {
        18: {'C': +0.2},        # C en pos18 aporta +0.2
        19: {'G': +0.2},        # G en pos19 aporta +0.2
        20: {'G': +0.3, 'T': -0.3}  # G en pos20 +0.3, T en pos20 -0.3
    }
    for pos in range(1, 21):  # Posiciones 1 a 20 de la guía
        base = secuencia_20nt[pos-1]  # Ajuste de índice
        if pos in peso_posicional and base in peso_posicional[pos]:
            puntaje += peso_posicional[pos][base]

    # 6. Evaluar auto-complementariedad (ejemplo: contar "GCGC" o "CGCG")
    posibles_palindromos = 0
    for i in range(0, 17):  # Hasta 17 para comparar tramos de 4 bases
        tramo = secuencia_20nt[i:i+4]
        comp = complemento_inverso(tramo)
        if comp in secuencia_20nt:
            posibles_palindromos += 1

    # Límite máximo de penalización por estructuras secundarias
    max_penalizacion_palindromos = 0.8  # No restar más de 0.8 en total
    penalizacion_real = min(0.3 * posibles_palindromos, max_penalizacion_palindromos)
    puntaje -= penalizacion_real

    # 7. Asegurar que el puntaje no sea negativo
    puntaje = max(0, puntaje)

    return puntaje


# Función para predecir la eficiencia de una nueva secuencia sgRNA
def predecir_eficiencia(seq):
    if len(seq) != 20:
        raise ValueError("La secuencia debe tener exactamente 20 nucleótidos.")

    # Codificar la secuencia con one-hot encoding
    X_new = np.array([one_hot_encode(seq)])

    # Hacer la predicción con el modelo entrenado
    prediccion = rf_model.predict(X_new)[0]

    prediccion =prediccion*100

    return prediccion



# Función de predicción con Red Neuronal
def predecir_eficiencia_nn(seq):
    if len(seq) != 20:
        raise ValueError("La secuencia debe tener exactamente 20 nucleótidos.")

    # Codificar la secuencia
    X_new = np.array([one_hot_encode(seq)])

    # Hacer la predicción
    prediccion = nn_model.predict(X_new)[0][0] * 100

    return prediccion



# Función de predicción con XGBoost
def predecir_eficiencia_xgb(seq):
    if len(seq) != 20:
        raise ValueError("La secuencia debe tener exactamente 20 nucleótidos.")

    # Codificar la secuencia
    X_new = np.array([one_hot_encode(seq)])

    # Hacer la predicción
    prediccion = xgb_model.predict(X_new)[0] * 100

    return prediccion


def predecir_eficiencia_combined(seq):
    pred_xgb = predecir_eficiencia_xgb(seq)
    pred_nn = predecir_eficiencia_nn(seq)

    # Promedio de ambas predicciones
    pred_final = (pred_xgb + pred_nn) / 2
    return pred_final


# Ejemplo de predicción con una nueva guía CRISPR
# nueva_sec = "AGCTGATCGATGCGTGCTAG"  # Ejemplo de secuencia de 20 nt
# eficiencia_predicha = predecir_eficiencia(nueva_sec)
# print(eficiencia_predicha)
