# sgRNA_package

## Descripción
**sgRNA_package** es una herramienta para el diseño y la predicción de eficiencia de secuencias sgRNA utilizadas en el sistema CRISPR-Cas9. Emplea modelos de aprendizaje automático como Random Forest, Redes Neuronales y XGBoost para estimar la eficacia de las guías generadas.

## Estructura del Proyecto
```
sgRNA_package/
├── sgRNA/                # Módulo principal del paquete
│   ├── __init__.py
│   ├── percent.py        # Carga y predicción con modelos
│   ├── soporte.py        # Funciones principales de diseño de sgRNAs
│   ├── pdf.py            # Generación de reportes en PDF
│   ├── models/           # Modelos entrenados
│   │   ├── modelo_sgRNA.pkl
│   │   ├── modelo_sgRNA_nn.h5
│   │   ├── modelo_sgRNA_xgb.pkl
├── scripts/
│   ├── realnow.py        # Script de ejemplo para ejecución
├── main.py               # Archivo principal de ejecución
├── setup.py              # Configuración del paquete
├── requirements.txt      # Dependencias del proyecto
├── README.md             # Documentación del paquete
```

## Instalación

Para instalar el paquete, ejecute:
```bash
pip install .
```
Esto instalará todas las dependencias necesarias definidas en `setup.py`.

## Uso
Una vez instalado el paquete, se puede ejecutar desde la terminal con el siguiente comando:
```bash
sgRNA-run genome.fasta target.fasta --modelo nn --output resultados
```

### Parámetros
| Parámetro        | Descripción |
|-----------------|-------------|
| `genome.fasta`  | Archivo FASTA con el genoma de referencia |
| `target.fasta`  | Archivo FASTA con la secuencia objetivo |
| `--modelo`      | Modelo a utilizar: `rf` (Random Forest), `nn` (Red Neuronal), `xgb` (XGBoost). Valor por defecto: `rf` |
| `--output`      | Nombre base del archivo de salida (sin extensión) |

## Resultados
El script generará dos archivos de salida:
- `resultados.csv`: Contiene las secuencias sgRNA diseñadas y su eficiencia estimada.
- `resultados.pdf`: Informe en formato PDF con la información de las guías generadas.

## Dependencias
El paquete requiere las siguientes bibliotecas:
- `biopython`
- `numpy`
- `pandas`
- `tensorflow`
- `joblib`
- `reportlab`
- `xgboost`
- `scikit-learn`

Para instalar manualmente las dependencias:
```bash
pip install -r requirements.txt
```

Para probar los cambios localmente:
```bash
python main.py genome.fasta target.fasta --modelo xgb --output test
```

## Licencia
Este proyecto está bajo la licencia **MIT**. Ver el archivo `LICENSE` para más detalles.

