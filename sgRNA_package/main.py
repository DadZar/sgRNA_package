import argparse
import pandas as pd
from sgRNA.paq1_soporte import design_sgRNAs_rf, design_sgRNAs_nn, design_sgRNAs_xgb
from sgRNA.crr_pdf import create_pdf

def main():
    parser = argparse.ArgumentParser(description='Ejecuta el diseño de sgRNAs para CRISPR-Cas9.')
    parser.add_argument('gen_file', type=str, help='Archivo FASTA del genoma de referencia')
    parser.add_argument('target_file', type=str, help='Archivo FASTA con la secuencia objetivo')
    parser.add_argument('--modelo', type=str, choices=['rf', 'nn', 'xgb'], default='rf', help='Modelo a utilizar: Random Forest (rf), Red Neuronal (nn) o XGBoost (xgb)')
    parser.add_argument('--output', type=str, default='sgRNA_results', help='Nombre base del archivo de salida (sin extensión)')
    
    args = parser.parse_args()
    
    print("Diseñando sgRNAs...")
    
    if args.modelo == 'rf':
        df_sgRNA = design_sgRNAs_rf(args.gen_file, args.target_file)
    elif args.modelo == 'nn':
        df_sgRNA = design_sgRNAs_nn(args.gen_file, args.target_file)
    else:
        df_sgRNA = design_sgRNAs_xgb(args.gen_file, args.target_file)
    
    csv_output = f"{args.output}.csv"
    pdf_output = f"{args.output}"
    
    df_sgRNA.to_csv(csv_output, index=False)
    create_pdf(df_sgRNA, pdf_output)
    
    print(f"Resultados guardados en {csv_output} y {pdf_output}.pdf")

if __name__ == "__main__":
    main()
