from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import pandas as pd

def create_pdf(df_sgRNA, name_file_pdf):
    """
    Genera un PDF con información de guías sgRNA usando reportlab.
    
    Parámetros:
        df_sgRNA (pd.DataFrame): DataFrame con las guías encontradas.
        name_file_pdf (str): Nombre del archivo de salida (sin extensión).
    """
    # Definir tamaño del PDF
    pdf_file = f"{name_file_pdf}.pdf"
    c = canvas.Canvas(pdf_file, pagesize=letter)
    
    # Título alineado a la izquierda
    c.setFont("Helvetica-Bold", 16)
    c.drawString(40, 750, "Guías de sgRNA")
    
    # Espacio antes de la sección
    c.setFont("Helvetica", 12)
    c.drawString(40, 720, f"Se encontraron {len(df_sgRNA)} guías para la secuencia objetivo.")
    c.drawString(40, 705, "Solo se muestran las guías cuyo GC% tiene valores entre 40%-80%.")
    
    # Posición inicial de la tabla
    x_start = 40
    y_start = 680
    row_height = 20

    # Definir anchos de columna ajustados
    col_widths = [160, 50, 50, 60, 50, 70]  # Más espacio para "gRNA"

    # Dibujar encabezado de la tabla
    headers = ["gRNA", "PAM", "GC%", "Posición", "Hebra", "Eficiencia"]
    c.setFont("Helvetica-Bold", 10)
    for col, (header, width) in enumerate(zip(headers, col_widths)):
        c.drawString(x_start + sum(col_widths[:col]), y_start, header)

    # Dibujar filas de la tabla
    c.setFont("Helvetica", 9)
    for i, row in df_sgRNA.iterrows():
        y_pos = y_start - (i + 1) * row_height
        if y_pos < 50:  # Cambiar de página si se acaba el espacio
            c.showPage()
            y_pos = 750  # Reiniciar la posición en la nueva página

        # Dibujar cada dato en la columna correspondiente
        c.drawString(x_start, y_pos, row["gRNA"])  # gRNA más ancho
        c.drawString(x_start + col_widths[0], y_pos, row["PAM"])
        c.drawString(x_start + sum(col_widths[:2]), y_pos, f"{row['GC_content']}%")
        c.drawString(x_start + sum(col_widths[:3]), y_pos, str(row["position"]))
        c.drawString(x_start + sum(col_widths[:4]), y_pos, row["hebra"])
        c.drawString(x_start + sum(col_widths[:5]), y_pos, f"{row['Eficiencia']:.2f}")

    # Guardar PDF
    c.save()
    print(f"PDF generado: {pdf_file}")
