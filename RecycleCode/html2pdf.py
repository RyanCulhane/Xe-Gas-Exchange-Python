import pdfkit
import os

cwd = os.getcwd()

options = {
    'page-width':140,
    'page-height':100,
    'margin-top': 1,
    'margin-right': 0.1,
    'margin-bottom': 0.1,
    'margin-left': 0.1,
    'encoding': "UTF-8",
    # 'dpi':400
    }

pdfkit.from_file('clinical_report.html', 'clinical_report.pdf',options=options)
