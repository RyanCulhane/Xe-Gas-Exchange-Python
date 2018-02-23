import pdfkit
import os

cwd = os.getcwd()

options = {
    'page-width':160,
    'page-height':80,
    'margin-top': 1,
    'margin-right': 0.1,
    'margin-bottom': 0.1,
    'margin-left': 0.1,
    'encoding': "UTF-8",
    # 'dpi':400
    }

pdfkit.from_file('report_clinical.html', 'report_clinical.pdf',options=options)
