

class ReportGenerator:
    """
    Génération rapports synthétiques multi-expériences.
    
    Formats:
        - HTML interactif (plotly.js)
        - PDF statique (matplotlib + reportlab)
        - JSON structuré (archivage)
    """
    
    def generate_html_report(self, pipeline_results: PipelineResults,
                            output_path: str = "report.html"):
        """
        Génère rapport HTML avec :
        - Tableau récapitulatif validations
        - Figures interactives (zoom, hover)
        - Logs exécution (temps, warnings)
        """
        
    def generate_pdf_report(self, pipeline_results: PipelineResults,
                            output_path: str = "report.pdf"):
        """
        PDF publication-ready avec :
        - Résumé exécutif
        - Figures haute résolution (300 dpi)
        - Annexes (paramètres, équations)
        """