"""
Exécute catalogue complet expériences prédéfinies.

Usage:
    python run_gallery.py --all               # Toutes expériences
    python run_gallery.py --tag 2d            # Seulement 2D
    python run_gallery.py --experiment double_slit
"""

import argparse
from quantum_simulation.experiments.gallery import EXPERIMENT_REGISTRY
from quantum_simulation.orchestration.pipeline import ExperimentPipeline
from quantum_simulation.orchestration.reports import ReportGenerator

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--all', action='store_true')
    parser.add_argument('--tag', type=str, help='2d, 3d, tunneling, etc.')
    parser.add_argument('--experiment', type=str)
    parser.add_argument('--parallel', action='store_true')
    parser.add_argument('--output', default='./gallery_results/')
    
    args = parser.parse_args()
    
    # Sélection expériences
    if args.all:
        experiments = EXPERIMENT_REGISTRY.get_all()
    elif args.tag:
        experiments = EXPERIMENT_REGISTRY.get_by_tag(args.tag)
    elif args.experiment:
        experiments = [EXPERIMENT_REGISTRY.get(args.experiment)]
    else:
        print("Erreur : spécifier --all, --tag ou --experiment")
        return
        
    # Exécution pipeline
    pipeline = ExperimentPipeline(experiments)
    results = pipeline.run(parallel=args.parallel, n_workers=4)
    
    # Rapport
    reporter = ReportGenerator()
    reporter.generate_html_report(results, f"{args.output}/gallery_report.html")
    
    print(f"✓ {len(experiments)} expériences terminées")
    print(f"  Rapport : {args.output}/gallery_report.html")

if __name__ == "__main__":
    main()