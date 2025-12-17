from quantum_simulation.experiments.gallery.tunneling_barrier import TunnelingBarrier
from quantum_simulation.experiments.gallery.double_slit_2d import DoubleSlit2D
from quantum_simulation.experiments.gallery.wave_packet_evolution import WavePacketEvolution
from quantum_simulation.experiments.gallery.measurement_statistics import MeasurementStatistics
from quantum_simulation.orchestration.pipeline import ExperimentPipeline

def test_pipeline_sequential():
    """Pipeline 2 expériences séquentielles."""
    exp1 = WavePacketEvolution(config1)
    exp2 = MeasurementStatistics(config2)
    
    pipeline = ExperimentPipeline([exp1, exp2])
    results = pipeline.run(parallel=False)
    
    assert results.n_experiments == 2
    assert results.all_passed
    
def test_pipeline_checkpoint():
    """Reprise calcul depuis checkpoint."""
    pipeline = ExperimentPipeline(experiments)
    pipeline.checkpoint('checkpoint.pkl')
    
    # Simuler crash
    pipeline2 = ExperimentPipeline([])
    pipeline2.load_checkpoint('checkpoint.pkl')
    
    results = pipeline2.run()
    assert results.n_experiments == len(experiments)