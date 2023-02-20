import meep as mp
import numpy as np
import ax



def objective(hyperparameters):
    # Define the particle geometry and the simulation parameters
    radius = hyperparameters['radius']
    refractive_index = hyperparameters['refractive_index']
    resolution = hyperparameters['resolution']
    wavelength = hyperparameters['wavelength']
    source_position = hyperparameters['source_position']

    geometry = [mp.Sphere(radius=radius, material=mp.Medium(index=refractive_index))]
    cell_size = mp.Vector3(2 * radius, 2 * radius, 2 * radius)
    pml_layers = [mp.PML(1.0)]
    source = mp.Source(mp.ContinuousSource(wavelength=wavelength),
                       component=mp.Ez,
                       center=mp.Vector3(source_position),
                       size=mp.Vector3(0, 0, 0))

    # Run the simulation and calculate the scattered power
    sim = mp.Simulation(cell_size=cell_size,
                        geometry=geometry,
                        sources=[source],
                        resolution=resolution,
                        boundary_layers=pml_layers)

    sim.run(until=200)

    flux_box = sim.get_flux_box(center=mp.Vector3(), size=mp.Vector3(2 * radius, 2 * radius, 2 * radius))
    scattered_power = np.sum(np.abs(flux_box.get_E() ** 2))

    return {'scattered_power': scattered_power}
	
parameters = [
    {"name": "radius", "type": "range", "bounds": [0.1, 1.0]},
    {"name": "refractive_index", "type": "range", "bounds": [1.0, 3.0]},
    {"name": "resolution", "type": "range", "bounds": [10, 50]},
    {"name": "wavelength", "type": "range", "bounds": [0.4, 0.8]},
    {"name": "source_position", "type": "range", "bounds": [-1.0, 1.0]}
]

experiment = ax.Experiment(
    name="particle_scattering",
    search_space=ax.SearchSpace(parameters),
    objective=ax.Objective(name="scattered_power", minimize=True),
)

for i in range(20):
    # Generate a new set of hyperparameters using the Bayesian optimization algorithm
    parameters, trial_index = experiment.new_trial()

    # Evaluate the objective function with the new hyperparameters
    hyperparameters = parameters.mapping
    evaluation = objective(hyperparameters)

    # Report the evaluation to Ax
    experiment.complete_trial(trial_index=trial_index, raw_data=evaluation)
