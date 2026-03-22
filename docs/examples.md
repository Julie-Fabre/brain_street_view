# Example usage

Brain Street View has been used in published neuroscience research. Below are examples showing how BSV was applied in each study.

---

## Mapping dopaminergic projections to the striatum

Pan-Vazquez et al. (2025) used BSV to map and visualize projection patterns from the ventral tegmental area (VTA) to the striatum, characterizing the spatial organization of dopaminergic inputs to the caudate putamen.

> Pan-Vazquez, A., Zimmerman, C. A., McMannon, B., Fabre, J. M. J., et al. *VTA dopamine neuron activity produces spatially organized value representations.* bioRxiv (2025). [DOI: 10.1101/2025.11.04.685995](https://doi.org/10.1101/2025.11.04.685995)

```python
import bsv

# Find VTA injection experiments
experiment_ids = bsv.find_connectivity_experiments(['VTA'])

# Fetch and normalize by injection intensity
experiment_imgs, injection_summary, _, _ = bsv.fetch_connectivity_data(
    experiment_ids, save_location, '',
    'injectionIntensity', False,
    allen_atlas_path=allen_atlas_path)

# Visualize projections to caudate putamen
bsv.plot_connectivity(experiment_imgs, allen_atlas_path, 'CP',
                       10, 15, 'coronal', True, 2, 'global', None,
                       'injectionIntensity')

# 3D rendering
bsv.plot_connectivity_3d(injection_summary, allen_atlas_path, 'CP',
                          plot_patch=True)
```

---

## Prefrontal-striatal projections in sensorimotor learning

Song & Peters (2025) used BSV to examine how prefrontal cortex projects to the striatum, comparing projection patterns across cortical regions involved in cross-modal sensorimotor learning.

> Song, D. & Peters, A. J. *Prefrontal cortex and striatum dissociate modality and learning order in cross-modal sensorimotor learning.* bioRxiv (2025). [DOI: 10.1101/2025.11.10.687671](https://doi.org/10.1101/2025.11.10.687671)

```python
import bsv

# Find prefrontal cortex injection experiments
experiment_ids = bsv.find_connectivity_experiments(['PL', 'ILA', 'ACAd', 'ACAv'])

experiment_imgs, injection_summary, _, _ = bsv.fetch_connectivity_data(
    experiment_ids, save_location, '',
    'injectionIntensity', False,
    allen_atlas_path=allen_atlas_path)

# Compare projections to striatum across prefrontal regions
bsv.plot_connectivity(experiment_imgs, allen_atlas_path, 'CP',
                       10, 15, 'coronal', True, 2, 'global', None,
                       'injectionIntensity')
```

---

## Amygdala connectivity in risk-based decision making

Piantadosi et al. (2025) used BSV to visualize amygdala projection patterns to understand the anatomical connectivity underlying risk-sensitive choice behavior.

> Piantadosi, P. T., Coden, K. M., Choi, H., et al. *Risk reshapes amygdala representation of choice.* bioRxiv (2025). [DOI: 10.1101/2025.10.06.680813](https://doi.org/10.1101/2025.10.06.680813)

```python
import bsv

# Find basolateral amygdala injection experiments
experiment_ids = bsv.find_connectivity_experiments(['BLA'])

experiment_imgs, injection_summary, _, _ = bsv.fetch_connectivity_data(
    experiment_ids, save_location, '',
    'injectionIntensity', False,
    allen_atlas_path=allen_atlas_path)

# Visualize amygdala projections
bsv.plot_connectivity(experiment_imgs, allen_atlas_path, 'CP',
                       10, 15, 'coronal', True, 2, 'global', None,
                       'injectionIntensity')
```

---

If you use Brain Street View in your research, please cite it using the [CITATION.cff](https://github.com/Julie-Fabre/brain_street_view/blob/main/CITATION.cff) file. To add your paper here, open a [GitHub issue](https://github.com/Julie-Fabre/brain_street_view/issues) or pull request.
