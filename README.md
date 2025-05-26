# Drone and Vehicle Based Quantum Key Distribution

This is a public version of the code used in *Drone and Vehicle Based Quantum Key Distribution* \[[Arxiv](https://arxiv.org/abs/2505.17587)\]. This was built for [v2.0.2](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/releases/tag/v2.0.2) of the Open QKD Security package.

For each operational mode, e.g. "Air to Air", run the main file in the corresponding folder. To switch between different experimental runs of one operational mode follow the instructions summarized in the table below.


| Operation mode    | Folder           | Selecting different experimental runs                                       |
| ----------------- | ---------------- | --------------------------------------------------------------------------- |
| Air to Air        | `AirtoAir`       | Toggle between different presets in `ThreeStateDataAirDifInt_Main.m`        |
| Air to Car        | `AirtoCar`       | Toggle between different presets in `ThreeStateDataAirtoCarDifInt_Main.m`   |
| Car to Car        | `CartoCar`       | Toggle between different data files in `ThreeStateDataCarDifInt_Preset.m`   |
| Ground to Ground  | `GroundtoGround` | Only one experimental run performed                                         |
| Table Top         | `TableTop`       | Toggle between different data files in `ThreeStateDataTableDifInt_Preset.m` |



## Installation instructions
> [!CAUTION]
> This repository is for archival and transparency purposes; we do not guarantee compatibility with other versions of the Open QKD Security package beyond the ones listed above.

### As zip
1. Download the linked version of the code from above and follow all [installation instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/v2.0.2).
2. Also follow the additional Mosek install instructions if you want an exact match.
3. Download the latest release on the side bar and unzip in your preferred directory and add this folder to the Matlab path.


### with git
1. Clone this repository and its exact submodules navigate to your desired directory and run,
```
git clone --recurse-submodules https://github.com/Optical-Quantum-Communication-Theory/Drone-and-Vehicile-based-QKD
```
2. Follow all further [installation instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/v2.0.2).
3. Also follow the additional Mosek install instructions if you want an exact match.
