# Millimeter Wave Blockage Modeling

Work towards my thesis at NYU. This work is published as conference paper at international teletraffic congress (ITC30) and the journal paper at JSAC special issue on URLLC applications.

We present a simplified model for the key QoS parameters such as blockage probability, frequency, and duration in mmWave cellular systems. Our model considers an open park-like area with dynamic blockage due to mobile blockers and self-blockage due to user's own body. A typical user is at the center and the BSs are distributed uniformly around the user. The user is considered blocked when all potential BSs around the UE are blocked simultaneously. 

## Understand the Code

#### Directory: [Simulations](Simulations)

The MATLAB simulations considers random waypoint mobility model for blockers. The main file is [SimulationLOS.m](Simulations/SimulationLOS.m).

It takes around 1 hour to run for a single iteration. We obtained results for 10,000 iterations using [NYU HPC](HPC_Matlab.md) (high performance computing)

Note: SimulationLOS.m uses the following the functions- [Generate_Mobility.m](Simulations/Generate_Mobility.m) (copied from [here)](https://www.mathworks.com/matlabcentral/fileexchange/30939-random-waypoint-mobility-model), [BlockageSimFn.m](Simulations/BlockageSimFn.m), [find_blockage_distance.m](Simulations/find_blockage_distance.m).



#### Directory: [Theory](Theory)

The [TheoryLOS.m](Theory/TheoryLOS.m) implements the theoretical results of LOS blockage and [TheoryNLOS.m](Theory/TheoryNLOS.m) for NLOS blockage.

#### Directory: [HexagonalCase](HexagonalCase)
[hexagonal.m](HexagonalCase/hexagonal.m) considers hexagonal cell deployment of BSs for open park scenario. 

#### Directory: [DataProcessing](DataProcessing)
The data obtained from NYU HPC by running SimulationLOS.m is analysed using [processData9.m](DataProcessing/processData9.m) 
Finally, [plotResults.m](DataProcessing/plotResults.m) takes the data from csv files and plots the nice figures comparing theory and simulation :)

#### Directory: [CaseStudy](CaseStudy)
It provides several code and plots for better understanding the paper. But the codes are not well maintained in this directory.

## Results
Our results tentatively show that the density of BS required to provide acceptable quality of experience to AR/VR applications is much higher than that obtained by capacity requirements alone. This suggests that the mmWave cellular networks may be blockage limited instead of capacity limited. 

## Slides and papers
Please look at the ITC slides [here](ITC_slides.pdf)

Jain, Ish Kumar, Rajeev Kumar, and Shiendra Panwar. "[Driven by capacity or blockage? a millimeter wave blockage analysis](https://ieeexplore.ieee.org/abstract/document/8493070)." 2018 30th International Teletraffic Congress (ITC 30). Vol. 1. IEEE, 2018.

Jain, Ish Kumar, Rajeev Kumar, and Shivendra Panwar. "[Can Millimeter Wave Cellular Systems provide High Reliability and Low Latency? An analysis of the impact of Mobile Blockers.](https://arxiv.org/pdf/1807.04388.pdf)" arXiv preprint arXiv:1807.04388 (2018).


## Future Work

The following extensions are planned for future work

* **Generalized blockage model:** We can improve the theoretical analysis of static, dynamic, and self-blockage. 
* **Capacity analysis:** The data rates of a typical user can be evaluated using the generalized blockage model. We are interested in evaluating whether 5G mmWave is capacity limited or blockage limited.
* **Fallback to 4G LTE:** We plan to explore the potential solution to blockages as switching to 4G LTE. Whether 4G would be able to handle the huge intermittent 5G traffic.
* **Deterministic networks:** We have considered a random deployment of BSs in our analysis. However, in most cases, the deployments of BSs are based on a deterministic hexagonal grid. Therefore, a blockage model for the deterministic networks is more practical.
* **Backhoul capacity analysis:** With the UEs switching between BSs in case of blockage events, the backhoul capacity requirements for the BS may have high fluctuations. It is interesting to study that random process.
* **Correlated blockage:** We assumed the dynamic blockages of all the BSs are independent (So the overall blockage probability is the product of the blockage probabilities of all the BS-UE links). However, when blockage of multiple BSs are correlated (depends on blocker location;) we might get a higher blockage probability.

Please cite [our paper](https://arxiv.org/pdf/1807.04388.pdf) if you are using any part of this code.
