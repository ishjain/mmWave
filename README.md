# Millimeter Wave Blockage Modeling

Work towards my thesis at NYU

We presented a simplified model for the key QoS parameters such as blockage probability, frequency, and duration in mmWave cellular systems. Our model considered an open park-like area with dynamic blockage due to mobile blockers and self-blockage due to user's own body. A typical user is at the center and the BSs are distributed uniformly around the user. The user is considered blocked when all potential BSs around the UE are blocked simultaneously. 

# Simulations
The MATLAB simulations considers random waypoint mobility model for blockers. The main file is Simulation.m
It takes around 1 hour to run for a single iteration. We obtained results for 10,000 iterations using [NYU HPC](HPC_Matlab.md) (high performance computing)

# Results
Our results tentatively show that the density of BS required to provide acceptable quality of experience to AR/VR applications is much higher than that obtained by capacity requirements alone. This suggests that the mmWave cellular networks may be blockage limited instead of capacity limited. 

Also, the optimal height of BSs would be lower as compared to microwave BSs. 

# Future Work

The following extensions are planned for future work

* **Generalized blockage model:** We can add a simple model for static blockage in our analysis of dynamic and self-blockage. 
\item Data rate analysis: The data rates of a typical user can be evaluated using the generalized blockage model. We are interested in evaluating whether 5G mmWave is capacity limited or blockage limited.
* **Fallback to 4G LTE:** We plan to explore the potential solution to blockages as switching to 4G LTE. Whether 4G would be able to handle the huge intermittent 5G traffic.
* **Deterministic networks:** We have considered a random deployment of BSs in our analysis. However, in most cases, the deployments of BSs are based on a deterministic hexagonal grid. Therefore, a blockage model for the deterministic networks is more practical.
* **Backhoul capacity analysis:** With the UEs switching between BSs in case of blockage events, the backhoul capacity requirements for the BS may have high fluctuations. It is interesting to study that random process.


