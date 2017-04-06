# Resonant-Acoustic-Mixing
I don't consider myself as the most elegant programmer, but I know enough to get around. The following GIT repository contains some critical coding files in development toward my masters study in dimensionalising the Resonant Acoustic Mixing problem.

# Dimensional Analysis
This file creates a DIMENSIONAL sweep of the design parameters in question the following list of parameters is swept to find the effect of each of the parameters:

The parameters take effect on a sinusoidal movement that is defined by:

y(t) = A sin(w t)

1. Global parameters:
    * Cup Height: H
    * Cup Diameter: D
    * Amplitude: A
    * Frequncy: f
    * Run Time: T
    * Gravity constant: g
    
2. Local parameters, for $n$ particle TYPES:

    1. Particle Properties
        * Particle Diameter: d (size: [n x 1])
        * Particle Number: N (size: [n x 1])
        * Particle Density: p (size: [n x 1])
        * Young's Modulus: E (size: [n x 1])
        * Poison's Ratio: v (size: [n x 1])
    2. Inter Particle Properties
        * Restitution Coefficient [normal]: R (size: [n x n], symmetrical)
        * Restitution Coefficient [tangential]: R' (size: [n x n], symmetrical)
        * Wall Restitution Coefficient [normal]: W (size: [n x n], symmetrical)
        * Wall Restitution Coefficient [tangential]: W' (size: [n x n], symmetrical)
        * Static friction Coefficient [normal]: u (size: [n x n], symmetrical)

3. Producing dependent variables:
    * Mean collision time: dt
    * Mean free path: y
    * Particle velocity: vi

By using these parameters we can write a java macro script as a input file to StarCCM+. The master.sim file will be used as a base file.

Addline