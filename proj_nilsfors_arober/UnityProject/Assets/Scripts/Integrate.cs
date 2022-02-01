using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class Integrate
{
    public static float timeStep = 0.02f;
    public static float slowDownFactor = 1f;
    public static void forwardEuler(Particle p, float dt)
    {
        p.velocity += dt * p.force / p.density;
        p.position += dt * p.velocity;
    }

    public static void LeapFrogPos(Particle p, float dt)
    {
        p.position += slowDownFactor * (p.velocity * dt + 0.5f * p.force * dt * dt /p.density);
        
    }
    public static void LeapFrogVel(Particle p, float dt)
    {
        p.velocity += slowDownFactor * (0.5f * (p.oldForce/p.oldDensity + p.force/p.density) * dt);
    }
}
