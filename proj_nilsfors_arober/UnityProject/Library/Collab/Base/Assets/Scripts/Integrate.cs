using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class Integrate
{
    public static float timeStep = 0.02f;
    public static void forwardEuler(Particle p, float dt)
    {
        p.velocity += dt * p.force / p.density;
        p.position += dt * p.velocity;
    }

    public static void LeapFrogPos(Particle p, float dt)
    {
        p.position += p.velocity * timeStep + 0.5f * p.force * dt * dt;
        
    }
    public static void LeapFrogVel(Particle p, float dt)
    {
        p.velocity += 0.5f * (p.oldForce + p.force) * dt / p.density;
    }
}
