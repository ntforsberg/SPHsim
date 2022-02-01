using System.Collections;
using System.Collections.Generic;
using System.Threading.Tasks;
using UnityEngine;
using UnityEngine.Profiling;

public class Fluid : MonoBehaviour
{
    [SerializeField]
    private GameObject particlePrefab = null;

    [SerializeField]
    private float density = 997; // density at rest, measured in kg/m^3

    [SerializeField]
    private float volume = 1; // total volume of fluid, measured in m^3

    //[SerializeField]
    public float particleMass = 1;

    [SerializeField]
    private int count = 125; // number of particles in simulation

    [SerializeField]
    private float h = 3f; // neighbor radius

    private float hSqr;

    [SerializeField]
    private int neighborCount = 32;

    [SerializeField]
    private float k = 50.0f; // gas stiffness constant, k = nRT

    [SerializeField]
    private float viscosityCoefficient = 0.8f; //approx viscosity of water at 30deg C

    [SerializeField]
    private float restitutionCoefficient = 0.0f;

    [SerializeField]
    private float tensionCoefficient = 0.0728f;

    [SerializeField]
    private float tensionThreshold = 7.0f;

    [SerializeField]
    public float viscosityForceMultiplier = 1;

    [SerializeField]
    public float pressureForceMultiplier = 1;

    [SerializeField]
    public float externalForceMultiplier = 1;

    [SerializeField]
    public float particleRadius = 1;

    [SerializeField]
    public float separationFactor = 1.4f;

    [SerializeField]
    private bool precompute = true;

    [SerializeField]
    private float time = 15.0f; // in seconds


    private float particleScale;

    private bool oddStep = true;

    private bool onlyOnce = true;

    private List<Particle> particles;
    private SpatialHash spatialHash;

    CustomSampler setNeighbors;
    CustomSampler updateSpatialHash;
    CustomSampler computeDensity;
    CustomSampler computePressure;
    CustomSampler computeForces;
    CustomSampler computeInternalForces;
    CustomSampler computeExternalForces;
    CustomSampler integration;

    // Start is called before the first frame update
    void Start()
    {
        setNeighbors = CustomSampler.Create("SetNeighbors");
        updateSpatialHash = CustomSampler.Create("UpdateSpatialHash");
        computeDensity = CustomSampler.Create("ComputeDensity");
        computePressure = CustomSampler.Create("ComputePressure");
        computeForces = CustomSampler.Create("ComputeForces");
        computeInternalForces = CustomSampler.Create("ComputeInternalForces");
        computeExternalForces = CustomSampler.Create("ComputeExternalForces");
        integration = CustomSampler.Create("Integration");

        initializeVariables();

        

        particles = GenerateParticlesPerturbedCube(volume, count);

        
        //particles = new List<Particle>();
        //particles.Add(AddSingleParticle());



        if (neighborCount != 0)
        {
            h = Mathf.Pow(3*volume*neighborCount / (4 * Mathf.PI * count), 1f/3f);
            //h = 2 * (separationFactor * particleRadius);
        }

        hSqr = h * h;

        spatialHash = new SpatialHash(count, Mathf.CeilToInt(h), h);
        spatialHash.Insert(particles);

        // Set the particle mass
        particleMass = density * (volume / count);
        //density = mass * count / volume;
        Debug.Log("mass: " + particleMass);
        foreach (Particle p in particles)
        {
            p.mass = particleMass;
            // NOTE This is only a temporary thing to test variable scale of particles, remove later
            p.density = density;
            p.particleScale = particleScale;
        }

        // NOTE Debug first particle
        Particle d = particles[0];
        d.debug = true;
        Debug.Log("START");
        Debug.Log("-----");

        Debug.Log("position: " + d.position);
        Debug.Log("velocity: " + d.velocity);

        Debug.Log("force: " + d.force);

        Debug.Log("mass: " + d.mass);
        Debug.Log("density: " + d.density);
        Debug.Log("pressure: " + d.pressure);

        Debug.Log("-----");


        // Precompute positions
        if (precompute)
        {
            Physics.autoSimulation = false;

            int steps = (int) (time / Time.fixedDeltaTime);
            Debug.Log("time: " + time);
            Debug.Log("fixedDeltaTime: " + Time.fixedDeltaTime);
            Debug.Log("steps: " + steps);

            // Setup position storage for particles
            for (int i = 0; i < particles.Count; i++)
            {
                particles[i].positions = new List<Vector3>(steps);
            }

            // Setup ParticleCollision
            for (int i = 0; i < particles.Count; i++)
            {
                particles[i].gameObject.SendMessage("Start");
            }

            for (int i = 0; i < steps; i++)
            {
                Physics.Simulate(Time.fixedDeltaTime);
                Simulate();
                // Store position
                for (int j = 0; j < particles.Count; j++)
                {
                    particles[j].positions.Add(particles[j].position);
                }
            }
        }
    }

    List<Particle> GenerateParticlesCube(float volume, int count)
    {
        Debug.Log("volume; " + volume);
        Debug.Log("count: " + count);

        float sideLength = Mathf.Pow(volume, 1f / 3f);
        int sideCount = (int) Mathf.Pow(count, 1f / 3f);
        int extra = count - (int) Mathf.Pow(sideCount, 3);

        Debug.Log("sideLength: " + sideLength);
        Debug.Log("sideCount: " + sideCount);
        Debug.Log("extra: " + extra);

        float delta = (float) sideLength / (float) sideCount;

        Debug.Log("delta: " + delta);

        List<Particle> l = new List<Particle>();
        Vector3 origin = new Vector3(-sideLength/2, 1.0f, -sideLength / 2);

        for (int i = 0; i < sideCount; i++)
        {
            for (int j = 0; j < sideCount; j++)
            {
                for (int k = 0; k < sideCount; k++)
                {
                    Vector3 pos = origin + Vector3.right * i * delta + Vector3.forward * j * delta + Vector3.up * k * delta;
                    GameObject o = Instantiate(particlePrefab, pos, particlePrefab.transform.rotation, this.transform);
                    Particle p = o.GetComponent<Particle>();
                    p.position = pos;
                    p.velocity = Vector3.zero;
                    l.Add(p);
                }
            }
        }

        // TODO Add the extra particles
        return l;
    }

    Particle AddSingleParticle()
    {
        Vector3 pos = new Vector3(Random.Range(-6, 6), 3f, Random.Range(-6, 6));
        GameObject o = Instantiate(particlePrefab, pos, particlePrefab.transform.rotation, this.transform);
        Particle p = o.GetComponent<Particle>();
        p.position = pos;
        p.velocity = Vector3.zero;
        return p;
    }

    List<Particle> GenerateParticlesPerturbedCube(float volume, int count)
    {
        List<Particle> particles = GenerateParticlesCube(volume, count);
        foreach (Particle p in particles)
        {
            p.position += Random.insideUnitSphere * 0.05f;
        }
        return particles;
    }

    void FixedUpdate()
    {
        if (!precompute)
            Simulate();
    }

    void Simulate()
    {
        //Debug.Log("FIXEDUPDATE");
        /*
        if(particles.Count < count)
        {
            particles.Add(AddSingleParticle());
        }
        */

        updateSpatialHash.Begin();
        spatialHash.Update();
        updateSpatialHash.End();

        setNeighbors.Begin();
        spatialHash.SetNeighbors();
        setNeighbors.End();

        //For all particles, compute density -> pressure -> forces
        Parallel.ForEach(particles, particle =>
        {
            ComputeDensity(particle);
            ComputePressure(particle);
        });

        Parallel.ForEach(particles, particle =>
        {
            ComputeForces(particle);
        });

        integration.Begin();
        float dt = Time.fixedDeltaTime;
        //Debug.Log("FDT:" + Time.fixedDeltaTime);
        if (oddStep)
        {
            foreach (Particle particle in particles)
            {
                Integrate.LeapFrogVel(particle, dt);
                //particle.positionUpdate();
            }
            oddStep = false;
        }
        else
        {
            foreach (Particle particle in particles)
            {
                Integrate.LeapFrogPos(particle, dt);
                particle.positionUpdate();
            }
            oddStep = true;
        }


        integration.End();
        
        Parallel.ForEach(particles, particle =>
        {
            velocityCorrection(particle);
        });

        
        //Note - need to first update density for all particles, then pressure for all particles etc.
        //So neighbors need to be saved away somehow to be reused in different fnk calls
        //Maybe a list of every particle of its neighbors? Or is that too inefficient memory wise?
    }

    void ComputeForces(Particle p)
    {
        computeForces.Begin();
        p.oldForce = p.force;
        p.force = Vector3.zero;
        ComputeInternalForces(p);
        ComputeExternalForces(p);
        computeForces.End();
    }

    void ComputeDensity(Particle i)
    {
        computeDensity.Begin();
        //float dens = i.mass;
        i.oldDensity = i.density;
        i.density = i.mass; //density appears in denominators later on, if i has no neighbors this will get bad
        foreach (Particle j in i.neighbors)
        {
            i.density += j.mass * Kernel.Poly6(i.position - j.position, h); 
        }
        //i.density = Mathf.Max(i.density, density);
        /*
        if (i.debug)
        {
            Debug.Log("density: " + i.density);
        }
        */
        computeDensity.End();
    }

    void ComputePressure(Particle i)
    {
        computePressure.Begin();
        float p = k * (i.density - density);
        i.pressure = p;
        /*
        if (i.debug)
            Debug.Log("pressure: " + p);
        */
        computePressure.End();
    }

    void ComputeInternalForces(Particle i)
    {
        computeInternalForces.Begin();

        ComputePressureForce(i);
        ComputeViscosityForce(i);

        computeInternalForces.End();
    }

    void ComputePressureForce(Particle i)
    {
        Vector3 pressureGrad = Vector3.zero;

        foreach (Particle j in i.neighbors)
        {
//            pressureGrad += j.mass * (i.pressure / Mathf.Pow(i.density, 2) + (j.pressure / Mathf.Pow(j.density, 2))) * Kernel.SpkiyGrad(i.position - j.position, h);
            pressureGrad += j.mass * (i.pressure + j.pressure) / (2 * j.density) * Kernel.Poly6Grad(i.position - j.position, h);
        }

        //Assign directly to force here?
        //Debug.Log("pressureGrad:" + pressureGrad);
//        pressureGrad *= i.mass;
        i.force -= pressureGrad * pressureForceMultiplier;
    }

    void ComputeViscosityForce(Particle i)
    {
        Vector3 viscosityForce = Vector3.zero;

        foreach (Particle j in i.neighbors)
        {
            viscosityForce += j.mass * (j.velocity - i.velocity) / j.density * Kernel.ViscosityLap(i.position - j.position, h);
        }

        i.force += viscosityCoefficient * viscosityForce * viscosityForceMultiplier;
    }

    void ComputeExternalForces(Particle i)
    {
        computeExternalForces.Begin();
        

        // TODO Calculate other external forces

        //ComputeCollision(i);
        //ComputeCollisionFlat(i);
        ComputeCollisionKruger(i);

        ComputeGravity(i);
        ComputeSurfaceTension(i);

        computeExternalForces.End();
    }

    void ComputeCollision(Particle i)
    {
        float k = 1000.0f;
        if (i.collision)
        {
            Vector3 f = Vector3.zero;
            f = -k * i.collisionSeparation * i.collisionNormal;
            i.force += f;
        }
    }

    void ComputeCollisionFlat(Particle i)
    {
        if (i.collision)
        {
            i.position -= i.collisionNormal * i.collisionSeparation;
            i.velocity = Vector3.ProjectOnPlane(i.collisionNormal, i.velocity);
            //i.force = Vector3.ProjectOnPlane(cp.normal, i.force);
        }
    }

    void ComputeCollision2point0(Particle p)
    {
        if (p.collision)
        {
            if(p.collisionSeparation < 0)
            {
                p.position = p.collisionPoint;
                p.velocity = Vector3.ProjectOnPlane(p.collisionNormal, p.velocity);
            }
            
            
        }
    }

    void ComputeCollisionKruger(Particle i)
    {
        if (i.collision)
        {
            if (i.collisionSeparation < 0)
            {
                //i.position = i.collisionPoint;
                i.position -= i.collisionSeparation * i.collisionNormal;
                i.velocity -= (1 + restitutionCoefficient) * Vector3.Dot(i.velocity, i.collisionNormal) * i.collisionNormal;
            }
        }
    }

    void ComputeSurfaceTension(Particle i)
    {
        float colorField = 0.0f;
        Vector3 colorGradient = Vector3.zero;
        float colorLaplacian = 0.0f;
        foreach (Particle j in i.neighbors)
        {
            colorField += j.mass / j.density * Kernel.Poly6(i.position - j.position, h);
        }

        foreach (Particle j in i.neighbors)
        {
            colorGradient += j.mass / j.density * Kernel.Poly6Grad(i.position - j.position, h);
        }

        foreach (Particle j in i.neighbors)
        {
            colorLaplacian += j.mass / j.density * Kernel.Poly6Lap(i.position - j.position, h);
        }

        if (colorGradient.sqrMagnitude < tensionThreshold * tensionThreshold)
        {
            // This particle isn't close to the surface
            return;
        }

        i.force += -tensionCoefficient * colorLaplacian * colorGradient / colorGradient.magnitude;
    }

    void ComputeGravity(Particle i)
    {
       i.force += i.mass * Physics.gravity * externalForceMultiplier; 
    }

    void initializeVariables()
    {
        particleScale = 2 * particleRadius;
        volume = Mathf.Pow(particleScale, 3) * count * separationFactor;
        //density = particleMass * count / volume;
        particleMass = density * volume / count;
        Physics.gravity = new Vector3(0f, -9.88f, 0f);
    }

    void velocityCorrection(Particle particle)
    {
        Vector3 correctionValue = Vector3.zero;
        foreach(Particle p in particle.neighbors)
        {
            correctionValue += (p.velocity - particle.velocity) * p.mass * Kernel.Poly6(particle.position - p.position, h) / (p.density + particle.density);
        }
        particle.velocity += correctionValue;
    }

    /*
    void ComputeSurfaceTension(Particle p)
    {
        float surfaceForce = 0.0f;
        foreach(Particle p in p.neighbors)
        {
        }
    }

    void ComputeSurfaceNormal()
    {

    }
    */

}
