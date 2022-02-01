using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Particle : MonoBehaviour
{
    public List<Vector3> positions = null;

    public float mass = 1.0f;
    public float density = 1f; //about 0.998
    public float pressure = 1.0f;

    public bool debug = false;

    //public Collision collision = null;
    public bool collision = false;
    public Vector3 collisionNormal = Vector3.zero;
    public Vector3 collisionPoint = Vector3.zero;
    public float collisionSeparation = 0.0f;

    public List<Particle> neighbors;

    public Vector3 position;
    public Vector3 velocity;

    public Vector3 force;

    public Vector3 oldForce = Vector3.zero;
    public float oldDensity = 1f;
    //public Vector3 oldPosition = Vector3.zero;
    //public Vector3 oldVelocity = Vector3.zero;

    public float particleScale = 1;

    public int timeStep = 0;

    void Start()
    {
        //Probably won't have anything here
    }

    void Update()
    {
        // The paper by Kruger (2006) suggests that particles can be rendered with varying size.
        // This will make it look more natural.
        // The basis for the calculation is that volume = mass / density.
        /*
        if (density != 0.0f)
        {
            float r = Mathf.Pow((3 * mass) / (4 * Mathf.PI * density), 1f / 3f);
            transform.localScale = new Vector3(r, r, r) * 2 * particleScale; // multiply by 2 to get diameter
        }
        else
        {
            transform.localScale = Vector3.one * particleScale;
        }
        */
        transform.localScale = Vector3.one * particleScale; //simplified for testing

        // Update position
        //transform.position = position;

        if (positions != null)
        {
            transform.position = positions[timeStep];
            timeStep++;
            timeStep = timeStep % positions.Count;
        }
    }

    public void positionUpdate()
    {
        transform.position = position;
    }
}
