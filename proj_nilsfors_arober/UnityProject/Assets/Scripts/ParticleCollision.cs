using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ParticleCollision : MonoBehaviour
{
    private Particle p;

    void Start()
    {
        p = this.GetComponent<Particle>();
    }

    private void OnCollisionEnter(Collision c)
    {
        //Debug.Log("coll c: " + c.contactCount);
        p.collision = true;
        int cc = c.contactCount;

        do {
            cc--;
            ContactPoint cp = c.GetContact(cc);
            if (Mathf.Abs(p.collisionSeparation) < Mathf.Abs(cp.separation))
            {
                p.collisionPoint = cp.point;
                p.collisionNormal = cp.normal;
                p.collisionSeparation = cp.separation;
            }
        } while (cc > 0);
    }

    private void OnCollisionStay(Collision c)
    {
        //Debug.Log("coll c: " + c.contactCount);
        p.collision = true;
        int cc = c.contactCount;

        do {
            cc--;
            ContactPoint cp = c.GetContact(cc);
            if (Mathf.Abs(p.collisionSeparation) < Mathf.Abs(cp.separation))
            {
                p.collisionPoint = cp.point;
                p.collisionNormal = cp.normal;
                p.collisionSeparation = cp.separation;
            }
        } while (cc > 0);
    }

    private void OnCollisionExit(Collision c)
    {
        p.collision = false;
        p.collisionSeparation = 0.0f;
    }
}
