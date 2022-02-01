using System.Collections;
using System.Collections.Generic;
using System.Threading.Tasks;
using UnityEngine;
using UnityEngine.Profiling;

public class SpatialHash
{
    private float cellSize;
    private Vector3 cellVec;

    private float h;
    private float hSqr;

    private int tableSize;

    private Dictionary<int, List<Particle>> table;

    CustomSampler getBoundingBox;

    public SpatialHash(int particleCount, int cellSize, float h)
    {
        getBoundingBox = CustomSampler.Create("GetBoundingBox");
        // Set up table with 2x the number of particles;
//        this.tableSize = 2 * particleCount;
        this.tableSize = 30*30*10 / (cellSize * cellSize * cellSize);
        this.table = new Dictionary<int, List<Particle>>(tableSize);

        Debug.Log("tableSize: " + tableSize);

        this.h = h;
        this.hSqr = h * h;

        this.cellSize = cellSize;
        this.cellVec = Vector3.one * cellSize;

        for (int i = 0; i < tableSize; i++)
        {
            table.Add(i, new List<Particle>());
        }
    }

    public void Clear()
    {
        table.Clear();
        for (int i = 0; i < tableSize; i++)
        {
            table.Add(i, new List<Particle>());
        }
    }

    public void Insert(List<Particle> particles)
    {
        for (int i = 0; i < particles.Count; i++)
        {
            Particle p = particles[i];
            table[Hash(Discretize(p.position))].Add(p);
        }
    }

    public void Update()
    {
        for (int i = 0; i < tableSize; i++)
        {
//        Parallel.For(0, tableSize, i => {
            List<Particle> cell = table[i];
            for (int j = 0; j < cell.Count; j++)
            {
                Particle p = cell[j];
                int hash = Hash(Discretize(p.position));
                if (hash != i)
                {
                    cell.RemoveAt(j);
                    table[hash].Add(p);
                }
            }
//        });
        }
    }

    public List<Particle> GetNeighbors(Particle p)
    {
        getBoundingBox.Begin();
        List<int> bb = GetBoundingBox(p);
        getBoundingBox.End();
        List<Particle> neighbors = new List<Particle>();
        for (int i = 0; i < bb.Count; i++)
        {
            List<Particle> cell = table[bb[i]];
            for (int j = 0; j < cell.Count; j++) {
                if ((p.position - cell[j].position).sqrMagnitude <= hSqr)
                    neighbors.Add(cell[j]);
            }
        }

        return neighbors;
    }

    public void SetNeighbors()
    {
 //       for (int i = 0; i < tableSize; i++)
//        {
        Parallel.For(0, tableSize, i =>
        {
            List<Particle> cell = table[i];
            if (cell.Count > 0)
            {
                // The bounding box will be the same for all particles in the same cell
                getBoundingBox.Begin();
                List<int> bb = GetBoundingBox(cell[0]);
                getBoundingBox.End();
                for (int j = 0; j < cell.Count; j++)
                {
                    Particle p = cell[j];
                    p.neighbors.Clear();
                    for (int k = 0; k < bb.Count; k++)
                    {
                        List<Particle> bbCell = table[bb[k]];
                        for (int l = 0; l < bbCell.Count; l++)
                        {
                            if (p != bbCell[l] && (p.position - bbCell[l].position).sqrMagnitude <= hSqr)
                                p.neighbors.Add(bbCell[l]);
                        }
                    }
                }
            }
        });
//        }
    }

    private Vector3 Discretize(Vector3 r)
    {
        return new Vector3(
            Mathf.Floor(r.x / cellSize),
            Mathf.Floor(r.y / cellSize),
            Mathf.Floor(r.z / cellSize)
        );
    }

    private List<int> GetBoundingBox(Particle p)
    {
        List<int> cellHashes = new List<int>();
        Vector3 min = Discretize(p.position - cellVec);
        Vector3 max = Discretize(p.position + cellVec);
/*
        if (p.debug)
        {
            Debug.Log("bb min: " + min);
            Debug.Log("bb max: " + max);
        }
*/
//        Debug.Log("bb size: " + ((max.x - min.x) * (max.y - min.y) * (max.z - min.z)));

        // Iterate from min to max
/*
        for (int dx = 0; dx <= max.x - min.x; dx++)
        {
            for (int dy = 0; dy <= max.y - min.y; dy++)
            {
                for (int dz = 0; dz <= max.z - min.z; dz++)
                {
                    int cellHash = Hash(min + new Vector3(dx, dy, dz));
                    cellHashes.Add(cellHash);
                }
            }
        }
*/
        for (int dx = 0; dx <= 2; dx++)
        {
            for (int dy = 0; dy <= 2; dy++)
            {
                for (int dz = 0; dz <= 2; dz++)
                {
                    int cellHash = Hash(min + new Vector3(dx, dy, dz));
                    cellHashes.Add(cellHash);
                }
            }
        }
        return cellHashes;
    }

    // r should already be discretized
    private int Hash(Vector3 r)
    {
        int x = (int) r.x;
        int y = (int) r.y;
        int z = (int) r.z;
        int p1 = 72856093;
        int p2 = 19349663;
        int p3 = 83492791;
        int hash = ((x * p1) ^ (y * p2) ^ (z * p3)) % tableSize;
        if (hash < 0)
            hash = hash + tableSize;
        return hash;
    }
}
