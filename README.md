# PDBCluster

K-medoids and K-centers clustering of PDB files (RNA/DNA/protein).

## Compiling
```bash
make
```

## Running
```
./PDBCluster -name <atomtype1> [-name <atomtype2> ..] -clusteringType <distSum|maxDistSum|maxDist> -clusters <integer> [list of pdb-files]
```

or: 

```
./PDBCluster -name <atomtype1> [-name <atomtype2> ..] -dumpDistances <filename>
```

## Example

Example of running k-centers for k=10:

```bash
./PDBCluster -name CA -name CB -clusteringType maxDist -clusters 10 sample1.pdb sample2.pdb sample3.pdb ...
```

 * Standard value for -name is C4' and CA
 * Standard value for -clusteringType is distSum
 * Standard value for -clusters is 10

`-clusteringType distSum` will perform standard k-medoids clustering where the
cost associated with a set of medoids is the average distance from each member
to its nearest medoid.

`-clusteringType maxDistSum` will perform a variant of k-centers clustering
where the cost associated with a set of medoids is the average distance from
each medoid to its furthest associated member.

`-clusteringType maxDist` will perform k-centers clustering where the cost
associated with a set of medoids is the largest distance from any site to the
nearest medoid.

For all clustering types, the dRMSD metric is used to indicate the distance
between two structures. If `-dumpDistances` is specified no clustering is
performed, the structure dRMSD distance matrix is just written to the specified
file.
