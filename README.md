# tets2vox
Converts tetrahedral meshes (MSH 2.2) to voxel meshes (MSH 2.2, VoxHenry).

## Usage

Input file format: MSH 2.2 Format tetrahedra mesh.
Output file format: MSH 2.2 if output extension is `.msh` or `.msh2`, VoxHenry if output extension is `.vhr`. 

```
Usage: tets2vox.exe --input <input> --output <output> --res <res>

Options:
  -i, --input <input>    Input file
  -o, --output <output>  Output file
  -r, --res <res>        Resolution
  -h, --help             Print help
  -V, --version          Print version
```
